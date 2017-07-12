
!/*****************************************************************************/! *
! *  Elmer, A Finite Element Software for Multiphysical Problems
! *
! *  Copyright 1st April 1995 - , CSC - Scientific Computing Ltd., Finland
! *
! *  This program is free software; you can redistribute it and/or
! *  modify it under the terms of the GNU General Public License
! *  as published by the Free Software Foundation; either version 2
! *  of the License, or (at your option) any later version.
! *
! *  This program is distributed in the hope that it will be useful,
! *  but WITHOUT ANY WARRANTY; without even the implied warranty of
! *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! *  GNU General Public License for more details.
! *
! *  You should have received a copy of the GNU General Public License
! *  along with this program (in file fem/GPL-2); if not, write to the
! *  Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor,
! *  Boston, MA 02110-1301, USA.
! *
! *****************************************************************************/
!
!/******************************************************************************
! *
! *  Module containing solvers and routines for 1D temperature equation
! *
! ******************************************************************************
! *
! *  Author: Hakime
! *  Derived from the original 3D temperature equation solver
! *
! *****************************************************************************/
!------------------------------------------------------------------------------
   RECURSIVE SUBROUTINE TemperateIce1D( Model,Solver,Timestep,TransientSimulation )
!DLLEXPORT TemprateIceSolver
!------------------------------------------------------------------------------
!******************************************************************************
!
!  Solve the convection diffusion equation with limiters!
!
!  ARGUMENTS:
!
!  TYPE(Model_t) :: Model,  
!     INPUT: All model information (mesh,materials,BCs,etc...)
!
!  TYPE(Solver_t) :: Solver
!     INPUT: Linear equation solver options
!
!  REAL(KIND=dp) :: Timestep
!     INPUT: Timestep size for time dependent simulations
!
!******************************************************************************
     USE DiffuseConvectiveGeneral
     USE Differentials
     USE MaterialModels
!     USE Adaptive
     USE DefUtils

!------------------------------------------------------------------------------
     IMPLICIT NONE
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!    External variables
!------------------------------------------------------------------------------
     TYPE(Model_t)  :: Model
     TYPE(Solver_t), TARGET :: Solver
     LOGICAL :: TransientSimulation
     REAL(KIND=dp) :: Timestep
!------------------------------------------------------------------------------
!    Local variables
!------------------------------------------------------------------------------
     TYPE(Solver_t), POINTER :: PointerToSolver
     TYPE(Matrix_t), POINTER :: Systemmatrix
     TYPE(Nodes_t) :: ElementNodes, BDElementNodes
     TYPE(Element_t),POINTER :: Element, BDElement
     TYPE(Variable_t), POINTER :: TempSol,FlowSol,CurrentSol, MeshSol,VarTempHom,VarTempResidual, DepthSol, HeightSol, asSol, &
                                  gasSol, fieldTempSol, fieldTempHomoSol, wProfileSol, FSSol
     TYPE(ValueList_t), POINTER :: Equation,Material,SolverParams,BodyForce,BC,Constants

     INTEGER :: i,j,k,l,m,n,nbd,t,tt,iter,body_id,eq_id,material_id, &
          istat, LocalNodes,bf_id, bc_id, bc_idd, DIM, &
          NSDOFs, NonlinearIter, NonlinearIterMin, DepthDOFs, HeightDOFs, asDOFs, gasDOFs, fieldTempDOFs, fieldTempHomoDOFs,&
          FSDOFs

     INTEGER, POINTER :: NodeIndexes(:), TempPerm(:),FlowPerm(:),CurrentPerm(:),MeshPerm(:), DepthPerm(:), HeightPerm(:), &
                         asPerm(:), gasPerm(:), fieldTempPerm(:), fieldTempHomoPerm(:), wProfilePerm(:), FSPerm(:)

     CHARACTER(LEN=MAX_NAME_LEN) :: ConvectionFlag, VariableName, SolverName

     LOGICAL :: Stabilize = .TRUE., &
          Found, FluxBC, Permeable=.TRUE., IsPeriodicBC=.FALSE.,&
          AllocationsDone = .FALSE.,  SubroutineVisited = .FALSE., FirstTime=.TRUE.,&
          LimitSolution, ApplyDirichlet, FlowSolutionFound
     LOGICAL, ALLOCATABLE ::  LimitedSolution(:), ActiveNode(:)
     LOGICAL :: strainHeating, foundNodes

     REAL(KIND=dp) :: NonlinearTol, LinearTol, Relax, &
          SaveRelax,dt,CumulativeTime, RelativeChange, &
          Norm,PrevNorm,S,C, &
          ReferencePressure=0.0d0, &
          HeatCapacityGradient(3), round = 0.0D00
     REAL(KIND=dp), POINTER :: Temp(:), FlowSolution(:), &
          ForceVector(:), PrevSolution(:), HC(:), Hwrk(:,:,:),&
          PointerToResidualVector(:),&
          ResidualVector(:), TempHomologous(:), DepthValues(:), HeightValues(:), asValues(:), gasValues(:), fieldTempValues(:), &
          fieldTempHomoValues(:), wProfileValues(:), FSValues(:)
     REAL(KIND=dp), ALLOCATABLE :: MASS(:,:), &
       STIFF(:,:), LOAD(:), HeatConductivity(:,:,:), &
         FORCE(:), Pressure(:),  MeshVelocity(:,:),&
         IceVeloU(:),IceVeloV(:),IceVeloW(:),TimeForce(:), &
         TransferCoeff(:), LocalTemp(:), Work(:), C1(:), C0(:), CT(:), Zero(:), Viscosity(:),&
         UpperLimit(:), HeatCapacity(:),  Density(:), TempExt(:), &
         StiffVector(:), OldValues(:), OldRHS(:)
     REAL(KIND=dp) :: at,at0,totat,st,totst,t1,CPUTime,RealTime
     REAL(KIND=dp) :: z_tildes, c_1, c_2, c_3, c_4, zstar, vz_b, as, iceThickness

     SAVE &
          OldValues,             &
          OldRHS,                &
          MeshVelocity,          &
          IceVeloU,             &
          IceVeloV,             &
          IceVeloW,             &
          Pressure,              &
          ElementNodes,          &
          BDElementNodes,        &
          Work,Zero,             &
          Viscosity,             &
          HeatCapacity,          &
          Density,               &
          TempExt,               &
          C1,                    &
          C0,                    &
          CT,                    &
          TransferCoeff,         &
          LocalTemp,             &
          HeatConductivity,      &
          MASS,                  &
          STIFF,LOAD,            &
          FORCE,                 &
          TimeForce,             &
          StiffVector,           &
          ResidualVector,        &
          UpperLimit,            &
          LimitedSolution,       &
          ActiveNode,            &
          AllocationsDone, FirstTime, Hwrk, VariableName, SolverName, NonLinearTol, M, round

!------------------------------------------------------------------------------
!    Get variables needed for solution
!------------------------------------------------------------------------------
     DIM = CoordinateSystemDimension()
     SolverName = 'TemperateIceSolver ('// TRIM(Solver % Variable % Name) // ')'
     VariableName = TRIM(Solver % Variable % Name)

     IF ( .NOT. ASSOCIATED( Solver % Matrix ) ) RETURN
     SystemMatrix => Solver % Matrix
     ForceVector => Solver % Matrix % RHS

     PointerToSolver => Solver

     TempSol => Solver % Variable
     TempPerm  => TempSol % Perm
     Temp => TempSol % Values
     
     LocalNodes = COUNT( TempPerm > 0 )
     IF ( LocalNodes <= 0 ) RETURN

!------------------------------------------------------------------------------
!    Allocate some permanent storage, this is done first time only
!------------------------------------------------------------------------------
     IF ( .NOT. AllocationsDone .OR. Solver % Mesh % Changed ) THEN
        N = Solver % Mesh % MaxElementNodes
        M = Model % Mesh % NumberOfNodes
        K = SIZE( SystemMatrix % Values )
        L = SIZE( SystemMatrix % RHS )

        IF ( AllocationsDone ) THEN
           DEALLOCATE(                    &
                OldValues,                &
                OldRHS,                   &
                MeshVelocity,             &
                IceVeloU,                 &
                IceVeloV,                 &
                IceVeloW,                 &
                Pressure,                 &
                ElementNodes % x,         &
                ElementNodes % y,         &
                ElementNodes % z,         &
                BDElementNodes % x,       &
                BDElementNodes % y,       &
                BDElementNodes % z,       &
                Work,Zero,                &
                Viscosity,                &
                HeatCapacity,             &
                Density,                  &
                TempExt,                  &
                C1,                       &
                C0,                       &
                CT,                       &
                TransferCoeff,            &
                LocalTemp,                &
                HeatConductivity,         &
                MASS,                     &
                STIFF,LOAD,               &
                FORCE,                    &
                TimeForce,                &
                StiffVector,              &
                ResidualVector,           &
                UpperLimit,               &
                LimitedSolution,          &
                ActiveNode)              
        END IF                           
        
        ALLOCATE(                                  &
             OldValues( K ), &
             OldRHS( L ), &
             MeshVelocity( 3,N ),                  &
             IceVeloU( N ),                        &
             IceVeloV( N ),                        &
             IceVeloW( N ),                        &
             Pressure( N ),                        &
             ElementNodes % x( N ),                &
             ElementNodes % y( N ),                &
             ElementNodes % z( N ),                &
             BDElementNodes % x ( N ),             &
             BDElementNodes % y ( N ),             &
             BDElementNodes % z ( N ),             &
             Work( N ), Zero( N ),                 &
             Viscosity( N ),                       &
             HeatCapacity( M ),                    &
             Density( N ),                         &
             TempExt( N ),                         &
             C1( N ),                              &
             C0( N ),                              &
             CT( N ),                              &
             TransferCoeff( N ),                   &
             LocalTemp( N ),                       &
             HeatConductivity( 3,3,N ),            &
             MASS(  2*N,2*N ),                     &
             STIFF( 2*N,2*N ),LOAD( N ),           &
             FORCE( 2*N ),                         &
             TimeForce( 2*N ),                     &
             StiffVector( M ),                     &
             ResidualVector(K),                    &
             UpperLimit( M ),                      &
             LimitedSolution( M ),                 &
             ActiveNode( M ),                      &
             STAT=istat )

        IF ( istat /= 0 ) THEN
           CALL FATAL( SolverName, 'Memory allocation error' )
        ELSE
           CALL INFO(SolverName, 'Memory allocation done', level=1 )
        END IF
        
        AllocationsDone = .TRUE.
     
     END IF

!------------------------------------------------------------------------------
!    Say hello
!------------------------------------------------------------------------------
     WRITE(Message,'(A,A)')&
          'Limited diffusion Solver for variable ', VariableName
     CALL INFO(SolverName,Message,Level=1)

!------------------------------------------------------------------------------
!    Read physical and numerical constants and initialize 
!------------------------------------------------------------------------------
     Constants => GetConstants()
     SolverParams => GetSolverParams()

     Stabilize = GetLogical( SolverParams,'Stabilize',Found )
     IF (.NOT. Found) Stabilize = .FALSE.

     LinearTol = GetConstReal( SolverParams, &
          'Linear System Convergence Tolerance',    Found )
     IF ( .NOT.Found ) THEN
        CALL FATAL(SolverName, 'No >Linear System Convergence Tolerance< found')
     END IF

     NonlinearIter = GetInteger(   SolverParams, &
                     'Nonlinear System Max Iterations', Found )
     IF ( .NOT.Found ) NonlinearIter = 2

     NonlinearIterMin = GetInteger(   SolverParams, &
                     'Nonlinear System Min Iterations', Found )
     IF ( .NOT.Found ) NonlinearIterMin = 2
     IF (NonlinearIterMin < 2) THEN
        NonlinearIterMin = 2
        CALL Info( SolverName, 'Given Nonlinear System Min Iterations smaller than 2. Needs and will use 2', Level=2 )
    END IF

     NonlinearTol  = GetConstReal( SolverParams, &
          'Nonlinear System Convergence Tolerance',    Found )

     Relax = GetConstReal( SolverParams, &
               'Nonlinear System Relaxation Factor',Found )

     IF ( .NOT.Found ) Relax = 1.0D00

     ApplyDirichlet = GetLogical( SolverParams, &
          'Apply Dirichlet', Found)
     IF ( .NOT.Found ) THEN
        ApplyDirichlet = .FALSE.
     END IF

     SaveRelax = Relax
     dt = Timestep
     CumulativeTime = 0.0d0
     ActiveNode = .FALSE.

!------------------------------------------------------------------------------
!   time stepping loop.
!------------------------------------------------------------------------------
     DO WHILE( CumulativeTime < Timestep-1.0d-12 .OR. .NOT. TransientSimulation )
        round = round +1.0D00
!------------------------------------------------------------------------------
!       The first time around this has been done by the caller...
!------------------------------------------------------------------------------
        IF ( TransientSimulation .AND. .NOT.FirstTime ) THEN
           CALL InitializeTimestep( Solver )
        END IF

!------------------------------------------------------------------------------
!       Save current solution
!------------------------------------------------------------------------------
        ALLOCATE( PrevSolution(LocalNodes) )
        PrevSolution = Temp(1:LocalNodes)

        totat = 0.0d0
        totst = 0.0d0

!------------------------------------------------------------------------------
!       Get externally declared DOFs
!------------------------------------------------------------------------------
        IF (.NOT.ApplyDirichlet) ActiveNode = .FALSE.
        VarTempHom => VariableGet( Model % Mesh % Variables, TRIM(Solver % Variable % Name) // ' Homologous' )
        IF (.NOT.ASSOCIATED(VarTempHom)) THEN
           WRITE(Message,'(A)') TRIM(Solver % Variable % Name) // ' Homologous not associated'
           CALL FATAL( SolverName, Message)
        END IF        

        VarTempResidual => VariableGet( Model % Mesh % Variables, TRIM(Solver % Variable % Name) // ' Residual' )
        IF (.NOT.ASSOCIATED(VarTempResidual)) THEN
           WRITE(Message,'(A)') '>' // TRIM(Solver % Variable % Name) // ' Residual< not associated'
           CALL FATAL( SolverName, Message)
        END IF
        PointerToResidualVector => VarTempResidual % Values

        DepthSol => VariableGet( Solver % Mesh % Variables, 'Depth' )
        IF ( ASSOCIATED( DepthSol ) ) THEN
            DepthPerm => DepthSol % Perm
            DepthValues => DepthSol % Values
            DepthDOFs = DepthSol % DOFs
        ELSE
            WRITE(Message,'(a,I6)') 'Could not find Depth field variable'
            CALL FATAL('InitAge',Message)
        END IF

        HeightSol => VariableGet( Solver % Mesh % Variables, 'Height' )
        IF ( ASSOCIATED( HeightSol ) ) THEN
            HeightPerm => HeightSol % Perm
            HeightValues => HeightSol % Values
            HeightDOFs = HeightSol % DOFs
        ELSE
            WRITE(Message,'(a,I6)') 'Could not find Height field variable'
            CALL FATAL('InitAge',Message)
        END IF

        asSol => VariableGet( Solver % Mesh % Variables, 'as' )
        IF ( ASSOCIATED( asSol ) ) THEN
            asPerm   => asSol % Perm
            asValues => asSol % Values
            asDOFs   =  asSol % DOFs
        ELSE
            WRITE(Message,'(a,I6)') 'Could not find as field variable'
            CALL FATAL('InitAge',Message)
        END IF

        gasSol => VariableGet( Solver % Mesh % Variables, 'gas' )
        IF ( ASSOCIATED( gasSol ) ) THEN
            gasPerm   => gasSol % Perm
            gasValues => gasSol % Values
            gasDOFs   =  gasSol % DOFs
        ELSE
            WRITE(Message,'(a,I6)') 'Could not find gas field variable'
            CALL FATAL('InitAge',Message)
        END IF

        wProfileSol => VariableGet( Solver % Mesh % Variables, 'wprofile' )
        IF ( ASSOCIATED( wProfileSol ) ) THEN
            wProfilePerm   => wProfileSol % Perm
            wProfileValues => wProfileSol % Values
        ELSE
            WRITE(Message,'(a,I6)') 'Could not find wprofile field variable'
            CALL FATAL('InitAge',Message)
        END IF

        FSSol => VariableGet( Solver % Mesh % Variables, 'FS' )
        IF ( ASSOCIATED( FSSol ) ) THEN
            FSPerm   => FSSol % Perm
            FSValues => FSSol % Values
            FSDOFs   =  FSSol % DOFs
        ELSE
            WRITE(Message,'(a,I6)') 'Could not find FS field variable'
            CALL FATAL('InitAge',Message)
        END IF

        ! Look at the node on the surface belonging to the same column and get the accumulation
        CALL Info( SolverName, 'Accumulation for each column of ice...', Level=2 )
        DO t=1,Solver % NumberOfActiveElements
            Element => GetActiveElement(t)
            n = GetElementNOFNodes()
            CALL GetElementNodes( ElementNodes )
            
            DO i=1,n

                foundNodes = .FALSE.
                DO tt=1, Solver % Mesh % NumberOfBoundaryElements
                    BDElement => GetBoundaryElement(tt)
                    IF ( .NOT.ActiveBoundaryElement(BDElement) ) CYCLE
                    nbd = GetElementNOFNodes(BDElement)
                    IF ( GetElementFamily(BDElement) == 1 ) CYCLE
                    bc_id = GetBCId( BDElement )
                    CALL GetElementNodes( BDElementNodes, BDElement )
                    IF (bc_id == 6) THEN
                        DO j=1,nbd
                            IF (ElementNodes % x(i) == BDElementNodes % x(j) .AND. & 
                                    ElementNodes % y(i) == BDElementNodes % y(j)) THEN
                                gasValues(gasPerm(Element % NodeIndexes(i))) = & 
                                        asValues(asPerm(BDElement % NodeIndexes(j)))
                                foundNodes = .True.
                                EXIT
                            END IF
                        END DO
                    END IF
                    IF (foundNodes) EXIT
                END DO

            END DO
        END DO
       
        zstar = 0.25 !(1.0/3.0)
        vz_b = -2.5 * 10.0**(-3.0) !0.0_dp

!------------------------------------------------------------------------------
!       non-linear system iteration loop
!------------------------------------------------------------------------------
        DO iter=1,NonlinearIter
           FirstTime = .FALSE.           
           !------------------------------------------------------------------------------
           ! print out some information
           !------------------------------------------------------------------------------
           at  = CPUTime()
           at0 = RealTime()

           CALL Info( SolverName, ' ', Level=4 )
           CALL Info( SolverName, ' ', Level=4 )
           CALL Info( SolverName, '-------------------------------------',Level=4 )
           WRITE( Message,'(A,A,I3,A,I3)') &
                TRIM(Solver % Variable % Name),  ' iteration no.', iter,' of ',NonlinearIter
           CALL Info( SolverName, Message, Level=4 )
           CALL Info( SolverName, '-------------------------------------',Level=4 )
           CALL Info( SolverName, ' ', Level=4 )
           CALL Info( SolverName, 'Starting Assembly...', Level=4 )
             
           !------------------------------------------------------------------------------
           ! lets start
           !------------------------------------------------------------------------------
           CALL DefaultInitialize()
           !-----------------------------------------------------------------------------
           ! Get lower and Upper limit:
           !-----------------------------------------------------------------------------
           DO t=1,Solver % NumberOfActiveElements
              Element => GetActiveElement(t)
              n = GetElementNOFNodes()
              CALL GetElementNodes( ElementNodes )
              Material => GetMaterial()
              ! upper limit
              !------------
              UpperLimit(Element % Nodeindexes(1:N)) = ListGetReal(Material,TRIM(VariableName) // & 
                   ' Upper Limit',n,Element % NodeIndexes, Found)
              IF (.NOT. Found) THEN
                 LimitedSolution(Element % Nodeindexes(1:N)) = .FALSE.
                 WRITE(Message,'(a,i10)') 'No upper limit of solution for element no. ', t
                 CALL INFO(SolverName, Message, level=10)
              ELSE
                 LimitedSolution(Element % Nodeindexes(1:N)) = .TRUE.
              END IF
           END DO
           !------------------------------------------------------------------------------
           ! write some info on max/min values
           !------------------------------------------------------------------------------
           WRITE(Message,'(a,e13.6,a,e13.6)') &
                'Max/min values Temperature:', MAXVAL( Temp(:)),'/',MINVAL( Temp(:))
           CALL INFO(SolverName,Message,Level=4)
           !-----------------------------------------------------------------------------
           body_id = -1
           NULLIFY(Material)
           !------------------------------------------------------------------------------
           ! Bulk elements
           !------------------------------------------------------------------------------
           DO t=1,Solver % NumberOfActiveElements
              !------------------------------------------------------------------------------
              ! write some info on status of assembly
              !------------------------------------------------------------------------------
              IF ( RealTime() - at0 > 1.0 ) THEN
                 WRITE(Message,'(a,i3,a)' ) '   Assembly: ', INT(100.0 - 100.0 * &
                      (Solver % NumberOfActiveElements-t) / &
                      (1.0*Solver % NumberOfActiveElements)), ' % done'

                 CALL Info( SolverName, Message, Level=5 )

                 at0 = RealTime()
              END IF
              !------------------------------------------------------------------------------
              ! Check if this element belongs to a body where scalar equation
              ! should be calculated
              !------------------------------------------------------------------------------
              Element => GetActiveElement(t,Solver)
              IF (.NOT.ASSOCIATED(Element)) CYCLE
              IF ( Element % BodyId /= body_id ) THEN
                 Equation => GetEquation()
                 IF (.NOT.ASSOCIATED(Equation)) THEN
                    WRITE (Message,'(A,I3)') 'No Equation  found for boundary element no. ', t
                    CALL FATAL(SolverName,Message)
                 END IF

              ConvectionFlag = GetString( Equation, 'Convection', Found )

              Material => GetMaterial()
                 IF (.NOT.ASSOCIATED(Material)) THEN
                    WRITE (Message,'(A,I3)') 'No Material found for boundary element no. ', t
                    CALL FATAL(SolverName,Message)
                 ELSE
                    material_id = GetMaterialId( Element, Found)
                    IF(.NOT.Found) THEN
                       WRITE (Message,'(A,I3)') 'No Material ID found for boundary element no. ', t
                       CALL FATAL(SolverName,Message)
                    END IF
                 END IF
              END IF
            
              !------------------------------------------------------------------------------
              ! Get element material parameters
              !------------------------------------------------------------------------------              
              N = GetElementNOFNodes(Element)
              CALL GetElementNodes( ElementNodes )
              CALL ListGetRealArray( Material,TRIM(Solver % Variable % Name) // &
                   ' Heat Conductivity',Hwrk,n, Element % NodeIndexes )
              HeatConductivity = 0.0d0
              IF ( SIZE(Hwrk,1) == 1 ) THEN
                 DO i=1,3
                    HeatConductivity( i,i,1:N ) = Hwrk( 1,1,1:N)
                 END DO
              ELSE IF ( SIZE(Hwrk,2) == 1 ) THEN
                 DO i=1,MIN(3,SIZE(Hwrk,1))
                    HeatConductivity(i,i,1:N) = Hwrk(i,1,1:N)
                 END DO
              ELSE
                 DO i=1,MIN(3,SIZE(Hwrk,1))
                    DO j=1,MIN(3,SIZE(Hwrk,2))
                       HeatConductivity(i,j,1:N) = Hwrk(i,j,1:N)
                    END DO
                 END DO
              END IF              
              HeatCapacity(1:N) =  ListGetReal( Material,  TRIM(Solver % Variable % Name) // &
                   ' Heat Capacity', n, Element % NodeIndexes, Found )
              IF (.NOT.Found) THEN
                 HeatCapacity = 0.0D00
                 WRITE(Message,'(a,a,a,i5,a,i5,a)') 'Keyword >', TRIM(Solver % Variable % Name) // &
                   ' Heat Capacity', '< not found for element ', t, ' material ', material_id
                 CALL INFO(SolverName,Message,Level=4)
              END IF

              Density(1:N) = ListGetReal( Material, 'Density',  N, Element % NodeIndexes, Found )
              IF (.NOT.Found) THEN
                 Density = 0.0D00
                 WRITE(Message,'(a,i5,a,i5,a)') 'Keyword >Density< not found for element ',&
                      t, ' material ', material_id
                 CALL INFO(SolverName,Message,Level=4)
              END IF

              !------------------------------------------
              ! Zero strain heating
              !------------------------------------------
              Viscosity = 0.0D00
            
              !------------------------------------------------------------------------------
              ! No mesh velocity
              !------------------------------------------------------------------------------
              MeshVelocity = 0.0d0
  
              !------------------------------------------------------------------------------         
              ! asuming convection or ALE mesh contribution by default
              !------------------------------------------------------------------------------         
              C1(1:N) = Density(1:N) * HeatCapacity(1:N)
              CT(1:N) = Density(1:N) * HeatCapacity(1:N)
              
              !------------------------------------------------------------------------------
              ! Get scalar velocity
              !------------------------------------------------------------------------------         
              IceVeloU = 0.0d00
              IceVeloV = 0.0d00
              IceVeloW = 0.0d00
              ! constant (i.e., in section Material given) velocity
              !---------------------------------------------------
              IF ( ConvectionFlag == 'constant' ) THEN
                 IceVeloU(1:N)= GetReal( Material, 'Convection Velocity 1', Found )
                 IceVeloV(1:N) = GetReal( Material, 'Convection Velocity 2', Found )
                 IceVeloW(1:N) = GetReal( Material, 'Convection Velocity 3', Found )                 
              ! computed velocity
              !------------------
              ELSE IF ( ConvectionFlag == 'computed' ) THEN
                 DO i=1,n

                    IceVeloU(i) = 0.0_dp
                    IceVeloV(i) = 0.0_dp
                
                    k = Element % NodeIndexes(i)
                    iceThickness = DepthValues(DepthPerm(k))+HeightValues(HeightPerm(k))
                    z_tildes = HeightValues(HeightPerm(k))/iceThickness

                    c_1 = (2.0*(1.0+vz_b))/(2.0-zstar)
                    c_2 = (zstar+2.0*vz_b)/(2.0-zstar)
                    c_3 = (1.0+vz_b)/(zstar*(2.0-zstar))
                    c_4 = -vz_b

                    IF (z_tildes >= zstar) THEN 
                        IceVeloW(i) = (-c_1*z_tildes + c_2) * gasValues(gasPerm(k)) 
                        wProfileValues(wProfilePerm(k)) = (-c_1*z_tildes + c_2) * gasValues(gasPerm(k))
                    ELSE IF (z_tildes <= zstar) THEN
                        IceVeloW(i) = (-c_3*z_tildes**2.0 - c_4) * gasValues(gasPerm(k))
                        wProfileValues(wProfilePerm(k)) = (-c_3*z_tildes**2.0 - c_4) * gasValues(gasPerm(k))
                    END IF

                 END DO
                 WRITE(Message,'(a,i5, a, i5)') 'Convection in element ', t, &
                      ' material ',  material_id
              ELSE  ! Conduction and ALE contribution only
                 IF (ANY( MeshVelocity /= 0.0d0 )) THEN
                    WRITE(Message,'(a,i5, a, i5)') 'Only mesh deformation in element ', t,&
                         ' material ',  material_id
                 ELSE ! neither convection nor ALE mesh deformation contribution -> all C1(1:N)=0
                    C1 = 0.0D0 
                    WRITE(Message,'(a,i5, a, i5)') 'No convection and mesh deformation in element ', t,&
                         ' material ',  material_id
                 END IF                 
              END IF
              CALL INFO(SolverName,Message,Level=10)
              !------------------------------------------------------------------------------
              ! no contribution proportional to temperature by default
              !------------------------------------------------------------------------------
              C0=0.0d00
              !------------------------------------------------------------------------------
              ! Add body forces
              !------------------------------------------------------------------------------
              LOAD = 0.0D00
              BodyForce => GetBodyForce()
              IF ( ASSOCIATED( BodyForce ) ) THEN
                 bf_id = GetBodyForceId()
                 LOAD(1:N) = LOAD(1:N) +   &
                      GetReal( BodyForce, TRIM(Solver % Variable % Name) // ' Volume Source', Found )
              END IF
  
              !------------------------------------------------------------------------------
              ! dummy input array for faking   heat capacity, density, temperature, 
              !                                enthalpy and viscosity
              !------------------------------------------------------------------------------
              Work = 1.0d00
              Zero = 0.0D00

              !------------------------------------------------------------------------------
              ! Get element local matrices, and RHS vectors
              !------------------------------------------------------------------------------
              MASS = 0.0d00
              STIFF = 0.0d00
              FORCE = 0.0D00
              ! cartesian coords
              !----------------
              IF ( CurrentCoordinateSystem() == Cartesian ) THEN
                 CALL DiffuseConvectiveCompose( &
                      MASS, STIFF, FORCE, LOAD, &
                      CT(1:N), C0, C1(1:N), HeatConductivity, &
                      .FALSE., Zero, Zero, IceVeloU, IceVeloV, IceVeloW, &
                      MeshVelocity(1,1:N),MeshVelocity(2,1:N),MeshVelocity(3,1:N),&
                      Viscosity, Density, Pressure, Zero, Zero,&
                      .FALSE., Stabilize, .FALSE., Element, n, ElementNodes )
              ! special coords (account for metric)
              !-----------------------------------
              ELSE
                 CALL FATAL(SolverName, 'Non cartesian system not supported.')                 
              END IF              
              !------------------------------------------------------------------------------
              ! If time dependent simulation add mass matrix to stiff matrix
              !------------------------------------------------------------------------------
              TimeForce  = FORCE
              IF ( TransientSimulation ) THEN
                 CALL Default1stOrderTime( MASS,STIFF,FORCE )
              END IF
              !------------------------------------------------------------------------------
              !  Update global matrices from local matrices
              !------------------------------------------------------------------------------
              CALL DefaultUpdateEquations( STIFF, FORCE )
           END DO     !  Bulk elements

           ! This was introduced to enable the computation of loads in the new system 
           CALL DefaultFinishBulkAssembly()
      
           !------------------------------------------------------------------------------
           ! Neumann & Newton boundary conditions
           !------------------------------------------------------------------------------
           DO t=1, Solver % Mesh % NumberOfBoundaryElements

              ! get element information
              Element => GetBoundaryElement(t)
              IF ( .NOT.ActiveBoundaryElement() ) CYCLE
              n = GetElementNOFNodes()
              IF ( GetElementFamily() == 1 ) CYCLE
              BC => GetBC()
              bc_id = GetBCId( Element )
              CALL GetElementNodes( ElementNodes )

              IF ( ASSOCIATED( BC ) ) THEN            
                 ! Check that we are on the correct boundary part!
                 STIFF=0.0D00
                 FORCE=0.0D00
                 MASS=0.0D00
                 LOAD=0.0D00
                 TransferCoeff = 0.0D00
                 TempExt = 0.0D00
                 FluxBC = .FALSE.
                 FluxBC =  GetLogical(BC,TRIM(Solver % Variable % Name) // ' Flux BC', Found)

                 IF (FluxBC) THEN
                    !------------------------------
                    !BC: -k@T/@n = \alpha(T - TempExt)
                    !------------------------------
                    TransferCoeff(1:N) = GetReal( BC, TRIM(Solver % Variable % Name) //  ' Transfer Coefficient',Found )
                    IF ( ANY(TransferCoeff(1:N) /= 0.0d0) ) THEN
                       TempExt(1:N) = GetReal( BC, TRIM(Solver % Variable % Name) // ' External Value',Found )   
                       DO j=1,n
                          LOAD(j) = LOAD(j) +  TransferCoeff(j) * TempExt(j)
                       END DO
                    END IF
                    !---------------
                    !BC: -k@T/@n = q
                    !---------------
                    LOAD(1:N)  = LOAD(1:N) + &
                         GetReal( BC, TRIM(Solver % Variable % Name) // ' Heat Flux', Found )

                    ! -------------------------------------
                    ! set boundary due to coordinate system
                    ! -------------------------------------
                    IF ( CurrentCoordinateSystem() == Cartesian ) THEN
                       CALL DiffuseConvectiveBoundary( STIFF,FORCE, &
                            LOAD,TransferCoeff,Element,n,ElementNodes )
                    ELSE
                       CALL FATAL(SolverName, 'Non cartesian system not supported.')                        
                    END IF
                 END IF
              END IF

              !------------------------------------------------------------------------------
              ! Update global matrices from local matrices
              !------------------------------------------------------------------------------
              IF ( TransientSimulation ) THEN
                 MASS = 0.d0
                 CALL Default1stOrderTime( MASS, STIFF, FORCE )
              END IF
          
              CALL DefaultUpdateEquations( STIFF, FORCE )
           END DO   ! Neumann & Newton BCs
           !------------------------------------------------------------------------------

           CALL DefaultFinishAssembly()
           CALL DefaultDirichletBCs()


           OldValues = SystemMatrix % Values
           OldRHS = ForceVector

           !------------------------------------------------------------------------------
           ! Dirichlet method - matrix and force-vector manipulation
           !------------------------------------------------------------------------------
           IF (ApplyDirichlet) THEN
              ! manipulation of the matrix
              !---------------------------
              DO i=1,Model % Mesh % NumberOfNodes
                 k = TempPerm(i)           
                 IF (ActiveNode(i) .AND. (k > 0)) THEN
                    CALL ZeroRow( SystemMatrix, k ) 
                    CALL SetMatrixElement( SystemMatrix, k, k, 1.0d0 )
                    SystemMatrix % RHS(k) = UpperLimit(i)
                 END IF
              END DO
           END IF

           CALL Info( SolverName, 'Assembly done', Level=4 )

           !------------------------------------------------------------------------------
           !     Solve the system and check for convergence
           !------------------------------------------------------------------------------
           at = CPUTime() - at
           st = CPUTime()

           PrevNorm = Solver % Variable % Norm

           Norm = DefaultSolve()

           st = CPUTIme()-st
           totat = totat + at
           totst = totst + st
           WRITE(Message,'(a,i4,a,F8.2,F8.2)') 'iter: ',iter,' Assembly: (s)', at, totat
           CALL Info( SolverName, Message, Level=4 )
           WRITE(Message,'(a,i4,a,F8.2,F8.2)') 'iter: ',iter,' Solve:    (s)', st, totst
           CALL Info( SolverName, Message, Level=4 )
   
           RelativeChange = Solver % Variable % NonlinChange

           !IF ( PrevNorm + Norm /= 0.0d0 ) THEN
           !   RelativeChange = 2.0d0 * ABS( PrevNorm-Norm ) / (PrevNorm + Norm)
           !ELSE
           !   RelativeChange = 0.0d0
           !END IF

           WRITE( Message, * ) 'Result Norm   : ',Norm
           CALL Info( SolverName, Message, Level=4 )
           WRITE( Message, * ) 'Relative Change : ',RelativeChange
           CALL Info( SolverName, Message, Level=4 )

           SystemMatrix % Values = OldValues
           ForceVector = OldRHS

           !------------------------------------------------------------------------------
           ! compute residual
           !------------------------------------------------------------------------------ 
           IF ( ParEnv % PEs > 1 ) THEN !!!!!!!!!!!!!!!!!!!!!! we have a parallel run
             CALL ParallelInitSolve( SystemMatrix, Temp, ForceVector, ResidualVector )
             CALL ParallelMatrixVector( SystemMatrix, Temp, StiffVector, .TRUE. )
             ResidualVector =  StiffVector - ForceVector
             CALL ParallelSumVector( SystemMatrix, ResidualVector )
           ELSE !!!!!!!!!!!!!!!!!!!!!! serial run 
             CALL CRS_MatrixVectorMultiply( SystemMatrix, Temp, StiffVector)
             ResidualVector =  StiffVector - ForceVector
           END IF

           !-----------------------------
           ! determine "active" nodes set
           !-----------------------------
           IF (ASSOCIATED(VarTempHom)) THEN
              TempHomologous => VarTempHom % Values
              DO i=1,Model % Mesh % NumberOfNodes
                 k = VarTempHom % Perm(i)
                 l= TempPerm(i)
                 TempHomologous(k) = Temp(l) - UpperLimit(i)
                 IF (ApplyDirichlet) THEN
                    !---------------------------------------------------------
                    ! if upper limit is exceeded, manipulate matrix in any case
                    !----------------------------------------------------------
                    IF (TempHomologous(k) >= 0.0 ) THEN
                       ActiveNode(i) = .TRUE.
                       TempHomologous(k) = LinearTol
                    END IF
                    !---------------------------------------------------
                    ! if there is "heating", don't manipulate the matrix
                    !---------------------------------------------------
                    IF (ResidualVector(l) > - LinearTol &
                         .AND. iter>1) ActiveNode(i) = .FALSE.
                 END IF
                 IF( .NOT.ActiveNode(i) ) THEN
                    PointerToResidualVector(VarTempResidual % Perm(i)) = 0.0D00
                 ELSE
                    PointerToResidualVector(VarTempResidual % Perm(i)) = ResidualVector(l)
                 END IF
              END DO
           ELSE
              WRITE(Message,'(A)') TRIM(Solver % Variable % Name) // ' Homologous not associated'
              CALL FATAL( SolverName, Message)
           END IF
           !------------------------------------------
           ! special treatment for periodic boundaries
           !------------------------------------------
           k=0
           DO t=1, Solver % Mesh % NumberOfBoundaryElements

              ! get element information
              Element => GetBoundaryElement(t)
              IF ( .NOT.ActiveBoundaryElement() ) CYCLE
              n = GetElementNOFNodes()
              IF ( GetElementFamily() == 1 ) CYCLE
              BC => GetBC()
              bc_id = GetBCId( Element )
              CALL GetElementNodes( ElementNodes )


              IF ( ASSOCIATED( BC ) ) THEN    
                 IsPeriodicBC = GetLogical(BC,'Periodic BC ' // TRIM(Solver % Variable % Name),Found)
                 IF (.NOT.Found) IsPeriodicBC = .FALSE.
                 IF (IsPeriodicBC) THEN 
                    DO i=1,N
                       IF  (ActiveNode(Element % NodeIndexes(i))) THEN
                          k = k + 1
                          ActiveNode(Element % NodeIndexes(i)) = .FALSE.
                       END IF
                    END DO
                 END IF
              END IF
           END DO
           !----------------------
           ! check for convergence
           !----------------------
           IF ( RelativeChange < NonlinearTol .AND. iter > NonlinearIterMin ) THEN
              EXIT
           ELSE
              IF (ApplyDirichlet) THEN
                 WRITE(Message,'(a,i10)') 'Deactivated Periodic BC nodes:', k
                 CALL INFO(SolverName,Message,Level=1)
                 WRITE(Message,'(a,i10)') 'Number of constrained points:', COUNT(ActiveNode)
                 CALL INFO(SolverName,Message,Level=1)
              END IF
           END IF
        END DO ! of the nonlinear iteration
        !------------------------------------------------------------------------------

        !------------------------------------------------------------------------------
        !   Compute cumulative time done by now and time remaining
        !------------------------------------------------------------------------------
        IF ( .NOT. TransientSimulation ) EXIT
        CumulativeTime = CumulativeTime + dt
        dt = Timestep - CumulativeTime
     END DO ! time interval
     !------------------------------------------------------------------------------

     DEALLOCATE( PrevSolution )

     !! Copy the 1D solution to our Temperature variable
     CALL Info( SolverName, 'Copy solution to Temperature field variable...', Level=4 )
     
     fieldTempSol => VariableGet( Solver % Mesh % Variables, 'Temp' )
        IF ( ASSOCIATED( fieldTempSol ) ) THEN
            fieldTempPerm   => fieldTempSol % Perm
            fieldTempValues => fieldTempSol % Values
            fieldTempDOFs   =  fieldTempSol % DOFs
        ELSE
            WRITE(Message,'(a,I6)') 'Could not find Temp field variable'
            CALL FATAL('InitAge',Message)
     END IF

     fieldTempHomoSol => VariableGet( Solver % Mesh % Variables, 'Temp Homologous' )
        IF ( ASSOCIATED( fieldTempHomoSol ) ) THEN
            fieldTempHomoPerm   => fieldTempHomoSol % Perm
            fieldTempHomoValues => fieldTempHomoSol % Values
            fieldTempHomoDOFs   =  fieldTempHomoSol % DOFs
        ELSE
            WRITE(Message,'(a,I6)') 'Could not find Temp Homologous field variable'
            CALL FATAL('InitAge',Message)
     END IF

     DO i=1,Solver % NumberOFActiveElements
        Element => GetActiveElement(i)
        NodeIndexes => Element % NodeIndexes

        DO k=1, GetElementNOFNodes(Element)
             j = fieldTempPerm(NodeIndexes(k))
             l = fieldTempHomoPerm(NodeIndexes(k))

             fieldTempValues(fieldTempDOFs*(j-1)+1) = &
                        TempSol % Values(1*(TempSol % Perm(NodeIndexes(k))-1)+1)

             fieldTempHomoValues(fieldTempHomoDOFs*(l-1)+1) = &
                        VarTempHom % Values(1*(VarTempHom % Perm(NodeIndexes(k))-1)+1)
        END DO
     END DO

     SubroutineVisited = .TRUE.

 CONTAINS

!------------------------------------------------------------------------------
!>  Return element local matrices and RHS vector for diffusion-convection
!>  equation: 
!------------------------------------------------------------------------------
   SUBROUTINE DiffuseConvectiveCompose( MassMatrix,StiffMatrix,ForceVector,  &
      LoadVector,NodalCT,NodalC0,NodalC1,NodalC2,PhaseChange,NodalTemperature, &
         Enthalpy,Ux,Uy,Uz,MUx,MUy,MUz,Nodalmu,Nodalrho,NodalPressure, &
            NodaldPressureDt, NodalPressureCoeff, Compressible, Stabilize, &
              UseBubbles, Element,n,Nodes )

     REAL(KIND=dp), DIMENSION(:)   :: ForceVector,UX,UY,UZ,MUX,MUY,MUZ,LoadVector
     REAL(KIND=dp), DIMENSION(:,:) :: MassMatrix,StiffMatrix
     REAL(KIND=dp) :: NodalTemperature(:),Enthalpy(:),Nodalmu(:), &
       NodaldPressureDt(:),NodalPressure(:),NodalPressureCoeff(:),Nodalrho(:)
     REAL(KIND=dp) :: NodalC0(:),NodalC1(:),NodalCT(:),NodalC2(:,:,:)

     LOGICAL :: UseBubbles,PhaseChange,Compressible,Stabilize

     INTEGER :: n

     TYPE(Nodes_t) :: Nodes
     TYPE(Element_t), POINTER :: Element

!------------------------------------------------------------------------------
!    Local variables
!------------------------------------------------------------------------------

     CHARACTER(LEN=MAX_NAME_LEN) :: StabilizeFlag
     REAL(KIND=dp) :: dBasisdx(2*n,3),detJ
     REAL(KIND=dp) :: Basis(2*n)
     REAL(KIND=dp) :: ddBasisddx(n,3,3),dNodalBasisdx(n,n,3)

     REAL(KIND=dp) :: Velo(3),Grad(3,3),Force

     REAL(KIND=dp) :: A,M
     REAL(KIND=dp) :: Load

     REAL(KIND=dp) :: VNorm,hK,mK
     REAL(KIND=dp) :: Lambda=1.0,Pe,Pe1,Pe2,Tau,x,y,z

     REAL(KIND=dp) :: Tau_M, Tau_C, Gmat(3,3),Gvec(3), dt=0._dp, LC1, NodalVelo(4,n), &
       RM(n),LC(3,n), gradP(n), VRM(3), GradNodal(n,3,3), PVelo(3), NodalPVelo(4,n), &
       Work(3,n), Grav(3), expc, reft, temperature

     REAL(KIND=dp), POINTER :: gWrk(:,:)

     INTEGER :: i,j,k,c,p,q,t,dim,N_Integ,NBasis,Order

     TYPE(GaussIntegrationPoints_t), TARGET :: IntegStuff
     REAL(KIND=dp) :: s,u,v,w,dEnth,dTemp,mu,DivVelo,Pressure,rho,&
                      Pcoeff, minl, maxl

     REAL(KIND=dp) :: C0,C00,C1,CT,C2(3,3),dC2dx(3,3,3),SU(n),SW(n)

     REAL(KIND=dp), DIMENSION(:), POINTER :: U_Integ,V_Integ,W_Integ,S_Integ

     LOGICAL :: Vms, Found, Transient, stat,Convection,ConvectAndStabilize,Bubbles
     TYPE(ValueList_t), POINTER :: BodyForce

!------------------------------------------------------------------------------

     StabilizeFlag = GetString( GetSolverParams(),'Stabilization Method',Found )
     Vms = .FALSE.

     Transient = GetString(GetSimulation(),'Simulation type',Found)=='transient'

     dim = 3

     ForceVector = 0.0D0
     StiffMatrix = 0.0D0
     MassMatrix  = 0.0D0
     Load = 0.0D0

     Convection =  ANY( NodalC1 /= 0.0d0 )
     NBasis = n
     Bubbles = .FALSE.
     IF ( Convection .AND. .NOT. (Vms .OR. Stabilize) .AND. UseBubbles ) THEN
        NBasis = 2*n
        Bubbles = .TRUE.
     END IF

!------------------------------------------------------------------------------
!    Integration stuff
!------------------------------------------------------------------------------
     IF ( Bubbles ) THEN
        IntegStuff = GaussPoints( element, Element % Type % GaussPoints2 )
     ELSE
        IntegStuff = GaussPoints( element )
     END IF
     U_Integ => IntegStuff % u
     V_Integ => IntegStuff % v
     W_Integ => IntegStuff % w
     S_Integ => IntegStuff % s
     N_Integ =  IntegStuff % n

!------------------------------------------------------------------------------
!    Stabilization parameters: hK, mK (take a look at Franca et.al.)
!    If there is no convection term we don t need stabilization.
!------------------------------------------------------------------------------
     hK = element % hK
     mK = element % StabilizationMK

     ConvectAndStabilize = .FALSE.
     IF  ( Vms ) THEN

         ! No ops!!!!

     ELSE IF ( Stabilize .AND. Convection ) THEN
       ConvectAndStabilize = .TRUE.
       hK = element % hK
       mK = element % StabilizationMK
       dNodalBasisdx = 0._dp
       DO p=1,n
         u = Element % Type % NodeU(p)
         v = Element % Type % NodeV(p)
         w = Element % Type % NodeW(p)
         stat = ElementInfo( Element, Nodes, u,v,w, detJ, Basis, dBasisdx )
         dNodalBasisdx(1:n,p,:) = dBasisdx(1:n,:)
       END DO
     END IF

     BodyForce => GetBodyForce()

!------------------------------------------------------------------------------
!    Now we start integrating
!------------------------------------------------------------------------------
     DO t=1,N_Integ

       u = U_Integ(t)
       v = V_Integ(t)
       w = W_Integ(t)

!------------------------------------------------------------------------------
!      Basis function values & derivatives at the integration point
!------------------------------------------------------------------------------
       stat = ElementInfo( Element,Nodes,u,v,w,detJ, &
             Basis,dBasisdx, Bubbles=Bubbles )

       s = detJ * S_Integ(t)
!------------------------------------------------------------------------------
!      Coefficient of the convection and time derivative terms
!      at the integration point
!------------------------------------------------------------------------------
       C0 = SUM( NodalC0(1:n) * Basis(1:n) )
       C1 = SUM( NodalC1(1:n) * Basis(1:n) )
       CT = SUM( NodalCT(1:n) * Basis(1:n) )
!------------------------------------------------------------------------------
!       Compute effective heatcapacity, if modelling phase change,
!       at the integration point.
!       NOTE: This is for heat equation only, not generally for diff.conv. equ.
!------------------------------------------------------------------------------
        IF ( PhaseChange ) THEN
            ! No ops!!!
        END IF
!------------------------------------------------------------------------------
!      Coefficient of the diffusion term & it s derivatives at the
!      integration point
!------------------------------------------------------------------------------
       rho = SUM( Nodalrho(1:n) * Basis(1:n) ) 

       DO i=3,dim
         DO j=3,dim
           C2(i,j) = SUM( NodalC2(i,j,1:n) * Basis(1:n) )
         END DO
       END DO

       DO i=3,dim
          C2(i,i) = EffectiveConductivity( C2(i,i), rho, Element, &
                 NodalTemperature, UX,UY,UZ, Nodes, n, n, u, v, w )
       END DO
!------------------------------------------------------------------------------
!      If there's no convection term we don't need the velocities, and
!      also no need for stabilization
!------------------------------------------------------------------------------
       Convection = .FALSE.
       IF ( C1 /= 0.0D0 ) THEN
          Convection = .TRUE.
          ! IF ( PhaseChange ) C1 = CT
!------------------------------------------------------------------------------
!         Velocity from previous iteration at the integration point
!------------------------------------------------------------------------------
          Velo = 0.0D0
          Velo(1) = SUM( (UX(1:n)-MUX(1:n))*Basis(1:n) )
          Velo(2) = SUM( (UY(1:n)-MUY(1:n))*Basis(1:n) )
          IF ( dim > 2 ) Velo(3) = SUM( (UZ(1:n)-MUZ(1:n))*Basis(1:n) )

          IF ( Compressible ) THEN
             ! No ops!!!
          END IF

          IF ( Vms ) THEN

             ! No ops!!!

          ELSE IF ( Stabilize ) THEN
!------------------------------------------------------------------------------
!           Stabilization parameter Tau
!------------------------------------------------------------------------------
            VNorm = SQRT( SUM(Velo(1:dim)**2) )

            Pe  = MIN( 1.0D0, mK*hK*C1*VNorm/(2*ABS(C2(1,1))) )
            Tau = 0.0D0
            IF ( VNorm /= 0.0 ) THEN
               Tau = hK * Pe / (2 * C1 * VNorm)
            END IF
!------------------------------------------------------------------------------

            DO i=3,dim
              DO j=3,dim
                DO k=3,dim
                  dC2dx(i,j,k) = SUM( NodalC2(i,j,1:n)*dBasisdx(1:n,k) )
                END DO
              END DO
            END DO

!------------------------------------------------------------------------------
!           Compute residual & stablization vectors
!------------------------------------------------------------------------------
            DO p=1,N
              SU(p) = C0 * Basis(p)
              DO i = 3,dim
                SU(p) = SU(p) + C1 * dBasisdx(p,i) * Velo(i)
                DO j=3,dim
                  SU(p) = SU(p) - dC2dx(i,j,j) * dBasisdx(p,i)
                  SU(p) = SU(p) - C2(i,j) * SUM(dNodalBasisdx(p,1:n,i)*dBasisdx(1:n,j))
                END DO
              END DO

              SW(p) = C0 * Basis(p)
              DO i = 3,dim
                SW(p) = SW(p) + C1 * dBasisdx(p,i) * Velo(i)
                DO j=3,dim
                  SW(p) = SW(p) - dC2dx(i,j,j) * dBasisdx(p,i)
                  SW(p) = SW(p) - C2(i,j) * SUM(dNodalBasisdx(p,1:n,i)*dBasisdx(1:n,j))
                END DO
              END DO
            END DO
          END IF
        END IF

!------------------------------------------------------------------------------
!       Loop over basis functions of both unknowns and weights
!------------------------------------------------------------------------------
        DO p=1,NBasis
        DO q=1,NBasis
!------------------------------------------------------------------------------
!         The diffusive-convective equation without stabilization
!------------------------------------------------------------------------------
          M = CT * Basis(q) * Basis(p)
          A = C0 * Basis(q) * Basis(p)
!------------------------------------------------------------------------------
!         The diffusion term
!------------------------------------------------------------------------------
          DO i=3,dim
            DO j=3,dim
              A = A + C2(i,j) * dBasisdx(q,i) * dBasisdx(p,j)
            END DO
          END DO

          IF ( Convection ) THEN
!------------------------------------------------------------------------------
!           The convection term
!------------------------------------------------------------------------------
            DO i=3,dim
              A = A + C1 * Velo(i) * dBasisdx(q,i) * Basis(p)
            END DO
!------------------------------------------------------------------------------
!           Next we add the stabilization...
!------------------------------------------------------------------------------
            IF ( Vms ) THEN

               ! No ops!!!

            ELSE IF ( Stabilize ) THEN
              A = A + Tau * SU(q) * SW(p)
              M = M + Tau * CT * Basis(q) * SW(p)
            END IF
          END IF

          StiffMatrix(p,q) = StiffMatrix(p,q) + s * A
          MassMatrix(p,q)  = MassMatrix(p,q)  + s * M
        END DO
        END DO

!------------------------------------------------------------------------------
!       The righthand side...
!------------------------------------------------------------------------------
!       Force at the integration point
!------------------------------------------------------------------------------
        Force = SUM( LoadVector(1:n)*Basis(1:n) ) + &
          JouleHeat( Element, Nodes, u, v, w, n )

        IF ( Convection ) THEN
!         IF ( Compressible ) Force = Force - Pressure * DivVelo

          Pcoeff = SUM(NodalPressureCoeff(1:n)*Basis(1:n))
          IF ( Pcoeff /= 0.0_dp ) THEN
            Force = Force + Pcoeff * SUM(NodalDPressureDt(1:n)*Basis(1:n))
            DO i=3,dim
              Force = Force + Pcoeff*Velo(i)*SUM(NodalPressure(1:n)*dBasisdx(1:n,i))
            END DO
          END IF

        END IF
!------------------------------------------------------------------------------
        DO p=1,NBasis
          Load = Force * Basis(p)
          IF ( Vms ) THEN
            DO i=3,dim
              Load = Load + C1 * Velo(i) * Tau  * Force * dBasisdx(p,i)
              Load = Load - C1 * Tau_M * VRM(i) * Tau * Force * dBasisdx(p,i)
            END DO
          ELSE IF ( ConvectAndStabilize ) THEN
            Load = Load + Tau * Force * SW(p)
          END IF
          ForceVector(p) = ForceVector(p) + s * Load
        END DO
      END DO
!------------------------------------------------------------------------------
   END SUBROUTINE DiffuseConvectiveCompose
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!>  Return element local matrices and RSH vector for boundary conditions
!>  of diffusion convection equation: 
!------------------------------------------------------------------------------
   SUBROUTINE DiffuseConvectiveBoundary( BoundaryMatrix,BoundaryVector, &
               LoadVector,NodalAlpha,Element,n,Nodes )
!------------------------------------------------------------------------------

     REAL(KIND=dp) :: BoundaryMatrix(:,:),BoundaryVector(:), &
                    LoadVector(:),NodalAlpha(:)

     TYPE(Nodes_t)   :: Nodes
     TYPE(Element_t) :: Element

     INTEGER :: n

     REAL(KIND=dp) :: ddBasisddx(n,3,3)
     REAL(KIND=dp) :: Basis(n)
     REAL(KIND=dp) :: dBasisdx(n,3),detJ

     REAL(KIND=dp) :: U,V,W,S
     REAL(KIND=dp) :: Force,Alpha
     REAL(KIND=dp), POINTER :: U_Integ(:),V_Integ(:),W_Integ(:),S_Integ(:)

     INTEGER :: i,t,q,p,N_Integ

     TYPE(GaussIntegrationPoints_t), TARGET :: IntegStuff

     LOGICAL :: stat
!------------------------------------------------------------------------------

     BoundaryVector = 0.0D0
     BoundaryMatrix = 0.0D0
!------------------------------------------------------------------------------
!    Integration stuff
!------------------------------------------------------------------------------
     IntegStuff = GaussPoints( Element )
     U_Integ => IntegStuff % u
     V_Integ => IntegStuff % v
     W_Integ => IntegStuff % w
     S_Integ => IntegStuff % s
     N_Integ =  IntegStuff % n

!------------------------------------------------------------------------------
!   Now we start integrating
!------------------------------------------------------------------------------
     DO t=1,N_Integ
       U = U_Integ(t)
       V = V_Integ(t)
       W = W_Integ(t)
!------------------------------------------------------------------------------
!     Basis function values & derivatives at the integration point
!------------------------------------------------------------------------------
       stat = ElementInfo( Element,Nodes,u,v,w,detJ, &
                  Basis,dBasisdx )

       S = detJ * S_Integ(t)
!------------------------------------------------------------------------------
       Force = SUM( LoadVector(1:n)*Basis )
       Alpha = SUM( NodalAlpha(1:n)*Basis )

       DO p=1,N
         DO q=1,N
           BoundaryMatrix(p,q) = BoundaryMatrix(p,q) + &
              s * Alpha * Basis(q) * Basis(p)
         END DO
       END DO

       DO q=1,N
         BoundaryVector(q) = BoundaryVector(q) + s * Basis(q) * Force
       END DO
     END DO
   END SUBROUTINE DiffuseConvectiveBoundary

!------------------------------------------------------------------------------
   END SUBROUTINE TemperateIce1D
!------------------------------------------------------------------------------



