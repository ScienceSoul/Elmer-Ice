
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
! *  Module containing solvers and routines for the SIA temperature equation
! *
! ******************************************************************************
! *
! *  Author: Hakime Seddik
! *  Derived from the original temperature equation solver
! *
! *****************************************************************************/
!------------------------------------------------------------------------------
   RECURSIVE SUBROUTINE TemperateIceSolverSIA( Model,Solver,Timestep,TransientSimulation )
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
     USE DiffuseConvective
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
     TYPE(Nodes_t) :: ElementNodes
     TYPE(Element_t),POINTER :: Element, ParentElement
     TYPE(Variable_t), POINTER :: TempSol,FlowSol,CurrentSol, MeshSol,VarTempHom,VarTempResidual, DepthVar, &
                                  FreeSurfGrad1, FreeSurfGrad2, NormalVar
     TYPE(ValueList_t), POINTER :: Equation,Material,SolverParams,BodyForce,BC,Constants, ParentMaterial

     INTEGER :: i,j,k,l,m,n,t,iter,body_id,eq_id,material_id, &
          istat, LocalNodes,bf_id, bc_id,  DIM, &
          NSDOFs, NonlinearIter, NonlinearIterMin, Ind(3,3)

     INTEGER, POINTER :: NodeIndexes(:), TempPerm(:),FlowPerm(:),CurrentPerm(:),MeshPerm(:), DepthPerm(:), &
                         FreeSurfGrad1Perm(:), FreeSurfGrad2Perm(:), NormalPerm(:)

     CHARACTER(LEN=MAX_NAME_LEN) :: ConvectionFlag, VariableName, SolverName, FlowSolName

     LOGICAL :: Stabilize = .TRUE., Bubbles = .TRUE., UseBubbles, &
          Found, FluxBC, Permeable=.TRUE., IsPeriodicBC=.FALSE.,&
          AllocationsDone = .FALSE.,  SubroutineVisited = .FALSE., FirstTime=.TRUE.,&
          LimitSolution, ApplyDirichlet, FlowSolutionFound
     LOGICAL, ALLOCATABLE ::  LimitedSolution(:), ActiveNode(:)
     LOGICAL :: strainHeating, IgnoreMeshVelocity, isSliding

     REAL(KIND=dp) :: NonlinearTol, LinearTol, Relax, &
          SaveRelax,dt,CumulativeTime, RelativeChange, &
          Norm,PrevNorm,S,C, &
          ReferencePressure=0.0d0, &
          HeatCapacityGradient(3), round = 0.0D00, vSn, un, rho, gravity, gradh, sigma
     REAL(KIND=dp), POINTER :: Temp(:), FlowSolution(:), &
          ForceVector(:), PrevSolution(:), HC(:), Hwrk(:,:,:),&
          PointerToResidualVector(:),&
          ResidualVector(:), TempHomologous(:), DepthValues(:), FreeSurfGrad1Values(:), FreeSurfGrad2Values(:), NormalValues(:)
     REAL(KIND=dp), ALLOCATABLE :: MASS(:,:), &
       STIFF(:,:), LOAD(:), HeatConductivity(:,:,:), &
         FORCE(:), Pressure(:),  MeshVelocity(:,:),&
         IceVeloU(:),IceVeloV(:),IceVeloW(:),TimeForce(:), &
         TransferCoeff(:), LocalTemp(:), Work(:), C1(:), C0(:), CT(:), Zero(:), Viscosity(:),&
         UpperLimit(:), HeatCapacity(:),  Density(:), TempExt(:), &
         StiffVector(:), OldValues(:), OldRHS(:), normal(:), velo(:), Sn(:), ut(:), Sig(:,:), Arrhenius(:)
     REAL(KIND=dp) :: at,at0,totat,st,totst,t1,CPUTime,RealTime

     SAVE &
          OldValues,             &
          OldRHS,                &
          MeshVelocity,          &
          IceVeloU,             &
          IceVeloV,             &
          IceVeloW,             &
          Pressure,              &
          ElementNodes    ,      &
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
          AllocationsDone, FirstTime, Hwrk, VariableName, SolverName, NonLinearTol, M, round, Ind, normal, velo, Sn, ut, Sig, &
          Arrhenius

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
                IceVeloU,                &
                IceVeloV,                &
                IceVeloW,                &
                Pressure,                 &
                ElementNodes % x,         &
                ElementNodes % y,         &
                ElementNodes % z,         &
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
                ActiveNode,               &
                normal,                   &
                velo,                     &
                Sn,                       &
                ut,                       &
                Sig,                      &
                Arrhenius)              
        END IF                           
        
        ALLOCATE(                                  &
             OldValues( K ), &
             OldRHS( L ), &
             MeshVelocity( 3,N ),                  &
             IceVeloU( N ),                       &
             IceVeloV( N ),                       &
             IceVeloW( N ),                       &
             Pressure( N ),                        &
             ElementNodes % x( N ),                &
             ElementNodes % y( N ),                &
             ElementNodes % z( N ),                &
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
             normal(DIM),                          &
             velo(DIM),                            &
             Sn(DIM),                              &
             ut(DIM),                              &
             Sig(DIM,DIM),                         &
             Arrhenius(N),                         &
             STAT=istat )

        IF ( istat /= 0 ) THEN
           CALL FATAL( SolverName, 'Memory allocation error' )
        ELSE
           CALL INFO(SolverName, 'Memory allocation done', level=1 )
        END IF
        
        AllocationsDone = .TRUE.

        Do i=1, 3
          Ind(i,i) = i
        END DO
        Ind(1,2) = 4
        Ind(2,1) = 4
        Ind(2,3) = 5
        Ind(3,2) = 5
        Ind(3,1) = 6
        Ind(1,3) = 6

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
     UseBubbles = GetLogical( SolverParams,'Bubbles',Found )
     IF ( .NOT.Found .AND. (.NOT.Stabilize)) UseBubbles = .TRUE.

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

              IgnoreMeshVelocity =  GetLogical( Equation, 'Ignore Mesh Velocity', Found)
                 IF (.NOT.Found) IgnoreMeshVelocity = .FALSE.

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


              k = ListGetInteger( Model % Bodies(Element % BodyId) % Values, 'Equation', &
                   minv=1, maxv=Model % NumberOFEquations )

              SELECT CASE( ListGetString( Model % Equations(k) % Values, &
                   'Convection', Found ) )

                 !-----------------
              CASE( 'computed' )
                 !-----------------

                 FlowSolName =  GetString( Model % Equations(k) % Values,'Flow Solution Name', Found)
                 IF(.NOT.Found) THEN        
                    CALL WARN(SolverName,'Keyword >Flow Solution Name< not found in section >Equation<')
                    CALL WARN(SolverName,'Taking default value >Flow Solution<')
                    WRITE(FlowSolName,'(A)') 'Flow Solution'
                 END IF


                 FlowSol => VariableGet( Solver % Mesh % Variables, FlowSolName )
                 IF ( ASSOCIATED( FlowSol ) ) THEN
                    FlowPerm     => FlowSol % Perm
                    NSDOFs       =  FlowSol % DOFs
                    FlowSolution => FlowSol % Values
                    FlowSolutionFound = .TRUE.
                 ELSE
                    CALL INFO(SolverName,'No Flow Solution associated',Level=1)
                    FlowSolutionFound = .FALSE.
                 END IF
              CASE( "none")
                 FlowSolutionFound = .FALSE.

              END SELECT

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
              ! NB.: viscosity needed for strain heating
              !      but Newtonian flow is assumed
              !------------------------------------------
              !Viscosity = 0.0D00
              BodyForce => GetBodyForce()
              IF ( ASSOCIATED( BodyForce ) ) THEN
                 strainHeating =  GetLogical( BodyForce, 'Friction Heat', Found)
                 IF(.NOT. Found) strainHeating = .FALSE.
              ELSE
                  CALL FATAL( SolverName, 'Body Force not found in model.')
              END IF
              
              IF (strainHeating) THEN
                 CALL INFO(SolverName,'Strain heating included in temperature model',Level=15)
                 FreeSurfGrad1 => VariableGet( Model % Variables, 'FreeSurfGrad1' )
                 IF ( ASSOCIATED( FreeSurfGrad1 ) ) THEN
                      FreeSurfGrad1Perm    =>  FreeSurfGrad1 % Perm
                      FreeSurfGrad1Values  =>  FreeSurfGrad1 % Values
                 ELSE
                      CALL FATAL(SolverName,'No variable FreeSurfGrad1 found')  
                 END IF

                 FreeSurfGrad2 => VariableGet( Model % Variables, 'FreeSurfGrad2' )
                 IF ( ASSOCIATED( FreeSurfGrad2 ) ) THEN
                      FreeSurfGrad2Perm    =>  FreeSurfGrad2 % Perm
                      FreeSurfGrad2Values  =>  FreeSurfGrad2 % Values
                 ELSE
                      CALL FATAL(SolverName,'No variable FreeSurfGrad2 found')  
                 END IF

                 DepthVar =>  VariableGet(Model % Variables,'Depth')
                 IF ( ASSOCIATED( DepthVar ) ) THEN
                      DepthPerm   => DepthVar % Perm
                      DepthValues => DepthVar % Values
                 ELSE
                      CALL FATAL(SolverName, 'Depth variable not found')
                 END IF
                 
                 Arrhenius(1:N) = GetReal( Material,'Temperature SIA Arrhenius', found)
                 IF(.NOT. Found) THEN 
                     CALL FATAL(SolverName, 'Temperature SIA Arrhenius not found')
                 END IF

                 body_id = Element % BodyId
                 bf_id = ListGetInteger( Model % Bodies(body_id) % Values, &
                    'Body Force', Found, 1, Model % NumberOfBodyForces )
                 gravity = GetConstReal(Model % BodyForces(bf_id) % Values, 'Flow Bodyforce 3', Found)
                 IF(.NOT. Found) THEN 
                     CALL FATAL(SolverName, 'Flow Body Force 3 not found')
                 END IF
                 gravity = gravity*(-1.0)

                 DO k=1,n
                    gradh =  ( FreeSurfGrad1Values(FreeSurfGrad1Perm(Element % NodeIndexes(k)))**2.0_dp & 
                             + FreeSurfGrad2Values(FreeSurfGrad2Perm(Element % NodeIndexes(k)))**2.0_dp )**0.5_dp
                    sigma =  Density(k) * gravity * DepthValues(DepthPerm(Element % NodeIndexes(k))) * gradh
                    Viscosity(k) = 2.0_dp * Arrhenius(k) * sigma**(3.0_dp+1.0_dp)
                 END DO
              ELSE
                   Viscosity = 0.0D00
              END IF

              NULLIFY(BodyForce)              
              
              !------------------------------------------------------------------------------
              ! Get mesh velocity
              !------------------------------------------------------------------------------
              MeshVelocity = 0.0d0
              IF (.NOT.IgnoreMeshVelocity) &
                 CALL GetVectorLocalSolution( MeshVelocity, 'Mesh Velocity')         
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
              ELSE IF (( ConvectionFlag == 'computed' ) .AND. FlowSolutionFound) THEN
                 DO i=1,n
                    k = FlowPerm(Element % NodeIndexes(i))
                    IF ( k > 0 ) THEN
                       Pressure(i) = FlowSolution(NSDOFs*k) + ReferencePressure
                       SELECT CASE( NSDOFs )
                       CASE(3)
                          IceVeloU(i) = FlowSolution( NSDOFs*k-2 )
                          IceVeloV(i) = FlowSolution( NSDOFs*k-1 )
                          IceVeloW(i) = 0.0D0
                       CASE(4)
                          IceVeloU(i) = FlowSolution( NSDOFs*k-3 )
                          IceVeloV(i) = FlowSolution( NSDOFs*k-2 )
                          IceVeloW(i) = FlowSolution( NSDOFs*k-1 )
                       END SELECT
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
              ! Do we really need residual free Bubbles
              !------------------------------------------------------------------------------
              Bubbles = UseBubbles  .AND. &
                   ( ConvectionFlag == 'computed' .OR. ConvectionFlag == 'constant' )         
              !------------------------------------------------------------------------------
              ! Get element local matrices, and RHS vectors
              !------------------------------------------------------------------------------
              MASS = 0.0d00
              STIFF = 0.0d00
              FORCE = 0.0D00
              ! cartesian coords
              !----------------
              IF ( CurrentCoordinateSystem() == Cartesian ) THEN
                 CALL DiffuseConvectiveComposeSIA( &
                      MASS, STIFF, FORCE, LOAD, &
                      CT(1:N), C0, C1(1:N), HeatConductivity, &
                      .FALSE., Zero, Zero, IceVeloU, IceVeloV, IceVeloW, &
                      MeshVelocity(1,1:N),MeshVelocity(2,1:N),MeshVelocity(3,1:N),&
                      Viscosity, Density, Pressure, Zero, Zero,&
                      .FALSE., Stabilize, Bubbles, Element, n, ElementNodes )
              ! special coords (account for metric)
              !-----------------------------------
              ELSE
                CALL FATAL(SolverName,'Matrix assembly not supported for non cartesian coordinate system.')                  
              END IF              
              !------------------------------------------------------------------------------
              ! If time dependent simulation add mass matrix to stiff matrix
              !------------------------------------------------------------------------------
              TimeForce  = FORCE
              IF ( TransientSimulation ) THEN
                 IF ( Bubbles ) FORCE = 0.0d0
                 CALL Default1stOrderTime( MASS,STIFF,FORCE )
              END IF
              !------------------------------------------------------------------------------
              !  Update global matrices from local matrices
              !------------------------------------------------------------------------------
              IF (  Bubbles ) THEN
                 CALL Condensate( N, STIFF, FORCE, TimeForce )
                 IF (TransientSimulation) CALL DefaultUpdateForce( TimeForce )
              END IF

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

              ParentElement => Element % BoundaryInfo % Left
              IF ( .NOT. ASSOCIATED( ParentElement ) ) &
                       ParentElement => Element % BoundaryInfo % Right
              body_id = ParentElement % BodyId
              material_id = ListGetInteger(Model % Bodies(body_id) % Values, 'Material', Found)
              ParentMaterial => Model % Materials(material_id) % Values
              IF (.NOT. ASSOCIATED(ParentMaterial)) THEN
                  CALL FATAL(SolverName,'Parent material not found')  
              END IF 
               
              !Get BofyForce ID
              !----------------
              bf_id = ListGetInteger( Model % Bodies(body_id) % Values, &
                    'Body Force', Found, 1, Model % NumberOfBodyForces )

              IF ( ASSOCIATED( BC ) ) THEN            
                 ! Check that we are on the correct boundary part!
                 STIFF=0.0D00
                 FORCE=0.0D00
                 MASS=0.0D00
                 LOAD=0.0D00
                 TransferCoeff = 0.0D00
                 TempExt = 0.0D00
                 FluxBC = .FALSE.
                 isSliding = .FALSE.
                 FluxBC =  GetLogical(BC,TRIM(Solver % Variable % Name) // ' Flux BC', Found)
                 isSliding =  GetLogical(BC,TRIM(Solver % Variable % Name) // ' Flux BC Sliding Heat', Found)

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

                    !---------------	
                    ! Heat sliding
                    !--------------- 
                    IF (isSliding) THEN

                      ! Get the surface gradiants
                      !------------------------
                      FreeSurfGrad1 => VariableGet( Model % Variables, 'FreeSurfGrad1' )
                      IF ( ASSOCIATED( FreeSurfGrad1 ) ) THEN
                        FreeSurfGrad1Perm    =>  FreeSurfGrad1 % Perm
                        FreeSurfGrad1Values  =>  FreeSurfGrad1 % Values
                      ELSE
                        CALL FATAL(SolverName,'No variable FreeSurfGrad1 found')  
                      END IF

                      FreeSurfGrad2 => VariableGet( Model % Variables, 'FreeSurfGrad2' )
                      IF ( ASSOCIATED( FreeSurfGrad2 ) ) THEN
                        FreeSurfGrad2Perm    =>  FreeSurfGrad2 % Perm
                        FreeSurfGrad2Values  =>  FreeSurfGrad2 % Values
                      ELSE
                        CALL FATAL(SolverName,'No variable FreeSurfGrad2 found')  
                      END IF

                      DepthVar =>  VariableGet(Model % Variables,'Depth')
                      IF ( ASSOCIATED( DepthVar ) ) THEN
                        DepthPerm   => DepthVar % Perm
                        DepthValues => DepthVar % Values
                      ELSE
                        CALL FATAL(SolverName, 'Depth variable not found')
                      END IF

                      NormalVar =>  VariableGet(Model % Variables,'Normal Vector')
                      IF ( ASSOCIATED( NormalVar ) ) THEN
                         NormalPerm => NormalVar % Perm
                         NormalValues => NormalVar % Values
                      ELSE
                         CALL FATAL(SolverName, 'Normal Vector variable not found')
                      END IF
                      
                      rho = GetConstReal(ParentMaterial, 'Density', Found)
                      IF(.NOT. Found) THEN 
                        CALL FATAL(SolverName, 'Density not found in Material')
                      END IF

                      gravity = GetConstReal(Model % BodyForces(bf_id) % Values, 'Flow Bodyforce 3', Found)
                      IF(.NOT. Found) THEN 
                        CALL FATAL(SolverName, 'Flow Body Force 3 not found')
                      END IF

                      gravity = gravity*(-1.0)

                      DO k=1,n
  
                        DO i=1, DIM
                            normal(i) = NormalValues(DIM*(NormalPerm(Element % NodeIndexes(k))-1) + i)
                            velo(i) = FlowSolution( (DIM+1)*(FlowPerm(Element % NodeIndexes(k))-1) + i )
                        END DO
   
                        !Tengantial velocity
                        un = SUM(velo(1:DIM)*(normal(1:DIM))) 
                        ut = velo(1:DIM)-un*(normal(1:DIM))
   
                        ! Shallow ice Cauchy stress tensor
                        Sig = 0.0_dp
                        IF(DIM == 2) THEN
                            Sig(1,1) = -rho * gravity * DepthValues(DepthPerm(Element % NodeIndexes(k)))
                            Sig(1,2) = -rho * gravity * DepthValues(DepthPerm(Element % NodeIndexes(k))) * &
                                       FreeSurfGrad1Values(FreeSurfGrad1Perm(Element % NodeIndexes(k)))
                            Sig(2,1) = -rho * gravity * DepthValues(DepthPerm(Element % NodeIndexes(k))) * &
                                       FreeSurfGrad1Values(FreeSurfGrad1Perm(Element % NodeIndexes(k)))
                            Sig(2,2) = -rho * gravity * DepthValues(DepthPerm(Element % NodeIndexes(k)))
                        ELSE IF (DIM == 3) THEN
                            Sig(1,1) = -rho * gravity * DepthValues(DepthPerm(Element % NodeIndexes(k)))
                            Sig(1,2) = 0.0_dp
                            Sig(1,3) = -rho * gravity * DepthValues(DepthPerm(Element % NodeIndexes(k))) * &
                                       FreeSurfGrad1Values(FreeSurfGrad1Perm(Element % NodeIndexes(k)))
                            Sig(2,1) = 0.0_dp
                            Sig(2,2) = -rho * gravity * DepthValues(DepthPerm(Element % NodeIndexes(k)))
                            Sig(2,3) = -rho * gravity * DepthValues(DepthPerm(Element % NodeIndexes(k))) * &
                                       FreeSurfGrad2Values(FreeSurfGrad2Perm(Element % NodeIndexes(k)))
                            Sig(3,1) = -rho * gravity * DepthValues(DepthPerm(Element % NodeIndexes(k))) * &
                                       FreeSurfGrad1Values(FreeSurfGrad1Perm(Element % NodeIndexes(k)))
                            Sig(3,2) = -rho * gravity * DepthValues(DepthPerm(Element % NodeIndexes(k))) * &
                                       FreeSurfGrad2Values(FreeSurfGrad2Perm(Element % NodeIndexes(k)))
                            Sig(3,3) = -rho * gravity * DepthValues(DepthPerm(Element % NodeIndexes(k)))
                        ELSE
                            CALL FATAL(SolverName, 'Unsupported dimmension in Cauchy stress tensor calculation')
                        END IF

                        DO i=1, DIM
                            Sn(i) = SUM(Sig(i,1:DIM)*normal(1:DIM)) 
                        END DO  
 
                        vSn = SUM(ut(1:DIM)*Sn(1:DIM))
   
                        LOAD(k) = LOAD(k) - MIN(vSn,0.0d0)
   
                      END DO

                    END IF

                    ! -------------------------------------
                    ! set boundary due to coordinate system
                    ! -------------------------------------
                    IF ( CurrentCoordinateSystem() == Cartesian ) THEN
                       CALL DiffuseConvectiveBoundary( STIFF,FORCE, &
                            LOAD,TransferCoeff,Element,n,ElementNodes )
                    ELSE
                       CALL DiffuseConvectiveGenBoundary(STIFF,FORCE,&
                            LOAD,TransferCoeff,Element,n,ElementNodes ) 
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

     SubroutineVisited = .TRUE.

CONTAINS

!------------------------------------------------------------------------------
!>  Return element local matrices and RHS vector for diffusion-convection
!>  equation: 
!------------------------------------------------------------------------------
   SUBROUTINE DiffuseConvectiveComposeSIA( MassMatrix,StiffMatrix,ForceVector,  &
      LoadVector,NodalCT,NodalC0,NodalC1,NodalC2,PhaseChange,NodalTemperature, &
         Enthalpy,Ux,Uy,Uz,MUx,MUy,MUz,Nodalmu,Nodalrho,NodalPressure, &
            NodaldPressureDt, NodalPressureCoeff, Compressible, Stabilize, &
              UseBubbles, Element,n,Nodes )
!------------------------------------------------------------------------------
!
!  REAL(KIND=dp) :: MassMatrix(:,:)
!     OUTPUT: time derivative coefficient matrix
!
!  REAL(KIND=dp) :: StiffMatrix(:,:)
!     OUTPUT: rest of the equation coefficients
!
!  REAL(KIND=dp) :: ForceVector(:)
!     OUTPUT: RHS vector
!
!  REAL(KIND=dp) :: LoadVector(:)
!     INPUT:
!
!  REAL(KIND=dp) :: NodalCT,NodalC0,NodalC1
!     INPUT: Coefficient of the time derivative term, 0 degree term, and
!            the convection term respectively
!
!  REAL(KIND=dp) :: NodalC2(:,:,:)
!     INPUT: Nodal values of the diffusion term coefficient tensor
!
!  LOGICAL :: PhaseChange
!     INPUT: Do we model phase change here...
!
!  REAL(KIND=dp) :: NodalTemperature
!     INPUT: NodalTemperature from previous iteration
!
!  REAL(KIND=dp) :: Enthalpy
!     INPUT: Enthalpy from previous iteration, needed if we model
!            phase change
!
!  REAL(KIND=dp) :: UX(:),UY(:),UZ(:)
!     INPUT: Nodal values of velocity components from previous iteration
!           used only if coefficient of the convection term (C1) is nonzero
!
!  REAL(KIND=dp) :: Nodalmu(:)
!     INPUT: Nodal values of the viscosity
!
!  LOGICAL :: Stabilize
!     INPUT: Should stabilzation be used ? Used only if coefficient of the
!            convection term (C1) is nonzero
!
!  TYPE(Element_t) :: Element
!       INPUT: Structure describing the element (dimension,nof nodes,
!               interpolation degree, etc...)
!
!  INTEGER :: n
!       INPUT: Number of element nodes
!
!  TYPE(Nodes_t) :: Nodes
!       INPUT: Element node coordinates
!
!------------------------------------------------------------------------------

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

     LOGICAL :: Vms, Found, Transient, stat,Convection,ConvectAndStabilize,Bubbles, &
          FrictionHeat
     TYPE(ValueList_t), POINTER :: BodyForce

!------------------------------------------------------------------------------

     StabilizeFlag = GetString( GetSolverParams(),'Stabilization Method',Found )
     Vms = StabilizeFlag == 'vms'

     Transient = GetString(GetSimulation(),'Simulation type',Found)=='transient'

     dim = CoordinateSystemDimension()
     c = dim + 1

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
        IntegStuff = GaussPoints( element, Element % TYPE % GaussPoints2 )
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
       NodalVelo(1,1:n) = Ux(1:n)
       NodalVelo(2,1:n) = Uy(1:n)
       NodalVelo(3,1:n) = Uz(1:n)
       NodalVelo(dim+1,1:n) = NodalPressure(1:n)

       dNodalBasisdx = 0._dp
       GradNodal = 0._dp
       DO p=1,n
          u = Element % TYPE % NodeU(p)
          v = Element % TYPE % NodeV(p)
          w = Element % TYPE % NodeW(p)
          stat = ElementInfo( Element, Nodes, u,v,w, detJ, Basis, dBasisdx )
 
          dNodalBasisdx(1:n,p,:) = dBasisdx(1:n,:)
          GradNodal(p,1:dim,1:dim) = MATMUL( NodalVelo(1:dim,1:n), dBasisdx(1:n,1:dim) )
          GradNodal(p,dim+1,1:dim) = MATMUL( NodalVelo(dim+1,1:n), dBasisdx(1:n,1:dim) )
       END DO

       NodalPvelo = 0._dp
       IF ( Transient ) THEN
         dt = CurrentModel % Solver % dt
         Order = MIN(CurrentModel % Solver % DoneTime,CurrentModel % Solver % Order)
 
         CALL GetVectorLocalSolution( NodalPVelo, 'Flow Solution', tStep=-1 )
 
         IF ( Order<2 ) THEN
           NodalPVelo(1:dim,1:n)=(NodalVelo(1:dim,1:n)-NodalPVelo(1:dim,1:n))/dt
         ELSE
           CALL GetVectorLocalSolution( Work, 'Flow Solution', tStep=-2 )
           NodalPVelo(1:dim,1:n)=(1.5_dp*NodalVelo(1:dim,1:n) - &
                  2._dp*NodalPVelo(1:dim,1:n)+0.5_dp*Work(1:dim,1:n) )/dt
         END IF
       END IF

       expc = GetCReal(GetMaterial(),'Heat Expansion Coefficient',Found)
       reft = GetCReal(GetMaterial(),'Reference Temperature',Found)
       CALL GetConstRealArray( GetConstants(), gWrk, 'Grav',Found )
       IF ( Found ) THEN
         grav = gWrk(1:3,1)*gWrk(4,1)
       ELSE
         grav    =  0.00_dp
         grav(2) = -9.81_dp
       END IF

       LC1 = 2._dp/mK
       LC(1,1:n) = Element % TYPE % NodeU(1:n)
       LC(2,1:n) = Element % TYPE % NodeV(1:n)
       LC(3,1:n) = Element % TYPE % NodeW(1:n)

       DO i=1,Element % TYPE % DIMENSION
         minl=MINVAL(LC(i,1:n))
         maxl=MAXVAL(LC(i,1:n))
         LC(i,1:n) = 2*(LC(i,1:n)-minl)/(maxl-minl)-1
       END DO
     ELSE IF ( Stabilize .AND. Convection ) THEN
       ConvectAndStabilize = .TRUE.
       hK = element % hK
       mK = element % StabilizationMK
       dNodalBasisdx = 0._dp
       DO p=1,n
         u = Element % TYPE % NodeU(p)
         v = Element % TYPE % NodeV(p)
         w = Element % TYPE % NodeW(p)
         stat = ElementInfo( Element, Nodes, u,v,w, detJ, Basis, dBasisdx )
         dNodalBasisdx(1:n,p,:) = dBasisdx(1:n,:)
       END DO
     END IF

     BodyForce => GetBodyForce()
     FrictionHeat = .FALSE.
     IF (ASSOCIATED(BodyForce))  &
       FrictionHeat = GetLogical( BodyForce, 'Friction Heat', Found )
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
          dEnth = 0.0D0
          dTemp = 0.0D0
          DO i=1,3
            dEnth = dEnth + SUM( Enthalpy(1:n) * dBasisdx(1:n,i) )**2
            dTemp = dTemp + SUM( NodalTemperature(1:n) * dBasisdx(1:n,i) )**2
          END DO

          IF( dTemp > TINY( dTemp ) ) THEN
            CT = SQRT( dEnth/dTemp )
          ELSE
            CALL Info('DiffuseConvectiveCompose',&
                'Temperature difference almost zero, cannot account for phase change!',Level=7)
            CT = 0.0
          END IF
        END IF
!------------------------------------------------------------------------------
!      Coefficient of the diffusion term & it s derivatives at the
!      integration point
!------------------------------------------------------------------------------
       rho = SUM( Nodalrho(1:n) * Basis(1:n) ) 

       DO i=1,dim
         DO j=1,dim
           C2(i,j) = SUM( NodalC2(i,j,1:n) * Basis(1:n) )
         END DO
       END DO

       DO i=1,dim
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
          IF ( PhaseChange ) C1 = CT
!------------------------------------------------------------------------------
!         Velocity from previous iteration at the integration point
!------------------------------------------------------------------------------
          Velo = 0.0D0
          Velo(1) = SUM( (UX(1:n)-MUX(1:n))*Basis(1:n) )
          Velo(2) = SUM( (UY(1:n)-MUY(1:n))*Basis(1:n) )
          IF ( dim > 2 ) Velo(3) = SUM( (UZ(1:n)-MUZ(1:n))*Basis(1:n) )

          IF ( Compressible ) THEN
            Grad = 0.0D0
            DO i=1,3
              Grad(1,i) = SUM( UX(1:n)*dBasisdx(1:n,i) )
              Grad(2,i) = SUM( UY(1:n)*dBasisdx(1:n,i) )
              IF ( dim > 2 ) Grad(3,i) = SUM( UZ(1:n)*dBasisdx(1:n,i) )
            END DO

            Pressure = SUM( NodalPressure(1:n)*Basis(1:n) )
            DivVelo = 0.0D0
            DO i=1,dim
              DivVelo = DivVelo + Grad(i,i)
            END DO
          END IF


          IF ( Vms ) THEN
            mu = GetCReal( GetMaterial(), 'Viscosity', Found )
            mu = EffectiveViscosity( mu, rho, Ux, Uy, Uz, &
                   Element, Nodes, n, n, u,v,w )

            Grad = 0.0D0
            DO i=1,3
              Grad(1,i) = SUM( Ux(1:n)*dBasisdx(1:n,i) )
              Grad(2,i) = SUM( Uy(1:n)*dBasisdx(1:n,i) )
              Grad(3,i) = SUM( Uz(1:n)*dBasisdx(1:n,i) )
            END DO
            VNorm = SQRT( SUM(Velo(1:dim)**2) )

            Temperature = SUM(Basis(1:n)*NodalTemperature(1:n))

            DO i=1,dim
              GradP(i) = SUM( NodalPressure(1:n)*dBasisdx(1:n,i) )
            END DO


            Gmat = 0._dp
            Gvec = 0._dp
            DO i=1,dim
              DO j=1,dim
                Gvec(i) = Gvec(i) + SUM(LC(j,1:n)*dBasisdx(1:n,i))
                DO k=1,dim
                  Gmat(i,j) = Gmat(i,j) + SUM(LC(k,1:n)*dBasisdx(1:n,i)) * &
                                          SUM(LC(k,1:n)*dBasisdx(1:n,j))
                END DO
              END DO
            END DO

            IF ( Transient ) THEN
              Tau_M = 1._dp / SQRT( SUM(Velo*MATMUL(Gmat,Velo)) + &
                    LC1**2 * (mu/rho)**2*SUM(Gmat*Gmat)/dim + 4/dt**2 )
            ELSE
              Tau_M = 1._dp / SQRT( SUM(Velo*MATMUL(Gmat,Velo)) + &
                    LC1**2 * (mu/rho)**2 * SUM(Gmat*Gmat)/dim )
            END IF

            Pe  = MIN( 1.0_dp, mK*hK*C1*VNorm/(2*ABS(C2(1,1))) )
            IF ( VNorm /= 0.0 ) THEN
               Tau = hK * Pe / (2 * C1 * VNorm)
            END IF

            RM = 0._dp
            DO p=1,n
              RM(p) = C0 * Basis(p)
              DO i=1,dim
                RM(p) = RM(p) + C1 * Velo(i) * dBasisdx(p,i)
                DO j=1,dim
                  RM(p) = RM(p) - C2(i,j)*SUM(dNodalBasisdx(p,1:n,i)*dBasisdx(1:n,j))
                END DO
              END DO
            END DO

            VRM = 0._dp
            DO i=1,dim
              VRM(i) = SUM(NodalPVelo(i,1:n)*Basis(1:n))
              DO j=1,dim
                VRM(i) = VRM(i) + Velo(j) * Grad(i,j)
                VRM(i) = VRM(i) - (mu/rho)*SUM(GradNodal(1:n,i,j)*dBasisdx(1:n,j))
              END DO
              VRM(i) = VRM(i) + GradP(i)
              VRM(i) = VRM(i) + Grav(i)*ExpC*( Temperature - RefT )
            END DO
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

            DO i=1,dim
              DO j=1,dim
                DO k=1,dim
                  dC2dx(i,j,k) = SUM( NodalC2(i,j,1:n)*dBasisdx(1:n,k) )
                END DO
              END DO
            END DO

!------------------------------------------------------------------------------
!           Compute residual & stablization vectors
!------------------------------------------------------------------------------
            DO p=1,N
              SU(p) = C0 * Basis(p)
              DO i = 1,dim
                SU(p) = SU(p) + C1 * dBasisdx(p,i) * Velo(i)
                DO j=1,dim
                  SU(p) = SU(p) - dC2dx(i,j,j) * dBasisdx(p,i)
                  SU(p) = SU(p) - C2(i,j) * SUM(dNodalBasisdx(p,1:n,i)*dBasisdx(1:n,j))
                END DO
              END DO

              SW(p) = C0 * Basis(p)
              DO i = 1,dim
                SW(p) = SW(p) + C1 * dBasisdx(p,i) * Velo(i)
                DO j=1,dim
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
!         The diffusion term, SIA only assumes the vertical component
!------------------------------------------------------------------------------
          DO i=dim,dim
            DO j=dim,dim
              A = A + C2(i,j) * dBasisdx(q,i) * dBasisdx(p,j)
            END DO
          END DO

          IF ( Convection ) THEN
!------------------------------------------------------------------------------
!           The convection term
!------------------------------------------------------------------------------
            DO i=1,dim
              A = A + C1 * Velo(i) * dBasisdx(q,i) * Basis(p)
            END DO
!------------------------------------------------------------------------------
!           Next we add the stabilization...
!------------------------------------------------------------------------------
            IF ( Vms ) THEN
              DO i=1,dim
                A = A - C1 * Tau_M * VRM(i) * dBasisdx(q,i) * Basis(p)

                A = A + C1 * Velo(i) * Tau * RM(q) * dBasisdx(p,i)
                M = M + C1 * Velo(i) * Tau * CT*Basis(q) * dBasisdx(p,i)

                A = A - C1 * Tau_M * VRM(i) * Tau * RM(q) * dBasisdx(p,i)
                M = M - C1 * Tau_M * VRM(i) * Tau * CT*Basis(q) * dBasisdx(p,i)
              END DO
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
            DO i=1,dim
              Force = Force + Pcoeff*Velo(i)*SUM(NodalPressure(1:n)*dBasisdx(1:n,i))
            END DO
          END IF

          IF ( FrictionHeat ) THEN
            mu = SUM( Nodalmu(1:n) * Basis(1:n) )
            IF ( mu > 0.0d0 ) THEN
               Force = Force + mu
            END IF
          END IF
        END IF
!------------------------------------------------------------------------------
        DO p=1,NBasis
          Load = Force * Basis(p)
          IF ( Vms ) THEN
            DO i=1,dim
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
   END SUBROUTINE DiffuseConvectiveComposeSIA
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
   END SUBROUTINE TemperateIceSolverSIA
!------------------------------------------------------------------------------



