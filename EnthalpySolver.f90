
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
! *  Module containing solvers and routines for standard ice dynamic problems
! *
! ******************************************************************************
! *
! *  Authors: Thomas Zwinger,  Juha Ruokolainen, Hakime Seddik
! *  Email:   Thomas.Zwinger@csc.fi, Juha.Ruokolainen@csc.fi 
! *  Web:     http://www.csc.fi/elmer
! *  Address: CSC - Scientific Computing Ltd.
! *           Keilaranta 14
! *           02101 Espoo, Finland 
! *
! *  Original Date: 14 May 2007
! *
! *****************************************************************************/
!------------------------------------------------------------------------------
   RECURSIVE SUBROUTINE EnthalpySolver( Model,Solver,Timestep,TransientSimulation )
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
     TYPE(Element_t),POINTER :: Element
     TYPE(Variable_t), POINTER :: EnthalpySol,FlowSol,MeshSol
     TYPE(ValueList_t), POINTER :: Equation,Material,SolverParams,BodyForce,BC,Constants

     INTEGER :: i,j,k,l,m,n,t,iter,body_id,eq_id,material_id, &
          istat, LocalNodes,bf_id, bc_id,  DIM, &
          NSDOFs, NonlinearIter
     INTEGER, POINTER :: NodeIndexes(:), EnthalpyPerm(:),FlowPerm(:),MeshPerm(:)

     CHARACTER(LEN=MAX_NAME_LEN) :: ConvectionFlag, VariableName, SolverName, FlowSolName

     LOGICAL :: Stabilize = .TRUE., Bubbles = .TRUE., UseBubbles, &
          Found, FluxBC, Permeable=.TRUE., IsPeriodicBC=.FALSE.,&
          AllocationsDone = .FALSE.,  SubroutineVisited = .FALSE., FirstTime=.TRUE.,&
          FlowSolutionFound
     LOGICAL :: strainHeating, IgnoreMeshVelocity

     REAL(KIND=dp) :: NonlinearTol, LinearTol, Relax, &
          SaveRelax,dt,CumulativeTime, RelativeChange, &
          Norm,PrevNorm,S,C, LatentHeat, EnthalpyRange, &
          ReferencePressure=0.0d0, VelocityScaling = 1.0D00,&
          round = 0.0D00
     REAL(KIND=dp), POINTER :: Enthalpy(:), FlowSolution(:), &
          ForceVector(:), PrevSolution(:), HC(:)
     REAL(KIND=dp), ALLOCATABLE :: & 
          MASS(:,:), STIFF(:,:), LOAD(:),&
          HeatDiffusivity(:,:,:),WaterDiffusivity(:,:,:),EnthalpyDiffusivity(:,:,:), &
          FORCE(:), Pressure(:),  MeshVelocity(:,:),&
          IceVeloU(:),IceVeloV(:),IceVeloW(:),TimeForce(:), &
          TransferCoeff(:), C1(:), C0(:), Zero(:), Unity(:),Viscosity(:),&
          IceHeatCapacity(:),  IceHeatConductivity(:), H2OHeatCapacity(:), H2ODiffusivity(:),&
          EffectiveHeatCapacity(:), Density(:), EnthalpyExt(:), PressureMeltingPoint(:), &
          StiffVector(:), OldValues(:), OldRHS(:), EnthalpyH2O(:), EnthalpyIce(:), &
          LocalEnthalpy(:), LocalEnthalpyFlux(:)
     REAL(KIND=dp) :: at,at0,totat,st,totst,t1,CPUTime,RealTime,WaterContent

     SAVE &
          OldValues,             &
          OldRHS,                &
          MeshVelocity,          &
          IceVeloU,              &
          IceVeloV,              &
          IceVeloW,              &
          Pressure,              &
          ElementNodes    ,      &
          Zero,Unity,            &
          Viscosity,             &
          EffectiveHeatCapacity, &
          IceHeatCapacity,       &
          IceHeatConductivity,   &
          H2OHeatCapacity,       &
          H2ODiffusivity,        &
          Density,               &
          EnthalpyExt,           &
          C1,                    &
          C0,                    &
          TransferCoeff,         &
          EnthalpyH2O,           &
          EnthalpyIce,           &
          LocalEnthalpy,         &
          LocalEnthalpyFlux,     &
          EnthalpyDiffusivity,   &
          HeatDiffusivity,       &
          WaterDiffusivity,      &
          MASS,                  &
          STIFF,LOAD,            &
          FORCE,                 &
          TimeForce,             &
          PressureMeltingPoint,  &
          AllocationsDone,       &
          FirstTime,             &
          VariableName,          &
          SolverName,            &
          NonLinearTol,          &
          M,                     &
          round

!------------------------------------------------------------------------------
!    Get variables needed for solution
!------------------------------------------------------------------------------
     DIM = CoordinateSystemDimension()
     SolverName = 'EnthalpySolver'
     VariableName = TRIM(Solver % Variable % Name)

     IF ( .NOT. ASSOCIATED( Solver % Matrix ) ) RETURN
     SystemMatrix => Solver % Matrix
     ForceVector => Solver % Matrix % RHS

     PointerToSolver => Solver

     EnthalpySol => Solver % Variable
     EnthalpyPerm  => EnthalpySol % Perm
     Enthalpy => EnthalpySol % Values
     
     LocalNodes = COUNT( EnthalpyPerm > 0 )
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
                Zero,Unity,               &
                Viscosity,                &
                EffectiveHeatCapacity,    &
                IceHeatCapacity,          &
                IceHeatConductivity,      &
                H2OHeatCapacity,          &
                H2ODiffusivity,           &
                Density,                  &
                EnthalpyExt,              &
                C1,                       &
                C0,                       &
                TransferCoeff,            &
                EnthalpyH2O,              &
                EnthalpyIce,              &
                LocalEnthalpy,            &
                LocalEnthalpyFlux,        &
                EnthalpyDiffusivity,      &
                HeatDiffusivity,          &
                WaterDiffusivity,         &
                MASS,                     &
                STIFF,LOAD,               &
                FORCE,                    &
                TimeForce,                &
                PressureMeltingPoint)
        END IF                           
        
        ALLOCATE(                                  &
             OldValues( K ),                       &
             OldRHS( L ),                          &
             MeshVelocity( 3,N ),                  &
             IceVeloU( N ),                        &
             IceVeloV( N ),                        &
             IceVeloW( N ),                        &
             Pressure( N ),                        &
             ElementNodes % x( N ),                &
             ElementNodes % y( N ),                &
             ElementNodes % z( N ),                &
             Zero( N ), Unity( N ),                &
             Viscosity( N ),                       &
             EffectiveHeatCapacity(N),             &
             IceHeatCapacity(N),                   &
             IceHeatConductivity(N),               &
             H2OHeatCapacity(N),                   &
             H2ODiffusivity(N),                    &
             Density( N ),                         &
             EnthalpyExt( N ),                     &
             C1( N ),                              &
             C0( N ),                              &
             TransferCoeff( N ),                   &
             EnthalpyH2O( N ),                     &
             EnthalpyIce( N ),                     &
             LocalEnthalpy( N ),                   &
             LocalEnthalpyFlux(N),                 &
             EnthalpyDiffusivity( 3,3,N ),         &
             HeatDiffusivity( 3,3,N ),             &
             WaterDiffusivity( 3,3,N ),            &
             MASS(  2*N,2*N ),                     &
             STIFF( 2*N,2*N ),LOAD( N ),           &
             FORCE( 2*N ),                         &
             TimeForce( 2*N ),                     &
             PressureMeltingPoint( N ),            &
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
          'Enthalpy Solver for variable ', VariableName
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
     IF ( .NOT.Found ) NonlinearIter = 1


     NonlinearTol  = GetConstReal( SolverParams, &
          'Nonlinear System Convergence Tolerance',    Found )

     Relax = GetConstReal( SolverParams, &
               'Nonlinear System Relaxation Factor',Found )

     IF ( .NOT.Found ) Relax = 1.0D00

  

     SaveRelax = Relax
     dt = Timestep
     CumulativeTime = 0.0d0

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
        PrevSolution = Enthalpy(1:LocalNodes)

        totat = 0.0d0
        totst = 0.0d0

        
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
           END DO
           !------------------------------------------------------------------------------
           ! write some info on max/min values
           !------------------------------------------------------------------------------
           WRITE(Message,'(a,e13.6,a,e13.6)') &
                'Max/min values Enthalpy:', MAXVAL( Enthalpy(:)),'/',MINVAL( Enthalpy(:))
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
              IceHeatCapacity(1:N) =  ListGetReal( Material,  &
                   'Ice Heat Capacity', n, Element % NodeIndexes, Found )
              IF (.NOT.Found) THEN
                 IceHeatCapacity = 0.0D00
                 WRITE(Message,'(a,i5,a,i5,a)') 'Keyword >Ice Heat Capacity< not found for element ',&
                      t, ' material ', material_id
                 CALL INFO(SolverName,Message,Level=4)
              END IF
              H2OHeatCapacity(1:N) =  ListGetReal( Material,  &
                   'Water Heat Capacity', n, Element % NodeIndexes, Found )
              IF (.NOT.Found) THEN
                 H2OHeatCapacity = 0.0D00
                 WRITE(Message,'(a,i5,a,i5,a)') 'Keyword >H2O Heat Capacity< not found for element ',&
                      t, ' material ', material_id
                 CALL INFO(SolverName,Message,Level=4)
              END IF
              Density(1:N) = ListGetReal( Material, 'Mixture Density',  N, Element % NodeIndexes, Found )
              IF (.NOT.Found) THEN
                 Density = 0.0D00
                 WRITE(Message,'(a,i5,a,i5,a)') 'Keyword >Mixture Density< not found for element ',&
                      t, ' material ', material_id
                 CALL INFO(SolverName,Message,Level=4)
              END IF
              LatentHeat = GetConstReal(Material, 'Enthalpy Latent Heat', Found)
              IF (.NOT. Found) THEN
                 WRITE(Message,'(a,I2,a,I2,a)') 'No constant >Enthalpy Latent Heat<  found in Material ', &
                      material_id
                 CALL FATAL(SolverName, Message)
              END IF

              EnthalpyRange = GetConstReal(Material, 'Enthalpy Range', Found)
              IF (.NOT. Found) THEN
                 WRITE(Message,'(a,I2,a,I2,a)') 'No constant >Enthalpy Range<  found in Material ', &
                      material_id
                 CALL FATAL(SolverName, Message)
              END IF

              PressureMeltingPoint(1:N) = ListGetReal( Material, 'Pressure Melting Point',&
                   N, Element % NodeIndexes, Found )
              IF (.NOT.Found) THEN
                 PressureMeltingPoint = 0.0D00
                 WRITE(Message,'(a,i5,a,i5,a)') 'Keyword >Pressure Melting Point< not found for element ',&
                      t, ' material ', material_id
                 CALL INFO(SolverName,Message,Level=4)
              END IF
              IceHeatConductivity(1:N) = ListGetReal( Material, 'Ice Heat Conductivity',  N, Element % NodeIndexes, Found )
              IF (.NOT.Found) THEN
                 WRITE(Message,'(a,i5,a,i5,a)') 'Keyword >Ice Heat Conductivity< not found for element ',&
                      t, ' material ', material_id
                 CALL FATAL(SolverName,Message)
              END IF
              H2ODiffusivity(1:N) = ListGetReal( Material, 'Water Diffusivity',  N, Element % NodeIndexes, Found )
              IF (.NOT.Found) THEN
                 WRITE(Message,'(a,i5,a,i5,a)') 'Keyword >Water Diffusivity< not found for element ',&
                      t, ' material ', material_id
                 CALL FATAL(SolverName,Message)
              END IF

              VelocityScaling = GetConstReal( Material, 'Enthalpy Velocity Scaling',  Found )
              IF (.NOT.Found) THEN
                 VelocityScaling = 1.0d00
              END IF

              !CALL ListGetRealArray( Material, 'Ice Heat Conductivity', Hwrk1, n, Element % NodeIndexes )
              !CALL ListGetRealArray( Material, 'Water Heat Conductivity', Hwrk2, n, Element % NodeIndexes )
 
              EnthalpyDiffusivity = 0.0d0 
              LocalEnthalpy(1:N) = Enthalpy(EnthalpyPerm(Element % NodeIndexes(1:N)))
              DO k=1,N
                 ! check whether we are in temperate ice or in cold ice or in between
                 EnthalpyH2O(k) = H2OHeatCapacity(k)*PressureMeltingPoint(k)
                 EnthalpyIce(k) = EnthalpyH2O(k) - LatentHeat

                 ! cold ice
                 IF ( LocalEnthalpy(k) <= (EnthalpyIce(k) - EnthalpyRange)) THEN
                    WaterContent = 0.0
                    DO i=1,3
                       EnthalpyDiffusivity( i,i,1:N ) = IceHeatConductivity(k)/iceHeatCapacity(k)
                    END DO
                 ! transition cold ice - temperate ice (regularization!)
                 ELSE IF (LocalEnthalpy(k) > (EnthalpyIce(k) - EnthalpyRange) &
                      .AND. LocalEnthalpy(k)  < EnthalpyIce(k)) THEN
                    WaterContent = 0.0
                    DO i=1,3
                       EnthalpyDiffusivity( i,i,1:N ) = ((EnthalpyIce(k) - LocalEnthalpy(k))/EnthalpyRange)&
                            * IceHeatConductivity(k)/IceHeatCapacity(k) + &
                            Density(k)*H2ODiffusivity(k)*(((LocalEnthalpy(k) - EnthalpyIce(k))/EnthalpyRange) + 1.0)
                    END DO
                 ! temperate ice 
                 ELSE IF (LocalEnthalpy(k) >= EnthalpyIce(k) .AND. LocalEnthalpy(k) < EnthalpyH2O(k) ) THEN
                    WaterContent = (LocalEnthalpy(k) - EnthalpyIce(k))/LatentHeat
                    DO i=1,3
                       EnthalpyDiffusivity( i,i,1:N ) = Density(k)*H2ODiffusivity(k)
                    END DO
                 ! pure water
                 ELSE IF (LocalEnthalpy(k) >= EnthalpyH2O(k)) THEN
                    WaterContent = 1.0
                    CALL FATAL(SolverName, 'Melted all ice - this case is not included')
                 END IF

                 EffectiveHeatCapacity(k)  = (1.0-WaterContent)*IceHeatCapacity(k) & 
                      + WaterContent*H2OHeatCapacity(k) 

              END DO


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
                 CALL INFO(SolverName,'Strain heating included in temperature model', &
                 Level=15)
                 Viscosity(1:N) = GetReal( Material,'Viscosity', found)
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
              MeshVelocity = &
                  MeshVelocity * VelocityScaling
              !------------------------------------------------------------------------------         
              ! asuming convection or ALE mesh contribution by default
              !------------------------------------------------------------------------------         
              DO i=1,N
                 C1(i) = Density(i)
!                 PRINT *, 'C1(',i,')=',C1(i)
              END DO
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
                          IceVeloU(i) = FlowSolution( NSDOFs*k-2 )*VelocityScaling
                          IceVeloV(i) = FlowSolution( NSDOFs*k-1 )*VelocityScaling
                          IceVeloW(i) = 0.0D0
                       CASE(4)
                          IceVeloU(i) = FlowSolution( NSDOFs*k-3 )*VelocityScaling
                          IceVeloV(i) = FlowSolution( NSDOFs*k-2 )*VelocityScaling
                          IceVeloW(i) = FlowSolution( NSDOFs*k-1 )*VelocityScaling
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
              !                                (the wrong) enthalpy and viscosity
              !------------------------------------------------------------------------------
              Unity = 1.0d00
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
                 CALL DiffuseConvectiveCompose( &
                      MASS, STIFF, FORCE, LOAD, &
                      Unity, C0, C1(1:N), EnthalpyDiffusivity, &
                      .FALSE., Zero, Zero, IceVeloU, IceVeloV, IceVeloW, &
                      MeshVelocity(1,1:N),MeshVelocity(2,1:N),MeshVelocity(3,1:N),&
                      Viscosity, Density, Pressure, Zero, Zero,&
                      .FALSE., Stabilize, Bubbles, Element, n, ElementNodes )
              ! special coords (account for metric)
              !-----------------------------------
              ELSE
                 CALL DiffuseConvectiveGenCompose( &
                      MASS, STIFF, FORCE, LOAD, &
                      IceHeatCapacity, C0, C1(1:N), EnthalpyDiffusivity, &
                      .FALSE., Zero, Zero, IceVeloU, IceVeloV, IceVeloW, &
                      MeshVelocity(1,1:N),MeshVelocity(2,1:N),MeshVelocity(3,1:N), Viscosity,&
                      Density, Pressure, Zero, Zero,.FALSE.,&
                      Stabilize, Element, n, ElementNodes )

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


           !------------------------------------------------------------------------------
           ! Neumann & Newton boundary conditions
           !------------------------------------------------------------------------------
           DO t=1, Solver % Mesh % NumberOfBoundaryElements
 
              ! get element information
              Element => GetBoundaryElement(t)
              IF ( .NOT.ActiveBoundaryElement() ) CYCLE
              n = GetElementNOFNodes()
              IF ( GetElementFamily() == 1 ) CYCLE

              IF (ParEnv % myPe .NE. Element % partIndex) CYCLE

              BC => GetBC()
              bc_id = GetBCId( Element )
              CALL GetElementNodes( ElementNodes )

              LocalEnthalpy(1:N) = Enthalpy(EnthalpyPerm(Element % NodeIndexes(1:N)))
 

              IF ( ASSOCIATED( BC ) ) THEN            
                 ! Check that we are on the correct boundary part!
                 STIFF=0.0D00
                 FORCE=0.0D00
                 MASS=0.0D00
                 LOAD=0.0D00
                 TransferCoeff = 0.0D00
                 EnthalpyExt = 0.0D00
                 FluxBC = .FALSE.
                 FluxBC =  GetLogical(BC,TRIM(Solver % Variable % Name) // ' Flux BC', Found)

                 IF (FluxBC) THEN
                    !---------------
                    !BC: -rho alpha @H/@n = q
                    !---------------
                    LocalEnthalpyFlux(1:N) = GetReal( BC, TRIM(Solver % Variable % Name) // ' Flux', Found )
                    IF (.NOT.Found) THEN
                       LocalEnthalpyFlux = 0.0d00
                       WRITE(Message,'(a,i5, a, i5)') 'No boundary value >',&
                            TRIM(Solver % Variable % Name) // ' Flux',&
                            '< found for boundary element no. ',t,'. Setting zero flux'
                       CALL WARN(SolverName,Message)
                    ELSE
                       DO k=1,N
                          ! check whether we are in temperate ice or in cold ice or in between
                          EnthalpyH2O(k) = H2OHeatCapacity(k)*PressureMeltingPoint(k)
                          EnthalpyIce(k) = EnthalpyH2O(k) - LatentHeat
                          IF ( LocalEnthalpy(k) < EnthalpyIce(k)) THEN
                             LOAD(k)  = LOAD(k) + LocalEnthalpyFlux(k)
                          END IF                         
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


           CALL Info( SolverName, 'Assembly done', Level=4 )

           !------------------------------------------------------------------------------
           !     Solve the system and check for convergence
           !------------------------------------------------------------------------------
           at = CPUTime() - at
           st = CPUTime()

           PrevNorm = Solver % Variable % Norm

           Norm = DefaultSolve()

           st = CPUTime()-st
           totat = totat + at
           totst = totst + st
           WRITE(Message,'(a,i4,a,F8.2,F8.2)') 'iter: ',iter,' Assembly: (s)', at, totat
           CALL Info( SolverName, Message, Level=4 )
           WRITE(Message,'(a,i4,a,F8.2,F8.2)') 'iter: ',iter,' Solve:    (s)', st, totst
           CALL Info( SolverName, Message, Level=4 )


           IF ( PrevNorm + Norm /= 0.0d0 ) THEN
              RelativeChange = 2.0d0 * ABS( PrevNorm-Norm ) / (PrevNorm + Norm)
           ELSE
              RelativeChange = 0.0d0
           END IF

           WRITE( Message, * ) 'Result Norm   : ',Norm
           CALL Info( SolverName, Message, Level=4 )
           WRITE( Message, * ) 'Relative Change : ',RelativeChange
           CALL Info( SolverName, Message, Level=4 )

           SystemMatrix % Values = OldValues
           ForceVector = OldRHS

           !----------------------
           ! check for convergence
           !----------------------
           IF ( RelativeChange < NonlinearTol ) THEN
              EXIT
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

!------------------------------------------------------------------------------
   END SUBROUTINE EnthalpySolver
!------------------------------------------------------------------------------

FUNCTION EnthalpyNeumannBoundary (Model, Node, Enthalpy ) RESULT(SurfaceEnthalpyFlux)
  USE Types
  USE CoordinateSystems
  USE SolverUtils
  USE ElementDescription
  USE DefUtils
!-----------------------------------------------------------
  IMPLICIT NONE
!------------ external variables ------------------------------------------
  TYPE(Model_t) :: Model
  INTEGER :: Node
  REAL (KIND=dp) :: Enthalpy
  REAL (KIND=dp) ::  SurfaceEnthalpyFlux
!------------ internal variables-------------------------------------------
  REAL (KIND=dp), ALLOCATABLE :: &
       PressureMeltingPoint(:), LatentHeat(:),H2OHeatCapacity(:), GeothermalHeatFlux(:)
  REAL(KIND=dp) :: EnthalpyH2O, EnthalpyIce
  INTEGER :: i,k,l,t,nMax,body_id,other_body_id,material_id,&
       NBoundary,NParent,istat,ParentElementNode,BoundaryElementNode
  LOGICAL :: FirstTime=.TRUE., GotIt
  CHARACTER(LEN=MAX_NAME_LEN) :: FunctionName
  TYPE(ValueList_t), POINTER :: ParentMaterial, BC
  TYPE(Element_t), POINTER :: BoundaryElement, ParentElement 

  SAVE FirstTime, H2OHeatCapacity,&
       PressureMeltingPoint, LatentHeat,GeothermalHeatFlux

  !-------------------------------------------
  ! Allocations 
  !------------------------------------------- 
  IF (FirstTime) THEN
     WRITE(FunctionName, '(A)') 'EnthalpySolver(EnthalpyNeumannBoundary)'
     nMax = Model % MaxElementNodes
     ALLOCATE(PressureMeltingPoint(nMax), H2OHeatCapacity(nMax), GeothermalHeatFlux(nMax),&
          LatentHeat(nMax), &
         STAT=istat)
     IF ( istat /= 0 ) THEN
        CALL FATAL(FunctionName,'Memory allocation error, Aborting.')
     END IF
     FirstTime = .FALSE.
     CALL INFO(FunctionName,'Memory allocation done', level=3)
  END IF
  !-------------------------------------------
  ! get boundary element and BC pointer
  !-------------------------------------------   
  IF ( .NOT. ASSOCIATED(Model % CurrentElement) ) THEN
     CALL FATAL(FunctionName, 'Model % CurrentElement not associated')
  ELSE
     BoundaryElement => Model % CurrentElement
  END IF
  IF ( .NOT. ASSOCIATED(BoundaryElement) ) THEN
     CALL FATAL(FunctionName,'No boundary element found')
  END IF
  BC => GetBC()
  !---------------------------------------------------------
  ! Do nothing if this is a halo-element of a parallel mesh
  !---------------------------------------------------------
  IF (ParEnv % myPe .NE. BoundaryElement % partIndex) RETURN
  !------------------------------------------------
  ! get corresponding point within boundary element
  !------------------------------------------------
  NBoundary = BoundaryElement % Type % NumberOfNodes
    DO BoundaryElementNode=1,NBoundary
     IF (Node .EQ. BoundaryElement % NodeIndexes(BoundaryElementNode)) THEN
        GotIt = .TRUE.
        EXIT
     END IF
  END DO
  IF (.NOT.GotIt) THEN
     CALL FATAL(FunctionName,'Node not found in Current Boundary Element')
  END IF
  !------------------------------------------------
  ! get parent element of boundary element
  !------------------------------------------------
  other_body_id = BoundaryElement % BoundaryInfo % outbody
  IF (other_body_id < 1) THEN ! only one body in calculation
     ParentElement => BoundaryElement % BoundaryInfo % Right
     IF ( .NOT. ASSOCIATED(ParentElement) ) ParentElement => BoundaryElement % BoundaryInfo % Left
  ELSE ! we are dealing with a body-body boundary and asume that the normal is pointing outwards
     ParentElement => BoundaryElement % BoundaryInfo % Right
     IF (ParentElement % BodyId == other_body_id) ParentElement => BoundaryElement % BoundaryInfo % Left
  END IF
  ! just to be on the save side, check again
  IF ( .NOT. ASSOCIATED(ParentElement) ) THEN
     WRITE(Message,'(A,I10,A)')&
          'Parent Element for Boundary element no. ',&
          BoundaryElement % ElementIndex, ' not found'
     CALL FATAL(FunctionName,Message)
  END IF  
  body_id = ParentElement % BodyId
  NParent = ParentElement % Type % NumberOfNodes
  DO ParentElementNode=1,NParent
     IF ( Node == ParentElement % NodeIndexes(ParentElementNode) ) EXIT
  END DO
  !-------------------------------------------
  ! get material properties of parent element
  !-------------------------------------------
  material_id = ListGetInteger(Model % Bodies(body_id) % Values, 'Material', GotIt)
  ParentMaterial => Model % Materials(material_id) % Values
  IF ((.NOT. ASSOCIATED(ParentMaterial)) .OR. (.NOT. GotIt)) THEN
     WRITE(Message,'(A,I10,A,I10)')&
          'No material values found for body no ', body_id,&
          ' under material id ', material_id
     CALL FATAL(FunctionName,Message)
  END IF 
  Model % CurrentElement => ParentElement ! to be on the save side, if this is needed by a function

  H2OHeatCapacity(1:NParent) = ListGetReal( ParentMaterial,'Water Heat Capacity', NParent, &
       ParentElement % NodeIndexes, GotIt )
  IF (.NOT. GotIt) THEN
     WRITE(Message,'(a,I2,a,I2,a)') 'No >Water Heat Capacity<  found in Material ', &
          material_id,' for node ', Node
     CALL FATAL(FunctionName, Message)
  END IF
  PressureMeltingPoint(1:NParent) = ListGetReal( ParentMaterial, 'Pressure Melting Point',&
        NParent, ParentElement % NodeIndexes, GotIt )
  IF (.NOT.GotIt) THEN
     WRITE(Message,'(a,i5,a,i5,a)') 'Keyword >Pressure Melting Point< not found for element ',&
          t, ' Material ', material_id
     CALL FATAL(FunctionName,Message)
  END IF
  LatentHeat(1:NParent) = ListGetReal(ParentMaterial, 'Enthalpy Latent Heat', &
       NParent, ParentElement % NodeIndexes, GotIt)
  IF (.NOT. GotIt) THEN
     CALL FATAL(FunctionName,'No value for >Enthalpy Latent Heat< found')
  END IF
  Model % CurrentElement => BoundaryElement
  !----------------------------------
  ! get the boundary condition values
  !----------------------------------
  IF (.NOT.ASSOCIATED(BC)) THEN
     CALL FATAL(FunctionName,'No Boundary Condition associated')
  END IF
  GeothermalHeatFlux(1:NBoundary) = GetReal(BC, 'Geothermal Heat Flux', GotIt)
    IF (.NOT.GotIt) THEN
     WRITE(Message,'(a,i5,a,i5,a)') 'Keyword >Geothermal Heat Flux< not found for node ',&
          Node
     CALL FATAL(FunctionName,Message)
  END IF

  !---------------------
  ! finally, do the math
  !---------------------
  EnthalpyH2O = H2OHeatCapacity(ParentElementNode) * PressureMeltingPoint(ParentElementNode)
  EnthalpyIce = EnthalpyH2O - LatentHeat(ParentElementNode)
  ! cold ice
  IF (Enthalpy < EnthalpyIce) THEN ! cold ice boundary
     SurfaceEnthalpyFlux =  GeothermalHeatFlux(BoundaryElementNode)
  ELSE IF ((Enthalpy .GE. EnthalpyIce) .AND. (EnthalpyIce < EnthalpyH2O)) THEN !temperate ice boundary
     SurfaceEnthalpyFlux = 0.0D00
  ELSE
     WRITE(Message,'(a,i5,a)') 'Node ',&
          Node, ' is already pure water. Aborting'
     CALL FATAL(FunctionName,Message)
  END IF

END FUNCTION EnthalpyNeumannBoundary

FUNCTION EnthalpyDirichletBoundary (Model, Node, Enthalpy ) RESULT(SurfaceEnthalpy)
  USE Types
  USE CoordinateSystems
  USE SolverUtils
  USE ElementDescription
  USE DefUtils
!-----------------------------------------------------------
  IMPLICIT NONE
!------------ external variables ------------------------------------------
  TYPE(Model_t) :: Model
  INTEGER :: Node
  REAL (KIND=dp) :: Enthalpy
  REAL (KIND=dp) ::  SurfaceEnthalpy
!------------ internal variables-------------------------------------------
  REAL (KIND=dp), ALLOCATABLE :: IceHeatConductivity(:), IceHeatCapacity(:),&
       PressureMeltingPoint(:), H2OHeatCapacity(:),LatentHeat(:),&
       SurfaceTemperature(:), SurfaceWaterContent(:)
  REAL(KIND=dp) :: EnthalpyH2O, EnthalpyIce
  INTEGER :: i,k,l,t,nMax,body_id,other_body_id,material_id,&
       NBoundary,NParent,istat,ParentElementNode,BoundaryElementNode
  LOGICAL :: FirstTime=.TRUE., GotIt
  CHARACTER(LEN=MAX_NAME_LEN) :: FunctionName
  TYPE(ValueList_t), POINTER :: ParentMaterial, BC
  TYPE(Element_t), POINTER :: BoundaryElement, ParentElement 

  SAVE FirstTime, IceHeatConductivity, H2OHeatCapacity,&
       IceHeatCapacity, PressureMeltingPoint, LatentHeat,&
       SurfaceTemperature, SurfaceWaterContent

  !-------------------------------------------
  ! Allocations 
  !------------------------------------------- 
  IF (FirstTime) THEN
     WRITE(FunctionName, '(A)') 'EnthalpySolver(EnthalpyDirichletBoundary)'
     nMax = Model % MaxElementNodes
     ALLOCATE(IceHeatConductivity(nMax), IceHeatCapacity(nMax),&
          PressureMeltingPoint(nMax), H2OHeatCapacity(nMax), &
          LatentHeat(nMax), SurfaceTemperature(nMax), SurfaceWaterContent(nMax), &
         STAT=istat)
     IF ( istat /= 0 ) THEN
        CALL FATAL(FunctionName,'Memory allocation error, Aborting.')
     END IF
     FirstTime = .FALSE.
     CALL INFO(FunctionName,'Memory allocation done', level=3)
  END IF
  !-------------------------------------------
  ! get boundary element and BC pointer
  !-------------------------------------------   
  IF ( .NOT. ASSOCIATED(Model % CurrentElement) ) THEN
     CALL FATAL(FunctionName, 'Model % CurrentElement not associated')
  ELSE
     BoundaryElement => Model % CurrentElement
  END IF
  IF ( .NOT. ASSOCIATED(BoundaryElement) ) THEN
     CALL FATAL(FunctionName,'No boundary element found')
  END IF
  BC => GetBC()
  !---------------------------------------------------------
  ! Do nothing if this is a halo-element of a parallel mesh
  !---------------------------------------------------------
  IF (ParEnv % myPe .NE. BoundaryElement % partIndex) RETURN
  !------------------------------------------------
  ! get corresponding point within boundary element
  !------------------------------------------------
  NBoundary = BoundaryElement % Type % NumberOfNodes
    DO BoundaryElementNode=1,NBoundary
     IF (Node .EQ. BoundaryElement % NodeIndexes(BoundaryElementNode)) THEN
        GotIt = .TRUE.
        EXIT
     END IF
  END DO
  IF (.NOT.GotIt) THEN
     CALL FATAL(FunctionName,'Node not found in Current Boundary Element')
  END IF
  !------------------------------------------------
  ! get parent element of boundary element
  !------------------------------------------------
  other_body_id = BoundaryElement % BoundaryInfo % outbody
  IF (other_body_id < 1) THEN ! only one body in calculation
     ParentElement => BoundaryElement % BoundaryInfo % Right
     IF ( .NOT. ASSOCIATED(ParentElement) ) ParentElement => BoundaryElement % BoundaryInfo % Left
  ELSE ! we are dealing with a body-body boundary and asume that the normal is pointing outwards
     ParentElement => BoundaryElement % BoundaryInfo % Right
     IF (ParentElement % BodyId == other_body_id) ParentElement => BoundaryElement % BoundaryInfo % Left
  END IF
  ! just to be on the save side, check again
  IF ( .NOT. ASSOCIATED(ParentElement) ) THEN
     WRITE(Message,'(A,I10,A)')&
          'Parent Element for Boundary element no. ',&
          BoundaryElement % ElementIndex, ' not found'
     CALL FATAL(FunctionName,Message)
  END IF  
  body_id = ParentElement % BodyId
  NParent = ParentElement % Type % NumberOfNodes
  DO ParentElementNode=1,NParent
     IF ( Node == ParentElement % NodeIndexes(ParentElementNode) ) EXIT
  END DO
  !-------------------------------------------
  ! get material properties of parent element
  !-------------------------------------------
  material_id = ListGetInteger(Model % Bodies(body_id) % Values, 'Material', GotIt)
  ParentMaterial => Model % Materials(material_id) % Values
  IF ((.NOT. ASSOCIATED(ParentMaterial)) .OR. (.NOT. GotIt)) THEN
     WRITE(Message,'(A,I10,A,I10)')&
          'No material values found for body no ', body_id,&
          ' under material id ', material_id
     CALL FATAL(FunctionName,Message)
  END IF 
  Model % CurrentElement => ParentElement ! to be on the save side, if this is needed by a function
  IceHeatConductivity(1:NParent) = ListGetReal( ParentMaterial,'Ice Heat Conductivity ', NParent, &
       ParentElement % NodeIndexes, GotIt )
  IF (.NOT. GotIt) THEN
     WRITE(Message,'(a,I2,a,I2,a)') 'No >Ice Heat Conductivity< found in Material ', &
          material_id,' for node ', Node
     CALL FATAL(FunctionName, Message)
  END IF
  IceHeatCapacity(1:NParent) = ListGetReal( ParentMaterial,'Ice Heat Capacity', NParent, &
       ParentElement % NodeIndexes, GotIt )
  IF (.NOT. GotIt) THEN
     WRITE(Message,'(a,I2,a,I2,a)') 'No >Ice Heat capacity found in Material< ', &
          material_id,' for node ', Node
     CALL FATAL(FunctionName, Message)
  END IF
  H2OHeatCapacity(1:NParent) = ListGetReal( ParentMaterial,'Water Heat Capacity', NParent, &
       ParentElement % NodeIndexes, GotIt )
  IF (.NOT. GotIt) THEN
     WRITE(Message,'(a,I2,a,I2,a)') 'No >Water Heat Capacity<  found in Material ', &
          material_id,' for node ', Node
     CALL FATAL(FunctionName, Message)
  END IF
  PressureMeltingPoint(1:NParent) = ListGetReal( ParentMaterial, 'Pressure Melting Point',&
        NParent, ParentElement % NodeIndexes, GotIt )
  IF (.NOT.GotIt) THEN
     WRITE(Message,'(a,i5,a,i5,a)') 'Keyword >Pressure Melting Point< not found for element ',&
          t, ' Material ', material_id
     CALL FATAL(FunctionName,Message)
  END IF
  LatentHeat(1:NParent) = ListGetReal(ParentMaterial, 'Enthalpy Latent Heat', &
       NParent, ParentElement % NodeIndexes, GotIt)
  IF (.NOT. GotIt) THEN
     CALL FATAL(FunctionName,'No value for >Enthalpy Latent Heat< found')
  END IF
  Model % CurrentElement => BoundaryElement
  !----------------------------------
  ! get the boundary condition values
  !----------------------------------
  IF (.NOT.ASSOCIATED(BC)) THEN
     CALL FATAL(FunctionName,'No Boundary Condition associated')
  END IF
  SurfaceTemperature(1:NBoundary) = GetReal(BC, 'Surface Temperature', GotIt)
    IF (.NOT.GotIt) THEN
     WRITE(Message,'(a,i5,a,i5,a)') 'Keyword >Surface Temperature< not found for node ',&
          Node
     CALL FATAL(FunctionName,Message)
  END IF
  SurfaceWaterContent(1:NBoundary) = GetReal(BC, 'Surface Water Content', GotIt)
    IF (.NOT.GotIt) THEN
     WRITE(Message,'(a,i5)') 'Keyword >Surface Water Content< not found for node ',&
          Node
     CALL FATAL(FunctionName,Message)
  END IF
  !---------------------
  ! finally, do the math
  !---------------------
  EnthalpyH2O = H2OHeatCapacity(ParentElementNode) * PressureMeltingPoint(ParentElementNode)
  EnthalpyIce = EnthalpyH2O - LatentHeat(ParentElementNode)
  ! cold ice
  IF (Enthalpy < EnthalpyIce) THEN ! cold ice boundary
     SurfaceEnthalpy = EnthalpyIce +&
          IceHeatCapacity(ParentElementNode) * (SurfaceTemperature(BoundaryElementNode) - PressureMeltingPoint(ParentElementNode))
  ELSE IF ((Enthalpy .GE. EnthalpyIce) .AND. (EnthalpyIce < EnthalpyH2O)) THEN !temperate ice boundary
     SurfaceEnthalpy = EnthalpyIce +&
          LatentHeat(ParentElementNode) * SurfaceWaterContent(BoundaryElementNode)
  ELSE
     WRITE(Message,'(a,i5,a)') 'Node ',&
          Node, ' is already pure water. Aborting'
     CALL FATAL(FunctionName,Message)
  END IF
END FUNCTION EnthalpyDirichletBoundary


!----------------------------------------------------------------------------------    
RECURSIVE SUBROUTINE TemperatureEnthalpy( Model,Solver,dt,TransientSimulation )
!---------------------------------------------------------------------------------- 

  USE DefUtils

  IMPLICIT NONE
!------------------------------------------------------------------------------
!******************************************************************************
!
!  Solve stress equations for one timestep
!
!  ARGUMENTS:
!
!  TYPE(Model_t) :: Model,  
!     INPUT: All model information (mesh,materials,BCs,etc...)
!
!  TYPE(Solver_t) :: Solver
!     INPUT: Linear equation solver options
!
!  REAL(KIND=dp) :: dt,
!     INPUT: Timestep size for time dependent simulations (NOTE: Not used
!            currently)
!
!******************************************************************************

   TYPE(Model_t)  :: Model
   TYPE(Solver_t), TARGET :: Solver

   LOGICAL ::  TransientSimulation
   REAL(KIND=dp) :: dt
!------------------------------------------------------------------------------
!    Local variables
!------------------------------------------------------------------------------

 INTEGER :: N, i, j, k, t, TempDOFs, HomoTempDOFs,  EnthalpyDOFs, DIM, elementNbNodes
 TYPE(Element_t),POINTER :: BoundaryElement, CurrentElement
 REAL(KIND=dp), POINTER :: TempValues(:), HomoTempValues(:), EnthalpyValues(:),&
      WaterContValues(:), EnthalpyStateValues(:)
 TYPE(Variable_t), POINTER :: TempSol, HomoTempSol, EnthalpySol, WaterContSol, EnthalpyStateSol
 TYPE(ValueList_t),POINTER :: Material, BodyForce
 INTEGER, POINTER :: TempPerm(:), HomoTempPerm(:),  EnthalpyPerm(:),&
      WaterContPerm(:),EnthalpyStatePerm(:),NodeIndexes(:)
 REAL(KIND=dp) :: enthalpyPureIce, enthalpyPureWater, latentHeat, enthalpy
 REAL(KIND=dp) :: enthalpyRange, mixtureCapacity, waterContent
 REAL (KIND=dp), ALLOCATABLE :: waterCapacity(:), T_m(:), iceCapacity(:)
 INTEGER :: material_id, body_id, bf_id, istat
 LOGICAL :: GotIt, FirstTime = .TRUE.
 CHARACTER(LEN=MAX_NAME_LEN) :: SolverName
!------------ remember this -------------------------------------------------
  Save FirstTime,waterCapacity, iceCapacity,T_m, SolverName
  
  !-------------------------------------------
  ! Allocations 
  !------------------------------------------- 
  IF (FirstTime) THEN
     WRITE(SolverName, '(A)') 'EnthalpySolver(TemperatureEnthalpy)'
     N = Solver % Mesh % MaxElementNodes
     ALLOCATE( waterCapacity(N), T_m(N), iceCapacity(N), &
          STAT=istat)
     IF ( istat /= 0 ) THEN
        CALL FATAL(SolverName,'Memory allocation error, Aborting.')
     END IF
     FirstTime = .FALSE.
     CALL INFO(SolverName,'Memory allocation done', level=3)
  END IF

  TempSol => VariableGet( Solver % Mesh % Variables, 'Enthalpy Temperature' )
  IF ( ASSOCIATED( TempSol ) ) THEN
       TempPerm    => TempSol % Perm
       TempValues  => TempSol % Values
       TempDOFs = TempSol % DOFs
       CALL Info(SolverName, 'Found variable >Enthalpy Temperature<  variable', Level=4)
  ELSE
       CALL Info(SolverName, 'Could not find >Enthalpy Temperature<  variable', Level=4)
  END IF

  HomoTempSol => VariableGet( Solver % Mesh % Variables, 'Enthalpy Temp Homologous' )
  IF ( ASSOCIATED( HomoTempSol ) ) THEN
       HomoTempPerm    => HomoTempSol % Perm
       HomoTempValues  => HomoTempSol % Values
       HomoTempDOFs = HomoTempSol % DOFs
       CALL Info(SolverName, 'Found variable >Enthalpy Temp Homologous<  variable', Level=4)
  ELSE
       CALL Info(SolverName, 'Could not find >Enthalpy Temp Homologous< variable. No output of this variable', Level=4)
  END IF

  EnthalpySol => VariableGet( Model % Solver % Mesh % Variables, "Enthalpy" )
  IF ( ASSOCIATED( EnthalpySol) ) THEN
        EnthalpyPerm    => EnthalpySol % Perm
        EnthalpyValues  => EnthalpySol % Values
        EnthalpyDOFs = EnthalpySol % DOFs   
        CALL Info(SolverName, 'Found variable >Enthalpy<  variable', Level=4)
  ELSE
        CALL Fatal(SolverName, 'No variable for Enthalpy associated.')
  END IF

  WaterContSol => VariableGet( Model % Solver % Mesh % Variables, "Water Content" )
  IF ( ASSOCIATED( WaterContSol) ) THEN
        WaterContPerm    => WaterContSol % Perm
        WaterContValues  => WaterContSol % Values
        CALL Info(SolverName, 'Found variable >Water Content<  variable', Level=4)
  ELSE
        CALL Info(SolverName, 'No variable for >Water Content< associated. Not output',level=3)
  END IF

  EnthalpyStateSol => VariableGet( Model % Solver % Mesh % Variables, "Enthalpy State" )
  IF ( ASSOCIATED( EnthalpyStateSol) ) THEN
        EnthalpyStatePerm    => EnthalpyStateSol % Perm
        EnthalpyStateValues  => EnthalpyStateSol % Values 
        CALL Info(SolverName, 'Found variable >Enthalpy State<  variable', Level=4)
  ELSE
     CALL Info(SolverName, 'No variable for >Enthalpy State< associated. Not output',level=3)
  END IF
  DIM = CoordinateSystemDimension()

  !---------------------------------------------------------------
  ! Compute the temperatures from the entropy for all nodes
  !---------------------------------------------------------------
   DO t=1,Solver % NumberOFActiveElements
       CurrentElement => GetActiveElement(t)
       IF (ParEnv % myPe .NE. CurrentElement % partIndex) CYCLE
       elementNbNodes = GetElementNOFNodes( CurrentElement)
       NodeIndexes => CurrentElement % NodeIndexes
       body_id = Model % CurrentElement % BodyId
       material_id = ListGetInteger(Model % Bodies(body_id) % Values, 'Material', GotIt)
       Material => Model % Materials(material_id) % Values
       IF ((.NOT. ASSOCIATED(Material)) .OR. (.NOT. GotIt)) THEN
              WRITE(Message,'(A)') 'No material found for model '
              CALL FATAL(SolverName,Message)
       END IF
       latentHeat = GetConstReal(Material, 'Enthalpy Latent Heat', GotIt)
       IF(.NOT. GotIt) THEN
           CALL FATAL(SolverName, 'Enthalpy Latent Heat not found in Material')
       END IF

       enthalpyRange = GetConstReal(Material, 'Enthalpy Range', GotIt)
       IF(.NOT. GotIt) THEN
           CALL FATAL(SolverName, 'Enthalpy Range not found in Material')
       END IF


       T_m(1:elementNbNodes) = ListGetReal( Material,'Pressure Melting Point', elementNbNodes, &
       Model % CurrentElement % NodeIndexes, GotIt )
       IF (.NOT. GotIt) THEN
          WRITE(Message,'(a,I2,a,I2,a)') 'No Temp Upper Limit  found in Material ', &
          material_id
          CALL FATAL(SolverName, Message)
       END IF

       iceCapacity(1:elementNbNodes) = ListGetReal( Material,'Ice Heat Capacity', elementNbNodes, &
       Model % CurrentElement % NodeIndexes, GotIt )
       IF (.NOT. GotIt) THEN
          WRITE(Message,'(a,I2,a,I2,a)') 'No Ice Heat Capacity specific heat  found in Material ', &
          material_id
          CALL FATAL(SolverName, Message)
       END IF

       Watercapacity(1:elementNbNodes) = ListGetReal( Material,'Water Heat Capacity', elementNbNodes, &
       Model % CurrentElement % NodeIndexes, GotIt )
       IF (.NOT. GotIt) THEN
          WRITE(Message,'(a,I2,a,I2,a)') 'No >Water Heat Capcity< found in Material ', &
          material_id
          CALL FATAL(SolverName, Message)
       END IF



       DO i=1,elementNbNodes

          k = Model % CurrentElement % NodeIndexes(i)

          !---------------------------------------------------------------
          ! compute the enthalpy of puce ice and pure water, H_s, and H_l
          !---------------------------------------------------------------

          enthalpyPureWater = waterCapacity(i)*T_m(i)

          enthalpyPureIce = -latentHeat+enthalpyPureWater

          

          !---------------------------------------------------------------
          ! enthalpy
          !---------------------------------------------------------------

          enthalpy = EnthalpyValues(EnthalpyPerm(k))

          !---------------------------------------------------------------
          ! compute the water content
          !---------------------------------------------------------------

          IF ( enthalpy <= (enthalpyPureIce-enthalpyRange)) THEN
             waterContent = 0.0
          ELSE IF (enthalpy > (enthalpyPureIce-enthalpyRange) .AND. enthalpy < enthalpyPureIce) THEN
             waterContent = 0.0
          ELSE IF (enthalpy >= enthalpyPureIce .AND. enthalpy < enthalpyPureWater ) THEN
             waterContent = (enthalpy-enthalpyPureIce)/latentHeat
          ELSE IF (enthalpy >= enthalpyPureWater) THEN
             waterContent = 1.0
          END IF

          IF (ASSOCIATED(WaterContValues)) &
               WaterContValues(WaterContPerm(k)) = waterContent

          !--------------------------------------------
          ! Compute the capacity of the mixture
          !--------------------------------------------

          mixtureCapacity  = (1.0-waterContent)*iceCapacity(i) + waterContent*waterCapacity(i)

 

          IF (enthalpy <= (enthalpyPureIce-enthalpyRange)) THEN    
             TempValues(TempPerm(k)) = (1.0/mixtureCapacity)*(enthalpy-enthalpyPureIce) + T_m(i)
             IF (ASSOCIATED(EnthalpyStateValues))  EnthalpyStateValues(EnthalpyStatePerm(k)) = -1.0
          ELSE IF (enthalpy > (enthalpyPureIce-enthalpyRange) .AND. enthalpy < enthalpyPureIce) THEN
             TempValues(TempPerm(k)) = T_m(i)
             IF (ASSOCIATED(EnthalpyStateValues))  EnthalpyStateValues(EnthalpyStatePerm(k)) = 0.0
          ELSE IF (enthalpy >= enthalpyPureIce .AND. enthalpy < enthalpyPureWater ) THEN
             TempValues(TempPerm(k)) = T_m(i)
             IF (ASSOCIATED(EnthalpyStateValues))  EnthalpyStateValues(EnthalpyStatePerm(k)) = 1.0
          ELSE IF (enthalpy >= enthalpyPureWater) THEN
             TempValues(TempPerm(k)) = (1.0/mixtureCapacity)*(enthalpy+enthalpyPureWater) + T_m(i)
             IF (ASSOCIATED(EnthalpyStateValues))  EnthalpyStateValues(EnthalpyStatePerm(k)) = 10.0
          END IF
       
          IF (ASSOCIATED(HomoTempValues)) &
               HomoTempValues(HomoTempPerm(k)) = TempValues(TempPerm(k)) - T_m(i)

       END DO ! end loop nodes of element
 
    END DO ! end loop elements

  WRITE(Message,*) 'solve done', minval(TempValues), maxval(TempValues)
      CALL Info( SolverName, Message, Level=4 )

  WRITE(Message,'(a)') '------------------------------------------------------'
  CALL Info(SolverName, Message, Level=3)
  CALL Info(SolverName, 'Compute Temperature:..........done', Level=3)
  WRITE(Message,'(a)') '------------------------------------------------------'
  CALL Info(SolverNAme, Message, Level=3)

!------------------------------------------------------------------------------
END SUBROUTINE TemperatureEnthalpy
!------------------------------------------------------------------------------
