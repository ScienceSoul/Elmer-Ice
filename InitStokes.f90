!----------------------------------------------------------------------------------    
RECURSIVE SUBROUTINE InitStokes( Model,Solver,dt,TransientSimulation )
!---------------------------------------------------------------------------------- 

  USE DefUtils

  IMPLICIT NONE
!------------------------------------------------------------------------------
!******************************************************************************
!
!  Initialize Stokes fields with SIA velocities and pressure
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
 
 INTEGER :: nvalue, i, j, k, l, m, n, o, FlowDOFs, SicoVxDOFS, SicoVyDOFS, SicoVzDOFS, SicoPressureDOFS, &
            DIM, istat
 TYPE(Element_t), POINTER :: CurrentElement
 TYPE(ValueList_t),POINTER :: SolverParams
 TYPE(ValueList_t), POINTER :: Material
 REAL(KIND=dp), POINTER :: FlowValues(:), SicoVxVal(:),  SicoVyVal(:),  SicoVzVal(:), SicoPressureVal(:)
 TYPE(Variable_t), POINTER :: FlowSol, SicoVxSol, SicoVySol, SicoVzSol, SicoPressureSol
 INTEGER, POINTER :: FlowPerm(:), SicoVxPerm(:), SicoVyPerm(:), SicoVzPerm(:), SicoPressurePerm(:), NodeIndexes(:)
 CHARACTER(LEN=MAX_NAME_LEN) :: SolverName, VarName
 LOGICAL :: AllocationsDone = .FALSE., Found, done = .FALSE.

 SAVE done

 SolverName = 'InitStokes'

 SolverParams => Solver % Values

 M = Model % MaxElementNodes
 dim = CoordinateSystemDimension()  

 VarName = GetString(SolverParams,'Flow Solver Name',Found )
 IF(.NOT. Found) THEN
    CALL Fatal(SolverName, 'Variable for flow solution not found')
 END IF

 FlowSol => VariableGet( Solver % Mesh % Variables, VarName )
 IF ( ASSOCIATED( FlowSol ) ) THEN
       FlowPerm    => FlowSol % Perm
       FlowValues  => FlowSol % Values
       FlowDOFs = FlowSol % DOFs
 ELSE
       CALL Fatal(SolverName, 'Could not find velocity field variable')
 END IF

 SicoVxSol => VariableGet( Solver % Mesh % Variables, 'SIAFlow 1' )
 IF ( ASSOCIATED( SicoVxSol ) ) THEN
       SicoVxPerm    => SicoVxSol % Perm
       SicoVxVal  => SicoVxSol % Values
       SicoVxDOFs = SicoVxSol % DOFs
 ELSE
       CALL Fatal(SolverName, 'Could not find variable for SIAFlow 1')
 END IF

 SicoVySol => VariableGet( Solver % Mesh % Variables, 'SIAFlow 2' )
 IF ( ASSOCIATED( SicoVySol ) ) THEN
       SicoVyPerm    => SicoVySol % Perm
       SicoVyVal  => SicoVySol % Values
       SicoVyDOFs = SicoVySol % DOFs
 ELSE
       CALL Fatal(SolverName, 'Could not find variable for SIAFlow 2')
 END IF

 SicoVzSol => VariableGet( Solver % Mesh % Variables, 'SIAFlow 3' )
 IF ( ASSOCIATED( SicoVzSol ) ) THEN
       SicoVzPerm    => SicoVzSol % Perm
       SicoVzVal  => SicoVzSol % Values
       SicoVzDOFs = SicoVzSol % DOFs
 ELSE
       CALL Fatal(SolverName, 'Could not find variable for SIAFlow 3')
 END IF

 SicoPressureSol => VariableGet( Solver % Mesh % Variables, 'SIAFlow 4' )
 IF ( ASSOCIATED( SicoPressureSol ) ) THEN
       SicoPressurePerm    => SicoPressureSol % Perm
       SicoPressureVal  => SicoPressureSol % Values
       SicoPressureDOFs = SicoPressureSol % DOFs
 ELSE
       CALL Fatal(SolverName, 'Could not find variable for SIAFlow 4')
 END IF

 IF (.NOT. done) THEN
    DO i=1,Solver % NumberOFActiveElements

      CurrentElement => GetActiveElement(i)
      NodeIndexes => CurrentElement % NodeIndexes

      DO k=1, GetElementNOFNodes(CurrentElement)
   
        j = FlowPerm(NodeIndexes(k))
        l = SicoVxPerm(NodeIndexes(k))
        m = SicoVyPerm(NodeIndexes(k))
        n = SicoVzPerm(NodeIndexes(k))
        o = SicoPressurePerm(NodeIndexes(k))

        FlowValues(FlowDOFs*(j-1)+1) = 0.0_dp !SicoVxVal(SicoVxDOFs*(l-1)+1)
        FlowValues(FlowDOFs*(j-1)+2) = 0.0_dp !SicoVyVal(SicoVyDOFs*(m-1)+1)
        FlowValues(FlowDOFs*(j-1)+3) = 0.0_dp !SicoVzVal(SicoVzDOFs*(n-1)+1)
        FlowValues(FlowDOFs*(j-1)+4) = 0.0_dp !SicoVzVal(SicoVzDOFs*(n-1)+1)

      END DO
    END DO
    WRITE(Message,'(a)') '---------------------------------------------------------------'
    CALL Info(SolverName,Message, Level=3)
    CALL Info(SolverName,'Initialize Stokes fields:..........done', Level=3)
    WRITE(Message,'(a)') '---------------------------------------------------------------'
    CALL Info(SolverName,Message, Level=3)

    done = .TRUE.
 ELSE
    CALL Info(SolverName, 'No ops....')
 END IF

!------------------------------------------------------------------------------
END SUBROUTINE InitStokes
!------------------------------------------------------------------------------
