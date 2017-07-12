!----------------------------------------------------------------------------------    
RECURSIVE SUBROUTINE createStore( Model,Solver,dt,TransientSimulation )
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
 
 INTEGER :: nvalue, i, j, k, l, m, n, o, DIM, istat, StoreDOFs, &
            TargetDOFs
 TYPE(Element_t), POINTER :: CurrentElement
 TYPE(ValueList_t),POINTER :: SolverParams
 TYPE(ValueList_t), POINTER :: Material
 REAL(KIND=dp), POINTER :: TargetValues(:), StoreValues(:)
 TYPE(Variable_t), POINTER :: TargetSol, StoreSol
 INTEGER, POINTER :: TargetPerm(:), StorePerm(:), NodeIndexes(:)
 CHARACTER(LEN=MAX_NAME_LEN) :: SolverName, VarStore, varTarget
 LOGICAL :: AllocationsDone = .FALSE., Found, done = .FALSE.

 SAVE done

 SolverName = 'createStore'

 SolverParams => Solver % Values

 M = Model % MaxElementNodes
 dim = CoordinateSystemDimension()  

 VarStore = GetString(SolverParams,'Name of variable to store',Found )
 IF(.NOT. Found) THEN
    CALL Fatal(SolverName, 'Name of variable to store not found')
 END IF

 VarTarget = GetString(SolverParams,'Name of target variable',Found )
 IF(.NOT. Found) THEN
    CALL Fatal(SolverName, 'Name of target variable not found')
 END IF

 StoreSol => VariableGet( Solver % Mesh % Variables, VarStore)
 IF ( ASSOCIATED( StoreSol ) ) THEN
       StorePerm    => StoreSol % Perm
       StoreValues  => StoreSol % Values
       StoreDOFs    = StoreSol % DOFs
 ELSE
       CALL Fatal(SolverName, 'Could not find storage variable')
 END IF

 TargetSol => VariableGet( Solver % Mesh % Variables, VarTarget )
 IF ( ASSOCIATED( TargetSol ) ) THEN
       TargetPerm    => TargetSol % Perm
       TargetValues  => TargetSol % Values
       TargetDOFs    = TargetSol % DOFs
 ELSE
       CALL Fatal(SolverName, 'Could not find target variable')
 END IF

 IF(StoreDOFs > 1 .OR. TargetDOFs > 1) THEN
    CALL Fatal(SolverName, 'Variables with DOF > 1 not yet supported')
 END IF

 IF (.NOT. done) THEN
    DO i=1,Solver % NumberOFActiveElements

      CurrentElement => GetActiveElement(i)
      NodeIndexes => CurrentElement % NodeIndexes

      DO k=1, GetElementNOFNodes(CurrentElement)
           TargetValues(TargetPerm(NodeIndexes(k))) =  StoreValues(StorePerm(NodeIndexes(k))) 
      END DO
    END DO
    WRITE(Message,'(a)') '---------------------------------------------------------------'
    CALL Info(SolverName,Message, Level=3)
    CALL Info(SolverName,'Create store:..........done', Level=3)
    WRITE(Message,'(a)') '---------------------------------------------------------------'
    CALL Info(SolverName,Message, Level=3)

    done = .TRUE.
 ELSE
    CALL Info(SolverName, 'No ops....')
 END IF

!------------------------------------------------------------------------------
END SUBROUTINE createStore
!------------------------------------------------------------------------------
