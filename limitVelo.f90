!----------------------------------------------------------------------------------    
RECURSIVE SUBROUTINE LimitVelo( Model,Solver,dt,TransientSimulation )
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

 INTEGER :: i, j, k, DIM, FlowDOFs, count
 TYPE(Element_t), POINTER :: CurrentElement
 TYPE(Variable_t), POINTER :: FlowSol
 INTEGER, POINTER :: FlowPerm(:), NodeIndexes(:)
 REAL(KIND=dp), POINTER :: FlowValues(:)
 REAL(KIND=dp) :: UpperLimit, LowerLimit
 CHARACTER(LEN=MAX_NAME_LEN) :: SolverName
 LOGICAL :: AllocationsDone = .FALSE., Found

 SAVE AllocationsDone, DIM

 SolverName = 'LimitVelo'

 IF ( .NOT. AllocationsDone .OR. Solver % Mesh % Changed ) THEN
   
    dim = CoordinateSystemDimension()  
        
    AllocationsDone = .TRUE.

 END IF

 FlowSol => VariableGet( Solver % Mesh % Variables, 'flow solution' )
 IF ( ASSOCIATED( FlowSol ) ) THEN
       FlowPerm    => FlowSol % Perm
       FlowValues  => FlowSol % Values
       FlowDOFs = FlowSol % DOFs
 ELSE
       CALL Fatal(SolverName, 'Could not find velocity field variable')
 END IF

 UpperLImit = GetConstReal(Solver % Values, 'Upper Limit', Found)
 IF (.NOT. Found) THEN
    CALL FATAL(SolverName,' Parameter Upper Limit not found .')
 END IF

 LowerLImit = GetConstReal(Solver % Values, 'Lower Limit', Found)
 IF (.NOT. Found) THEN
    CALL FATAL(SolverName,' Parameter Lower Limit not found .')
 END IF

 count = 0

 DO i=1,Solver % NumberOFActiveElements

   CurrentElement => GetActiveElement(i)
   NodeIndexes => CurrentElement % NodeIndexes

   DO k=1, GetElementNOFNodes(CurrentElement)

      j = FlowPerm(NodeIndexes(k))

      !vx
      IF (FlowValues(FlowDOFs*(j-1)+1) > UpperLimit) THEN
         FlowValues(FlowDOFs*(j-1)+1) = UpperLimit
         count = count+1
      ELSE IF (FlowValues(FlowDOFs*(j-1)+1) < LowerLimit) THEN
         FlowValues(FlowDOFs*(j-1)+1) = LowerLimit
         count = count+1
      END IF

      !vy
      IF (FlowValues(FlowDOFs*(j-1)+2) > UpperLimit) THEN
         FlowValues(FlowDOFs*(j-1)+2) = UpperLImit
         count = count+1
      ELSE IF (FlowValues(FlowDOFs*(j-1)+2) < LowerLimit) THEN
         FlowValues(FlowDOFs*(j-1)+2) = LowerLImit
         count = count+1
      END IF

      !vz
      IF (FlowValues(FlowDOFs*(j-1)+3) > UpperLimit) THEN
         FlowValues(FlowDOFs*(j-1)+3) = UpperLimit
         count = count+1
      ELSE IF (FlowValues(FlowDOFs*(j-1)+3) < LowerLimit) THEN
         FlowValues(FlowDOFs*(j-1)+3) = LowerLimit
         count = count+1
      END IF

   END DO

 END DO

WRITE(Message,'(a)') '----------------------------------------------------------'
CALL Info(SolverName,Message, Level=3)
WRITE(Message,'(a,I4,a)') 'Limit ', count, ' nodes for Velocities'
CALL Info(SolverName,Message, Level=3)
WRITE(Message,'(a)') '---------------------------------------------------------------'
CALL Info(SolverName,Message, Level=3)

!------------------------------------------------------------------------------
END SUBROUTINE LimitVelo
!------------------------------------------------------------------------------
