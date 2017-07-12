!----------------------------------------------------------------------------------    
RECURSIVE SUBROUTINE LimitTemp( Model,Solver,dt,TransientSimulation )
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

 INTEGER :: i, j, k, l, m, DIM, TempDOFs, TempHomoDOFs, FlowDOFs, count
 TYPE(Element_t), POINTER :: CurrentElement
 TYPE(Variable_t), POINTER :: TempSol, TempHomoSol, FlowSol
 INTEGER, POINTER :: TempPerm(:), TempHomoPerm(:), FlowPerm(:), NodeIndexes(:)
 REAL(KIND=dp), POINTER :: TempValues(:), TempHomoValues(:), FlowValues(:)
 REAL(KIND=dp) :: beta, pressureMelting, pressure
 CHARACTER(LEN=MAX_NAME_LEN) :: SolverName
 LOGICAL :: AllocationsDone = .FALSE., Found

 SAVE AllocationsDone, DIM

 SolverName = 'LimitTemp'

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

 TempSol => VariableGet( Solver % Mesh % Variables, 'Temp' )
 IF ( ASSOCIATED( TempSol ) ) THEN
       TempPerm    => TempSol % Perm
       TempValues  => TempSol % Values
       TempDOFs = TempSol % DOFs
 ELSE
       CALL Fatal(SolverName, 'Could not find temp field variable')
 END IF

 TempHomoSol => VariableGet( Solver % Mesh % Variables, 'Temp Homologous' )
 IF ( ASSOCIATED( TempHomoSol ) ) THEN
       TempHomoPerm    => TempHomoSol % Perm
       TempHomoValues  => TempHomoSol % Values
       TempHomoDOFs = TempHomoSol % DOFs
 ELSE
       CALL Fatal(SolverName, 'Could not find temp homologous field variable')
 END IF

 beta = 9.8E-08*1.0E06;

 count = 0

 DO i=1,Solver % NumberOFActiveElements

   CurrentElement => GetActiveElement(i)
   NodeIndexes => CurrentElement % NodeIndexes

   DO k=1, GetElementNOFNodes(CurrentElement)

      j = TempPerm(NodeIndexes(k))
      l =  FlowPerm(NodeIndexes(k))
      m =  TempHomoPerm(NodeIndexes(k))

      pressure = FlowValues(FlowDOFs*(l-1)+4)
      IF (pressure < 0.0) pressure = 0.0

      pressureMelting = 273.15 - ( beta*pressure)

      IF (TempValues(TempDOFs*(j-1)+1) > 273.15) THEN
         TempValues(TempDOFs*(j-1)+1) = 273.15
         TempHomoValues(TempHomoDOFs*(m-1)+1) = TempValues(TempDOFs*(j-1)+1) - pressureMelting
         count = count+1
      ELSE IF (TempValues(TempDOFs*(j-1)+1) < 210.0) THEN
         TempValues(TempDOFs*(j-1)+1) = 210.0
         TempHomoValues(TempHomoDOFs*(m-1)+1) =  TempValues(TempDOFs*(j-1)+1) - pressureMelting
         count = count+1
      END IF

   END DO

 END DO

WRITE(Message,'(a)') '----------------------------------------------------------'
CALL Info(SolverName,Message, Level=3)
WRITE(Message,'(a,I4,a)') 'Limited ', count, ' nodes for temperature'
CALL Info(SolverName,Message, Level=3)
WRITE(Message,'(a)') '---------------------------------------------------------------'
CALL Info(SolverName,Message, Level=3)

!------------------------------------------------------------------------------
END SUBROUTINE LimitTemp
!------------------------------------------------------------------------------

