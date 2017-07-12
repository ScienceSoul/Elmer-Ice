!----------------------------------------------------------------------------------    
RECURSIVE SUBROUTINE LineOutputs( Model,Solver,dt,TransientSimulation )
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
 
 INTEGER :: nvalue, i, j, k, DIM, istat, FlowDOFs
 TYPE(Element_t), POINTER :: CurrentElement
 TYPE(ValueList_t),POINTER :: SolverParams
 TYPE(ValueList_t), POINTER :: Material
 REAL(KIND=dp), POINTER :: FlowValues(:), DepthValues(:), HeightValues(:)
 REAL(KIND=dp) :: velocity
 TYPE(Variable_t), POINTER :: FlowSol, DepthSol, HeightSol
 INTEGER, POINTER :: FlowPerm(:), DepthPerm(:), HeightPerm(:), NodeIndexes(:)
 TYPE(Nodes_t) :: ElementNodes
 CHARACTER(LEN=MAX_NAME_LEN) :: SolverName, FileNameSurf, FileNameBase
 LOGICAL :: AllocationsDone = .FALSE., Found, done = .FALSE.

 SAVE done, ElementNodes

 SolverName = 'LimitLines'

 SolverParams => Solver % Values

 dim = CoordinateSystemDimension()  

 FileNameSurf = GetString(SolverParams,'LineOut Surface File Name',Found )
 IF(.NOT. Found) THEN
    CALL Fatal(SolverName, 'LineOut Surface File name not found')
 END IF

 FileNameBase = GetString(SolverParams,'LineOut Base File Name',Found )
 IF(.NOT. Found) THEN
    CALL Fatal(SolverName, 'LineOut Base File name not found')
 END IF


 OPEN (12, FILE=FileNameSurf,POSITION='APPEND')
 OPEN (13, FILE=FileNameBase,POSITION='APPEND')

 FlowSol => VariableGet( Solver % Mesh % Variables, 'Flow Solution' )
 IF ( ASSOCIATED( FlowSol ) ) THEN
       FlowPerm    => FlowSol % Perm
       FlowValues  => FlowSol % Values
       FlowDOFs = FlowSol % DOFs
 ELSE
       CALL Fatal(SolverName, 'Could not find velocity field variable')
 END IF

 DepthSol => VariableGet( Solver % Mesh % Variables, 'Depth' )
 IF ( ASSOCIATED( DepthSol ) ) THEN
       DepthPerm    => DepthSol % Perm
       DepthValues  => DepthSol % Values
 ELSE
       CALL Fatal(SolverName, 'Could not find Depth field variable')
 END IF

 HeightSol => VariableGet( Solver % Mesh % Variables, 'Height' )
 IF ( ASSOCIATED(  HeightSol ) ) THEN
       HeightPerm    =>  HeightSol % Perm
       HeightValues  =>  HeightSol % Values
 ELSE
       CALL Fatal(SolverName, 'Could not find Height field variable')
 END IF

 DO i=1,Solver % NumberOFActiveElements
    CurrentElement => GetActiveElement(i)
    NodeIndexes => CurrentElement % NodeIndexes
    CALL GetElementNodes( ElementNodes )
    DO k=1, GetElementNOFNodes(CurrentElement)
       j = FlowPerm(NodeIndexes(k))
       velocity = SQRT(FlowValues(FlowDOFs*(j-1)+1)**2.0_dp + FlowValues(FlowDOFs*(j-1)+2)**2.0_dp + &
                       FlowValues(FlowDOFs*(j-1)+3)**2.0_dp)
       IF(DepthValues((DepthPerm(NodeIndexes(k))-1)+1) == 0.0_dp) THEN
          write(12,'(e13.5,2x,e15.8,2x,e15.8,2x,e15.8)') ElementNodes % x(k), ElementNodes % y(k), &
                                                         ElementNodes % z(k), velocity
       ELSE IF(HeightValues((HeightPerm(NodeIndexes(k))-1)+1) == 0.0_dp) THEN
          write(13,'(e13.5,2x,e15.8,2x,e15.8,2x,e15.8)') ElementNodes % x(k), ElementNodes % y(k), &
                                                         ElementNodes % y(k), velocity
       END IF
    END DO
 END DO

 close(12)
 close(13)
!------------------------------------------------------------------------------
END SUBROUTINE LineOutputs
!------------------------------------------------------------------------------
