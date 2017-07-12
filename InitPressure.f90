!----------------------------------------------------------------------------------    
RECURSIVE SUBROUTINE InitPressure( Model,Solver,dt,TransientSimulation )
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

 INTEGER :: nvalue, i, j, k, FlowDOFs, DepthDOFs, DIM
 REAL(KIND=dp), POINTER :: FlowValues(:), DepthValues(:)
 TYPE(Variable_t), POINTER :: FlowSol, DepthSol
 INTEGER, POINTER :: FlowPerm(:), DepthPerm(:)
 REAL(KIND=dp) :: rho, g

 rho = 918.0
 g = 9.81 

 FlowSol => VariableGet( Solver % Mesh % Variables, 'flow solution' )
 IF ( ASSOCIATED( FlowSol ) ) THEN
       FlowPerm    => FlowSol % Perm
       FlowValues  => FlowSol % Values
       FlowDOFs = FlowSol % DOFs
 ELSE
       CALL Info('InitPressure', 'Could not find velocity field variable', Level=4)
 END IF

 DepthSol => VariableGet( Solver % Mesh % Variables, 'Depth' )
 IF ( ASSOCIATED( DepthSol ) ) THEN
        DepthPerm => DepthSol % Perm
        DepthValues => DepthSol % Values
        DepthDOFs = DepthSol % DOFs
 ELSE
 WRITE(Message,'(a,I6)') 'Could not find Depth field variable'
 CALL FATAL('InitAge',Message)
 END IF

 dim = CoordinateSystemDimension()  

 DO i=1,Solver % Mesh % NumberOfNodes
   
   j = FlowPerm(i)
   k = DepthPerm(i)
   FlowValues(FlowDOFs*(j-1)+(DIM+1)) = rho*g*DepthValues(DepthDOFs*(k-1)+1)
 END DO

WRITE(Message,'(a)') '------------------------------------------------------'
CALL Info('InitAge',Message, Level=3)
CALL Info('InitAge','Initialize values for pressure:..........done', Level=3)
WRITE(Message,'(a)') '------------------------------------------------------'
CALL Info('InitAge',Message, Level=3)

!------------------------------------------------------------------------------
END SUBROUTINE InitPressure
!------------------------------------------------------------------------------

