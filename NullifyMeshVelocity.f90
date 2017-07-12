!---------------------------------------------------------------------------------------------------------    
RECURSIVE SUBROUTINE NullifyMeshVelocity( Model,Solver,dt,TransientSimulation )
!---------------------------------------------------------------------------------------------------------- 

  USE DefUtils
  USE Types

  IMPLICIT NONE
!------------------------------------------------------------------------------
!******************************************************************************
!
!  Nullify Mesh Velicity variable
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
  TYPE(Element_t),POINTER :: CurrentElement
  INTEGER, POINTER :: NodeIndexes(:)
  TYPE(ValueList_t), POINTER :: SolverParams

  INTEGER :: nvalue, i, j, k, N, MeshVeloDOFs, DIM
  REAL(KIND=dp), POINTER :: MeshVeloValues(:)
  TYPE(Variable_t), POINTER :: MeshVeloSol
  INTEGER, POINTER :: MeshVeloPerm(:)
  LOGICAL :: AllocationsDone = .FALSE., Found

  SAVE AllocationsDone
 

  !------------------------------------------------------------------------------
  !     Allocate some permanent storage, this is done first time only
  !------------------------------------------------------------------------------
  IF ( .NOT. AllocationsDone .OR. Solver % Mesh % Changed) THEN
      N = Model % MaxElementNodes
      DIM = CoordinateSystemDimension()
      AllocationsDone = .TRUE. 
  END IF


  MeshVeloSol => VariableGet( Solver % Mesh % Variables, &
            'Mesh Velocity' )

  IF ( ASSOCIATED( MeshVeloSol ) ) THEN
       MeshVeloPerm    => MeshVeloSol % Perm
       MeshVeloValues  => MeshVeloSol % Values
       MeshVeloDOFs = MeshVeloSol % DOFs
  END IF

  DO i=1,Solver % Mesh % NumberOfNodes

     k = MeshVeloDOFs
     DO j=1,k
       MeshVeloValues(k*(MeshVeloPerm(i)-1)+j) = 0.0
     END DO
  
  END DO

  WRITE(Message,'(a)') '------------------------------------------------------'
  CALL Info('NullifyMeshVelocity',Message, Level=3)
  CALL Info('NullifyMeshVelocity','Nullify Mesh Velocity:..........done', Level=3)
  WRITE(Message,'(a)') '------------------------------------------------------'
  CALL Info('NullifyMeshVelocity',Message, Level=3)

!------------------------------------------------------------------------------
END SUBROUTINE NullifyMeshVelocity
!------------------------------------------------------------------------------
