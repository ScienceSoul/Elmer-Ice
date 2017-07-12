SUBROUTINE ComputeNormalSolver( Model, Solver, dt, TransientSimulation )
!DEC$ATTRIBUTES DLLEXPORT :: ComputeNormalSolver
!------------------------------------------------------------------------------
!******************************************************************************
!
!  ARGUMENTS:
!
!  TYPE(Model_t) :: Model,  
!     INPUT: All model information (mesh, materials, BCs, etc...)
!
!  TYPE(Solver_t) :: Solver
!     INPUT: Linear & nonlinear equation solver options
!
!  REAL(KIND=dp) :: dt,
!     INPUT: Timestep size for time dependent simulations
!
!  LOGICAL :: TransientSimulation
!     INPUT: Steady state or transient simulation
!
!******************************************************************************
  USE DefUtils

  IMPLICIT NONE
!------------------------------------------------------------------------------
  TYPE(Solver_t), TARGET :: Solver
  TYPE(Model_t) :: Model

  REAL(KIND=dp) :: dt
  LOGICAL :: TransientSimulation

!------------------------------------------------------------------------------
! Local variables
!------------------------------------------------------------------------------
  TYPE(Element_t), POINTER :: Element
  TYPE(Variable_t), POINTER :: NormalSolution

  INTEGER :: i, j, k, n, t, DIM
  REAL(KIND=dp) :: u, v, w, s 
  REAL(KIND=dp), POINTER :: Nvector(:)
  INTEGER, POINTER :: Permutation(:)

  TYPE(Nodes_t) :: Nodes
  REAL(KIND=dp) :: Bu, Bv, Normal(3) 

  SAVE Nodes

!---------------------------------------
! Setup pointer to the current solution:
!---------------------------------------
  NormalSolution => VariableGet( Solver % Mesh % Variables, 'Normal Vector' ) 
  IF ( ASSOCIATED( NormalSolution ) ) THEN 
     Nvector => NormalSolution % Values
     Permutation => NormalSolution % Perm
  ELSE
     PRINT *,'FATAL: Unable to set pointer to the current solution'
     STOP
  END IF

  PRINT *,'-----------------------------------'
  PRINT *,' Computing Normal Vector for Nodes'
  PRINT *,'-----------------------------------'


 DIM = CoordinateSystemDimension()

  DO t=1,Solver % Mesh % NumberOfBoundaryElements
     Element => GetBoundaryElement( t )
     n = GetElementNOFNodes( Element )
     IF (n==1) CYCLE

     CALL GetElementNodes( Nodes , Element )

     DO i = 1,n
        j = Element % NodeIndexes( i )
        k = Permutation(j)

              Bu = Element % Type % NodeU(i)
              IF ( Element % Type % Dimension > 1 ) THEN
                Bv = Element % Type % NodeV(i)
              ELSE
                Bv = 0.0D0
              END IF
        Normal = NormalVector(Element, Nodes, Bu, Bv, .TRUE.)
        Nvector(DIM*(k-1)+1:DIM*k) = Nvector(DIM*(k-1)+1:DIM*k) +& 
                               Normal(1:DIM) 
     END DO 
  END DO

  DO i=1,Model % NumberOfNodes
      k = Permutation(i)
      IF ( k > 0 ) THEN
        s = SQRT( SUM( Nvector(DIM*(k-1)+1:DIM*k)**2 ) )

        IF ( s /= 0.0D0 ) THEN
        Nvector(DIM*(k-1)+1:DIM*k) = Nvector(DIM*(k-1)+1:DIM*k)/s 
        END IF
      END IF
  END DO


  PRINT *,'Done!'

!------------------------------------------------------------------------------
END SUBROUTINE ComputeNormalSolver
!------------------------------------------------------------------------------
