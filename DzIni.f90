FUNCTION DzIni ( Model, nodenumber, x) RESULT(Dz)
   USE types
   USE CoordinateSystems
   USE SolverUtils
   USE ElementDescription
   USE DefUtils
   IMPLICIT NONE
   TYPE(Model_t) :: Model
   TYPE(Solver_t), TARGET :: Solver
   INTEGER :: nodenumber,  NMAX, i
   REAL(KIND=dp) :: x,   Dz
   REAL(KIND=dp), ALLOCATABLE :: Dz0(:)
   LOGICAL :: FirstTime=.True.

   SAVE FirstTime
   SAVE Dz0

   IF (FirstTime) THEN
          FirstTime=.False.
          NMAX = Model % NumberOfNodes
          ALLOCATE(Dz0(NMAX))
          Do i = 1, NMAX
          Dz0(i) = Model % Nodes % z (i)
          EndDo
   END IF

       Dz = Dz0(nodenumber)

END FUNCTION DzIni
