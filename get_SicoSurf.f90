FUNCTION get_SicoSurf( Model, nodenumber, x) RESULT(surfSico)
  USE Types
  USE CoordinateSystems
  USE SolverUtils
  USe ElementDescription
!-------------------------------------------------------------------
IMPLICIT NONE
!-------------------------external variables-----------------------
TYPE(Model_t) :: Model 
INTEGER :: nodenumber, NMAX
INTEGER, POINTER :: Surf_AnomaliesPermutation(:)
REAL(KIND=dp) :: x, surfSico
TYPE(Variable_t), POINTER :: Surf_AnomaliesSolution 

!Find variable for Surface Anomalies 
Surf_AnomaliesSolution  => VariableGet(Model % Variables, 'SurfAnomalies', .TRUE.)
IF ( ASSOCIATED( Surf_AnomaliesSolution) ) THEN 
Surf_AnomaliesPermutation  => Surf_AnomaliesSolution % Perm
ELSE 
CALL Fatal('get_SicoSurf', 'No variable SurfAnomalies found') 
END IF 

surfSico = Surf_AnomaliesSolution % Values(Surf_AnomaliesPermutation(nodenumber))

!WRITE(*,*) surfSico

END FUNCTION get_SicoSurf
