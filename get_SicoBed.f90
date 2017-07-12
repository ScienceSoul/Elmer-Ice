FUNCTION get_SicoBed( Model, nodenumber, x) RESULT(bedSico)
  USE Types
  USE CoordinateSystems
  USE SolverUtils
  USe ElementDescription
!-------------------------------------------------------------------
IMPLICIT NONE
!-------------------------external variables-----------------------
TYPE(Model_t) :: Model 
INTEGER :: nodenumber
INTEGER, POINTER :: Bed_AnomaliesPermutation(:)
REAL(KIND=dp) :: x, bedSico
TYPE(Variable_t), POINTER :: Bed_AnomaliesSolution 

!Find variable for Surface Anomalies 
Bed_AnomaliesSolution  => VariableGet(Model % Variables, 'BedAnomalies', .TRUE.)
IF ( ASSOCIATED( Bed_AnomaliesSolution) ) THEN 
Bed_AnomaliesPermutation  => Bed_AnomaliesSolution % Perm
ELSE 
CALL Fatal('get_SicoBed', 'No variable BedAnomalies found') 
END IF 

bedSico = Bed_AnomaliesSolution % Values(Bed_AnomaliesPermutation(nodenumber))

!WRITE(*,*) bedSico

END FUNCTION get_SicoBed
