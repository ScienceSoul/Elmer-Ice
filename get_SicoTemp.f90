FUNCTION get_SicoTemp( Model, nodenumber, x) RESULT(tempSico)
  USE Types
  USE CoordinateSystems
  USE SolverUtils
  USe ElementDescription
!-------------------------------------------------------------------
IMPLICIT NONE
!-------------------------external variables-----------------------
TYPE(Model_t) :: Model 
INTEGER :: nodenumber
INTEGER, POINTER :: Temp_AnomaliesPermutation(:)
REAL(KIND=dp) :: x, tempSico, variableStore
TYPE(Variable_t), POINTER :: Temp_AnomaliesSolution 

!Find variable for Temperature Anomalies 
Temp_AnomaliesSolution  => VariableGet(Model % Variables, 'TempAnomalies', .TRUE.)
IF ( ASSOCIATED( Temp_AnomaliesSolution) ) THEN 
Temp_AnomaliesPermutation  => Temp_AnomaliesSolution % Perm
ELSE 
CALL Fatal('get_SicoTemp', 'No variable TempAnomalies found') 
END IF 

 
variableStore = -54.3+Temp_AnomaliesSolution % Values(Temp_AnomaliesPermutation(nodenumber))

tempSico = variableStore+273.15

!WRITE(*,*) tempSico

END FUNCTION get_SicoTemp
