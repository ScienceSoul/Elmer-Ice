!*********************************************************************************************/
FUNCTION getViscosityFactor (Model, nodenumber, y) RESULT(viscFact)

   USE types
   USE CoordinateSystems
   USE SolverUtils
   USE ElementDescription
   USE DefUtils

   IMPLICIT NONE

   TYPE(Model_t) :: Model
   REAL (KIND=dp) :: y             
   INTEGER :: nodenumber

!------------------------------------------------------------------------------
!    Local variables
!------------------------------------------------------------------------------

 INTEGER :: nvalue, i, j, k, ViscDOFs, DIM
 REAL(KIND=dp), POINTER :: ViscValues(:)
 TYPE(Variable_t), POINTER :: ViscSol
 INTEGER, POINTER :: ViscPerm(:)
 REAL(KIND=dp) :: viscFact

 ViscSol => VariableGet( Model % Variables, 'viscosity factor' )
 IF ( ASSOCIATED( ViscSol ) ) THEN
      ViscPerm    => ViscSol % Perm
      ViscValues  => ViscSol % Values
      ViscDOFs = ViscSol % DOFs
 ELSE
     CALL FATAL('getViscosityFactor', 'Could not find viscosity factor variable')
 END IF

 viscFact = ViscValues(ViscPerm(nodenumber))

 END FUNCTION getViscosityFactor
