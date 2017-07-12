FUNCTION fortuin_temp_range( Model, node, zs) RESULT(temp)
    USE Types
    USE CoordinateSystems
    USE SolverUtils
    USe ElementDescription

    !-------------------------------------------------------------------
    IMPLICIT NONE
    !-------------------------external variables-----------------------
    TYPE(Model_t) :: Model 
    INTEGER :: node
    TYPE(Variable_t), POINTER :: phiSol
    REAL(KIND=dp), POINTER :: phii(:)
    INTEGER, POINTER :: phiPerm(:)
    REAL(KIND=dp) :: temp, zs
    REAL(KIND=dp) :: pi_180_inv = 180.0_dp/pi
    INTEGER ::  phiDOFs
    REAL(KIND=dp) :: theta_ma, c_ma, gamma_ma, &
            theta_ma_1, c_ma_1, gamma_ma_1, &
            theta_ma_2, c_ma_2, gamma_ma_2, &
            theta_ma_3, c_ma_3, gamma_ma_3, &
            zs_sep_1, zs_sep_2

    ! Parameterisation by Fortuin and Oerlemans
    !  (1990), separately for three different
    !  elevation ranges

    zs_sep_1   =  200.0_dp
    zs_sep_2   = 1500.0_dp

    theta_ma_1 =  49.642_dp
    gamma_ma_1 =   0.0_dp
    c_ma_1     =  -0.943_dp

    theta_ma_2 =  36.689_dp
    gamma_ma_2 =  -5.102e-03_dp
    c_ma_2     =  -0.725_dp

    theta_ma_3 =   7.405_dp
    gamma_ma_3 = -14.285e-03_dp
    c_ma_3     =  -0.180_dp

    phiSol => VariableGet( Model % Variables, 'Latitude' )
    IF ( ASSOCIATED( phiSol ) ) THEN
       phiPerm    => phiSol % Perm
       phii       => phiSol % Values
       phiDOFs = phiSol % DOFs
    ELSE
       CALL FATAL('fausto_temp', ' Could not find Latitude field variable')
    END IF

    IF ( zs <= zs_sep_1 ) THEN
      temp = theta_ma_1 + gamma_ma_1*zs &
                        + c_ma_1*abs(phii(phiDOFs*(phiPerm(node)-1)+1)) !*pi_180_inv
    ELSE IF ( zs <= zs_sep_2 ) THEN
      temp = theta_ma_2 + gamma_ma_2*zs &
                        + c_ma_2*abs(phii(phiDOFs*(phiPerm(node)-1)+1)) !*pi_180_inv
    ELSE
      temp = theta_ma_3 + gamma_ma_3*zs &
                        + c_ma_3*abs(phii(phiDOFs*(phiPerm(node)-1)+1)) !*pi_180_inv
    END IF 

    temp = temp + 273.15

END FUNCTION fortuin_temp_range

