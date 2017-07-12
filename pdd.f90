RECURSIVE SUBROUTINE pdd( Model,Solver,Timestep,TransientSimulation )

  USE DefUtils
  USE MaterialModels
  USE Types
  IMPLICIT NONE

  !------------------------------------------------------------------------------
  !    External variables
  !------------------------------------------------------------------------------
  TYPE(Model_t)  :: Model
  TYPE(Solver_t), TARGET :: Solver
  LOGICAL :: TransientSimulation
  REAL(KIND=dp) :: Timestep
  !------------------------------------------------------------------------------
  !    Local variables
  !------------------------------------------------------------------------------
  INTEGER :: i, j, n, m, k, l, DIM, TimeDOFs, zsDOFs, precip_ma_presDOFs, zs_refDOFs, phiDOFs, lambdaDOFs, &
           asDOFs, tempDOFs, annualtempDOFs, julytempDOFs, ar4_precip_DOFS, ar4_annualtemp_DOFs, &
           ar4_julytemp_DOFs, istat, ios, body_id, old_body = -1

  INTEGER :: i_gr, i_kl
  INTEGER :: TSURFACE, ABLSURFACE, SOLID_PRECIP, ACCSURFACE, PRECIP_ANOM_INTERPOL, TEMP_PRESENT_PARA, GRID
  INTEGER :: grip_time_min, grip_time_stp, grip_time_max, ndata_grip
  INTEGER :: gi_time_min, gi_time_stp, gi_time_max, ndata_gi
  TYPE(ValueList_t), POINTER :: Material
  TYPE(Element_t), POINTER :: CurrentElement
  
  TYPE(Variable_t), POINTER :: TimeSol, zsSol, precip_ma_presSol, phiSol, lambdaSol, zs_refSol, asSol, tempSol, annualTempSol, &
                        julytempSol, ar4_precipSol, ar4_annualtempSol, ar4_julytempSol

  REAL(KIND=dp), POINTER :: Time(:), zs(:), precip_ma_present(:), phii(:), lambda(:), zs_ref(:), as(:), temp(:), annualtemp(:), &
                         julytemp(:), ar4_precip(:), ar4_annualtemp(:), ar4_julytemp(:)

  INTEGER, POINTER :: TimePerm(:), zsPerm(:), precip_ma_presPerm(:), phiPerm(:), lambdaPerm(:), zs_refPerm(:), &
                      NodeIndexes(:), asPerm(:), tempPerm(:), annualtempPerm(:), julytempPerm(:), ar4_precip_Perm(:), &
                      ar4_annualtemp_Perm(:), ar4_julytemp_Perm(:)

  REAL(KIND=dp) :: delta_ts, glac_index
  REAL(KIND=dp) :: time_gr, time_kl
  REAL(KIND=dp), ALLOCATABLE :: temp_ma(:), temp_mm(:,:), temp_ampl(:), precip(:,:), precip_present(:,:), &
                            temp_ma_lgm_anom(:), temp_mj_lgm_anom(:), precip_lgm_anom(:,:), &
                            gamma_precip_lgm_anom(:,:), accum(:), evap(:), runoff(:), as_perp(:), &
                            temp_s(:), snowfall(:), rainfall(:), ET(:), melt(:), melt_star(:), &
                            glacial_index(:), ar4_time1_precip(:), ar4_time1_annualtemp(:), ar4_time1_julytemp(:), &
                            anomalies_precip(:), anomaliesPrecipAtT94(:), anomalies_annualtemp(:), anomalies_julytemp(:), & 
                            anamaliesAnnualTempAtT94(:), anamaliesJulyTempAtT94(:)

  REAL(KIND=dp) :: theta_ma, c_ma, kappa_ma, gamma_ma, &
                   theta_mj, c_mj, kappa_mj, gamma_mj, zs_sep_1, zs_sep_2, theta_ma_1, gamma_ma_1, &
                   c_ma_1, theta_ma_2, gamma_ma_2, c_ma_2, c_ma_3, theta_ma_3, gamma_ma_3
  REAL(KIND=dp), PARAMETER :: GRIP_TEMP_FACT = 1.0d0
  REAL(KIND=dp), PARAMETER :: DELTA_TS0 = 0.0d0
  REAL(KIND=dp), PARAMETER :: SINE_AMPLIT = 10.0d0
  REAL(KIND=dp), PARAMETER :: pi_180 = pi/180.0_dp
  REAL(KIND=dp), PARAMETER :: pi_180_inv = 180.0_dp/pi
  REAL(KIND=dp), PARAMETER :: eps = 1.0e-05_dp
  REAL(KIND=dp) :: sine_factor
  REAL(KIND=dp), DIMENSION(12) :: temp_mm_help
  REAL(KIND=dp) :: temp_jja_help
  Real(KIND=dp) :: griptemp(0:10000)
  REAL(KIND=dp) :: gamma_p, zs_thresh, temp_rain, temp_snow, &
                  inv_delta_temp_rain_snow, coeff(0:5), inv_sqrt2_s_stat, &
                  precip_fact, frac_solid
  REAL(KIND=dp) :: s_stat, phi_sep, temp_lt, temp_ht, inv_delta_temp_ht_lt, &
                  beta1_lt, beta1_ht, beta2_lt, beta2_ht, beta1, beta2, beta_ice, beta_snow, Pmax, &
                  mu, lambda_lti, temp_lti, sTime
  REAL(KIND=dp) :: S_STAT_0      ! Standard deviation of the air termperature for the degree-day model
  REAL(KIND=dp) :: PHI_SEP_0     ! Separation latitude for the computation of the degree-day
                                ! factors beta1 and beta2: Equatorward of phi_sep, only the
                                ! high-temperature values are used, whereas poleward of phi_sep,
                                ! beta1 and beta2 are temperature-dependent
  REAL(KIND=dp) :: PMAX_0        ! Saturation factor for the formation of superimposed ice
  REAL(KIND=dp) :: MU_0          ! Firn-warming correction
  REAL(KIND=dp) :: RHO_W         ! Density of pure water
  REAL(KIND=dp) :: RHO           ! Density of ice
  REAL(KIND=dp) :: BETA1_LT_0    ! Degree-day factor for snow at low summer temperatures
  REAL(KIND=dp) :: BETA1_HT_0    ! Degree-day factor for snow at high summer temperatures
  REAL(KIND=dp) :: BETA2_LT_0    ! Degree-day factor for ice at low summer temperatures
  REAL(KIND=dp) :: BETA2_HT_0    ! Degree-day factor for ice at high summer temperatures
  REAL(KIND=dp) :: BETA_SNOW_0   ! Degree-day factor for snow
  REAL(KIND=dp) :: BETA_ICE_0    ! Degree-day factor for ice
  REAL(KIND=dp), PARAMETER :: LAMBDA_LTI_C = 500.0d0
  REAL(KIND=dp), PARAMETER :: TEMP_LTI_C = 2.0d0
  REAL(KIND=dp), PARAMETER :: ACCFACT = 1.0d0
  REAL(KIND=dp), PARAMETER :: GAMMA_S = 0.070458d0
  REAL(KIND=dp), parameter :: &
          inv_twelve = 1.0_dp/12.0_dp, one_third = 1.0_dp/3.0_dp
  REAL(KIND=dp) :: YEAR_SEC
  REAL(KIND=dp) :: d_dummy
  REAL(KIND=dp) :: anomaliesScaling
  LOGICAL :: AllocationsDone = .FALSE., Found
  LOGICAL, ALLOCATABLE :: watchdog(:)
  LOGICAL :: Anomalies, anomaliesParam, lapsrate, convert, TarasovPDD = .FALSE., HuybrechtsPDD = .FALSE.
  CHARACTER(LEN=MAX_NAME_LEN) :: SolverName, ch_dummy, INPATH, GRIP_TEMP_FILE, GLAC_IND_FILE, &
                              TEMP_MA_ANOM_FILE, TEMP_MJ_ANOM_FILE, PrecipMode, TempParam, continent
   
  SAVE AllocationsDone, SolverName, old_body, TSURFACE, ABLSURFACE, SOLID_PRECIP, &
      ACCSURFACE, PRECIP_ANOM_INTERPOL, TEMP_PRESENT_PARA, GRID, YEAR_SEC, ndata_grip, grip_time_min, grip_time_stp, &
      grip_time_max, griptemp, gi_time_min, gi_time_stp, gi_time_max, ndata_gi, glacial_index, &
      S_STAT_0, PHI_SEP_0, PMAX_0, MU_0, RHO_W, RHO, BETA1_LT_0, BETA1_HT_0, BETA2_LT_0, BETA2_HT_0, &
      BETA_SNOW_0, BETA_ICE_0, temp_ma, temp_mm, temp_ampl, precip, precip_present, temp_ma_lgm_anom, temp_mj_lgm_anom, &
      precip_lgm_anom, gamma_precip_lgm_anom, accum, evap, runoff, as_perp, temp_s, snowfall, rainfall, &
      ET, melt, melt_star, watchdog, PrecipMode, continent, TempParam, Anomalies, anomaliesParam, anomaliesScaling, & 
      ar4_time1_precip, ar4_time1_annualtemp, ar4_time1_julytemp, &
      anomalies_precip, anomaliesPrecipAtT94, anomalies_annualtemp, anomalies_julytemp, anamaliesAnnualTempAtT94, & 
      anamaliesJulyTempAtT94, lapsrate, convert, sTime

  SolverName = 'PDDSolver' 
 
  IF ( .NOT. AllocationsDone .OR. Solver % Mesh % Changed ) THEN

        DIM = CoordinateSystemDimension()
        M = Model % Mesh % NumberOfNodes

        IF ( AllocationsDone ) THEN
           DEALLOCATE(                     &
               temp_ma,                    &
               temp_mm,                    &
               temp_ampl,                  &
               precip,                     &
               precip_present,             &
               temp_ma_lgm_anom,           &
               temp_mj_lgm_anom,           &
               precip_lgm_anom,            &
               gamma_precip_lgm_anom,      &
               accum,                      &
               evap,                       &
               runoff,                     &
               as_perp,                    &
               temp_s,                     &
               snowfall,                   &
               rainfall,                   &
               ET,                         &
               melt,                       &
               melt_star,                  &
               watchdog,                   &
               ar4_time1_precip,           &
               ar4_time1_annualtemp,       &
               ar4_time1_julytemp,         &
               anomalies_precip,           &
               anomaliesPrecipAtT94,       &
               anomalies_annualtemp,       &
               anomalies_julytemp,         &
               anamaliesAnnualTempAtT94,   &
               anamaliesJulyTempAtT94)              
        END IF                           
        
        ALLOCATE(                             &
               temp_ma(M),                  &
               temp_mm(M,12),               &
               temp_ampl(M),                &
               precip(M,0:12),              &
               precip_present(M,12),        &
               temp_ma_lgm_anom(M),         &
               temp_mj_lgm_anom(M),         &
               precip_lgm_anom(M,12),       &
               gamma_precip_lgm_anom(M,12), & 
               accum(M),                    &
               evap(M),                     &
               runoff(M),                   &
               as_perp(M),                  &
               temp_s(M),                   &      
               snowfall(M),                 &
               rainfall(M),                 &
               ET(M),                       &
               melt(M),                     &
               melt_star(M),                &
               watchdog(M),                 &
               ar4_time1_precip(M),         &
               ar4_time1_annualtemp(M),     &
               ar4_time1_julytemp(M),       &
               anomalies_precip(M),         &
               anomaliesPrecipAtT94(M),     &
               anomalies_annualtemp(M),     &
               anomalies_julytemp(M),       &
               anamaliesAnnualTempAtT94(M), &
               anamaliesJulyTempAtT94(M),   &
               STAT=istat )

        IF ( istat /= 0 ) THEN
           CALL FATAL( SolverName, 'Memory allocation error' )
        ELSE
           CALL INFO(SolverName, 'Memory allocation done', level=1 )
        END IF

        TSURFACE = GetInteger(Solver % Values, 'TSURFACE', Found)
        IF (.NOT. Found) THEN
            CALL FATAL(SolverName,' Parameter TSURFACE not found .')
        END IF

        ABLSURFACE = GetInteger(Solver % Values, 'ABLSURFACE', Found)
        IF (.NOT. Found) THEN
            CALL FATAL(SolverName,' Parameter ABLSURFACE not found .')
        END IF

        SOLID_PRECIP = GetInteger(Solver % Values, 'SOLID_PRECIP', Found)
        IF (.NOT. Found) THEN
            CALL FATAL(SolverName,' Parameter SOLID_PRECIP not found .')
        END IF

        ACCSURFACE = GetInteger(Solver % Values, 'ACCSURFACE', Found)
        IF (.NOT. Found) THEN
            CALL FATAL(SolverName,' Parameter ACCSURFACE not found .')
        END IF

        PRECIP_ANOM_INTERPOL = GetInteger(Solver % Values, 'PRECIP_ANOM_INTERPOL', Found)
        IF (.NOT. Found) THEN
            CALL FATAL(SolverName,' Parameter PRECIP_ANOM_INTERPOL not found .')
        END IF

        TEMP_PRESENT_PARA = GetInteger(Solver % Values, 'TEMP_PRESENT_PARA', Found)
        IF (.NOT. Found) THEN
            CALL FATAL(SolverName,' Parameter TEMP_PRESENT_PARA not found .')
        END IF

        GRID = GetInteger(Solver % Values, 'GRID', Found)
        IF (.NOT. Found) THEN
            CALL FATAL(SolverName,' Parameter GRID not found .')
        END IF
        
        YEAR_SEC = GetConstReal(Solver % Values, 'YEAR_SEC', Found)
        IF (.NOT. Found) THEN
            CALL FATAL(SolverName,' Parameter YEAR_SEC not found .')
        END IF

        PrecipMode = GetString(Solver % Values, 'Precipitation input mode', Found)
        IF (.NOT. Found) THEN
            CALL FATAL(SolverName,' Parameter Precipitation input mode  not found .')
        END IF

        TempParam = GetString(Solver % Values, 'Temperature input mode', Found)
        IF (.NOT. Found) THEN
            CALL FATAL(SolverName,' Parameter Temperature input mode not found .')
        END IF

        continent = GetString(Solver % Values, 'Continent', Found)
        IF (.NOT. Found) THEN
            CALL FATAL(SolverName,' Parameter Continent not found .')
        END IF

        Anomalies = GetLogical(Solver % Values, 'Use Anomalies', Found)
        IF (.NOT. Found) THEN
            CALL FATAL(SolverName, 'Parameter: Use Anomalies not found')
        END IF
        IF(Anomalies) THEN

            CALL INFO(SolverName, 'Compute anomalies', level=3) 

            anomaliesParam =  GetLogical(Solver % Values, 'Anomalies parametrization', Found)
            IF (.NOT. Found) THEN
                anomaliesParam = .FALSE.
            END IF

            IF (anomaliesParam) THEN
               CALL INFO(SolverName, 'Use anomalies parametrization', level=3)  
            END IF
 
            anomaliesScaling = GetConstReal(Solver % Values, 'Anomalies scaling factor', Found)
            IF (.NOT. Found) THEN
                anomaliesScaling = 1.0
            END IF

            ar4_precipSol => VariableGet( Solver % Mesh % Variables, 'ar4_precipitation' )
            IF ( ASSOCIATED( ar4_precipSol ) ) THEN
                ar4_precip_Perm => ar4_precipSol % Perm
                ar4_precip      => ar4_precipSol % Values
                ar4_precip_DOFS = ar4_precipSol % DOFs
            ELSE
                CALL FATAL(SolverName, ' Could not find ar4_precipitation field variable')
            END IF

            ar4_annualtempSol => VariableGet( Solver % Mesh % Variables, 'ar4_annualtemp' )
            IF ( ASSOCIATED( ar4_annualtempSol ) ) THEN
                ar4_annualtemp_Perm   => ar4_annualtempSol % Perm
                ar4_annualtemp      => ar4_annualtempSol % Values
                ar4_annualtemp_DOFs   =  ar4_annualtempSol % DOFs
            ELSE
                CALL FATAL(SolverName, ' Could not find ar4_annualtemp field variable')
            END IF

            SELECT CASE(continent)
    
            CASE('greenland')


                ar4_julytempSol => VariableGet( Solver % Mesh % Variables, 'ar4_julytemp' )
                IF ( ASSOCIATED( ar4_julytempSol ) ) THEN
                    ar4_julytemp_Perm   => ar4_julytempSol % Perm
                    ar4_julytemp      => ar4_julytempSol % Values
                    ar4_julytemp_DOFs   =  ar4_julytempSol % DOFs
                ELSE
                    CALL FATAL(SolverName, ' Could not find ar4_julytemp field variable')
                END IF

            CASE('antarctica')

                ar4_julytempSol => VariableGet( Solver % Mesh % Variables, 'ar4_annualtemp' )
                IF ( ASSOCIATED( ar4_julytempSol ) ) THEN
                    ar4_julytemp_Perm   => ar4_julytempSol % Perm
                    ar4_julytemp      => ar4_julytempSol % Values
                    ar4_julytemp_DOFs   =  ar4_julytempSol % DOFs
                ELSE
                    CALL FATAL(SolverName, ' Could not find ar4_annualtemp field variable')
                END IF

            END SELECT


            watchdog = .FALSE.

            DO i=1,Solver % NumberOFActiveElements
   
                CurrentElement => GetActiveElement(i)
                NodeIndexes => CurrentElement % NodeIndexes
    
                DO k=1, GetElementNOFNodes(CurrentElement)

                    IF( .NOT. watchdog(NodeIndexes(k)) ) THEN
    
                        ar4_time1_precip(NodeIndexes(k)) = ar4_precip(ar4_precip_DOFS*(ar4_precip_Perm(NodeIndexes(k))-1)+1)
                        ar4_time1_annualtemp(NodeIndexes(k)) = & 
                            ar4_annualtemp(ar4_annualtemp_DOFs*(ar4_annualtemp_Perm(NodeIndexes(k))-1)+1)
                        ar4_time1_julytemp(NodeIndexes(k)) = ar4_julytemp(ar4_julytemp_DOFs*(ar4_julytemp_Perm(NodeIndexes(k))-1)+1)

                        watchdog(NodeIndexes(k)) = .TRUE.

                    END IF

                END DO

            END DO

            NULLIFY(ar4_precipSol, ar4_precip, ar4_precip_Perm)
            NULLIFY(ar4_annualtempSol, ar4_annualtemp, ar4_annualtemp_Perm)
            NULLIFY(ar4_julytempSol, ar4_julytemp, ar4_julytemp_Perm)

            WRITE(Message,'(a,F8.2)') 'Anomalies used with scaling factor: ', anomaliesScaling   
            CALL INFO(SolverName, Message, level=3)   
    
        END IF

        lapsrate = GetLogical(Solver % Values, 'Atmospheric laps rate', Found)
        IF (.NOT. Found) THEN
            CALL FATAL(SolverName, 'Parameter: Atmospheric laps rate not found')
        ELSE 
            IF (lapsrate) THEN  
                CALL INFO(SolverName, 'Atmospheric laps rate included in temperature parametrization', level=3)   
            ELSE
                CALL INFO(SolverName, 'Atmospheric laps rate not included in temperature parametrization', level=3)   
            END IF
        END IF

        convert = GetLogical(Solver % Values, 'Ice equivalent conversion', Found)
        IF (.NOT. Found) THEN
            CALL FATAL(SolverName, 'Parameter: Ice equivalent conversion not found')
        END IF 

        IF (TSURFACE == 4) THEN
            
           !-------- Read data for delta_ts --------
           INPATH = GetString(Solver % Values,'PATH',Found )
           IF (.NOT. Found) THEN
              CALL FATAL(SolverName,' PATH for accessing files not found.')
           END IF

           GRIP_TEMP_FILE = GetString(Solver % Values,'Grip Temperature file',Found )
           IF (.NOT. Found) THEN
              CALL FATAL(SolverName,' Name for grip temperature file not found.')
           END IF
                     
           griptemp = 0.0_dp 

           OPEN(21, iostat=ios, &
            file=GRIP_TEMP_FILE, &                   !INPATH//'/general/'//GRIP_TEMP_FILE, &
            status='old')
           IF (ios /= 0) &
           CALL FATAL(SolverName, ' Error when opening the data file for delta_ts!')

           READ(21,*) ch_dummy, grip_time_min, grip_time_stp, grip_time_max

           IF (ch_dummy /= '#') THEN
              CALL FATAL(SolverName, &
               ' grip_time_min, grip_time_stp, grip_time_max not defined in data file for delta_ts!')
           END IF

           ndata_grip = (grip_time_max-grip_time_min)/grip_time_stp

           IF (ndata_grip .GT. 10000) &
           CALL FATAL(SolverName, ' Too many data in data file for delta_ts!')

           DO n=0, ndata_grip
            READ(21,*) d_dummy, griptemp(n)
           END DO

           CLOSE(21, status='keep')
       
        END IF

        IF (TSURFACE == 5) THEN

           !-------- Read data for the glacial index --------
           GLAC_IND_FILE = GetString(Solver % Values,'Glacial-index file',Found )
           IF (.NOT. Found) THEN
              CALL FATAL(SolverName,' Name for glaciel-index file not found.')
           END IF

           OPEN(21, iostat=ios, &
            file=INPATH//'/general/'//GLAC_IND_FILE, status='old')
           IF (ios /= 0) CALL FATAL(SolverName, ' Error when opening the glacial-index file!')

           READ(21,*) ch_dummy, gi_time_min, gi_time_stp, gi_time_max

           IF (ch_dummy /= '#') THEN
              CALL FATAL(SolverName, & 
                 ' gi_time_min, gi_time_stp, gi_time_max not defined in glacial-index file!')
           END IF

           ndata_gi = (gi_time_max-gi_time_min)/gi_time_stp

           IF ( ALLOCATED(glacial_index) ) DEALLOCATE (glacial_index)
           ALLOCATE(glacial_index(0:ndata_gi))
 
           DO n=0, ndata_gi
             READ(21,*) d_dummy, glacial_index(n)
           END DO

           CLOSE(21, status='keep')

           !-------- Reading of LGM mean-annual and mean-July surface-temperature 
           !         anomalies --------

           TEMP_MA_ANOM_FILE = GetString(Solver % Values,'LLGM mean-annual anomalies',Found )
           IF (.NOT. Found) THEN
              CALL FATAL(SolverName,' Name for LLGM mean-annual anomalies file not found.')
           END IF

           TEMP_MJ_ANOM_FILE = GetString(Solver % Values,'LLGM mean-July anomalies',Found )
           IF (.NOT. Found) THEN
              CALL FATAL(SolverName,' Name for LLGM mean-July anomalies file not found.')
           END IF

           IF (GRID == 0 .OR. GRID == 1) THEN
             OPEN(21, iostat=ios, &
              file=INPATH//'/grl/'//TEMP_MA_ANOM_FILE, &
              recl=8192, status='old')
           END IF

           IF (ios /= 0) CALL FATAL(SolverName, ' Error when opening the temp_ma anomaly file!')

           IF (GRID == 0 .OR. GRID == 1) THEN
             OPEN(22, iostat=ios, &
              file=INPATH//'/grl/'//TEMP_MJ_ANOM_FILE, &
              recl=8192, status='old')
           END IF

           if (ios /= 0) CALL FATAL(SolverName, ' Error when opening the temp_mj anomaly file!')

           DO n=1, 6; READ(21,'(a)') ch_dummy; END DO
           Do n=1, 6; READ(22,'(a)') ch_dummy; END DO

           !The rest of the code which gets the data is not implemented yet as it needs to interpolate from the sicopilis grid to the finite-element mesh

           CLOSE(21, status='keep')
           CLOSE(22, status='keep')

        END IF
  
        sTime = 0.0_dp

        AllocationsDone = .TRUE.

  END IF

  TimeSol => VariableGet( Solver % Mesh % Variables, 'Timestep' )
  IF ( ASSOCIATED( TimeSol ) ) THEN
       TimePerm    => TimeSol % Perm
       Time        => TimeSol % Values
       TimeDOFs = TimeSol % DOFs
  ELSE
       CALL FATAL(SolverName, ' Could not find Time field variable')
  END IF

  zsSol => VariableGet( Solver % Mesh % Variables, 'FS' )
  If ( ASSOCIATED( zsSol ) ) THEN
       zsPerm    => zsSol % Perm
       zs        => zsSol % Values
       zsDOFs = zsSol % DOFs
  ELSE
       CALL FATAL(SolverName, ' Could not find FS field variable')
  END IF

  phiSol => VariableGet( Solver % Mesh % Variables, 'Latitude' )
  IF ( ASSOCIATED( phiSol ) ) THEN
       phiPerm    => phiSol % Perm
       phii       => phiSol % Values
       phiDOFs = phiSol % DOFs
  ELSE
       CALL FATAL(SolverName, ' Could not find Latitude field variable')
  END IF

  lambdaSol => VariableGet( Solver % Mesh % Variables, 'Longitude' )
  IF ( ASSOCIATED( lambdaSol ) ) THEN
       lambdaPerm    => lambdaSol % Perm
       lambda        => lambdaSol % Values
       lambdaDOFs = lambdaSol % DOFs
  ELSE
       CALL FATAL(SolverName, ' Could not find Longitude field variable')
  END IF

  zs_refSol => VariableGet( Solver % Mesh % Variables, 'ReferenceFS' )
  IF ( ASSOCIATED( zs_refSol ) ) THEN
       zs_refPerm    => zs_refSol % Perm
       zs_ref        => zs_refSol % Values
       zs_refDOFs = zs_refSol % DOFs
  ELSE
       CALL FATAL(SolverName, ' Could not find ReferenceFS field variable')
  END IF

  asSol => VariableGet( Solver % Mesh % Variables, 'as' )
  IF ( ASSOCIATED( asSol ) ) THEN
       asPerm    => asSol % Perm
       as        => asSol % Values
       asDOFs = asSol % DOFs
  ELSE
       CALL FATAL(SolverName, ' Could not find as field variable')
  END IF

  tempSol => VariableGet( Solver % Mesh % Variables, 'temp_surf' )
  IF ( ASSOCIATED( tempSol ) ) THEN
       tempPerm    => tempSol % Perm
       temp      => tempSol % Values
       tempDOFs = tempSol % DOFs
  ELSE
       CALL FATAL(SolverName, ' Could not find temp_surf field variable')
  END IF
       
  sTime = sTime + TimeStep
  !WRITE(*,*) 'Debug info: ', TimeStep, Time(1), sTime

  SELECT CASE(PrecipMode)
   CASE('constant')

    precip_ma_presSol => VariableGet( Solver % Mesh % Variables, 'precipitation' )
    IF ( ASSOCIATED( precip_ma_presSol ) ) THEN
           precip_ma_presPerm    => precip_ma_presSol % Perm
           precip_ma_present     => precip_ma_presSol % Values
           precip_ma_presDOFs = precip_ma_presSol % DOFs
    ELSE
           CALL FATAL(SolverName, ' Could not find precipitation field variable')
    END IF

    IF (Anomalies) THEN

        ar4_precipSol => VariableGet( Solver % Mesh % Variables, 'ar4_precipitation' )
            IF ( ASSOCIATED( ar4_precipSol ) ) THEN
                ar4_precip_Perm => ar4_precipSol % Perm
                ar4_precip      => ar4_precipSol % Values
                ar4_precip_DOFS = ar4_precipSol % DOFs
            ELSE
                CALL FATAL(SolverName, ' Could not find ar4_precipitation field variable')
            END IF

    END IF
 
    IF (ACCSURFACE <= 4) THEN
    !-------- Present monthly precipitation rates
    !         (assumed to be equal to the mean annual precipitation rate) --------
   
    DO n=1, 12

       watchdog = .FALSE.

        DO i=1,Solver % NumberOFActiveElements
   
            CurrentElement => GetActiveElement(i)
            NodeIndexes => CurrentElement % NodeIndexes

            Material => GetMaterial()
            IF (.NOT. ASSOCIATED(Material)) THEN
                WRITE(Message,'(A,I5,A)') 'No material found.'
                CALL FATAL(SolverName,Message)
            END IF

            RHO_W = GetConstReal( Material, 'Density of pure water', Found )
            IF (.NOT. Found) THEN
                CALL FATAL(SolverName,' Density of pure water not found', Found )
            END IF

            RHO = GetConstReal( Material, 'Density of ice', Found )
            IF (.NOT. Found) THEN
                CALL FATAL(SolverName,' Density of ice not found', Found )
            END IF
    
            DO k=1, GetElementNOFNodes(CurrentElement)

                IF( .NOT. watchdog(NodeIndexes(k)) ) THEN
      
                  precip_present(NodeIndexes(k),n) = &
                    (precip_ma_present(precip_ma_presDOFs*(precip_ma_presPerm(NodeIndexes(k))-1)+1)/YEAR_SEC)

                    ! m/a water equiv. -> m/a ice equiv.
                    IF(convert) THEN
                         precip_present(NodeIndexes(k),n) = precip_present(NodeIndexes(k),n) * (RHO_W/RHO)
                    END IF

                  watchdog(NodeIndexes(k)) = .TRUE.

                END IF   

            END DO
        END DO
    END DO
   
    NULLIFY(Material)
   
    END IF


   CASE('scenarios')

    precip_ma_presSol => VariableGet( Solver % Mesh % Variables, 'precipitation' )
    IF ( ASSOCIATED( precip_ma_presSol ) ) THEN
         precip_ma_presPerm    => precip_ma_presSol % Perm
         precip_ma_present     => precip_ma_presSol % Values
         precip_ma_presDOFs = precip_ma_presSol % DOFs
    ELSE
         CALL FATAL(SolverName, ' Could not find Precipitation field variable')
    END IF

    IF (ACCSURFACE <= 4) THEN
    !-------- Present monthly precipitation rates
    !         (assumed to be equal to the mean annual precipitation rate) --------
   
    DO n=1, 12

       watchdog = .FALSE.

       DO i=1,Solver % NumberOFActiveElements
   
          CurrentElement => GetActiveElement(i)
          NodeIndexes => CurrentElement % NodeIndexes

          Material => GetMaterial()
          IF (.NOT. ASSOCIATED(Material)) THEN
              WRITE(Message,'(A,I5,A)') 'No material found.'
              CALL FATAL(SolverName,Message)
          END IF

          RHO_W = GetConstReal( Material, 'Density of pure water', Found )
          IF (.NOT. Found) THEN
              CALL FATAL(SolverName,' Density of pure water not found', Found )
          END IF

          RHO = GetConstReal( Material, 'Density of ice', Found )
          IF (.NOT. Found) THEN
              CALL FATAL(SolverName,' Density of ice not found', Found )
          END IF
    
          DO k=1, GetElementNOFNodes(CurrentElement)

             IF( .NOT. watchdog(NodeIndexes(k)) ) THEN
      
                 precip_present(NodeIndexes(k),n) = &
                   (precip_ma_present(precip_ma_presDOFs*(precip_ma_presPerm(NodeIndexes(k))-1)+1)/YEAR_SEC)

                 ! m/a water equiv. -> m/a ice equiv.
                 IF(convert) THEN
                    precip_present(NodeIndexes(k),n) = precip_present(NodeIndexes(k),n) * (RHO_W/RHO)
                 END IF

                 watchdog(NodeIndexes(k)) = .TRUE.

             END IF   

          END DO
       END DO

   END DO

   NULLIFY(Material)
   
   END IF

   CASE DEFAULT
    
    CALL FATAL(SolverName,' Given Precipitation input mode not supported.', Found )   

  END SELECT

  delta_ts   = 0.0_dp
  glac_index = 0.0_dp

  !-------- Surface-temperature deviation from present values --------
  
  IF (TSURFACE == 1) THEN 
     delta_ts = DELTA_TS0

  ELSE IF (TSURFACE == 3) THEN
     delta_ts = SINE_AMPLIT*COS(2.0_dp*pi*Time(1)/(100000.0d0*YEAR_SEC))-SINE_AMPLIT 

  ELSE IF (TSURFACE == 4) THEN

    !  ------ delta_ts from ice-core record
   
    IF (Time(1)/YEAR_SEC .LT. REAL(grip_time_min,dp)) THEN
       delta_ts = griptemp(0)
    ELSE IF (Time(1)/YEAR_SEC .LT. REAL(grip_time_max,dp)) THEN

       i_kl = FLOOR(((Time(1)/YEAR_SEC) &
          -REAL(grip_time_min,dp))/REAL(grip_time_stp,dp))
       i_kl = MAX(i_kl, 0)

       i_gr = CEILING(((Time(1)/YEAR_SEC) &
          -REAL(grip_time_min,dp))/REAL(grip_time_stp,dp))
       i_gr = MIN(i_gr, ndata_grip)

       IF (i_kl .EQ. i_gr) THEN

          delta_ts = griptemp(i_kl)

       ELSE
 
          time_kl = (grip_time_min + i_kl*grip_time_stp) *YEAR_SEC
          time_gr = (grip_time_min + i_gr*grip_time_stp) *YEAR_SEC
          
          delta_ts = griptemp(i_kl) &
                  +(griptemp(i_gr)-griptemp(i_kl)) &
                  *(Time(1)-time_kl)/(time_gr-time_kl)
                  ! linear interpolation of the ice-core data

       END IF

    ELSE
       delta_ts  = griptemp(ndata_grip)
    END IF

    delta_ts = delta_ts * GRIP_TEMP_FACT
     !                ! modification by constant factor

  !-------- Glacial index --------

  ELSE IF (TSURFACE == 5) THEN

    IF (Time(1)/YEAR_SEC < REAL(gi_time_min,dp)) THEN
       glac_index = glacial_index(0)
    ELSE IF (Time(1)/YEAR_SEC < REAL(gi_time_max,dp)) THEN
       
       i_kl = FLOOR(((Time(1)/YEAR_SEC) &
          -REAL(gi_time_min,dp))/REAL(gi_time_stp,dp))
       i_kl = MAX(i_kl, 0)

       i_gr = CEILING(((Time(1)/YEAR_SEC) &
          -REAL(gi_time_min,dp))/REAL(gi_time_stp,dp))
       i_gr = MIN(i_gr, ndata_gi)

       IF (i_kl == i_gr) THEN

          glac_index = glacial_index(i_kl)

       ELSE
      
          time_kl = (gi_time_min + i_kl*gi_time_stp) *YEAR_SEC
          time_gr = (gi_time_min + i_gr*gi_time_stp) *YEAR_SEC

          glac_index = glacial_index(i_kl) &
                +(glacial_index(i_gr)-glacial_index(i_kl)) &
                *(Time(1)-time_kl)/(time_gr-time_kl)
                 ! linear interpolation of the glacial-index data

       END IF

    ELSE
       glac_index  = glacial_index(ndata_gi)
    END IF
  
  END IF

  SELECT CASE(TempParam)

   CASE('parametrization')

    SELECT CASE(continent)
    
     CASE('greenland')

      !-------- Surface air temperatures --------
      IF (TEMP_PRESENT_PARA == 1) THEN   ! Parameterization by Ritz et al. (1997)
       theta_ma = 49.13_dp
       gamma_ma = -7.992e-03_dp
       c_ma     = -0.7576_dp
       kappa_ma =  0.0_dp

       theta_mj = 30.38_dp
       gamma_mj = -6.277e-03_dp
       c_mj     = -0.3262_dp
       kappa_mj =  0.0_dp
      ELSE IF (TEMP_PRESENT_PARA == 2) THEN ! Parameterization by Fausto et al. (2009)
       theta_ma = 41.83_dp
       gamma_ma = -6.309e-03_dp
       c_ma     = -0.7189_dp
       kappa_ma = -0.0672_dp

       theta_mj = 14.70_dp
       gamma_mj = -5.426e-03_dp
       c_mj     = -0.1585_dp
       kappa_mj = -0.0518_dp
      ELSE 
       CALL FATAL(SolverName, ' Parameter TEMP_PRESENT_PARA must be either 1 or 2')
      END IF

     CASE('antarctica')

      !-------- Surface air temperatures --------
      IF (TEMP_PRESENT_PARA == 1) THEN    !  Parameterisation by Fortuin and Oerlemans
                                          !  (1990) for the whole ice sheet   
       theta_ma = 34.461_dp
       gamma_ma = -9.14e-03_dp
       c_ma     = -0.688_dp

       theta_mj = 16.81_dp
       gamma_mj = -6.92e-03_dp
       c_mj     = -0.27973_dp
      ELSE IF (TEMP_PRESENT_PARA == 2) THEN  ! Parameterisation by Fortuin and Oerlemans
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

       theta_mj   =  16.81_dp
       gamma_mj   =  -6.92e-03_dp
       c_mj       =  -0.27937_dp
      ELSE 
       CALL FATAL(SolverName, ' Parameter TEMP_PRESENT_PARA must be either 1 or 2')
      END IF

     CASE DEFAULT
      CALL FATAL(SolverName, ' Continent type not supported!!')

    END SELECT

    IF (Anomalies) THEN

        SELECT CASE(continent)
    
        CASE('greenland')
 

            ar4_annualtempSol => VariableGet( Solver % Mesh % Variables, 'ar4_annualtemp' )
            IF ( ASSOCIATED( ar4_annualtempSol ) ) THEN
                ar4_annualtemp_Perm   => ar4_annualtempSol % Perm
                ar4_annualtemp        => ar4_annualtempSol % Values
                ar4_annualtemp_DOFs   =  ar4_annualtempSol % DOFs
            ELSE
                CALL FATAL(SolverName, ' Could not find ar4_annualtemp field variable')
            END IF

            ar4_julytempSol => VariableGet( Solver % Mesh % Variables, 'ar4_julytemp' )
            IF ( ASSOCIATED( ar4_julytempSol ) ) THEN
                ar4_julytemp_Perm   => ar4_julytempSol % Perm
                ar4_julytemp        => ar4_julytempSol % Values
                ar4_julytemp_DOFs   =  ar4_julytempSol % DOFs
            ELSE
                CALL FATAL(SolverName, ' Could not find ar4_julytemp field variable')
            END IF

        CASE('antarctica')

            ar4_annualtempSol => VariableGet( Solver % Mesh % Variables, 'ar4_annualtemp' )
            IF ( ASSOCIATED( ar4_annualtempSol ) ) THEN
                ar4_annualtemp_Perm   => ar4_annualtempSol % Perm
                ar4_annualtemp        => ar4_annualtempSol % Values
                ar4_annualtemp_DOFs   =  ar4_annualtempSol % DOFs
            ELSE
                CALL FATAL(SolverName, ' Could not find ar4_annualtemp field variable')
            END IF

            ar4_julytempSol => VariableGet( Solver % Mesh % Variables, 'ar4_annualtemp' )
            IF ( ASSOCIATED( ar4_julytempSol ) ) THEN
                ar4_julytemp_Perm   => ar4_julytempSol % Perm
                ar4_julytemp        => ar4_julytempSol % Values
                ar4_julytemp_DOFs   =  ar4_julytempSol % DOFs
            ELSE
                CALL FATAL(SolverName, ' Could not find ar4_annualtemp field variable')
            END IF

        END SELECT

    END IF

    watchdog = .FALSE.

    DO i=1,Solver % NumberOFActiveElements
   
     CurrentElement => GetActiveElement(i)
     NodeIndexes => CurrentElement % NodeIndexes
    
     DO k=1, GetElementNOFNodes(CurrentElement)

       IF( .NOT. watchdog(NodeIndexes(k)) ) THEN

          IF( lapsrate ) THEN

              SELECT CASE(continent)
               CASE('greenland')

                !  ------ Present-day mean-annual air temperature
                temp_ma(NodeIndexes(k)) = theta_ma &
                  + gamma_ma*zs(zsDOFs*(zsPerm(NodeIndexes(k))-1)+1) &
                  + c_ma*phii(phiDOFs*(phiPerm(NodeIndexes(k))-1)+1)*pi_180_inv &
                  + kappa_ma*(modulo(lambda(lambdaDOFS*(lambdaPerm(NodeIndexes(k))-1)+1)+pi, 2.0_dp*pi)-pi) &
                    *pi_180_inv
                  ! west longitudes counted negatively
                    
                !  ------ Present-day mean-July air temperature
                temp_mm(NodeIndexes(k),7) = theta_mj &
                    + gamma_mj*zs(zsDOFs*(zsPerm(NodeIndexes(k))-1)+1) &
                    + c_mj*phii(phiDOFs*(phiPerm(NodeIndexes(k))-1)+1)*pi_180_inv &
                    + kappa_mj*(modulo(lambda(lambdaDOFS*(lambdaPerm(NodeIndexes(k))-1)+1)+pi, 2.0_dp*pi)-pi) &
                      *pi_180_inv
                    ! west longitudes counted negatively
                
               CASE('antarctica')
                !  ------ Present-day mean-annual air temperature
                IF (TEMP_PRESENT_PARA == 1) THEN
                   temp_ma(NodeIndexes(k)) =  theta_ma + gamma_ma*zs(zsDOFs*(zsPerm(NodeIndexes(k))-1)+1)  &
                                   + c_ma*abs(phii(phiDOFs*(phiPerm(NodeIndexes(k))-1)+1))*pi_180_inv
                ELSE IF (TEMP_PRESENT_PARA == 2) THEN 
                   IF ( zs(zsDOFs*(zsPerm(NodeIndexes(k))-1)+1) <= zs_sep_1 ) THEN
                       temp_ma(NodeIndexes(k)) =  theta_ma_1 + gamma_ma_1*zs(zsDOFs*(zsPerm(NodeIndexes(k))-1)+1)  &
                                        + c_ma_1*abs(phii(phiDOFs*(phiPerm(NodeIndexes(k))-1)+1))*pi_180_inv
                   ELSE IF ( zs(zsDOFs*(zsPerm(NodeIndexes(k))-1)+1) <= zs_sep_2 ) then
                       temp_ma(NodeIndexes(k)) =  theta_ma_2 + gamma_ma_2*zs(zsDOFs*(zsPerm(NodeIndexes(k))-1)+1)  &
                                        + c_ma_2*abs(phii(phiDOFs*(phiPerm(NodeIndexes(k))-1)+1))*pi_180_inv
                   ELSE
                       temp_ma(NodeIndexes(k)) = theta_ma_3 + gamma_ma_3*zs(zsDOFs*(zsPerm(NodeIndexes(k))-1)+1)  &
                                        + c_ma_3*abs(phii(phiDOFs*(phiPerm(NodeIndexes(k))-1)+1))*pi_180_inv
                   END IF
                END IF
                !  ------ Present-day mean-July air temperature
                temp_mm(NodeIndexes(k),7) =  theta_mj + gamma_mj*zs(zsDOFs*(zsPerm(NodeIndexes(k))-1)+1) &
                                   + c_mj*abs(phii(phiDOFs*(phiPerm(NodeIndexes(k))-1)+1))*pi_180_inv
              END SELECT

          ELSE 

             SELECT CASE(continent)
               CASE('greenland')

                !  ------ Present-day mean-annual air temperature
                temp_ma(NodeIndexes(k)) = theta_ma &
                  + gamma_ma*zs_ref(zs_refDOFs*(zs_refPerm(NodeIndexes(k))-1)+1) &
                  + c_ma*phii(phiDOFs*(phiPerm(NodeIndexes(k))-1)+1)*pi_180_inv &
                  + kappa_ma*(modulo(lambda(lambdaDOFS*(lambdaPerm(NodeIndexes(k))-1)+1)+pi, 2.0_dp*pi)-pi) &
                    *pi_180_inv
                  ! west longitudes counted negatively

                !  ------ Present-day mean-July air temperature
                temp_mm(NodeIndexes(k),7) = theta_mj &
                    + gamma_mj*zs_ref(zs_refDOFs*(zs_refPerm(NodeIndexes(k))-1)+1) &
                    + c_mj*phii(phiDOFs*(phiPerm(NodeIndexes(k))-1)+1)*pi_180_inv &
                    + kappa_mj*(modulo(lambda(lambdaDOFS*(lambdaPerm(NodeIndexes(k))-1)+1)+pi, 2.0_dp*pi)-pi) &
                      *pi_180_inv
                    ! west longitudes counted negatively

               CASE('antarctica')
                 !  ------ Present-day mean-annual air temperature
                IF (TEMP_PRESENT_PARA == 1) THEN
                   temp_ma(NodeIndexes(k)) =  theta_ma + gamma_ma*zs_ref(zs_refDOFs*(zs_refPerm(NodeIndexes(k))-1)+1) &
                                   + c_ma*abs(phii(phiDOFs*(phiPerm(NodeIndexes(k))-1)+1))*pi_180_inv
                ELSE IF (TEMP_PRESENT_PARA == 2) THEN 
                   IF ( zs_ref(zs_refDOFs*(zs_refPerm(NodeIndexes(k))-1)+1) <= zs_sep_1 ) THEN
                       temp_ma(NodeIndexes(k)) =  theta_ma_1 + gamma_ma_1*zs_ref(zs_refDOFs*(zs_refPerm(NodeIndexes(k))-1)+1) &
                                        + c_ma_1*abs(phii(phiDOFs*(phiPerm(NodeIndexes(k))-1)+1))*pi_180_inv
                   ELSE IF ( zs_ref(zs_refDOFs*(zs_refPerm(NodeIndexes(k))-1)+1)<= zs_sep_2 ) then
                       temp_ma(NodeIndexes(k)) =  theta_ma_2 + gamma_ma_2*zs_ref(zs_refDOFs*(zs_refPerm(NodeIndexes(k))-1)+1) &
                                        + c_ma_2*abs(phii(phiDOFs*(phiPerm(NodeIndexes(k))-1)+1))*pi_180_inv
                   ELSE
                       temp_ma(NodeIndexes(k)) = theta_ma_3 + gamma_ma_3*zs_ref(zs_refDOFs*(zs_refPerm(NodeIndexes(k))-1)+1) &
                                        + c_ma_3*abs(phii(phiDOFs*(phiPerm(NodeIndexes(k))-1)+1))*pi_180_inv
                   END IF
                END IF
                !  ------ Present-day mean-July air temperature
                temp_mm(NodeIndexes(k),7) =  theta_mj + gamma_mj*zs_ref(zs_refDOFs*(zs_refPerm(NodeIndexes(k))-1)+1) &
                                   + c_mj*abs(phii(phiDOFs*(phiPerm(NodeIndexes(k))-1)+1))*pi_180_inv
               END SELECT

          END IF
     
          IF (TSURFACE <= 4) THEN
             !  ------ Correction with deviation delta_ts
             temp_ma(NodeIndexes(k)) = temp_ma(NodeIndexes(k)) + delta_ts
             temp_mm(NodeIndexes(k),7) = temp_mm(NodeIndexes(k),7) + delta_ts

          ELSE IF (TSURFACE == 5) THEN
             !  ------ Correction with LGM anomaly and glacial index
             !WARNING: The program does not get the data for the anomalies correctly, so using TSURFACE == 5 won't work
             temp_ma(NodeIndexes(k)) = temp_ma(NodeIndexes(k)) + glac_index*temp_ma_lgm_anom(NodeIndexes(k))
             temp_mm(NodeIndexes(k),7) = temp_mm(NodeIndexes(k),7) + glac_index*temp_mj_lgm_anom(NodeIndexes(k))
 
          END IF   

          ! ------ Add anomalies if required
          IF (Anomalies) THEN 

            anomalies_annualtemp(NodeIndexes(k)) = &
            ar4_annualtemp(ar4_annualtemp_DOFs*(ar4_annualtemp_Perm(NodeIndexes(k))-1)+1) - ar4_time1_annualtemp(NodeIndexes(k))

            anomalies_julytemp(NodeIndexes(k)) = &
            ar4_julytemp(ar4_julytemp_DOFs*(ar4_julytemp_Perm(NodeIndexes(k))-1)+1) - ar4_time1_julytemp(NodeIndexes(k))

            IF (anomaliesParam) THEN
                
                IF (TimeStep < 1.0) THEN
                    IF ( sTime < (95.0/YEAR_SEC)) THEN

                        ! Remember the anomalies, the last time is before t = 95 a
                        anamaliesAnnualTempAtT94(NodeIndexes(k)) = anomalies_annualtemp(NodeIndexes(k)) 
                        anamaliesJulyTempAtT94(NodeIndexes(k))  =  anomalies_julytemp(NodeIndexes(k)) 

                        anomalies_annualtemp(NodeIndexes(k)) = anomalies_annualtemp(NodeIndexes(k)) * 1.5
                        anomalies_julytemp(NodeIndexes(k)) =  anomalies_julytemp(NodeIndexes(k)) * 1.5 
                   
                    ELSE IF (sTime >= (95.0/YEAR_SEC) .AND. sTime < (200.0/YEAR_SEC)) THEN

                        anomalies_annualtemp(NodeIndexes(k)) = &
                            1.5 *  ( anamaliesAnnualTempAtT94(NodeIndexes(k)) / (94.0/YEAR_SEC) ) * sTime;
                        anomalies_julytemp(NodeIndexes(k)) = & 
                            1.5 *  ( anamaliesJulyTempAtT94(NodeIndexes(k)) / (94.0/YEAR_SEC) ) * sTime;

                    ELSE IF (sTime >=  (200.0/YEAR_SEC) .AND. sTime <= (500.0/YEAR_SEC)) THEN

                        anomalies_annualtemp(NodeIndexes(k)) = &
                        1.5 * ( (200.0/YEAR_SEC) * anamaliesAnnualTempAtT94(NodeIndexes(k)) / (94.0/YEAR_SEC) ) &
                        + 0.07*( anamaliesAnnualTempAtT94(NodeIndexes(k)) / (94.0/YEAR_SEC) ) * (stime-(200.0/YEAR_SEC)) 

                        anomalies_julytemp(NodeIndexes(k)) = & 
                        1.5 * ( (200.0/YEAR_SEC) * anamaliesJulyTempAtT94(NodeIndexes(k)) / (94.0/YEAR_SEC) ) &
                        + 0.07*( anamaliesJulyTempAtT94(NodeIndexes(k)) / (94.0/YEAR_SEC) ) * (sTime-(200.0/YEAR_SEC)) 

                    END IF

                ELSE IF (TimeStep >= 1) THEN
                    CALL FATAL(SolverName, ' Time step greater than or equal to 1 not yet supported in anomalies param.')
                END IF
            
            ElSE

                 anomalies_annualtemp(NodeIndexes(k)) = anomalies_annualtemp(NodeIndexes(k))*anomaliesScaling 
                 anomalies_julytemp(NodeIndexes(k)) =  anomalies_julytemp(NodeIndexes(k))*anomaliesScaling 
           
            END IF

            temp_ma(NodeIndexes(k)) = temp_ma(NodeIndexes(k)) + anomalies_annualtemp(NodeIndexes(k))
            temp_mm(NodeIndexes(k),7) = temp_mm(NodeIndexes(k),7) + anomalies_julytemp(NodeIndexes(k))

          END IF

          !------ Amplitude of the annual cycle
          temp_ampl(NodeIndexes(k)) = temp_mm(NodeIndexes(k),7) - temp_ma(NodeIndexes(k))

          IF (temp_ampl(NodeIndexes(k)) < eps ) THEN
             temp_ampl(NodeIndexes(k)) = eps  ! Correction of amplitude, if required
          END IF

          watchdog(NodeIndexes(k)) = .TRUE.

       END IF

     END DO ! End number of element nodes
    END DO ! End number of elements

   CASE('scenarios')

    annualtempSol => VariableGet( Solver % Mesh % Variables, 'annualtemp' )
    IF ( ASSOCIATED( annualtempSol ) ) THEN
       annualtempPerm    => annualtempSol % Perm
       annualtemp        => annualtempSol % Values
       annualtempDOFs = annualtempSol % DOFs
    ELSE
       CALL FATAL(SolverName, ' Could not find annualtemp field variable')
    END IF

    julytempSol => VariableGet( Solver % Mesh % Variables, 'julytemp' )
    IF ( ASSOCIATED( julytempSol ) ) THEN
       julytempPerm    => julytempSol % Perm
       julytemp        => julytempSol % Values
       julytempDOFs = julytempSol % DOFs
    ELSE
       CALL FATAL(SolverName, ' Could not find julytemp field variable')
    END IF

    watchdog = .FALSE.

    DO i=1,Solver % NumberOFActiveElements
   
     CurrentElement => GetActiveElement(i)
     NodeIndexes => CurrentElement % NodeIndexes
    
     DO k=1, GetElementNOFNodes(CurrentElement)

       IF( .NOT. watchdog(NodeIndexes(k)) ) THEN

       temp_ma(NodeIndexes(k)) = annualtemp(annualTempDOFs*(annualTempPerm(NodeIndexes(k))-1)+1)

       temp_mm(NodeIndexes(k),7) = julytemp(julyTempDOFs*(julyTempPerm(NodeIndexes(k))-1)+1)

       IF (TSURFACE <= 4) THEN
       !  ------ Correction with deviation delta_ts
       temp_ma(NodeIndexes(k)) = temp_ma(NodeIndexes(k)) + delta_ts
       temp_mm(NodeIndexes(k),7) = temp_mm(NodeIndexes(k),7) + delta_ts

       ELSE IF (TSURFACE == 5) THEN
       !  ------ Correction with LGM anomaly and glacial index
       !WARNING: The program does not get the data for the anomalies correctly, so using TSURFACE == 5 won't work
       temp_ma(NodeIndexes(k)) = temp_ma(NodeIndexes(k)) + glac_index*temp_ma_lgm_anom(NodeIndexes(k))
       temp_mm(NodeIndexes(k),7) = temp_mm(NodeIndexes(k),7) + glac_index*temp_mj_lgm_anom(NodeIndexes(k))
     
       END IF
  
       !------ Amplitude of the annual cycle
       temp_ampl(NodeIndexes(k)) = temp_mm(NodeIndexes(k),7) - temp_ma(NodeIndexes(k))

       IF (temp_ampl(NodeIndexes(k)) < eps ) THEN
        temp_ampl(NodeIndexes(k)) = eps  ! Correction of amplitude, if required
       END IF

       watchdog(NodeIndexes(k)) = .TRUE.

       END IF

     END DO ! End number of element nodes
    END DO ! End number of elements

   CASE DEFAULT
    
    CALL FATAL(SolverName,' Given temperature input mode not supported.', Found )

  END SELECT


  !  ------ Monthly temperatures
  DO n=1, 12   ! month counter

     sine_factor = sin((real(n,dp)-4.0_dp)*pi/6.0_dp)

     DO i=1,Solver % NumberOFActiveElements
      
        CurrentElement => GetActiveElement(i)
        NodeIndexes => CurrentElement % NodeIndexes
        
         DO k=1, GetElementNOFNodes(CurrentElement)

           temp_mm(NodeIndexes(k),n) = temp_ma(NodeIndexes(k)) + sine_factor*temp_ampl(NodeIndexes(k))

         END DO
     END DO

  END DO 

  watchdog = .FALSE.

  DO i=1,Solver % NumberOFActiveElements

   CurrentElement => GetActiveElement(i)
   NodeIndexes => CurrentElement % NodeIndexes

   body_id = CurrentElement % BodyId

   IF (body_id /= old_body) Then 
     
      old_body = body_id

      !--------Get Material parameters--------
 
      Material => GetMaterial()
      IF (.NOT. ASSOCIATED(Material)) THEN
           WRITE(Message,'(A,I5,A)') 'No material found.'
           CALL FATAL(SolverName,Message)
      END IF

      S_STAT_0 = GetConstReal( Material, 'Standard deviation of the air temperature', Found )
      IF (.NOT. Found) THEN
         CALL FATAL(SolverName,' Standard deviation of the air temperature not found.')
      END IF

      PHI_SEP_0 = GetConstReal( Material, 'Separation latitude', Found )
      IF (.NOT. Found) THEN
         CALL FATAL(SolverName,' Separation latitude not found.')
      END IF

      PMAX_0 = GetConstReal( Material, 'Saturation factor', Found )
      IF (.NOT. Found) THEN
         CALL FATAL(SolverName,' Separation factor not found.')
      END IF

      MU_0 = GetConstReal( Material, 'Firn-warming correction', Found )
      IF (.NOT. Found) THEN
         CALL FATAL(SolverName,' Firn-warming correction not found.')
      END IF

      RHO_W = GetConstReal( Material, 'Density of pure water', Found )
      IF (.NOT. Found) THEN
         CALL FATAL(SolverName,' Density of pure water not found', Found )
      END IF

      RHO = GetConstReal( Material, 'Density of ice', Found )
      IF (.NOT. Found) THEN
         CALL FATAL(SolverName,' Density of ice not found', Found )
      END IF

     TarasovPDD = GetLogical(Material, 'Tarasav PDD factors', Found)
     HuybrechtsPDD = GetLogical(Material, 'Huybrechts PDD factors', Found)

     IF (TarasovPDD) THEN 

        BETA1_LT_0 = &
          GetConstReal( Material, 'Degree-day factor for snow at low summer temperature', Found )
        IF (.NOT. Found) THEN
         CALL FATAL(SolverName,' Degree-day factor for snow at low summer temperature not found', Found )
        END IF

        BETA1_HT_0 = &
          GetConstReal( Material, 'Degree-day factor for snow at high summer temperature', Found )
        IF (.NOT. Found) THEN
         CALL FATAL(SolverName,' Degree-day factor for snow at high summer temperature not found', Found )
        END IF

        BETA2_LT_0 = &
          GetConstReal( Material, 'Degree-day factor for ice at low summer temperature', Found )
        IF (.NOT. Found) THEN
         CALL FATAL(SolverName,' Degree-day factor for ice at low summer temperature not found', Found )
        END IF

        BETA2_HT_0 = &
          GetConstReal( Material, 'Degree-day factor for ice at high summer temperature', Found )
        IF (.NOT. Found) THEN
         CALL FATAL(SolverName,' Degree-day factor for ice at high summer temperature not found', Found )
        END IF

    ELSE IF (HuybrechtsPDD) THEN

        BETA_SNOW_0 =  GetConstReal( Material, 'Degree-day factor for snow', Found )
        IF (.NOT. Found) THEN
         CALL FATAL(SolverName,' Degree-day factor for snow not found', Found )
        END IF

        BETA_ICE_0 =  GetConstReal( Material, 'Degree-day factor for ice', Found )
        IF (.NOT. Found) THEN
         CALL FATAL(SolverName,' Degree-day factor for ice not found', Found )
        END IF

    ELSE 
    
        CALL FATAL(SolverName,'Method for PDD factors not specified', Found )

    END IF
        
    NULLIFY(Material)

   END IF

   !-------- Accumulation-ablation function as_perp --------
   gamma_p   = -0.001_dp*log(2.0_dp)     ! Precipitation lapse rate,
                                        ! in m^(-1)
   zs_thresh = 2000.0_dp                ! Threshold for elevation desertification, in m
  
   IF (SOLID_PRECIP == 1) THEN           ! Marsiat (1994)
      temp_rain =    7.0_dp             ! Threshold monthly mean temperature for
                                        ! precipitation = 100% rain, in deg C
      temp_snow =  -10.0_dp             ! Threshold monthly mean temperature for &
                                        ! precipitation = 100% snow, in deg C
      inv_delta_temp_rain_snow = 1.0_dp/(temp_rain-temp_snow)

   ELSE IF (SOLID_PRECIP == 2) THEN       ! Bales et al. (2009)
      temp_rain =    7.2_dp             ! Threshold monthly mean temperature for &
                                        ! precipitation = 100% rain, in deg C
      temp_snow =  -11.6_dp             ! Threshold monthly mean temperature for &
                                        ! precipitation = 100% snow, in deg C

      coeff(0) =  5.4714e-01_dp         ! Coefficients
      coeff(1) = -9.1603e-02_dp         ! of
      coeff(2) = -3.314e-03_dp          ! the
      coeff(3) =  4.66e-04_dp           ! fifth-order
      coeff(4) =  3.8e-05_dp            ! polynomial
      coeff(5) =  6.0e-07_dp            ! fit

   ELSE IF (SOLID_PRECIP == 3) THEN       ! Huybrechts and de Wolde (1999)
      temp_rain = 2.0_dp                ! Threshold instantaneous temperature for &
                                        ! precipitation = 100% rain, in deg C
      temp_snow = temp_rain             ! Threshold instantaneous temperature for &
                                        ! precipitation = 100% snow, in deg C

      s_stat    = S_STAT_0              ! Standard deviation of the air termperature
                                        ! (same parameter as in the PDD model)

      inv_sqrt2_s_stat = 1.0_dp/(sqrt(2.0_dp)*s_stat)
  
   END IF

   IF (ABLSURFACE == 1 .OR. ABLSURFACE == 2) THEN
    
      s_stat   = S_STAT_0
     
      phi_sep  = PHI_SEP_0*pi_180       ! separates different domains for computation of
                                        ! degree-day factors beta1 and beta2
                                        ! deg N --> rad

      temp_lt  = -1.0_dp
      temp_ht  = 10.0_dp
      inv_delta_temp_ht_lt = 1.0_dp/(temp_ht-temp_lt)

      IF (TarasovPDD) THEN

        beta1_lt = BETA1_LT_0  *(0.001_dp*365.0_dp)*(RHO_W/RHO)
                           ! (mm WE)/(d*deg C) --> (m IE)/(a*deg C)
        beta1_ht = BETA1_HT_0  *(0.001_dp*365.0_dp)*(RHO_W/RHO)
                           ! (mm WE)/(d*deg C) --> (m IE)/(a*deg C)
        beta2_lt = BETA2_LT_0  *(0.001_dp*365.0_dp)*(RHO_W/RHO)
                           ! (mm WE)/(d*deg C) --> (m IE)/(a*deg C)
        beta2_ht = BETA2_HT_0  *(0.001_dp*365.0_dp)*(RHO_W/RHO)
                           ! (mm WE)/(d*deg C) --> (m IE)/(a*deg C)

      ELSE IF (HuybrechtsPDD) THEN

        beta_snow = BETA_SNOW_0 *(0.001_dp*365.0_dp)*(RHO_W/RHO)
                           ! (mm WE)/(d*deg C) --> (m IE)/(a*deg C)
        beta_ice =  BETA_ICE_0 *(0.001_dp*365.0_dp)*(RHO_W/RHO)
                           ! (mm WE)/(d*deg C) --> (m IE)/(a*deg C)

      END IF

      Pmax     = PMAX_0
      mu       = MU_0        *(1000.0_dp/365.0_dp)*(RHO/RHO_W)
                           ! (d*deg C)/(mm WE) --> (a*deg C)/(m IE)

   ELSE IF (ABLSURFACE == 3) THEN
     
      lambda_lti = LAMBDA_LTI_C *(0.001_dp/YEAR_SEC)*(RHO_W/RHO)
                                  ! (mm WE)/(a*deg C) --> depending on YEAR_SEC value (m IE)/(s*deg C) or (m IE)/(a*deg C)
      temp_lti   = TEMP_LTI_C
      mu         = 0.0_dp       ! no superimposed ice considered

   END IF
   
   DO k=1,GetElementNOFNodes(CurrentElement)

    IF( .NOT. watchdog(NodeIndexes(k)) ) THEN

     IF (ABLSURFACE == 1 .OR. ABLSURFACE == 2) THEN
      
      IF (TarasovPDD) THEN

        IF (phii(phiDOFs*(phiPerm(NodeIndexes(k))-1)+1) <= phi_sep) THEN
          beta1 = beta1_ht
          beta2 = beta2_ht
        ELSE
          IF (temp_mm(NodeIndexes(k),7) >= temp_ht) THEN
             beta1 = beta1_ht
             beta2 = beta2_ht
          ELSE IF (temp_mm(NodeIndexes(k),7) <= temp_lt) THEN
             beta1 = beta1_lt
             beta2 = beta2_lt
          ELSE
             beta1 = beta1_lt &
                 + (beta1_ht-beta1_lt) &
                   *inv_delta_temp_ht_lt*(temp_mm(NodeIndexes(k),7)-temp_lt)
             beta2 = beta2_ht &
                 + (beta2_lt-beta2_ht) &
                   *(inv_delta_temp_ht_lt*(temp_ht-temp_mm(NodeIndexes(k),7)))**3
          END IF
        END IF
 
      ELSE IF(HuybrechtsPDD) THEN

        beta1 = beta_snow
        beta2 = beta_ice

      END IF

     END IF

     !  ------ Accumulation
 
     IF (ACCSURFACE <= 4) THEN

      !    ---- Elevation desertification of precipitation

      IF (zs_ref(zs_refDOFs*(zs_refPerm(NodeIndexes(k))-1)+1) < zs_thresh) THEN
         precip_fact = exp(gamma_p*(max(zs(zsDOFs*(zsPerm(NodeIndexes(k))-1)+1),zs_thresh)-zs_thresh))
      ELSE 
         precip_fact = exp(gamma_p*(max(zs(zsDOFs*(zsPerm(NodeIndexes(k))-1)+1),zs_thresh) &
                     -zs_ref(zs_refDOFs*(zs_refPerm(NodeIndexes(k))-1)+1)))
      END IF

      DO n=1, 12   ! month counter
         precip(NodeIndexes(k),n) = precip_present(NodeIndexes(k),n)*precip_fact
      END DO
     END IF

     !    ---- Precipitation change related to changing climate
   
     IF (ACCSURFACE == 1) THEN
      precip_fact = ACCFACT
     ELSE IF (ACCSURFACE == 2) THEN
      precip_fact = 1.0_dp + GAMMA_S*delta_ts
     ELSE IF (ACCSURFACE == 3) THEN
      precip_fact = exp(GAMMA_S*delta_ts)
     END IF

     IF (ACCSURFACE <= 4) THEN

      precip(NodeIndexes(k),0) = 0.0_dp   ! initialisation value for mean annual precip

      DO n=1, 12   ! month counter
         precip(NodeIndexes(k),n) = precip(NodeIndexes(k),n)*precip_fact    ! monthly precip
         precip(NodeIndexes(k),0) = precip(NodeIndexes(k),0) + precip(NodeIndexes(k),n)*inv_twelve
                                                                            ! mean annual precip
      END DO
   
     ELSE IF (ACCSURFACE == 5) THEN
     ! WARNING: data for the precipitation anomalies are not correctly read by the program as it needs interpolation from sicopilis grid to finite-element
     ! mesh. The option ACCSURFACE == 5 won't work.
 
       precip(NodeIndexes(k),0) = 0.0_dp   ! initialisation value for mean annual precip

       DO n=1, 12   ! month counter
        IF (PRECIP_ANOM_INTERPOL == 1) THEN
            precip_fact = 1.0_dp-glac_index+glac_index*precip_lgm_anom(NodeIndexes(k),n)
                    ! interpolation with a linear function
        ELSE IF (PRECIP_ANOM_INTERPOL == 2) THEN 
            precip_fact = exp(-glac_index*gamma_precip_lgm_anom(NodeIndexes(k),n))
                    ! interpolation with an exponential function
        END IF

        precip(NodeIndexes(k),n) = precip_present(NodeIndexes(k),n)*precip_fact  ! monthly precip
        precip(NodeIndexes(k),0) = precip(NodeIndexes(k),0) + precip(NodeIndexes(k),n)*inv_twelve
                                                                                 ! mean annual precip
       END DO
 
     ELSE IF (ACCSURFACE == 6) THEN 
 
       precip(NodeIndexes(k),0) = 0.0_dp   ! initialisation value for mean annual precip

       !      -- Mean-annual precipitation already read from data (above) --------

       IF (Anomalies) THEN
 
            anomalies_precip(NodeIndexes(k)) = &
            ar4_precip(ar4_precip_DOFS*(ar4_precip_Perm(NodeIndexes(k))-1)+1) - ar4_time1_precip(NodeIndexes(k))
  
            IF (anomaliesParam) THEN

                IF (TimeStep < 1.0) THEN 
                    IF (sTime < (95.0/YEAR_SEC))THEN
                    
                        ! Remember the anomalies, the last time is before t = 95 a
                        anomaliesPrecipAtT94(NodeIndexes(k)) = anomalies_precip(NodeIndexes(k)) 

                        anomalies_precip(NodeIndexes(k)) = anomalies_precip(NodeIndexes(k)) * 1.5
                   
                    ELSE IF (sTime >= (95.0/YEAR_SEC) .AND. sTime < (200.0/YEAR_SEC)) THEN

                        anomalies_precip(NodeIndexes(k)) = &
                            1.5 *  ( anomaliesPrecipAtT94(NodeIndexes(k)) / (94.0/YEAR_SEC) ) * sTime;

                    ELSE IF (sTime >=  (200.0/YEAR_SEC) .AND. sTime <= (500.0/YEAR_SEC)) THEN

                        anomalies_precip(NodeIndexes(k)) = &
                            1.5 * ( (200.0/YEAR_SEC) * anomaliesPrecipAtT94(NodeIndexes(k)) / (94.0/YEAR_SEC) ) + &
                            0.07*( anomaliesPrecipAtT94(NodeIndexes(k)) / (94.0/YEAR_SEC) ) * (sTime-(200.0/YEAR_SEC)) 

                    END IF

                ELSE IF (TimeStep >= 1) THEN
                    CALL FATAL(SolverName, ' Time step greater than or equal to 1 not yet supported in anomalies param.')
                END IF
            
            ElSE
                
                anomalies_precip(NodeIndexes(k)) = anomalies_precip(NodeIndexes(k)) * anomaliesScaling
           
            END IF

            IF(.NOT. convert) THEN
                anomalies_precip(NodeIndexes(k)) =  anomalies_precip(NodeIndexes(k)) * (RHO_W/RHO) 
            END IF
        
            precip(NodeIndexes(k),0) = ( precip_ma_present(precip_ma_presDOFs*(precip_ma_presPerm(NodeIndexes(k))-1)+1) + &
                anomalies_precip(NodeIndexes(k)) )/YEAR_SEC ! mean annual precip

       ELSE
            precip(NodeIndexes(k),0) = precip_ma_present(precip_ma_presDOFs*(precip_ma_presPerm(NodeIndexes(k))-1)+1)/YEAR_SEC 
                 ! mean annual precip
   
       END IF

       IF(convert) THEN
          precip(NodeIndexes(k),0) = precip(NodeIndexes(k),0) * (RHO_W/RHO) 
       END IF

       !    ---- Make sure no negative precip

       IF (  precip(NodeIndexes(k),0) < -eps) THEN
        
            precip(NodeIndexes(k),0) = 0.0_dp

       END IF 

       DO n=1, 12   ! month counter
  
          precip(NodeIndexes(k),n) = precip(NodeIndexes(k),0)   ! monthly precip, assumed to be equal
                                                                ! to the mean annual precip
       END DO

     END IF

          
     !    ---- Annual accumulation, snowfall and rainfall rates

     accum(NodeIndexes(k)) = precip(NodeIndexes(k),0)

     snowfall(NodeIndexes(k)) = 0.0_dp   ! initialisation value

     DO n=1, 12   ! month counter   
      
      IF (SOLID_PRECIP == 1) THEN  ! Marsiat (1994) 
         IF (temp_mm(NodeIndexes(k),n) >= temp_rain) THEN
            frac_solid = 0.0_dp
         ELSE IF (temp_mm(NodeIndexes(k),n) <= temp_snow) THEN
            frac_solid = 1.0_dp
         ELSE
            frac_solid = (temp_rain-temp_mm(NodeIndexes(k),n))*inv_delta_temp_rain_snow
         END IF
      ELSE IF (SOLID_PRECIP == 2) THEN ! Bales et al. (2009) 
         IF (temp_mm(NodeIndexes(k),n) >= temp_rain) THEN
            frac_solid = 0.0_dp
         ELSE IF (temp_mm(NodeIndexes(k),n) <= temp_snow) THEN
            frac_solid = 1.0_dp
         ELSE
            frac_solid = coeff(0) + temp_mm(NodeIndexes(k),n) * ( coeff(1) &
                                  + temp_mm(NodeIndexes(k),n) * ( coeff(2) &
                                  + temp_mm(NodeIndexes(k),n) * ( coeff(3) &
                                  + temp_mm(NodeIndexes(k),n) * ( coeff(4) &
                                  + temp_mm(NodeIndexes(k),n) *   coeff(5) ) ) ) )
               ! evaluation of 5th-order polynomial by Horner scheme
          END IF
       ELSE IF (SOLID_PRECIP == 3) THEN ! Huybrechts and de Wolde (1999)

          frac_solid = 1.0_dp &
                   - 0.5_dp*erfcc((temp_rain-temp_mm(NodeIndexes(k),n))*inv_sqrt2_s_stat)

       END IF

       snowfall(NodeIndexes(k)) = snowfall(NodeIndexes(k)) + precip(NodeIndexes(k),n)*frac_solid*inv_twelve
    
     END DO

     rainfall(NodeIndexes(k)) = precip(NodeIndexes(k),0) - snowfall(NodeIndexes(k))
    
     IF (snowfall(NodeIndexes(k)) < -eps) CALL FATAL(SolverName, ' boundary: Negative snowfall rate')
     IF (rainfall(NodeIndexes(k)) < -eps) CALL FATAL(SolverName, ' boundary: Negative rainfall rate!')

     !  ------ Ablation

     !    ---- Runoff
    
     IF (ABLSURFACE == 1 .OR. ABLSURFACE == 2) THEN
 
       !      -- Temperature excess ET

       Do n=1, 12   ! month counter
         temp_mm_help(n) = temp_mm(NodeIndexes(k),n)
       end do

       call pdd_f(temp_mm_help, s_stat, ET(NodeIndexes(k)))

       !      -- Formation rate of superimposed ice (melt_star), melt rate (melt)
       !         and runoff rate (runoff)

       IF (ABLSURFACE == 1) THEN

          IF ((beta1*ET(NodeIndexes(k))) <= (Pmax*snowfall(NodeIndexes(k)))) THEN
             melt_star(NodeIndexes(k)) = beta1*ET(NodeIndexes(k))
             melt(NodeIndexes(k))      = 0.0_dp
             runoff(NodeIndexes(k))    = melt(NodeIndexes(k))+rainfall(NodeIndexes(k))
          ELSE
             melt_star(NodeIndexes(k)) = Pmax*snowfall(NodeIndexes(k))
             melt(NodeIndexes(k))      = beta2*(ET(NodeIndexes(k))-melt_star(NodeIndexes(k))/beta1)
             runoff(NodeIndexes(k))    = melt(NodeIndexes(k))+rainfall(NodeIndexes(k))
          END IF
     
       ELSE IF (ABLSURFACE == 2) THEN
             IF ( rainfall(NodeIndexes(k)) <= (Pmax*snowfall(NodeIndexes(k))) ) THEN
                IF ( (rainfall(NodeIndexes(k))+beta1*ET(NodeIndexes(k))) &
                                        <= (Pmax*snowfall(NodeIndexes(k))) ) THEN
                   melt_star(NodeIndexes(k)) = rainfall(NodeIndexes(k))+beta1*ET(NodeIndexes(k))
                   melt(NodeIndexes(k))      = 0.0_dp
                   runoff(NodeIndexes(k))    = melt(NodeIndexes(k))
                ELSE
                   melt_star(NodeIndexes(k)) = Pmax*snowfall(NodeIndexes(k))
                   melt(NodeIndexes(k))      = beta2 &
                          *(ET(NodeIndexes(k))-(melt_star(NodeIndexes(k))-rainfall(NodeIndexes(k)))/beta1)
                   runoff(NodeIndexes(k))    = melt(NodeIndexes(k))
                END IF
             ELSE
                melt_star(NodeIndexes(k)) = Pmax*snowfall(NodeIndexes(k))
                melt(NodeIndexes(k))      = beta2*ET(NodeIndexes(k))
                runoff(NodeIndexes(k))    = melt(NodeIndexes(k)) + & 
                                 rainfall(NodeIndexes(k))-Pmax*snowfall(NodeIndexes(k))
             END IF
        END IF
     ELSE IF (ABLSURFACE == 3) THEN

        temp_jja_help  = & 
              one_third*(temp_mm(NodeIndexes(k),6)+temp_mm(NodeIndexes(k),7)+temp_mm(NodeIndexes(k),8))
        melt_star(NodeIndexes(k)) = 0.0_dp   ! no superimposed ice considered
        melt(NodeIndexes(k))      = lambda_lti*max((temp_jja_help-temp_lti), 0.0_dp)
        runoff(NodeIndexes(k))    = melt(NodeIndexes(k)) + rainfall(NodeIndexes(k))
     END IF

     !    ---- Evaporation

     evap(NodeIndexes(k)) = 0.0_dp

     !  ------ Accumulation minus ablation

     as_perp(NodeIndexes(k)) = accum(NodeIndexes(k)) - evap(NodeIndexes(k)) - runoff(NodeIndexes(k))

     !  ------ Ice-surface temperature (10-m firn temperature) temp_s,
     !         including empirical firn-warming correction due to
     !         refreezing meltwater when superimposed ice is formed
 
     IF (melt_star(NodeIndexes(k)) .GE. melt(NodeIndexes(k))) THEN
      temp_s(NodeIndexes(k)) = temp_ma(NodeIndexes(k)) &
                    +mu*(melt_star(NodeIndexes(k))-melt(NodeIndexes(k)))
     ELSE
      temp_s(NodeIndexes(k)) = temp_ma(NodeIndexes(k))
     END IF

     IF (temp_s(NodeIndexes(k)) > -0.001_dp) temp_s(NodeIndexes(k)) = -0.001_dp
                               ! Cut-off of positive air temperatures
     
     watchdog(NodeIndexes(k)) = .TRUE.
    
    END IF
 
   END DO ! Number of element nodes

  END DO ! Number of elements

  DO i=1,Solver % NumberOFActiveElements

    CurrentElement => GetActiveElement(i)
    NodeIndexes => CurrentElement % NodeIndexes

    DO k=1,GetElementNOFNodes(CurrentElement)

     as(asDOFs*(asPerm(NodeIndexes(k))-1)+1) = as_perp(NodeIndexes(k))
     temp(tempDOFs*(tempPerm(NodeIndexes(k))-1)+1) = temp_s(NodeIndexes(k))

    END DO

  END DO 

  WRITE(Message,'(a)') 'Solver pdd done.'
  CALL INFO(SolverName,Message,Level=2)


CONTAINS

!-------------------------------------------------------------------------------
!> Computation of the positive degree days (PDD) with statistical temperature
!! fluctuations; based on semi-analytical solution by Calov and Greve (2005).
!! Note that this subroutine uses years as time unit
!! (as opposed to the seconds which are usually used by SICOPOLIS).
!<------------------------------------------------------------------------------
SUBROUTINE pdd_f(temp_mm, s_stat, ET)

 IMPLICIT NONE

 REAL(KIND=dp), DIMENSION(12), INTENT(IN) :: temp_mm
 REAL(KIND=dp), INTENT(IN) :: s_stat

 REAL(KIND=dp), INTENT(OUT) :: ET

 INTEGER :: n
 REAL(KIND=dp) :: inv_sqrt2pi, inv_s_stat, inv_sqrt2
 REAL(KIND=dp) :: pdd_sum

 REAL(KIND=dp), PARAMETER :: time_year     = 1.0_dp, &        ! period 1 year
                          time_year_inv = 1.0_dp/time_year, &
                          d_time        = 1.0_dp/12.0_dp   ! time-step 1 month

 inv_sqrt2pi = 1.0_dp/SQRT(2.0_dp*pi)
 inv_s_stat  = 1.0_dp/s_stat
 inv_sqrt2   = 1.0_dp/SQRT(2.0_dp)

 pdd_sum = 0.0_dp

 DO n=1, 12   ! month counter
   pdd_sum = pdd_sum &
             + ( s_stat*inv_sqrt2pi*EXP(-0.5_dp*(temp_mm(n)*inv_s_stat)**2) &
             + 0.5_dp*temp_mm(n)*erfcc(-temp_mm(n)*inv_s_stat*inv_sqrt2) ) &
             *d_time   ! positive degree days (in a * deg C)
 END DO

 ET   = pdd_sum*time_year_inv   ! temperature excess   (in deg C)

END SUBROUTINE pdd_f


!-------------------------------------------------------------------------------
!> Computation of the complementary error function erfc(x) = 1-erf(x)
!! with a fractional error everywhere less than 1.2 x 10^(-7)
!! (formula by Press et al., 'Numerical Recipes in Fortran 77').
!<------------------------------------------------------------------------------
REAL(KIND=dp) function erfcc(x)

IMPLICIT NONE

REAL(KIND=dp), INTENT(IN) :: x

REAL(KIND=dp) :: t, z

z = ABS(x)
t = 1.0_dp/(1.0_dp+0.5_dp*z)

erfcc = t * EXP( -z*z     -1.26551223_dp &
                 + t  * (  1.00002368_dp &
                 + t  * (  0.37409196_dp &
                 + t  * (  0.09678418_dp &
                 + t  * ( -0.18628806_dp &
                 + t  * (  0.27886807_dp &
                 + t  * ( -1.13520398_dp &
                 + t  * (  1.48851587_dp &
                 + t  * ( -0.82215223_dp &
                 + t  *    0.17087277_dp ) ) ) ) ) ) ) ) )

if (x < 0.0_dp) erfcc = 2.0_dp-erfcc

END FUNCTION erfcc


!-------------------------------------------------------------------------------
!> Computation of longitude lambda and latitude phi for position (x,y) in the
!! numerical domain.
!<------------------------------------------------------------------------------
SUBROUTINE geo_coord(phi_val, lambda_val, x_val, y_val, GRID)

IMPLICIT NONE

REAL(KIND=dp), INTENT(IN) :: x_val, y_val

REAL(KIND=dp), INTENT(OUT) :: lambda_val, phi_val

REAL(KIND=dp) :: K

REAL(KIND=dp), PARAMETER :: phi0 = 1.239184d+00       ! Standard parallel +71 deg (71 deg N) in rad
REAL(KIND=dp), PARAMETER :: lambda0 = -0.6806784d+00  ! Reference longitude -39 deg (39 deg W) in rad

REAL(KIND=dp), PARAMETER :: A = 6.378137d+06          ! Semi-major axis of the Earth = 6378137 m
REAL(KIND=dp), PARAMETER :: B = 6.3567523142d+06      ! Semi-minor axis of the Earth = 6356752.3142 m  

INTEGER :: GRID
 
IF (GRID == 0 .OR.  GRID == 1) THEN

!-------- Inverse stereographic projection --------

call stereo_inv_ellipsoid(x_val, y_val, A, B, &
                          lambda0, phi0, lambda_val, phi_val)

ELSE IF (GRID == 2) THEN

!-------- Identity --------

lambda_val = x_val
phi_val    = y_val

END IF

END SUBROUTINE geo_coord

!-------------------------------------------------------------------------------
!> Inverse stereographic projection for an ellipsoidal planet.
!<------------------------------------------------------------------------------
SUBROUTINE stereo_inv_ellipsoid(x_val, y_val, A, B, &
                                lambda0, phi0, lambda_val, phi_val)

IMPLICIT NONE

REAL(kind=dp), INTENT(IN)  :: x_val, y_val, A, B, lambda0, phi0
REAL(KIND=dp), INTENT(OUT) :: lambda_val, phi_val

INTEGER :: l
INTEGER :: sign_phi0
REAL(KIND=dp) :: phi_aux, phi0_aux
REAL(KIND=dp) :: e, mc, t, tc, kp, rho, phi_p, residual
REAL(KIND=dp) :: sinphi0, sinlambda0, cosphi0, coslambda0

REAL(KIND=dp), PARAMETER :: eps = 1.0e-05_dp, &
                         eps_residual = 1.0e-09_dp

IF (phi0 > eps) THEN   ! for northern hemisphere
   sign_phi0 =  1
ELSE IF (phi0 < (-eps)) THEN  ! for southern hemisphere
   sign_phi0 = -1
ELSE
   stop ' geo_coord: phi0 must be different from zero!'
END IF

phi0_aux = phi0 * sign_phi0

e=SQRT((A**2-B**2)/(A**2))

sinphi0    = SIN(phi0_aux)
sinlambda0 = SIN(lambda0)
cosphi0    = COS(phi0_aux)
coslambda0 = COS(lambda0)
  
tc=SQRT(((1.0_dp-sinphi0)/(1.0_dp+sinphi0))* &
       ((1.0_dp+e*sinphi0)/(1.0_dp-e*sinphi0))**e)
mc=cosphi0/SQRT(1.0_dp-e*e*sinphi0*sinphi0)
rho=SQRT(x_val*x_val+y_val*y_val)
t=rho*tc/(A*mc)

lambda_val = lambda0 + sign_phi0*ATAN2(y_val,x_val) + 0.5_dp*pi

!  fix point iteration

phi_p=0.5_dp*pi-2.0_dp*ATAN(t)
l=0
residual=3600.0_dp
DO WHILE(residual >= eps_residual)
   l=l+1
   phi_aux=0.5_dp*pi-2.0_dp*ATAN(t*((1.0_dp-e*SIN(phi_p))/ &
           (1.0_dp+e*SIN(phi_p)))**(0.5_dp*e))
   residual=ABS(phi_aux-phi_p)
   phi_p=phi_aux
END DO

phi_val = phi_aux * sign_phi0

IF (lambda_val < 0.0_dp) THEN
   lambda_val = lambda_val + 2.0_dp*pi
ELSE IF (lambda_val >= (2.0_dp*pi)) THEN
   lambda_val = lambda_val - 2.0_dp*pi
END IF

END SUBROUTINE stereo_inv_ellipsoid



END SUBROUTINE pdd

