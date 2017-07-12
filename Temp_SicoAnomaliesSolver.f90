!----------------------------------------------------------------------------------    
RECURSIVE SUBROUTINE Temp_SicoAnomalies( Model,Solver,dt,TransientSimulation )
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

   INTEGER :: nvalue, i, j, k, indx, Temp_AnomaliesDOFs, DIM, N, M, countNodes
   REAL(KIND=dp) :: buffer
   REAL(KIND=dp), POINTER :: Temp_AnomaliesValues(:)
   REAL(KIND=dp), ALLOCATABLE :: got_temp_anomalies(:)
   TYPE(Variable_t), POINTER :: Temp_AnomaliesSol
   INTEGER, POINTER :: Temp_AnomaliesPerm(:)
   CHARACTER(len=50) :: filename1, filename2
   INTEGER :: status
   LOGICAL :: FirstTime=.True.

   SAVE FirstTime
   SAVE got_temp_anomalies
   SAVE k

   filename1 ='nb_input.dat'
   filename2 = 'surface_temperature_anomalies.dat'

!-----------------------------------------------------------------------
! get solver variable
!-----------------------------------------------------------------------
   Temp_AnomaliesSol => Solver % Variable
   IF (.NOT.ASSOCIATED(Temp_AnomaliesSol)) THEN
     CALL FATAL('Temp_SicoAnomaliesSolver','No variable associated')
   END IF
   Temp_AnomaliesPerm  => Temp_AnomaliesSol % Perm
   Temp_AnomaliesValues => Temp_AnomaliesSol % Values
   Temp_AnomaliesDOFs =  Temp_AnomaliesSol % DOFs
   Temp_AnomaliesValues = 0.0d00
   
   IF (FirstTime) THEN
      FirstTime=.False.
      OPEN(UNIT=1, FILE=filename1, STATUS='OLD', ACTION='READ', &
               IOSTAT=status)
      IF (status == 0) THEN
         READ(1, *) nvalue
      ELSE
         CALL FATAL('Temp_SicoAnomaliesSolver','Error opening file')
      END IF

      CLOSE (UNIT=1)    

      ALLOCATE(got_temp_anomalies(nvalue), STAT=status)

      IF (status /= 0) THEN
        WRITE(Message,'(a,I6)') 'Memory allocation error. Stat:', status 
        CALL FATAL('Temp_SicoAnomaliesSolver', Message)
      END IF
      
      indx = 0

      OPEN(UNIT=2, FILE=filename2, STATUS='OLD', ACTION='READ', &
               IOSTAT=status)
      IF (status == 0) THEN
      DO
        READ(2,*,IOSTAT=status) buffer
        IF(status /= 0) EXIT
        indx = indx+1
        got_temp_anomalies(indx) = buffer
      END DO
      IF (status > 0) THEN
          WRITE(Message,'(a,I6)') 'An error occured reading line', indx
          CALL FATAL('Temp_SicoAnomaliesSolver', Message)
      END IF
      ELSE
         CALL FATAL('Temp_SicoAnomaliesSolver','Error opening input file')
      END IF

      CLOSE (UNIT=2)
      
      k=1
      
      DIM = CoordinateSystemDimension()
      N = Solver % Mesh % MaxElementNodes
      M = Model % Mesh % NumberOfNodes

     WRITE(Message,'(a,I6,a)') 'The file nb_input contains', nvalue, ' values'
     CALL Info('Temp_SicoAnomaliesSolver',Message, Level=3)
     WRITE(Message,'(a,I6,a)') 'The table got_temp_anomalies has', indx, ' values'
     CALL Info('Temp_SicoAnomaliesSolver',Message, Level=3)        
   END IF


   !-----------------------------------------------------------------------
   ! loop over all nodes in mesh
   !-----------------------------------------------------------------------
   countNodes = 0
   DO i=1,Solver % Mesh % NumberOfNodes
      
      j = Temp_AnomaliesPerm(i)
      Temp_AnomaliesValues(Temp_AnomaliesDOFs*(j-1)+1) = got_temp_anomalies(k)
      countNodes = countNodes+1
   END DO

   k = k+1
   
   IF(countNodes == Solver % Mesh % NumberOfNodes) THEN
      CALL Info('Temp_SicoAnomaliesSolver','Got the Sicopolis temperatures anomalies', Level=3)
   ELSE
      CALL FATAL('Temp_SicoAnomaliesSolver','Error getting the Sicopolis temperature anomalies')
   END IF
   
   WRITE(Message,'(a)') '-----------------------------------------------'
   CALL Info('Temp_SicoAnomaliesSolver',Message, Level=3)

!------------------------------------------------------------------------------
 END SUBROUTINE  Temp_SicoAnomalies
!------------------------------------------------------------------------------
