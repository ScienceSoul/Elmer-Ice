RECURSIVE SUBROUTINE Bed_SicoAnomalies( Model,Solver,dt,TransientSimulation )
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

   INTEGER :: nvalue, i, j, k, indx, Bed_AnomaliesDOFs, DIM, N, M, countNodes
   REAL(KIND=dp) :: buffer
   REAL(KIND=dp), POINTER :: Bed_AnomaliesValues(:)
   REAL(KIND=dp), ALLOCATABLE :: got_bed_anomalies(:)
   TYPE(Variable_t), POINTER :: Bed_AnomaliesSol
   INTEGER, POINTER :: Bed_AnomaliesPerm(:)
   CHARACTER(len=50) :: filename1, filename2
   INTEGER :: status
   LOGICAL :: FirstTime=.True.

   SAVE FirstTime
   SAVE got_bed_anomalies
   SAVE k

   filename1 ='nb_input.dat'
   filename2 = 'bed_elevation_anomalies.dat'

!-----------------------------------------------------------------------
! get solver variable
!-----------------------------------------------------------------------
   Bed_AnomaliesSol => Solver % Variable
   IF (.NOT.ASSOCIATED(Bed_AnomaliesSol)) THEN
     CALL FATAL('Bed_SicoAnomaliesSolver','No variable associated')
   END IF
   Bed_AnomaliesPerm  => Bed_AnomaliesSol % Perm
   Bed_AnomaliesValues => Bed_AnomaliesSol % Values
   Bed_AnomaliesDOFs =  Bed_AnomaliesSol % DOFs
   Bed_AnomaliesValues = 0.0d00

IF (FirstTime) THEN
      FirstTime=.False.
      OPEN(UNIT=1, FILE=filename1, STATUS='OLD', ACTION='READ', &
               IOSTAT=status)
      IF (status == 0) THEN
         READ(1, *) nvalue
      ELSE
         CALL FATAL('Bed_SicoAnomaliesSolver','Error opening file')
      END IF

      CLOSE (UNIT=1)    

      ALLOCATE(got_bed_anomalies(nvalue), STAT=status)

      IF (status /= 0) THEN
        WRITE(Message,'(a,I6)') 'Memory allocation error. Stat:', status 
        CALL FATAL('Bed_SicoAnomaliesSolver', Message)
      END IF
      
      indx = 0

      OPEN(UNIT=2, FILE=filename2, STATUS='OLD', ACTION='READ', &
               IOSTAT=status)
      IF (status == 0) THEN
      DO
        READ(2,*,IOSTAT=status) buffer
        IF(status /= 0) EXIT
        indx = indx+1
        got_bed_anomalies(indx) = buffer
      END DO
      IF (status > 0) THEN
          WRITE(Message,'(a,I6)') 'An error occured reading line', indx
          CALL FATAL('Bed_SicoAnomaliesSolver', Message)
      END IF
      ELSE
         CALL FATAL('Bed_SicoAnomaliesSolver','Error opening input file')
      END IF

      CLOSE (UNIT=2)
      
      k=1
      
      DIM = CoordinateSystemDimension()
      N = Solver % Mesh % MaxElementNodes
      M = Model % Mesh % NumberOfNodes

      WRITE(Message,'(a,I6,a)') 'The file nb_input contains', nvalue, ' values'
      CALL Info('Bed_SicoAnomaliesSolver',Message, Level=3)
      WRITE(Message,'(a,I6,a)') 'The table got_bed_anomalies has', indx, ' values'
      CALL Info('Bed_SicoAnomaliesSolver',Message, Level=3)        
   END IF

 !-----------------------------------------------------------------------
 ! loop over all nodes in mesh
 !-----------------------------------------------------------------------
   countNodes = 0
   DO i=1,Solver % Mesh % NumberOfNodes
      
      j = Bed_AnomaliesPerm(i)
      Bed_AnomaliesValues(Bed_AnomaliesDOFs*(j-1)+1) = got_bed_anomalies(k)
      countNodes = countNodes+1
   END DO

   k = k+1

   IF(countNodes == Solver % Mesh % NumberOfNodes) THEN
     CALL Info('Bed_SicoAnomaliesSolver','Got the Sicopolis bed elevation anomalies', Level=3)
   ELSE
     CALL FATAL('Bed_SicoAnomaliesSolver','Error getting the Sicopolis bed elevation anomalies')
   END IF

   WRITE(Message,'(a)') '-----------------------------------------------'
   CALL Info('Bed_SicoAnomaliesSolver',Message, Level=3)

!------------------------------------------------------------------------------
 END SUBROUTINE  Bed_SicoAnomalies
!------------------------------------------------------------------------------
