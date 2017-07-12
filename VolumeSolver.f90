RECURSIVE SUBROUTINE VolumeSolver ( Model,Solver,Timestep,TransientSimulation )
  USE DefUtils
    
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
 INTEGER :: i,j,k,t,n,m,DIM,istat,incr,TimeDOFs,TimeSimul,ierror
 TYPE(Element_t), POINTER :: CurrentElement
 TYPE(Nodes_t) :: Nodes 
 TYPE(GaussIntegrationPoints_t) :: IP 
 INTEGER, POINTER :: NodeIndexes(:)
 TYPE(Variable_t), POINTER :: TimeSol
 REAL(KIND=dp), POINTER :: Time(:)
 INTEGER, POINTER :: TimePerm(:)
 REAL(KIND=dp) :: ElementVolume, Volume
 REAL(KIND=dp) :: DetJ
 REAL(KIND=dp), ALLOCATABLE :: Basis(:), dBasisdx(:,:), ddBasisddx(:,:,:), buffer(:)
 CHARACTER(LEN=MAX_NAME_LEN) :: SolverName, string, FILENAME
 LOGICAL :: AllocationsDone = .FALSE., stat, GotIt

 SAVE SolverName, Basis, dBasisdx, ddBasisddx, AllocationsDone, buffer, incr, TimeSimul, FILENAME

 IF ( .NOT. AllocationsDone .OR. Solver % Mesh % Changed ) THEN

        SolverName = 'VolumeSolver'

        DIM = CoordinateSystemDimension()
        N = Solver % Mesh % MaxElementNodes
        M = Model % Mesh % NumberOfNodes

        IF ( AllocationsDone ) THEN
           DEALLOCATE(          &
               Basis,         &
               dBasisdx,      &
               ddBasisddx,    &
               buffer)              
        END IF                           
        
        ALLOCATE(                   &
               Basis(N),          &
               dBasisdx(N,3),     &
               ddBasisddx(N,3,3), &
               buffer(M),         &
               STAT=istat )

        IF ( istat /= 0 ) THEN
           CALL FATAL( SolverName, 'Memory allocation error' )
        ELSE
           CALL INFO(SolverName, 'Memory allocation done', level=1 )
        END IF 

        incr = 1

        IF ( TransientSimulation ) THEN

         TimeSimul = GetInteger( GetSimulation(),'Timestep Intervals',GotIt)
         IF (.NOT. GotIt) THEN
           CALL FATAL(SolverName,'Could not find parameter Timestep Intervals')
         END IF

        END IF

        FILENAME = GetString( Solver % Values,'Output Filename',GotIt)
        IF (.NOT. GotIt) THEN
           CALL FATAL(SolverName,'Could not find parameter FILENAME')
        END IF

        AllocationsDone = .TRUE.
 END IF
 
 Volume = 0.0_dp

 DO t=1,Solver % NumberOFActiveElements

    CurrentElement => GetActiveElement(t)
    NodeIndexes => CurrentElement % NodeIndexes

    CALL GetElementNodes( Nodes, CurrentElement ) 
    IP = GaussPoints( CurrentElement )

    ElementVolume = 0.0_dp

    !---------------------------------------------------------------- 
    ! Loop over Gauss-points (element Integration) 
    !---------------------------------------------------------------- 
    DO i=1,IP % n

       stat = ElementInfo( CurrentElement, Nodes, IP % U(i), IP % V(i), & 
                 IP % W(i), DetJ, Basis, dBasisdx, ddBasisddx, .FALSE. )

       ElementVolume = ElementVolume + IP % s(i) * DetJ
    
    END DO
   
    Volume = Volume+ElementVolume

 END DO

 buffer(incr) = Volume/1000000000.0
 !WRITE(*,*) buffer(1:incr)
 
IF ( TransientSimulation ) THEN

    TimeSol => VariableGet( Solver % Mesh % Variables, 'Time' )
    IF ( ASSOCIATED( TimeSol ) ) THEN
       TimePerm    => TimeSol % Perm
       Time        => TimeSol % Values
       TimeDOFs = TimeSol % DOFs
    ELSE
       CALL FATAL(SolverName, ' Could not find Time field variable')
    END IF
    
    IF ( incr == TimeSimul) THEN
       ! serial run 
       IF ( ParEnv % PEs <= 1 ) THEN 
    
          OPEN (UNIT=10, FILE=FILENAME, STATUS='REPLACE', ACTION='WRITE', IOSTAT=ierror)
          IF( ierror == 0) THEN 
            DO i=1,incr
               WRITE(10,*) buffer(i)
            END DO
            CLOSE(10)
          ELSE 
            CALL FATAL(SolverName, ' Error opening file for writing volume')
          END IF

       ELSE ! parallel run 
          
          WRITE (string,'(I0)') ParEnv % myPe
          OPEN (UNIT=10, FILE=TRIM(string) // FILENAME, STATUS='REPLACE', ACTION='WRITE', &
               IOSTAT=ierror)
          IF( ierror == 0) THEN 
            DO i=1,incr
               WRITE(10,*) buffer(i)
            END DO
            CLOSE(10)
          ELSE 
            CALL FATAL(SolverName, ' Error opening file for writing volume')
          END IF
  
       ENDIF

    END IF 

 ELSE IF ( .NOT. TransientSimulation ) THEN

    IF ( ParEnv % PEs <= 1 ) THEN 
    
        OPEN (UNIT=10, FILE=FILENAME, STATUS='REPLACE', ACTION='WRITE', IOSTAT=ierror)
        IF( ierror == 0) THEN 
           DO i=1,incr
              WRITE(10,*) buffer(i)
           END DO
           CLOSE(10)
        ELSE 
           CALL FATAL(SolverName, ' Error opening file for writing volume')
        END IF

    ELSE ! parallel run 
        
        WRITE (string,'(I0)') ParEnv % myPe
        OPEN (UNIT=10, FILE=TRIM(string) // FILENAME, STATUS='REPLACE', ACTION='WRITE', &
               IOSTAT=ierror)
        IF( ierror == 0) THEN 
           DO i=1,incr
              WRITE(10,*) buffer(i)
           END DO
           CLOSE(10)
        ELSE 
           CALL FATAL(SolverName, ' Error opening file for writing volume')
        END IF
  
    ENDIF

 END IF

 incr = incr+1

 WRITE(Message,'(a)') 'Volume Solver done.'
 CALL INFO(SolverName,Message,Level=2)
   

END SUBROUTINE VolumeSolver
