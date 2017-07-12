!----------------------------------------------------------------------------------    
RECURSIVE SUBROUTINE getVeloAtProfiles( Model,Solver,dt,TransientSimulation )
!---------------------------------------------------------------------------------- 

  USE DefUtils

  IMPLICIT NONE
!------------------------------------------------------------------------------

   TYPE(Model_t)  :: Model
   TYPE(Solver_t), TARGET :: Solver

   LOGICAL ::  TransientSimulation
   REAL(KIND=dp) :: dt
!------------------------------------------------------------------------------
!    Local variables
!------------------------------------------------------------------------------
 
 INTEGER :: nvalue, i, j, k, l, n, t, DIM, istat, FLowDOFs, bc_id, &
            Npes, ierr, nbd1, nbd2, ni1, ni2,  tot1, tot2
 INTEGER :: status(MPI_STATUS_SIZE)
 TYPE(Element_t), POINTER :: Element
 TYPE(ValueList_t),POINTER :: SolverParams, BC
 TYPE(ValueList_t), POINTER :: Material
 TYPE(Nodes_t) :: ElementNodes
 REAL(KIND=dp), POINTER :: Flow(:), Depth(:), Height(:)
 TYPE(Variable_t), POINTER :: FlowSol, DepthSol, HeightSol
 REAL(KIND=dp) :: velocity, vector(3), rotVector(3)
 REAL(KIND=dp), ALLOCATABLE :: gatherVeloSurf(:), gatherSurfX(:), &
                               gatherVeloBed(:), gatherBedX(:), &
                               VeloSurf(:), XSurf(:), VeloBed(:), XBed(:)
 REAL(KIND=dp), ALLOCATABLE :: SortedVeloSurf(:), SortedSurfX(:), SortedVeloBed(:), &
                               SortedBedX(:), R(:,:)
 INTEGER, ALLOCATABLE :: indexes(:), NodePerPe1(:), NodePerPe2(:)
 INTEGER, POINTER :: FlowPerm(:), DepthPerm(:), NodeIndexes(:), HeightPerm(:)
 CHARACTER(LEN=MAX_NAME_LEN) :: SolverName, OutFile1, OutFile2, BCName, &
                                transetName
 LOGICAL :: AllocationsDone = .FALSE., Found, done = .FALSE., rotate

 SAVE done

 SolverName = 'getVeloAtProfiles'

 SolverParams => Solver % Values

 dim = CoordinateSystemDimension()  

 FlowSol => VariableGet( Solver % Mesh % Variables, 'Flow Solution')
 IF ( ASSOCIATED( FlowSol ) ) THEN
       FlowPerm    => FlowSol % Perm
       Flow        => FlowSol % Values
       FLowDOFs    =  FlowSol % DOFs
 ELSE
       CALL Fatal(SolverName, 'Could not find Flow Solution variable')
 END IF

 DepthSol => VariableGet( Solver % Mesh % Variables, 'Depth')
 IF ( ASSOCIATED( DepthSol ) ) THEN
       DepthPerm    => DepthSol % Perm
       Depth        => DepthSol % Values
 ELSE
       CALL Fatal(SolverName, 'Could not find Depth variable')
 END IF

 HeightSol => VariableGet( Solver % Mesh % Variables, 'Height')
 IF ( ASSOCIATED( HeightSol ) ) THEN
       HeightPerm    => HeightSol % Perm
       Height        => HeightSol % Values
 ELSE
       CALL Fatal(SolverName, 'Could not find Height variable')
 END IF

 transetName = GetString( Solver % Values, 'transect name', Found )
 IF(.NOT. Found) THEN
    CALL Fatal(SolverName, 'Could not find < transect name > parameter')
 END IF
 WRITE(Message,'(a,a)') 'Output surface and bed nodes for transect: ', transetName
 CALL Info(SolverName, Message, Level=4 )

 rotate = GetLogical(Solver % Values, 'rotate velocity vector', Found)
 IF(.NOT. Found) THEN
    CALL Fatal(SolverName, 'Could not find < rotate velocity vector > parameter')
 END IF

 Npes = ParEnv % PEs

 ALLOCATE(NodePerPe1(Npes), NodePerPe2(NPes))
 IF ( istat /= 0 ) THEN
        CALL FATAL( SolverName, 'Memory allocation error' )
 END IF
       
 nbd1 = 0
 nbd2 = 0
 DO t=1, Solver % Mesh % NumberOfBoundaryElements
    Element => GetBoundaryElement(t)
    
    BC => GetBC()
    IF ( .NOT. ASSOCIATED(BC) ) CYCLE    
    BCName =  ListGetString( BC,'Name', Found)

    CALL GetElementNodes( ElementNodes, Element )
    n = GetElementNOFNodes(Element)
    NodeIndexes => Element % NodeIndexes
    IF (BCName == transetName) THEN
        DO j=1,n
           k = FlowPerm(NodeIndexes(j))
           IF(ABS(Depth(DepthPerm(NodeIndexes(j)))-0.0_dp) <= EPSILON(0.0_dp)) THEN
                 nbd1 = nbd1 + 1
           END IF
           IF(ABS(Height(HeightPerm(NodeIndexes(j)))-0.0_dp) <= EPSILON(0.0_dp)) THEN
                 nbd2 = nbd2 + 1
           END IF 
        END DO
    END IF
 END DO

 CALL MPI_Gather(nbd1,1,MPI_Integer,NodePerPe1,1,MPI_Integer,0,MPI_COMM_WORLD,ierr)
 CALL MPI_Gather(nbd2,1,MPI_Integer,NodePerPe2,1,MPI_Integer,0,MPI_COMM_WORLD,ierr)

 IF (ParEnv % MyPe == 0) then
     DO i=1,Npes
        WRITE(Message,'(a,i4,a,i4)') 'Process rank: ',i ,' found number of surface nodes: ', NodePerPe1(i)
        CALL Info(SolverName, Message, Level=4 )
    END DO
    DO i=1,Npes
        WRITE(Message,'(a,i4,a,i4)') 'Process rank: ',i ,' found number of bed nodes: ', NodePerPe2(i)
        CALL Info(SolverName, Message, Level=4 )
    END DO
 END IF

 IF (ParEnv % MyPe == 0) then

    tot1 = SUM(NodePerPe1)
    tot2 = SUM(NodePerPe2)
    ALLOCATE(gatherVeloSurf(tot1), gatherSurfX(tot1), STAT=istat)
    ALLOCATE(gatherVeloBed(tot2), gatherBedX(tot2), STAT=istat)
    IF ( istat /= 0 ) THEN
        CALL FATAL( SolverName, 'Memory allocation error' )
    END IF

    gatherVeloSurf = 0.0_dp
    gatherVeloBed = 0.0_dp

    WRITE(Message,'(a,i4)') 'Total surface nodes found: ', tot1
    CALL Info(SolverName, Message, Level=4 )
    WRITE(Message,'(a,i4)') 'Total bed nodes found: ', tot2
    CALL Info(SolverName, Message, Level=4 )
 END IF

 IF (nbd1 /= 0 .AND. nbd2 /= 0) THEN
    ALLOCATE(VeloSurf(nbd1), XSurf(nbd1), STAT=istat )
    ALLOCATE(VeloBed(nbd2), XBed(nbd2), STAT=istat )
    IF ( istat /= 0 ) THEN
        CALL FATAL( SolverName, 'Memory allocation error' )
    END IF

    ALLOCATE(R(3,3))
    IF ( istat /= 0 ) THEN
        CALL FATAL( SolverName, 'Memory allocation error' )
    END IF
       
    IF(transetName == 'transverse5') THEN
        R(1,1) = COS(230_dp*PI/180.0_dp)
        R(1,2) = -SIN(230_dp*PI/180.0_dp)
        R(1,3) = 0.0_dp
        R(2,1) = SIN(230_dp*PI/180.0_dp)
        R(2,2) = COS(230_dp*PI/180.0_dp)
        R(2,3) = 0.0_dp
        R(3,1) = 0.0_dp
        R(3,2) = 0.0_dp
        R(3,3) = 1.0_dp
    ELSE IF (transetName == 'transverse6') THEN
        R(1,1) = COS(205_dp*PI/180.0_dp)
        R(1,2) = -SIN(205_dp*PI/180.0_dp)
        R(1,3) = 0.0_dp
        R(2,1) = SIN(205_dp*PI/180.0_dp)
        R(2,2) = COS(205_dp*PI/180.0_dp)
        R(2,3) = 0.0_dp
        R(3,1) = 0.0_dp
        R(3,2) = 0.0_dp
        R(3,3) = 1.0_dp
    END IF

    VeloSurf = 0.0_dp
    XSurf = 0.0_dp
    VeloBed = 0.0_dp
    XBed = 0.0_dp

    i = 0
    l = 0
    DO t=1, Solver % Mesh % NumberOfBoundaryElements
        Element => GetBoundaryElement(t)
    
        BC => GetBC()
        IF ( .NOT. ASSOCIATED(BC) ) CYCLE    
        BCName =  ListGetString( BC,'Name', Found)

        CALL GetElementNodes( ElementNodes, Element )
        n = GetElementNOFNodes(Element)
        NodeIndexes => Element % NodeIndexes
        IF (BCName == transetName) THEN
            DO j=1,n
                k = FlowPerm(NodeIndexes(j))
                 IF(ABS(Depth(DepthPerm(NodeIndexes(j)))-0.0_dp) <= EPSILON(0.0_dp)) THEN
                    IF (rotate) THEN
                        vector(1) = Flow(FlowDOFs*(k-1)+1)
                        vector(2) = Flow(FlowDOFs*(k-1)+2)
                        vector(3) = Flow(FlowDOFs*(k-1)+3)
                        rotVector = MATMUL(R,vector)
                        velocity = SQRT(rotVector(1)**2.0_dp + &
                                        rotVector(2)**2.0_dp + &
                                        rotVector(3)**2.0_dp )
                    ELSE 
                        velocity = SQRT(Flow(FlowDOFs*(k-1)+1)**2.0_dp + &
                                        Flow(FlowDOFs*(k-1)+2)**2.0_dp + &
                                        Flow(FlowDOFs*(k-1)+3)**2.0_dp )
                    END IF
                    i = i + 1
                    VeloSurf(i) = velocity
                    XSurf(i) =  ElementNodes % x(j)
                ELSE IF(ABS(Height(HeightPerm(NodeIndexes(j)))-0.0_dp) <= EPSILON(0.0_dp)) THEN
                    IF (rotate) THEN
                        vector(1) = Flow(FlowDOFs*(k-1)+1)
                        vector(2) = Flow(FlowDOFs*(k-1)+2)
                        vector(3) = Flow(FlowDOFs*(k-1)+3)
                        rotVector = MATMUL(R,vector)
                        velocity = SQRT(rotVector(1)**2.0_dp + &
                                        rotVector(2)**2.0_dp + &
                                        rotVector(3)**2.0_dp )
                    ELSE 
                        velocity = SQRT(Flow(FlowDOFs*(k-1)+1)**2.0_dp + &
                                        Flow(FlowDOFs*(k-1)+2)**2.0_dp + &
                                        Flow(FlowDOFs*(k-1)+3)**2.0_dp )
                    END IF
                    l = l + 1
                    VeloBed(l) = velocity
                    XBed(l) =  ElementNodes % x(j)
                END IF
            END DO
        END IF
    END DO

 END IF

 if (ParEnv % MyPe /= 0) then

    IF (nbd1 /= 0 .AND. nbd2 /= 0) THEN
        CALL MPI_SEND(VeloSurf(1),nbd1,MPI_DOUBLE_PRECISION,0,8003,MPI_COMM_WORLD,ierr)
        CALL MPI_SEND(XSurf(1),nbd1,MPI_DOUBLE_PRECISION,0,8004,MPI_COMM_WORLD,ierr)

        CALL MPI_SEND(VeloBed(1),nbd2,MPI_DOUBLE_PRECISION,0,8005,MPI_COMM_WORLD,ierr)
        CALL MPI_SEND(XBed(1),nbd2,MPI_DOUBLE_PRECISION,0,8006,MPI_COMM_WORLD,ierr)
    END IF

 ELSE
    
    IF(nbd1 /= 0 .AND. nbd2 /= 0) THEN
       gatherVeloSurf(1:nbd1) = VeloSurf(1:nbd1)
       gatherSurfX(1:nbd1) = XSurf(1:nbd1)
       ni1 = 1 + nbd1

       gatherVeloBed(1:nbd2) = VeloBed(1:nbd2)
       gatherBedX(1:nbd2) = XBed(1:nbd2)
       ni2 = 1 + nbd2
    ELSE
       ni1 = 1
       ni2 = 1
    END IF
    DO i=2,Npes
        IF(NodePerPe1(i) /= 0) THEN
            call MPI_RECV(gatherVeloSurf(ni1),NodePerPe1(i),MPI_DOUBLE_PRECISION,i-1,8003,MPI_COMM_WORLD, status, ierr )
            call MPI_RECV(gatherSurfX(ni1),NodePerPe1(i),MPI_DOUBLE_PRECISION,i-1,8004,MPI_COMM_WORLD, status, ierr )
            ni1 = ni1 + NodePerPe1(i)
        END IF
    END DO

    DO i=2,Npes
        IF(NodePerPe2(i) /= 0) THEN
            call MPI_RECV(gatherVeloBed(ni2),NodePerPe2(i),MPI_DOUBLE_PRECISION,i-1,8005,MPI_COMM_WORLD, status, ierr )
            call MPI_RECV(gatherBedX(ni2),NodePerPe2(i),MPI_DOUBLE_PRECISION,i-1,8006,MPI_COMM_WORLD, status, ierr )
            ni2 = ni2 + NodePerPe2(i)
        END IF
    END DO

 END IF

 IF (ParEnv % MyPe == 0) THEN
    ALLOCATE(SortedVeloSurf(tot1), SortedSurfX(tot1), STAT=istat)
    ALLOCATE(SortedVeloBed(tot2), SortedBedX(tot2), STAT=istat)
    IF ( istat /= 0 ) THEN
        CALL FATAL( SolverName, 'Memory allocation error' )
    END IF

    SortedVeloSurf = 0.0_dp
    SortedSurfX = 0.0_dp
    SortedVeloBed = 0.0_dp
    SortedBedX = 0.0_dp

    k = 1
    DO i=1, tot1
        Found = .FALSE.
        Do j=1, tot1
           IF(SortedSurfX(j) == gatherSurfX(i)) THEN
              Found = .TRUE.
              EXIT
           END IF
        END DO
        IF (.NOT. Found) THEN
            SortedVeloSurf(k) = gatherVeloSurf(i)
            SortedSurfX(k) = gatherSurfX(i)
             k = k + 1
        END IF
    END DO

    k = 1
    DO i=1, tot2
        Found = .FALSE.
        Do j=1, tot2
           IF(SortedBedX(j) == gatherBedX(i)) THEN
              Found = .TRUE.
              EXIT
           END IF
        END DO
        IF (.NOT. Found) THEN
            SortedVeloBed(k) = gatherVeloBed(i)
            SortedBedX(k) = gatherBedX(i)
             k = k + 1
        END IF
    END DO

   CALL SortArray( tot1,SortedSurfX,SortedVeloSurf)
   CALL SortArray( tot2,SortedBedX,SortedVeloBed)

   OutFile1 = ListGetString(Solver % Values,'Surface transect file name',Found )
   IF(.NOT. Found) THEN
        CALL Fatal(SolverName, 'Could not find < Surface transect file name >')
   END IF
   OutFile2 = ListGetString(Solver % Values,'Bed transect file name',Found )
   IF(.NOT. Found) THEN
        CALL Fatal(SolverName, 'Could not find < Bed transect file name >')
   END IF
   OPEN (12, FILE=OutFile1)
   OPEN (13, FILE=OutFile2)
   DO i=1, tot1
      IF(SortedSurfX(i) /= 0.0_dp) THEN
         WRITE(12,'(e15.8,e15.8)') SortedSurfX(i), SortedVeloSurf(i)
      END IF
   END DO
   DO i=1, tot2
      IF(SortedBedX(i) /= 0.0_dp) THEN
         WRITE(13,'(e15.8,e15.8)') SortedBedX(i), SortedVeloBed(i)
      END IF
   END DO
   CLOSE(12)
   CLOSE(13)
 END IF

 DEALLOCATE(NodePerPe1)
 DEALLOCATE(NodePerPe2)

 IF (ParEnv % MyPe == 0) then
    DEALLOCATE(gatherVeloSurf, gatherSurfX, &
               gatherVeloBed, gatherBedX)
    DEALLOCATE(SortedVeloSurf, SortedSurfX, &
               SortedVeloBed, SortedBedX)
 END IF

 IF (nbd1 /= 0 .AND. nbd2 /= 0) THEN
    DEALLOCATE(VeloSurf, XSurf, &
               VeloBed, XBed)
    DEALLOCATE(R)
 END IF

 CONTAINS

!------------------------------------------------------------------------------
!> Sort an real array, and change the order of an index array accordingly.
!------------------------------------------------------------------------------
 SUBROUTINE SortArray( n,a,b )
!------------------------------------------------------------------------------
     INTEGER :: n
     REAL(KIND=dp) :: a(:), b(:)
!------------------------------------------------------------------------------

     INTEGER :: i,j,l,ir
     REAL(KIND=dp) :: ra, rb
!------------------------------------------------------------------------------

      IF ( n <= 1 ) RETURN
 
      l = n / 2 + 1
      ir = n
      DO WHILE( .TRUE. )

        IF ( l > 1 ) THEN
          l = l - 1
          ra = a(l)
          rb = b(l)
        ELSE
          ra = a(ir)
          rb = b(ir)
          a(ir) = a(1)
          b(ir) = b(1)
          ir = ir - 1
          IF ( ir == 1 ) THEN
            a(1) = ra
            b(1) = rb
            RETURN
          END IF
        END IF
        i = l
        j = l + l
        DO WHILE( j <= ir )
          IF ( j<ir  ) THEN
            IF ( a(j)<a(j+1) ) j = j+1
          END IF
          IF ( ra<a(j) ) THEN
            a(i) = a(j)
            b(i) = b(j)
            i = j
            j = j + i
          ELSE
            j = ir + 1
          END IF
          a(i) = ra
          b(i) = rb
       END DO
     END DO

 END SUBROUTINE SortArray

!------------------------------------------------------------------------------
END SUBROUTINE getVeloAtProfiles
!------------------------------------------------------------------------------
