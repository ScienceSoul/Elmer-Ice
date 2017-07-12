!----------------------------------------------------------------------------------    
RECURSIVE SUBROUTINE getB1301Velo( Model,Solver,dt,TransientSimulation )
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
 
 INTEGER :: nvalue, i, j, k, n, t, DIM, istat, FLowDOFs, bc_id, &
            Npes, ierr, indx1, indx2
 TYPE(Element_t), POINTER :: Element
 TYPE(ValueList_t),POINTER :: SolverParams
 TYPE(ValueList_t), POINTER :: Material
 TYPE(Nodes_t) :: ElementNodes
 REAL(KIND=dp), POINTER :: Flow(:)
 TYPE(Variable_t), POINTER :: FlowSol
 REAL(KIND=dp) :: B1301_x, B1301_y, velocity, dist, minValue
 REAL(KIND=dp), ALLOCATABLE :: gatherVelo(:), gatherDistance(:), distance(:,:), &
                               storedTriangle(:,:)
 INTEGER, ALLOCATABLE :: indexes(:)
 INTEGER, POINTER :: FlowPerm(:), NodeIndexes(:)
 CHARACTER(LEN=MAX_NAME_LEN) :: SolverName, OutFile, Method
 LOGICAL :: AllocationsDone = .FALSE., Found, done = .FALSE.

 SAVE done

 SolverName = 'getB1301'

 CALL INFO(SolverName,'Getting B1301....')

 SolverParams => Solver % Values

 dim = CoordinateSystemDimension()  

 B1301_x = -663809.8062_dp
 B1301_y = -1168195.263_dp

 FlowSol => VariableGet( Solver % Mesh % Variables, 'Flow Solution')
 IF ( ASSOCIATED( FlowSol ) ) THEN
       FlowPerm    => FlowSol % Perm
       Flow        => FlowSol % Values
       FLowDOFs    =  FlowSol % DOFs
 ELSE
       CALL Fatal(SolverName, 'Could not find Flow Solution variable')
 END IF

 Npes = ParEnv % PEs

 Method = GetString( Solver % Values, 'Point look-up method', Found )
 IF(.NOT. Found) THEN
     CALL Fatal(SolverName, 'Could not find < Point look-up method > parameter')
 END IF

 SELECT CASE(Method)
   CASE('nearest')
       ALLOCATE(gatherVelo(Npes), gatherDistance(Npes), STAT=istat )
       gatherVelo = 0.0_dp
       gatherDistance = 0.0;
       velocity = 0.0_dp

       ALLOCATE(distance(Solver % NumberOfActiveElements, Solver % Mesh % MaxElementNodes), & 
         STAT=istat )
       IF ( istat /= 0 ) THEN
          CALL FATAL( SolverName, 'Memory allocation error' )
       ELSE
          CALL INFO(SolverName, 'Memory allocation done', level=1 )
       END IF

       DO t=1, Solver % NumberOfActiveElements
          Element => GetActiveElement(t)
          n = GetElementNOFNodes(Element)
          CALL GetElementNodes( ElementNodes, Element )
          DO j=1,n
             dist = SQRT( (B1301_x - ElementNodes % x(j))**2.0 + &
                     (B1301_y - ElementNodes % y(j))**2.0 )
             distance(t, j) = dist
          END DO
       END DO

       minValue  = HUGE(minValue)
       DO t=1, Solver % NumberOfActiveElements
          Element => GetActiveElement(t)
          n = GetElementNOFNodes(Element)
          NodeIndexes => Element % NodeIndexes
          DO j=1,n
            IF(distance(t,j) < minValue) THEN
              k = FlowPerm(NodeIndexes(j))
              minValue = distance(t,j)
              velocity = SQRT( Flow(FlowDOFs*(k-1)+1)**2.0_dp + &
                               Flow(FlowDOFs*(k-1)+2)**2.0_dp + &
                               Flow(FlowDOFs*(k-1)+3)**2.0_dp )
            END IF
          END DO
      END DO

      CALL MPI_Gather(velocity,1,MPI_DOUBLE,gatherVelo,1,MPI_DOUBLE,0,MPI_COMM_WORLD,ierr)
      CALL MPI_Gather(minValue,1,MPI_DOUBLE,gatherDistance,1,MPI_DOUBLE,0,MPI_COMM_WORLD,ierr)
 
      DEALLOCATE(distance)

   CASE('triangulate')
       ALLOCATE(gatherVelo(Npes), STAT=istat )
       ALLOCATE(storedTriangle(3,3),  STAT=istat )
       IF ( istat /= 0 ) THEN
          CALL FATAL( SolverName, 'Memory allocation error' )
       ELSE
          CALL INFO(SolverName, 'Memory allocation done', level=1 )
       END IF

       gatherVelo = 0.0_dp
       storedTriangle = 0.0_dp
       velocity = 0.0_dp

       DO t=1, Solver % NumberOfActiveElements
         Element => GetActiveElement(t)
         n = GetElementNOFNodes(Element)
         NodeIndexes => Element % NodeIndexes
         CALL GetElementNodes( ElementNodes, Element ) 
         CALL pointTriangle(B1301_x, B1301_y, ElementNodes % x(1),  ElementNodes % y(1), &
                           ElementNodes % x(2),  ElementNodes % y(2), &
                           ElementNodes % x(3),  ElementNodes % y(3), j)
         IF(j == 1) THEN
            storedTriangle(1,1) = ElementNodes % x(1)
            storedTriangle(1,2) = ElementNodes % y(1)
            k = FlowPerm(NodeIndexes(1))
            storedTriangle(1,3) = SQRT( Flow(FlowDOFs*(k-1)+1)**2.0_dp + &
                                        Flow(FlowDOFs*(k-1)+2)**2.0_dp + &
                                        Flow(FlowDOFs*(k-1)+3)**2.0_dp )

            storedTriangle(2,1) = ElementNodes % x(2)
            storedTriangle(2,2) = ElementNodes % y(2)
            k = FlowPerm(NodeIndexes(2))
            storedTriangle(2,3) = SQRT( Flow(FlowDOFs*(k-1)+1)**2.0_dp + &
                                        Flow(FlowDOFs*(k-1)+2)**2.0_dp + &
                                        Flow(FlowDOFs*(k-1)+3)**2.0_dp )

            storedTriangle(3,1) = ElementNodes % x(3)
            storedTriangle(3,2) = ElementNodes % y(3)
            k = FlowPerm(NodeIndexes(3))
            storedTriangle(3,3) = SQRT( Flow(FlowDOFs*(k-1)+1)**2.0_dp + &
                                        Flow(FlowDOFs*(k-1)+2)**2.0_dp + &
                                        Flow(FlowDOFs*(k-1)+3)**2.0_dp )

         END IF
       END DO
       CALL interpolate(storedTriangle(1,1), storedTriangle(1,2), storedTriangle(2,1), storedTriangle(2,2), storedTriangle(3,1), &
             storedTriangle(3,2), B1301_x, B1301_y,  storedTriangle(1,3),  storedTriangle(2,3) , storedTriangle(3,3), &
             velocity, 2.0_dp)
       CALL MPI_Gather(velocity,1,MPI_DOUBLE,gatherVelo,1,MPI_DOUBLE,0,MPI_COMM_WORLD,ierr)

       DEALLOCATE(storedTriangle)

   CASE DEFAULT
       CALL INFO(SolverName, 'Method look-up not supported')
 END SELECT
        

 IF (ParEnv % MyPe == 0) then
     OutFile = ListGetString(Solver % Values,'Output file name',Found )
     IF(.NOT. Found) THEN
          CALL Fatal(SolverName, 'Could not find name for output file')
     END IF
     OPEN (12, FILE=OutFile)

     SELECT CASE(Method)
        CASE('nearest')
           minValue = MINVAL(gatherDistance)

           DO i=1,Npes
              IF (gatherDistance(i) == minValue) THEN
                WRITE(12,'(e15.8)') gatherVelo(i)
              END IF
           END DO
           DEALLOCATE(gatherVelo, gatherDistance)

        CASE('triangulate')
           DO i=1,Npes
              IF (gatherVelo(i) /= 0.0_dp) THEN
                WRITE(12,'(e15.8)') gatherVelo(i)
              END IF
           END DO
           DEALLOCATE(gatherVelo)

        CASE DEFAULT
           CALL INFO(SolverName, 'Method look-up not supported')
        END SELECT

     CLOSE(12)
 END IF

 CONTAINS

SUBROUTINE pointTriangle(p1, p2, a1, a2, b1, b2, c1, c2, state)

 INTEGER :: state
 REAL(KIND=dp) :: p1, p2, a1, a2, b1, b2, c1, c2
 REAL(KIND=dp) :: v0(2), v1(2), v2(2) 
 REAL(KIND=dp) :: dot00, dot01, dot02, dot11, dot12, &
                  invDenom, u, v

 u = 0.0_dp
 v = 0.0_dp
 invDenom = 0.0_dp

 v0(1) = c1 - a1
 v0(2) = c2 - a2

 v1(1) = b1 - a1
 v1(2) = b2 - a2
 
 v2(1) = p1 - a1 
 v2(2) = p2 - a2

 ! Compute dot products
 dot00 = 0.0_dp
 dot00 = SUM( v0(1:2) * v0(1:2) )

 dot01 = 0.0_dp
 dot01 = SUM( v0(1:2) * v1(1:2) )

 dot02 = 0.0_dp
 dot02 = SUM( v0(1:2) * v2(1:2) )

 dot11 = 0.0_dp
 dot11 = SUM( v1(1:2) * v1(1:2) )

 dot12 = 0.0_dp
 dot12 = SUM( v1(1:2) * v2(1:2) )

 ! Compute barycentric coordinates
 invDenom = 1.0_dp / (dot00 * dot11 - dot01 * dot01)
 u = (dot11 * dot02 - dot01 * dot12) * invDenom
 v = (dot00 * dot12 - dot01 * dot02) * invDenom

 ! Check if point is in triangle
 IF(u > 0.0_dp .AND. v > 0.0_dp .AND. (u+v) < 1.0_dp) THEN
   state = 1
 ELSE 
   state = -1
 END IF
 
END SUBROUTINE pointTriangle

SUBROUTINE interpolate(t11, t12, t21, t22, t31, t32, p1, p2, val1, val2, val3, interVal, exponentVal)

 REAL(KIND=dp) :: t11, t12, t21, t22, t31, t32, p1, p2, &
                  val1, val2, val3
 REAL(KIND=dp) :: distanceToPoint(3), exponentVal, interVal, weight, weightsum

 distanceToPoint = 0.0_dp
 distanceToPoint(1) = SQRT( (p1 - t11)**2.0 + (p2 - t12)**2.0 )
 distanceToPoint(2) = SQRT( (p1 - t21)**2.0 + (p2 - t22)**2.0 )
 distanceToPoint(3) = SQRT( (p1 - t31)**2.0 + (p2 - t32)**2.0 )

 weight = 0.0_dp
 weightsum = 0.0_dp
 interVal = 0.0_dp

 weight = distanceToPoint(1)**(-exponentVal)
 interVal = interVal + weight * val1
 weightsum = weightsum + weight

 weight = distanceToPoint(2)**(-exponentVal)
 interVal = interVal + weight * val2
 weightsum = weightsum + weight

 weight = distanceToPoint(3)**(-exponentVal)
 interVal = interVal + weight * val3
 weightsum = weightsum + weight

 interVal = interVal/weightsum
 
END

!------------------------------------------------------------------------------
END SUBROUTINE getB1301Velo
!------------------------------------------------------------------------------
