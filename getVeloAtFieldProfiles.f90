!----------------------------------------------------------------------------------    
RECURSIVE SUBROUTINE getVeloAtFieldProfiles( Model,Solver,dt,TransientSimulation )
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
 
 INTEGER :: nvalue, i, ii, j, k, n, t, tt, DIM, istat, FLowDOFs, bc_id, &
            Npes, ierr, nbNodes
 TYPE(Element_t), POINTER :: Element
 TYPE(ValueList_t),POINTER :: SolverParams
 TYPE(ValueList_t), POINTER :: Material
 TYPE(Nodes_t) :: ElementNodes
 REAL(KIND=dp), POINTER :: Flow(:)
 TYPE(Variable_t), POINTER :: FlowSol
 REAL(KIND=dp) :: T1_X(17), T1_Y(17), T2_X(22), T2_Y(22), T3_X(31), T3_Y(31), &
                  L_X(30), L_Y(30), velocity, dist, minValue, x, y
 REAL(KIND=dp), ALLOCATABLE :: gatherVelo(:), gatherDistance(:), distance(:,:), &
                               storedTriangle(:,:)
 INTEGER, ALLOCATABLE :: indexes(:)
 INTEGER, POINTER :: FlowPerm(:), NodeIndexes(:)
 CHARACTER(LEN=MAX_NAME_LEN) :: SolverName, OutFile
 LOGICAL :: AllocationsDone = .FALSE., Found, done = .FALSE.

 SAVE done

 SolverName = 'getVeloAtFieldProfiles'

 CALL INFO(SolverName,'Getting Velocities....')

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

 Npes = ParEnv % PEs

 T1_X(1)  = -663966.2921_dp
 T1_X(2)  = -663879.4264_dp
 T1_X(3)  = -663809.8062_dp
 T1_X(4)  = -663743.0009_dp
 T1_X(5)  = -663654.8937_dp
 T1_X(6)  = -663605.4159_dp
 T1_X(7)  = -663507.8737_dp
 T1_X(8)  = -663449.0433_dp
 T1_X(9)  = -663397.5884_dp
 T1_X(10) = -663341.884_dp
 T1_X(11) = -663284.3923_dp
 T1_X(12) = -663226.9605_dp
 T1_X(13) = -663161.9942_dp
 T1_X(14) = -663127.2479_dp
 T1_X(15) = -663070.7474_dp
 T1_X(16) = -663008.4851_dp
 T1_X(17) = -662971.3888_dp
 
 T1_Y(1)  = -1167986.661_dp
 T1_Y(2)  = -1168097.439_dp
 T1_Y(3)  = -1168195.263_dp
 T1_Y(4)  = -1168285.12_dp
 T1_Y(5)  = -1168326.826_dp
 T1_Y(6)  = -1168394.045_dp
 T1_Y(7)  = -1168536.414_dp
 T1_Y(8)  = -1168646.08_dp
 T1_Y(9)  = -1168710.948_dp
 T1_Y(10) = -1168807.375_dp
 T1_Y(11) = -1168862.986_dp
 T1_Y(12) = -1168954.012_dp
 T1_Y(13) = -1169033.821_dp
 T1_Y(14) = -1169124.136_dp
 T1_Y(15) = -1169208.732_dp
 T1_Y(16) = -1169302.461_dp
 T1_Y(17) = -1169371.681_dp

 T2_X(1)  = -663344.1996_dp
 T2_X(2)  = -663300.2985_dp
 T2_X(3)  = -663263.3585_dp
 T2_X(4)  = -663197.1934_dp
 T2_X(5)  = -663090.8422_dp
 T2_X(6)  = -663023.1106_dp
 T2_X(7)  = -662948.8732_dp
 T2_X(8)  = -662861.6233_dp
 T2_X(9)  = -662774.0779_dp
 T2_X(10) = -662689.7394_dp
 T2_X(11) = -662619.4076_dp
 T2_X(12) = -662552.2841_dp
 T2_X(13) = -662426.8616_dp
 T2_X(14) = -662061.667_dp
 T2_X(15) = -661929.9573_dp
 T2_X(16) = -661837.9219_dp
 T2_X(17) = -661765.6161_dp
 T2_X(18) = -661703.1267_dp
 T2_X(19) = -661618.8854_dp
 T2_X(20) = -661572.554_dp
 T2_X(21) = -661511.6536_dp
 T2_X(22) = -661380.0387_dp

 T2_Y(1)  = -1166130.806_dp
 T2_Y(2)  = -1166198.599_dp
 T2_Y(3)  = -1166301.062_dp
 T2_Y(4)  = -1166368.024_dp
 T2_Y(5)  = -1166444.66_dp
 T2_Y(6)  = -1166509.778_dp
 T2_Y(7)  = -1166564.402_dp
 T2_Y(8)  = -1166645.485_dp
 T2_Y(9)  = -1166722.177_dp
 T2_Y(10) = -1166785.383_dp
 T2_Y(11) = -1166851.14_dp
 T2_Y(12) = -1166915.637_dp
 T2_Y(13) = -1167029.919_dp
 T2_Y(14) = -1167361.951_dp
 T2_Y(15) = -1167482.461_dp
 T2_Y(16) = -1167558.816_dp
 T2_Y(17) = -1167607.368_dp
 T2_Y(18) = -1167699.416_dp
 T2_Y(19) = -1167759.183_dp
 T2_Y(20) = -1167834.586_dp
 T2_Y(21) = -1167917.436_dp
 T2_Y(22) = -1167998.331_dp

 T3_X(1)  = -662890.1396_dp
 T3_X(2)  = -662802.3219_dp
 T3_X(3)  = -662741.9164_dp
 T3_X(4)  = -662663.2705_dp
 T3_X(5)  = -662585.2473_dp
 T3_X(6)  = -662495.9469_dp
 T3_X(7)  = -662413.2646_dp
 T3_X(8)  = -662295.8688_dp
 T3_X(9)  = -662180.4386_dp
 T3_X(10) = -662048.5109_dp
 T3_X(11) = -661902.9063_dp
 T3_X(12) = -661721.7629_dp
 T3_X(13) = -661650.3932_dp
 T3_X(14) = -661587.336_dp
 T3_X(15) = -661482.3696_dp
 T3_X(16) = -661307.7577_dp
 T3_X(17) = -661245.1281_dp
 T3_X(18) = -661102.2408_dp
 T3_X(19) = -660938.88_dp
 T3_X(20) = -660856.3239_dp
 T3_X(21) = -660787.3484_dp
 T3_X(22) = -660695.8378_dp
 T3_X(23) = -660634.6415_dp
 T3_X(24) = -660552.0778_dp
 T3_X(25) = -660477.4099_dp
 T3_X(26) = -660393.7655_dp
 T3_X(27) = -660321.4734_dp
 T3_X(28) = -660247.1083_dp
 T3_X(29) = -660181.7881_dp
 T3_X(30) = -660110.912_dp
 T3_X(31) = -660053.6477_dp

 T3_Y(1)  = -1164929.342_dp
 T3_Y(2)  = -1164999.325_dp
 T3_Y(3)  = -1165053.665_dp
 T3_Y(4)  = -1165121.692_dp
 T3_Y(5)  = -1165181.713_dp
 T3_Y(6)  = -1165251.491_dp
 T3_Y(7)  = -1165317.963_dp
 T3_Y(8)  = -1165431.754_dp
 T3_Y(9)  = -1165513.441_dp
 T3_Y(10) = -1165629.556_dp
 T3_Y(11) = -1165747.89_dp
 T3_Y(12) = -1165862.081_dp
 T3_Y(13) = -1165936.809_dp
 T3_Y(14) = -1166000.112_dp
 T3_Y(15) = -1166087.718_dp
 T3_Y(16) = -1166260.711_dp
 T3_Y(17) = -1166290.969_dp
 T3_Y(18) = -1166384.098_dp
 T3_Y(19) = -1166517.847_dp
 T3_Y(20) = -1166574.052_dp
 T3_Y(21) = -1166631.936_dp
 T3_Y(22) = -1166682.742_dp
 T3_Y(23) = -1166723.112_dp
 T3_Y(24) = -1166836.567_dp
 T3_Y(25) = -1166879.031_dp
 T3_Y(26) = -1166999.411_dp
 T3_Y(27) = -1167114.243_dp
 T3_Y(28) = -1167161.9_dp
 T3_Y(29) = -1167236.917_dp
 T3_Y(30) = -1167310.288_dp
 T3_Y(31) = -1167369.831_dp

 L_X(1)  = -663968.6933_dp
 L_X(2)  = -663586.955_dp
 L_X(3)  = -663303.9103_dp
 L_X(4)  = -663131.1088_dp
 L_X(5)  = -662982.1971_dp
 L_X(6)  = -662848.1675_dp
 L_X(7)  = -662633.2881_dp
 L_X(8)  = -662449.5498_dp
 L_X(9)  = -662293.8343_dp
 L_X(10) = -662232.0748_dp
 L_X(11) = -661926.1007_dp
 L_X(12) = -661807.1384_dp
 L_X(13) = -661718.8471_dp
 L_X(14) = -661583.8865_dp
 L_X(15) = -661413.6382_dp
 L_X(16) = -661307.7577_dp
 L_X(17) = -661202.063_dp
 L_X(18) = -660996.7597_dp
 L_X(19) = -660930.9652_dp
 L_X(20) = -660841.4055_dp
 L_X(21) = -660784.2046_dp
 L_X(22) = -660705.1557_dp
 L_X(23) = -660651.6858_dp
 L_X(24) = -660500.7636_dp
 L_X(25) = -660440.2045_dp
 L_X(26) = -660361.4857_dp
 L_X(27) = -660299.0993_dp
 L_X(28) = -660236.752_dp
 L_X(29) = -660171.3885_dp
 L_X(30) = -659899.2529_dp

 L_Y(1)  = -1168363.18_dp
 L_Y(2)  = -1168129.177_dp
 L_Y(3)  = -1167984.271_dp
 L_Y(4)  = -1167878.329_dp
 L_Y(5)  = -1167768.591_dp
 L_Y(6)  = -1167816.349_dp
 L_Y(7)  = -1167683.871_dp
 L_Y(8)  = -1167533.783_dp
 L_Y(9)  = -1167485.344_dp
 L_Y(10) = -1167382.354_dp
 L_Y(11) = -1167078.56_dp
 L_Y(12) = -1166836.033_dp
 L_Y(13) = -1166711.573_dp
 L_Y(14) = -1166562.919_dp
 L_Y(15) = -1166391.943_dp
 L_Y(16) = -1166260.711_dp
 L_Y(17) = -1166176.854_dp
 L_Y(18) = -1166001.181_dp
 L_Y(19) = -1165917.307_dp
 L_Y(20) = -1165806.572_dp
 L_Y(21) = -1165746.609_dp
 L_Y(22) = -1165663.224_dp
 L_Y(23) = -1165593.379_dp
 L_Y(24) = -1165433.154_dp
 L_Y(25) = -1165370.128_dp
 L_Y(26) = -1165291.214_dp
 L_Y(27) = -1165203.102_dp
 L_Y(28) = -1165111.217_dp
 L_Y(29) = -1165035.849_dp
 L_Y(30) = -1164756.907_dp


 ALLOCATE(gatherVelo(Npes), STAT=istat )
 ALLOCATE(storedTriangle(3,3),  STAT=istat )
 IF ( istat /= 0 ) THEN
    CALL FATAL( SolverName, 'Memory allocation error' )
 ELSE
    CALL INFO(SolverName, 'Memory allocation done', level=1 )
 END IF

 ! Loop over all profiles
 DO tt=1, 4
   IF (tt == 1) THEN
       nbNodes = 17
   ELSE IF (tt == 2) THEN
       nbNodes = 22
   ELSE IF (tt == 3) THEN
       nbNodes = 31
   ELSE IF (tt == 4) THEN
       nbNodes = 30
   END IF

   IF (ParEnv % MyPe == 0) THEN
      IF(tt == 1) THEN
        OutFile = ListGetString(Solver % Values,'T1 profile file name',Found )
        IF(.NOT. Found) THEN
          CALL Fatal(SolverName, 'Could not find file name for T1')
        END IF
      ELSE IF (tt == 2) THEN
        OutFile = ListGetString(Solver % Values,'T2 profile file name',Found )
        IF(.NOT. Found) THEN
           CALL Fatal(SolverName, 'Could not find file name for T2')
        END IF
      ELSE IF (tt == 3) THEN
        OutFile = ListGetString(Solver % Values,'T3 profile file name',Found )
        IF(.NOT. Found) THEN
          CALL Fatal(SolverName, 'Could not find file name for T3')
        END IF
      ELSE IF (tt == 4) THEN
        OutFile = ListGetString(Solver % Values,'L profile file name',Found )
        IF(.NOT. Found) THEN
          CALL Fatal(SolverName, 'Could not find file name for L')
        END IF        
      END IF
      
      OPEN (12, FILE=OutFile)
   END IF
   
   ! Loop of over nodes for given profile
   DO i=1, nbNodes

      IF (tt == 1) THEN
           x = T1_X(i)
           y = T1_Y(i)
      ELSE IF (tt == 2) THEN
           x = T2_X(i)
           y = T2_Y(i)
      ELSE IF (tt == 3) THEN
           x = T3_X(i)
           y = T3_Y(i) 
      ELSE IF (tt == 4) THEN
           x = L_X(i)
           y = L_Y(i)
      END IF

      gatherVelo = 0.0_dp
      storedTriangle = 0.0_dp
      velocity = 0.0_dp

      DO t=1, Solver % NumberOfActiveElements
         Element => GetActiveElement(t)
         n = GetElementNOFNodes(Element)
         NodeIndexes => Element % NodeIndexes
         CALL GetElementNodes( ElementNodes, Element )
         CALL pointTriangle(x, y, ElementNodes % x(1),  ElementNodes % y(1), &
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
             storedTriangle(3,2), x, y,  storedTriangle(1,3),  storedTriangle(2,3) , storedTriangle(3,3), &
             velocity, 2.0_dp)
      CALL MPI_Gather(velocity,1,MPI_DOUBLE,gatherVelo,1,MPI_DOUBLE,0,MPI_COMM_WORLD,ierr)

     IF (ParEnv % MyPe == 0) then
        DO ii=1,Npes
              IF (gatherVelo(ii) /= 0.0_dp) THEN
                WRITE(12,'(e15.8,1x,e15.8,e15.8)') x, y, gatherVelo(ii)
              END IF
        END DO
     END IF

   END DO
    
   IF (ParEnv % MyPe == 0) then
      CLOSE(12)
   END IF
 
 END DO
 
 DEALLOCATE(storedTriangle)
 DEALLOCATE(gatherVelo)

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
END SUBROUTINE getVeloAtFieldProfiles
!------------------------------------------------------------------------------
