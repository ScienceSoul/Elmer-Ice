RECURSIVE SUBROUTINE InterpolatePointValue( Model,Solver,Timestep,TransientSimulation )
  USE DefUtils
! USE realloc_mod

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
  INTEGER :: istat, NoVariables, VariableNo, LocalNodes, DataChannel,&
       timearound=1, AllocationIncrement, DIM, dataread, elementnumber,&
       i, j, k, l,ll, SupportingPoints(99), NoDim(99),VariableDirections(99,3), indx, indx2, &
       DepthDOFs, HeightDOFs, counter, nx, ny
  INTEGER, POINTER :: Permutation(:),VarPerm(:),NodeIndexes(:), DepthPerm(:), HeightPerm(:)
  REAL(KIND=dp), ALLOCATABLE :: sigma(:), sigmaMesh(:)
  REAL(KIND=dp), POINTER :: Field(:), VarVal(:), PInData(:,:), DepthValues(:), HeightValues(:)
  REAL(KIND=dp) :: x, y, z, bed, thick, delta, dy, x1(151), x2(281), ya(151,281)
  LOGICAL, ALLOCATABLE:: IsToBeInterpolated(:)
  LOGICAL :: AllocationsDone=.FALSE., FirstTime=.TRUE., Found, GotVar, VariablesExist, & 
           InterpolateVariable, reinitiate, interpolateByLayers(99), isThere
  CHARACTER(LEN=MAX_NAME_LEN) :: SolverName, Name, DataName, VariableName(99) , VariableDataName(99), &
                              InterpolationMethod(99), PATH(99)
  CHARACTER(LEN=MAX_NAME_LEN), ALLOCATABLE :: layersName(:)
  INTEGER :: NbLayers(99), NbMeshLayers(99)
  TYPE(Variable_t), POINTER :: Var, DepthSol, HeightSol
  TYPE(Solver_t), POINTER :: PointerToSolver
  TYPE(ValueList_t), POINTER :: BC, Equation
  TYPE(Element_t), POINTER :: CurrentElement


  SAVE AllocationsDone, FirstTime, SolverName,&
       DIM,NoVariables, VariablesExist, LocalNodes,Permutation,&
       IsToBeInterpolated, SupportingPoints, VariableName,&
       VariableDataName
  
  IF (FirstTime) WRITE(SolverName,'(A)') 'InterpolatePointValue'
  PointerToSolver => Solver
  IF ( .NOT. ASSOCIATED( PointerToSolver ) ) THEN
    CALL FATAL(SolverName, ' No Solver Pointer associated')
  END IF

  !----------------------------------------
  ! Do these things for the first time only
  !------------------------------------------------------------------------------------------------------------------
  IF (FirstTime) THEN
     LocalNodes = Model % Mesh % NumberOfNodes
     DIM = CoordinateSystemDimension()
     FirstTime = .FALSE.
     CALL INFO(SolverName,'(Re-)Initialization started.',Level=3)

     !-------------------------------------------------------------
     ! Find out how many variables should be read and interpolated
     !-------------------------------------------------------------
     NoVariables = 0
     GotVar = .TRUE.
     DO WHILE(GotVar) 
        NoVariables = NoVariables + 1
        IF (NoVariables > 99) &
             CALL FATAL(SolverName,'Number of parameters cannot exceed 99')
        IF(NoVariables < 10) THEN
           WRITE (Name,'(A,I2)') 'Variable',NoVariables
           WRITE (DataName,'(A,I2)') 'Variable Data',NoVariables
        ELSE
           WRITE (Name,'(A,I3)') 'Variable',NoVariables
           WRITE (DataName,'(A,I3)') 'Variable Data',NoVariables
        END IF
        VariableName(NoVariables) = ListGetString( Solver % Values, TRIM(Name), GotVar)
        IF(GotVar) THEN
           WRITE(Message,'(A,A,A)') TRIM(Name),': ', VariableName(NoVariables)
           CALL INFO(SolverName,Message,Level=3)
        ELSE
           EXIT
        END IF

        VariableDataName(NoVariables) = ListGetString( Solver % Values, TRIM(DataName), Found)
        IF (Found) THEN
           WRITE(Message,'(A,A,A,A,A)') TRIM(Name),&
                ': ', TRIM(VariableName(NoVariables)), ' Data: ',&
                TRIM(VariableDataName(NoVariables))
           CALL INFO(SolverName,Message,Level=3)
        ELSE
           WRITE(Message,'(A,A,A,A)') TRIM(Name),' Data file: ',&
                VariableDataName(NoVariables), ' not found.'
           CALL FATAL(SolverName,Message)
        END IF
        SupportingPoints(NoVariables) = &
             ListGetInteger( Solver % Values, TRIM(Name) // ' Supporting Points', Found)
        IF (.NOT.Found) THEN
           WRITE(Message,'(A,A)') TRIM(Name),&
                ' Number of supporting points not found - setting to 2'
           CALL INFO(SolverName,Message,Level=3)
           SupportingPoints(NoVariables) = 2
        ELSE
           WRITE(Message,'(A,A,I6)') TRIM(Name),&
                ' Number of supporting points: ', SupportingPoints(NoVariables)
           CALL INFO(SolverName,Message,Level=4)
        END IF
        NoDIM(NoVariables) = &
             ListGetInteger( Solver % Values, TRIM(Name) // ' Dimensions', Found)
        VariableDirections(NoVariables,1:3) = 0
        IF (.NOT.Found) THEN
           NoDIM(NoVariables) = DIM
           DO i=1,DIM
              VariableDirections(NoVariables,i) = i
           END DO
        ELSE
           !WRITE(Message,'(A,A,I1)') TRIM(Name),&
           !     ' Dimensions set to ', NoDIM
           VariableDirections(NoVariables,1:NoDIM(NoVariables)) = &
                ListGetIntegerArray( Solver % Values,TRIM(Name) // ' Directions',Found)
           IF (Found) THEN
              WRITE(Message,'(A,A,AI1,I1,I1,A)')&
                   'Directions for Variable ', TRIM(Name), ': ',VariableDirections(NoVariables,1:3)
              CALL INFO(SolverName,Message,Level=4)
           ELSE
              WRITE(Message,'(A,A,A)') &
                   TRIM(Name) // ' Dimensions', ' found, but no keyword ', TRIM(Name) // ' Directions'
              CALL FATAL(SolverName,Message)
           END IF
        END IF

        interpolateByLayers(NoVariables) = GetLogical( Solver % Values, &
                   TRIM(Name) // ' Interpolate By Layers', Found )
        IF (Found) THEN
         IF (interpolateByLayers(NoVariables)) THEN
           WRITE(Message,'(A,A)') TRIM(Name), ' will be interpolated with layers'
           CALL INFO(SolverName,Message,Level=3)
           NbLayers(NoVariables) = ListGetInteger( Solver % Values, &
                TRIM(Name) // ' Number Of Layers', Found)
           IF (.NOT. Found) THEN
              CALL FATAL(SolverName,'Interpolation by layers but no number of layers specified.')
           END IF
           NbMeshLayers(NoVariables) = ListGetInteger( Solver % Values, &
                TRIM(Name) // ' Mesh Number Of Layers', Found)
           IF (.NOT. Found) THEN
              CALL FATAL(SolverName,'Interpolation by layers but no number of mesh layers specified.')
           END IF
           InterpolationMethod(NoVariables) = ListGetString( Solver % Values, &
             TRIM(Name) // ' Interpolation Method', Found)
           IF (.NOT. Found) THEN
              CALL FATAL(SolverName,'Interpolation by layers but interpolation method not found.')
           END IF
           PATH(NoVariables) = ListGetString( Solver % Values, &
             TRIM(Name) // ' PATH', Found)
           IF (.NOT. Found) THEN
              CALL FATAL(SolverName,'Interpolation by layers but PATH parameter not found.')
           END IF
         ELSE
           NbLayers(NoVariables) = 0
           NbMeshLayers(NoVariables) = 0
         END IF
        ELSE
         interpolateByLayers(NoVariables) = .FALSE.
         NbLayers(NoVariables) = 0
         NbMeshLayers(NoVariables) = 0         
       END IF
     END DO
     NoVariables = NoVariables-1

     ! --------------------------------------------
     ! Allocate space for new variables to be added
     ! and add them
     ! --------------------------------------------
     IF(NoVariables > 0) VariablesExist = .TRUE.
     IF (VariablesExist) THEN 
        ALLOCATE(Permutation(LocalNodes))
        IF (.NOT.ASSOCIATED(Permutation)) &
             CALL FATAL(Solvername,'Failed to allocate permutation vector\n')
        DO i=1,LocalNodes
           Permutation(i) = i
        END DO
        DO VariableNo = 1, NoVariables
           Var => VariableGet( Model % Variables, TRIM(VariableName(VariableNo)), .TRUE.)     
           IF(.NOT. ASSOCIATED( Var ) ) THEN               
              ALLOCATE(Field(LocalNodes))
              Field(1:LocalNodes) = 1._dp * Permutation(1:LocalNodes) * VariableNo
              CALL VariableAdd( Solver % Mesh % Variables, Solver % Mesh, PointerToSolver, &
                   TRIM(VariableName(VariableNo)), 1, Field, Permutation )          
              WRITE(Message,'(A,I2,A,A,A)') 'Variable no. ',&
                   VariableNo, ' (',TRIM(VariableName(VariableNo)),') added.'
              CALL INFO(SolverName,Message,Level=3)
              NULLIFY( Field ) 
           END IF
        END DO
     ELSE
        CALL FATAL(SolverName, 'No valid parameters found.')
     END IF
     CALL INFO(Solvername,'(Re-)Initialization done',Level=1)
  END IF ! FirstTime---------------------------------------------------------------------------------------------

  !----------------------
  ! Read in data and 
  ! interpolate variables
  !----------------------
  
  AllocationIncrement = MAX(LocalNodes/10,1)
  IF (VariablesExist) THEN
     DO VariableNo = 1,NoVariables
        WRITE(Message,'(A,I3,A,I3)') 'Processing variable ',VariableNo,'/',NoVariables
        CALL INFO(SolverName,Message,Level=3)
        NULLIFY(Var,VarVal,VarPerm)
        Var => VariableGet( Model % Variables, TRIM(VariableName(VariableNo)), .TRUE.) 
        IF (.NOT.ASSOCIATED(Var)) THEN
           WRITE(Message,'(A A)') 'Variable ', TRIM(VariableName(VariableNo)), ' not associated'
           CALL FATAL(SolverName,Message)
        END IF
        VarPerm  => Var % Perm
        VarVal => Var % Values

        IF(interpolateByLayers(VariableNo)) THEN
          
          !Get Depth and Height
          DepthSol => VariableGet( Solver % Mesh % Variables, 'Depth' )
          IF ( ASSOCIATED( DepthSol ) ) THEN
             DepthPerm    => DepthSol % Perm
             DepthValues  => DepthSol % Values
             DepthDOFs = DepthSol % DOFs
          ELSE
             CALL Info(SolverName, 'Could not find Depth field variable', Level=4)
          END IF

          HeightSol => VariableGet( Solver % Mesh % Variables, 'Height' )
          IF ( ASSOCIATED( HeightSol ) ) THEN
             HeightPerm    => HeightSol % Perm
             HeightValues  => HeightSol % Values
             HeightDOFs = HeightSol % DOFs
          ELSE
             CALL Info(SolverName, 'Could not find Height field variable', Level=4)
          END IF

          ! We also store the interpolation result for point x-y in XY(:,:)
          ALLOCATE(sigma(NbLayers(VariableNo)),sigmaMesh(NbMeshLayers(VariableNo)),&
                 layersName(NbLayers(VariableNo)), STAT=istat)
          IF ( istat /= 0 ) THEN
           CALL FATAL( SolverName, 'Memory allocation error for XY table' )
          END IF

          OPEN(16, FILE=TRIM(PATH(VariableNo)) // 'layers.dat', STATUS='OLD', ACTION='READ', IOSTAT=Istat)
          IF (ISTAT /= 0) CALL FATAL(SolverName,'Error opening layers.dat')
          dataread = 0
          DO
           dataread = dataread + 1
           100 FORMAT (A30)
           READ(16, FMT=100, IOSTAT=Istat) layersName(dataread)
           IF (Istat /= 0) THEN 
              EXIT
           END IF
          END DO
          CLOSE(16)
     
          OPEN(18, FILE=TRIM(PATH(VariableNo)) // 'sigma.dat', STATUS='OLD', ACTION='READ', IOSTAT=ISTAT)
          IF (ISTAT /= 0) CALL FATAL(SolverName,'Error opening sigma  file')
          dataread = 0
          DO
           dataread = dataread + 1
           READ (18, *, IOSTAT=Istat) sigma(dataread)
           IF (Istat /= 0) THEN 
              EXIT
           END IF
          END DO
          CLOSE(18)
         
          delta = 1.0_dp/REAL((NbMeshLayers(VariableNo)-1),dp)
          sigmaMesh(1) = 0.0_dp
          DO i=2,NbMeshLayers(VariableNo) 
             sigmaMesh(i) = sigmaMesh(i-1)+delta
             !WRITE(*,*) sigmaMesh(i)
          END DO

          DO l=1,NbMeshLayers(VariableNo)

           DO ll=1, NbLayers(VariableNo)
            IF(ABS(sigmaMesh(l)-sigma(ll)) < 1.0E-04) THEN

             IF (ASSOCIATED(PInData)) DEALLOCATE(PIndata)
             ALLOCATE(PInData(AllocationIncrement, NoDIM(VariableNo)+1), STAT=istat)
             IF (istat /= 0) &
             CALL FATAL(SolverName, 'Allocation error of initial pointer')

             OPEN(17, FILE=layersName(ll), STATUS='OLD', ACTION='READ', IOSTAT=ISTAT)
             IF (ISTAT /= 0) CALL FATAL(SolverName,'Error opening layer corresponding data file')
             dataread = 0
             DO 
              dataread = dataread + 1
              IF (dataread > SIZE(PInData,1)) PInData => reallocate2d(PInData,AllocationIncrement,0)
              READ (17 ,*, IOSTAT=Istat) PInData(dataread,1:NoDIM(VariableNo)+1)
              IF (Istat /= 0) THEN 
                  EXIT
              END IF
              !WRITE(*,*) Istat, dataread, PInData(dataread,1:NoDIM(VariableNo)+1)
             END DO
            CLOSE(17)

             dataread = dataread - 1
             WRITE(Message, '(A,I10,A,A,A,I3,A,A,A)') &
             'Data read  (',dataread,' datasets) from file ', TRIM(layersName(ll)),&
             ' for variable no. ', VariableNo,' (', TRIM(VariableName(VariableNo)),')'
             CALL INFO(SolverName,Message,Level=1)
             
             SELECT CASE (InterpolationMethod(VariableNo))
             CASE('inverse-distance')

             reinitiate = .TRUE.

             CASE('bilinear')

             ! Reorder data so that they are compatible to the input format for the routine locate
             ! Assume that the y coordinate is the fastest varying
             nx = 1
             ny = 1
             DO j=2, dataread
                IF (PInData(j,1) /= PInData(j-1,1)) THEN
                   nx = nx+1
                END IF
             END DO  
             DO j=2, dataread
                IF (PInData(j,1) /= PInData(j-1,1)) THEN
                   ny = j
                   EXIT
                END IF
             END DO  
             ny = ny-1
     
             x1(1) = PInData(1,1)
             indx2 = 2
             k = 0
             DO j=1, dataread
                k = k+1
                IF (k == ny ) THEN
                   x1(indx2) = PInData(j+1,1)
                   indx2 = indx2+1
                   k = 0
                   END IF
             END DO
             DO j=1, dataread
                IF (j == ny+1) EXIT
                   x2(j) = PInData(j,2)
             END DO
             k = 0
             indx2 = 1
             DO j=1, dataread
                k = k+1
                ya(indx2,k) = PInData(j,3)
                IF (k == ny) THEN
                   indx2 = indx2+1
                   k = 0
                END IF
             END DO
             END SELECT

             DO i=1,Solver % NumberOFActiveElements
   
               CurrentElement => GetActiveElement(i)
               IF (ParEnv % myPe .NE. CurrentElement % partIndex) CYCLE
               NodeIndexes => CurrentElement % NodeIndexes
 
               DO k=1, GetElementNOFNodes(CurrentElement)
 
                  thick = DepthValues(DepthDOFs*(DepthPerm(NodeIndexes(k))-1)+1) + &
                      HeightValues(HeightDOFs*(HeightPerm(NodeIndexes(k))-1)+1)
                  bed = Model % Nodes % z(NodeIndexes(k))-HeightValues(HeightDOFs*(HeightPerm(NodeIndexes(k))-1)+1)
                  z = (Model % Nodes % z(NodeIndexes(k))-bed)/thick

                  IF (ABS(z-sigmaMesh(l)) < 1.0E-04 ) THEN
                
                    !------------------------------------------------------
                    ! interpolate the values
                    !------------------------------------------------------
                    SELECT CASE (InterpolationMethod(VariableNo))
                    CASE('inverse-distance')

                    VarVal(VarPerm(NodeIndexes(k))) =  &
                     GetRadiallyInterpolatedValue(PInData,&
                       Model % Nodes % x(NodeIndexes(k)), Model % Nodes % y(NodeIndexes(k)), 0.0_dp, &
                       dataread,&
                       NoDIM(VariableNo),&
                       VariableDirections(VariableNo,1:3), &
                       2.0_dp, &
                       reinitiate,&
                       SupportingPoints(VariableNo),&
                       SolverName)
                    IF (j == 1) reinitiate = .FALSE. 

                    CASE('bilinear')
                    
                    VarVal(VarPerm(NodeIndexes(k))) = &
                         biliearInterpolate(x1, x2, ya, Model % Nodes % x(NodeIndexes(k)), &
                         Model % Nodes % y(NodeIndexes(k)))
                    END SELECT

                  END IF
               END DO
             END DO                                       
             EXIT
            END IF
           END DO ! Number of layers
         END DO   !Number of mesh layers
         DEALLOCATE(sigma)
         DEALLOCATE(sigmaMesh)
         DEALLOCATE(layersName)
         NULLIFY(DepthSol)
         NULLIFY(HeightSol)
         WRITE(*,*) 'Done Variable'
        ELSE
         OPEN (15, FILE=VariableDataName(VariableNo), STATUS="OLD", ACTION='READ', IOSTAT=Istat)
         IF (Istat /= 0) THEN 
           WRITE(Message,'(A A)') 'Error in opening file ', VariableDataName(VariableNo)
           CALL FATAL(SolverName,Message)
         END IF

         IF (ASSOCIATED(PInData)) DEALLOCATE(PIndata)
         ALLOCATE(PInData(AllocationIncrement, NoDIM(VariableNo)+1), STAT=istat)
         IF (istat /= 0) &
             CALL FATAL(SolverName, 'Allocation error of initial pointer')
         dataread = 0

         DO 
           dataread = dataread + 1
           IF (dataread > SIZE(PInData,1)) PInData => reallocate2d(PInData,AllocationIncrement,0)
           READ (15, *, END=10, IOSTAT=Istat) PInData(dataread,1:NoDIM(VariableNo)+1) 
           IF (Istat /= 0) THEN 
              EXIT
           END IF      
         END DO        
10       CLOSE(15)
         dataread = dataread - 1
         WRITE(Message, '(A,I10,A,A,A,I3,A,A,A)') &
             'Data read  (',dataread,' datasets) from file ', TRIM(VariableDataName(VariableNo)),&
             ' for variable no. ', VariableNo,' (',&
             TRIM(VariableName(VariableNo)),')'
         CALL INFO(SolverName,Message,Level=1)
      
        !------------------------------------------------------
        ! interpolate the values
        !------------------------------------------------------
        reinitiate = .TRUE.
        !------------------------------------------------------
        ! Loop all active elements of solver
        !------------------------------------------------------
         DO elementNumber=1,Solver % NumberOFActiveElements
           CurrentElement => GetActiveElement(elementNumber)
           IF (.NOT.ASSOCIATED(CurrentElement)) CALL FATAL(SolverName,'Element pointer not associated')
           Equation => GetEquation()
           IF ( ASSOCIATED( Equation ) ) THEN           
              InterpolateVariable = ListGetLogical( Equation, &
                   'Interpolate' // TRIM(VariableName(VariableNo)), Found)
              IF (.NOT.  Found) &
                   InterpolateVariable = .FALSE.
           ELSE
              WRITE(Message,'(A,I3,A)') 'Equation for body no ', &
                   CurrentElement % BodyId, ' not found'
              CALL FATAL(SolverName,Message)
           END IF          
           DO i=1,GetElementNOFNodes(CurrentElement)  
              VarVal(VarPerm(CurrentElement % NodeIndexes(i))) = &
                   GetRadiallyInterpolatedValue(PInData,&
                   Model % Nodes % x( CurrentElement % NodeIndexes( i ) ),&
                   Model % Nodes % y( CurrentElement % NodeIndexes( i ) ),&
                   Model % Nodes % z( CurrentElement % NodeIndexes( i ) ),&
                   dataread,&
                   NoDIM(VariableNo),&
                   VariableDirections(VariableNo,1:3), &
                   2.0_dp, &
                   reinitiate,&
                   SupportingPoints(VariableNo),&
                   SolverName)
           END DO
           IF (elementNumber==1) reinitiate = .FALSE.
         END DO ! DO elementNumber
        END IF
        WRITE(*,*) 'Done variable'
     END DO !DO VariableNo
  END IF
  WRITE(*,*) 'ALL DONE.'   

  RETURN

!20 CLOSE(15)
!  WRITE(Message,'(A A)') 'Error in reading file ', VariableDataName(VariableNo)
!  CALL FATAL(SolverName,Message)

CONTAINS
!-----------------------------------------------------------------------------------------------------------

  FUNCTION reallocate2d(inpointer, incr1, incr2)               ! reallocate REAL(:,:)
    REAL(KIND=dp), POINTER, DIMENSION(:,:) :: inpointer, reallocate2d
!    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:,:) :: temporarystorage
    INTEGER, INTENT(IN) :: incr1, incr2
    INTEGER :: nold(2), ierr

    IF(.NOT. ASSOCIATED(inpointer)) RETURN
    nold(1) = SIZE(inpointer,1)
    nold(2) = SIZE(inpointer,2)
 
    IF((incr1 < 0) .OR. (incr2 < 0)) THEN 
       CALL FATAL('reallocate2d', 'Increment cannot be smaler than zero!')
    END IF     
    ALLOCATE(reallocate2d(1:(nold(1) + incr1), 1:(nold(2) + incr2)), STAT=ierr)
    IF(ierr /= 0) STOP "re-allocation error"
    reallocate2d(1:nold(1),1:nold(2)) = inpointer(1:nold(1),1:nold(2))
    reallocate2d(nold(1): nold(1) + incr1,nold(2): nold(2) + incr2) = 0.0_dp
    IF (.NOT.ASSOCIATED(reallocate2d)) THEN 
       CALL FATAL('reallocate2d','Failed to reallocate array')
    END IF
    DEALLOCATE(inpointer)
  END FUNCTION reallocate2d

 FUNCTION reallocate1d(inpointer, incr1)               ! reallocate REAL
    REAL(KIND=dp), POINTER, DIMENSION(:) :: inpointer, reallocate1d
    INTEGER, INTENT(IN) :: incr1
    INTEGER :: nold(2), ierr

    IF(.NOT. ASSOCIATED(inpointer)) RETURN
    nold(1) = SIZE(inpointer)
    IF(incr1 < 0) THEN 
       WRITE(*,*) "Increment cannot be smaler than zero!"
       STOP
    END IF
    ALLOCATE(reallocate1d(1:(nold(1) + incr1)), STAT=ierr)
    IF(ierr /= 0) STOP "allocate error"
    reallocate1d(1:nold(1)+ incr1) = 0.0_dp
    IF (.NOT.ASSOCIATED(reallocate1d)) CALL FATAL('reallocate1d','Failed to reallocate vector')
    reallocate1d(1:nold(1)) = inpointer(1:nold(1))
    DEALLOCATE(inpointer) 
  END FUNCTION reallocate1d

 FUNCTION GetRadiallyInterpolatedValue(InData,XI,YI,ZI,&
     NoInData,DIM,Directions,exponent,reinitiate,numberOfSupportingPoints,SolverName)&
     RESULT(theInterpolatedValue)
  USE DefUtils
  REAL(KIND=dp) :: theInterpolatedValue
  REAL(KIND=dp), INTENT(IN)  :: InData(:,:), XI,YI,ZI,exponent
  INTEGER, INTENT(IN)  :: DIM, NoInData, numberOfSupportingPoints, Directions(3)
  LOGICAL, INTENT(IN) :: reinitiate
  CHARACTER(LEN=MAX_NAME_LEN), INTENT(IN) :: SolverName


  REAL(KIND=dp) :: maxdistance, difference, actualdifference, minmaxxy(2,3), &
       radius, weightsum, weight, datasum, X(3)
  REAL(KIND=dp), ALLOCATABLE ::  distanceToPoint(:)
  INTEGER :: i,j,k,usedSupportingPoints,actualpoint
  INTEGER, ALLOCATABLE :: closestPoints(:)
  LOGICAL :: isSmallerThanAny

  SAVE distanceToPoint, closestPoints, maxdistance
  
  IF (numberOfSupportingPoints < 2) &
     CALL FATAL(SolverName // TRIM('(GetRadiallyInterpolatedValue)'),'The number of supporting points must be at least 2')
  IF (exponent < 1) &
       CALL  FATAL(SolverName // TRIM('(GetRadiallyInterpolatedValue)'),'The exponent should be larger or equal to unity')

  !PRINT *, Directions(1:3)
  X(1:3) = 0.0_dp
  DO i=1,DIM
     SELECT CASE (Directions(i))
     CASE(1)
        X(i) = XI
     CASE(2)
        X(i) = YI
     CASE(3)
        X(i) = ZI
     CASE DEFAULT
        X(i) = 0.0_dp
     END SELECT
  END DO


  !--------------------
  ! (re)initialization
  !--------------------
  IF (reinitiate) THEN
     !------------
     ! allocations
     !------------
     IF (ALLOCATED(distanceToPoint)) DEALLOCATE(distanceToPoint)
     IF (ALLOCATED(closestPoints)) DEALLOCATE(closestPoints)
     ALLOCATE(distanceToPoint(numberOfSupportingPoints),closestPoints(numberOfSupportingPoints))
     !-----------------
     ! get bounding box
     !-----------------
     minmaxxy(1,1:DIM) = InData(1,1:DIM)
     minmaxxy(2,1:DIM) = InData(1,1:DIM)
     DO i=1,NoInData
        DO j=1,DIM
           IF (InData(i,j) <  minmaxxy(1,j)) minmaxxy(1,j) = InData(i,j)
           IF (InData(i,j) >  minmaxxy(2,j)) minmaxxy(2,j) = InData(i,j)
        END DO
     END DO
     maxdistance = 0.0_dp
     DO j=1,DIM
        maxdistance = maxdistance + (minmaxxy(2,j) - minmaxxy(1,j))**2.0_dp
     END DO
     maxdistance = 2.0_dp * sqrt(maxdistance)     
  END IF
  !-------------------

  !-------------------
  ! get the supporting
  ! points for the 
  ! interpolation
  !-------------------
  distanceToPoint = maxdistance
  closestPoints = 0  
  DO i=1,NoInData
     !-------------------------------
     ! get radius to input data point
     !-------------------------------
     radius = 0.0_dp
     DO j=1,DIM
        radius = radius + (X(j) -  InData(i,j))**2.0_dp
     END DO
     IF (radius > 1.0D-06) THEN 
        radius = SQRT(radius)
     ELSE
        radius = 0.0_dp
     END IF
     !----------------------------------------
     ! check whether one and if so which one 
     ! of the current supporting
     ! points has a longer distance
     !----------------------------------------
     actualdifference = distanceToPoint(numberOfSupportingPoints)
     actualpoint = 0
     isSmallerThanAny = .FALSE.
     DO j =1,numberOfSupportingPoints,1        
 !       PRINT *,'J=',j,' R=', radius,' D=',distanceToPoint(j)
        difference = distanceToPoint(j) - radius
        IF (difference > 0.0_dp) THEN
           isSmallerThanAny = .TRUE.
           IF (difference < actualdifference) THEN
              actualpoint = j
              actualdifference = difference
           END IF
        END IF
     END DO
     !-------------
     ! if so, swap
     ! and reorder
     !-------------
     IF (isSmallerThanAny) THEN
        DO k=numberOfSupportingPoints-1,actualpoint,-1
           distanceToPoint(k+1) =  distanceToPoint(k)
           closestPoints(k+1) = closestPoints(k)
        END DO
        distanceToPoint(actualpoint) = radius
        closestPoints(actualpoint) = i
     END IF
!     IF (ANY(closestPoints(1:numberOfSupportingPoints) == 0)) THEN        
!        CALL WARN(SolverName,'Less than the reqested supporting points found')
!     END IF
        
!     PRINT *,'N=',closestPoints(1:numberOfSupportingPoints)
!     PRINT *,'D=',distanceToPoint(1:numberOfSupportingPoints)
    
  END DO  
  !-----------------


  !--------------------------
  ! do we have a bull's eye?
  !--------------------------
  IF (distanceToPoint(1) < 1.0D-12) THEN
     theInterpolatedValue = InData(closestPoints(1),DIM+1)
  ELSE 
     !-----------
     ! inerpolate
     !-----------
     weightsum = 0.0_dp
     theInterpolatedValue = 0.0_dp
     usedSupportingPoints = 0
     DO k=1,numberOfSupportingPoints
        IF (closestPoints(k) /= 0) THEN
           usedSupportingPoints = usedSupportingPoints + 1
           weight = (distanceToPoint(k))**(-exponent)
!           PRINT *,'w=',weight,distanceToPoint(k),InData(closestPoints(k),DIM+1) 
           theInterpolatedValue = theInterpolatedValue + weight * InData(closestPoints(k),DIM+1)
           weightsum = weightsum + weight
        END IF
     END DO
!     PRINT *,'--------'
     IF (usedSupportingPoints < numberOfSupportingPoints) THEN
        WRITE(Message,'(A,F10.3,F10.3,F10.3,A,I3,A,I3)')&
             'Number of supporting points used for point (',&
             XI,YI,ZI,') =', usedSupportingPoints,&
             ' smaller than requested ', numberOfSupportingPoints
        CALL WARN(TRIM(Solvername) // '(GetRadiallyInterpolatedValue)',&
             Message)
     END IF
     IF (usedSupportingPoints == 0) THEN
        WRITE(Message,'(A,F10.3,F10.3,F10.3,A)') &
             'No supporting point for point (',&
                XI,YI,ZI,') found'
        CALL FATAL(TRIM(Solvername) // '(GetRadiallyInterpolatedValue)',&
             Message)
     END IF
     theInterpolatedValue = theInterpolatedValue/weightsum
  END IF
  RETURN

END FUNCTION GetRadiallyInterpolatedValue

FUNCTION biliearInterpolate(x1, x2, ya, x, y) RESULT(theInterpolatedValue)
USE DefUtils
IMPLICIT NONE

INTEGER :: i, j, k
REAL(KIND=dp), DIMENSION(:), INTENT(IN):: x1
REAL(KIND=dp), DIMENSION(:), INTENT(IN) :: x2
REAL(KIND=dp), DIMENSION(:,:), INTENT(IN) :: ya
REAL(KIND=dp), INTENT(IN) ::  x, y
REAL(KIND=dp) :: y1, y2, y3, y4, t, u
REAL(KIND=dp) :: theInterpolatedValue

j = locate(x1, x)
k = locate(x2, y)

y1 = ya(j,k)
y2 = ya(j+1,k)
y3 = ya(j+1,k+1)
y4 = ya(j,k+1)

t = (x-x1(j))/(x1(j+1)-x1(j))
u = (y-x2(k))/(x2(k+1)-x2(k))

theInterpolatedValue = (1.0_dp-t)*(1.0_dp-u)*y1+t*(1.0_dp-u)*y2+t*u*y3+(1.0_dp-t)*u*y4

END FUNCTION biliearInterpolate

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   Function locate
!   Given an array xx(1:N), and given a value x, returns a value j such that x is netween xx(j) and xx(j+1). 
!   xx muste be monotonic , either increasing or decreasing. j=0 orj=N is returned to indicate that x is out of range.
!   Courtesy of Numerical Recipies in Fortran 90.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
FUNCTION locate(xx,x)
USE DefUtils
IMPLICIT NONE

REAl(KIND=dp), DIMENSION(:), INTENT(IN) :: xx
REAL(KIND=dp), INTENT(IN) :: x
INTEGER :: locate
INTEGER :: n, jl, jm, ju
LOGICAL :: ascnd

n = size(xx)
ascnd = (xx(n) >= xx(1))
jl=0
ju=n+1
DO 
  IF (ju-jl <= 1) EXIT
  jm=(ju+jl)/2
  IF (ascnd .eqv. (x >= xx(jm))) THEN
    jl=jm
  ELSE
    ju=jm
  END IF
  END DO
IF (x == xx(1)) THEN
   locate=1
ELSE IF (x == xx(n)) THEN
   locate=n-1
ELSE 
   locate=jl
END IF
END FUNCTION locate 

FUNCTION BOOL(x,y) RESULT(TRUE)
REAL(KIND=dp), INTENT(IN) :: x, y
LOGICAL:: TRUE

IF (ABS(x-y) < 1.0E-04) THEN
  TRUE=.TRUE.
ELSE
  TRUE=.FALSE.
END IF
END FUNCTION BOOL

END SUBROUTINE InterpolatePointValue


! sets variable entries to variable of this solver
!--------------------------------------------------
RECURSIVE SUBROUTINE PointwiseVariableSet( Model,Solver,Timestep,TransientSimulation )
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
  REAL(KIND=dp), ALLOCATABLE :: InData(:)

  TYPE(Variable_t), POINTER :: PointerVar

  INTEGER :: i,j

  PointerVar => Solver % Variable

  IF ( ASSOCIATED(PointerVar) ) THEN
     ! sets the array for data I/O to the size of your mesh points
     ! i.e., we are dealing with a scalar variable
     ! modify, if vector/tensor parameters are needed
     ALLOCATE(InData(Model % Mesh % NumberOfNodes))  
     
     ! PETE: here you have to write yourself the I/O part, that fills InData with the point-wise information

     DO i=1,SIZE(PointerVar % Perm)
        j = PointerVar % Perm(i)
        IF ( j>0 ) THEN
           PointerVar % Values(j) = InData(i)
        END IF
     END DO
     DEALLOCATE(InData)
  ELSE 
     WRITE(Message,'(A,A)')&
          'Could not find variable',TRIM(Solver % Variable % Name)
     CALL FATAL("PointwiseVariableSet",Message);
  END IF
END SUBROUTINE PointwiseVariableSet
