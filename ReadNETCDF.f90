RECURSIVE SUBROUTINE ReadNETCDF( Model,Solver,Timestep,TransientSimulation )
    USE netcdf
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
 INTEGER, POINTER :: Permutation(:),VarPerm(:),NodeIndexes(:)
 REAL(KIND=dp), POINTER :: Field(:), VarVal(:)
 REAL(KIND=dp), ALLOCATABLE :: PInData(:,:)
 TYPE(Variable_t), POINTER :: Var
 TYPE(Solver_t), POINTER :: PointerToSolver
 TYPE(ValueList_t), POINTER :: BC, Equation
 TYPE(Element_t), POINTER :: CurrentElement
 INTEGER :: istat, NoVariables, VariableNo, LocalNodes, DataChannel, &
           timearound=1, AllocationIncrement, DIM, dataread, elementnumber,&
           i, j, k, SupportingPoints(99), NoDim(99),NetCDFDIM(99), VariableDirections(99,3), &
           AssocVariableDirections(99,3), VariableSizeInDirections(99,3), incrTime(99)
 CHARACTER(LEN=MAX_NAME_LEN) :: SolverName, Name, DataName, VariableName(99) , VariableDataName(99), &
                             AssociatedVariableName(99), XName(99), YName(99)
 LOGICAL :: AllocationsDone=.FALSE., FirstTime=.TRUE., Found, GotVar, VariablesExist, & 
          InterpolateVariable, reinitiate, IsInteger(99)
 !------------------------------------------------------------------------------
 !     Sone stuff needed for NETCDF reading
 !------------------------------------------------------------------------------
 INTEGER :: N1, N2, N3
 REAL(KIND=dp), ALLOCATABLE :: ddata_in2D(:,:), ddata_in3D(:,:,:), Xcoord(:), Ycoord(:), Zcoord(:)
 INTEGER, ALLOCATABLE :: idata_in2D(:,:), idata_in3D(:,:,:)
 INTEGER :: ncid, varid
 INTEGER :: x, y, z

 SAVE   AllocationsDone, FirstTime, SolverName,&
       DIM,NoVariables, VariablesExist, LocalNodes,Permutation,&
       SupportingPoints, VariableName, VariableDirections, AssocVariableDirections, &
       VariableDataName, VariableSizeInDirections, AssociatedVariableName, XName, YName, &
       NoDim, NetCDFDIM, IsInteger, incrTime
 
 IF (FirstTime) WRITE(SolverName,'(A)') 'ReadNETCDF'
 PointerToSolver => Solver
 IF ( .NOT. ASSOCIATED( PointerToSolver ) ) THEN
    CALL FATAL(SolverName, ' No Solver Pointer associated')
 END IF

 !------------------------------------------------------------------
 ! Do these things for the first time only
 !------------------------------------------------------------------
 IF (FirstTime) THEN
    LocalNodes = Model % NumberOfNodes
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
           WRITE (Name,'(A,I2)') 'NETCDF Variable',NoVariables
           WRITE (DataName,'(A,I2)') 'NETCDF Variable Data',NoVariables
       ELSE
           WRITE (Name,'(A,I3)') 'NETCDF Variable',NoVariables
           WRITE (DataName,'(A,I3)') 'NETCDF Variable Data',NoVariables
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

       AssociatedVariableName(NoVariables) = ListGetString( Solver % Values, & 
                   TRIM(Name) // ' Associated Variable', Found)
       IF (Found) THEN
           WRITE(Message,'(A,A,A,A,A)') TRIM(Name),&
                ': ', TRIM(VariableName(NoVariables)), ' Associated Variable: ',&
                TRIM(AssociatedVariableName(NoVariables))
           CALL INFO(SolverName,Message,Level=3)
       ELSE
           WRITE(Message,'(A,A,A,A)') TRIM(Name),' Associated Variable: ',&
                VariableDataName(NoVariables), ' not found.'
           CALL FATAL(SolverName,Message)
       END IF  

       NetCDFDIM(NoVariables) = &
             ListGetInteger( Solver % Values, TRIM(Name) // ' Dimensions', Found)
       VariableDirections(NoVariables,1:3) = 0      
       IF (.NOT.Found) THEN
           NetCDFDIM(NoVariables) = DIM
           DO i=1,DIM
              VariableDirections(NoVariables,i) = i
           END DO
       ELSE
           VariableDirections(NoVariables,1:NetCDFDIM(NoVariables)) = &
                ListGetIntegerArray( Solver % Values,TRIM(Name) // ' Directions',Found)
           IF (Found) THEN
              WRITE(Message,'(A,A,A,I1,I1,I1)')&
                   'Directions for Variable ', TRIM(Name), ': ',VariableDirections(NoVariables,1:3)
              CALL INFO(SolverName,Message,Level=4)
           ELSE
              WRITE(Message,'(A,A,A)') &
                   TRIM(Name) // ' Dimensions', ' found, but no keyword ', TRIM(Name) // ' Directions'
              CALL FATAL(SolverName,Message)
           END IF
       END IF

       IsInteger(NoVariables) = GetLogical( Solver % Values , TRIM(Name) // ' Integer Type', Found )
       IF (Found) THEN
         IF (IsInteger(NoVariables)) THEN
           WRITE(Message,'(A,A)') TRIM(Name), ' is of type integer'
           CALL INFO(SolverName,Message,Level=3)
         ELSE 
           WRITE(Message,'(A,A)') TRIM(Name), ' is of type float'
           CALL INFO(SolverName,Message,Level=3)
         END IF
       ELSE
           IsInteger(NoVariables) = .FALSE.
           WRITE(Message,'(A,A)') TRIM(Name), ' is of type float'
           CALL INFO(SolverName,Message,Level=3)
       END IF

       VariableSizeInDirections(NoVariables,1:NetCDFDIM(NoVariables)) = &
                ListGetIntegerArray( Solver % Values,TRIM(Name) // ' size in directions',Found)
       IF (Found) THEN
          WRITE(Message,'(A,A,I6,I6,I6)') TRIM(Name), ': ',VariableSizeInDirections(NoVariables,1:3)
          CALL INFO(SolverName,Message,Level=4)
       ELSE
          WRITE(Message,'(A,A)') &
                   TRIM(Name) // ' size in directions', ' not found'
          CALL FATAL(SolverName,Message)
       END IF

       XName(NoVariables) = ListGetString( Solver % Values, TRIM(Name) // ' Coordinate 1 name', Found)
       IF (Found) THEN
           WRITE(Message,'(A,A,A,A,A)') TRIM(Name),&
                ': ', TRIM(VariableName(NoVariables)), ' Coordinate 1 name: ', XName(NoVariables)
           CALL INFO(SolverName,Message,Level=3)
       ELSE
           WRITE(Message,'(A,A,A,A,A)') TRIM(Name),&
                ': ', TRIM(VariableName(NoVariables)), ' Coordinate 1 name: ', 'not found'
           CALL FATAL(SolverName,Message)
       END IF

       YName(NoVariables) = ListGetString( Solver % Values, TRIM(Name) // ' Coordinate 2 name', Found)
       IF (Found) THEN
           WRITE(Message,'(A,A,A,A,A)') TRIM(Name),&
                ': ', TRIM(VariableName(NoVariables)), ' Coordinate 2 name: ', YName(NoVariables)
           CALL INFO(SolverName,Message,Level=3)
       ELSE
           WRITE(Message,'(A,A,A,A,A)') TRIM(Name),&
                ': ', TRIM(VariableName(NoVariables)), ' Coordinate 2 name: ', 'not found'
           CALL FATAL(SolverName,Message)
       END IF
   
       SupportingPoints(NoVariables) = &
             ListGetInteger( Solver % Values, TRIM(Name) // ' Associated Variable Supporting Points', Found)
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
             ListGetInteger( Solver % Values, TRIM(Name) // ' Associated Variable Dimensions', Found)
       AssocVariableDirections(NoVariables,1:3) = 0
       IF (.NOT.Found) THEN
           NoDIM(NoVariables) = DIM
           DO i=1,DIM
              AssocVariableDirections(NoVariables,i) = i
           END DO
       ELSE
           !WRITE(Message,'(A,A,I1)') TRIM(Name),&
           !     ' Dimensions set to ', NoDIM
           AssocVariableDirections(NoVariables,1:NoDIM(NoVariables)) = &
                ListGetIntegerArray( Solver % Values,TRIM(Name) // ' Associated Variable Directions',Found)
           IF (Found) THEN
              WRITE(Message,'(A,A,A,I1,I1,I1)')&
                   'Directions for Variable ', TRIM(Name), ': ',AssocVariableDirections(NoVariables,1:3)
              CALL INFO(SolverName,Message,Level=4)
           ELSE
              WRITE(Message,'(A,A,A)') &
                   TRIM(Name) // ' Associaced Variable Dimensions', ' found, but no keyword ', TRIM(Name) // &
                                ' Associated Variable Directions'
              CALL FATAL(SolverName,Message)
           END IF
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
           Var => VariableGet( Model % Variables, TRIM(AssociatedVariableName(VariableNo)), .TRUE.)     
           IF(.NOT. ASSOCIATED( Var ) ) THEN               
              ALLOCATE(Field(LocalNodes))
              Field(1:LocalNodes) = 1._dp * Permutation(1:LocalNodes) * VariableNo
              CALL VariableAdd( Solver % Mesh % Variables, Solver % Mesh, PointerToSolver, &
                   TRIM(AssociatedVariableName(VariableNo)), 1, Field, Permutation )          
              WRITE(Message,'(A,I2,A,A,A)') 'Variable no. ',&
                   VariableNo, ' (',TRIM(AssociatedVariableName(VariableNo)),') added.'
              CALL INFO(SolverName,Message,Level=3)
              NULLIFY( Field ) 
           END IF
        END DO
     ELSE
        CALL FATAL(SolverName, 'No valid parameters found.')
     END IF
     CALL INFO(Solvername,'(Re-)Initialization done',Level=1)

     !Starts at Time indx=1 for all variables
     incrTime = 1
 END IF ! FirstTime---------------------------------------------------------------------------------------------


 !-------------------------
 ! Read in NETCDF data and 
 ! interpolate variables
 !------------------------
 AllocationIncrement = MAX(LocalNodes/10,1)
 IF (VariablesExist) THEN
    DO VariableNo = 1,NoVariables
       WRITE(Message,'(A,I3,A,I3)') 'Processing variable ',VariableNo,'/',NoVariables
       CALL INFO(SolverName,Message,Level=3)
       NULLIFY(Var,VarVal,VarPerm)
       Var => VariableGet( Model % Variables, TRIM(AssociatedVariableName(VariableNo)), .TRUE.) 
       IF (.NOT.ASSOCIATED(Var)) THEN
           WRITE(Message,'(A,A)') 'Variable ', TRIM(AssociatedVariableName(VariableNo)), ' not associated'
           CALL FATAL(SolverName,Message)
       END IF
       VarPerm  => Var % Perm
       VarVal => Var % Values

       IF (NetCDFDIM(VariableNo) == 2) THEN
           N1 = VariableSizeInDirections(VariableNo,1)
           N2 = VariableSizeInDirections(VariableNo,2)
           N3 = 0
       ELSE IF (NetCDFDIM(VariableNo) == 3) THEN
           N1 = VariableSizeInDirections(VariableNo,1)
           N2 = VariableSizeInDirections(VariableNo,2)
           N3 = VariableSizeInDirections(VariableNo,3)
       ELSE 
          CALL FATAL(SolverName, 'Unsupported dimension for input NETCDF variable.')
       END IF

       IF (IsInteger(VariableNo)) THEN
          IF (NetCDFDIM(VariableNo) == 2) THEN         
             ALLOCATE(idata_in2D(N1,N2),STAT=istat) 
             IF ( istat /= 0 ) THEN
                CALL FATAL( SolverName, 'Memory allocation error' )
             ELSE
                CALL INFO(SolverName, 'Memory allocation done', level=1 )
             END IF       
          ELSE IF (NetCDFDIM(VariableNo) == 3) THEN
                 ALLOCATE(idata_in3D(N1,N2,N3),STAT=istat)
                 IF ( istat /= 0 ) THEN
                    CALL FATAL( SolverName, 'Memory allocation error' )
                 ELSE
                    CALL INFO(SolverName, 'Memory allocation done', level=1 )
                 END IF
          END IF
       ELSE
          IF (NetCDFDIM(VariableNo) == 2) THEN
             ALLOCATE(ddata_in2D(N1,N2),STAT=istat)
                IF ( istat /= 0 ) THEN
                    CALL FATAL( SolverName, 'Memory allocation error' )
                ELSE
                    CALL INFO(SolverName, 'Memory allocation done', level=1 )
                END IF            
          ELSE IF (NetCDFDIM(VariableNo) == 3) THEN
                 ALLOCATE(ddata_in3D(N1,N2,N3),STAT=istat)
                 IF ( istat /= 0 ) THEN
                    CALL FATAL( SolverName, 'Memory allocation error' )
                 ELSE
                    CALL INFO(SolverName, 'Memory allocation done', level=1 )
                 END IF
          END IF
       END IF
       
       ! Open the file. NF90_NOWRITE tells netCDF we want read-only access to the file.
       istat  = nf90_open(TRIM(VariableDataName(VariableNo)), NF90_NOWRITE, ncid)
       if (istat /= nf90_noerr) CALL FATAL(SolverName, 'No NETCDF input file  found.')
       ! Get the varid of the data variable, based on its name.
       call check( nf90_inq_varid(ncid, TRIM(VariableName(VariableNo)), varid) )
       ! Read the data.
       IF (IsInteger(VariableNo)) THEN
        IF (NetCDFDIM(VariableNo) == 2) THEN
           call check( nf90_get_var(ncid, varid, idata_in2D) )
        ELSE IF  (NetCDFDIM(VariableNo) == 3) THEN
           call check( nf90_get_var(ncid, varid, idata_in3D) )
        END IF
       ELSE
        IF (NetCDFDIM(VariableNo) == 2) THEN
           call check( nf90_get_var(ncid, varid, ddata_in2D) )
        ELSE IF  (NetCDFDIM(VariableNo) == 3) THEN
           call check( nf90_get_var(ncid, varid, ddata_in3D) )
        END IF
       END IF

       IF (NoDIM(VariableNo) == 2) THEN
          If (VariableDirections(VariableNo,1) == 1) THEN
             ALLOCATE(Xcoord(VariableSizeInDirections(VariableNo,1)), & 
                    Ycoord(VariableSizeInDirections(VariableNo,2)), STAT=istat)
             IF ( istat /= 0 ) THEN
                CALL FATAL( SolverName, 'Memory allocation error' )
             ELSE
                CALL INFO(SolverName, 'Memory allocation done', level=1 )
             END IF
          ELSE IF (VariableDirections(VariableNo,1) == 2) THEN
             ALLOCATE(Xcoord(VariableSizeInDirections(VariableNo,2)), & 
                    Ycoord(VariableSizeInDirections(VariableNo,1)), STAT=istat)
             IF ( istat /= 0 ) THEN
                CALL FATAL( SolverName, 'Memory allocation error' )
             ELSE
                CALL INFO(SolverName, 'Memory allocation done', level=1 )
             END IF       
          END IF
       ELSE
          CALL FATAL(SolverName, 'Interpolation in 3D space (x,y,z) variable not implemented yet .')
       END IF

       ! Get the varid of the data variable for the X coordinate and read the data
       call check( nf90_inq_varid(ncid, XName(VariableNo), varid) )   
       call check( nf90_get_var(ncid, varid, Xcoord) )
       ! Get the varid of the data variable for the Y coordinate and read the data
       call check( nf90_inq_varid(ncid, YName(VariableNo), varid) )   
       call check( nf90_get_var(ncid, varid, Ycoord) )

       ! Close the file, freeing all resources.
       call check( nf90_close(ncid) )
       WRITE(Message,'(A,A,A)') '*** SUCCESS reading netcdf file', VariableDataName(VariableNo), '!'
       CALL INFO(SolverName,Message,Level=3)

       !Reorganise the netcdf data, x-y-variableToInterpolate in one single  table that will passed to the interpoler routine. 
       !Data are organized as x-y-variableToInterpolate
       dataread=0
       IF (NetCDFDIM(VariableNo) == 2) THEN
          DO i=1, N1
            DO j=1,N2
              dataread = dataread+1
            END DO
          END DO
       ELSE IF (NetCDFDIM(VariableNo) == 3) THEN

          IF (NetCDFDIM(VariableNo) == 3 .AND. NoDIM(VariableNo) == 2) THEN

             DO i=1, N1
              DO j=1,N2
                dataread = dataread+1
              END DO
             END DO

          ELSE IF (NetCDFDIM(VariableNo) == 3 .AND. NoDIM(VariableNo) == 3) THEN

             DO i=1, N1
              DO j=1,N2
                DO k=1,N3
                  dataread = dataread+1
                END DO
              END DO
             END DO

          ELSE

             CALL FATAL(SolverName, 'Unsupported dimensions in computation of size of interpolants table.')

          END IF
      
       ELSE
 
          CALL FATAL(SolverName, 'NetCDF data with dimension higher than three are not supported.')   

       END IF

       Allocate(PInData(dataread,DIM+1),STAT=istat)
       IF ( istat /= 0 ) THEN
          CALL FATAL( SolverName, 'Memory allocation error' )
       ELSE
          CALL INFO(SolverName, 'Memory allocation done', level=1 )
       END IF  
       dataread=1
       !Check associated variable dimension. If 2D and netcdf variable is 3D, only copy two dimensions from netcdf variable.
       !If 3D and netcdf variable is 3D, copy all dimensions
       !Check directions in netcdf variable then copy to PInData
       IF (NetCDFDIM(VariableNo) == 2 .AND. NoDIM(VariableNo) == 2) THEN
           If (VariableDirections(VariableNo,1) == 1) THEN
              DO i=1,N2
                DO j=1,N1
                   PInData(dataread, 1) = Xcoord(j)
                   PInData(dataread, 2) = Ycoord(i)
                   PInData(dataread, 3) = ddata_in2D(j,i)
                   dataread = dataread +1
                END DO
              END DO
            Else IF (VariableDirections(VariableNo,1) == 2) THEN
              DO i=1,N2
                DO j=1,N1
                   PInData(dataread, 1) = Ycoord(j)
                   PInData(dataread, 2) = Xcoord(i)
                   PInData(dataread, 3) = ddata_in2D(j,i)
                   dataread = dataread +1
                END DO
              END DO
           END IF
        ELSE IF (NetCDFDIM(VariableNo) == 3 .AND. NoDIM(VariableNo) == 2) THEN
           IF (VariableDirections(VariableNo,1) == 1 .AND. VariableDirections(VariableNo,2) == 2 .AND. &
                 VariableDirections(VariableNo,3) == 3 ) THEN
              IF (VariableDirections(VariableNo,3) > 1) THEN
               DO j=1,N2
                  DO k=1,N1
                     PInData(dataread, 1) = Xcoord(k)
                     PInData(dataread, 2) = Ycoord(j)
                     PInData(dataread, 3) = ddata_in3D(k,j,incrTime(VariableNo))
                     dataread = dataread+1
                  END DO
               END DO
              ELSE
               DO i=1,N3
                DO j=1,N2
                  DO k=1,N1
                     PInData(dataread, 1) = Xcoord(k)
                     PInData(dataread, 2) = Ycoord(j)
                     PInData(dataread, 3) = ddata_in3D(k,j,i)
                     dataread = dataread+1
                  END DO
                END DO
               END DO
              END IF
           ELSE IF (VariableDirections(VariableNo,1) == 1 .AND. VariableDirections(VariableNo,2) == 3 .AND. &
                 VariableDirections(VariableNo,3) == 2 ) THEN  
              IF (VariableDirections(VariableNo,2) > 1) THEN
               DO i=1,N3
                  DO k=1,N1
                     PInData(dataread, 1) = Xcoord(k)
                     PInData(dataread, 2) = Ycoord(i)
                     PInData(dataread, 3) = ddata_in3D(k,incrTime(VariableNo),i)
                     dataread = dataread+1
                  END DO
               END DO
              ELSE
               DO i=1,N3
                DO j=1,N2
                  DO k=1,N1
                     PInData(dataread, 1) = Xcoord(k)
                     PInData(dataread, 2) = Ycoord(i)
                     PInData(dataread, 3) = ddata_in3D(k,j,i)
                     dataread = dataread+1
                  END DO
                END DO
               END DO
              END IF
           ELSE IF (VariableDirections(VariableNo,1) == 2 .AND. VariableDirections(VariableNo,2) == 1 .AND. &
                 VariableDirections(VariableNo,3) == 3 ) THEN
              IF (VariableDirections(VariableNo,3) > 1) THEN
               DO j=1,N2
                  DO k=1,N1
                     PInData(dataread, 1) = Ycoord(k)
                     PInData(dataread, 2) = Xcoord(j)
                     PInData(dataread, 3) = ddata_in3D(k,j,incrTime(VariableNo))
                     dataread = dataread+1
                  END DO
               END DO
              ELSE
               DO i=1,N3
                DO j=1,N2
                  DO k=1,N1
                     PInData(dataread, 1) = Ycoord(k)
                     PInData(dataread, 2) = Xcoord(j)
                     PInData(dataread, 3) = ddata_in3D(k,j,i)
                     dataread = dataread+1
                  END DO
                END DO
               END DO
              END IF
           ELSE IF (VariableDirections(VariableNo,1) == 2 .AND. VariableDirections(VariableNo,2) == 3 .AND. &
                 VariableDirections(VariableNo,3) == 1 ) THEN
              IF (VariableDirections(VariableNo,2) > 1) THEN
               DO i=1,N3
                  DO k=1,N1
                     PInData(dataread, 1) = Ycoord(k)
                     PInData(dataread, 2) = Xcoord(i)
                     PInData(dataread, 3) = ddata_in3D(k,incrTime(VariableNo),i)
                     dataread = dataread+1
                  END DO
               END DO
              ELSE
               DO i=1,N3
                DO j=1,N2
                  DO k=1,N1
                     PInData(dataread, 1) = Ycoord(k)
                     PInData(dataread, 2) = Xcoord(i)
                     PInData(dataread, 3) = ddata_in3D(k,j,i)
                     dataread = dataread+1
                  END DO
                END DO
               END DO
              END IF
           ELSE IF (VariableDirections(VariableNo,1) == 3 .AND. VariableDirections(VariableNo,2) == 1 .AND. &
                 VariableDirections(VariableNo,3) == 2 ) THEN
              IF (VariableDirections(VariableNo,1) > 1) THEN
               DO i=1,N3
                  DO j=1,N2   
                     PInData(dataread, 1) = Xcoord(j)
                     PInData(dataread, 2) = Ycoord(i)
                     PInData(dataread, 3) = ddata_in3D(incrTime(VariableNo),j,i)
                     dataread = dataread+1
                  END DO
                END DO
              ELSE
               DO i=1,N3
                DO j=1,N2
                  DO k=1,N1
                     PInData(dataread, 1) = Xcoord(j)
                     PInData(dataread, 2) = Ycoord(i)
                     PInData(dataread, 3) = ddata_in3D(k,j,i)
                     dataread = dataread+1
                  END DO
                END DO
               END DO
              END IF
           ELSE IF (VariableDirections(VariableNo,1) == 3 .AND. VariableDirections(VariableNo,2) == 2 .AND. &
                 VariableDirections(VariableNo,3) == 1 ) THEN
              IF (VariableDirections(VariableNo,1) > 1) THEN
               DO i=1,N3
                  DO j=1,N2
                     PInData(dataread, 1) = Ycoord(j)
                     PInData(dataread, 2) = Xcoord(i)
                     PInData(dataread, 3) = ddata_in3D(incrTime(VariableNo),j,i)
                     dataread = dataread+1
                  END DO
               END DO
              ELSE
               DO i=1,N3
                DO j=1,N2
                  DO k=1,N1
                     PInData(dataread, 1) = Ycoord(j)
                     PInData(dataread, 2) = Xcoord(i)
                     PInData(dataread, 3) = ddata_in3D(k,j,i)
                     dataread = dataread+1
                  END DO
                END DO
               END DO
              END IF
           END IF
        ELSE IF (NetCDFDIM(VariableNo) == 3 .AND. NoDIM(VariableNo) == 3) THEN
           IF (VariableDirections(VariableNo,1) == 1 .AND. VariableDirections(VariableNo,2) == 2 .AND. &
                 VariableDirections(VariableNo,3) == 3 ) THEN
              DO i=1,N3
                DO j=1,N2
                  DO k=1,N1
                     PInData(dataread, 1) = Xcoord(k)
                     PInData(dataread, 2) = Ycoord(j)
                     PInData(dataread, 3) = Zcoord(i)
                     PInData(dataread, 4) = ddata_in3D(k,j,i)
                     dataread = dataread+1
                  END DO
                END DO
              END DO
           ELSE IF (VariableDirections(VariableNo,1) == 1 .AND. VariableDirections(VariableNo,2) == 3 .AND. &
                 VariableDirections(VariableNo,3) == 2 ) THEN  
              DO i=1,N3
                DO j=1,N2
                  DO k=1,N1
                     PInData(dataread, 1) = Xcoord(k)
                     PInData(dataread, 2) = Zcoord(j)
                     PInData(dataread, 3) = Ycoord(i)
                     PInData(dataread, 4) = ddata_in3D(k,j,i)
                     dataread = dataread+1
                  END DO
                END DO
              END DO
           ELSE IF (VariableDirections(VariableNo,1) == 2 .AND. VariableDirections(VariableNo,2) == 1 .AND. &
                 VariableDirections(VariableNo,3) == 3 ) THEN  
              DO i=1,N3
                DO j=1,N2
                  DO k=1,N1
                     PInData(dataread, 1) = Ycoord(k)
                     PInData(dataread, 2) = Xcoord(j)
                     PInData(dataread, 3) = Zcoord(i)
                     PInData(dataread, 4) = ddata_in3D(k,j,i)
                     dataread = dataread+1
                  END DO
                END DO
              END DO
           ELSE IF (VariableDirections(VariableNo,1) == 2 .AND. VariableDirections(VariableNo,2) == 3 .AND. &
                 VariableDirections(VariableNo,3) == 1 ) THEN
              DO i=1,N3
                DO j=1,N2
                  DO k=1,N1
                     PInData(dataread, 1) = Ycoord(k)
                     PInData(dataread, 2) = Zcoord(j)
                     PInData(dataread, 3) = Xcoord(i)
                     PInData(dataread, 4) = ddata_in3D(k,j,i)
                     dataread = dataread+1
                  END DO
                END DO
              END DO
           ELSE IF (VariableDirections(VariableNo,1) == 3 .AND. VariableDirections(VariableNo,2) == 1 .AND. &
                 VariableDirections(VariableNo,3) == 2 ) THEN
              DO i=1,N3
                DO j=1,N2
                  DO k=1,N1
                     PInData(dataread, 1) = Zcoord(k)
                     PInData(dataread, 2) = Xcoord(j)
                     PInData(dataread, 3) = Ycoord(i)
                     PInData(dataread, 4) = ddata_in3D(k,j,i)
                     dataread = dataread+1
                  END DO
                END DO
              END DO
           ELSE IF (VariableDirections(VariableNo,1) == 3 .AND. VariableDirections(VariableNo,2) == 2 .AND. &
                 VariableDirections(VariableNo,3) == 1 ) THEN
              DO i=1,N3
                DO j=1,N2
                  DO k=1,N1
                     PInData(dataread, 1) = Zcoord(k)
                     PInData(dataread, 2) = Ycoord(j)
                     PInData(dataread, 3) = Xcoord(i)
                     PInData(dataread, 4) = ddata_in3D(k,j,i)
                     dataread = dataread+1
                  END DO
                END DO
              END DO
            END IF
        ELSE
            CALL FATAL(SolverName, 'Unsupported dimensions in data copy.')
        END IF

       dataread = dataread-1

       !------------------------------------------------------
       ! interpolate the values
       !------------------------------------------------------
       reinitiate = .TRUE.
       !------------------------------------------------------
       ! Loop all active elements of solver
       !------------------------------------------------------
       DO elementNumber=1,Solver % NumberOFActiveElements
         ! PRINT *,'dataread=',dataread  
         !PRINT *,'e=',elementNumber,'/',Solver % Mesh % NumberOfBulkElements
          CurrentElement => GetActiveElement(elementNumber)
          IF (.NOT.ASSOCIATED(CurrentElement)) CALL FATAL(SolverName,'Element pointer not associated')
          Equation => GetEquation()
          IF ( ASSOCIATED( Equation ) ) THEN           
              InterpolateVariable = ListGetLogical( Equation, &
                   'Interpolate ' // TRIM(AssociatedVariableName(VariableNo)), Found)
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
                   AssocVariableDirections(VariableNo,1:3), &
                   2.0_dp, &
                   reinitiate,&
                   SupportingPoints(VariableNo),&
                   SolverName)
              !WRITE(*,*) VarVal(VarPerm(CurrentElement % NodeIndexes(i)))   
          END DO
          IF (elementNumber==1) reinitiate = .FALSE.
       END DO ! DO elementNumber
       WRITE(*,*) 'Done variable'

       DEALLOCATE(PInData, Xcoord, YCoord)
       IF (IsInteger(VariableNo)) THEN
         IF (NetCDFDIM(VariableNo) == 2) THEN  
            DEALLOCATE(idata_in2D)       
         ELSE IF (NetCDFDIM(VariableNo) == 3) THEN
            DEALLOCATE(idata_in3D)
         END IF
       ELSE
         IF (NetCDFDIM(VariableNo) == 2) THEN
            DEALLOCATE(ddata_in2D)
         ELSE IF (NetCDFDIM(VariableNo) == 3) THEN
            DEALLOCATE(ddata_in3D)      
         END IF
       END IF

    IF (NetCDFDIM(VariableNo) == 3 .AND. NoDIM(VariableNo) == 2) THEN
      IF (VariableDirections(VariableNo,1) == 1 .AND. VariableDirections(VariableNo,2) == 2 .AND. &
                 VariableDirections(VariableNo,3) == 3 ) THEN
         IF(VariableSizeInDirections(VariableNo,3) > 1) THEN
           incrTime(VariableNo) = incrTime(VariableNo)+1
           IF (incrTime(VariableNo) > VariableSizeInDirections(VariableNo,3)) THEN
              incrTime(VariableNo) = VariableSizeInDirections(VariableNo,3)
           END IF
           WRITE(*,*) '********************',incrTime(VariableNo) 
         END IF
      ELSE IF (VariableDirections(VariableNo,1) == 1 .AND. VariableDirections(VariableNo,2) == 3 .AND. &
                 VariableDirections(VariableNo,3) == 2 ) THEN
         IF(VariableSizeInDirections(VariableNo,2) > 1) THEN
           incrTime(VariableNo) = incrTime(VariableNo)+1
           IF (incrTime(VariableNo) > VariableSizeInDirections(VariableNo,2)) THEN
              incrTime(VariableNo) = VariableSizeInDirections(VariableNo,2)
           END IF
         END IF
      ELSE IF (VariableDirections(VariableNo,1) == 2 .AND. VariableDirections(VariableNo,2) == 1 .AND. &
                 VariableDirections(VariableNo,3) == 3 ) THEN
         IF(VariableSizeInDirections(VariableNo,3) > 1) THEN
           incrTime(VariableNo) = incrTime(VariableNo)+1
           IF (incrTime(VariableNo) > VariableSizeInDirections(VariableNo,3)) THEN
              incrTime(VariableNo) = VariableSizeInDirections(VariableNo,3)
           END IF
         END IF
      ELSE IF (VariableDirections(VariableNo,1) == 2 .AND. VariableDirections(VariableNo,2) == 3 .AND. &
                 VariableDirections(VariableNo,3) == 1 ) THEN
         IF(VariableSizeInDirections(VariableNo,2) > 1) THEN
           incrTime(VariableNo) = incrTime(VariableNo)+1
           IF (incrTime(VariableNo) > VariableSizeInDirections(VariableNo,2)) THEN
              incrTime(VariableNo) = VariableSizeInDirections(VariableNo,2)
           END IF
         END IF
      ELSE IF (VariableDirections(VariableNo,1) == 3 .AND. VariableDirections(VariableNo,2) == 1 .AND. &
                 VariableDirections(VariableNo,3) == 2 ) THEN
         IF(VariableSizeInDirections(VariableNo,1) > 1) THEN
           incrTime(VariableNo) = incrTime(VariableNo)+1
           IF (incrTime(VariableNo) > VariableSizeInDirections(VariableNo,1)) THEN
              incrTime(VariableNo) = VariableSizeInDirections(VariableNo,1)
           END IF
         END IF
      ELSE IF (VariableDirections(VariableNo,1) == 3 .AND. VariableDirections(VariableNo,2) == 2 .AND. &
                 VariableDirections(VariableNo,3) == 1 ) THEN
         IF(VariableSizeInDirections(VariableNo,1) > 1) THEN
           incrTime(VariableNo) = incrTime(VariableNo)+1
           IF (incrTime(VariableNo) > VariableSizeInDirections(VariableNo,1)) THEN
              incrTime(VariableNo) = VariableSizeInDirections(VariableNo,1)
           END IF
         END IF
      END IF
    ELSE IF (NetCDFDIM(VariableNo) == 2 .AND. NoDIM(VariableNo) == 2 .AND. TransientSimulation) THEN
       !CALL FATAL(SolverName, 'Cant figure out what to interpolate in time.')
       !CALL FATAL(SolverName, 'NetCDF variable and target variable are 2D, then only 2D spatial data')
    ELSE IF (NetCDFDIM(VariableNo) == 3 .AND. NoDIM(VariableNo) == 3 .AND. TransientSimulation) THEN
       !CALL FATAL(SolverName, 'Cant figure out what to interpolate in time.')
       !CALL FATAL(SolverName, 'NetCDF variable and target variable are 3D, then only 3D spatial data')
    END IF
     
    END DO !DO VariableNo
 END IF
                  
 CONTAINS
!-----------------------------------------------------------------------------------------------------------
 SUBROUTINE check(status)
    INTEGER, INTENT (IN) :: status
    
    IF(status /= nf90_noerr) THEN
      print *, trim(nf90_strerror(status))
      stop "Stopped"
    END IF
 END SUBROUTINE check  

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
     !theInterpolatedValue = closestPoints(1)
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


END SUBROUTINE ReadNETCDF

