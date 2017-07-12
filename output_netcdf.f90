RECURSIVE SUBROUTINE output_nc( Model,Solver,Timestep,TransientSimulation )
    use pnetcdf
    USE DefUtils
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
 INTEGER :: i,j,n,t, LocalNodes, DIM, istat, FlowDOFs, SIAFlowDOFs, TempDOFs, HomoTempDOFs, DepthDOFs, FSDOFs, FSRefDOFs, &
            BasalMeltDOFS, incr, Nsize, asDOFs, ObsFlowxDOFs, ObsFlowyDOFs, BetaDOFs, StressDOFs, indx
 INTEGER(8) :: Nsize2, TimeSimul2, MeshNumberOFNodes, Reduced
 INTEGER :: TimeSimul, nbNodesOnSurface
 TYPE(Solver_t), POINTER :: PointerToSolver
 TYPE(ValueList_t),POINTER :: SolverParams
 TYPE(Element_t),POINTER :: Element
 TYPE(ValueList_t), POINTER :: BC
 TYPE(Nodes_t) :: ElementNodes
 TYPE(Variable_t), POINTER :: FlowSol, SIAFlowSol, TempSol, HomoTempSol, DepthSol, FSSol, FSRefSol, BasalMeltSol, asSol, &
                              ObsFlowxSol, ObsFlowySol, BetaSol, StressSol
 REAL(KIND=dp), POINTER :: Flow(:), SIAFlow(:), Temp(:), HomoTemp(:), Depth(:), FS(:), FSRef(:), BasalMelt(:), as(:), ObsFlowx(:), &
                           ObsFlowy(:), Beta(:), Stress(:)
 INTEGER, POINTER :: FlowPerm(:), SIAFlowPerm(:), TempPerm(:), HomoTempPerm(:), DepthPerm(:), FSPerm(:), FSRefPerm(:), &
                     BasalMeltPerm(:), NodeIndexes(:), asPerm(:), ObsFlowxPerm(:), ObsFlowyPerm(:), BetaPerm(:), StressPerm(:)
 REAL(KIND=dp), ALLOCATABLE:: XCoord(:), YCoord(:), DataBuffer(:)

 INTEGER, ALLOCATABLE :: mask(:), elements(:)
 CHARACTER(LEN=MAX_NAME_LEN) :: SolverName, FILE_NAME, TYPEOfOutput, VarName
 LOGICAL :: AllocationsDone=.FALSE., FirstTime=.TRUE., Found, countNodes
 !------------------------------------------------------------------------------
 !     Sone stuff needed for NETCDF reading
 !------------------------------------------------------------------------------
 INTEGER :: NDIMS
 INTEGER :: ncid, varid1, varid2, varid3, varid4, varid5, varid6
 INTEGER :: varid7, varid8, varid9, varid10, varid11, varid12, varid13, varid14, varid15, varid16, varid17, varid18, varid19, &
            varid20, varid21, varid22, varid23, varid24, varid25, varid26, varid27, varid28
 INTEGER, ALLOCATABLE ::  dimids(:)
 INTEGER :: x_dimid, y_dimid, z_dimid, t_dimid
 INTEGER(KIND=MPI_OFFSET_KIND), ALLOCATABLE :: start(:), count(:)
 INTEGER(KIND=MPI_OFFSET_KIND) :: unit=1

 !------------------------------------------------------------------------------
 !     MPI stuff: number of processors, rank of this processor, and error 
 !     code.
 !------------------------------------------------------------------------------
 INTEGER(8) :: pp
 INTEGER :: p, my_rank, ierr
 
 SAVE AllocationsDone, FirstTime, DIM, Nsize, Nsize2, XCoord, YCoord, start, count, dimids, nbNodesOnSurface, NDIMS, FILE_NAME
 SAVE SolverName, TYPEOfOutput, TimeSimul, TimeSimul2, incr, my_rank, p, ncid, varid1, varid2, varid3, varid4, varid5, varid6
 SAVE varid7, varid8, varid9, varid10, varid11, varid12, varid13, varid14, varid15, varid16, varid17, varid18, varid19, varid20 
 SAVE varid21, varid22, varid23, varid24, varid25, varid26, varid27, varid28, mask, elements, DataBuffer, VarName
 
 SolverParams => Solver % Values

 IF (FirstTime) WRITE(SolverName,'(A)') 'output_netcdf'
 PointerToSolver => Solver
 IF ( .NOT. ASSOCIATED( PointerToSolver ) ) THEN
    CALL FATAL(SolverName, ' No Solver Pointer associated')
 END IF
    
 !------------------------------------------------------------------
 ! Do these things for the first time only
 !------------------------------------------------------------------
 IF ( .NOT. AllocationsDone .OR. Solver % Mesh % Changed) THEN
 
    DIM = CoordinateSystemDimension()
    Nsize = Model % Mesh % NumberOfNodes

    indx = 0
    Nsize2 = 0
    DO t = 1, Solver % NumberOfActiveElements

        Element => GetActiveElement(t)
        n = GetElementNOFNodes()
        NodeIndexes => Element % NodeIndexes

        DO i=1,n
            indx = indx+1
            Nsize2 = Nsize2+1
        END DO
    END DO

    IF ( AllocationsDone ) THEN
        DEALLOCATE( XCoord, YCoord, mask, elements, DataBuffer)
    END IF

    ALLOCATE( XCoord(Nsize), YCoord(indx), mask(indx), elements(indx), DataBuffer(indx), STAT=istat)
    IF ( istat /= 0 ) THEN
          CALL Fatal(SolverName, 'Memory allocation error.' )
    END IF
    CALL Info(SolverName,'Memory allocations done', Level=3)

    VarName = GetString(SolverParams,'Flow Solver Name',Found )
    IF(.NOT. Found) THEN
          CALL Fatal(SolverName, 'Variable for flow solution not found')
    END IF

    AllocationsDone = .TRUE.

  END IF
  
  IF (FirstTime) THEN
  
    FirstTime = .FALSE.
    CALL INFO(SolverName,'(Re-)Initialization started.',Level=3)

    ! GET dimensions
    NDIMS = GetInteger( Solver % Values, 'Dimension of variables', Found)
    IF (.NOT. Found) THEN 
           Call FATAL(SolverName, 'Cant find dimensions for variable: ')
    END IF
  
    ! Get the NETCDF file name to write on disk
    FILE_NAME = GetString( Solver % Values, 'File Name', Found)
    IF(.NOT. Found) THEN
       Call FATAL(SolverName, 'Cant find a name to use for the NETCDF file.')
    END IF

    TYPEOfOutput = GetString( Solver % Values, 'Type of netcdf output')
    IF(.NOT. Found) THEN
       Call FATAL(SolverName, 'Cant find a type of netcdf output to use.')
    END IF


    IF ( TransientSimulation ) THEN

        TimeSimul = GetInteger( GetSimulation(),'NETCDF Timestep Intervals', Found)
        IF (.NOT. Found) THEN
            CALL FATAL(SolverName,'Could not find parameter NETCDF Timestep Intervals')
        END IF

    END IF


    IF ( TransientSimulation ) THEN

        TimeSimul2 = GetInteger( GetSimulation(), 'NETCDF Timestep Intervals', Found)
        IF (.NOT. Found) THEN
            CALL FATAL(SolverName,'Could not find parameter NETCDF Timestep Intervals')
        END IF

    ELSE

    TimeSimul2 = 1

    END IF

    SELECT CASE(TYPEOfOutput)

    CASE('searise')
  
        IF (NDIMS < 3) THEN 
            CALL FATAL(SolverName, 'Given variable dimension for netcdf not consistent with searise type. Should be three.')
        END IF
  
    CASE DEFAULT
  
        CALL FATAL(SolverName, 'Type of netcdf output not supported in initialization')

    END SELECT

    ! These will tell where in the data file this processor should write.
    ALLOCATE(start(NDIMS), count(NDIMS), STAT=istat)
    IF(istat .NE. 0) THEN
      CALL FATAL(SolverName, 'Error in allocation')
    END IF

    ALLOCATE(dimids(NDIMS), STAT=istat)
    IF(istat .NE. 0) THEN
      CALL FATAL(SolverName, 'Error in allocation')
    END IF

    CALL MPI_Comm_rank(MPI_COMM_WORLD, my_rank, ierr) 
    CALL MPI_Comm_size(MPI_COMM_WORLD, p, ierr)

    ! Set the netcdf staff
    ! Create the netCDF file. The comm and info parameters cause parallel I/O to be enabled.
    CALL check( nfmpi_create(MPI_COMM_WORLD, FILE_NAME, IOR(NF_CLOBBER, NF_64BIT_OFFSET), MPI_INFO_NULL, ncid) )


    ! The number of columns in variable which is the max value of nodes accros all partitions
    CALL MPI_ALLREDUCE(Nsize2,Reduced,1,MPI_INT,MPI_MAX,MPI_COMM_WORLD,ierr)  
    WRITE(*,*) 'Nb. nodes, nb. elements, part. nb, tot. nb. nodes:', Nsize2,  Solver % NumberOfActiveElements, &
               ParEnv % myPe, Reduced

    MeshNumberOFNodes = Reduced

    ! Define the dimensions. NetCDF will hand back an ID for 
    ! each. Metadata operations must take place on all processors.
    SELECT CASE(TYPEOfOutput)
 
    CASE('searise')
    
        pp = p
        CALL check( nfmpi_def_dim(ncid, "x", pp, x_dimid) )
        CALL check( nfmpi_def_dim(ncid, "y", MeshNumberOFNodes, y_dimid) )
        CALL check( nfmpi_def_dim(ncid, "time", TimeSimul2, t_dimid) )

        dimids = (/ y_dimid, x_dimid, t_dimid /)

    CASE DEFAULT
     
        CALL Fatal(SolverName, 'Type of netcdf output not supported in defining dimensions')

    END SELECT


    call check( nfmpi_def_var(ncid, "zs", NF_DOUBLE, 3, dimids, varid1) )
    call check( nfmpi_def_var(ncid, "xcoord", NF_DOUBLE,3, dimids, varid2) )
    call check( nfmpi_def_var(ncid, "ycoord", NF_DOUBLE, 3, dimids, varid3) )
    call check( nfmpi_def_var(ncid, "flowx", NF_DOUBLE, 3, dimids, varid4) )
    call check( nfmpi_def_var(ncid, "flowy", NF_DOUBLE, 3, dimids, varid5) )
    call check( nfmpi_def_var(ncid, "flowz", NF_DOUBLE, 3, dimids, varid6) )
    call check( nfmpi_def_var(ncid, "siaflowx", NF_DOUBLE, 3, dimids, varid7) )
    call check( nfmpi_def_var(ncid, "siaflowy", NF_DOUBLE, 3, dimids, varid8) )
    call check( nfmpi_def_var(ncid, "siaflowz", NF_DOUBLE, 3, dimids, varid9) )
    call check( nfmpi_def_var(ncid, "temp", NF_DOUBLE, 3, dimids, varid10) )
    call check( nfmpi_def_var(ncid, "thick", NF_DOUBLE, 3, dimids, varid11) )
    call check( nfmpi_def_var(ncid, "zcoord", NF_DOUBLE, 3, dimids, varid12) )
    call check( nfmpi_def_var(ncid, "meltrate", NF_DOUBLE, 3, dimids, varid13) )
    call check( nfmpi_def_var(ncid, "as", NF_DOUBLE, 3, dimids, varid14) )
    call check( nfmpi_def_var(ncid, "refzs", NF_DOUBLE, 3, dimids, varid15) )
    call check( nfmpi_def_var(ncid, "mask", NF_DOUBLE, 3, dimids, varid16) )
    call check( nfmpi_def_var(ncid, "elements", NF_DOUBLE, 3, dimids, varid17) )
    call check( nfmpi_def_var(ncid, "temphomo", NF_DOUBLE, 3, dimids, varid18) )
    call check( nfmpi_def_var(ncid, "pressure", NF_DOUBLE, 3, dimids, varid19) )
    call check( nfmpi_def_var(ncid, "obsflowx", NF_DOUBLE, 3, dimids, varid20) )
    call check( nfmpi_def_var(ncid, "obsflowy", NF_DOUBLE, 3, dimids, varid21) )
    call check( nfmpi_def_var(ncid, "beta", NF_DOUBLE, 3, dimids, varid22) )
    call check( nfmpi_def_var(ncid, "sxx", NF_DOUBLE, 3, dimids, varid23) )
    call check( nfmpi_def_var(ncid, "syy", NF_DOUBLE, 3, dimids, varid24) )
    call check( nfmpi_def_var(ncid, "szz", NF_DOUBLE, 3, dimids, varid25) )
    call check( nfmpi_def_var(ncid, "sxy", NF_DOUBLE, 3, dimids, varid26) )
    call check( nfmpi_def_var(ncid, "syz", NF_DOUBLE, 3, dimids, varid27) )
    call check( nfmpi_def_var(ncid, "sxz", NF_DOUBLE, 3, dimids, varid28) )

    call check( nfmpi_enddef(ncid) )

    incr = 1

 END IF ! FirstTime---------------------------------------------------------------------------------------------
 
 ! Get field variables
 SELECT CASE(TYPEOfOutput)
 
    CASE('searise')
    
        FlowSol => VariableGet( Solver % Mesh % Variables, VarName)
        IF ( ASSOCIATED( FlowSol ) ) THEN
            FlowPerm    => FlowSol % Perm
            Flow        => FlowSol % Values
            FlowDOFs    =  FlowSol % DOFs
        ELSE
            CALL Fatal(SolverName, 'Could not find velocity field variable')
        END IF

        SIAFlowSol => VariableGet( Solver % Mesh % Variables, 'SIAFlow' )
        IF ( ASSOCIATED( SIAFlowSol ) ) THEN
            SIAFlowPerm     => SIAFlowSol % Perm
            SIAFlow         => SIAFlowSol % Values
            SIAFlowDOFs     =  SIAFlowSol % DOFs
        ELSE
            CALL Fatal(SolverName, 'Could not find SIAFlow variable')
        END IF

        TempSol => VariableGet( Solver % Mesh % Variables, 'Temp' )
        IF ( ASSOCIATED( TempSol ) ) THEN
            TempPerm    => TempSol % Perm
            Temp        => TempSol % Values
            TempDOFs    =  TempSol % DOFs
        ELSE
            CALL Fatal(SolverName, 'Could not find Temp field variable')
        END IF

        HomoTempSol => VariableGet( Solver % Mesh % Variables, 'Temp Homologous' )
        IF ( ASSOCIATED( HomoTempSol ) ) THEN
            HomoTempPerm   => HomoTempSol % Perm
            HomoTemp       => HomoTempSol % Values
            HomoTempDOFs   =  HomoTempSol % DOFs
        ELSE
            CALL FATAL(SolverName, 'Could not find Temp Homologous field variable')
        END IF

        DepthSol => VariableGet( Solver % Mesh % Variables, 'Depth' ) 
        IF ( ASSOCIATED( DepthSol ) ) THEN
            DepthPerm   => DepthSol % Perm
            Depth       => DepthSol % Values
            DepthDOFs   =  DepthSol % DOFs
        ELSE
            CALL FATAL(SolverName, 'Could not find Depth field variable')
        END IF

        FSSol => VariableGet( Solver % Mesh % Variables, 'FS' )
        IF ( ASSOCIATED( FSSol ) ) THEN
            FSPerm   => FSSol % Perm
            FS       => FSSol % Values
            FSDOFs   =  FSSol % DOFs
        ELSE
            CALL Fatal(SolverName, 'Could not find FS field variable')
        END IF

        FSRefSol => VariableGet( Solver % Mesh % Variables, 'ReferenceFS' )
        IF ( ASSOCIATED( FSRefSol ) ) THEn
            FSRefPerm   => FSRefSol % Perm
            FSRef       => FSRefSol % Values
            FSRefDOFs   = FSRefSol % DOFs
        ELSE
            CALL Fatal(SolverName, 'Could not find ReferenceFS field variable')
        END IF

        BasalMeltSol => VariableGet( Solver % Mesh % Variables, 'BasalMelting' )
        IF ( ASSOCIATED( BasalMeltSol ) ) THEN
            BasalMeltPerm   => BasalMeltSol % Perm
            BasalMelt       => BasalMeltSol % Values
            BasalMeltDOFS   =  BasalMeltSol % DOFs
        ELSE
            CALL Fatal(SolverName, 'Could not find BasalMelting field variable')
        END IF

        asSol => VariableGet( Solver % Mesh % Variables, 'as' )
        IF ( ASSOCIATED( asSol ) ) THEN
            asPerm   => asSol % Perm
            as       => asSol % Values
            asDOFs   =  asSol % DOFs
        ELSE
            CALL FATAL(SolverName, 'Could not find as field variable')
        END IF

        ObsFlowxSol => VariableGet( Solver % Mesh % Variables, 'vxdata' )
        IF ( ASSOCIATED( ObsFlowxSol ) ) THEN
            ObsFlowxPerm   => ObsFlowxSol % Perm
            ObsFlowx       => ObsFlowxSol % Values
            ObsFlowxDOFs   =  ObsFlowxSol % DOFs
        ELSE
            CALL FATAL(SolverName, 'Could not find vxdata field variable')
        END IF

        ObsFlowySol => VariableGet( Solver % Mesh % Variables, 'vydata' )
        IF ( ASSOCIATED( ObsFlowySol ) ) THEN
            ObsFlowyPerm   => ObsFlowySol % Perm
            ObsFlowy       => ObsFlowySol % Values
            ObsFlowyDOFs   =  ObsFlowySol % DOFs
        ELSE
            CALL FATAL(SolverName, 'Could not find vydata field variable')
        END IF

        BetaSol => VariableGet( Solver % Mesh % Variables, 'Beta' )
        IF ( ASSOCIATED( BetaSol ) ) THEN
            BetaPerm   => BetaSol % Perm
            Beta       => BetaSol % Values
            BetaDOFs   =  BetaSol % DOFs
        ELSE
            CALL FATAL(SolverName, 'Could not find Beta field variable')
        END IF

        StressSol => VariableGet( Solver % Mesh % Variables, 'Stress' )
        IF ( ASSOCIATED( StressSol ) ) THEN
            StressPerm   => StressSol % Perm
            Stress       => StressSol % Values
            StressDOFs   =  StressSol % DOFs
        ELSE
            CALL FATAL(SolverName, 'Could not find Stress field variable')
        END IF

    CASE DEFAULT

        CALL Fatal(SolverName, 'Type of netcdf output not supported in retrieving field variables')
    
 END SELECT
   
 !Copy field variable to buffers, define variable and write to netcdf file
 SELECT CASE(TYPEOfOutput)
 
 CASE('searise')
 
    DataBuffer = 0.0_dp
    mask = 0
 
    ! zs variable
    indx = 1
    DO t = 1, Solver % NumberOfActiveElements
    
        Element => GetActiveElement(t)
        n = GetElementNOFNodes()
        NodeIndexes => Element % NodeIndexes
 
        CALL GetElementNodes( ElementNodes, Element, Solver)
   
        DO i=1,n
            XCoord(indx) = ElementNodes % x(i)
            YCoord(indx) = ElementNodes % y(i)
            DataBuffer(indx) = FS(FSDOFs*(FSPerm(NodeIndexes(i))-1)+1)
            mask(indx) = 1
            elements(indx) = t
            indx = indx+1
        END DO
   
    END DO

    start = (/ 1, my_rank + 1, incr/)
    count = (/ Nsize2, unit, unit /)

    call check( nfmpi_put_vara_double_all(ncid, varid1, start, count, DataBuffer) )
    call check( nfmpi_put_vara_double_all(ncid, varid2, start, count, XCoord) )
    CALL check( nfmpi_put_vara_double_all(ncid, varid3, start, count, Ycoord) )
    CALL check( nfmpi_put_vara_int_all(ncid, varid16, start, count, mask) )
    CALL check( nfmpi_put_vara_int_all(ncid, varid17, start, count, elements) )
    
    DataBuffer = 0.0_dp
    
    ! flowx varible
    indx = 1
    DO t = 1, Solver % NumberOfActiveElements

        Element => GetActiveElement(t)
        n = GetElementNOFNodes()
        NodeIndexes => Element % NodeIndexes
  
        CALL GetElementNodes( ElementNodes, Element, Solver)
   
        DO i=1,n
            DataBuffer(indx) = Flow(FlowDOFs*(FlowPerm(NodeIndexes(i))-1)+1)
            indx = indx+1
        END DO
   
    END DO

    call check( nfmpi_put_vara_double_all(ncid, varid4, start, count, DataBuffer) )

    DataBuffer = 0.0_dp

    ! flowy variable
    indx = 1
    DO t = 1, Solver % NumberOfActiveElements

        Element => GetActiveElement(t)
        n = GetElementNOFNodes()
        NodeIndexes => Element % NodeIndexes
  
        CALL GetElementNodes( ElementNodes, Element, Solver)
   
        DO i=1,n    
            DataBuffer(indx) = Flow(FlowDOFs*(FlowPerm(NodeIndexes(i))-1)+2)
            indx = indx+1
        END DO
   
    END DO

    call check( nfmpi_put_vara_double_all(ncid, varid5, start, count, DataBuffer) )

    DataBuffer = 0.0_dp

    ! flowz variable
    indx = 1
    DO t = 1, Solver % NumberOfActiveElements

        Element => GetActiveElement(t)
        n = GetElementNOFNodes()
        NodeIndexes => Element % NodeIndexes
  
        CALL GetElementNodes( ElementNodes, Element, Solver)
   
        DO i=1,n
            DataBuffer(indx) = Flow(FlowDOFs*(FlowPerm(NodeIndexes(i))-1)+3)
            indx = indx+1
        END DO
   
    END DO

    call check( nfmpi_put_vara_double_all(ncid, varid6, start, count, DataBuffer) )

    DataBuffer = 0.0_dp
    
    ! SIAFlowx varible
    indx = 1
    DO t = 1, Solver % NumberOfActiveElements

        Element => GetActiveElement(t)
        n = GetElementNOFNodes()
        NodeIndexes => Element % NodeIndexes
  
        CALL GetElementNodes( ElementNodes, Element, Solver)
   
        DO i=1,n
            DataBuffer(indx) = SIAFlow(SIAFlowDOFs*(SIAFlowPerm(NodeIndexes(i))-1)+1)
            indx = indx+1
        END DO
   
    END DO

    call check( nfmpi_put_vara_double_all(ncid, varid7, start, count, DataBuffer) )

    DataBuffer = 0.0_dp

    ! SIAFlowy variable
    indx = 1
    DO t = 1, Solver % NumberOfActiveElements

        Element => GetActiveElement(t)
        n = GetElementNOFNodes()
        NodeIndexes => Element % NodeIndexes
  
        CALL GetElementNodes( ElementNodes, Element, Solver)
   
        DO i=1,n
            DataBuffer(indx) = SIAFlow(SIAFlowDOFs*(SIAFlowPerm(NodeIndexes(i))-1)+2)
            indx = indx+1
        END DO
   
    END DO

    call check( nfmpi_put_vara_double_all(ncid, varid8, start, count, DataBuffer) )

    DataBuffer = 0.0_dp

    ! SIAFlowz variable
    indx = 1
    DO t = 1, Solver % NumberOfActiveElements

        Element => GetActiveElement(t)
        n = GetElementNOFNodes()
        NodeIndexes => Element % NodeIndexes
  
        CALL GetElementNodes( ElementNodes, Element, Solver)
   
        DO i=1,n
            DataBuffer(indx) = SIAFlow(SIAFlowDOFs*(SIAFlowPerm(NodeIndexes(i))-1)+3)
            indx = indx+1  
        END DO
   
    END DO

    call check( nfmpi_put_vara_double_all(ncid, varid9, start, count, DataBuffer) )

    DataBuffer = 0.0_dp

    ! pressure variable
    indx = 1
    DO t = 1, Solver % NumberOfActiveElements

        Element => GetActiveElement(t)
        n = GetElementNOFNodes()
        NodeIndexes => Element % NodeIndexes
  
        CALL GetElementNodes( ElementNodes, Element, Solver)
  
        DO i=1,n
            DataBuffer(indx) = Flow(FlowDOFs*(FlowPerm(NodeIndexes(i))-1)+4)
            indx = indx+1
        END DO
   
    END DO

    call check( nfmpi_put_vara_double_all(ncid, varid19, start, count, DataBuffer) )

    DataBuffer = 0.0_dp

    ! temp variable
    indx = 1
    DO t = 1, Solver % NumberOfActiveElements

        Element => GetActiveElement(t)
        n = GetElementNOFNodes()
        NodeIndexes => Element % NodeIndexes
 
        CALL GetElementNodes( ElementNodes, Element, Solver)
   
        DO i=1,n
            DataBuffer(indx) = Temp(TempDOFs*(TempPerm(NodeIndexes(i))-1)+1)
            indx = indx+1
        END DO
   
    END DO

    call check( nfmpi_put_vara_double_all(ncid, varid10, start, count, DataBuffer) )

    DataBuffer = 0.0_dp

    ! temp homo variable
    indx = 1
    DO t = 1, Solver % NumberOfActiveElements

        Element => GetActiveElement(t)
        n = GetElementNOFNodes()
        NodeIndexes => Element % NodeIndexes
  
        CALL GetElementNodes( ElementNodes, Element, Solver)
   
        DO i=1,n
            DataBuffer(indx) = HomoTemp(HomoTempDOFs*(HomoTempPerm(NodeIndexes(i))-1)+1)
            indx = indx+1
        END DO
   
    END DO

    call check( nfmpi_put_vara_double_all(ncid, varid18, start, count, DataBuffer) )

    DataBuffer = 0.0_dp

    ! thick variable
    indx = 1
    DO t = 1, Solver % NumberOfActiveElements

        Element => GetActiveElement(t)
        n = GetElementNOFNodes()
        NodeIndexes => Element % NodeIndexes

        CALL GetElementNodes( ElementNodes, Element, Solver)
   
        DO i=1,n
            DataBuffer(indx) = Depth(DepthDOFs*(DepthPerm(NodeIndexes(i))-1)+1)
            indx = indx+1
      END DO
   
    END DO

    call check( nfmpi_put_vara_double_all(ncid, varid11, start, count, DataBuffer) )

    DataBuffer = 0.0_dp

    ! z variable
    indx = 1
    DO t = 1, Solver % NumberOfActiveElements

        Element => GetActiveElement(t)
        n = GetElementNOFNodes()
        NodeIndexes => Element % NodeIndexes

        CALL GetElementNodes( ElementNodes, Element, Solver)
   
        DO i=1,n
            DataBuffer(indx) = ElementNodes % z(i)
            indx = indx+1
        END DO
   
    END DO

    call check( nfmpi_put_vara_double_all(ncid, varid12, start, count, DataBuffer) )

    DataBuffer = 0.0_dp

    ! basalmelt variable
    indx = 1
    DO t = 1, Solver % NumberOfActiveElements

        Element => GetActiveElement(t)
        n = GetElementNOFNodes()
        NodeIndexes => Element % NodeIndexes
  
        CALL GetElementNodes( ElementNodes, Element, Solver)
   
        DO i=1,n
            DataBuffer(indx) = Basalmelt(BasalMeltDOFS*(BasalMeltPerm(NodeIndexes(i))-1)+1)
            indx = indx+1
        END DO
   
    END DO

    call check( nfmpi_put_vara_double_all(ncid, varid13, start, count, DataBuffer) )

    DataBuffer = 0.0_dp

    ! as variable
    indx = 1
    DO t = 1, Solver % NumberOfActiveElements

        Element => GetActiveElement(t)
        n = GetElementNOFNodes()
        NodeIndexes => Element % NodeIndexes
  
        CALL GetElementNodes( ElementNodes, Element, Solver)
  
        DO i=1,n
            DataBuffer(indx) = as(asDOFs*(asPerm(NodeIndexes(i))-1)+1)
            indx = indx+1
        END DO
   
    END DO

    call check( nfmpi_put_vara_double_all(ncid, varid14, start, count, DataBuffer) )

    DataBuffer = 0.0_dp

    ! refzs variable
    indx = 1
    DO t = 1, Solver % NumberOfActiveElements

        Element => GetActiveElement(t)
        n = GetElementNOFNodes()
        NodeIndexes => Element % NodeIndexes
  
        CALL GetElementNodes( ElementNodes, Element, Solver)
   
        DO i=1,n
            DataBuffer(indx) = FSRef(FSRefDOFs*(FSRefPerm(NodeIndexes(i))-1)+1)
            indx = indx+1
        END DO
   
    END DO

    call check( nfmpi_put_vara_double_all(ncid, varid15, start, count, DataBuffer) )

     DataBuffer = 0.0_dp

    ! Observed velo x variable
    indx = 1
    DO t = 1, Solver % NumberOfActiveElements

        Element => GetActiveElement(t)
        n = GetElementNOFNodes()
        NodeIndexes => Element % NodeIndexes
  
        CALL GetElementNodes( ElementNodes, Element, Solver)
   
        DO i=1,n
            DataBuffer(indx) = ObsFlowx(ObsFlowxDOFs*(ObsFlowxPerm(NodeIndexes(i))-1)+1)
            indx = indx+1
        END DO
   
    END DO

    call check( nfmpi_put_vara_double_all(ncid, varid20, start, count, DataBuffer) )

     DataBuffer = 0.0_dp

    ! Observed velo y variable
    indx = 1
    DO t = 1, Solver % NumberOfActiveElements

        Element => GetActiveElement(t)
        n = GetElementNOFNodes()
        NodeIndexes => Element % NodeIndexes
  
        CALL GetElementNodes( ElementNodes, Element, Solver)
   
        DO i=1,n
            DataBuffer(indx) = ObsFlowy(ObsFlowyDOFs*(ObsFlowyPerm(NodeIndexes(i))-1)+1)
            indx = indx+1
        END DO
   
    END DO

    call check( nfmpi_put_vara_double_all(ncid, varid21, start, count, DataBuffer) )

    DataBuffer = 0.0_dp

    ! Beta variable
    indx = 1
    DO t = 1, Solver % NumberOfActiveElements

        Element => GetActiveElement(t)
        n = GetElementNOFNodes()
        NodeIndexes => Element % NodeIndexes
  
        CALL GetElementNodes( ElementNodes, Element, Solver)
  
        DO i=1,n
            DataBuffer(indx) = 10.0_dp**Beta(BetaDOFs*(BetaPerm(NodeIndexes(i))-1)+1)
            indx = indx+1
        END DO
   
    END DO

    call check( nfmpi_put_vara_double_all(ncid, varid22, start, count, DataBuffer) )

    DataBuffer = 0.0_dp
    
    ! Sxx variable
    indx = 1
    DO t = 1, Solver % NumberOfActiveElements

        Element => GetActiveElement(t)
        n = GetElementNOFNodes()
        NodeIndexes => Element % NodeIndexes
  
        CALL GetElementNodes( ElementNodes, Element, Solver)
   
        DO i=1,n
            DataBuffer(indx) = Stress(StressDOFs*(StressPerm(NodeIndexes(i))-1)+1)
            indx = indx+1
        END DO
   
    END DO

    call check( nfmpi_put_vara_double_all(ncid, varid23, start, count, DataBuffer) )

    DataBuffer = 0.0_dp
    
    ! Syy variable
    indx = 1
    DO t = 1, Solver % NumberOfActiveElements

        Element => GetActiveElement(t)
        n = GetElementNOFNodes()
        NodeIndexes => Element % NodeIndexes
  
        CALL GetElementNodes( ElementNodes, Element, Solver)
   
        DO i=1,n
            DataBuffer(indx) = Stress(StressDOFs*(StressPerm(NodeIndexes(i))-1)+2)
            indx = indx+1
        END DO
   
    END DO

    call check( nfmpi_put_vara_double_all(ncid, varid24, start, count, DataBuffer) )

    DataBuffer = 0.0_dp
    
    ! Szz variable
    indx = 1
    DO t = 1, Solver % NumberOfActiveElements

        Element => GetActiveElement(t)
        n = GetElementNOFNodes()
        NodeIndexes => Element % NodeIndexes
  
        CALL GetElementNodes( ElementNodes, Element, Solver)
   
        DO i=1,n
            DataBuffer(indx) = Stress(StressDOFs*(StressPerm(NodeIndexes(i))-1)+3)
            indx = indx+1
        END DO
   
    END DO

    call check( nfmpi_put_vara_double_all(ncid, varid25, start, count, DataBuffer) )

    DataBuffer = 0.0_dp
    
    ! Sxy variable
    indx = 1
    DO t = 1, Solver % NumberOfActiveElements

        Element => GetActiveElement(t)
        n = GetElementNOFNodes()
        NodeIndexes => Element % NodeIndexes
  
        CALL GetElementNodes( ElementNodes, Element, Solver)
   
        DO i=1,n
            DataBuffer(indx) = Stress(StressDOFs*(StressPerm(NodeIndexes(i))-1)+4)
            indx = indx+1
        END DO
   
    END DO

    call check( nfmpi_put_vara_double_all(ncid, varid26, start, count, DataBuffer) )

    DataBuffer = 0.0_dp
    
    ! Syz variable
    indx = 1
    DO t = 1, Solver % NumberOfActiveElements

        Element => GetActiveElement(t)
        n = GetElementNOFNodes()
        NodeIndexes => Element % NodeIndexes
  
        CALL GetElementNodes( ElementNodes, Element, Solver)
   
        DO i=1,n
            DataBuffer(indx) = Stress(StressDOFs*(StressPerm(NodeIndexes(i))-1)+5)
            indx = indx+1
        END DO
   
    END DO

    call check( nfmpi_put_vara_double_all(ncid, varid27, start, count, DataBuffer) )

    DataBuffer = 0.0_dp
    
    ! Sxz variable
    indx = 1
    DO t = 1, Solver % NumberOfActiveElements

        Element => GetActiveElement(t)
        n = GetElementNOFNodes()
        NodeIndexes => Element % NodeIndexes
  
        CALL GetElementNodes( ElementNodes, Element, Solver)
   
        DO i=1,n
            DataBuffer(indx) = Stress(StressDOFs*(StressPerm(NodeIndexes(i))-1)+6)
            indx = indx+1
        END DO
   
    END DO

    call check( nfmpi_put_vara_double_all(ncid, varid28, start, count, DataBuffer) )

 CASE DEFAULT
     
    CALL Fatal(SolverName, 'Type of netcdf output not supported in defining dimensions')

 END SELECT
 
 Call INFO(SolverName, 'Write netcdf data:.....done.')
 
 IF ( TransientSimulation ) THEN
 
    IF( incr == TimeSimul ) THEN
 
        ! Close the file. This frees up any internal netCDF resources 
        ! associated with the file, and flushes any buffers.
        CALL Info(SolverName,'Closing netcdf file used for output', Level=3)
        call check( nfmpi_close(ncid) )
 
    ELSE
        incr = incr+1
  
    END IF
 ELSE
    CALL Info(SolverName,'Closing netcdf file used for output', Level=3)
    call check( nfmpi_close(ncid) )
 END IF
 
 
 CONTAINS
 !-----------------------------------------------------------------------------------------------------------
 
  SUBROUTINE check(status)
   INTEGER, INTENT (IN) :: status
   
   IF(status /= nf_noerr) THEN
        PRINT *, trim(nfmpi_strerror(status))
        STOP 2
   END IF
  END SUBROUTINE check

END SUBROUTINE output_nc



