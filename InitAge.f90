!------------------------------------------------------------------------------
   RECURSIVE SUBROUTINE InitAge( Model,Solver,dt,TransientSimulation )
!------------------------------------------------------------------------------

    USE DefUtils

    IMPLICIT NONE
!------------------------------------------------------------------------------
!******************************************************************************
!
!  Solve dating equation for one timestep (da/dt + da/dxi.ui = 1)
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

 INTEGER :: i,j,k,l,n,t, n1
 
 TYPE(Element_t),POINTER :: Element  

 TYPE(Variable_t), POINTER :: DepthSol, HeightSol, AgeSol, CompAgeSol

 TYPE(ValueList_t), POINTER :: Material

 INTEGER, POINTER :: DepthPerm(:), HeightPerm(:), AgePerm(:), CompAgePerm(:)

 REAL(KIND=dp), POINTER :: DepthValues(:), HeightValues(:), AgeValues(:), CompAgeValues(:), Ref(:)

 INTEGER :: DepthDOFs, HeightDOFs, AgeDOFs, CompAgeDOFs, Indexes(128), NDOFS

 LOGICAL ::  GotIt, AllocationsDone = .FALSE.

 REAL(KIND=dp) :: z_tildes, c1, c2, c3, c4, zstar, vz_b, AgeDimless, AgeDim

 REAL(KIND=dp) :: iceThickness, as

 INTEGER :: dim
 
!------------------------------------------------------------------------------
!    Get variables needed for solution
!------------------------------------------------------------------------------

  dim = CoordinateSystemDimension()   

  DepthSol => VariableGet( Solver % Mesh % Variables, 'Depth' )
  IF ( ASSOCIATED( DepthSol ) ) THEN
        DepthPerm => DepthSol % Perm
        DepthValues => DepthSol % Values
        DepthDOFs = DepthSol % DOFs
  ELSE
  WRITE(Message,'(a,I6)') 'Could not find Depth field variable'
  CALL FATAL('InitAge',Message)
  END IF

  HeightSol => VariableGet( Solver % Mesh % Variables, 'Height' )
  IF ( ASSOCIATED( HeightSol ) ) THEN
        HeightPerm => HeightSol % Perm
        HeightValues => HeightSol % Values
        HeightDOFs = HeightSol % DOFs
  ELSE
  WRITE(Message,'(a,I6)') 'Could not find Height field variable'
  CALL FATAL('InitAge',Message)
  END IF

  CompAgeSol =>  VariableGet( Solver % Mesh % Variables, 'CompAge' )
  IF ( ASSOCIATED( CompAgeSol ) ) THEN
        CompAgePerm => CompAgeSol % Perm
        CompAgeValues => CompAgeSol % Values
        CompAgeDOFs = CompAgeSol % DOFs
  ELSE
  WRITE(Message,'(a,I6)') 'Could not find CompAge field variable'
  CALL FATAL('InitAge',Message)
  END IF

  AgeSol =>  VariableGet( Solver % Mesh % Variables, 'Age' )
  IF ( ASSOCIATED( AgeSol ) ) THEN
        AgePerm => AgeSol % Perm
        AgeValues => AgeSol % Values
        AgeDOFs = AgeSol % DOFs
  ELSE
  WRITE(Message,'(a,I6)') 'Could not find Age field variable'
  CALL FATAL('InitAge',Message)
  END IF

!--------------------------------------------------------------------------------
! We work in dimensionless quantitites, cf Greve and others, An Glaciol, 35, 2002
!--------------------------------------------------------------------------------

as = ListGetConstReal( Model % Constants,'Accumulation rate',GotIt)
IF (.NOT. GotIt) THEN
   WRITE(Message,'(A)') 'No accumulation rate found in Constant '
   CALL Fatal('InitAge)', Message)
ELSE 
   WRITE(Message,'(A,F10.4)') 'Accumulation = ', as
   CALL INFO('InitAge', Message, Level = 25)
END IF

vz_b = ListGetConstReal( Model % Constants,'Basal velocity offset',GotIt)
IF (.NOT. GotIt) THEN
   WRITE(Message,'(A)') 'No basal velocity offset found in Constant '
   CALL Fatal('InitAge)', Message)
   ELSE 
   WRITE(Message,'(A,F10.4)') 'Basal velocity offset = ', vz_b
   CALL INFO('InitAge', Message, Level = 25)
END IF

zstar = (1.0/3.0)

!------------------------------------------------------------------------------
DO t=1,Solver % NumberOFActiveElements
!------------------------------------------------------------------------------

 Element => GetActiveElement(t) 
 n = GetElementNOFNodes(Element) 
 NDOFS = GetElementDOFs(Indexes, Element)

     DO i=1,n

        k = Element % NodeIndexes(i)


        iceThickness = DepthValues(DepthPerm(k))+HeightValues(HeightPerm(k))
        z_tildes = HeightValues(HeightPerm(k))/iceThickness

        c1 = (2.0*(1.0+vz_b))/(2.0-zstar)
        c2 = (zstar+2.0*vz_b)/(2.0-zstar)
        c3 = (1.0+vz_b)/(zstar*(2.0-zstar))
        c4 = -vz_b

        IF (z_tildes >= zstar) THEN 
        AgeDimless = (1.0/c1)*LOG(1.0/(c1*z_tildes-c2))
        ELSE 
           IF (z_tildes <= zstar) THEN
           AgeDimless = 1.0/(SQRT(c3*c4))*(ATAN(SQRT(c3*zstar)/SQRT(c4)) &
-ATAN(SQRT(c3*z_tildes)/SQRT(c4)))+1.0/c1*LOG(1.0/(c1*zstar-c2))
          END IF
        END IF

        AgeDim = (iceThickness/as)*AgeDimless
        CompAgeValues(CompAgePerm(Indexes(i))) = AgeDim
     END DO
END DO

! Average the elemental results to nodal values for continuous variable Age:
!---------------------------------------------------------------------------
n1 = Solver % Mesh % NumberOfNodes
      ALLOCATE( Ref(n1) )
      Ref = 0.
      AgeValues = 0.
      DO t=1,Solver % NumberOfActiveElements
         Element => GetActiveElement(t) 
         NDOFS = GetElementDOFs(Indexes, Element)
         n = GetElementNOFNodes(Element)
         DO i=1,n
            k = Element % NodeIndexes(i)
            AgeValues( AgePerm(k) ) =    & 
               AgeValues( AgePerm(k) ) + &
                  CompAgeValues(CompAgePerm(Indexes(i)) )
           Ref(k) = Ref(k) + 1
        END DO
      END DO

      DO i=1,n1
         j = AgePerm(i)
         IF (j < 1) CYCLE
         IF ( Ref(i) > 0 ) THEN
            AgeValues( AgePerm(i) ) = &
                AgeValues( AgePerm(i) ) / Ref(i)
         END IF
      END DO
  DEALLOCATE( Ref )

WRITE(Message,'(a)') '-----------------------------------------------'
CALL Info('InitAge',Message, Level=3)
CALL Info('InitAge','Initialize values for age:..........done', Level=3)
WRITE(Message,'(a)') '-----------------------------------------------'
CALL Info('InitAge',Message, Level=3)

!------------------------------------------------------------------------------
 END SUBROUTINE InitAge
!------------------------------------------------------------------------------
