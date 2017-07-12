!----------------------------------------------------------------------------------    
RECURSIVE SUBROUTINE ComputeIntegralTerm( Model,Solver,dt,TransientSimulation )
!---------------------------------------------------------------------------------- 

  USE DefUtils

  IMPLICIT NONE
!------------------------------------------------------------------------------
!******************************************************************************
!
!  Initialize Stokes fields with SIA velocities and pressure
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
 
 INTEGER :: nvalue, i, j, k, l, m, n, o, FlowDOFs, SicoVxDOFS, SicoVyDOFS, SicoVzDOFS, SicoPressureDOFS, &
            DIM, istat
 TYPE(Element_t), POINTER :: CurrentElement
 TYPE(ValueList_t),POINTER :: SolverParams
 TYPE(ValueList_t), POINTER :: Material
 REAL(KIND=dp), POINTER :: DepthValues(:), IntegralValues(:)
 REAL(KIND=dp) :: arrheniusFactor, enh
 TYPE(Variable_t), POINTER :: Depth, Integral
 INTEGER, POINTER :: DepthPerm(:), IntegralPerm(:), NodeIndexes(:)
 CHARACTER(LEN=MAX_NAME_LEN) :: SolverName
 LOGICAL :: AllocationsDone = .FALSE., Found

 SolverName = 'ComputeIntegralTerm'

 SolverParams => Solver % Values

 M = Model % MaxElementNodes
 dim = CoordinateSystemDimension()  

 Depth => VariableGet( Solver % Mesh % Variables, 'Depth' )
 IF ( ASSOCIATED( Depth ) ) THEN
       DepthPerm    => Depth % Perm
       DepthValues  => Depth % Values
 ELSE
       CALL Fatal(SolverName, 'Could not find Depth field variable')
 END IF

  Integral => VariableGet( Model % Variables, 'IntegralTerm' )
  IF ( ASSOCIATED( Integral ) ) THEN
        IntegralPerm    => Integral % Perm
        IntegralValues  =>  Integral % Values
  ELSE
        CALL FATAL(SolverName,'No variable Integral term found')  
  END IF

  DO i=1,Solver % NumberOFActiveElements

     CurrentElement => GetActiveElement(i)
     NodeIndexes => CurrentElement % NodeIndexes

     Material => GetMaterial()

     arrheniusFactor = GetConstReal(Material, 'Integral Arrhenius factor', Found)
     IF(.NOT. Found) THEN 
        CALL FATAL(SolverName, 'Integral Arrhenius factor not found in Material')
     END IF
     arrheniusFactor = arrheniusFactor * 31556926.0_dp * 10.0_dp**18.0_dp

     enh = GetConstReal(Material, 'Integral Enhancement factor', Found)
     IF(.NOT. Found) THEN 
        CALL FATAL(SolverName, 'Integral Enhancement factor not found in Material')
     END IF

     DO k=1, GetElementNOFNodes(CurrentElement)
   
        IntegralValues(IntegralPerm(NodeIndexes(k))) =  &
              enh * arrheniusFactor * (DepthValues(DepthPerm(NodeIndexes(k))))**(3.0_dp+1.0_dp)
     END DO
  END DO
  WRITE(Message,'(a)') '---------------------------------------------------------------'
  CALL Info(SolverName,Message, Level=3)
  CALL Info(SolverName,'Compute Integral term:..........done', Level=3)
  WRITE(Message,'(a)') '---------------------------------------------------------------'
  CALL Info(SolverName,Message, Level=3)

!------------------------------------------------------------------------------
END SUBROUTINE ComputeIntegralTerm
!------------------------------------------------------------------------------
