!----------------------------------------------------------------------------------    
RECURSIVE SUBROUTINE InitThickness( Model,Solver,dt,TransientSimulation )
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
 REAL(KIND=dp), POINTER :: SurfValues(:), BedValues(:), ThicknessValues(:)
 TYPE(Variable_t), POINTER :: SurfSol, BedSol, ThicknessSol
 INTEGER, POINTER :: SurfPerm(:), BedPerm(:), ThicknessPerm(:), NodeIndexes(:)
 CHARACTER(LEN=MAX_NAME_LEN) :: SolverName
 LOGICAL :: AllocationsDone = .FALSE., Found, done = .FALSE.

 SAVE done

 SolverName = 'InitThickness'

 SolverParams => Solver % Values

 M = Model % MaxElementNodes
 dim = CoordinateSystemDimension()  

 SurfSol => VariableGet( Solver % Mesh % Variables, 'surfData' )
 IF ( ASSOCIATED( SurfSol ) ) THEN
       SurfPerm    => SurfSol % Perm
       SurfValues  => SurfSol % Values
 ELSE
       CALL Fatal(SolverName, 'Could not find surfData field variable')
 END IF

 BedSol => VariableGet( Solver % Mesh % Variables, 'bedData' )
 IF ( ASSOCIATED( BedSol ) ) THEN
       BedPerm    => BedSol % Perm
       BedValues  => BedSol % Values
 ELSE
       CALL Fatal(SolverName, 'Could not find bedData field variable')
 END IF

 ThicknessSol => VariableGet( Solver % Mesh % Variables, 'H' )
 IF ( ASSOCIATED( ThicknessSol ) ) THEN
       ThicknessPerm    => ThicknessSol % Perm
       ThicknessValues  => ThicknessSol % Values
 ELSE
       CALL Fatal(SolverName, 'Could not find H field variable')
 END IF


 IF (.NOT. done) THEN
    DO i=1,Solver % NumberOFActiveElements

      CurrentElement => GetActiveElement(i)
      NodeIndexes => CurrentElement % NodeIndexes

      DO k=1, GetElementNOFNodes(CurrentElement)
   
        ThicknessValues(ThicknessPerm(NodeIndexes(k))) =  &
              SurfValues(SurfPerm(NodeIndexes(k))) - BedValues(BedPerm(NodeIndexes(k))) 

      END DO
    END DO
    WRITE(Message,'(a)') '---------------------------------------------------------------'
    CALL Info(SolverName,Message, Level=3)
    CALL Info(SolverName,'Initialize Thickness fields:..........done', Level=3)
    WRITE(Message,'(a)') '---------------------------------------------------------------'
    CALL Info(SolverName,Message, Level=3)

    done = .TRUE.
 ELSE
    CALL Info(SolverName, 'No ops....')
 END IF

!------------------------------------------------------------------------------
END SUBROUTINE InitThickness
!------------------------------------------------------------------------------
