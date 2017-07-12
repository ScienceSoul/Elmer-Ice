!----------------------------------------------------------------------------------    
RECURSIVE SUBROUTINE SmoothViscosity( Model,Solver,dt,TransientSimulation )
!---------------------------------------------------------------------------------- 

  USE DefUtils
  USE Types

  IMPLICIT NONE
!------------------------------------------------------------------------------
!******************************************************************************
!
!  Solve stress equations for one timestep
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
  TYPE(Element_t),POINTER :: CurrentElement
  INTEGER, POINTER :: NodeIndexes(:)
  TYPE(ValueList_t), POINTER :: SolverParams

  INTEGER :: nvalue, i, j, k, N, t, iter, elementNbNodes, ViscDOFs, TempHomoDOFs, DIM
  INTEGER :: NBSmoothingIter
  REAL(KIND=dp), POINTER :: ViscValues(:), TempHomoValues(:)
  TYPE(Variable_t), POINTER :: ViscSol, TempHomoSol
  INTEGER, POINTER :: ViscPerm(:), TempHomoPerm(:)
  REAL(KIND=dp) :: getArrheniusFactor, EF, homoTemp, viscFact, smoothedViscFact
  LOGICAL :: AllocationsDone = .FALSE., ifSmooth, Found

  SAVE AllocationsDone
 

  !------------------------------------------------------------------------------
  !     Allocate some permanent storage, this is done first time only
  !------------------------------------------------------------------------------
      IF ( .NOT. AllocationsDone .OR. Solver % Mesh % Changed) THEN
        N = Model % MaxElementNodes
        dim = CoordinateSystemDimension()
       AllocationsDone = .TRUE. 
      END IF
   
 ViscSol => VariableGet( Solver % Mesh % Variables, 'viscosity factor' )
 IF ( ASSOCIATED( ViscSol ) ) THEN
       ViscPerm    => ViscSol % Perm
       ViscValues  => ViscSol % Values
       ViscDOFs = ViscSol % DOFs
 ELSE
       CALL FATAL('SmoothViscosity', 'Could not find viscosity factor variable')
 END IF

 TempHomoSol => VariableGet( Solver % Mesh % Variables, 'Temp Homologous' )
 IF ( ASSOCIATED( TempHomoSol ) ) THEN
        TempHomoPerm => TempHomoSol % Perm
        TempHomoValues => TempHomoSol % Values
        TempHomoDOFs = TempHomoSol % DOFs
 ELSE
 WRITE(Message,'(a,I6)') 'Could not find Temp Homologous variable'
 CALL FATAL('SmoothViscosity',Message)
 END IF

 SolverParams => GetSolverParams()
 IF (.NOT. ASSOCIATED(SolverParams)) THEN
       CALL FATAL('SmoothViscosity','No Solver section found')
 END IF

 ifSmooth = GetLogical(SolverParams, 'Smooth Viscosity', Found)
 IF(.NOT.Found) THEN
    ifSmooth = .FALSE.
     CALL Info('SmoothViscosity','No smoothing is done to viscosity factor', Level=3)
 END IF 

 NBSmoothingIter =  ListGetInteger(SolverParams, 'Smoothing Max Iterations', Found)
 IF(.NOT.Found) THEN
   CALL Info('SmoothViscosity','Smoothing Max Iterations not found', Level=3)
   CALL Info('SmoothViscosity','Set it to 1', Level=3)
   NBSmoothingIter = 1  
 END IF 

 EF = GetConstReal(SolverParams, 'Enhancement Factor', Found)
 IF(.NOT.Found) THEN
   CALL Info('SmoothViscosity','Enhancement factor not found', Level=3)
   CALL Info('SmoothViscosity','Set EF=1.0', Level=3)
   EF = 1.0
 ELSE
   WRITE(Message,'(a,F4.1,a)') 'Use value ', EF, ' as Enhancement factor'
   CALL Info('SmoothViscosity', Message, Level=3)
 END IF


WRITE(Message,'(a)') '------------------------------------------------------'
CALL Info('SmoothViscosity',Message, Level=3)
CALL Info('SmoothViscosity','Computed viscosity factor for all mesh', Level=3)
WRITE(Message,'(a)') '------------------------------------------------------'
CALL Info('SmoothViscosity',Message, Level=3)

 DO i=1,Solver % Mesh % NumberOfNodes

    j = ViscPerm(i)
    k = TempHomoPerm(i)

    homoTemp = TempHomoValues(TempHomoDOFs*(k-1)+1)

    IF(homoTemp < -10.0) THEN
        getArrheniusFactor = 3.985E-13*EXP(-60.0E03/(8.314*(273.16+homoTemp)))
    ELSE IF(homoTemp > 0.0) THEN
         getArrheniusFactor = 1.916E03*EXP(-139.0E03/(8.314*(273.16)))
    ELSE
          getArrheniusFactor = 1.916E03*EXP(-139.0E03/(8.314*(273.16+homoTemp)))
    END IF

   viscFact = (2.0*EF*getArrheniusFactor)**(-1.0/3.0)
   viscFact = viscFact*31556926.0**(-1.0/3.0)*1.0E-06

   ViscValues(ViscDOFs*(j-1)+1) = viscFact
 END DO

IF(ifSmooth) THEN

 Do iter=1, NBSmoothingIter

  DO t=1,Solver % NumberOFActiveElements

     CurrentElement => GetActiveElement(t)
     elementNbNodes = GetElementNOFNodes( CurrentElement)
     NodeIndexes => CurrentElement % NodeIndexes
     smoothedViscFact = 0.0
     
     DO  i=1,elementNbNodes
       smoothedViscFact = smoothedViscFact+ViscValues(ViscDOFs*(ViscPerm(NodeIndexes(i))-1)+1)    
     END DO
     smoothedViscFact = smoothedViscFact/elementNbNodes  
     DO  i=1,elementNbNodes
        ViscValues(ViscDOFs*(ViscPerm(NodeIndexes(i))-1)+1) = smoothedViscFact
     END DO
 
  END DO


!  DO t=1, Solver % Mesh % NumberOfBoundaryElements 

!    CurrentElement => GetBoundaryElement(t)
!    IF ( .NOT.ActiveBoundaryElement() ) CYCLE 
!    elementNbNodes = GetElementNOFNodes(CurrentElement)
!    NodeIndexes => CurrentElement % NodeIndexes
!    IF ( GetElementFamily() == 1 ) CYCLE 
!    smoothedViscFact = 0.0

!     DO  i=1,elementNbNodes
!       smoothedViscFact = smoothedViscFact+ViscValues(ViscDOFs*(ViscPerm(NodeIndexes(i))-1)+1)    
!     END DO
!     smoothedViscFact = smoothedViscFact/elementNbNodes 
!     DO  i=1,elementNbNodes
!        ViscValues(ViscDOFs*(ViscPerm(NodeIndexes(i))-1)+1) = smoothedViscFact
!     END DO

!  END Do

 END DO

 WRITE(Message,'(a)') '------------------------------------------------------'
 CALL Info('SmoothViscosity',Message, Level=3)
 CALL Info('SmoothViscosity','Viscosity factor smoothing on all nodes:..........done', Level=3)
 WRITE(Message,'(a,I4,a)') 'Did it: ', NBSmoothingIter, ' time(s)'
 CALL Info('SmoothViscosity',Message, Level=3)
 WRITE(Message,'(a)') '------------------------------------------------------'
 CALL Info('SmoothViscosity',Message, Level=3)  

END IF

!------------------------------------------------------------------------------
END SUBROUTINE SmoothViscosity
!------------------------------------------------------------------------------
