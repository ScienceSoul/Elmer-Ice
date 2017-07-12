!---------------------------------------------------------------------------------------------------    
RECURSIVE SUBROUTINE SmoothTempHomo( Model,Solver,dt,TransientSimulation )
!--------------------------------------------------------------------------------------------------- 

  USE DefUtils
  USE Types

  IMPLICIT NONE
!------------------------------------------------------------------------------
!******************************************************************************
!
!  Smooth the temperature relative to melting point
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

  INTEGER :: nvalue, i, j, k, N, t, iter, elementNbNodes, SmoothHomoDOFs, & 
  TempHomoDOFs, DIM
  INTEGER :: NBSmoothingIter
  REAL(KIND=dp), POINTER :: SmoothHomoValues(:), TempHomoValues(:)
  TYPE(Variable_t), POINTER :: SmoothHomoSol, TempHomoSol
  INTEGER, POINTER :: SmoothHomoPerm(:), TempHomoPerm(:)
  REAL(KIND=dp) :: HomoTemp, smoothedHomoTemp
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
   
 SmoothHomoSol => VariableGet( Solver % Mesh % Variables, & 
     'Smoothed Homologous Temperature' )
 
IF ( ASSOCIATED( SmoothHomoSol ) ) THEN
       SmoothHomoPerm    => SmoothHomoSol % Perm
       SmoothHomoValues  => SmoothHomoSol% Values
       SmoothHomoDOFs = SmoothHomoSol % DOFs
 ELSE
       CALL FATAL('SmoothHomoTemp', 'Could not find Smoothed Homologous Temperature  variable')
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
       CALL FATAL('SmoothHomoTemp','No Solver section found')
 END IF

 ifSmooth = GetLogical(SolverParams, 'Smooth Homologous', Found)
 IF(.NOT.Found) THEN
    ifSmooth = .FALSE.
    CALL Info('SmoothHomoTemp','No smoothing is done to Homologous Temperature',&
     Level=3)
 END IF 

 NBSmoothingIter =  ListGetInteger(SolverParams, 'Smoothing Max Iterations', Found)
 IF(.NOT.Found) THEN
   CALL Info('SmoothHomo','Smoothing Max Iterations not found', Level=3)
   CALL Info('SmoothHomo','Set it to 1', Level=3)
   NBSmoothingIter = 1  
 END IF 

IF(ifSmooth) THEN

 Do iter=1, NBSmoothingIter

  DO t=1,Solver % NumberOFActiveElements

     CurrentElement => GetActiveElement(t)
     elementNbNodes = GetElementNOFNodes( CurrentElement)
     NodeIndexes => CurrentElement % NodeIndexes
     smoothedHomoTemp = 0.0
     
     DO  i=1,elementNbNodes
       smoothedHomoTemp = smoothedHomoTemp+TempHomoValues(TempHomoDOFs* &
       (TempHomoPerm(NodeIndexes(i))-1)+1)    
     END DO
     smoothedHomoTemp = smoothedHomoTemp/elementNbNodes  
     DO  i=1,elementNbNodes
        SmoothHomoValues(SmoothHomoDOFs*(SmoothHomoPerm(NodeIndexes(i))-1)+1)&
           = smoothedHomoTemp
     END DO
 
  END DO


 END DO

 WRITE(Message,'(a)') '------------------------------------------------------'
 CALL Info('SmoothHomo',Message, Level=3)
 CALL Info('SmoothHomo','Homologous temperature  smoothing on all nodes:..........done', Level=3)
 WRITE(Message,'(a,I4,a)') 'Did it: ', NBSmoothingIter, ' time(s)'
 CALL Info('SmoothHomo',Message, Level=3)
 WRITE(Message,'(a)') '------------------------------------------------------'
 CALL Info('SmoothHomo',Message, Level=3)  

END IF

!------------------------------------------------------------------------------
END SUBROUTINE SmoothTempHomo
!------------------------------------------------------------------------------
