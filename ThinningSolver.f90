! ******************************************************************************
! *
! *                    Author:       Hakime Seddik
! *
! *                    Address: Institute of Low Temperature Science
! *                             Hokkaido University
! *                             email:hakime@pop.lowtem.hokudai.ac.jp
! *
! *                    Date: May 1 2009
! *
! *                    Modified from: AdvectionReaction
! *
! *****************************************************************************/

!----------------------------------------------------------------------------------    
RECURSIVE SUBROUTINE ThinningSolver( Model,Solver,dt,TransientSimulation )
!---------------------------------------------------------------------------------- 

  USE DefUtils

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

 INTEGER :: nvalue, i, j, k, t, ThinFunctionDOFs, DepthDOFs, ThinningDOFs, DIM, elementNbNodes
 TYPE(Element_t),POINTER :: BoundaryElement, CurrentElement
 REAL(KIND=dp), POINTER :: ThinFunctionValues(:), DepthValues(:), ThinningValues(:)
 TYPE(Variable_t), POINTER :: ThinFunctionSol, DepthSol, ThinningSol
 TYPE(ValueList_t),POINTER :: Material, BodyForce
 INTEGER, POINTER :: ThinFunctionPerm(:), DepthPerm(:), ThinningPerm(:), NodeIndexes(:)
 REAL(KIND=dp) :: volumefraction
 INTEGER :: material_id, body_id, bf_id 
 LOGICAL :: GotIt

 ThinFunctionSol => VariableGet( Solver % Mesh % Variables, 'ThinFunction' )
 IF ( ASSOCIATED( ThinFunctionSol ) ) THEN
       ThinFunctionPerm    => ThinFunctionSol % Perm
       ThinFunctionValues  => ThinFunctionSol % Values
       ThinFunctionDOFs = ThinFunctionSol % DOFs
 ELSE
       CALL Info('Thinning:', 'Could not find Thinning Function variable', Level=4)
 END IF

 DepthSol => VariableGet( Solver % Mesh % Variables, 'Depth' )
 IF ( ASSOCIATED( DepthSol ) ) THEN
        DepthPerm => DepthSol % Perm
        DepthValues => DepthSol % Values
        DepthDOFs = DepthSol % DOFs
 ELSE
 WRITE(Message,'(a,I6)') 'Could not find Depth field variable'
 CALL FATAL('Thinning:',Message)
 END IF

 ThinningSol => VariableGet( Solver % Mesh % Variables, 'Thinning' )
 IF ( ASSOCIATED( ThinningSol ) ) THEN
        ThinningPerm => ThinningSol % Perm
        ThinningValues => ThinningSol % Values
        ThinningDOFs = ThinningSol % DOFs
 ELSE
 WRITE(Message,'(a,I6)') 'Could not find Thinning field variable'
 CALL FATAL('Thinning:',Message)
 END IF

 dim = CoordinateSystemDimension()  

 DO i=1,Solver % Mesh % NumberOfNodes
  
  IF (DepthValues(DepthDOFs*(DepthPerm(i)-1)+1) < 0.0D00) THEN    
     DepthValues(DepthDOFs*(DepthPerm(i)-1)+1) = 0.0D00  
  END IF  

  volumefraction = 910.0*( 1.0D00 - 0.55 * EXP(-3.8D-02*DepthValues(DepthDOFs*(DepthPerm(i)-1)+1)) )

  ThinningValues(ThinningDOFs*(ThinningPerm(i)-1)+1) = (volumefraction/910.0)* &
                   ThinFunctionValues(ThinFunctionDOFs*(ThinFunctionPerm(i)-1)+1)

 END DO

 WRITE(Message,*) 'solve done', minval(ThinningValues), maxval(ThinningValues)
      CALL Info( 'ThinningSolver', Message, Level=4 )

WRITE(Message,'(a)') '------------------------------------------------------'
CALL Info('Thinning:',Message, Level=3)
CALL Info('Thinning:','Compute Thinning:..........done', Level=3)
WRITE(Message,'(a)') '------------------------------------------------------'
CALL Info('Thinning:',Message, Level=3)

!------------------------------------------------------------------------------
END SUBROUTINE ThinningSolver
!------------------------------------------------------------------------------

