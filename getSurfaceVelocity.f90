FUNCTION getSurfaceVelocity1(Model, nodenumber, adummy) RESULT(velocity)
  USE Types
  USE CoordinateSystems
  USE SolverUtils
  USE ElementDescription
  USE DefUtils
!-------------------------------------------------------------------
  IMPLICIT NONE
!-------------------------external variables-----------------------
  TYPE(Model_t) :: Model 
  TYPE(Variable_t), POINTER :: VeloSol
  TYPE(Element_t), POINTER :: BoundaryElement, ParentElement, BDElement
  TYPE(ValueList_t), POINTER :: BC, ParentMaterial, BodyForce
  TYPE(Nodes_t) :: BoundaryElementNodes, BDElementNodes
  INTEGER, POINTER :: VeloPerm(:)
  REAL(KIND=dp), POINTER :: Velo(:)
  INTEGER :: material_id, body_id, other_body_id, bf_id
  INTEGER :: i, j, n, nb1, nb2, t, bc_id, istat, nodenumber, NBoundary, BoundaryElementNode
  REAL(KIND=dp) :: adummy, velocity
  LOGICAL :: GotIt, foundNodes

!----------------! get some information upon active boundary element and its parent  !----------------------------------------

  BoundaryElement => Model % CurrentElement

!-----------------------Do nothing if this is a halo-element of a parallel mesh?  !------------------------------------------

  IF (ParEnv % myPe .NE. BoundaryElement % partIndex) RETURN

  IF ( .NOT. ASSOCIATED(BoundaryElement) ) THEN     
    CALL FATAL('getSurfaceVelocity','No boundary element found')  
  END IF

  NBoundary = BoundaryElement % Type % NumberOfNodes  
  DO BoundaryElementNode=1,NBoundary     
   IF (nodenumber .EQ. BoundaryElement % NodeIndexes(BoundaryElementNode)) THEN        
      GotIt = .TRUE.        
      EXIT     
   END IF  
  END DO

  IF (.NOT.GotIt) THEN
    CALL WARN('getSurfaceVelocity','Node not found in Current Element')     
    velo = 0.0D00
    RETURN
  END IF

  other_body_id = BoundaryElement % BoundaryInfo % outbody
  IF (other_body_id < 1) THEN ! only one body in calculation
     ParentElement => BoundaryElement % BoundaryInfo % Right
     IF ( .NOT. ASSOCIATED(ParentElement) ) ParentElement => BoundaryElement % BoundaryInfo % Left
  ELSE ! we are dealing with a body-body boundary and asume that the normal is pointing outwards
     ParentElement => BoundaryElement % BoundaryInfo % Right
     IF (ParentElement % BodyId == other_body_id) ParentElement => BoundaryElement % BoundaryInfo % Left
  END IF

! just to be on the save side, check again
!-----------------------------------------
  IF ( .NOT. ASSOCIATED(ParentElement) ) THEN
     WRITE(Message,'(A,I10,A)')&
          'Parent Element for Boundary element no. ',&
          BoundaryElement % ElementIndex, ' not found'
     !CALL FATAL('getSurfaceVelocity',Message)
  END IF

  VeloSol => VariableGet( Model % Variables, "vxdata" )
  IF ( ASSOCIATED( VeloSol ) ) THEN
       VeloPerm => VeloSol % Perm
       Velo     => VeloSol % Values
  ELSE
       WRITE(Message, '(A)') 'Could not find vxdata variable '
       CALL FATAL('getSurfaceVelocity',Message)
  END IF

  N = Model % Mesh % MaxElementNodes

  ALLOCATE(BoundaryElementNodes % x( N ), BoundaryElementNodes % y( N ), BoundaryElementNodes % z( N ), STAT=istat)
  IF ( istat /= 0 ) THEN
      CALL FATAL('getSurfaceVelocity', 'Memory allocation error' )
  END IF

  ALLOCATE(BDElementNodes % x( N ), BDElementNodes % y( N ), BDElementNodes % z( N ), STAT=istat)
  IF ( istat /= 0 ) THEN
      CALL FATAL('getSurfaceVelocity', 'Memory allocation error' )
  END IF


  CALL GetElementNodes( BoundaryElementNodes, BoundaryElement )
  nb1 =  GetElementNOFNodes(BoundaryElement)

  foundNodes = .FALSE.
  DO t=1,  Model % Mesh % NumberOfBoundaryElements
     BDElement => GetBoundaryElement(t)
     IF ( .NOT.ActiveBoundaryElement(BDElement) ) CYCLE
     IF ( GetElementFamily(BDElement) == 1 ) CYCLE
     bc_id = GetBCId( BDElement )
     IF (bc_id == 6) THEN
         nb2 = GetElementNOFNodes(BDElement)
         CALL GetElementNodes( BDElementNodes, BDElement )
         DO i=1,nb2
           DO j=1, nb1
              IF ((BoundaryElementNodes % x(j) == BDElementNodes % x(i) .AND. &
                  BoundaryElementNodes % y(j) == BDElementNodes % y(i) ) .AND. BoundaryElement % NodeIndexes(j) == nodenumber) THEN
                  velocity = Velo(VeloPerm(BDElement % NodeIndexes(i)))
                  foundNodes = .True.
                  EXIT
              END IF
           END DO
           IF (foundNodes) EXIT
         END DO
         IF (foundNodes) EXIT
     END IF
  END DO

  DEALLOCATE(BoundaryElementNodes % x, BoundaryElementNodes % y, BoundaryElementNodes % z)
  DEALLOCATE(BDElementNodes % x, BDElementNodes % y, BDElementNodes % z)

END FUNCTION getSurfaceVelocity1


FUNCTION getSurfaceVelocity2(Model, nodenumber, adummy) RESULT(velocity)
  USE Types
  USE CoordinateSystems
  USE SolverUtils
  USE ElementDescription
  USE DefUtils
!-------------------------------------------------------------------
  IMPLICIT NONE
!-------------------------external variables-----------------------
  TYPE(Model_t) :: Model 
  TYPE(Variable_t), POINTER :: VeloSol
  TYPE(Element_t), POINTER :: BoundaryElement, ParentElement, BDElement
  TYPE(ValueList_t), POINTER :: BC, ParentMaterial, BodyForce
  TYPE(Nodes_t) :: BoundaryElementNodes, BDElementNodes
  INTEGER, POINTER :: VeloPerm(:)
  REAL(KIND=dp), POINTER :: Velo(:)
  INTEGER :: material_id, body_id, other_body_id, bf_id
  INTEGER :: i, j, n, nb1, nb2, t, bc_id, istat, nodenumber, NBoundary, BoundaryElementNode
  REAL(KIND=dp) :: adummy, velocity
  LOGICAL :: GotIt, foundNodes

!----------------! get some information upon active boundary element and its parent  !----------------------------------------

  BoundaryElement => Model % CurrentElement

!-----------------------Do nothing if this is a halo-element of a parallel mesh?  !------------------------------------------

  IF (ParEnv % myPe .NE. BoundaryElement % partIndex) RETURN

  IF ( .NOT. ASSOCIATED(BoundaryElement) ) THEN     
    CALL FATAL('getSurfaceVelocity','No boundary element found')  
  END IF

  NBoundary = BoundaryElement % Type % NumberOfNodes  
  DO BoundaryElementNode=1,NBoundary     
   IF (nodenumber .EQ. BoundaryElement % NodeIndexes(BoundaryElementNode)) THEN        
      GotIt = .TRUE.        
      EXIT     
   END IF  
  END DO

  IF (.NOT.GotIt) THEN
    CALL WARN('getSurfaceVelocity','Node not found in Current Element')     
    velo = 0.0D00
    RETURN
  END IF

  other_body_id = BoundaryElement % BoundaryInfo % outbody
  IF (other_body_id < 1) THEN ! only one body in calculation
     ParentElement => BoundaryElement % BoundaryInfo % Right
     IF ( .NOT. ASSOCIATED(ParentElement) ) ParentElement => BoundaryElement % BoundaryInfo % Left
  ELSE ! we are dealing with a body-body boundary and asume that the normal is pointing outwards
     ParentElement => BoundaryElement % BoundaryInfo % Right
     IF (ParentElement % BodyId == other_body_id) ParentElement => BoundaryElement % BoundaryInfo % Left
  END IF

! just to be on the save side, check again
!-----------------------------------------
  IF ( .NOT. ASSOCIATED(ParentElement) ) THEN
     WRITE(Message,'(A,I10,A)')&
          'Parent Element for Boundary element no. ',&
          BoundaryElement % ElementIndex, ' not found'
     !CALL FATAL('getSurfaceVelocity',Message)
  END IF

  VeloSol => VariableGet( Model % Variables, "vydata" )
  IF ( ASSOCIATED( VeloSol ) ) THEN
       VeloPerm => VeloSol % Perm
       Velo     => VeloSol % Values
  ELSE
       WRITE(Message, '(A)') 'Could not find vydata variable '
       CALL FATAL('getSurfaceVelocity',Message)
  END IF

  N = Model % Mesh % MaxElementNodes

  ALLOCATE(BoundaryElementNodes % x( N ), BoundaryElementNodes % y( N ), BoundaryElementNodes % z( N ), STAT=istat)
  IF ( istat /= 0 ) THEN
      CALL FATAL('getSurfaceVelocity', 'Memory allocation error' )
  END IF

  ALLOCATE(BDElementNodes % x( N ), BDElementNodes % y( N ), BDElementNodes % z( N ), STAT=istat)
  IF ( istat /= 0 ) THEN
      CALL FATAL('getSurfaceVelocity', 'Memory allocation error' )
  END IF


  CALL GetElementNodes( BoundaryElementNodes, BoundaryElement )
  nb1 =  GetElementNOFNodes(BoundaryElement)

  foundNodes = .FALSE.
  DO t=1,  Model % Mesh % NumberOfBoundaryElements
     BDElement => GetBoundaryElement(t)
     IF ( .NOT.ActiveBoundaryElement(BDElement) ) CYCLE
     IF ( GetElementFamily(BDElement) == 1 ) CYCLE
     bc_id = GetBCId( BDElement )
     IF (bc_id == 6) THEN
         nb2 = GetElementNOFNodes(BDElement)
         CALL GetElementNodes( BDElementNodes, BDElement )
         DO i=1,nb2
           DO j=1, nb1
              IF ((BoundaryElementNodes % x(j) == BDElementNodes % x(i) .AND. &
                  BoundaryElementNodes % y(j) == BDElementNodes % y(i) ) .AND. BoundaryElement % NodeIndexes(j) == nodenumber) THEN
                  velocity = Velo(VeloPerm(BDElement % NodeIndexes(i)))
                  foundNodes = .True.
                  EXIT
              END IF
           END DO
           IF (foundNodes) EXIT
         END DO
         IF (foundNodes) EXIT
     END IF
  END DO

  DEALLOCATE(BoundaryElementNodes % x, BoundaryElementNodes % y, BoundaryElementNodes % z)
  DEALLOCATE(BDElementNodes % x, BDElementNodes % y, BDElementNodes % z)

END FUNCTION getSurfaceVelocity2

