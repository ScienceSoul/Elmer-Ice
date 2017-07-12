FUNCTION getOceanWaterPressure( Model, nodenumber, elevation) RESULT(pressure)
  USE Types
  USE CoordinateSystems
  USE SolverUtils
  USE ElementDescription
  USE DefUtils
 !-------------------------------------------------------------------
 IMPLICIT NONE
 !-------------------------external variables-----------------------
 TYPE(Model_t) :: Model 
 TYPE(Element_t), POINTER :: BoundaryElement, ParentElement
 TYPE(ValueList_t), POINTER :: BC, ParentMaterial, BodyForce
 INTEGER :: material_id, body_id, other_body_id, bf_id
 INTEGER :: nodenumber, NBoundary, BoundaryElementNode
 REAL(KIND=dp) :: elevation, rhow, gravity, pressure
 LOGICAL :: GotIt

 !----------------! get some information upon active boundary element and its parent  !----------------------------------------

 BoundaryElement => Model % CurrentElement

 !-----------------------Do nothing if this is a halo-element of a parallel mesh?  !------------------------------------------

 IF (ParEnv % myPe .NE. BoundaryElement % partIndex) RETURN

 IF ( .NOT. ASSOCIATED(BoundaryElement) ) THEN     
 CALL FATAL('getOceanWaterPressure','No boundary element found')  
 END IF

 NBoundary = BoundaryElement % Type % NumberOfNodes  
 DO BoundaryElementNode=1,NBoundary     
   IF (nodenumber .EQ. BoundaryElement % NodeIndexes(BoundaryElementNode)) THEN        
      GotIt = .TRUE.        
      EXIT     
   END IF  
 END DO

IF (.NOT.GotIt) THEN
   CALL WARN('getOceanWaterPressure','Node not found in Current Element')     
   pressure = 0.0D00
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
     CALL FATAL('getOceanWaterPressure',Message)
 END IF

 body_id = ParentElement % BodyId
 material_id = ListGetInteger(Model % Bodies(body_id) % Values, 'Material', GotIt)
 ParentMaterial => Model % Materials(material_id) % Values
 IF ((.NOT. ASSOCIATED(ParentMaterial)) .OR. (.NOT. GotIt)) THEN
     WRITE(Message,'(A,I10,A,I10)')&
          'No material values found for body no ', body_id,&
          ' under material id ', material_id
     CALL FATAL('getOceanWaterPressure',Message)
 END IF

 !Get BofyForce ID
 !----------------
 bf_id = ListGetInteger( Model % Bodies(body_id) % Values, &
              'Body Force', gotIt, 1, Model % NumberOfBodyForces )

 rhow = GetConstReal(ParentMaterial, 'Sea water density', GotIt)
 IF(.NOT. GotIt) THEN 
    CALL FATAL('getOceanWaterPressure', 'Density not found in Material')
 END IF

 gravity = GetConstReal(Model % BodyForces(bf_id) % Values, 'Flow Bodyforce 3', GotIt)
 IF(.NOT. GotIt) THEN 
    CALL FATAL('USF_Sliding(Sliding_Weertman)', 'Flow Body Force 3 not found')
 END IF

 !We need a positive value for the gravity. It was orginally negative as it was acquired from the body force value in the model setting.
 gravity = gravity*(-1.0)

 pressure = -MAX(gravity * rhow * (0.0 - elevation), 0.0)
 !write(*,*) pressure

END FUNCTION getOceanWaterPressure
