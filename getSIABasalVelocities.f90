FUNCTION getSIABasalVelocitiesX (Model, nodenumber, depth) RESULT(Bvelo)

   USE types
   USE CoordinateSystems
   USE SolverUtils
   USE ElementDescription
   USE DefUtils

   IMPLICIT NONE

   TYPE(Model_t) :: Model
   REAL (KIND=dp) :: depth , x              
   INTEGER :: nodenumber

   TYPE(Element_t), POINTER :: BoundaryElement, ParentElement
   TYPE(ValueList_t), POINTER :: BC, ParentMaterial, BodyForce
   TYPE(Variable_t), POINTER :: FreeSurfGrad1, FreeSurfGrad2, Beta
   REAL(KIND=dp), POINTER :: FreeSurfGrad1Values(:), FreeSurfGrad2Values(:), BetaValues(:)
   INTEGER, POINTER :: FreeSurfGrad1Perm(:), FreeSurfGrad2Perm(:), BetaPerm(:)
   INTEGER :: DIM, i, j
   INTEGER :: NBoundary, NParent, BoundaryElementNode
   INTEGER :: material_id, body_id, other_body_id, bf_id
   REAL (KIND=dp) :: C, Bvelo, scaleFactor, gradh, dragBeta
   REAL (KIND=dp) :: rho, gravity, homoTemp
   REAL (KIND=dp) :: C_b, gamma, exponent_p, exponent_q
   REAL (KIND=dp) :: tau_b
   LOGICAL :: GotIt, FirstTime = .TRUE., scaling = .FALSE.

   SAVE :: DIM, FirstTime, scaling, scaleFactor

   BC => GetBC(Model % CurrentElement)

   IF (FirstTime) THEN
        FirstTime = .FALSE.  
        DIM = CoordinateSystemDimension()

        CALL INFO('getSIABasalVelocitiesX', & 
                'Slip condition will be calculated on bedrock', level=3)        
        
        scaling = GetLogical(BC, 'SIA basal sliding scaling', GotIt)
        IF(GotIt .AND. scaling) THEN

            scaleFactor = GetConstReal( BC, 'SIA basal scaling factor', GotIt )       
            IF (.NOT. GotIt) THEN
                CALL FATAL('getSIABasalVelocitiesX', 'SIA basal scaling factor not found')
            END IF
 
            WRITE(Message,'(a,F8.2)') 'SIA: Scaling of the basal sliding by: ', scaleFactor   
            CALL INFO('getSIABasalVelocitiesX', Message, level=3)   

        END IF
    
   END IF

   !----------------! get some information upon active boundary element and its parent  !----------------------------------------

   BoundaryElement => Model % CurrentElement

   !-----------------------Do nothing if this is a halo-element of a parallel mesh?  !------------------------------------------

   IF (ParEnv % myPe .NE. BoundaryElement % partIndex) RETURN

   IF ( .NOT. ASSOCIATED(BoundaryElement) ) THEN     
        CALL FATAL('getSIABasalVelocities','No boundary element found')  
   END IF

   NBoundary = BoundaryElement % Type % NumberOfNodes  
   DO BoundaryElementNode=1,NBoundary     
    IF (nodenumber .EQ. BoundaryElement % NodeIndexes(BoundaryElementNode)) THEN        
        GotIt = .TRUE.        
        EXIT     
        END IF  
   END DO

   IF (.NOT.GotIt) THEN
        CALL WARN('getSIABasalVelocities','Node not found in Current Element')     
        Bvelo = 0.0D00
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
        CALL FATAL('getSIABasalVelocitiesX',Message)
   END IF

   body_id = ParentElement % BodyId
   material_id = ListGetInteger(Model % Bodies(body_id) % Values, 'Material', GotIt)
   ParentMaterial => Model % Materials(material_id) % Values
   IF ((.NOT. ASSOCIATED(ParentMaterial)) .OR. (.NOT. GotIt)) THEN
        WRITE(Message,'(A,I10,A,I10)')&
          'No material values found for body no ', body_id,&
          ' under material id ', material_id
        CALL FATAL('getSIABasalVelocitiesX',Message)
   END IF

   !Get BofyForce ID
   !----------------
    bf_id = ListGetInteger( Model % Bodies(body_id) % Values, &
              'Body Force', gotIt, 1, Model % NumberOfBodyForces )

  
  FreeSurfGrad1 => VariableGet( Model % Variables, 'FreeSurfGrad1' )
  IF ( ASSOCIATED( FreeSurfGrad1 ) ) THEN
        FreeSurfGrad1Perm    =>  FreeSurfGrad1 % Perm
        FreeSurfGrad1Values  =>  FreeSurfGrad1 % Values
  ELSE
        CALL FATAL('getSIABasalVelocitiesX','No variable FreeSurfGrad1 found')  
  END IF

  FreeSurfGrad2 => VariableGet( Model % Variables, 'FreeSurfGrad2' )
  IF ( ASSOCIATED( FreeSurfGrad2 ) ) THEN
        FreeSurfGrad2Perm    =>  FreeSurfGrad2 % Perm
        FreeSurfGrad2Values  =>  FreeSurfGrad2 % Values
  ELSE
        CALL FATAL('getSIABasalVelocitiesX','No variable FreeSurfGrad2 found')  
  END IF

  Beta => VariableGet( Model % Variables, 'Beta' )
  IF ( ASSOCIATED( Beta ) ) THEN
        BetaPerm    =>  Beta % Perm
        BetaValues  =>  Beta % Values
  ELSE
        CALL FATAL('getSIABasalVelocitiesX','No variable Beta found')  
  END IF

  rho = GetConstReal(ParentMaterial, 'Density', GotIt)
  IF(.NOT. GotIt) THEN 
        CALL FATAL('getSIABasalVelocitiesX', 'Density not found in Material')
  END IF

  gravity = GetConstReal(Model % BodyForces(bf_id) % Values, 'Flow Bodyforce 3', GotIt)
  IF(.NOT. GotIt) THEN 
        CALL FATAL('getSIABasalVelocitiesX', 'Flow Body Force 3 not found')
  END IF

  gravity = gravity*(-1.0_dp)

  gradh =  ( FreeSurfGrad1Values(FreeSurfGrad1Perm(nodenumber))**2.0_dp &
            + FreeSurfGrad2Values(FreeSurfGrad2Perm(nodenumber))**2.0_dp )**0.5_dp

  dragBeta = 10.0_dp**BetaValues(BetaPerm(nodenumber))
  IF (scaling) THEN
      dragBeta = dragBeta / scaleFactor
  END IF

  Bvelo = -( (rho * gravity * depth)*FreeSurfGrad1Values(FreeSurfGrad1Perm(nodenumber)) ) / dragBeta

END FUNCTION getSIABasalVelocitiesX

!---------------------------------------------------------------------------------------------------------------------------------

FUNCTION getSIABasalVelocitiesY (Model, nodenumber, depth) RESULT(Bvelo)

   USE types
   USE CoordinateSystems
   USE SolverUtils
   USE ElementDescription
   USE DefUtils

   IMPLICIT NONE

   TYPE(Model_t) :: Model
   REAL (KIND=dp) :: depth , x              
   INTEGER :: nodenumber

   TYPE(Element_t), POINTER :: BoundaryElement, ParentElement
   TYPE(ValueList_t), POINTER :: BC, ParentMaterial, BodyForce
   TYPE(Variable_t), POINTER :: FreeSurfGrad1, FreeSurfGrad2, Beta
   REAL(KIND=dp), POINTER :: FreeSurfGrad1Values(:), FreeSurfGrad2Values(:), BetaValues(:)
   INTEGER, POINTER :: FreeSurfGrad1Perm(:), FreeSurfGrad2Perm(:), BetaPerm(:)
   INTEGER :: DIM, i, j
   INTEGER :: NBoundary, NParent, BoundaryElementNode
   INTEGER :: material_id, body_id, other_body_id, bf_id
   REAL (KIND=dp) :: C, Bvelo, scaleFactor, gradh, dragBeta
   REAL (KIND=dp) :: rho, gravity, homoTemp
   REAL (KIND=dp) :: C_b, gamma, exponent_p, exponent_q
   REAL (KIND=dp) :: tau_b
   LOGICAL :: GotIt, FirstTime = .TRUE., scaling = .FALSE.

   SAVE :: DIM, FirstTime, scaling, scaleFactor

   BC => GetBC(Model % CurrentElement)

   IF (FirstTime) THEN
        FirstTime = .FALSE.  
        DIM = CoordinateSystemDimension()

        CALL INFO('getSIABasalVelocitiesY', & 
                'Slip condition will be calculated on bedrock', level=3)        
        
        scaling = GetLogical(BC, 'SIA basal sliding scaling', GotIt)
        IF(GotIt .AND. scaling) THEN

            scaleFactor = GetConstReal( BC, 'SIA basal scaling factor', GotIt )       
            IF (.NOT. GotIt) THEN
                CALL FATAL('getSIABasalVelocitiesY', 'SIA basal scaling factor not found')
            END IF
 
            WRITE(Message,'(a,F8.2)') 'SIA: Scaling of the basal sliding by: ', scaleFactor   
            CALL INFO('getSIABasalVelocitiesY', Message, level=3)   

        END IF
    
   END IF

   !----------------! get some information upon active boundary element and its parent  !----------------------------------------

   BoundaryElement => Model % CurrentElement

   !-----------------------Do nothing if this is a halo-element of a parallel mesh?  !------------------------------------------

   IF (ParEnv % myPe .NE. BoundaryElement % partIndex) RETURN

   IF ( .NOT. ASSOCIATED(BoundaryElement) ) THEN     
        CALL FATAL('getSIABasalVelocitiesY','No boundary element found')  
   END IF

   NBoundary = BoundaryElement % Type % NumberOfNodes  
   DO BoundaryElementNode=1,NBoundary     
    IF (nodenumber .EQ. BoundaryElement % NodeIndexes(BoundaryElementNode)) THEN        
        GotIt = .TRUE.        
        EXIT     
        END IF  
   END DO

   IF (.NOT.GotIt) THEN
        CALL WARN('getSIABasalVelocitiesY','Node not found in Current Element')     
        Bvelo = 0.0D00
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
        CALL FATAL('getSIABasalVelocitiesY',Message)
   END IF

   body_id = ParentElement % BodyId
   material_id = ListGetInteger(Model % Bodies(body_id) % Values, 'Material', GotIt)
   ParentMaterial => Model % Materials(material_id) % Values
   IF ((.NOT. ASSOCIATED(ParentMaterial)) .OR. (.NOT. GotIt)) THEN
        WRITE(Message,'(A,I10,A,I10)')&
          'No material values found for body no ', body_id,&
          ' under material id ', material_id
        CALL FATAL('getSIABasalVelocitiesY',Message)
   END IF

   !Get BofyForce ID
   !----------------
    bf_id = ListGetInteger( Model % Bodies(body_id) % Values, &
              'Body Force', gotIt, 1, Model % NumberOfBodyForces )

  
  FreeSurfGrad1 => VariableGet( Model % Variables, 'FreeSurfGrad1' )
  IF ( ASSOCIATED( FreeSurfGrad1 ) ) THEN
        FreeSurfGrad1Perm    =>  FreeSurfGrad1 % Perm
        FreeSurfGrad1Values  =>  FreeSurfGrad1 % Values
  ELSE
        CALL FATAL('getSIABasalVelocitiesY','No variable FreeSurfGrad1 found')  
  END IF

  FreeSurfGrad2 => VariableGet( Model % Variables, 'FreeSurfGrad2' )
  IF ( ASSOCIATED( FreeSurfGrad2 ) ) THEN
        FreeSurfGrad2Perm    =>  FreeSurfGrad2 % Perm
        FreeSurfGrad2Values  =>  FreeSurfGrad2 % Values
  ELSE
        CALL FATAL('getSIABasalVelocitiesY','No variable FreeSurfGrad2 found')  
  END IF

  Beta => VariableGet( Model % Variables, 'Beta' )
  IF ( ASSOCIATED( Beta ) ) THEN
        BetaPerm    =>  Beta % Perm
        BetaValues  =>  Beta % Values
  ELSE
        CALL FATAL('getSIABasalVelocitiesY','No variable Beta found')  
  END IF

  rho = GetConstReal(ParentMaterial, 'Density', GotIt)
  IF(.NOT. GotIt) THEN 
        CALL FATAL('getSIABasalVelocitiesY', 'Density not found in Material')
  END IF

  gravity = GetConstReal(Model % BodyForces(bf_id) % Values, 'Flow Bodyforce 3', GotIt)
  IF(.NOT. GotIt) THEN 
        CALL FATAL('getSIABasalVelocitiesY', 'Flow Body Force 3 not found')
  END IF

  gravity = gravity*(-1.0_dp)

  gradh =  ( FreeSurfGrad1Values(FreeSurfGrad1Perm(nodenumber))**2.0_dp &
            + FreeSurfGrad2Values(FreeSurfGrad2Perm(nodenumber))**2.0_dp )**0.5_dp

  dragBeta =  10.0_dp**BetaValues(BetaPerm(nodenumber))
  IF (scaling) THEN
      dragBeta = dragBeta / scaleFactor
  END IF

  Bvelo = -( (rho * gravity * depth)*FreeSurfGrad2Values(FreeSurfGrad2Perm(nodenumber)) ) / dragBeta

END FUNCTION getSIABasalVelocitiesY

!---------------------------------------------------------------------------------------------------------------------------------

!**************************************************************************
!*
!*  basal melting velocity  as a function of effectiveheat flux 
!* (provided externaly)
!*
!**************************************************************************
FUNCTION getBasalMeltingVelocity( Model, Node, HeatFlux ) RESULT(basalMeltingvelocity)

!-----------------------------------------------------------
  USE DefUtils
  USE SolverUtils
!-----------------------------------------------------------
  IMPLICIT NONE
!-----------------------------------------------------------
  !external variables
  TYPE(Model_t), TARGET :: Model
  INTEGER :: Node
  REAL(KIND=dp) :: HeatFlux, basalMeltingvelocity

  !internal variables
  TYPE(ValueList_t), POINTER :: ParentMaterial, BC
  TYPE(Element_t), POINTER :: BoundaryElement, ParentElement 
  TYPE(Variable_t), POINTER :: VarTemp, VarResidual, FlowVariable, NormalVar
  TYPE(Variable_t), POINTER :: basalMeltingSol
  TYPE(Variable_t), POINTER :: FreeSurfGrad1, FreeSurfGrad2, DepthSol
  REAL(KIND=dp), POINTER :: FreeSurfGrad1Values(:), FreeSurfGrad2Values(:), DepthValues(:)
  INTEGER, POINTER :: FreeSurfGrad1Perm(:), FreeSurfGrad2Perm(:), DepthPerm(:)
  REAL(KIND=dp), POINTER :: basalMeltingValues(:), FlowValues(:)
  REAL(KIND=dp), POINTER :: NormalValues(:)
  INTEGER, POINTER :: basalMeltingPerm(:), FlowPerm(:), NormalPerm(:)
  INTEGER :: basalMeltingDOFs
  INTEGER :: i, j, N, DIM, NBoundary, NParent, BoundaryElementNode, Ind(3,3), &
       ParentElementNode, body_id, other_body_id, material_id, equation_id, bf_id, istat
  REAL(KIND=dp), ALLOCATABLE :: LatentHeat(:), Density(:),ExternalHF(:),PressureMeltingPoint(:)
  REAL(KIND=dp), ALLOCATABLE :: Sig(:,:), normal(:), velo(:), Sn(:), ut(:)
  REAL (KIND=dp) :: rho, gravity, gradh  
  REAL (KIND=dp) :: HomTemp, BasalMeltingOffset, vSn, un, limitValue
  CHARACTER(LEN=MAX_NAME_LEN) :: TempName, SolverName
  LOGICAL ::  FirstTime = .TRUE., GotIt, stat, isSliding, isLimited

  SAVE FirstTime, LatentHeat, Density, DIM, N, ExternalHF, PressureMeltingPoint, Sig, &
  normal, velo, Ind, Sn, ut

  WRITE(SolverName, '(A)') 'IceFlowProperties (getBasalMeltingVelocity))'

  !--------------------------------
  ! Allocations
  !--------------------------------
  IF ( FirstTime .OR. Model % Mesh % Changed) THEN
     DIM = CoordinateSystemDimension()
     N = Model % MaxElementNodes 

     IF ((DIM == 2).OR.(DIM == 3))  THEN
         ALLOCATE( Sig(DIM,DIM),normal(DIM), velo(DIM), Sn(DIM), ut(DIM) )
     ELSE
         CALL FATAL('IceFlowPoperties(getBasalMeltingVelocity)', 'Bad dimension of the problem')
     END IF
      
     Do i=1, 3
        Ind(i,i) = i
     END DO
     Ind(1,2) = 4
     Ind(2,1) = 4
     Ind(2,3) = 5
     Ind(3,2) = 5
     Ind(3,1) = 6
     Ind(1,3) = 6
     CALL INFO('IceFlowPoperties(getBasalMeltingVelocity)', 'Melting condition will be calculated on bedrock', level=3)

      IF (.NOT.FirstTime) THEN
         DEALLOCATE( &
              LatentHeat,&        
              Density,&
              ExternalHF,&
			  PressureMeltingPoint)
      END IF
      ALLOCATE( &
           LatentHeat( N ),&
           Density( N ),&
           ExternalHF( N ),&
		   PressureMeltingPoint( N ),&
           STAT = istat)
      IF (istat /= 0) THEN
         CALL FATAL(SolverName,'Allocations failed')
      ELSE 
         CALL INFO(SolverName,'Allocations done',Level=1)
      END IF

      FirstTime = .FALSE.
   END IF

   basalMeltingSol => VariableGet( Model % Variables, 'BasalMelting')
   IF ( ASSOCIATED(basalMeltingSol) ) THEN
      basalMeltingPerm => basalMeltingSol % Perm
      basalMeltingValues => basalMeltingSol % Values
      basalMeltingDOFs = basalMeltingSol % DOFs
   ELSE
      CALL INFO(' IceFlowPoperties', 'No variable associated for metling velocity')
   END IF

  !-----------------------------------------------------------------
  ! get some information upon active boundary element and its parent
  !-----------------------------------------------------------------
  BoundaryElement => Model % CurrentElement
  !---------------------------------------------------------
  ! Do nothing if this is a halo-element of a parallel mesh?
  !---------------------------------------------------------
  IF (ParEnv % myPe .NE. BoundaryElement % partIndex) RETURN


  IF ( .NOT. ASSOCIATED(BoundaryElement) ) THEN
     CALL FATAL(SolverName,'No boundary element found')
  END IF  
  NBoundary = BoundaryElement % Type % NumberOfNodes
  DO BoundaryElementNode=1,NBoundary
     IF (Node .EQ. BoundaryElement % NodeIndexes(BoundaryElementNode)) THEN
        GotIt = .TRUE.
        EXIT
     END IF
  END DO
  IF (.NOT.GotIt) THEN
     CALL WARN(SolverName,'Node not found in Current Element')
     basalMeltingVelocity = 0.0D00
  IF ( ASSOCIATED(basalMeltingSol) ) THEN
     basalMeltingValues(basalMeltingDOFs*(basalMeltingPerm(Node)-1)+1) = basalMeltingVelocity
  END IF
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
  IF ( .NOT. ASSOCIATED(ParentElement) ) THEN
     WRITE(Message,'(A,I10,A)')&
          'Parent Element for Boundary element no. ',&
          BoundaryElement % ElementIndex, ' not found'
     CALL FATAL(SolverName,Message)
  END IF  

  body_id = ParentElement % BodyId
  material_id = ListGetInteger(Model % Bodies(body_id) % Values, 'Material', GotIt)
  ParentMaterial => Model % Materials(material_id) % Values
  IF ((.NOT. ASSOCIATED(ParentMaterial)) .OR. (.NOT. GotIt)) THEN
     WRITE(Message,'(A,I10,A,I10)')&
          'No material values found for body no ', body_id,&
          ' under material id ', material_id
     CALL FATAL(SolverName,Message)
  END IF 
  NParent = ParentElement % Type % NumberOfNodes
  DO ParentElementNode=1,NParent
     IF ( Node == ParentElement % NodeIndexes(ParentElementNode) ) EXIT
  END DO

   !Get BofyForce ID
   !----------------
    bf_id = ListGetInteger( Model % Bodies(body_id) % Values, &
              'Body Force', gotIt, 1, Model % NumberOfBodyForces )

  !-----------------------------------------------------------------------
  ! get temperature variable
  !-----------------------------------------------------------------------
  TempName =  GetString(ParentMaterial ,'Temperature Name', GotIt)
  IF (.NOT.GotIt) THEN
     WRITE(Message,'(A, I10)') 'No Temperature Name found for parent material no. ', material_id
     CALL FATAL(SolverName, Message)
  ELSE
     WRITE(Message,'(a,a)') 'Variable Name for temperature: ', TempName
     CALL INFO(SolverName,Message,Level=21)
  END IF   
  VarTemp => VariableGet( Model % Variables, TRIM(TempName), .TRUE. )
  IF ( .NOT.ASSOCIATED( VarTemp ) ) THEN
     WRITE(Message,'(A,A,A)') 'Variable >', TRIM(TempName) , '< not found'
     CALL FATAL(SolverName,Message)
  END IF
  !-----------------------------------------------------------------------
  ! get residual variable (if 0 indicates no melt, if /= 0 melting)
  !-----------------------------------------------------------------------
  VarResidual  => VariableGet( Model % Variables, TRIM(TempName) // ' Residual', .TRUE. )
  IF ( .NOT.ASSOCIATED( VarResidual ) ) THEN
     WRITE(Message,'(A,A,A)') 'Variable >', TRIM(TempName) // ' Residual', '< not found'
     CALL FATAL(SolverName,Message)
  END IF
  !-------------------------
  ! Get material parameters
  !-------------------------
  Model % CurrentElement => ParentElement
  LatentHeat(1:NParent) = ListGetReal(ParentMaterial, 'Latent Heat', NParent, ParentElement % NodeIndexes, GotIt)
  IF (.NOT. GotIt) THEN
     CALL FATAL(SolverName,'No value for >Latent Heat< found')
  END IF
  Density(1:NParent) = ListGetReal( ParentMaterial, 'Density', NParent, ParentElement % NodeIndexes)
  IF (.NOT. GotIt) THEN
     CALL FATAL(SolverName,'No value for >Density< found')
  END IF
  
  PressureMeltingPoint(1:NParent)  = ListGetReal(ParentMaterial, TRIM(TempName) // ' Upper Limit', &
           NParent, ParentElement % NodeIndexes, GotIt)
  IF(.NOT. GotIt) THEN
      CALL FATAL(SolverName, 'Temp Upper Limit not found in Material')
  END IF

  !-------------------------
  ! Get boundary parameters
  !-------------------------
  Model % CurrentElement => BoundaryElement
  BC => GetBC()
  IF (.NOT.ASSOCIATED(BC)) THEN
     CALL FATAL(SolverName,'No Boundary Condition associated')
  ELSE
     ExternalHF(1:NBoundary) = GetReal(BC, TRIM(TempName) // ' Heat Flux', GotIt)
     IF (.NOT. GotIt) THEN
        WRITE(Message,'(a,a,a)') 'Keyword >', TRIM(TempName) // ' Heat Flux','< not found'
        CALL WARN(SolverName,Message)
        ExternalHF(1:NBoundary) = 0.0D00
     END IF
  END IF

  BasalMeltingOffset = GetConstReal(BC, 'Basal Melting Offset', GotIt)
  IF (GotIt) THEN
     WRITE(Message,'(A,F8.2)') 'Setting Basal melting offset to', BasalMeltingOffset
     CALL INFO(SolverName,Message,level=10)
  ELSE
     BasalMeltingOffset = 0.0D00
  END IF

  isSliding = GetLogical(BC, 'Sliding in melting equation', GotIt)
  IF (.NOT. GotIt) THEN
     WRITE(Message,'(A)') 'No basal sliding in basal melting equation found'
     CALL INFO(SolverName,Message,level=12)
     isSliding = .FALSE. 
  END IF
  
  isLimited = GetLogical(BC, 'Limit basal melting', GotIt)
  IF (.NOT. GotIt) THEN
     isLimited = .FALSE. 
  ELSE IF (GotIt .AND. isLimited) THEN
	 limitValue = GetConstReal(BC, 'Basal melting limit value', GotIt)
	 IF (.NOT. GotIt) THEN
	    CALL FATAL(SolverName, 'Limitation but cant find value to use')
	 END IF
  END IF

IF(isSliding) THEN

   equation_id =  ListGetInteger( Model % Bodies(body_id) % Values, 'Equation', &
                   minv=1, maxv=Model % NumberOFEquations )

  ! Get the variable velocity
  !---------------------------
  FlowVariable => VariableGet( Model % Variables, 'SIAFlow' )
  IF ( ASSOCIATED( FlowVariable ) ) THEN
      FlowPerm    => FlowVariable % Perm
      FlowValues  => FlowVariable % Values
  ELSE
      CALL FATAL('IceFlowPoperties(getBasalMeltingVelocity)', 'Need NS Solver, Flow Solution not found')
  END IF

  FreeSurfGrad1 => VariableGet( Model % Variables, 'FreeSurfGrad1' )
  IF ( ASSOCIATED( FreeSurfGrad1 ) ) THEN
        FreeSurfGrad1Perm    =>  FreeSurfGrad1 % Perm
        FreeSurfGrad1Values  =>  FreeSurfGrad1 % Values
  ELSE
        CALL FATAL('getSIABasalVelocities','No variable FreeSurfGrad1 found')  
  END IF

  FreeSurfGrad2 => VariableGet( Model % Variables, 'FreeSurfGrad2' )
  IF ( ASSOCIATED( FreeSurfGrad2 ) ) THEN
        FreeSurfGrad2Perm    =>  FreeSurfGrad2 % Perm
        FreeSurfGrad2Values  =>  FreeSurfGrad2 % Values
  ELSE
        CALL FATAL('getSIABasalVelocities','No variable FreeSurfGrad2 found')  
  END IF

  DepthSol => VariableGet( Model % Variables, 'Depth' )
  IF ( ASSOCIATED( DepthSol ) ) THEN
        DepthPerm    =>  DepthSol % Perm
        DepthValues  =>  DepthSol % Values
  ELSE
        CALL FATAL('getSIABasalVelocities','No variable Depth found')  
  END IF


  rho = GetConstReal(ParentMaterial, 'Density', GotIt)
  IF(.NOT. GotIt) THEN 
        CALL FATAL('getSIABasalVelocities', 'Density not found in Material')
  END IF

  gravity = GetConstReal(Model % BodyForces(bf_id) % Values, 'Flow Bodyforce 3', GotIt)
  IF(.NOT. GotIt) THEN 
        CALL FATAL('getSIABasalVelocities', 'Flow Body Force 3 not found')
  END IF

  gravity = gravity*(-1.0_dp)

  gradh =  ( FreeSurfGrad1Values(FreeSurfGrad1Perm(Node))**2.0_dp &
            + FreeSurfGrad2Values(FreeSurfGrad2Perm(Node))**2.0_dp )**0.5_dp

  ! Get the variable for normal vector
  NormalVar =>  VariableGet(Model % Variables,'Normal Vector')
  IF ( ASSOCIATED( NormalVar ) ) THEN
      NormalPerm => NormalVar % Perm
      NormalValues => NormalVar % Values
  ELSE
      CALL FATAL('IceFlowPoperties(getBasalMeltingVelocity)', 'Normal Vector variable not found')
  END IF

  DO i=1, DIM
     normal(i) = NormalValues(DIM*(NormalPerm(Node)-1) + i)      
     velo(i) = FlowValues( (DIM+1)*(FlowPerm(Node)-1) + i )
  END DO
  
  !Tengantial velocity
  un = SUM(velo(1:DIM)*(normal(1:DIM))) 
  ut = velo(1:DIM)-un*(normal(1:DIM))

  Sig = 0.0_dp
  
  Sig(1,1) = -rho * gravity * DepthValues(DepthPerm(Node))
  Sig(1,2) = 0.0_dp
  Sig(1,3) = -rho * gravity * DepthValues(DepthPerm(Node)) * FreeSurfGrad1Values(FreeSurfGrad1Perm(Node))
  Sig(2,1) = 0.0_dp
  Sig(2,2) = -rho * gravity * DepthValues(DepthPerm(Node))
  Sig(2,3) = -rho * gravity * DepthValues(DepthPerm(Node)) * FreeSurfGrad2Values(FreeSurfGrad2Perm(Node))
  Sig(3,1) = -rho * gravity * DepthValues(DepthPerm(Node)) * FreeSurfGrad1Values(FreeSurfGrad1Perm(Node))
  Sig(3,2) = -rho * gravity * DepthValues(DepthPerm(Node)) * FreeSurfGrad2Values(FreeSurfGrad2Perm(Node))
  Sig(3,3) = -rho * gravity * DepthValues(DepthPerm(Node))

  DO i=1, DIM
     Sn(i) = SUM(Sig(i,1:DIM)*normal(1:DIM)) 
  END DO  

  vSn = SUM(ut(1:DIM)*Sn(1:DIM)) 

  !WRITE(*,*) vSn
  
END IF

  !------------------------------------
  ! Get basal melting rate and velocity
  !------------------------------------
! IF (VarResidual % Values(VarResidual % Perm(Node)) < 0.0D00) THEN
  If (VarTemp % Values(VarTemp % Perm(Node)) >= PressureMeltingPoint(BoundaryElementNode)) THEN
     IF(isSliding) THEN
       basalMeltingvelocity = MAX((ExternalHF(BoundaryElementNode)+HeatFlux-MIN(vSn,0.0d0)) &
          /(LatentHeat(BoundaryElementNode) * Density(BoundaryElementNode)), &
          BasalMeltingOffset)
	   IF (isLimited) THEN
	       IF (basalMeltingvelocity > limitValue) THEN
		       basalMeltingvelocity = limitValue
		   END IF
	   END IF
     ELSE
       basalMeltingvelocity = MAX((ExternalHF(BoundaryElementNode)+HeatFlux) &
          /(LatentHeat(BoundaryElementNode) * Density(BoundaryElementNode)), &
          BasalMeltingOffset)
	   IF (isLimited) THEN
	       IF (basalMeltingvelocity > limitValue) THEN
		       basalMeltingvelocity = limitValue
		   END IF
	   END IF
     END IF
 !    PRINT *, Node, basalMeltingvelocity, ExternalHF(BoundaryElementNode), HeatFlux
     IF ( ASSOCIATED(basalMeltingSol) ) THEN   
        basalMeltingValues(basalMeltingDOFs*(basalMeltingPerm(Node)-1)+1) = & 
        basalMeltingVelocity
     END IF
  ELSE
     basalMeltingvelocity = BasalMeltingOffset
     IF ( ASSOCIATED(basalMeltingSol) ) THEN
        basalMeltingValues(basalMeltingDOFs*(basalMeltingPerm(Node)-1)+1) = & 
        BasalMeltingOffset
   END IF
  END IF
END FUNCTION getBasalMeltingVelocity

