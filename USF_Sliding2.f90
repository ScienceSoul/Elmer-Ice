! *********************************************************************************************/
! *
! *                USF_Sliding.f90
! *
! *         Gives the basal drag for different sliding law
! *
! *   (1) Sliding_Weertman
! *   Need some inputs in the sif file.
! *   Parameters: Sliding Coefficient      -> C 
!             Sliding Exponent         -> m
!             Sliding Linear Velocity  -> ut0
!
! *   Compute the Bdrag coefficient such as tau_b = Bdrag ub
! *   for the non-linear Weertman law tau_b = C ub^m
! *   To linearize the Weertman law, we can distinguish 4 cases:
! *     1/ ut=0 , tau_b=0     => C= infinity (no sliding, first step)
! *     2/ ut=0 , tau_b =/0   => C=C^1/m tau_b^(1-1/m)
! *     3/ ut=/0 , tau_b=0    => C=Cub^(m-1)
! *     4/ ut=/0 , tau_b=/0   => case 3 
! *   For cases 3 and 4, if ut < ut0, Bdrag = C ut0^{m-1}
! *
! ********************************************************************************************/
! *
! *  2008 Hakime Seddik
! *  2007/10/25. Gael Durand
! *  2008/04/06 OG 2D -> 3D
! *
! *********************************************************************************************/


!*********************************************************************************************/
FUNCTION Sliding_Weertman (Model, nodenumber, y) RESULT(Bdrag)


   USE types
   USE CoordinateSystems
   USE SolverUtils
   USE ElementDescription
   USE DefUtils

   IMPLICIT NONE

   TYPE(Model_t) :: Model
   REAL (KIND=dp) :: y , x              
   INTEGER :: nodenumber

   TYPE(Element_t), POINTER :: BoundaryElement, ParentElement
   TYPE(ValueList_t), POINTER :: BC, ParentMaterial, BodyForce
   TYPE(Variable_t), POINTER :: StressVariable, NormalVar, FlowVariable, TempHomologous, Temp, DepthVariable
   REAL(KIND=dp), POINTER :: StressValues(:), NormalValues(:), FlowValues(:), TempHomologousValues(:), & 
                             TempValues(:), DepthValues(:)
   INTEGER, POINTER :: StressPerm(:), NormalPerm(:), FlowPerm(:), TempHomologousPerm(:), TempPerm(:), DepthPerm(:)
   INTEGER :: DIM, i, j, Ind(3,3)
   INTEGER :: NBoundary, NParent, BoundaryElementNode, FlowDOFs
   INTEGER :: material_id, body_id, other_body_id, bf_id
   REAL (KIND=dp) :: C, m, Bdrag, scaleFactor 
   REAL (KIND=dp) :: rho, gravity, homoTemp, T_m
   REAL (KIND=dp) :: C_b, gamma, exponent_p, exponent_q
   REAL (KIND=dp) :: Snt, Snn, ut, un, ut0, normalStress
   REAL (KIND=dp), ALLOCATABLE :: Sig(:,:), normal(:), velo(:), Sn(:) 
   CHARACTER(LEN=MAX_NAME_LEN) :: slidingType 
   LOGICAL :: GotIt, FirstTime = .TRUE., scaling = .FALSE., constantTemp = .FALSE.

   SAVE :: Sig, normal, velo, DIM, Ind, Sn, FirstTime, slidingType, scaling, scaleFactor, constantTemp

!   body_id = -1  

! Get BC
!-------
    BC => GetBC(Model % CurrentElement)

    IF (FirstTime) THEN
        FirstTime = .FALSE.  
        DIM = CoordinateSystemDimension()
        IF ((DIM == 2).OR.(DIM == 3))  THEN
            ALLOCATE(Sig(DIM,DIM),normal(DIM), velo(DIM), Sn(DIM))
        ELSE
            CALL FATAL('USF_sliding(Sliding_Weertman)', 'Bad dimension of the problem')
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
        CALL INFO('USF_Sliding(Sliding_Weertman)', & 
                'Slip condition will be calculated on bedrock', level=3)
        slidingType = GetString(BC, 'Sliding law', GotIt)   
        IF(GotIt) THEN     
            WRITE(Message,'(a,a)') 'Sliding law: ', slidingType     
            CALL INFO('USF_Sliding(Sliding_Weertman)', Message, level=3)   
        ELSE     
            slidingType = 'mechanical'     
            CALL INFO('USF_Sliding(Sliding_Weertman)', &
            'Type of sliding law not found in BC. Setting to purely Mechanical', level=3)   
        END IF

        scaling = GetLogical(BC, 'SeaRise basal drag scaling', GotIt)
        IF(GotIt .AND. scaling) THEN

            scaleFactor = GetConstReal( BC, 'SeaRise scaling factor', GotIt )       
            IF (.NOT. GotIt) THEN
                CALL FATAL('USF_sliding(Sliding_Weertman)', 'SeaRise scaling factor not found')
            END IF
 
            WRITE(Message,'(a,F8.2)') 'SeaRise: Scaling of the basal drag by: ', scaleFactor   
            CALL INFO('USF_Sliding(Sliding_Weertman)', Message, level=3)   

        END IF
  
        constantTemp = GetLogical(BC, 'Constant Temp in Sliding', GotIt)
        IF(GotIt .AND. constantTemp) THEN     
            WRITE(Message,'(a,a)') 'Sliding computed with initial temperature for all times steps.'     
            CALL INFO('USF_Sliding(Sliding_Weertman)', Message, level=3)   
        END IF
  
    END IF

!----------------! get some information upon active boundary element and its parent  !----------------------------------------

 BoundaryElement => Model % CurrentElement

!-----------------------Do nothing if this is a halo-element of a parallel mesh?  !------------------------------------------

 IF (ParEnv % myPe .NE. BoundaryElement % partIndex) RETURN


 IF ( .NOT. ASSOCIATED(BoundaryElement) ) THEN     
 CALL FATAL('USF_Sliding(Sliding_Weertman)','No boundary element found')  
 END IF


 NBoundary = BoundaryElement % Type % NumberOfNodes  
 DO BoundaryElementNode=1,NBoundary     
   IF (nodenumber .EQ. BoundaryElement % NodeIndexes(BoundaryElementNode)) THEN        
      GotIt = .TRUE.        
      EXIT     
   END IF  
 END DO

 IF (.NOT.GotIt) THEN
   CALL WARN('USF_Sliding(Sliding_Weertman)','Node not found in Current Element')     
   Bdrag = 0.0D00
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
     CALL FATAL('USF_Sliding(Sliding_Weertman)',Message)
  END IF

 body_id = ParentElement % BodyId
  material_id = ListGetInteger(Model % Bodies(body_id) % Values, 'Material', GotIt)
  ParentMaterial => Model % Materials(material_id) % Values
  IF ((.NOT. ASSOCIATED(ParentMaterial)) .OR. (.NOT. GotIt)) THEN
     WRITE(Message,'(A,I10,A,I10)')&
          'No material values found for body no ', body_id,&
          ' under material id ', material_id
     CALL FATAL('USF_Sliding(Sliding_Weertman)',Message)
  END IF

!Get BofyForce ID
!----------------
 bf_id = ListGetInteger( Model % Bodies(body_id) % Values, &
              'Body Force', gotIt, 1, Model % NumberOfBodyForces )


! Read the needed variables for Sicopilis sliding law type in sif file
!---------------------------------------------------------------------

!----------------------------------------------------------------------
 SELECT CASE(slidingType)
!----------------------------------------------------------------------
    CASE('thermo-mechanical')


    IF (constantTemp) THEN

        TempHomologous => VariableGet( Model % Variables, 'sicotemp' )
        IF ( ASSOCIATED( TempHomologous ) ) THEN
            TempHomologousPerm    =>  TempHomologous % Perm
            TempHomologousValues  =>  TempHomologous % Values
        ELSE
            CALL FATAL('USF_Sliding(Sliding_Weertman)','No variable sicotemp found')  
        END IF
  
    ELSE  
        !TempHomologous => VariableGet( Model % Variables, 'Temp Homologous' )
        TempHomologous => VariableGet( Model % Variables, 'Temp Homologous' )
        IF ( ASSOCIATED( TempHomologous ) ) THEN
            TempHomologousPerm    =>  TempHomologous % Perm
            TempHomologousValues  =>  TempHomologous % Values
        ELSE
            CALL FATAL('USF_Sliding(Sliding_Weertman)','No variable Temp Homologous found')  
        END IF

    END IF
  
      DepthVariable => VariableGet( Model % Variables, 'Depth' )
      IF ( ASSOCIATED( DepthVariable ) ) THEN
         DepthPerm    =>  DepthVariable % Perm
         DepthValues  =>  DepthVariable % Values
      ELSE
         CALL FATAL('USF_Sliding(Sliding_Weertman)','No variable Depth found')  
      END IF

      rho = GetConstReal(ParentMaterial, 'Density', GotIt)
      IF(.NOT. GotIt) THEN 
         CALL FATAL('USF_Sliding(Sliding_Weertman)', 'Density not found in Material')
      END IF

      gravity = GetConstReal(Model % BodyForces(bf_id) % Values, 'Flow Bodyforce 3', GotIt)
      IF(.NOT. GotIt) THEN 
         CALL FATAL('USF_Sliding(Sliding_Weertman)', 'Flow Body Force 3 not found')
      END IF

      C_b = GetConstReal(BC, 'Thermo-Mechanical Sliding Coefficient', GotIt)
      IF(.NOT. GotIt) THEN 
         C_b = 1.0E05
         WRITE(Message,'(a,I2,a)') 'Thermo-Mechanical Sliding Coefficient not found for node ', &
          nodenumber, '. Setting C_b=1.OE05'
      CALL INFO('USF_Sliding(Sliding_Weertman)', Message, level=3)
      END IF

      gamma = GetConstReal(BC, 'Thermo-Mechanical Sliding Gamma', GotIt)
      IF(.NOT. GotIt) THEN 
         gamma = 1.0
      WRITE(Message,'(a,I2,a)') 'Value for Thermo-Mechanical Sliding Gamma not found for node ', & 
      nodenumber, '. Setting gamma=1.O'
      CALL INFO('USF_Sliding(Sliding_Weertman)', Message, level=3)
      END IF
 
      exponent_p = GetConstReal(BC, 'Thermo-Mechanical Sliding Exponent p', GotIt)
      IF(.NOT. GotIt) THEN 
         exponent_p  = 3.0
      WRITE(Message,'(a,I2,a)') 'Value for Thermo-Mechanical Sliding Exponent p not found for node ', & 
      nodenumber, '. Setting p=3.O'
      CALL INFO('USF_Sliding(Sliding_Weertman)', Message, level=3)
      END IF

      exponent_q = GetConstReal(BC, 'Thermo-Mechanical Sliding Exponent q', GotIt)
      IF(.NOT. GotIt) THEN 
         exponent_q  = 2.0
      WRITE(Message,'(a,I2,a)') 'Value for Thermo-Mechanical Sliding Exponent q not found for node ', & 
      nodenumber, '. Setting q=2.O'
      CALL INFO('USF_Sliding(Sliding_Weertman)', Message, level=3)
      END IF

   CASE('mechanical')

   ! Read the coefficients C and m for default law in the sif file 
   !--------------------------------------------------------------
   C = GetConstReal( BC, 'Mechanical Sliding Coefficient', GotIt )
   IF (.NOT.GotIt) THEN
      CALL FATAL('USF_sliding(Sliding_Weertman)', 'Need a Mechanical Sliding Coefficient for the choosen sliding law (weertman)') 
   END IF

   m = GetConstReal( BC, 'Mechanical Sliding Exponent', GotIt )
   IF (.NOT.GotIt) THEN
      CALL FATAL('USF_sliding(Sliding_Weertman)', 'Need a Mechanical Sliding Exponent for the choosen sliding law (weertman)')
   END IF 

   CASE('t_m-based')
  
   Temp => VariableGet( Model % Variables, 'Temp' )
   IF ( ASSOCIATED( Temp ) ) THEN
       TempPerm    =>  Temp % Perm
       TempValues  =>  Temp % Values
   ELSE
       CALL FATAL('USF_Sliding(Sliding_Weertman)','No variable Temp found')
   END IF

   C = GetConstReal( BC, 'Mechanical Sliding Coefficient', GotIt )
   IF (.NOT.GotIt) THEN
      CALL FATAL('USF_sliding(Sliding_Weertman)', 'Need a Mechanical Sliding Coefficient for the choosen sliding law (weertman)')
   END IF

   m = GetConstReal( BC, 'Mechanical Sliding Exponent', GotIt )
   IF (.NOT.GotIt) THEN
      CALL FATAL('USF_sliding(Sliding_Weertman)', 'Need a Mechanical Sliding Exponent for the choosen sliding law (weertman)')
   END IF

   T_m = GetConstReal(ParentMaterial, TRIM('Temp') // ' Upper Limit', GotIt)
   IF(.NOT. GotIt) THEN
      CALL FATAL('USF_Sliding(Sliding_Weertman)', 'Temp Upper Limit not found in Material')
   END IF

!----------------------------------------------------------------------
 END SELECT
!----------------------------------------------------------------------

   ut0 = GetConstReal( BC, 'Sliding Linear Velocity', GotIt )
   IF (.NOT.GotIt) THEN
      CALL FATAL('USF_sliding(Sliding_Weertman)', 'Need a Sliding Linear Velocity for the choosen sliding law (weertman)')
   END IF

! Get the variables to compute tau_b
!-----------------------------------
   StressVariable => VariableGet( Model % Variables, 'Stress' )
   IF ( ASSOCIATED( StressVariable ) ) THEN
      StressPerm    => StressVariable % Perm
      StressValues  => StressVariable % Values
   ELSE
      CALL FATAL('USF_Sliding(Sliding_Weertman)','No variable Stress found')   
   END IF

! Get the variables to compute ut
!--------------------------------
   FlowVariable => VariableGet( Model % Variables, 'Flow Solution' )
   IF ( ASSOCIATED( FlowVariable ) ) THEN
      FlowPerm    => FlowVariable % Perm
      FlowValues  => FlowVariable % Values
      FlowDOFs = FlowVariable % DOFs
   ELSE
      CALL FATAL('USF_sliding(Sliding_Weertman)', 'Need NS Solver, Flow Solution not found')
   END IF

! Get the variable to compute the normal
!---------------------------------------
   NormalVar =>  VariableGet(Model % Variables,'Normal Vector')
   IF ( ASSOCIATED( NormalVar ) ) THEN
      NormalPerm => NormalVar % Perm
      NormalValues => NormalVar % Values
   ELSE
      CALL FATAL('USF_sliding(Sliding_Weertman)', 'Need ComputeNormal Solver, Normal Vector not found')
   END IF

   DO i=1, DIM
     normal(i) = -NormalValues(DIM*(NormalPerm(Nodenumber)-1) + i)
     velo(i) = FlowValues( (DIM+1)*(FlowPerm(Nodenumber)-1) + i )
   END DO
   
   un = SUM(velo(1:DIM)*normal(1:DIM)) 
   ut = SQRT( SUM( (velo(1:DIM)-un*normal(1:DIM))**2.0 ) )

   IF ( ASSOCIATED( StressVariable ) ) THEN

       DO i=1, DIM
         DO j= 1, DIM
            Sig(i,j) =  &
              StressValues( 2*DIM *(StressPerm(Nodenumber)-1) + Ind(i,j) )
         END DO
       END DO

! Stress vector Sn
!-----------------       
       DO i=1, DIM
         Sn(i) = SUM(Sig(i,1:DIM)*normal(1:DIM)) 
       END DO  

       Snn = SUM( Sn(1:DIM) * normal(1:DIM) ) 
       Snt = SQRT( SUM((Sn(1:DIM) - Snn*normal(1:DIM))**2.0 ))
   END IF
 


!Basal Normal Stress (use pressure solution from Stokes equation)
!---------------------------------------------------------
IF (slidingType == 'thermo-mechanical') THEN
   !We need a positive value for the gravity. It was orginally negative as it was acquired from the body force value in the model setting.
   gravity = gravity*(-1.0)
   
   normalStress  =  FlowValues(FlowDOFs*(FlowPerm(nodenumber)-1)+(DIM+1))
   !we should not have negative pressure at the bedrock but who knows. 
   IF (normalStress < 0.0) THEN
      normalStress = normalStress*(-1.0)
   END IF
!-----------------------------------------------------------

   homoTemp = TempHomologousValues(TempHomologousPerm(nodenumber))
   IF (constantTemp) THEN
        homoTemp = homoTemp+273.15_dp
        homoTemp = homoTemp-(273.15_dp-(9.8E-08_dp*(910.0_dp*9.81_dp*DepthValues(DepthPerm(nodenumber)))))
   END IF
   
END IF


! Compute the basal drag (slip coefficient)
!------------------------------------------

!----------------------------------------------------------------------
 SELECT CASE(slidingType)
!----------------------------------------------------------------------
  CASE('thermo-mechanical')
   IF ( ut > 1.0e-6_dp ) THEN
      ! - case 3 and 4
      IF (ut > ut0) THEN
       Bdrag = ut**(1.0/exponent_p-1.0)* &
       ((rho*gravity)/(C_b*EXP(homoTemp/gamma)))** &
       (1.0/exponent_p)*normalStress**(exponent_q/exponent_p)
       !WRITE(*,*) '3)', Bdrag
      ELSE
       Bdrag = ut0**(1.0/exponent_p-1.0)* &
       ((rho*gravity)/(C_b*EXP(homoTemp/gamma)))** &
       (1.0/exponent_p)*normalStress**(exponent_q/exponent_p)
       !WRITE(*,*) '4)', Bdrag
      END IF
   ELSE
      IF (Snt > 1.0e-6_dp) THEN
      ! - case 2           
      Bdrag = Snt**(1.0-exponent_p)* &
      ((rho*gravity)/(C_b*EXP(homoTemp/gamma)))*normalStress**exponent_q
      !WRITE(*,*) '2)', Bdrag
      ELSE
      ! - case 1
      Bdrag = 1.0e20 
      !WRITE(*,*) '1)', Bdrag
      END IF
   END IF

  CASE('mechanical')
   IF ( ut > 1.0e-6_dp ) THEN
      ! - case 3 and 4           
      IF (ut > ut0) THEN
       Bdrag = C * ut**(m-1.0)
      ELSE
       Bdrag = C * ut0**(m-1.0)
      END IF
   ELSE 
      IF (Snt > 1.0e-6_dp) THEN
      ! - case 2           
       Bdrag = C**(1.0/m) * Snt**(1.0 - 1.0/m)
      ELSE
      ! - case 1           
      Bdrag = 1.0e20 
      END IF        
   END IF

  CASE('t_m-based')
   IF ( TempValues(TempPerm(nodenumber)) >= T_m ) THEN
     IF ( ut > 1.0e-6_dp ) THEN
      IF (ut > ut0) THEN
       Bdrag = C * ut**(m-1.0)
      ELSE
       Bdrag = C * ut0**(m-1.0)
      END IF
     ELSE 
       IF (Snt > 1.0e-6_dp) THEN
        Bdrag = C**(1.0/m) * Snt**(1.0 - 1.0/m)
       ELSE
        Bdrag = 1.0e20 
       END IF
     END IF
   ELSE 
      Bdrag = 1.0e20
   END IF
   
  CASE DEFAULT
     CALL FATAL('USF_sliding(Sliding_Weertman)', 'Sliding law not supported')
   
!----------------------------------------------------------------------
 END SELECT
!----------------------------------------------------------------------

 IF (scaling) THEN
    Bdrag = Bdrag/scaleFactor
 END IF

 Bdrag = MIN(Bdrag,1.0e20)


 END FUNCTION Sliding_Weertman 
