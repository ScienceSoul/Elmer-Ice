! *  Elmer, A Finite Element Software for Multiphysical Problems
! *
! *  Copyright 1st April 1995 - , CSC - Scientific Computing Ltd., Finland
! *
! *  This program is free software; you can redistribute it and/or
! *  modify it under the terms of the GNU General Public License
! *  as published by the Free Software Foundation; either version 2
! *  of the License, or (at your option) any later version.
! *
! *  This program is distributed in the hope that it will be useful,
! *  but WITHOUT ANY WARRANTY; without even the implied warranty of
! *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! *  GNU General Public License for more details.
! *
! *  You should have received a copy of the GNU General Public License
! *  along with this program (in file fem/GPL-2); if not, write to the
! *  Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor,
! *  Boston, MA 02110-1301, USA.
! *
! *****************************************************************************/
!
!/******************************************************************************
! *
! *  Module containing solvers and routines for standard ice dynamic problems
! *
! ******************************************************************************
! *
! *  Authors: Thomas Zwinger, Olivier Gagliardini, Juha.Ruokolainen
! *  Email:   Thomas.Zwinger@csc.fi, Juha.Ruokolainen@csc.fi 
! *  Web:     http://www.csc.fi/elmer
! *  Address: CSC - Scientific Computing Ltd.
! *           Keilaranta 14
! *           02101 Espoo, Finland 
! *
! *  Original Date: 14 May 2007
! *
! *****************************************************************************/

!*********************************************************************************************************************************
!*
!* heat conductivity of ice as a function of temperature (K):  k = c_1 * exp(c_2 * T[K]); c_2 < 0 
!*
!*********************************************************************************************************************************
FUNCTION getHeatConductivity( Model, N, temperature ) RESULT(conductivity)
  USE types
  USE CoordinateSystems
  USE SolverUtils
  USE ElementDescription
!-----------------------------------------------------------
  IMPLICIT NONE
!------------ external variables ---------------------------
  TYPE(Model_t) :: Model
  INTEGER :: N
  REAL(KIND=dp) :: temperature, conductivity
!------------ internal variables----------------------------
  TYPE(ValueList_t), POINTER :: Material
  INTEGER :: nMax,i,j,body_id,material_id,elementNodes,nodeInElement,istat
  REAL (KIND=dp), ALLOCATABLE :: conductivityExponentFactor(:), conductivityFactor(:)
  LOGICAL :: FirstTime = .TRUE., GotIt
  CHARACTER(LEN=MAX_NAME_LEN) :: SolverName
!------------ remember this -------------------------------
  Save FirstTime, conductivityExponentFactor, conductivityFactor, SolverName
  !-------------------------------------------
  ! Allocations 
  !------------------------------------------- 
  IF (FirstTime) THEN
     WRITE(SolverName, '(A)') 'IceFlowProperties (getHeatConductivity))'
     nMax = Model % MaxElementNodes
     ALLOCATE(conductivityExponentFactor(nMax),&
          conductivityFactor(nMax),&
          STAT=istat)
     IF ( istat /= 0 ) THEN
        CALL FATAL(SolverName,'Memory allocation error, Aborting.')
     END IF
     FirstTime = .FALSE.
     CALL INFO(SolverName,'Memory allocation done', level=3)
  END IF
  !-------------------------------------------
  ! get element properties
  !-------------------------------------------   
  IF ( .NOT. ASSOCIATED(Model % CurrentElement) ) THEN
     CALL FATAL(SolverName, 'Model % CurrentElement not associated')
  END IF
  body_id = Model % CurrentElement % BodyId
  material_id = ListGetInteger(Model % Bodies(body_id) % Values, 'Material', GotIt)
  elementNodes = Model % CurrentElement % Type % NumberOfNodes
  IF (.NOT. GotIt) THEN
     WRITE(Message,'(a,I2,a,I2,a)') 'No material id for current element of node ',n,', body ',body_id,' found'
     CALL FATAL(SolverName, Message)
  END IF
  DO nodeInElement=1,elementNodes
     IF ( N == Model % CurrentElement % NodeIndexes(nodeInElement) ) EXIT
  END DO
  Material => Model % Materials(material_id) % Values
  !-------------------------------------------
  ! get material properties
  !-------------------------------------------
  conductivityExponentFactor(1:elementNodes) = ListGetReal( Material,'Conductivity Exponent Factor', elementNodes, &
       Model % CurrentElement % NodeIndexes, GotIt )
  IF (.NOT. GotIt) THEN
     WRITE(Message,'(a,I2,a,I2,a)') 'No Conductivity Exponent Factor found in Material ', &
          material_id,' for node ', n, '.setting E=1'
     CALL FATAL(SolverName, Message)
  END IF
  conductivityFactor(1:elementNodes) = ListGetReal( Material,'Conductivity Factor', elementNodes, &
       Model % CurrentElement % NodeIndexes, GotIt )
  IF (.NOT. GotIt) THEN
     WRITE(Message,'(a,I2,a,I2,a)') 'No Conductivity Factor found in Material ', material_id,' for node ', n, '.setting E=1'
     CALL FATAL(SolverName, Message)
  END IF
  !-------------------------------------------
  ! compute heat conductivity
  !-------------------------------------------
  conductivity = conductivityFactor(nodeInElement)*EXP(conductivityExponentFactor(nodeInElement)*temperature)
END FUNCTION getHeatConductivity

!*********************************************************************************************************************************
!*
!* heat capacity of ice as a function of temperature [K] (internaly transformed into Celsius):  k = c_1  + c_2 * T[C];
!*
!*********************************************************************************************************************************
FUNCTION getHeatCapacity( Model, N, temperature ) RESULT(capacity)
  USE types
  USE CoordinateSystems
  USE SolverUtils
  USE ElementDescription
  USE DefUtils
!-----------------------------------------------------------
  IMPLICIT NONE
!------------ external variables ---------------------------
  TYPE(Model_t) :: Model
  INTEGER :: N
  REAL(KIND=dp) :: temperature, capacity
!------------ internal variables----------------------------
  REAL(KIND=dp) :: celsius
  INTEGER :: body_id, istat, i
  INTEGER :: nodeInElement(1), elementNodes
  TYPE(ValueList_t), POINTER :: Material
  REAL (KIND=dp) :: C(2)
  LOGICAL :: FirstTime = .TRUE., GotIt
  CHARACTER(LEN=MAX_NAME_LEN) :: SolverName
!------------ remember -------------------------------------
  SAVE SolverName, FirstTime

  IF (FirstTime) THEN
     WRITE(SolverName, '(A)') 'IceFlowProperties (getHeatCapacity)'
     FirstTime = .FALSE.
     CALL INFO(SolverName,'Initialization done', level=1)
  END IF
  !-------------------------------------------
  ! get element properties
  !-------------------------------------------   
  IF ( .NOT. ASSOCIATED(Model % CurrentElement) ) THEN
     CALL FATAL(SolverName, 'Model % CurrentElement not associated')
  END IF
  body_id = Model % CurrentElement % BodyId
  Material => GetMaterial()
  elementNodes = Model % CurrentElement % Type % NumberOfNodes


  DO i=1,elementNodes
     IF ( N == Model % CurrentElement % NodeIndexes(i) ) EXIT
  END DO
  nodeInElement(1) = i

  IF (.NOT.ASSOCIATED(Material)) CALL FATAL(SolverName,'No Material found')

  C(1:1) = ListGetReal( Material,'Reference Capacity', 1, &
       Model % CurrentElement % NodeIndexes, GotIt )
  IF (.NOT.GotIt) THEN
     WRITE(Message,'(A,I10)') 'Keyword >Reference Capacity< not found for node no. ', N
     CALL INFO(SolverName, Message, Level = 4)
     CALL INFO(SolverName, 'Using default value of 2.1275x10^03', Level = 4)
     C(1) = 2.1275D03
  END IF
  C(2:2) = ListGetReal( Material,'Capacity Temperature Gradient', 1, &
       NodeInElement, GotIt )
  IF (.NOT.GotIt) THEN
     WRITE(Message,'(A,I10)') 'Keyword >Capacity Temperature Gradient< not found for node no. ', N
     CALL INFO(SolverName, Message, Level = 4)
     CALL INFO(SolverName, 'Using default value of 7.253', Level = 4)
     C(2) = 7.253D00
  END IF
  !-------------------------------------------
  ! compute celsius temperature and limit it 
  ! to 0 deg
  !------------------------------------------- 
  celsius = MIN(temperature - 2.7316D02,0.0d00)
  !-------------------------------------------
  ! compute heat capacity
  !-------------------------------------------  
  capacity = C(1) + C(2) * celsius
END FUNCTION getHeatCapacity

!*********************************************************************************************************************************
!*
!* pressure melting point: T_m = 273.16 - beta * Pressure
!*
!*********************************************************************************************************************************
FUNCTION getPressureMeltingPoint( Model, n, Pressure ) RESULT(PressureMeltingPoint)
   USE types
   USE CoordinateSystems
   USE SolverUtils
   USE ElementDescription
   USE DefUtils
!-----------------------------------------------------------
   IMPLICIT NONE
!------------ external variables ---------------------------
   TYPE(Model_t) :: Model
   INTEGER :: n
   REAL(KIND=dp) :: PressureMeltingPoint, Pressure
!------------ internal variables----------------------------
   TYPE(Element_t), POINTER :: Element
   TYPE(ValueList_t),POINTER :: Material
   INTEGER :: body_id, material_id, nodeInElement, istat, elementNodes
   REAL(KIND=dp), ALLOCATABLE :: ClausiusClapeyron(:)
   LOGICAL :: GotIt, FirstTime = .TRUE., UseCelsius = .FALSE.
   CHARACTER(LEN=MAX_NAME_LEN) :: SolverName
!------------ remember this -------------------------------
   Save FirstTime, ClausiusClapeyron, SolverName
  !-------------------------------------------
  ! Allocations 
  !------------------------------------------- 
  IF (FirstTime) THEN
     WRITE(SolverName, '(A)') 'IceFlowProperties (getPressureMeltingPoint)'
     ALLOCATE(ClausiusClapeyron(Model % MaxElementNodes),&
          STAT=istat)
     IF ( istat /= 0 ) THEN
        CALL FATAL(SolverName,'Memory allocation error, Aborting.')
     END IF
     FirstTime = .FALSE.
     CALL INFO(SolverName,'Memory allocation done', level=3)
  END IF
  !-------------------------------------------
  ! get element properties
  !-------------------------------------------   
  IF ( .NOT. ASSOCIATED(Model % CurrentElement) ) THEN
     CALL FATAL(SolverName, 'Model % CurrentElement not associated')
  ELSE
     Element => Model % CurrentElement
  END IF
  body_id = Element % BodyId
  material_id = ListGetInteger(Model % Bodies(body_id) % Values, 'Material', GotIt)
  elementNodes = Element % Type % NumberOfNodes
  IF (.NOT. GotIt) THEN
     WRITE(Message,'(a,I2,a,I2,a)') 'No material id for current element of node ',elementNodes,', body ',body_id,' found'
     CALL FATAL(SolverName, Message)
  END IF
  DO nodeInElement=1,elementNodes
     IF ( N == Model % CurrentElement % NodeIndexes(nodeInElement) ) EXIT
  END DO
  Material => Model % Materials(material_id) % Values
  !-------------------------
  ! get material parameters
  !-------------------------
  ClausiusClapeyron(1:elementNodes) = &
       ListGetReal( Material, 'Clausius Clapeyron', elementNodes, Element % NodeIndexes, GotIt)
  IF (.NOT. GotIt) THEN
     CALL FATAL(SolverName,'No value for Clausius Clapeyron parameter found')
  END IF
  !-------------------------------
  ! compute pressure melting point
  !-------------------------------
  PressureMeltingPoint = 2.7316D02 - ClausiusClapeyron(nodeInElement)*MAX(Pressure,0.0d00)
END FUNCTION getPressureMeltingPoint



!**************************************************************************
!*
!*  computes effective heat flux
!*
!**************************************************************************
RECURSIVE SUBROUTINE getNetBoundaryHeatflux( Model,Solver,Timestep,TransientSimulation)
  USE DefUtils
  USE Materialmodels
  !-----------------------------------------------------------
  IMPLICIT NONE
  !------------ external variables ---------------------------
  TYPE(Model_t)  :: Model
  TYPE(Solver_t), TARGET :: Solver
  LOGICAL :: TransientSimulation
  REAL(KIND=dp) :: Timestep
  !------------------------------------------------------------------------------
  !    Local variables
  !------------------------------------------------------------------------------
  TYPE(Solver_t), POINTER :: PointerToSolver
  TYPE(Nodes_t) :: Nodes
  TYPE(Element_t),POINTER :: Element
  TYPE(Variable_t), POINTER :: HFSol, TempSol
  TYPE(ValueList_t), POINTER :: Equation,Material,SolverParams,BodyForce,BC,Constants
  CHARACTER(LEN=MAX_NAME_LEN) :: TempName, SolverName
  INTEGER :: i, j, k, t, N, BN,  M, DIM, istat, material_id
  INTEGER, POINTER :: TempPerm(:), HFPerm(:),BoundaryReorder(:)
  INTEGER, ALLOCATABLE :: contributedElements(:)
  REAL(KIND=dp) ::  U, V, W, SqrtElementMetric
  REAL(KIND=dp), ALLOCATABLE ::  Basis(:),dBasisdx(:,:), ddBasisddx(:,:,:),&
       LatentHeat(:), HeatConductivity(:), Density(:), heatFlux(:,:)
  REAL(KIND=dp), POINTER :: HF(:), Temp(:), BoundaryNormals(:,:), BoundaryTangent1(:,:), BoundaryTangent2(:,:)
  LOGICAL :: GotIt, AllocationsDone=.FALSE., stat

  SAVE AllocationsDone, SolverName, DIM,  heatFlux, &
       Basis, dBasisdx, ddBasisddx, Nodes, &
       LatentHeat, HeatConductivity, Density, contributedElements, &
       BoundaryNormals, BoundaryTangent1, BoundaryTangent2, BoundaryReorder, BN


  WRITE(SolverName, '(A)') 'IceFlowProperties (getNetBoundaryHeatflux))'
  !-----------------------------------------------------------------------
  ! get solver variable
  !-----------------------------------------------------------------------
  HFSol => Solver % Variable
  IF (.NOT.ASSOCIATED(HFSol)) THEN
     CALL FATAL(SolverName,'No variable associated')
  END IF
  HFPerm  => HFSol % Perm
  HF => HFSol % Values  
  HF = 0.0D00
  !-----------------------------------------------------------------------
  ! Allocations
  !-----------------------------------------------------------------------
  IF ( .NOT.AllocationsDone .OR. Solver % Mesh % Changed ) THEN
     DIM = CoordinateSystemDimension()
     N = Solver % Mesh % MaxElementNodes
     M = Model % Mesh % NumberOfNodes        
     IF ( AllocationsDone ) &
          DEALLOCATE( &
          Basis, &
          dBasisdx, &
          ddBasisddx, &
          LatentHeat, &
          HeatConductivity, &
          Density, &
          heatFlux,&
          Nodes % x, &
          Nodes % y, &
          Nodes % z, &
          contributedElements, &
          BoundaryReorder, &
          BoundaryNormals, &
          BoundaryTangent1, &
          BoundaryTangent2)
     CALL CheckNormalTangentialBoundary( Model, &
          'Basal Melting', BN, &
          BoundaryReorder, BoundaryNormals, BoundaryTangent1, &
          BoundaryTangent2, DIM )
     WRITE(Message,'(A,i6)') &
          'Number of boundary nodes on boundaries associated with melting:',&
          BN
     CALL INFO(SolverName,Message,Level=3)
     IF (BN > 0) THEN
        CALL AverageBoundaryNormals(Model, &
             'Basal Melting', BN, &
             BoundaryReorder, BoundaryNormals, BoundaryTangent1, &
             BoundaryTangent2, DIM )
     END IF
     ALLOCATE( &
          Basis(N), &
          dBasisdx(N,3), &
          ddBasisddx(N,3,3), &
          LatentHeat(N), &
          HeatConductivity(N), &
          Density(N), &
          heatFlux(M,3),& 
          Nodes % x(N), &
          Nodes % y(N), &
          Nodes % z(N), &
          contributedElements(M), &
          STAT = istat)
     AllocationsDone = .TRUE.
  END IF
  !-------------------------------------------------------------------------
  ! Calculate element-wise heat flux contributions for each point of element
  !-------------------------------------------------------------------------
  contributedElements = 0
  heatFlux = 0.0D00  

  IF (BN > 0) THEN
     DO t=1,Solver % NumberOfActiveElements
        Element => GetActiveElement(t)
        Model % CurrentElement => Element
        N = Element % Type % NumberOfNodes
        CALL GetElementNodes( Nodes )
        !-------------
        ! get material
        !-------------
        Material => GetMaterial()
        IF (.NOT.ASSOCIATED(Material)) THEN
           WRITE (Message,'(A,I3)') 'No Material found for boundary element no. ', t
           CALL FATAL(SolverName,Message)
        ELSE
           material_id = GetMaterialId( Element, GotIt)
           IF(.NOT.GotIt) THEN
              WRITE (Message,'(A,I3)') 'No Material ID found for boundary element no. ', t
              CALL FATAL(SolverName,Message)
           END IF
        END IF
        !-----------------------------------------------------------------------
        ! get temperature variable
        !-----------------------------------------------------------------------
        TempName =  GetString(Material ,'Temperature Name', GotIt)
        IF (.NOT.GotIt) THEN
           CALL FATAL(SolverName,'No Temperature Name found in Material section')
        ELSE
           WRITE(Message,'(a,a)') 'Variable Name for temperature: ', TempName
           CALL INFO(SolverName,Message,Level=21)
        END IF
        TempSol => VariableGet( Solver % Mesh % Variables, TRIM(TempName) )
        IF ( ASSOCIATED( TempSol ) ) THEN
           TempPerm => TempSol % Perm
           Temp => TempSol % Values
        ELSE
           WRITE(Message, '(A,A)') 'Could not find temperature field variable ',  TRIM(TempName)
           CALL FATAL(SolverName,Message)
        END IF
        !-----------------------
        ! get material parameter
        !-----------------------
        HeatConductivity(1:N) =  ListGetReal( Material,  TRIM(TempName) // &
             ' Heat Conductivity', n, Element % NodeIndexes, GotIt )
        IF (.NOT.GotIt) THEN
           HeatConductivity = 0.0D00
           WRITE(Message,'(a,a,a,i5,a,i5,a)') 'Keyword >', TRIM(TempName) // &
                ' Heat Conductivity', '< not found for element ', t, ' material ', material_id
           CALL FATAL(SolverName,Message)
        END IF
        !--------------------------
        ! loop all nodes in element
        !--------------------------
        DO i=1,N
           U = Element % Type % NodeU(i)
           V = Element % Type % NodeV(i)
           W =Element % Type % NodeW(i)
           stat = ElementInfo( Element,Nodes,U,V,W,SqrtElementMetric, &
                Basis,dBasisdx,ddBasisddx,.FALSE.,.FALSE. )
           !--------------------
           ! loop all dimensions
           !--------------------
           DO j=1,DIM
              DO k=1,N
                 heatFlux(Element % NodeIndexes(i),j) = & 
                      heatFlux(Element % NodeIndexes(i),j) - HeatConductivity(k)*dBasisdx(k,j)*Temp(TempPerm(k))
              END DO
           END DO
           contributedElements(Element % NodeIndexes(i)) = contributedElements(Element % NodeIndexes(i)) + 1
        END DO
     END DO

     !----------------------------------------------------------
     ! averaged internal flux for each point of element
     !----------------------------------------------------------
     DO i=1,Solver % Mesh % NumberOfNodes
        IF (contributedElements(i) > 0) THEN
           heatFlux(i,1:DIM) = heatFlux(i,1:DIM)/contributedElements(i)           
           j = BoundaryReorder(i)
           IF (j>0) THEN
              HF(HFPerm(i)) = SUM(heatFlux(i,1:DIM)*BoundaryNormals(j,1:DIM))
           ELSE
              HF(HFPerm(i)) = 0.0D00
           END IF
        ELSE
           WRITE(Message,'(A,i6,A)') 'Node no. ',i,' appears not to be part of any element.'
           CALL WARN(SolverName,Message)
           HF(HFPerm(i)) = 0.0
        END IF
     END DO
  END IF

END SUBROUTINE getNetBoundaryHeatflux

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
  TYPE(Variable_t), POINTER :: VarTemp, VarResidual, FlowVariable, StressVariable, NormalVar
  TYPE(Variable_t), POINTER :: basalMeltingSol
  REAL(KIND=dp), POINTER :: basalMeltingValues(:), FlowValues(:), StressValues(:)
  REAL(KIND=dp), POINTER :: NormalValues(:)
  INTEGER, POINTER :: basalMeltingPerm(:), FlowPerm(:), StressPerm(:), NormalPerm(:)
  INTEGER :: basalMeltingDOFs
  INTEGER :: i, j, N, DIM, NBoundary, NParent, BoundaryElementNode, Ind(3,3), &
       ParentElementNode, body_id, other_body_id, material_id, equation_id, istat
  REAL(KIND=dp), ALLOCATABLE :: LatentHeat(:), Density(:),ExternalHF(:),PressureMeltingPoint(:)
  REAL(KIND=dp), ALLOCATABLE :: Sig(:,:), normal(:), velo(:), Sn(:), ut(:)  
  REAL (KIND=dp) :: HomTemp, BasalMeltingOffset, vSn, un, limitValue
  CHARACTER(LEN=MAX_NAME_LEN) :: TempName, FlowSolName, SolverName
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
   FlowSolName =  GetString( Model % Equations(equation_id) % Values,'Flow Solution Name', GotIt)

  ! Get the variable velocity
  !---------------------------
  FlowVariable => VariableGet( Model % Variables, FlowSolName )
  IF ( ASSOCIATED( FlowVariable ) ) THEN
      FlowPerm    => FlowVariable % Perm
      FlowValues  => FlowVariable % Values
  ELSE
      CALL FATAL('IceFlowPoperties(getBasalMeltingVelocity)', 'Need NS Solver, Flow Solution not found')
  END IF

  ! Get the stress variable
  !------------------------
  StressVariable => VariableGet( Model % Variables, 'Stress' )
  IF ( ASSOCIATED( StressVariable ) ) THEN
      StressPerm    => StressVariable % Perm
      StressValues  => StressVariable % Values
  ELSE
      CALL FATAL('IceFlowPoperties(getBasalMeltingVelocity)','No variable Stress found')   
  END IF

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
  
  !WRITE(*,*) velo
  
  DO i=1, DIM
     DO j= 1, DIM
        Sig(i,j) =  &
           StressValues( 2*DIM *(StressPerm(Node)-1) + Ind(i,j) )
     END DO
  END DO

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



!*********************************************************************************************************************************
!*
!*  basal slip coefficient as a function of temperature
!*
!*********************************************************************************************************************************

FUNCTION basalSlip( Model, Node, Temperature ) RESULT(basalSlipCoefficient)
!-----------------------------------------------------------
  USE DefUtils
!-----------------------------------------------------------
  IMPLICIT NONE
!-----------------------------------------------------------
  !external variables
  TYPE(Model_t), TARGET :: Model
  INTEGER :: Node
  REAL(KIND=dp) :: Temperature, basalSlipCoefficient
  !internal variables
  TYPE(Element_t), POINTER :: CurrentElementAtBeginning, BoundaryElement, ParentElement
  TYPE(ValueList_t), POINTER :: ParentMaterial, BC
  TYPE(Variable_t), POINTER :: varTemperature, varPressure
  INTEGER :: N, NBoundary, NParent, BoundaryElementNode, ParentElementNode, &
       i, DIM, other_body_id, body_id, material_id, istat, NSDOFs
  REAL(KIND=dp) :: TempHom, ThermalCoefficient, SlipCoefficientMax
  REAL(KIND=dp), ALLOCATABLE :: PressureMeltingPoint(:), TemperateSlipCoefficient(:)
  CHARACTER(LEN=MAX_NAME_LEN) :: TempName
  LOGICAL ::  GotIt, stat, Jump=.FALSE.

  !---------------
  ! Initialization
  !---------------
  CurrentElementAtBeginning => Model % CurrentElement 
  N = Model % MaxElementNodes 
  ALLOCATE(TemperateSlipCoefficient(N),&  
       PressureMeltingPoint(N),&
       STAT = istat)
  IF (istat /= 0) THEN
     CALL FATAL('iceproperties (basalSlip)','Allocations failed')
  END IF

  TemperateSlipCoefficient = 1.0D30 ! high value - > no slip by default
  PressureMeltingPoint = 273.16D00

  !-----------------------------------------------------------------
  ! get some information upon active boundary element and its parent
  !-----------------------------------------------------------------
  BoundaryElement => Model % CurrentElement
  IF ( .NOT. ASSOCIATED(BoundaryElement) ) THEN
     CALL FATAL('iceproperties (basalMelting)','No boundary element found')
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
     CALL FATAL('iceproperties (basalMelting)',Message)
  END IF  
  Model % CurrentElement => ParentElement
  body_id = ParentElement % BodyId
  material_id = ListGetInteger(Model % Bodies(body_id) % Values, 'Material', GotIt)
  ParentMaterial => Model % Materials(material_id) % Values
  IF ((.NOT. ASSOCIATED(ParentMaterial)) .OR. (.NOT. GotIt)) THEN
     WRITE(Message,'(A,I10,A,I10)')&
          'No material values found for body no ', body_id,&
          ' under material id ', material_id
     CALL FATAL('iceproperties (basalMelting)',Message)
  END IF
  ! number of nodes and node in elements
  NBoundary = BoundaryElement % Type % NumberOfNodes
  NParent = ParentElement % Type % NumberOfNodes
  DO BoundaryElementNode=1,Nboundary
     IF ( Node == BoundaryElement % NodeIndexes(BoundaryElementNode) ) EXIT
  END DO
  DO ParentElementNode=1,NParent
     IF ( Node == ParentElement % NodeIndexes(ParentElementNode) ) EXIT
  END DO
  !-------------------------
  ! Get Pressure Melting Point
  !-------------------------
  TempName =  GetString(Model % Constants ,'Temperature Name', GotIt)
  PressureMeltingPoint(1:NParent) =&
       ListGetReal( ParentMaterial, TRIM(TempName) // ' Upper Limit',&
       NParent, ParentElement % NodeIndexes, GotIt)
  IF (.NOT. GotIt) THEN
     WRITE(Message,'(a,a,a)') 'No entry for ', TRIM(TempName) // ' Upper Limit', ' found'
     CALL FATAL('iceproperties (basalMelting)',Message)
  END IF
  !---------------------------------------------------------------------
  ! get temperate slip coefficient and upper limit for slip coefficient
  !---------------------------------------------------------------------
  Model % CurrentElement => BoundaryElement
  BC => GetBC()

  IF (.NOT.ASSOCIATED(BC)) THEN
     CALL FATAL('iceproperties (basalSlip)','No Boundary Condition associated')
  ELSE
     TemperateSlipCoefficient(1:NBoundary) = GetReal(BC, 'Temperate Slip Coefficient', GotIt)
     IF (.NOT. GotIt) THEN
        CALL WARN('iceproperties (basalSlip)','Keyword >Temperate Slip Coefficient< not found')
        CALL WARN('iceproperties (basalSlip)','Asuming Default 5.63D08 [kg /(s m^2)]')
        TemperateSlipCoefficient(1:NBoundary) = 5.63D08
     END IF
     ThermalCoefficient = GetConstReal(BC, 'Thermal Coefficient', GotIt)
     IF (.NOT. GotIt) THEN
        CALL WARN('iceproperties (basalSlip)','Keyword >Thermal Coefficient< not found')
        CALL WARN('iceproperties (basalSlip)','Asuming Default 1 [1/K]')
        ThermalCoefficient =  1.0D00
     END IF
     SlipCoefficientMax = GetConstReal(BC, 'Max Slip Coefficient', GotIt)
     IF (.NOT. GotIt) THEN
        CALL WARN('iceproperties (basalSlip)','Keyword >Max Slip Coefficient< not found')
        CALL WARN('iceproperties (basalSlip)','Asuming Default 1.0E+16 [kg /(s m^2)]')
        SlipCoefficientMax = 1.0E+16
     END IF
  END IF
  !------------------------------
  ! check homologous temperature
  !------------------------------
  TempHom = MIN(Temperature - PressureMeltingPoint(ParentElementNode),0.0D00)
  !------------------------------
  ! get the result and check it
  !------------------------------
  basalSlipCoefficient = TemperateSlipCoefficient(BoundaryElementNode)*EXP(-1.0D00*TempHom*ThermalCoefficient)
  IF (basalSlipCoefficient < TemperateSlipCoefficient(BoundaryElementNode)) THEN
     WRITE(Message,'(A,e10.4)') 'Applied lower threshold of ', &
          TemperateSlipCoefficient(BoundaryElementNode)
     CALL INFO('iceproperties (basalSlip)',Message,Level=9)
     basalSlipCoefficient =  TemperateSlipCoefficient(BoundaryElementNode)
  ELSE IF(basalSlipCoefficient > SlipCoefficientMax) THEN
     WRITE(Message,'(A,e10.4)') 'Applied upper threshold of ', SlipCoefficientMax
     CALL INFO('iceproperties (basalSlip)',Message,Level=9)
     basalSlipCoefficient = SlipCoefficientMax
  END IF
  !------------------------------
  ! clean up
  !------------------------------
  DEALLOCATE(TemperateSlipCoefficient, PressureMeltingPoint)
END FUNCTION basalSlip


!**************************************************************************
!*
!*  computes normal heat flux for melting boudnary
!*
!**************************************************************************
RECURSIVE SUBROUTINE getNormalHeatFlux( Model,Solver,Timestep,TransientSimulation)
  USE DefUtils
  USE Materialmodels
  !-----------------------------------------------------------
  IMPLICIT NONE
  !------------ external variables ---------------------------
  TYPE(Model_t)  :: Model
  TYPE(Solver_t), TARGET :: Solver
  LOGICAL :: TransientSimulation
  REAL(KIND=dp) :: Timestep
  !------------------------------------------------------------------------------
  !    Local variables
  !------------------------------------------------------------------------------
  TYPE(Solver_t), POINTER :: PointerToSolver
  TYPE(Nodes_t) :: Nodes
  TYPE(Element_t),POINTER :: Element
  TYPE(Variable_t), POINTER :: HFSol, HeatFluxSol
  TYPE(ValueList_t), POINTER :: Equation,Material,SolverParams,BodyForce,BC,Constants
  CHARACTER(LEN=MAX_NAME_LEN) :: TempName, SolverName
  INTEGER :: i, j, k, l, t, N, BN,  M, DIM, istat, material_id, HeatFluxDOFs
  INTEGER, POINTER :: HeatFluxPerm(:), HFPerm(:),BoundaryReorder(:)
  INTEGER, ALLOCATABLE :: contributedElements(:)
  REAL(KIND=dp) ::  U, V, W, SqrtElementMetric
  REAL(KIND=dp), ALLOCATABLE ::  Basis(:),dBasisdx(:,:), ddBasisddx(:,:,:), heatFlux(:,:)
  REAL(KIND=dp), POINTER :: HF(:), HeatFluxVar(:), BoundaryNormals(:,:), BoundaryTangent1(:,:), BoundaryTangent2(:,:)
  LOGICAL :: GotIt, AllocationsDone=.FALSE., stat

  SAVE AllocationsDone, SolverName, DIM,  heatFlux, &
       Basis, dBasisdx, ddBasisddx, Nodes, &
       BoundaryNormals, BoundaryTangent1, BoundaryTangent2, BoundaryReorder, BN


  WRITE(SolverName, '(A)') 'IceFlowProperties (getNormalHeatFlux))'
  !-----------------------------------------------------------------------
  ! get solver variable
  !-----------------------------------------------------------------------
  HFSol => Solver % Variable
  IF (.NOT.ASSOCIATED(HFSol)) THEN
     CALL FATAL(SolverName,'No variable associated')
  END IF
  HFPerm  => HFSol % Perm
  HF => HFSol % Values  
  HF = 0.0D00
  !-----------------------------------------------------------------------
  ! Allocations
  !-----------------------------------------------------------------------
  IF ( .NOT.AllocationsDone .OR. Solver % Mesh % Changed ) THEN
     DIM = CoordinateSystemDimension()
     N = Solver % Mesh % MaxElementNodes
     M = Model % Mesh % NumberOfNodes        
     IF ( AllocationsDone ) &
          DEALLOCATE( &
          Basis, &
          dBasisdx, &
          ddBasisddx, &
          heatFlux,&
          Nodes % x, &
          Nodes % y, &
          Nodes % z, &
          BoundaryReorder, &
          BoundaryNormals, &
          BoundaryTangent1, &
          BoundaryTangent2)
     CALL CheckNormalTangentialBoundary( Model, &
          'Basal Melting', BN, &
          BoundaryReorder, BoundaryNormals, BoundaryTangent1, &
          BoundaryTangent2, DIM )
     WRITE(Message,'(A,i6)') &
          'Number of boundary nodes on boundaries associated with melting:',&
          BN
     CALL INFO(SolverName,Message,Level=3)
     IF (BN > 0) THEN
        CALL AverageBoundaryNormals(Model, &
             'Basal Melting', BN, &
             BoundaryReorder, BoundaryNormals, BoundaryTangent1, &
             BoundaryTangent2, DIM )
     END IF
     ALLOCATE( &
          Basis(N), &
          dBasisdx(N,3), &
          ddBasisddx(N,3,3), &
          heatFlux(M,3),& 
          Nodes % x(N), &
          Nodes % y(N), &
          Nodes % z(N), &
          STAT = istat)
     AllocationsDone = .TRUE.
  END IF

  !-------------------------------------------------------------------------
  ! Calculate element-wise heat flux contributions for each point of element
  !-------------------------------------------------------------------------
  heatFlux = 0.0D00  
  IF (BN > 0) THEN
     DO t=1,Solver % NumberOfActiveElements
        Element => GetActiveElement(t)
        Model % CurrentElement => Element
        N = Element % Type % NumberOfNodes
        CALL GetElementNodes( Nodes )
        !-------------
        ! get material
        !-------------
        Material => GetMaterial()
        IF (.NOT.ASSOCIATED(Material)) THEN
           WRITE (Message,'(A,I3)') 'No Material found for boundary element no. ', t
           CALL FATAL(SolverName,Message)
        ELSE
           material_id = GetMaterialId( Element, GotIt)
           IF(.NOT.GotIt) THEN
              WRITE (Message,'(A,I3)') 'No Material ID found for boundary element no. ', t
              CALL FATAL(SolverName,Message)
           END IF
        END IF
        !-----------------------------------------------------------------------
        ! get temperature variable
        !-----------------------------------------------------------------------
        TempName =  GetString(Material ,'Temperature Name', GotIt)
        IF (.NOT.GotIt) THEN
           CALL FATAL(SolverName,'No Temperature Name found in Material section')
        ELSE
           WRITE(Message,'(a,a)') 'Variable Name for temperature: ', TempName
           CALL INFO(SolverName,Message,Level=21)
        END IF
        !-----------------------------------------------------------------------
        ! get flux
        !-----------------------------------------------------------------------
        HeatFluxSol => VariableGet( Solver % Mesh % Variables,  TRIM(TempName) // ' Internal Heat Flux')
        IF ( ASSOCIATED( HeatFluxSol ) ) THEN
           HeatFluxPerm => HeatFluxSol % Perm
           HeatFluxDOFs =  HeatFluxSol % DOFs
           HeatFluxVar  => HeatFluxSol % Values
        ELSE
           CALL FATAL(SolverName,'No >>' //  TRIM(TempName) // ' Internal Heat Flux<< associated')
        END IF
        DO i=1,N
           k = (Element % NodeIndexes(i))
           l = HeatFluxPerm(k)
           IF (l>0) THEN
              SELECT CASE( HeatFluxDOFs )
              CASE(3)
                 heatFlux(k,1) =  HeatFluxVar(HeatFluxDOFs*l-2)
                 heatFlux(k,2) =  HeatFluxVar(HeatFluxDOFs*l-1)
                 heatFlux(k,3) =  HeatFluxVar(HeatFluxDOFs*l)
              CASE(2)
                 heatFlux(k,1) =  HeatFluxVar(HeatFluxDOFs*l-1)
                 heatFlux(k,2) =  HeatFluxVar(HeatFluxDOFs*l)
                 heatFlux(k,3) =  0.0D00
              END SELECT
           END IF
        END DO
     END DO

     !----------------------------------------------------------
     ! averaged internal flux for each point at active boundaries
     !----------------------------------------------------------
     DO i=1,Solver % Mesh % NumberOfNodes
        j = BoundaryReorder(i)        
        IF (j>0) THEN
            HF(HFPerm(i)) = SUM(heatFlux(i,1:DIM)*BoundaryNormals(j,1:DIM))
        ELSE
           HF(HFPerm(i)) = 0.0D00
        END IF
     END DO
  END IF

END SUBROUTINE getNormalHeatFlux

!****************************************************************************************************************
!*
!* viscosity factor as a function of temperature
!*
!****************************************************************************************************************
FUNCTION getViscosityFactor( Model, n, temperature ) RESULT(visFact)
  USE types
  USE CoordinateSystems
  USE SolverUtils
  USE ElementDescription
  USE DefUtils
!-----------------------------------------------------------
  IMPLICIT NONE
!------------ external variables ---------------------------
  TYPE(Model_t) :: Model
  INTEGER :: n
  REAL(KIND=dp) :: temperature, visFact
!------------ internal variables----------------------------
  TYPE(ValueList_t), POINTER :: Material
  INTEGER :: DIM,nMax,i,j,body_id,material_id,elementNodes,nodeInElement,istat
  REAL(KIND=dp) ::&
       rateFactor, aToMinusOneThird, gasconst, temphom
  REAL(KIND=dp), POINTER :: Hwrk(:,:,:)
  REAL (KIND=dp), ALLOCATABLE :: activationEnergy(:,:), arrheniusFactor(:,:),&
       enhancementFactor(:), viscosityExponent(:), PressureMeltingPoint(:),&
       Ux(:), Uy(:), Uz(:), LimitTemp(:)
  LOGICAL :: FirstTime = .TRUE., GotIt
  CHARACTER(LEN=MAX_NAME_LEN) :: TempName, SolverName
!------------ remember this -------------------------------
  Save DIM, FirstTime, gasconst, activationEnergy, arrheniusFactor,&
       enhancementFactor, viscosityExponent, Hwrk, PressureMeltingPoint, &
       Ux, Uy, Uz, LimitTemp, SolverName

  !-----------------------------------------------------------
  ! Read in constants from SIF file and do some allocations
  !-----------------------------------------------------------
  IF (FirstTime) THEN
     WRITE(SolverName, '(A)') 'IceFlowProperties (getViscosityFactor)'
     ! inquire coordinate system dimensions  and degrees of freedom from NS-Solver
     ! ---------------------------------------------------------------------------
     DIM = CoordinateSystemDimension()
     ! inquire minimum temperature
     !------------------------- 
     gasconst = ListGetConstReal( Model % Constants,'Gas Constant',GotIt)
     IF (.NOT. GotIt) THEN
        gasconst = 8.314D00 ! m-k-s
        WRITE(Message,'(a,e10.4,a)') 'No entry for Gas Constant (Constants) in input file found. Setting to ',&
             gasconst,' (J/mol)'
        CALL INFO(SolverName, Message, level=4)
     END IF
     nMax = Model % MaxElementNodes
     ALLOCATE(activationEnergy(2,nMax),&
          arrheniusFactor(2,nMax),&
          enhancementFactor(nMax),&
          LimitTemp( nMax),&
          PressureMeltingPoint( nMax ),&
          viscosityExponent(nMax),&
          Ux(nMax),&
          Uy(nMax),&
          Uz(nMax),&
          STAT=istat)
     IF ( istat /= 0 ) THEN
        CALL Fatal(SolverName,'Memory allocation error, Aborting.')
     END IF
     NULLIFY( Hwrk )
     FirstTime = .FALSE.
     CALL Info(SolverName,'Memory allocations done', Level=3)
  END IF
  !-------------------------------------------
  ! get element properties
  !-------------------------------------------   
  body_id = Model % CurrentElement % BodyId
  material_id = ListGetInteger(Model % Bodies(body_id) % Values, 'Material', GotIt)
  IF (.NOT.GotIt) CALL FATAL(SolverName,'No Material ID found')
  elementNodes = Model % CurrentElement % Type % NumberOfNodes
  IF (.NOT. GotIt) THEN
     WRITE(Message,'(a,I2,a,I2,a)') 'No material id for current element of node ',n,', body ',body_id,' found'
     CALL FATAL(SolverName, Message)
  END IF
  DO nodeInElement=1,elementNodes
     IF ( N == Model % CurrentElement % NodeIndexes(nodeInElement) ) EXIT
  END DO
  Material => Model % Materials(material_id) % Values
  IF (.NOT.ASSOCIATED(Material)) THEN 
     WRITE(Message,'(a,I2,a,I2,a)') 'No Material for current element of node ',n,', body ',body_id,' found'
     CALL FATAL(SolverName,Message)
  END IF
  !-------------------------------------------
  ! get material properties
  !-------------------------------------------
  ! activation energies
  !--------------------
  CALL ListGetRealArray( Material,'Activation Energies',Hwrk,elementNodes, &
       Model % CurrentElement % NodeIndexes, GotIt )
  IF (.NOT. GotIt) THEN
     WRITE(Message,'(a,I2,a,I2)') 'No Value for Activation Energy  found in Material ', material_id,' for node ', n
     CALL FATAL(SolverName,Message)
  END IF
  IF ( SIZE(Hwrk,2) == 1 ) THEN
     DO i=1,MIN(3,SIZE(Hwrk,1))
        activationEnergy(i,1:elementNodes) = Hwrk(i,1,1:elementNodes)
     END DO
  ELSE
     WRITE(Message,'(a,I2,a,I2)') 'Incorrect array size for Activation Energy in Material ', material_id,' for node ', n
     CALL FATAL(SolverName,Message)
  END IF
  ! Arrhenius Factors
  !------------------
  CALL ListGetRealArray( Material,'Arrhenius Factors',Hwrk,elementNodes, &
       Model % CurrentElement % NodeIndexes, GotIt )
  IF (.NOT. GotIt) THEN
     WRITE(Message,'(a,I2,a,I2)') 'No Value for Arrhenius Factors  found in Material ', material_id,' for node ', n
     CALL FATAL(SolverName,Message)
  END IF
  IF ( SIZE(Hwrk,2) == 1 ) THEN
     DO i=1,MIN(3,SIZE(Hwrk,1))
        arrheniusFactor(i,1:elementNodes) = Hwrk(i,1,1:elementNodes)
     END DO
  ELSE
     WRITE(Message,'(a,I2,a,I2)') 'Incorrect array size for Arrhenius Factors in Material ', material_id,' for node ', n
     CALL FATAL(SolverName,Message)
  END IF
  ! Enhancement Factor
  !-------------------
  enhancementFactor(1:elementNodes) = ListGetReal( Material,'Enhancement Factor', elementNodes, &
       Model % CurrentElement % NodeIndexes, GotIt )
  IF (.NOT. GotIt) THEN
     enhancementFactor(1:elementNodes) = 1.0D00
     WRITE(Message,'(a,I2,a,I2,a)') 'No Enhancement Factor found in Material ', material_id,' for node ', n, '.setting E=1'
     CALL INFO(SolverName, Message, level=9)
  END IF
  ! Threshold temperature for switching activation energies and Arrhenius factors
  !------------------------------------------------------------------------------
  LimitTemp(1:elementNodes) = ListGetReal( Material,'Limit Temperature', elementNodes, &
       Model % CurrentElement % NodeIndexes, GotIt )
  IF (.NOT. GotIt) THEN
     LimitTemp(1:elementNodes) = -1.0D01
     WRITE(Message,'(a,I2,a,I2,a)') 'No keyword >Limit Temperature< found in Material ',&
          material_id,' for node ', n, '.setting to -10'
     CALL INFO(SolverName, Message, level=9)
  END IF
  ! Viscosity Exponent
  !-------------------
  viscosityExponent(1:elementNodes) = ListGetReal( Material,'Viscosity Exponent', elementNodes, &
       Model % CurrentElement % NodeIndexes, GotIt )
  IF (.NOT. GotIt) THEN
     viscosityExponent(1:elementNodes) = 1.0D00/3.0D00
     WRITE(Message,'(a,I2,a,I2,a)') 'No Viscosity Exponent found in Material ', material_id,' for node ', n, '.setting k=1/3'
     CALL INFO(SolverName, Message, level=9)
  END IF
  ! Pressure Melting Point and homologous temperature
  !--------------------------------------------------
  TempName =  GetString(Material  ,'Temperature Name', GotIt)
  IF (.NOT.GotIt) CALL FATAL(SolverName,'No Temperature Name found')
  PressureMeltingPoint(1:elementNodes) =&
       ListGetReal( Material, TRIM(TempName) // ' Upper Limit',&
       elementNodes, Model % CurrentElement % NodeIndexes, GotIt )
  IF (.NOT.GotIt) THEN
     temphom = 0.0d00
     WRITE(Message,'(A,A,A,i3,A)') 'No entry for ',TRIM(TempName) // ' Upper Limit',&
          ' found in material no. ', material_id,'. Using 273.16 K.'
     CALL WARN(SolverName,Message)
  ELSE
     temphom = MIN(temperature - PressureMeltingPoint(nodeInElement), 0.0d00)
  END IF
  !-------------------------------------------
  ! homologous Temperature is below 10 degrees
  !-------------------------------------------
  IF (temphom < LimitTemp(nodeInElement)) THEN
     i=1
     !-------------------------------------------
     ! homologous Temperature is above 10 degrees
     !-------------------------------------------
  ELSE
     i=2
  END IF
  rateFactor =&
       arrheniusFactor(i,nodeInElement)*exp(-1.0D00*activationEnergy(i,nodeInElement)/(gasconst*(2.7316D02 + temphom)))
  visFact =(2.0D00*enhancementFactor(nodeInElement)*rateFactor)**(-1.0D00*viscosityExponent(nodeInElement))
END FUNCTION getViscosityFactor



!****************************************************************************************************************
!*
!* fluidity as a function of homologous temperature
!*
!****************************************************************************************************************
FUNCTION getFluidity( Model, n, temperature ) RESULT(fluidity)
  USE types
  USE CoordinateSystems
  USE SolverUtils
  USE ElementDescription
  USE DefUtils
!-----------------------------------------------------------
  IMPLICIT NONE
!------------ external variables ---------------------------
  TYPE(Model_t) :: Model
  INTEGER :: n
  REAL(KIND=dp) :: temperature, fluidity
!------------ internal variables----------------------------
  TYPE(ValueList_t), POINTER :: Material
  INTEGER :: DIM,nMax,i,j,body_id,material_id,elementNodes,nodeInElement,istat
  REAL(KIND=dp) ::&
       rateFactor, aToMinusOneThird, gasconst, temphom
  REAL(KIND=dp), POINTER :: Hwrk(:,:,:)
  REAL (KIND=dp), ALLOCATABLE :: activationEnergy(:,:), arrheniusFactor(:,:),&
       enhancementFactor(:), viscosityExponent(:), PressureMeltingPoint(:),&
       LimitTemp(:)
  LOGICAL :: FirstTime = .TRUE., GotIt
  CHARACTER(LEN=MAX_NAME_LEN) :: TempName
  CHARACTER(LEN=MAX_NAME_LEN) :: SolverName
!------------ remember this -------------------------------
  Save DIM, FirstTime, gasconst, activationEnergy, arrheniusFactor,&
       enhancementFactor,  Hwrk, PressureMeltingPoint, LimitTemp, Solvername
!-----------------------------------------------------------
  !-----------------------------------------------------------
  ! Read in constants from SIF file and do some allocations
  !-----------------------------------------------------------
  IF (FirstTime) THEN
     WRITE(SolverName, '(A)') 'IceFlowProperties (getFluidity)'
     ! inquire coordinate system dimensions  and degrees of freedom from NS-Solver
     ! ---------------------------------------------------------------------------
     DIM = CoordinateSystemDimension()
     ! inquire minimum temperature
     !------------------------- 
     gasconst = ListGetConstReal( Model % Constants,'Gas Constant',GotIt)
     IF (.NOT. GotIt) THEN
        gasconst = 8.314D00 ! m-k-s
        WRITE(Message,'(a,e10.4,a)') 'No entry for Gas Constant (Constants) in input file found. Setting to ',&
             gasconst,' (J/mol)'
        CALL INFO(SolverName, Message, level=4)
     END IF
     nMax = Model % MaxElementNodes
     ALLOCATE(activationEnergy(2,nMax),&
          arrheniusFactor(2,nMax),&
          enhancementFactor(nMax),&
          PressureMeltingPoint( nMax ),&
          LimitTemp( nMax ), &
          STAT=istat)
     IF ( istat /= 0 ) THEN
        CALL Fatal(SolverName,'Memory allocation error, Aborting.')
     END IF
     NULLIFY( Hwrk )
     FirstTime = .FALSE.
     CALL Info(SolverName,'Memory allocations done', Level=3)
  END IF
  !-------------------------------------------
  ! get element properties
  !-------------------------------------------   
  body_id = Model % CurrentElement % BodyId
  material_id = ListGetInteger(Model % Bodies(body_id) % Values, 'Material', GotIt)
  IF (.NOT.GotIt) CALL FATAL(SolverName,'No Material ID found')
  elementNodes = Model % CurrentElement % Type % NumberOfNodes
  IF (.NOT. GotIt) THEN
     WRITE(Message,'(a,I2,a,I2,a)') 'No material id for current element of node ',n,', body ',body_id,' found'
     CALL FATAL(SolverName, Message)
  END IF
  DO nodeInElement=1,elementNodes
     IF ( N == Model % CurrentElement % NodeIndexes(nodeInElement) ) EXIT
  END DO
  Material => Model % Materials(material_id) % Values
  IF (.NOT.ASSOCIATED(Material)) THEN 
     WRITE(Message,'(a,I2,a,I2,a)') 'No Material for current element of node ',n,', body ',body_id,' found'
     CALL FATAL(SolverName,Message)
  END IF
  !-------------------------------------------
  ! get material properties
  !-------------------------------------------
  ! activation energies
  !--------------------
  CALL ListGetRealArray( Material,'Activation Energies',Hwrk,elementNodes, &
       Model % CurrentElement % NodeIndexes, GotIt )
  IF (.NOT. GotIt) THEN
     WRITE(Message,'(a,I2,a,I2)') 'No Value for Activation Energy  found in Material ', material_id,' for node ', n
     CALL FATAL(SolverName,Message)
  END IF
  IF ( SIZE(Hwrk,2) == 1 ) THEN
     DO i=1,MIN(3,SIZE(Hwrk,1))
        activationEnergy(i,1:elementNodes) = Hwrk(i,1,1:elementNodes)
     END DO
  ELSE
     WRITE(Message,'(a,I2,a,I2)') 'Incorrect array size for Activation Energy in Material ', material_id,' for node ', n
     CALL FATAL(SolverName,Message)
  END IF
  ! Arrhenius Factors
  !------------------
  CALL ListGetRealArray( Material,'Arrhenius Factors',Hwrk,elementNodes, &
       Model % CurrentElement % NodeIndexes, GotIt )
  IF (.NOT. GotIt) THEN
     WRITE(Message,'(a,I2,a,I2)') 'No Value for Arrhenius Factors  found in Material ', material_id,' for node ', n
     CALL FATAL(SolverName,Message)
  END IF
  IF ( SIZE(Hwrk,2) == 1 ) THEN
     DO i=1,MIN(3,SIZE(Hwrk,1))
        arrheniusFactor(i,1:elementNodes) = Hwrk(i,1,1:elementNodes)
     END DO
  ELSE
     WRITE(Message,'(a,I2,a,I2)') 'Incorrect array size for Arrhenius Factors in Material ', material_id,' for node ', n
     CALL FATAL(SolverName,Message)
  END IF
  ! Enhancement Factor
  !-------------------
  enhancementFactor(1:elementNodes) = ListGetReal( Material,'Enhancement Factor', elementNodes, &
       Model % CurrentElement % NodeIndexes, GotIt )
  IF (.NOT. GotIt) THEN
     enhancementFactor(1:elementNodes) = 1.0D00
     WRITE(Message,'(a,I2,a,I2,a)') 'No Enhancement Factor found in Material ', material_id,' for node ', n, '.setting E=1'
     CALL INFO(SolverName, Message, level=9)
  END IF
  ! Threshold temperature for switching activation energies and Arrhenius factors
  !------------------------------------------------------------------------------
  LimitTemp(1:elementNodes) = ListGetReal( Material,'Limit Temperature', elementNodes, &
       Model % CurrentElement % NodeIndexes, GotIt )
  IF (.NOT. GotIt) THEN
     LimitTemp(1:elementNodes) = -1.0D01
     WRITE(Message,'(a,I2,a,I2,a)') 'No keyword >Limit Temperature< found in Material ',&
          material_id,' for node ', n, '.setting to -10'
     CALL INFO(SolverName, Message, level=9)
  END IF
  ! Pressure Melting Point and homologous temperature
  !--------------------------------------------------
  TempName =  GetString(Material,'Temperature Name', GotIt)
  IF (.NOT.GotIt) CALL FATAL(SolverName,'No Temperature Name found')
  PressureMeltingPoint(1:elementNodes) =&
       ListGetReal( Material, TRIM(TempName) // ' Upper Limit',&
       elementNodes, Model % CurrentElement % NodeIndexes, GotIt )
  IF (.NOT.GotIt) THEN
     temphom = 0.0d00
     WRITE(Message,'(A,A,A,i3,A)') 'No entry for ',TRIM(TempName) // ' Upper Limit',&
          ' found in material no. ', material_id,'. Using 273.16 K.'
     CALL WARN(SolverName,Message)
  ELSE
     temphom = MIN(temperature - PressureMeltingPoint(nodeInElement), 0.0d00)
  END IF
  !-----------------------------------------------------
  ! homologous Temperature is below temperature treshold
  !----------------------------------------------------
  IF (temphom < LimitTemp(nodeInElement))THEN
     i=1
     !-----------------------------------------------------
     ! homologous Temperature is above temperature treshold
     !-----------------------------------------------------
  ELSE
     i=2
  END IF
  rateFactor =&
       arrheniusFactor(i,nodeInElement)*exp(-1.0D00*activationEnergy(i,nodeInElement)/(gasconst*(2.7316D02 + temphom)))
  fluidity = 2.0D00*enhancementFactor(nodeInElement)*rateFactor
END FUNCTION getFluidity


RECURSIVE SUBROUTINE computeMassBalance( Model,Solver,Timestep,TransientSimulation)
  USE DefUtils
  USE Materialmodels
!-----------------------------------------------------------
  IMPLICIT NONE
!------------ external variables ---------------------------
  TYPE(Model_t)  :: Model
  TYPE(Solver_t), TARGET :: Solver
  LOGICAL :: TransientSimulation
  REAL(KIND=dp) :: Timestep
!------------------------------------------------------------------------------
!    Local variables
!------------------------------------------------------------------------------
  TYPE(Variable_t), POINTER :: MBSol,FlowSol,SurfGrad1,SurfGrad2
  TYPE(Element_t),POINTER :: BoundaryElement, ParentElement 
  TYPE(ValueList_t), POINTER :: ParentMaterial
  INTEGER, POINTER :: BoundaryReorder(:), MBPerm(:),FlowPerm(:)
  CHARACTER(LEN=MAX_NAME_LEN) :: SolverName,FlowSolName,ConvectionFlag
  REAL(KIND=dp), POINTER :: MB(:), FlowSolution(:),BoundaryNormals(:,:), BoundaryTangent1(:,:),&
       BoundaryTangent2(:,:)
  REAL(KIND=dp), ALLOCATABLE :: VeloU(:),VeloV(:),VeloW(:)
  INTEGER :: i,j,k,l,N,t,NSDOFs,istat,DIM,body_id,other_body_id,material_id,BN
  LOGICAL, ALLOCATABLE :: Visited(:)
  LOGICAL :: Found, AllocationsDone

  SAVE AllocationsDone,SolverName, DIM,BN,&
       BoundaryNormals,BoundaryTangent1,BoundaryTangent2,BoundaryReorder,&
       VeloU,VeloV,VeloW

  !-----------------------------------------------------------------------
  ! get solver variable
  !-----------------------------------------------------------------------
  WRITE(SolverName, '(A)') 'IceFlowProperties (computeMassBalance))'
  MBSol => Solver % Variable
  IF (.NOT.ASSOCIATED(MBSol)) THEN
     CALL FATAL(SolverName,'No variable associated')
  END IF
  MBPerm  => MBSol % Perm
  MB => MBSol % Values  
  MB = 0.0D00
  !-----------------------------------------------------------------------
  ! get needed external solver variable
  !-----------------------------------------------------------------------
  SurfGrad1 => VariableGet( Model % Variables, 'FreeSurfGrad1', .TRUE. )
  IF (.NOT.ASSOCIATED(SurfGrad1)) THEN
     CALL FATAL(Solvername,' Variable >FreeSurfGrad1< not found')
  END IF
  SurfGrad2 => VariableGet( Model % Variables, 'FreeSurfGrad2', .TRUE. )
  IF (.NOT.ASSOCIATED(SurfGrad2)) THEN
     CALL FATAL(Solvername,' Variable >FreeSurfGrad2< not found')
  END IF  

  !-----------------------------------------------------------------------
  ! Allocations
  !-----------------------------------------------------------------------
  IF ( .NOT.AllocationsDone .OR. Solver % Mesh % Changed ) THEN
     DIM = CoordinateSystemDimension()
     N = Solver % Mesh % MaxElementNodes
     IF ( AllocationsDone ) &
          DEALLOCATE( &
          Visited, &
          BoundaryReorder, &
          BoundaryNormals, &
          BoundaryTangent1, &
          BoundaryTangent2)
     CALL CheckNormalTangentialBoundary( Model, &
          'Mass Balance', BN, &
          BoundaryReorder, BoundaryNormals, BoundaryTangent1, &
          BoundaryTangent2, DIM )
     WRITE(Message,'(A,i6)') &
          'Number of boundary nodes on boundaries associated with mass balance computation:',&
          BN
     CALL INFO(SolverName,Message,Level=3)
     IF (BN > 0) THEN
        CALL AverageBoundaryNormals(Model, &
             'Basal Melting', BN, &
             BoundaryReorder, BoundaryNormals, BoundaryTangent1, &
             BoundaryTangent2, DIM )
        ALLOCATE( &
             Visited(BN), &
             VeloU(N), &
             VeloV(N), &
             VeloW(N), &
             STAT = istat)   
        IF (istat /= 0) &
           CALL FATAL(SolverName,'Allocations failed')
     END IF
     AllocationsDone = .TRUE.
  END IF

  Visited = .FALSE.

  IF (BN > 0) THEN
     DO t=1, Solver % Mesh % NumberOfBoundaryElements
        ! get element information
        !------------------------
        BoundaryElement => GetBoundaryElement(t)

        IF ( (.NOT.ActiveBoundaryElement()) .OR. (GetElementFamily() == 1 )) CYCLE
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
        material_id = ListGetInteger(Model % Bodies(body_id) % Values, 'Material', Found)
        ParentMaterial => Model % Materials(material_id) % Values
        IF (.NOT.ASSOCIATED(ParentMaterial)) THEN
           WRITE(Message,'(A,I10,A)')&
                'No material for Parent Element no. ',&
                ParentElement % ElementIndex, ' found'
           CALL FATAL(SolverName,Message)
        END IF

        N = GetElementNOFNodes(ParentElement)
        VeloU = 0.0d00
        VeloV = 0.0d00
        VeloW = 0.0d00
        k = ListGetInteger( Model % Bodies(ParentElement % BodyId) % Values, 'Equation', &
             minv=1, maxv=Model % NumberOFEquations )
        ConvectionFlag = GetString( Model % Equations(k) % Values, 'Convection', Found ) 
        IF (.NOT.Found) THEN
           WRITE(Message,'(A,I6)') 'No keyword >Convection< found in section Equation ', k
           CALL FATAL(SolverName,Message)
        END IF
        ! constant (i.e., in section Material given) velocity
        !---------------------------------------------------
        IF ( ConvectionFlag == 'constant' ) THEN
           Model % CurrentElement => ParentElement
           VeloU(1:N) = GetReal( ParentMaterial, 'Convection Velocity 1', Found )
           VeloV(1:N) = GetReal( ParentMaterial, 'Convection Velocity 2', Found )
           VeloW(1:N) = GetReal( ParentMaterial, 'Convection Velocity 3', Found )                 
           ! computed velocity
           !------------------
        ELSE IF (( ConvectionFlag == 'computed' )) THEN

           FlowSolName =  GetString( Model % Equations(k) % Values,'Flow Solution Name', Found)
           IF(.NOT.Found) THEN        
              CALL WARN(SolverName,'Keyword >Flow Solution Name< not found in section >Equation<')
              CALL WARN(SolverName,'Taking default value >Flow Solution<')
              WRITE(FlowSolName,'(A)') 'Flow Solution'
           END IF


           FlowSol => VariableGet( Solver % Mesh % Variables, FlowSolName )
           IF ( ASSOCIATED( FlowSol ) ) THEN
              FlowPerm     => FlowSol % Perm
              NSDOFs       =  FlowSol % DOFs
              FlowSolution => FlowSol % Values
           ELSE
              CALL FATAL(SolverName,'No >Flow Solution< associated')
           END IF

           DO i=1,N
              k = FlowPerm(ParentElement % Nodeindexes(i))
              IF ( k > 0 ) THEN
                 SELECT CASE( NSDOFs )
                 CASE(3)
                    VeloU(i) = FlowSolution( NSDOFs*k-2 )
                    VeloV(i) = FlowSolution( NSDOFs*k-1 )
                    VeloW(i) = 0.0D0
                 CASE(4)
                    VeloU(i) = FlowSolution( NSDOFs*k-3 )
                    VeloV(i) = FlowSolution( NSDOFs*k-2 )
                    VeloW(i) = FlowSolution( NSDOFs*k-1 )
                 END SELECT
              END IF

              WRITE(Message,'(a,i5, a, i5)') 'Convection in element ', t, &
                   ' material ',  material_id
           END DO
        END IF


        DO i=1,N
           j = ParentElement % NodeIndexes(i)
           l = BoundaryReorder(j)           
           IF (l <= 0) CYCLE
           IF ( Visited(l) ) CYCLE
           Visited(l) = .TRUE.

           MB(MBPerm(j))  & 
                = SurfGrad1 % Values(SurfGrad1 % Perm(j)) * VeloU(i) &
                + SurfGrad2 % Values(SurfGrad2 % Perm(j)) * VeloV(i) &
                - VeloW(i)
        END DO
     END DO
  ELSE
     CALL WARN(SolverName,'No keyword >Mass Balance< found in any of the boundary conditions')
     CALL WARN(SolverName,'No points selected for mass abalance computation')
     CALL WARN(SolverName,'Are you sure you need that solver?')
  END IF

END SUBROUTINE computeMassBalance

! *****************************************************************************
! calculates the SIA stress components on the free boundaries
!
! *****************************************************************************
RECURSIVE SUBROUTINE getSIAstress( Model,Solver,Timestep,TransientSimulation)
  USE DefUtils
  USE Materialmodels
!-----------------------------------------------------------
  IMPLICIT NONE
!------------ external variables ---------------------------
  TYPE(Model_t)  :: Model
  TYPE(Solver_t), TARGET :: Solver
  LOGICAL :: TransientSimulation
  REAL(KIND=dp) :: Timestep
!------------------------------------------------------------------------------
!    Local variables
!------------------------------------------------------------------------------
  TYPE(Solver_t), POINTER :: PointerToSolver
  TYPE(Nodes_t) :: Nodes
  TYPE(Element_t),POINTER :: BoundaryElement, CurrentElement
  TYPE(Variable_t), POINTER :: SIASol, DepthSol, SurfSol1, SurfSol2
  TYPE(ValueList_t), POINTER :: Equation,Material,SolverParams,BodyForce,BC,Constants
  CHARACTER(LEN=MAX_NAME_LEN) :: TempName, SolverName
  INTEGER :: i, j, k, l, t, N, BN,  M, DIM, SIADOFs, istat, elementNbNodes
  INTEGER, POINTER :: DepthSolPerm(:), SurfSol1Perm(:), SurfSol2Perm(:),&
       SIAPerm(:),BoundaryReorder(:), NodeIndexes(:)
  INTEGER :: material_id, body_id, bf_id
  INTEGER :: NBoundary
  REAL(KIND=dp) ::  U, V, W, SqrtElementMetric, SIAstress(2), HydrostaticPressure, &
                    density, gravity
  REAL(KIND=dp), POINTER :: SIA(:), Depth(:), Surf1(:), Surf2(:), &
       BoundaryNormals(:,:), BoundaryTangent1(:,:), BoundaryTangent2(:,:)
  LOGICAL :: GotIt, FirstTimeAround=.TRUE., stat

  SAVE FirstTimeAround, SolverName, DIM, &
       BoundaryNormals, BoundaryTangent1, BoundaryTangent2, BoundaryReorder, BN

  ! assign solver name for communicative output
  !-----------------------------------------------------------------------
  WRITE(SolverName, '(A)') 'IceFlowproperties (getSIAstress))'
  !-----------------------------------------------------------------------
  ! get solver variable
  !-----------------------------------------------------------------------
  SIASol => Solver % Variable
  IF (.NOT.ASSOCIATED(SIASol)) THEN
     CALL FATAL(SolverName,'No variable associated')
  END IF
  SIAPerm  => SIASol % Perm
  SIA => SIASol % Values
  SIADOFs =  SIASol % DOFs
  SIA = 0.0d00
    
  !-----------------------------------------------------------------------
  ! 
  !-----------------------------------------------------------------------
  IF ( FirstTimeAround .OR. Solver % Mesh % Changed ) THEN
     DIM = CoordinateSystemDimension()
     N = Solver % Mesh % MaxElementNodes
     M = Model % Mesh % NumberOfNodes        
     IF ( .NOT.FirstTimeAround ) &
          DEALLOCATE( &
          BoundaryReorder, &
          BoundaryNormals, &
          BoundaryTangent1, &
          BoundaryTangent2)
     ! check boundaries for calculation of SIA-stress and allocate necessary 
     ! space for averaged Normal and Tangentials
     !-----------------------------------------------------------------------
     CALL CheckNormalTangentialBoundary( Model, &
          'Calc SIA', BN, &
          BoundaryReorder, BoundaryNormals, BoundaryTangent1, &
          BoundaryTangent2, DIM )
     WRITE(Message,'(A,i6)') &
          'Number of boundary nodes on boundaries associated with SIA stresses:',&
          BN
     CALL INFO(SolverName,Message,Level=3)
     ! compute averaged normals and tangentials for  boundaries designated for
     ! SIA-stress boundary elements
     !-----------------------------------------------------------------------
     CALL AverageBoundaryNormals(Model, &
          'Calc SIA', BN, &
          BoundaryReorder, BoundaryNormals, BoundaryTangent1, &
          BoundaryTangent2, DIM )

     ! read in variables for flow depth and for free surface gradients
     !-----------------------------------------------------------------------
     DepthSol => VariableGet( Solver % Mesh % Variables, "Depth" )
     IF ( ASSOCIATED( DepthSol ) ) THEN
        DepthSolPerm => DepthSol % Perm
        Depth => DepthSol % Values
     ELSE
        WRITE(Message, '(A)') 'Could not find surface Gradient 1 field variable '
        CALL FATAL(SolverName,Message)
     END IF

     SurfSol1 => VariableGet( Solver % Mesh % Variables, "FreeSurfGrad1" )
     IF ( ASSOCIATED( SurfSol1 ) ) THEN
        SurfSol1Perm => SurfSol1 % Perm
        Surf1 => SurfSol1 % Values
     ELSE
        WRITE(Message, '(A)') 'Could not find Surface Gradient 1 field variable '
        CALL FATAL(SolverName,Message)
     END IF
     IF (DIM > 2) THEN
        SurfSol2 => VariableGet( Solver % Mesh % Variables, "FreeSurfGrad2" )
        IF ( ASSOCIATED( SurfSol2 ) ) THEN
           SurfSol2Perm => SurfSol2 % Perm
           Surf2 => SurfSol2 % Values
        ELSE
           WRITE(Message, '(A)') 'Could not find Surface Gradient 2 field variable '
           CALL FATAL(SolverName,Message)
        END IF
     END IF

     body_id = -1  
     
     !-----------------------------------------------------------------------
     ! loop over all elements in mesh
     !-----------------------------------------------------------------------
     DO t=1,Solver % NumberOFActiveElements
 
        CurrentElement => GetActiveElement(t)
        IF (ParEnv % myPe .NE. CurrentElement % partIndex) CYCLE
        elementNbNodes = GetElementNOFNodes( CurrentElement)
        NodeIndexes => CurrentElement % NodeIndexes
  
        IF ( CurrentElement % BodyId /= body_id ) THEN
         body_id = CurrentElement % BodyId
       
         !Get BofyForce ID!---------------- 
         bf_id = ListGetInteger( Model % Bodies(body_id) % Values, &              
                  'Body Force', gotIt, 1, Model % NumberOfBodyForces )
 
         material_id = ListGetInteger(Model % Bodies(body_id) % Values, 'Material', GotIt)
         Material => Model % Materials(material_id) % Values
         IF ((.NOT. ASSOCIATED(Material)) .OR. (.NOT. GotIt)) THEN
              WRITE(Message,'(A)') 'No material found for model '
              CALL FATAL('IceFlowproperties (getSIAstress)',Message)
         END IF
   
         density = GetConstReal(Material, 'Density', GotIt)      
         IF(.NOT. GotIt) THEN         
             CALL FATAL('IceFlowproperties (getSIAstress)', & 
                  'Density not found in Material')      
         END IF

         gravity = GetConstReal(Model % BodyForces(bf_id) % Values, & 
                  'Flow Bodyforce 3', GotIt)      
         IF(.NOT. GotIt) THEN         
             CALL FATAL('IceFlowproperties (getSIAstress)', 'Flow Body Force 3 not found')
         END IF
        END IF

        DO i=1,elementNbNodes
        j = BoundaryReorder(NodeIndexes(i)) ! projection from real space to SIA-stress 
        ! boundary space
        k = SIAPerm(NodeIndexes(i)) ! projection from real space to solver matrix coordinate

        ! if boundary element node with SIA-stress condition enabled
        !-----------------------------------------------------------------------
         IF (j > 0) THEN
           HydrostaticPressure = density * gravity * Depth(DepthSolPerm(NodeIndexes(i)))        
           SIAstress(1) =  HydrostaticPressure * Surf1(SurfSol1Perm(NodeIndexes(i)))
           IF (DIM > 2) THEN
              SIAstress(2) = HydrostaticPressure * Surf2(SurfSol2Perm(NodeIndexes(i)))
           ELSE
              SIAstress(2) = 0.0D00
           END IF

           ! vector product between Cauchy-stress tensor and surface normal
           !  = stress vector
           !-----------------------------------------------------------------------
           !Compute first element of the stress vector P*n(1)+txy*n(2) For  DIM=2 (two dimensions)
           !Compute first element of the stress vector P*n(1)+txz*n(3) for  Dim=3 (three-dimensions)
           SIA(SIADOFs*(k-1)+1)= HydrostaticPressure * BoundaryNormals(BoundaryReorder(NodeIndexes(i)),1) &
                + SIAstress(1) * BoundaryNormals(BoundaryReorder(NodeIndexes(i)),DIM)
           !Compute the second element of the stress vector in three-dimension P*n(2)+tyz*n(3)
           SIA(SIADOFs*(k-1)+2)= HydrostaticPressure * BoundaryNormals(BoundaryReorder(NodeIndexes(i)),2) &
                + SIAstress(2) * BoundaryNormals(BoundaryReorder(NodeIndexes(i)),3) ! this line doesn't contribute if 2d
           !If DIM=2, compute the second element of the stress vector in two-dimensions P*n(2)+txy*n(1)
           IF (DIM == 2) &
                SIA(SIADOFs*(k-1)+2)= SIA(SIADOFs*(k-1)+2) &
                + SIAstress(1) * BoundaryNormals(BoundaryReorder(NodeIndexes(i)),1)
           !If DIM=3, compute the third element of the stress vector in three-dimensions P*n(3)+txz*n(1)+tyz*n(2)
           IF (DIM > 2) &
                SIA(SIADOFs*(k-1)+3)= HydrostaticPressure * BoundaryNormals(BoundaryReorder(NodeIndexes(i)),3) &
                + SIAstress(1) * BoundaryNormals(BoundaryReorder(NodeIndexes(i)),1) &
                + SIAstress(2) * BoundaryNormals(BoundaryReorder(NodeIndexes(i)),2)   
         END IF
        END DO
     END DO
     FirstTimeAround = .FALSE.
  END IF
END SUBROUTINE getSIAstress

!*********************************************************************************************************************************
!* Computation of component of the SIA stress vector. This code is  used for a rectangular type geometry where the automatic solution "getSIAstres" does not provide good results at the corner. Here we compute the components Pressure i in the sif file independantly for each side, North, Est, South and West. 
!*
!*********************************************************************************************************************************

!----------------------------------------------------------------------------------------------
!StressSIA: Compute the corresponding shallow-ice shear stress
!----------------------------------------------------------------------------------------------
FUNCTION StressSIA( Model, Node, surfaceGrad) RESULT(shearstress)
  USE Types
  USE CoordinateSystems
  USE SolverUtils
  USe ElementDescription
!-------------------------------------------------------------------
IMPLICIT NONE
!-------------------------external variables-----------------------
TYPE(Model_t) :: Model 
INTEGER :: Node
INTEGER, POINTER :: DepthPermutation(:)
REAL(KIND=dp) :: surfaceGrad, shearstress
TYPE(Variable_t), POINTER :: DepthSolution 

!Find variable for flow depth 
DepthSolution  => VariableGet(Model % Variables, 'Depth', .TRUE.)
DepthPermutation  => DepthSolution % Perm 
shearstress = -918.0D00*9.81D00*DepthSolution % Values(DepthPermutation(Node))*surfaceGrad
END FUNCTION StressSIA

!--------------------------------------------------------------------------------------------

FUNCTION invStressSIA( Model, Node, surfaceGrad) RESULT(shearstress)
  USE Types
  USE CoordinateSystems
  USE SolverUtils
  USe ElementDescription
!-------------------------------------------------------------------
IMPLICIT NONE
!-------------------------external variables-----------------------
TYPE(Model_t) :: Model 
INTEGER :: Node
INTEGER, POINTER :: DepthPermutation(:)
REAL(KIND=dp) :: surfaceGrad, shearstress
TYPE(Variable_t), POINTER :: DepthSolution 

!Find variable for flow depth 
DepthSolution  => VariableGet(Model % Variables, 'Depth', .TRUE.)
DepthPermutation  => DepthSolution % Perm 
shearstress = 918.0D00*9.81D00*DepthSolution % Values(DepthPermutation(Node))*surfaceGrad
END FUNCTION invStressSIA
