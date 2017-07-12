RECURSIVE SUBROUTINE getSIAStressBC( Model,Solver,Timestep,TransientSimulation)

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
  TYPE(Variable_t), POINTER :: SIASol, DepthSol, SurfSol1, SurfSol2, NormalVar
  TYPE(ValueList_t), POINTER :: Equation,Material,SolverParams,BodyForce,BC,Constants
  CHARACTER(LEN=MAX_NAME_LEN) :: TempName, SolverName
  INTEGER :: i, j, k, l, t, N,  BN, M, DIM, SIADOFs, istat, elementNbNodes
  INTEGER, POINTER :: DepthSolPerm(:), SurfSol1Perm(:), SurfSol2Perm(:),&
       SIAPerm(:), NodeIndexes(:), NormalPerm(:), BoundaryReorder(:)
  INTEGER :: material_id, body_id, bf_id
  INTEGER :: NBoundary
  REAL(KIND=dp) ::  U, V, W, SqrtElementMetric, SIAstress(2), HydrostaticPressure, &
                    density, gravity
  REAL(KIND=dp), POINTER :: SIA(:), Depth(:), Surf1(:), Surf2(:), NormalValues(:), BoundaryNormals(:,:), &
                            BoundaryTangent1(:,:), BoundaryTangent2(:,:)
  REAL (KIND=dp), ALLOCATABLE :: normal(:)
  LOGICAL :: GotIt, FirstTimeAround=.TRUE., stat

  SAVE FirstTimeAround, SolverName, DIM, normal, BoundaryNormals, &
       BoundaryTangent1, BoundaryTangent2, BoundaryReorder, BN


  ! assign solver name for communicative output
  !-----------------------------------------------------------------------
  WRITE(SolverName, '(A)') 'getSIAStressBC'
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
          BoundaryTangent2, &
          normal)

     IF ((DIM == 2).OR.(DIM == 3))  THEN
            ALLOCATE(normal(DIM))
     ELSE
            CALL FATAL(SolverName, 'Bad dimension of the problem')
     END IF

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

     ! Read in variables for flow depth and for free surface gradients
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

     ! Read normals
     !-----------------------------------------------------------------------
     NormalVar =>  VariableGet(Model % Variables,'Normal Vector')
     IF ( ASSOCIATED( NormalVar ) ) THEN
        NormalPerm => NormalVar % Perm
        NormalValues => NormalVar % Values
     ELSE
        CALL FATAL(SolverName, 'Need ComputeNormal Solver, Normal Vector not found')
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
              CALL FATAL(SolverName,Message)
         END IF
   
         density = GetConstReal(Material, 'Density', GotIt)      
         IF(.NOT. GotIt) THEN         
             CALL FATAL('IceFlowproperties (getSIAstress)', & 
                  'Density not found in Material')      
         END IF

         gravity = GetConstReal(Model % BodyForces(bf_id) % Values, & 
                  'Flow Bodyforce 3', GotIt)      
         IF(.NOT. GotIt) THEN         
             CALL FATAL(SolverName, 'Flow Body Force 3 not found')
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

           DO l=1, DIM
            normal(l) = NormalValues(DIM*(NormalPerm((NodeIndexes(i)))-1) + l)
           END DO

           ! vector product between Cauchy-stress tensor and surface normal
           !  = stress vector
           !-----------------------------------------------------------------------
           !Compute first element of the stress vector P*n(1)+txy*n(2) For  DIM=2 (two dimensions)
           !Compute first element of the stress vector P*n(1)+txz*n(3) for  Dim=3 (three-dimensions)
           SIA(SIADOFs*(k-1)+1)= HydrostaticPressure * normal(1) &
                + SIAstress(1) * normal(DIM)
           !Compute the second element of the stress vector in three-dimension P*n(2)+tyz*n(3)
           SIA(SIADOFs*(k-1)+2)= HydrostaticPressure * normal(2) &
                + SIAstress(2) * normal(3)! this line doesn't contribute if 2d
           !If DIM=2, compute the second element of the stress vector in two-dimensions P*n(2)+txy*n(1)
           IF (DIM == 2) &
                SIA(SIADOFs*(k-1)+2)= SIA(SIADOFs*(k-1)+2) &
                + SIAstress(1) * normal(1)
           !If DIM=3, compute the third element of the stress vector in three-dimensions P*n(3)+txz*n(1)+tyz*n(2)
           IF (DIM > 2) &
                SIA(SIADOFs*(k-1)+3)= HydrostaticPressure * normal(3) &
                + SIAstress(1) * normal(1) &
                + SIAstress(2) * normal(2)  
         END IF
        END DO
     END DO
     FirstTimeAround = .FALSE.
  END IF

END SUBROUTINE getSIAStressBC 
