!/*****************************************************************************/
! *
! *  Elmer, A Finite Element Software for Multiphysical Problems
! *
! *  Copyright 1st April 1995 - , CSC - IT Center for Science Ltd., Finland
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
! *************************************************************************************/
!
!/**************************************************************************************
! *
! *  Authors: Thomas Zwinger, Peter Råback, Juha Ruokolainen, Mikko Lyly, Hakime Seddik
! *  Email:   Thomas.Zwinger@csc.fi
! *  Web:     http://www.csc.fi/elmer
! *  Address: CSC - IT Center for Science Ltd.
! *           Keilaranta 14
! *           02101 Espoo, Finland 
! *
! *  Original Date: 17 May 2002
! *  Added limiters: 09 Jan 2008
! *
! *************************************************************************************/


!-----------------------------------------------------------------------------
!>  Solver for free surface evolution in 2d and 3d flows
!>  with or without surface flux, and upper and lower limiters.
!> \ingroup Solvers
!-----------------------------------------------------------------------------
SUBROUTINE FreeSurfaceSolverSIA( Model,Solver,dt,TransientSimulation )
  USE DefUtils
  USE Differentials
  USE MaterialModels
  IMPLICIT NONE

  !------------------------------------------------------------------------------
  !    external variables
  !------------------------------------------------------------------------------
  TYPE(Model_t) :: Model
  TYPE(Solver_t):: Solver
  REAL(KIND=dp) :: dt
  LOGICAL :: TransientSimulation
  !------------------------------------------------------------------------------
  !    Local variables
  !------------------------------------------------------------------------------

  LOGICAL ::&
       firstTime=.TRUE., Found, AllocationsDone = .FALSE., stat, &
       NeedOldValues, LimitDisp, &
       NormalFlux = .TRUE., SubstantialSurface = .TRUE.,&
       UseBodyForce = .TRUE., ApplyDirichlet=.FALSE.,  ALEFormulation=.FALSE., &
       scaling=.FALSE.
  LOGICAL, ALLOCATABLE ::  LimitedSolution(:,:), ActiveNode(:,:)

  INTEGER :: & 
       i,j,K,L, p, q, R, t,N,NMAX,MMAX,nfamily, deg, Nmatrix,&
       edge, bf_id,DIM,istat,LocalNodes,nocorr,&
       NSDOFs,NonlinearIter,iter, numberofsurfacenodes
  INTEGER, POINTER ::&
       FreeSurfPerm(:), NodeIndexes(:), EdgeMap(:,:), FreeSurfGrad1Perm(:), &
       FreeSurfGrad2Perm(:), BetaPerm(:), HeightPerm(:), &
       IntegralPerm(:), BasalMeltingPerm(:)


  REAL(KIND=dp) :: &
       at,st,totat,totst,CPUTime,Norm,PrevNorm,LocalBottom, cv, &
       Relax, MaxDisp, maxdh,LinearTol,NonlinearTol,RelativeChange,&
       smallestpossiblenumber, rr, ss, dragBeta, &
       scaleFactor

  REAL(KIND=dp), POINTER :: ForceVector(:), FreeSurf(:), PreFreeSurf(:,:), &
       PrevFlowSol(:,:), PointerToResidualVector(:), FreeSurfGrad1Values(:), &
       FreeSurfGrad2Values(:), BetaValues(:), HeightValues(:), &
       IntegralValues(:), BasalMeltingValues(:)

  REAL(KIND=dp), ALLOCATABLE :: ResidualVector(:), OldFreeSurf(:), &
       STIFF(:,:),SourceFunc(:),FORCE(:), TimeForce(:), &
       MASS(:,:), Flux(:,:), LowerLimit(:), UpperLimit(:), &
       OldValues(:), OldRHS(:),StiffVector(:),MeshVelocity(:,:), ElemFreeSurf(:), &
       Density(:), Gravity(:), Thickness(:), Gradh(:), IntegrationPart(:), &
       Drag(:), Melting(:)

  CHARACTER(LEN=MAX_NAME_LEN)  :: SolverName, VariableName, EquationName

  TYPE(Nodes_t)   :: ElementNodes
  TYPE(Element_t),POINTER :: CurrentElement
  TYPE(Variable_t), POINTER :: VarSurfResidual, FreeSurfGrad1, FreeSurfGrad2, Beta, &
                               Height, Integral, BasalMelting
  TYPE(ValueList_t), POINTER :: BodyForce, SolverParams, Material, Equation
  TYPE(Matrix_t), POINTER :: Systemmatrix
  !-----------------------------------------------------------------------------
  !      remember these variables
  !----------------------------------------------------------------------------- 
  SAVE STIFF, MASS, SourceFunc, FORCE, &
       ElementNodes, AllocationsDone, OldFreeSurf, TimeForce, &
       ElemFreeSurf, Flux, SubstantialSurface, &
       NormalFlux, UseBodyForce, LimitedSolution, LowerLimit, &
       UpperLimit, ActiveNode, OldValues, OldRHS, &
       ResidualVector, StiffVector, MeshVelocity, Density, Gravity, Thickness, &
       Gradh, IntegrationPart, Drag, Melting
  !------------------------------------------------------------------------------
  !    Get variables for the solution
  !------------------------------------------------------------------------------
  FreeSurf     => Solver % Variable % Values     ! Nodal values for free surface displacement
  IF (.NOT.ASSOCIATED(FreeSurf)) CALL FATAL(SolverName,'Variable values not associated')
  FreeSurfPerm => Solver % Variable % Perm       ! Permutations for free surface displacement
  PreFreeSurf  => Solver % Variable % PrevValues ! Nodal values for free surface displacement

  !------------------------------------------------------------------------------
  !    Get variabel/solver name
  !------------------------------------------------------------------------------
  VariableName = TRIM(Solver % Variable % Name)
  SolverName = 'FreeSurfaceSolver ('// TRIM(Solver % Variable % Name) // ')'
  !------------------------------------------------------------------------------
  !    if this partition (or the serial problem) has no free surface,
  !    then nothing to be doneGet variabel/solver name
  !------------------------------------------------------------------------------
  IF ( COUNT(FreeSurfPerm/=0)==0) THEN
     IF (ParEnv % PEs > 1) THEN
        WRITE(Message,'(A,I6,A)'), 'Partition ', ParEnv % myPE, 'has no free surface'
        CALL WARN(SolverName,Message)
        CALL WARN(SolverName,'This is not good for load balance!')
     ELSE
        CALL WARN(SolverName,'A serial run without a free surface, but the solver switched in - weird!')
     END IF
     RETURN
  END IF
  !------------------------------------------------------------------------------
  !    Get constants and solver params
  !------------------------------------------------------------------------------
  DIM = CoordinateSystemDimension()
  smallestpossiblenumber = TINY(smallestpossiblenumber)
  SolverParams => GetSolverParams()

  SystemMatrix => Solver % Matrix
  ForceVector => Solver % Matrix % RHS

  cv = GetConstReal( SolverParams, 'Velocity Implicity', Found)
  IF(.NOT. Found) cv = 1.0_dp 
  WRITE(Message,'(a,F8.2)') 'Velocity implicity (1=fully implicit)=', cv
  CALL INFO(SolverName, Message, Level=4)

  LinearTol = GetConstReal( SolverParams, &
       'Linear System Convergence Tolerance',    Found )
  IF ( .NOT.Found ) THEN
     CALL FATAL(SolverName, 'No >Linear System Convergence Tolerance< found')
  END IF
  NonlinearTol  = GetConstReal( SolverParams, &
       'Nonlinear System Convergence Tolerance',    Found )
  NonlinearIter = GetInteger(   SolverParams, &
       'Nonlinear System Max Iterations', Found )
  IF ( .NOT.Found ) NonlinearIter = 1

  MaxDisp = GetConstReal( SolverParams, 'Maximum Displacement', LimitDisp)

  Relax = GetCReal( SolverParams, 'Relaxation Factor', Found)
  IF(.NOT. Found) Relax = 1.0_dp
  NeedOldValues = Found .OR. LimitDisp 

  ApplyDirichlet = GetLogical( SolverParams, &
       'Apply Dirichlet', Found)
  IF ( .NOT.Found ) THEN
     ApplyDirichlet = .FALSE.
     CALL INFO(SolverName, 'No keyword >Apply Dirichlet< found. No limitation of solution',Level=1)
  ELSE
     IF (ApplyDirichlet) THEN
        CALL INFO(SolverName, 'Using Dirichlet method for limitation',Level=1)
        IF (NonlinearIter < 2) THEN
           CALL WARN(SolverName, 'Keyword >Apply Dirichlet< set, but >Nonlinear System Max Iterations< set to lower than 2')
        END IF
     ELSE
        CALL INFO(SolverName, 'No limitation of solution',Level=1)
     END IF
  END IF

  ALEFormulation = GetLogical( SolverParams, &
       'ALE Formulation', Found)
  IF ( .NOT.Found ) THEN
     ALEFormulation = .FALSE.
  END IF
  IF (ALEFormulation) THEN 
     CALL INFO(SolverName, 'Using horizontal ALE Formulation',Level=1)
  ELSE
     CALL INFO(SolverName, 'Using horizontal Eulerian Formulation',Level=1)     
  END IF

  WRITE(Message,'(a,i1,a,i1)') 'DIM=', DIM
  CALL INFO( SolverName, Message, level=4) 

  !------------------------------------------------------------------------------
  !    Allocate some permanent storage, this is done first time only
  !------------------------------------------------------------------------------

  IF ( .NOT. AllocationsDone .OR. Solver % Mesh % Changed) THEN
     NMAX = Model % MaxElementNodes
     MMAX = Model % Mesh % NumberOfNodes 
     K = SIZE( SystemMatrix % Values )
     L = SIZE( SystemMatrix % RHS )

     IF ( AllocationsDone ) THEN
        DEALLOCATE( ElementNodes % x,    &
             ElementNodes % y,    &
             ElementNodes % z,    &
             TimeForce,        &
             FORCE,    &
             STIFF, &
             MASS,  &
             MeshVelocity, &
             Flux, &
             ElemFreeSurf,&
             SourceFunc, &
             LowerLimit,                      &
             UpperLimit, &
             LimitedSolution,  &
             ActiveNode,                      & 
             ResidualVector, &
             StiffVector,  &
             OldValues, &
             OldRHS, &
             Density, &
             Gravity, &
             Thickness, &
             Gradh, &
             IntegrationPart, &
             Drag, &
             Melting)
     END IF

    Nmatrix = NMAX

     ALLOCATE( ElementNodes % x( NMAX ),    &
          ElementNodes % y( NMAX ),    &
          ElementNodes % z( NMAX ),    &
          TimeForce( Nmatrix ),        &
          FORCE( Nmatrix ),    &
          STIFF( Nmatrix, Nmatrix ), &
          MASS( Nmatrix, Nmatrix ),  &
          MeshVelocity( 3,NMAX ), &
          Flux( 3, NMAX), &
          ElemFreeSurf( NMAX ),&
          SourceFunc( NMAX ), &
          LowerLimit( MMAX ), &
          UpperLimit( MMAX ), &
          LimitedSolution( MMAX, 2 ),  &
          ActiveNode( MMAX, 2 ),                      &  
          ResidualVector( L ),                    &
          StiffVector( L ), &
          OldValues( K ), &
          OldRHS( L ), &
          Density(NMax), &
          Gravity(NMax), &
          Thickness(NMax), & 
          Gradh(NMax), &
          IntegrationPart(NMAx), &
          Drag(NMAx), &
          Melting(NMax), &
          STAT=istat )
     IF ( istat /= 0 ) THEN
        CALL Fatal(SolverName,'Memory allocation error, Aborting.')
     END IF

     IF(NeedOldValues) THEN
        ALLOCATE(OldFreeSurf(SIZE(FreeSurf)), STAT=istat)

        IF ( istat /= 0 ) THEN
           CALL Fatal(SolverName,'Memory allocation error, Aborting.')
        END IF
     END IF

     CALL INFO(SolverName,'Memory allocations done', Level=4)
     AllocationsDone = .TRUE.
     ActiveNode = .FALSE.
     ResidualVector = 0.0_dp
  END IF

  !                     from previous timestep
  IF( NeedOldValues) THEN
     OldFreeSurf = FreeSurf
  END IF

  !------------------------------------------------------------------------------
  !    Get variables for the residual
  !------------------------------------------------------------------------------
  VarSurfResidual => VariableGet( Model % Mesh % Variables, TRIM(VariableName) // ' Residual' )
  IF (.NOT.ASSOCIATED(VarSurfResidual)) THEN
     WRITE(Message,'(A)') '>' // TRIM(VariableName) // ' Residual< not associated'
     CALL FATAL( SolverName, Message)
  END IF
  PointerToResidualVector => VarSurfResidual % Values

  FreeSurfGrad1 => VariableGet( Model % Variables, 'FreeSurfGrad1' )
  IF ( ASSOCIATED( FreeSurfGrad1 ) ) THEN
        FreeSurfGrad1Perm    =>  FreeSurfGrad1 % Perm
        FreeSurfGrad1Values  =>  FreeSurfGrad1 % Values
  ELSE
        CALL FATAL(SolverName,'No variable FreeSurfGrad1 found')  
  END IF

  FreeSurfGrad2 => VariableGet( Model % Variables, 'FreeSurfGrad2' )
  IF ( ASSOCIATED( FreeSurfGrad2 ) ) THEN
        FreeSurfGrad2Perm    =>  FreeSurfGrad2 % Perm
        FreeSurfGrad2Values  =>  FreeSurfGrad2 % Values
  ELSE
        CALL FATAL(SolverName,'No variable FreeSurfGrad2 found')  
  END IF

  Beta => VariableGet( Model % Variables, 'bottom Beta' )
  IF ( ASSOCIATED( Beta ) ) THEN
        BetaPerm    =>  Beta % Perm
        BetaValues  =>  Beta % Values
  ELSE
        CALL FATAL(SolverName,'No variable bottom Beta found')  
  END IF

  Height => VariableGet( Model % Variables, 'Height' )
  IF ( ASSOCIATED( Height ) ) THEN
        HeightPerm    =>  Height % Perm
        HeightValues  =>  Height % Values
  ELSE
        CALL FATAL(SolverName,'No variable Height found')  
  END IF

  Integral => VariableGet( Model % Variables, 'int IntegralTerm' )
  IF ( ASSOCIATED( Integral ) ) THEN
        IntegralPerm    => Integral % Perm
        IntegralValues  =>  Integral % Values
  ELSE
        CALL FATAL(SolverName,'No variable int IntegralTerm found')  
  END IF

  BasalMelting => VariableGet( Model % Variables, 'bottom BasalMelting' )
  IF ( ASSOCIATED( BasalMelting ) ) THEN
        BasalMeltingPerm    => BasalMelting % Perm
        BasalMeltingValues  =>  BasalMelting % Values
  ELSE
        CALL FATAL(SolverName,'No variable bottom BasalMelting found')  
  END IF


  !------------------------------------------------------------------------------
  ! Non-linear iteration loop
  !------------------------------------------------------------------------------
  DO iter=1,NonlinearIter
     !------------------------------------------------------------------------------
     !    assign matrices
     !------------------------------------------------------------------------------
     LocalNodes = Model % NumberOfNodes
     !Norm = Solver % Variable % Norm     
     WRITE(Message,'(a,I4,a,I4)') 'Non-linear Iteration ', iter,' out of max. ',NonlinearIter
     CALL INFO( SolverName, Message, Level=4)
     !------------------------------------------------------------------------------
     !    Do some additional initialization, and go for it
     !------------------------------------------------------------------------------
     totat = 0.0_dp
     totst = 0.0_dp
     at = CPUTime()
     CALL INFO( SolverName, 'start assembly', Level=4 )
     CALL DefaultInitialize()
     !------------------------------------------------------------------------------
     !    Do the assembly
     !------------------------------------------------------------------------------
     DO t=1,Solver % NumberOfActiveElements
        CurrentElement => GetActiveElement(t)
        n = GetElementNOFNodes()
        NodeIndexes => CurrentElement % NodeIndexes


        ! set coords of highest occuring dimension to zero (to get correct path element)
        !-------------------------------------------------------------------------------
        ElementNodes % x(1:n) = Solver % Mesh % Nodes % x(NodeIndexes)
        IF (DIM == 2) THEN
           ElementNodes % y(1:n) = 0.0
           ElementNodes % z(1:n) = 0.0
        ELSE IF (DIM == 3) THEN
           ElementNodes % y(1:n) = Solver % Mesh % Nodes % y(NodeIndexes)
           ElementNodes % z(1:n) = 0.0_dp
        ELSE
           WRITE(Message,'(a,i1,a)')&
                'It is not possible to compute free-surface problems in DIM=',&
                DIM, ' dimensions. Aborting'
           CALL Fatal( SolverName, Message) 
           STOP   
        END IF


        ! get pointers on Equation, Material and body-Force section input
        !----------------------------------------------------------------
        Equation => GetEquation()
        Material => GetMaterial()
        BodyForce => GetBodyForce()

        Density(1:n) = GetReal(Material, 'Density', Found)
        IF(.NOT. Found) THEN 
            CALL FATAL(SolverName, 'Density not found in Material')
        END IF

        Gravity(1:n) = GetReal(BodyForce, 'Flow Bodyforce 3', Found)
        IF(.NOT. Found) THEN 
            CALL FATAL(SolverName, 'Flow Body Force 3 not found')
        END IF

        Gravity(1:n) = Gravity(1:n)*(-1.0_dp)

        DO i=1,n
            Thickness(i) = HeightValues(HeightPerm(NodeIndexes(i)))
        END DO

        DO i=1,n
           Gradh(i) =  ( FreeSurfGrad1Values(FreeSurfGrad1Perm(NodeIndexes(i)))**2.0_dp &
                       + FreeSurfGrad2Values(FreeSurfGrad2Perm(NodeIndexes(i)))**2.0_dp )**0.5_dp
        END DO

        DO i=1,n
           IntegrationPart(i) = IntegralValues(IntegralPerm(NodeIndexes(i)))
        END DO

        scaling = GetLogical(Material, 'Free surface basal sliding scaling', Found)
        IF (.NOT. Found) THEN
           CALL FATAL(SolverName, 'Free surface basal sliding scaling not found')
        END IF
        IF(scaling) THEN
            scaleFactor = GetConstReal(Material, 'Free surface basal scaling factor', Found )       
            IF (.NOT. Found) THEN
                CALL FATAL(SolverName, 'Free surface basal scaling factor not found')
            END IF
        END IF

        DO i=1,n
           Drag(i) =  10.0_dp**BetaValues(BetaPerm(NodeIndexes(i)))
           IF (scaling) THEN
              Drag(i) = Drag(i) / scaleFactor
           END IF
        END DO

        ! Basal melting
        !---------------
        DO i=1,n
           Melting(i) =  BasalMeltingValues(BasalMeltingPerm(NodeIndexes(i)))
        END DO

        ! get lower limit for solution 
        !-----------------------------
        LowerLimit(CurrentElement % Nodeindexes(1:N)) = &
             ListGetReal(Material,'Min ' // TRIM(VariableName),n,CurrentElement % NodeIndexes, Found) 
        LimitedSolution(CurrentElement % Nodeindexes(1:N), 1) = Found
        ! get upper limit for solution 
        !-----------------------------
        UpperLimit(CurrentElement % Nodeindexes(1:N)) = &
             ListGetReal(Material,'Max ' // TRIM(VariableName),n,CurrentElement % NodeIndexes, Found)              
        LimitedSolution(CurrentElement % Nodeindexes(1:N), 2) = Found

        !------------------------------------------------------------------------------
        ! Get mesh velocity
        !------------------------------------------------------------------------------
        MeshVelocity = 0.0_dp
        CALL GetVectorLocalSolution( MeshVelocity, 'Mesh Velocity',CurrentElement)
        !------------------------------------------------------------------------------
        !      get the accumulation/ablation rate (i.e. normal surface flux)
        !      from the body force section
        !------------------------------------------------------------------------------
        SourceFunc = 0.0_dp
        Flux  = 0.0_dp
        SubstantialSurface = .TRUE.

        IF (ASSOCIATED( BodyForce ) ) THEN
           SubstantialSurface = .FALSE.
           ! Accumulation/ablation is given in normal direction of surface:
           !---------------------------------------------------------------
           SourceFunc(1:n) = GetReal( BodyForce, &
                TRIM(VariableName) // ' Accumulation', NormalFlux ) 
           ! Accumulation/ablation has to be computed from given flux:
           !----------------------------------------------------------
           IF (.NOT.NormalFlux) THEN
              Flux(1,1:n) = GetReal( BodyForce, TRIM(VariableName) // ' Accumulation Flux 1',Found)
              IF (.NOT.Found) Flux(1,1:n) = 0.0_dp
              IF (DIM >= 2) THEN
                 Flux(2,1:n) = GetReal( BodyForce, TRIM(VariableName) // ' Accumulation Flux 2',Found )
                 IF (.NOT.Found) Flux(2,1:n) = 0.0_dp
              ELSE
                 Flux(2,1:n) = 0.0_dp
              END IF
              IF (DIM == 3) THEN
                 Flux(3,1:n) = GetReal( BodyForce, TRIM(VariableName) // ' Accumulation Flux 3',Found )
                 IF (.NOT.Found) Flux(3,1:n) = 0.0_dp
              ELSE
                 Flux(3,1:n) = 0.0_dp
              END IF
              SourceFunc = 0.0_dp
           END IF
        END IF

        IF( TransientSimulation) THEN
           ElemFreeSurf(1:n) = PreFreeSurf(FreeSurfPerm(NodeIndexes),1)
        END IF

        !------------------------------------------------------------------------------
        !      Get element local matrix, and rhs vector
        !------------------------------------------------------------------------------
        CALL LocalMatrix( STIFF, MASS, FORCE,&
             SourceFunc, ElemFreeSurf, MeshVelocity, CurrentElement,&
             n, ElementNodes, NodeIndexes, TransientSimulation,&
             Flux, Density, Gravity, Thickness, Gradh, IntegrationPart, Drag, Melting, &
             NormalFlux, SubstantialSurface, ALEFormulation)

        !------------------------------------------------------------------------------
        !      If time dependent simulation add mass matrix to stiff matrix
        !------------------------------------------------------------------------------
        TimeForce = 0.0_dp
        IF ( TransientSimulation ) THEN
           !------------------------------------------------------------------------------
           !        NOTE: This will replace STIFF and LocalForce with the
           !              combined information...
           !------------------------------------------------------------------------------
           CALL Default1stOrderTime( MASS, STIFF, FORCE )
        END IF

        !------------------------------------------------------------------------------
        !      Update global matrix and rhs vector from local matrix & vector
        !------------------------------------------------------------------------------
        CALL DefaultUpdateEquations( STIFF, FORCE )
        !------------------------------------------------------------------------------
     END DO ! End loop bulk elements
     CALL DefaultFinishBulkAssembly()

     !------------------------------------------------------------------------------
     !     Neumann & Newton boundary conditions
     !------------------------------------------------------------------------------
     !
     ! MIND: In weak formulation it is not possible to prescribe a contact angle on
     !       a boundary in this solver. This has to be taken care of in the boundary
     !       condition for the stress tensor in the Navier-Stokes Solver. Thus, in
     !       generally it does not make sense to prescribe a Neumann type of
     !       condition here.

     !------------------------------------------------------------------------------
     !    FinishAssemebly must be called after all other assembly steps, but before
     !    Dirichlet boundary settings. Actually no need to call it except for
     !    transient simulations.
     !------------------------------------------------------------------------------
     CALL DefaultFinishAssembly()
     CALL DefaultDirichletBCs()

     !------------------------------------------------------------------------------
     !    Manipulation of the assembled matrix due to limits
     !------------------------------------------------------------------------------
     OldValues = SystemMatrix % Values
     OldRHS = ForceVector

     IF (ApplyDirichlet) THEN
        ! manipulation of the matrix
        !---------------------------
        DO i=1,Model % Mesh % NumberOfNodes
           k = FreeSurfPerm(i)           
           IF ((ActiveNode(i,1) .OR. ActiveNode(i,2)) .AND. (k > 0)) THEN
              CALL ZeroRow( SystemMatrix, k ) 
              CALL SetMatrixElement( SystemMatrix, k, k, 1.0_dp ) 
              IF(ActiveNode(i,1)) THEN
                 SystemMatrix % RHS(k) = LowerLimit(i)
              ELSE
                 SystemMatrix % RHS(k) = UpperLimit(i)
              END IF
           END IF
        END DO
     END IF

        CALL INFO( SolverName, 'Assembly done', Level=4 )
        !------------------------------------------------------------------------------
        !    Solve System  and check for convergence
        !------------------------------------------------------------------------------
        at = CPUTime() - at
        st = CPUTime() 

        PrevNorm = Solver % Variable % Norm

        Norm = DefaultSolve()

        IF ( PrevNorm + Norm /= 0.0_dp ) THEN
           RelativeChange = 2.0_dp * ABS( PrevNorm-Norm ) / (PrevNorm + Norm)
        ELSE
           RelativeChange = 0.0_dp
        END IF

        WRITE( Message, * ) 'Result Norm   : ',Norm
        CALL INFO( SolverName, Message, Level=4 )
        WRITE( Message, * ) 'Relative Change : ',RelativeChange
        CALL INFO( SolverName, Message, Level=4 )

        !------------------------------------------------------------------------------
        ! compute residual
        !------------------------------------------------------------------------------ 
        SystemMatrix % Values = OldValues
        ForceVector = OldRHS

        IF ( ParEnv % PEs > 1 ) THEN !!!!!!!!!!!!!!!!!!!!!! we have a parallel run
           CALL ParallelInitSolve( SystemMatrix, FreeSurf, ForceVector, ResidualVector )
           CALL ParallelMatrixVector( SystemMatrix, FreeSurf, StiffVector, .TRUE. )
           ResidualVector =  StiffVector - ForceVector
           CALL ParallelSumVector( SystemMatrix, ResidualVector )
        ELSE !!!!!!!!!!!!!!!!!!!!!! serial run 
           CALL CRS_MatrixVectorMultiply( SystemMatrix, FreeSurf, StiffVector)
           ResidualVector =  StiffVector - ForceVector
        END IF
        !-----------------------------
        ! determine "active" nodes set
        !-----------------------------
        IF (ApplyDirichlet) THEN
           numberofsurfacenodes = 0
           DO i=1,Model % NumberOfNodes
              l= FreeSurfPerm(i)  
              IF (l<1) CYCLE
              numberofsurfacenodes = numberofsurfacenodes + 1
              !---------------------------------------------------------
              ! if upper limit is exceeded, manipulate matrix in any case
              !----------------------------------------------------------
              IF ((LimitedSolution(i,1)).AND.(FreeSurf(l)-LowerLimit(i)<0.0_dp )) THEN
                 ActiveNode(i,1) = .TRUE.
              END IF
              IF ((LimitedSolution(i,2)).AND.(FreeSurf(l)-UpperLimit(i)>0.0_dp )) THEN
                 ActiveNode(i,2) = .TRUE.
              END IF

              IF ( LimitedSolution(i,1) .AND. ResidualVector(l) < -LinearTol & 
                       .AND. iter>1 ) ActiveNode(i,1) = .FALSE.
              IF ( LimitedSolution(i,2) .AND. ResidualVector(l) >  LinearTol & 
                       .AND. iter>1 ) ActiveNode(i,2) = .FALSE.

              IF( .NOT.ActiveNode(i,1) .AND. .NOT.ActiveNode(i,2) ) THEN
                 PointerToResidualVector(VarSurfResidual % Perm(i)) = 0.0_dp
              ELSE
                 PointerToResidualVector(VarSurfResidual % Perm(i)) = ResidualVector(l)
              END IF
            END DO
        END IF
        !------------------------------------------
        ! special treatment for periodic boundaries
        !------------------------------------------

        !------------------------------------------------------------------------------
        ! Relaxation
        !------------------------------------------------------------------------------
        IF(NeedOldValues) THEN
           IF(LimitDisp) THEN 
              maxdh = -HUGE(maxdh)         
              DO i=1, Model % NumberOfNodes
                 j = FreeSurfPerm(i)
                 IF(j > 0) THEN
                    maxdh = MAX(maxdh, ABS(FreeSurf(j)-OldFreeSurf(j)))
                 END IF
              END DO
              IF(maxdh > MaxDisp) THEN
                 Relax = Relax * MaxDisp/maxdh
              END IF
              WRITE(Message,'(a,E8.2)') 'Maximum displacement ',maxdh
              CALL INFO( SolverName, Message, Level=4 )
           END IF
           WRITE(Message,'(a,F8.2)') 'pp Relaxation factor',Relax
           CALL INFO( SolverName, Message, Level=4 )
           DO i=1, Model % NumberOfNodes
              j = FreeSurfPerm(i)
              IF(j > 0) THEN
                 FreeSurf(j) = Relax * FreeSurf(j) + (1-Relax) * OldFreeSurf(j)
              END IF
           END DO
        END IF

        st = CPUTIme()-st
        totat = totat + at
        totst = totst + st

        WRITE(Message,'(a,F8.2,F8.2)') 'Assembly: (s)', at, totat
        CALL INFO( SolverName, Message, Level=4 )
        WRITE(Message,'(a,F8.2,F8.2)') ' Solve:    (s)', st, totst
        CALL INFO( SolverName, Message, Level=4 )
        !------------------------------------------------------------------------------
        ! write some info on max/min values
        !------------------------------------------------------------------------------
        WRITE(Message,'(a,e13.6,a,e13.6)') &
             'Max/min values surface:', MAXVAL(FreeSurf(:)),'/',MINVAL( FreeSurf(:))
        CALL INFO(SolverName,Message,Level=4)
        IF (ApplyDirichlet) THEN
           !           WRITE(Message,'(a,i10)') 'Deactivated Periodic BC nodes:', k
           !          CALL INFO(SolverName,Message,Level=1)
           WRITE(Message,'(a,i10)') 'Number of surface nodes: ', numberofsurfacenodes
           CALL INFO(SolverName,Message,Level=1)
           WRITE(Message,'(a,i10)') 'Number of constrained points (lower limit):', COUNT(ActiveNode(:,1))
           CALL INFO(SolverName,Message,Level=1)
           WRITE(Message,'(a,i10)') 'Number of constrained points (upper limit):', COUNT(ActiveNode(:,2))
           CALL INFO(SolverName,Message,Level=1)
        END IF
        !----------------------
        ! check for convergence
        !----------------------
        IF ( RelativeChange < NonlinearTol ) THEN
           WRITE(Message,'(a,i10,a)') 'Converged after', iter, ' iterations'
           CALL INFO(SolverName,Message,Level=1)
           EXIT
        ELSE

        END IF
     END DO ! End loop non-linear iterations
     !------------------------------------------------------------------------------
   CONTAINS

     !------------------------------------------------------------------------------
     !==============================================================================
     SUBROUTINE LocalMatrix( STIFF, MASS, FORCE,&
          SourceFunc, OldFreeSurf, MeshVelo, &
          Element, n, Nodes, NodeIndexes, TransientSimulation,&
          Flux, Density, Gravity, Thickness, Gradh, IntegrationPart, Drag, Melting, & 
          NormalFlux, SubstantialSurface, ALEFormulation)
       !------------------------------------------------------------------------------
       !    INPUT:  SourceFunc(:)   nodal values of the accumulation/ablation function
       !            
       !            Element         current element
       !            n               number of nodes
       !            Nodes           current node points
       !
       !    OUTPUT: STIFF(:,:)
       !            MASS(:,:)
       !            FORCE(:)
       !------------------------------------------------------------------------------
       !      external variables:
       !      ------------------------------------------------------------------------
       REAL(KIND=dp) ::&
            STIFF(:,:), MASS(:,:), FORCE(:), SourceFunc(:), &
            MeshVelo(:,:), OldFreeSurf(:), Flux(:,:), Density(:), &
            Gravity(:), Thickness(:), Gradh(:), IntegrationPart(:), &
            Drag(:), Melting(:)

       INTEGER :: NodeIndexes(:)
       TYPE(Nodes_t) :: Nodes
       TYPE(Element_t), POINTER :: Element
       LOGICAL :: TransientSimulation,NormalFlux,SubstantialSurface,ALEFormulation
       !------------------------------------------------------------------------------
       !      internal variables:
       !      ------------------------------------------------------------------------
       REAL(KIND=dp) :: &
            Basis(2*n),dBasisdx(2*n,3), &
            Source, gradFreeSurf(3), normGradFreeSurf,&
            FluxGauss(3),X,Y,Z,U,V,W,S,detJ, hk, Load
       REAL(KIND=dp) :: LocalDensity, LocalGravity, LocalThickness, LocalGradh, &
                        LocalIntegration, LocalDrag, LocalMelting
       TYPE(ElementType_t), POINTER :: SaveElementType
       INTEGER :: LinType(2:4) = [202,303,404]

       LOGICAL :: Stat, UseLinear
       INTEGER :: i,j,t,p,q, n
       TYPE(GaussIntegrationPoints_t) :: IntegStuff
       !------------------------------------------------------------------------------

       FORCE = 0.0_dp
       STIFF = 0.0_dp
       MASS  = 0.0_dp

       UseLinear = GetLogical( GetSolverParams(), 'Use linear elements', Stat )
       UseLinear = UseLinear .OR. ANY(ActiveNode(NodeIndexes,:))
       UseLinear = UseLinear .AND. Element % TYPE % BasisFunctionDegree==2

       IF ( UseLinear ) THEN
         SaveElementType => Element % TYPE
         Element % TYPE => GetElementType(LinType(GetElementFamily()))
       END IF

       hK = ElementDiameter( Element, Nodes )

       !
       !      Numerical integration:
       !      ----------------------
       IntegStuff = GaussPoints( Element )

       DO t = 1,IntegStuff % n
          U = IntegStuff % u(t)
          V = IntegStuff % v(t)
          W = IntegStuff % w(t)
          S = IntegStuff % s(t)
          !
          !        Basis function values & derivatives at the integration point:
          !        -------------------------------------------------------------
          stat = ElementInfo( Element,Nodes, U, V, W, detJ, &
                              Basis, dBasisdx)

          !        Correction from metric
          !        ----------------------
          S = S * detJ

          IF ( CurrentCoordinateSystem() /= Cartesian ) THEN
             X = SUM( Nodes % x(1:n) * Basis(1:n) )
             Y = SUM( Nodes % y(1:n) * Basis(1:n) )
             Z = SUM( Nodes % z(1:n) * Basis(1:n) )
             S = S * X
          END IF
          !
          !        Velocities and (norm of) gradient of free surface and source function 
          !        at Gauss point
          !        ---------------------------------------------------------------------

          gradFreeSurf=0.0_dp

          LocalDensity = SUM( Density(1:n) * Basis(1:n) )
          LocalGravity = SUM( Gravity(1:n) * Basis(1:n) )
          LocalThickness = SUM( Thickness(1:n) * Basis(1:n) )
          LocalGradh = SUM( Gradh(1:n) * Basis(1:n) )
          LocalIntegration = SUM( IntegrationPart(1:n) * Basis(1:n) )
          LocalDrag = SUM( Drag(1:n) * Basis(1:n) )
          LocalMelting = SUM( Melting(1:n) * Basis(1:n) )

          DO i=1,DIM-1
             gradFreeSurf(i) = SUM(dBasisdx(1:n,i)*OldFreeSurf(1:n))
          END DO

          gradFreeSurf(DIM) = 1.0_dp
       
          IF (DIM==3) THEN
             normGradFreeSurf = SQRT(1.0_dp + gradFreeSurf(1)**2 + &
                  gradFreeSurf(2)**2)
          ELSE
             normGradFreeSurf = SQRT(1.0_dp + gradFreeSurf(1)**2)
          END IF

          DO p=1,n
             DO q=1,n
               ! The mass matrix, if needed
               !--------------------------- 
               IF (TransientSimulation) THEN
                  MASS(p,q) = MASS(p,q) + s * Basis(q)*Basis(p)
               END IF

               ! Stiffness matrix:
               !-----------------
               STIFF(p,q) =  STIFF(p,q) + s * &
                             (LocalDensity*LocalGravity*LocalThickness**2.0_dp/LocalDrag + &
                             2.0_dp*(LocalDensity*LocalGravity)**3.0_dp*LocalGradh**(3.0_dp-1.0_dp)*LocalIntegration) * &
                             SUM(dBasisdx(q,1:DIM-1)*dBasisdx(p,1:DIM-1)) 
             END DO
          END DO

          !        Get accumulation/ablation function if flux input is given
          !        (i.e., calculate vector product between flux and normal)
          !        --------------------------------------------------------- 
          IF (.NOT.(SubstantialSurface)) THEN
             IF (NormalFlux) THEN 
                Source = normGradFreeSurf * SUM( SourceFunc(1:n) * Basis(1:n) )
             ELSE
                DO i=1,dim
                   FluxGauss(i) = SUM(Basis(1:n)*Flux(i,1:n))
                END DO
                Source = SUM(FluxGauss(1:DIM)*gradFreeSurf(1:DIM))
             END IF
          ELSE
             Source = 0.0_dp
          END IF

          !        Assemble force vector:
          !        ---------------------
          DO p=1,n
             Load = (Source-LocalMelting)*Basis(p)
             FORCE(p) = FORCE(p) + s * Load
          END DO

       END DO

       IF (UseLinear) THEN
         EdgeMap => GetEdgeMap(GetElementFamily())
         n = ELement % TYPE % NumberOfNodes
         DO i=n+1,n+SIZE(EdgeMap,1)
           j=EdgeMap(i-n,1)
           k=EdgeMap(i-n,2)
           STIFF(i,:) =  0._dp
           STIFF(:,i) =  0._dp
           MASS(i,:)  =  0._dp
           MASS(:,i)  =  0._dp
           STIFF(i,i) =  1._dp
           STIFF(i,j) = -0.5_dp
           STIFF(i,k) = -0.5_dp
           FORCE(i) = 0._dp
           Element % TYPE => SaveElementType
         END DO
       END IF

       !------------------------------------------------------------------------------
     END SUBROUTINE LocalMatrix

     !------------------------------------------------------------------------------
   END SUBROUTINE FreeSurfaceSolverSIA
!------------------------------------------------------------------------------
