!/******************************************************************************
! * 
! *       ELMER, A Computational Fluid Dynamics Program.
! *
! *       Copyright 1st April 1995 - , Center for Scientific Computing,
! *                                    Finland.
! *
! *       All rights reserved. No part of this program may be used,
! *       reproduced or transmitted in any form or by any means
! *       without the written permission of CSC.
! *
! *****************************************************************************/
!
!/******************************************************************************
! *
! * Module containing a solver for fabric parameter equations  
! *
! ******************************************************************************
! *
! *                    Author:       Juha Ruokolainen, Hakime Seddik
! *
! *                    Address: Center for Scientific Computing
! *                            Tietotie 6, P.O. Box 405
! *                              02101 Espoo, Finland
! *                              Tel. +358 0 457 2723
! *                            Telefax: +358 0 457 2302
! *                          EMail: Juha.Ruokolainen@csc.fi
! *
! *                       Date: 08 Jun 1997
! *
! *                Modified by: Fabien Gillet
! *
! *       Date of modification: apr 2004
! *
! *****************************************************************************/

!------------------------------------------------------------------------------
   RECURSIVE SUBROUTINE GradTransFSolver( Model,Solver,dt,TransientSimulation )
!------------------------------------------------------------------------------

      USE DefUtils

    IMPLICIT NONE
!------------------------------------------------------------------------------
!******************************************************************************
!
!  Solve equation for the evolution of the tensor gradient of
!  transformation : DF/Dt=L.F  L=velocity gradient tensor.
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
     TYPE(Solver_t), POINTER :: PSolver

     TYPE(Matrix_t),POINTER :: StiffMatrix

     INTEGER :: dim,n1,n2,i,j,k,l,n,t,iter,NDeg,STDOFs,LocalNodes,istat

     TYPE(ValueList_t),POINTER :: Material, BC
     TYPE(Nodes_t) :: ElementNodes
     TYPE(Element_t),POINTER :: CurrentElement, Element, &
              ParentElement, LeftParent, RightParent, Edge

     REAL(KIND=dp) :: RelativeChange,UNorm,PrevUNorm,Gravity(3), &
         Tdiff,Normal(3),NewtonTol,NonlinearTol,s


     INTEGER :: NewtonIter,NonlinearIter

     TYPE(Variable_t), POINTER :: FSol,FVariable,FlowVariable,MeshVeloVariable

     REAL(KIND=dp), POINTER ::  FValues(:), FlowValues(:), Solution(:), Ref(:),MeshVeloValues(:)

     INTEGER, POINTER :: FPerm(:),NodeIndexes(:),FlowPerm(:), MeshVeloPerm(:)


     LOGICAL :: GotForceBC,GotIt,NewtonLinearization = .FALSE.

     INTEGER :: body_id,bf_id,eq_id, comp, Indexes(128)
!
     INTEGER :: old_body = -1

     INTEGER :: NumberOfEdges

                        
     LOGICAL :: AllocationsDone = .FALSE., FreeSurface

     TYPE(Variable_t), POINTER :: TimeVar

     REAL(KIND=dp),ALLOCATABLE:: MASS(:,:),STIFF(:,:),&
       LOAD(:,:),Force(:), LocalTemperature(:), Alpha(:,:),Beta(:), &
       Velocity(:,:),MeshVelocity(:,:), F11(:),F22(:),F12(:),F21(:),F33(:),F13(:),F31(:),F23(:),F32(:)

     SAVE MASS,STIFF,LOAD,F11,F22,F12,F21,F33,F13,F31,F23,F32, &
       Force,ElementNodes,Alpha,Beta, LocalTemperature,  AllocationsDone, &
       Velocity,MeshVelocity, old_body,dim
! *
!------------------------------------------------------------------------------
     REAL(KIND=dp) :: Bu,Bv,Bw,RM(3,3), SaveTime = -1
     REAL(KIND=dp), POINTER :: PrevF(:),CurrF(:),PostF(:)

     SAVE  PrevF, CurrF,PostF
     REAL(KIND=dp) :: at, at0, CPUTime, RealTime

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!    Get variables needed for solution
!------------------------------------------------------------------------------
     IF ( .NOT. ASSOCIATED( Solver % Matrix ) ) RETURN

     Solution => Solver % Variable % Values
     STDOFs   =  Solver % Variable % DOFs

     FSol    => VariableGet( Solver % Mesh % Variables, 'GradTransF' )
     FPerm   => FSol % Perm
     FValues => FSol % Values


     FlowVariable => VariableGet( Solver % Mesh % Variables, 'flow solution' )
     IF ( ASSOCIATED( FlowVariable ) ) THEN
       FlowPerm    => FlowVariable % Perm
       FlowValues  => FlowVariable % Values
     END IF

     !!!!! Mesh Velo
     MeshVeloVariable => VariableGet( Solver % Mesh % Variables, &
            'Mesh Velocity' )

     IF ( ASSOCIATED( MeshVeloVariable ) ) THEN
       MeshVeloPerm    => MeshVeloVariable % Perm
       MeshVeloValues  => MeshVeloVariable % Values
     END IF

                                  
! *
     StiffMatrix => Solver % Matrix
     UNorm = Solver % Variable % Norm
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!     Allocate some permanent storage, this is done first time only
!------------------------------------------------------------------------------
     IF ( .NOT. AllocationsDone .OR. Solver % Mesh % Changed) THEN
       N = Model % MaxElementNodes
! * Fab
       dim = CoordinateSystemDimension()
! *
       IF ( AllocationsDone ) THEN
         DEALLOCATE( LocalTemperature,     &
                     Force, Velocity, MeshVelocity, &
                     MASS, STIFF,      &
                     LOAD, Alpha, Beta )
       END IF

       ALLOCATE( LocalTemperature( N ), &
                 Force( 2*STDOFs*N ), &
                 Velocity(4, N ), MeshVelocity(3,N), &
                 MASS( 2*STDOFs*N,2*STDOFs*N ),  &
                 STIFF( 2*STDOFs*N,2*STDOFs*N ),  &
                 LOAD( 4,N ), Alpha( 3,N ), Beta( N ),STAT=istat )


       IF ( istat /= 0 ) THEN
          CALL Fatal( 'GradTransFSolve', 'Memory allocation error.' )
       END IF

       ALLOCATE( CurrF( 9*SIZE(Solver % Variable % Values)/dim ))
       CurrF = 0
! Fab

       ALLOCATE( PostF( 9*SIZE(Solver % Variable % Values)/dim ))
       PostF = 0

! Fab
       IF ( TransientSimulation ) THEN
          ALLOCATE( PrevF( 9*SIZE(Solver % Variable % Values)/dim ))
         PrevF = 0
       END IF

       DO i=1,Solver % NumberOFActiveElements
          CurrentElement => GetActiveElement(i)   
          NodeIndexes => CurrentElement % NodeIndexes
          n = GetElementDOFs( Indexes )
          Indexes(1:n) = Solver % Variable % Perm( Indexes(1:n) )
          DO COMP=1,9
            CurrF(9*(Indexes(1:n)-1)+COMP) = &
                        FValues(9*(FPerm(NodeIndexes)-1)+COMP)
          END DO
       END DO

       AllocationsDone = .TRUE.
     END IF

     IF( TransientSimulation ) THEN
        TimeVar => VariableGet( Solver % Mesh % Variables, 'Time' )
        IF ( SaveTime /= TimeVar % Values(1) ) THEN
           SaveTime = TimeVar % Values(1)
           PrevF = CurrF
        END IF
     END IF

!------------------------------------------------------------------------------
!    Do some additional initialization, and go for it
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
     NonlinearTol = ListGetConstReal( Solver % Values, &
        'Nonlinear System Convergence Tolerance' )

     NewtonTol = ListGetConstReal( Solver % Values, &
        'Nonlinear System Newton After Tolerance' )

     NewtonIter = ListGetInteger( Solver % Values, &
        'Nonlinear System Newton After Iterations' )

     NonlinearIter = ListGetInteger( Solver % Values, &
         'Nonlinear System Max Iterations',GotIt )

     IF ( .NOT.GotIt ) NonlinearIter = 1
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------

     old_body = -1
!------------------------------------------------------------------------------
     DO iter=1,NonlinearIter
!------------------------------------------------------------------------------
       at  = CPUTime()
       at0 = RealTime()


       CALL Info( 'GradTransFSolve', ' ', Level=4 )
       CALL Info( 'GradTransFSolve', ' ', Level=4 )
       CALL Info( 'GradTransFSolve', '--------------------------------',Level=4 )
       WRITE( Message, * ) 'GradTransF solver  iteration', iter
       CALL Info( 'GradTransFSolve', Message,Level=4 )
       CALL Info( 'GradTransFSolve', '--------------------------------',Level=4 )
       CALL Info( 'GradTransFSolve', ' ', Level=4 )
       CALL Info( 'GradTransFSolve', 'Starting assembly...',Level=4 )

!------------------------------------------------------------------------------


       PrevUNorm = UNorm
! * Fab
       DO COMP=1,dim
! *
       Do i=1,size(Solver % Variable % Values),dim
        Do j=1,dim
         Solver % Variable % Values(i+j-1) = CurrF( 9*(i-1)/dim+(COMP-1)*dim+j )
        End do
       End do

       IF ( TransientSimulation ) THEN
       Do i=1,size(Solver % Variable % Values),dim
        Do j=1,dim
          Solver % Variable % PrevValues(i+j-1,1) = PrevF( 9*(i-1)/dim+(COMP-1)*dim+j )
        End do
       End do
       END IF

       CALL DefaultInitialize()
!------------------------------------------------------------------------------
       DO t=1,Solver % NumberOFActiveElements
!------------------------------------------------------------------------------

         IF ( RealTime() - at0 > 1.0 ) THEN
           WRITE(Message,'(a,i3,a)' ) '   Assembly: ', INT(100.0 - 100.0 * &
            (Solver % NumberOfActiveElements-t) / &
               (1.0*Solver % NumberOfActiveElements)), ' % done'
                       
           CALL Info( 'GradTransFSolve', Message, Level=5 )
           at0 = RealTime()
         END IF

         CurrentElement => GetActiveElement(t)
         CALL GetElementNodes( ElementNodes )
         n = GetElementNOFNodes()
         n = GetElementDOFs( Indexes )
         NodeIndexes => CurrentElement % NodeIndexes
        
         Material => GetMaterial()
         body_id = CurrentElement % BodyId
 
         k = FlowVariable % DOFs
         Velocity = 0.0d0
         DO i=1,k
            Velocity(i,1:n) = FlowValues(k*(FlowPerm(NodeIndexes)-1)+i)
         END DO
 
	 !-------------------mesh velo
         MeshVelocity=0._dp
         IF (ASSOCIATED(MeshVeloVariable)) Then
           k = MeshVeloVariable % DOFs
           DO i=1,k
              MeshVelocity(i,1:n) = MeshVeloValues(k*(MeshVeloPerm(NodeIndexes)-1)+i)
           END DO
         EndIF
!--------------------------


         CALL LocalMatrix( COMP, MASS, &
              STIFF, FORCE, LOAD, &
              Velocity,MeshVelocity, CurrentElement, n, ElementNodes)

!------------------------------------------------------------------------------
!        Update global matrices from local matrices 
!------------------------------------------------------------------------------
         IF ( TransientSimulation )  CALL Default1stOrderTime(MASS,STIFF,FORCE)
         CALL DefaultUpdateEquations( STIFF, FORCE )
!------------------------------------------------------------------------------
      END DO
      CALL Info( 'GradTransFSolve', 'Assembly done', Level=4 )
!------------------------------------------------------------------------------
   
!------------------------------------------------------------------------------
!     Assembly of the edge terms
!------------------------------------------------------------------------------
      
       If (dim.eq.3) then
         NumberOfEdges=Solver % Mesh % NumberOfFaces
       ELSE
         NumberOfEdges=Solver % Mesh % NumberOfEdges
      ENDIF
      
      DO t=1,NumberOfEdges

       If (dim.eq.3) then
         Edge => Solver % Mesh % Faces(t)
       ELSE
         Edge => Solver % Mesh % Edges(t)
      ENDIF
      
         IF ( .NOT. ActiveBoundaryElement(Edge) ) CYCLE
       
         LeftParent  => Edge % BoundaryInfo % Left
         RightParent => Edge % BoundaryInfo % Right
         IF ( ASSOCIATED(RightParent) .AND. ASSOCIATED(LeftParent) ) THEN
            n  = GetElementNOFNodes( Edge )
            n1 = GetElementNOFNodes( LeftParent )
            n2 = GetElementNOFNodes( RightParent )

            k = FlowVariable % DOFs
            Velocity = 0.0d0
            DO i=1,k
               Velocity(i,1:n) = FlowValues(k*(FlowPerm(Edge % NodeIndexes)-1)+i)
            END DO
	    
	    !-------------------mesh velo
         MeshVelocity=0._dp
      IF ( ASSOCIATED( MeshVeloVariable ) ) THEN
            k = MeshVeloVariable % DOFs
            DO i=1,k
               MeshVelocity(i,1:n) = MeshVeloValues(k*(MeshVeloPerm(Edge % NodeIndexes)-1)+i)
            END DO
      END IF
    !--------------------------

            FORCE = 0.0d0
            MASS  = 0.0d0
            CALL LocalJumps( STIFF,Edge,n,LeftParent,n1,RightParent,n2,Velocity,MeshVelocity  )
            IF ( TransientSimulation )  CALL Default1stOrderTime(MASS, STIFF, FORCE)
            CALL DefaultUpdateEquations( STIFF, FORCE, Edge )
         END IF
      END DO

!------------------------------------------------------------------------------
!     Loop over the boundary elements
!------------------------------------------------------------------------------
      DO t = 1, Solver % Mesh % NumberOfBoundaryElements
!------------------------------------------------------------------------------

         Element => GetBoundaryElement(t)
         IF( .NOT. ActiveBoundaryElement() )  CYCLE
         IF( GetElementFamily(Element) == 1 ) CYCLE

         ParentElement => Element % BoundaryInfo % Left
         IF ( .NOT. ASSOCIATED( ParentElement ) ) &
            ParentElement => Element % BoundaryInfo % Right
          
         n  = GetElementNOFNodes( Element )
         n1 = GetElementNOFnodes( ParentElement )
       
         k = FlowVariable % DOFs
         Velocity = 0.0d0
         DO i=1,k
            Velocity(i,1:n) = FlowValues(k*(FlowPerm(Element % NodeIndexes)-1)+i)
         END DO
	!-------------------mesh velo
         MeshVelocity=0._dp
      IF ( ASSOCIATED( MeshVeloVariable ) ) THEN
        k = MeshVeloVariable % DOFs
        DO i=1,k
         MeshVelocity(i,1:n) = MeshVeloValues(k*(MeshVeloPerm(Element % NodeIndexes)-1)+i)
        End do
      END IF
    !--------------------------
         BC => GetBC()
         LOAD = 0.0d0
         GotIt = .FALSE.
         IF ( ASSOCIATED(BC) ) THEN
           Do i=1,dim
            LOAD(i,1:n) = GetReal( BC, ComponentName('GradTransF', dim*(comp-1)+i ) , GotIt )
           End do
         END IF

         MASS = 0.0d0
         CALL LocalMatrixBoundary(  STIFF, FORCE, LOAD, &
                              Element, n, ParentElement, n1, Velocity, MeshVelocity, GotIt )

         IF ( TransientSimulation )  CALL Default1stOrderTime(MASS, STIFF, FORCE)
         CALL DefaultUpdateEquations( STIFF, FORCE )
      END DO

      CALL DefaultFinishAssembly()
!------------------------------------------------------------------------------
      CALL Info( 'GradTransFSolve', 'Set boundaries done', Level=4 )

!------------------------------------------------------------------------------
!     Solve the system and check for convergence
!------------------------------------------------------------------------------
      Unorm = DefaultSolve()
      !CurrF( COMP::9 ) = Solver % Variable % Values
print*,'solve done', minval( solver % variable % values ), maxval( Solver % variable % values)

      
      n1 = Solver % Mesh % NumberOfNodes
      ALLOCATE( Ref(n1) )
      Ref = 0
      
      Do i=1,dim
           FValues((comp-1)*dim+i::9 ) = 0.
      End do
      
      DO t=1,Solver % NumberOfActiveElements
         Element => GetActiveElement(t) 
         n = GetElementNOFNodes()
         n = GetElementDOFs( Indexes )
         
         DO i=1,n
            k = Element % NodeIndexes(i)
            Do j=1,STDOFs
            FValues( 9*(FPerm(k)-1) + (COMP-1)*dim+j ) =    & 
             FValues( 9*(FPerm(k)-1) + (COMP-1)*dim+j ) + &
              Solver % Variable % Values( STDOFs*(Solver % Variable % Perm(Indexes(i))-1)+j )
            End do
            Ref(k) = Ref(k) + 1
         END DO
      END DO

      DO i=1,n1
         k=FPerm(i)
         IF (k < 1) CYCLE
         Do j=1,STDOFs
          IF ( Ref(i) > 0 ) THEN
            FValues( 9*(FPerm(i)-1)+(COMP-1)*dim+j ) = &
                   FValues( 9*(FPerm(i)-1)+(COMP-1)*dim+j ) / Ref(i)
          END IF
         End do
      END DO

      DEALLOCATE( Ref )


!!!!*
      END DO ! End DO Comp


       DO i=1,Solver % NumberOFActiveElements
          CurrentElement => GetActiveElement(i)   
          NodeIndexes => CurrentElement % NodeIndexes
          n = GetElementDOFs( Indexes )
          Indexes(1:n) = Solver % Variable % Perm( Indexes(1:n) )
          DO COMP=1,9
            CurrF(9*(Indexes(1:n)-1)+COMP) = &
                        FValues(9*(FPerm(NodeIndexes)-1)+COMP)
          END DO
       END DO

! fab
      Unorm = SQRT( SUM( FValues**2 ) / SIZE(FValues) )
      Solver % Variable % Norm = Unorm  

      IF ( PrevUNorm + UNorm /= 0.0d0 ) THEN
         RelativeChange = 2.0d0 * ABS( PrevUNorm - UNorm) / ( PrevUnorm + UNorm)
      ELSE
         RelativeChange = 0.0d0
      END IF

      WRITE( Message, * ) 'Result Norm   : ',UNorm
      CALL Info( 'GradTransFSolve', Message, Level=4 )
      WRITE( Message, * ) 'Relative Change : ',RelativeChange
      CALL Info( 'GradTransFSolve', Message, Level=4 )

      
!------------------------------------------------------------------------------
      IF ( RelativeChange < NewtonTol .OR. &
            iter > NewtonIter ) NewtonLinearization = .TRUE.

      IF ( RelativeChange < NonLinearTol ) EXIT
!------------------------------------------------------------------------------
    END DO ! of nonlinear iter
!------------------------------------------------------------------------------

CONTAINS

!------------------------------------------------------------------------------
      SUBROUTINE LocalMatrix( Comp, MASS, STIFF, FORCE, LOAD, &
          NodalVelo,NodMeshVelo, Element, n, Nodes )
!------------------------------------------------------------------------------

     REAL(KIND=dp) :: STIFF(:,:),MASS(:,:)
     REAL(KIND=dp) :: LOAD(:,:), NodalVelo(:,:),NodMeshVelo(:,:)
     REAL(KIND=dp), DIMENSION(:) :: FORCE

     TYPE(Nodes_t) :: Nodes
     TYPE(Element_t) :: Element

     INTEGER :: n, Comp
!------------------------------------------------------------------------------
!
     REAL(KIND=dp) :: Basis(2*n),ddBasisddx(1,1,1)
     REAL(KIND=dp) :: dBasisdx(2*n,3),SqrtElementMetric


     REAL(KIND=dp) :: A(3,3),M(3,3), hK,tau,pe1,pe2,unorm,C0, SU(n), SW(n)
     REAL(KIND=dp) :: LoadAtIp, Temperature
! * Fab
    

     INTEGER :: i,j,k,p,q,t,dim,NBasis,ind(3), DOFs = 1

     REAL(KIND=dp) :: s,u,v,w, Radius, B(6,3), G(3,6), C44,C55,C66
     REAL(KIND=dp) :: Velo(3)
     REAL(KIND=dp) :: F11,F22,F12,F21,F33,F13,F31,F32,F23

     REAL(KIND=dp) :: LGrad(3,3),LGradp(3,3)
     REAL(KIND=dp) :: hmax,xmoy,socleB
     INTEGER :: N_Integ
     REAL(KIND=dp), DIMENSION(:), POINTER :: U_Integ,V_Integ,W_Integ,S_Integ

     LOGICAL :: stat, CSymmetry,Fond
     TYPE(GaussIntegrationPoints_t), TARGET :: IntegStuff

!------------------------------------------------------------------------------
     dim = CoordinateSystemDimension()

      Fond=.False.
      hmax=maxval(Nodes % y(1:n))
      xmoy=SUM(Nodes % x(1:n) * Basis(1:n))
      !IF (hmax.LT.(socleB(xmoy)+120.))  Fond=.True.
     

     FORCE = 0.0D0
     MASS  = 0.0D0
     STIFF = 0.0D0
!    
!    Integration stuff:
!    ------------------
     NBasis = n
     IntegStuff = GaussPoints( Element  )

     U_Integ => IntegStuff % u
     V_Integ => IntegStuff % v
     W_Integ => IntegStuff % w
     S_Integ => IntegStuff % s
     N_Integ =  IntegStuff % n

     hk = ElementDiameter( Element, Nodes )
!
!   Now we start integrating:
!   -------------------------
    DO t=1,N_Integ

      u = U_Integ(t)
      v = V_Integ(t)
      w = W_Integ(t)

!------------------------------------------------------------------------------
!     Basis function values & derivatives at the integration point
!------------------------------------------------------------------------------
      stat = ElementInfo( Element,Nodes,u,v,w,SqrtElementMetric, &
               Basis, dBasisdx, ddBasisddx, .FALSE. )

      s = SqrtElementMetric * S_Integ(t)
!------------------------------------------------------------------------------
!
!
!    Velocity Gradient
!    -----------------
      LGrad = MATMUL( NodalVelo(:,1:n), dBasisdx(1:n,:) )

!
!     Velocity :
!     ----------
      Velo = 0.0d0
      DO i=1,dim
	 Velo(i) = SUM( Basis(1:n) * (NodalVelo(i,1:n)-NodMeshVelo(i,1:n)) )
      END DO
      Unorm = SQRT( SUM( Velo**2 ) )



!     Loop over basis functions (of both unknowns and weights):
!     ---------------------------------------------------------
      DO p=1,NBasis
         DO q=1,NBasis
            A = 0.0d0
            M = 0.0d0

            Do i=1,dim
             M(i,i) = Basis(p) * Basis(q)
             Do j=1,dim
              A(i,i) = A(i,i) - Velo(j) * Basis(q) * dBasisdx(p,j)
              A(i,j) = A(i,j) - LGrad(i,j) * Basis(q) * Basis(p)
             End do
            End do

!           Add nodal matrix to element matrix:
!           -----------------------------------
           Do i=1,dim
            Do j=1,dim
             MASS(dim*(p-1)+i ,dim*(q-1)+j )  = MASS(dim*(p-1)+i,dim*(q-1)+j )  + s * M(i,j)
             STIFF( dim*(p-1)+i ,dim*(q-1)+j ) = STIFF( dim*(p-1)+i ,dim*(q-1)+j ) + s * A(i,j)
            End do
           End do
         END DO
      END DO

              
   END DO 

!------------------------------------------------------------------------------
 END SUBROUTINE LocalMatrix
!------------------------------------------------------------------------------



!------------------------------------------------------------------------------
    SUBROUTINE FindParentUVW( Edge, nEdge, Parent, nParent, U, V, W, Basis )
!------------------------------------------------------------------------------
      IMPLICIT NONE
      TYPE(Element_t), POINTER :: Edge, Parent
      INTEGER :: nEdge, nParent
      REAL( KIND=dp ) :: U, V, W, Basis(:)
!------------------------------------------------------------------------------
      INTEGER :: i, j,l
      REAL(KIND=dp) :: NodalParentU(nEdge),NodalParentV(nEdge),NodalParentW(nEdge)
!------------------------------------------------------------------------------
      DO i = 1,nEdge
        DO j = 1,nParent
          IF ( Edge % NodeIndexes(i) == Parent % NodeIndexes(j) ) THEN
            NodalParentU(i) = Parent % Type % NodeU(j)
            NodalParentV(i) = Parent % Type % NodeV(j)
            NodalParentW(i) = Parent % Type % NodeW(j)
            EXIT
          END IF
        END DO
      END DO
      U = SUM( Basis(1:nEdge) * NodalParentU(1:nEdge) )
      V = SUM( Basis(1:nEdge) * NodalParentV(1:nEdge) )
      W = SUM( Basis(1:nEdge) * NodalParentW(1:nEdge) )
!------------------------------------------------------------------------------      
    END SUBROUTINE FindParentUVW
!------------------------------------------------------------------------------      


!------------------------------------------------------------------------------
    SUBROUTINE LocalJumps( STIFF,Edge,n,LeftParent,n1,RightParent,n2,Velo,MeshVelo )
!------------------------------------------------------------------------------
      REAL(KIND=dp) :: STIFF(:,:), Velo(:,:),MeshVelo(:,:)
      INTEGER :: n,n1,n2
      TYPE(Element_t), POINTER :: Edge, LeftParent, RightParent
!------------------------------------------------------------------------------
      REAL(KIND=dp) :: EdgeBasis(n), EdgedBasisdx(n,3), EdgeddBasisddx(n,3,3)
      REAL(KIND=dp) :: LeftBasis(n1), LeftdBasisdx(n1,3), LeftddBasisddx(n1,3,3)
      REAL(KIND=dp) :: RightBasis(n2), RightdBasisdx(n2,3), RightddBasisddx(n2,3,3)
      REAL(KIND=dp) :: Jump(n1+n2), Average(n1+n2)
      REAL(KIND=dp) :: detJ, U, V, W, S, Udotn, xx, yy
      LOGICAL :: Stat
      INTEGER :: i, j, p, q, dim, t, nEdge, nParent
      TYPE(GaussIntegrationPoints_t) :: IntegStuff
      REAL(KIND=dp) :: hE, Normal(3), cu(3), LeftOut(3)
      Real(KIND=dp),dimension(3,3),parameter :: id2=reshape((/1.,0.,0.,0.,1.,0.,0.,0.,1./), (/3,3/) )

      TYPE(Nodes_t) :: EdgeNodes, LeftParentNodes, RightParentNodes
      SAVE EdgeNodes, LeftParentNodes, RightParentNodes
!------------------------------------------------------------------------------
      dim = CoordinateSystemDimension()
      STIFF = 0.0d0

      CALL GetElementNodes( EdgeNodes, Edge )
      CALL GetElementNodes( LeftParentNodes,  LeftParent )
      CALL GetElementNodes( RightParentNodes, RightParent )
!------------------------------------------------------------------------------
!     Numerical integration over the edge
!------------------------------------------------------------------------------
      IntegStuff = GaussPoints( Edge )

      LeftOut(1) = SUM( LeftParentNodes % x(1:n1) ) / n1
      LeftOut(2) = SUM( LeftParentNodes % y(1:n1) ) / n1
      LeftOut(3) = SUM( LeftParentNodes % z(1:n1) ) / n1
      LeftOut(1) = SUM( EdgeNodes % x(1:n) ) / n - LeftOut(1)
      LeftOut(2) = SUM( EdgeNodes % y(1:n) ) / n - LeftOut(2)
      LeftOut(3) = SUM( EdgeNodes % z(1:n) ) / n - LeftOut(3)
      
      DO t=1,IntegStuff % n
        U = IntegStuff % u(t)
        V = IntegStuff % v(t)
        W = IntegStuff % w(t)
        S = IntegStuff % s(t)

        ! Basis function values & derivatives at the integration point:
        !--------------------------------------------------------------
        stat = ElementInfo( Edge, EdgeNodes, U, V, W, detJ, &
             EdgeBasis, EdgedBasisdx, EdgeddBasisddx, .FALSE. )

        S = S * detJ

        Normal = NormalVector( Edge, EdgeNodes, U, V, .FALSE. )
        IF ( SUM( LeftOut*Normal ) < 0 ) Normal = -Normal

        ! Find basis functions for the parent elements:
        ! ---------------------------------------------
        CALL FindParentUVW( Edge,n,LeftParent,n1,U,V,W,EdgeBasis )
        stat = ElementInfo( LeftParent, LeftParentNodes, U, V, W, detJ, &
                LeftBasis, LeftdBasisdx, LeftddBasisddx, .FALSE. )

        CALL FindParentUVW( Edge,n,RightParent,n2,U,V,W,EdgeBasis )
        stat = ElementInfo( RightParent, RightParentNodes, U, V, W, detJ, &
              RightBasis, RightdBasisdx, RightddBasisddx, .FALSE. )

        ! Integrate jump terms:
        ! ---------------------
        Jump(1:n1) = LeftBasis(1:n1)
        Jump(n1+1:n1+n2) = -RightBasis(1:n2)

        Average(1:n1) = LeftBasis(1:n1) / 2
        Average(n1+1:n1+n2) = RightBasis(1:n2) / 2

        cu = 0.0d0
        DO i=1,dim
	  cu(i) = SUM( (Velo(i,1:n)-MeshVelo(i,1:n)) * EdgeBasis(1:n) )
        END DO
        Udotn = SUM( Normal * cu )

        DO p=1,n1+n2
          DO q=1,n1+n2
           Do i=1,dim
            Do j=1,dim
            STIFF(dim*(p-1)+i ,dim*(q-1)+j) = STIFF(dim*(p-1)+i ,dim*(q-1)+j) + s * Udotn * Average(q) * Jump(p)*id2(i,j)
            STIFF(dim*(p-1)+i ,dim*(q-1)+j) = STIFF(dim*(p-1)+i ,dim*(q-1)+j) + s * ABS(Udotn)/2 * Jump(q) * Jump(p)*id2(i,j)
           End do
          End do
          END DO
        END DO
     END DO
!------------------------------------------------------------------------------
   END SUBROUTINE LocalJumps
!------------------------------------------------------------------------------



!------------------------------------------------------------------------------
   SUBROUTINE LocalMatrixBoundary( STIFF, FORCE, LOAD, &
        Element, n, ParentElement, np, Velo, MeshVelo, InFlowBC )
!------------------------------------------------------------------------------
     REAL(KIND=dp) :: STIFF(:,:),  FORCE(:), LOAD(:,:), Velo(:,:),MeshVelo(:,:)
     INTEGER :: n, np
     LOGICAL :: InFlowBC
     TYPE(Element_t), POINTER :: Element, ParentElement
!------------------------------------------------------------------------------
     REAL(KIND=dp) :: Basis(n), dBasisdx(n,3), ddBasisddx(n,3,3)
     REAL(KIND=dp) :: ParentBasis(np), ParentdBasisdx(np,3), ParentddBasisddx(np,3,3)
     INTEGER :: i,j,p,q,t,dim

     REAL(KIND=dp) :: Normal(3), g, L(3), Udotn, UdotnA, cu(3), detJ,U,V,W,S
     LOGICAL :: Stat
     TYPE(GaussIntegrationPoints_t) :: IntegStuff
      Real(KIND=dp),dimension(3,3),parameter :: id2=reshape((/1.,0.,0.,0.,1.,0.,0.,0.,1./), (/3,3/) )

     TYPE(Nodes_t) :: Nodes, ParentNodes
     SAVE Nodes, ParentNodes
!------------------------------------------------------------------------------
     dim = CoordinateSystemDimension()
     FORCE = 0.0d0
     STIFF = 0.0d0

     CALL GetElementNodes( Nodes, Element )
     CALL GetElementNodes( ParentNodes, ParentElement )

     ! Numerical integration:
     !-----------------------
     IntegStuff = GaussPoints( Element )
!
! Compute the average velocity.dot.Normal        
!        
     UdotnA = 0.0   
     DO t=1,IntegStuff % n
       U = IntegStuff % u(t)
       V = IntegStuff % v(t)
       W = IntegStuff % w(t)
       S = IntegStuff % s(t)

       Normal = NormalVector( Element, Nodes, U, V, .TRUE. ) 

       ! Basis function values & derivatives at the integration point:
       ! -------------------------------------------------------------
       stat = ElementInfo( Element, Nodes, U, V, W, detJ, &
               Basis, dBasisdx, ddBasisddx, .FALSE. )
       S = S * detJ
       cu = 0.0d0
       DO i=1,dim
          
	  cu(i) = SUM( (Velo(i,1:n)-MeshVelo(i,1:n)) * Basis(1:n) )
       END DO
       UdotnA = UdotnA + s*SUM( Normal * cu )

     END DO

     DO t=1,IntegStuff % n
       U = IntegStuff % u(t)
       V = IntegStuff % v(t)
       W = IntegStuff % w(t)
       S = IntegStuff % s(t)

       Normal = NormalVector( Element, Nodes, U, V, .TRUE. ) 

       ! Basis function values & derivatives at the integration point:
       ! -------------------------------------------------------------
       stat = ElementInfo( Element, Nodes, U, V, W, detJ, &
               Basis, dBasisdx, ddBasisddx, .FALSE. )
       S = S * detJ

       CALL FindParentUVW( Element, n, ParentElement, np, U, V, W, Basis )
       stat = ElementInfo( ParentElement, ParentNodes, U, V, W, &
           detJ,  ParentBasis, ParentdBasisdx, ParentddBasisddx, .FALSE. )

      Do i=1,dim
       L(i) = SUM( LOAD(i,1:n) * Basis(1:n) )
      End do
       cu = 0.0d0
       DO i=1,dim
	  cu(i) = SUM( (Velo(i,1:n)-MeshVelo(i,1:n)) * Basis(1:n) )
       END DO
       Udotn = SUM( Normal * cu )

       DO p = 1,np
        IF (InFlowBC .And. (UdotnA < 0.) ) THEN
           Do i=1,dim
            FORCE(dim*(p-1)+i) = FORCE(dim*(p-1)+i) - s * Udotn*L(i)*ParentBasis(p)
           End do
         ELSE
           DO q=1,np
            Do i=1,dim
             Do j=1,dim
             STIFF(dim*(p-1)+i ,dim*(q-1)+j) = STIFF(dim*(p-1)+i ,dim*(q-1)+j) + s*Udotn*ParentBasis(q)*ParentBasis(p)*id2(i,j)
             End do
            End do
           END DO
         END IF
       END DO
     END DO
!------------------------------------------------------------------------------
   END SUBROUTINE LocalMatrixBoundary
!------------------------------------------------------------------------------



!------------------------------------------------------------------------------
 END SUBROUTINE GradTransFSolver
!------------------------------------------------------------------------------
