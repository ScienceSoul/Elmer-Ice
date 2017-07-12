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
! *                     Author:       Juha Ruokolainen
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
! *                Modified by: Hakime Seddik
! *
! *       Date of modification: 06/07
! *
! *****************************************************************************/

!------------------------------------------------------------------------------
      RECURSIVE SUBROUTINE CaffeSolveFabric( Model,Solver,dt,TransientSimulation )
!------------------------------------------------------------------------------

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

     TYPE(Matrix_t),POINTER :: StiffMatrix

     INTEGER :: dim,n1,n2,i,j,k,l,n,t,iter,STDOFs,LocalNodes,istat

     TYPE(ValueList_t),POINTER :: Material, BC, Equation
     TYPE(Nodes_t) :: ElementNodes
     TYPE(Element_t),POINTER :: CurrentElement, Element, &
              ParentElement, LeftParent, RightParent, Edge

     REAL(KIND=dp) :: RelativeChange,UNorm,PrevUNorm,  &
         Tdiff,Normal(3),NewtonTol,NonlinearTol,s


     INTEGER :: NewtonIter,NonlinearIter

     TYPE(Variable_t), POINTER :: FabricSol, FabricVariable, FlowVariable, &
                                  MeshVeloVariable

     REAL(KIND=dp), POINTER :: Fabric(:), &
           FabricValues(:), FlowValues(:), &
           MeshVeloValues(:), Solution(:), Ref(:)

     INTEGER, POINTER :: FabricPerm(:),NodeIndexes(:), &
                         FlowPerm(:), MeshVeloPerm(:)

     REAL(KIND=dp) :: iota
     REAL(KIND=dp) :: A1plusA2
     Real(KIND=dp), parameter :: Rad2deg=180._dp/Pi
     REAL(KIND=dp) :: a2(6)
     REAL(KIND=dp) :: ai(3), Angle(3)

     LOGICAL :: GotIt,NewtonLinearization = .FALSE.

     INTEGER :: body_id,bf_id,eq_id, comp, Indexes(128)
!
     INTEGER :: old_body = -1
                        
     LOGICAL :: AllocationsDone = .FALSE., FreeSurface

     TYPE(Variable_t), POINTER :: TimeVar

     REAL(KIND=dp), ALLOCATABLE:: MASS(:,:), STIFF(:,:),  &
       LocalFluidity(:), LOAD(:,:),Force(:), &
       K1(:), K2(:), E1(:), E2(:), E3(:), &
       Velocity(:,:), MeshVelocity(:,:)

     SAVE MASS, STIFF, LOAD, Force,ElementNodes,& 
          LocalFluidity,  AllocationsDone, K1, K2, &
          E1, E2, E3, iota, Velocity, MeshVelocity, old_body, dim, comp
!------------------------------------------------------------------------------------------

     REAL(KIND=dp) :: SaveTime = -1
     REAL(KIND=dp), POINTER :: PrevFabric(:),CurrFabric(:),TempFabVal(:)

     SAVE  PrevFabric, CurrFabric,TempFabVal
     REAL(KIND=dp) :: at, at0, CPUTime, RealTime
!------------------------------------------------------------------------------
                                                
!------------------------------------------------------------------------------
 CALL Info( 'CaffeSolveFabric', 'version 2.0', Level=4 )
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!    Get variables needed for solution
!------------------------------------------------------------------------------
      IF ( .NOT. ASSOCIATED( Solver % Matrix ) ) RETURN

      Solution => Solver % Variable % Values
      STDOFs   =  Solver % Variable % DOFs

      FabricSol    => VariableGet( Solver % Mesh % Variables, 'Fabric' )
      IF ( ASSOCIATED( FabricSol) ) THEN
        FabricPerm   => FabricSol % Perm
        FabricValues => FabricSol % Values
      ELSE
        CALL Fatal( 'CAFFE (Fabric solver)', 'No Fabric associated!.' ) 
      END IF

      FlowVariable => VariableGet( Solver % Mesh % Variables, 'flow solution' )
      IF ( ASSOCIATED( FlowVariable ) ) THEN
       FlowPerm    => FlowVariable % Perm
       FlowValues  => FlowVariable % Values
      ELSE
       CALL Fatal( 'CAFFE (Fabric solver)', 'No Flow Solution associated!.' )
      END IF
      
!!!!! Mesh Velo
     MeshVeloVariable => VariableGet( Solver % Mesh % Variables, &
            'Mesh Velocity' )

     IF ( ASSOCIATED( MeshVeloVariable ) ) THEN
       MeshVeloPerm    => MeshVeloVariable % Perm
       MeshVeloValues  => MeshVeloVariable % Values
     END IF
       
                                       
      StiffMatrix => Solver % Matrix
      ! UNorm = Solver % Variable % Norm
      Unorm = SQRT( SUM( FabricValues**2 ) / SIZE(FabricValues) )
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!     Allocate some permanent storage, this is done first time only
!------------------------------------------------------------------------------
      IF ( .NOT. AllocationsDone .OR. Solver % Mesh % Changed) THEN
        N = Model % MaxElementNodes

       dim = CoordinateSystemDimension()
       
       IF ( AllocationsDone ) THEN
         DEALLOCATE( K1,K2,E1,E2,E3, &      
                     Force, LocalFluidity, &
                     Velocity,MeshVelocity, &
                     MASS,STIFF,      &
                     LOAD, CurrFabric, TempFabVal )
       END IF

       ALLOCATE( LocalFluidity( N ), &
                 K1( N ), K2( N ), E1( N ), E2( N ), E3( N ), &
                 Force( 2*STDOFs*N ), &
                 Velocity(4, N ),MeshVelocity(3,N), &
                 MASS( 2*STDOFs*N,2*STDOFs*N ),  &
                 STIFF( 2*STDOFs*N,2*STDOFs*N ),  &
                 LOAD( 4,N ), &
                 CurrFabric( 5*SIZE(Solver % Variable % Values)), &
                 TempFabVal( 5*SIZE(Solver % Variable % Values)), &
                 STAT=istat )


       IF ( istat /= 0 ) THEN
          CALL Fatal( 'CAFFE (Fabric solver)', 'Memory allocation error.' )
       END IF
       CALL Info('CAFFE (Fabric solver)','Memory allocations done', Level=3)

       CurrFabric = 0.
       TempFabVal = 0.
       IF ( TransientSimulation ) THEN
          IF (AllocationsDone ) DEALLOCATE (PrevFabric)
          ALLOCATE( PrevFabric( 5*SIZE(Solver % Variable % Values)) )
         PrevFabric = 0.
       END IF

       DO i=1,Solver % NumberOFActiveElements
          CurrentElement => GetActiveElement(i)   
          n = GetElementDOFs( Indexes )
          n = GetElementNOFNodes()
          NodeIndexes => CurrentElement % NodeIndexes
          Indexes(1:n) = Solver % Variable % Perm( Indexes(1:n) )
          DO COMP=1,5
            IF ( TransientSimulation ) THEN
               PrevFabric(5*(Indexes(1:n)-1)+COMP) = &
                   FabricValues(5*(FabricPerm(NodeIndexes(1:n))-1)+COMP)
            END IF
               CurrFabric(5*(Indexes(1:n)-1)+COMP) = &
                   FabricValues(5*(FabricPerm(NodeIndexes(1:n))-1)+COMP)
          END DO
       END DO

       AllocationsDone = .TRUE.
      END IF

      IF( TransientSimulation ) THEN
        TimeVar => VariableGet( Solver % Mesh % Variables, 'Time' )
        IF ( SaveTime /= TimeVar % Values(1) ) THEN
           SaveTime = TimeVar % Values(1)
           PrevFabric = CurrFabric
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
	 
!------------------------------------------------------------------------------
      DO iter=1,NonlinearIter
!------------------------------------------------------------------------------
       at  = CPUTime()
       at0 = RealTime()

       CALL Info( 'CaffeSolveFabric', ' ', Level=4 )
       CALL Info( 'CaffeSolveFabric', ' ', Level=4 )
       CALL Info( 'CaffeSolveFabric', &
                    '-------------------------------------',Level=4 )
       WRITE( Message, * ) 'Caffe Fabric solver  iteration', iter
       CALL Info( 'CaffeSolveFabric', Message,Level=4 )
       CALL Info( 'CaffeSolveFabric', &
                     '-------------------------------------',Level=4 )
!------------------------------------------------------------------------------

       PrevUNorm = UNorm
       
!Loop over the number of fabric equations to solve
       DO COMP=1,2*dim-1
       
       Solver % Variable % Values = CurrFabric( COMP::5 )
       IF ( TransientSimulation ) THEN
          Solver % Variable % PrevValues(:,1) = PrevFabric( COMP::5 )
       END IF

       CALL DefaultInitialize()

       CALL Info( 'CaffeSolveFabric', ' ', Level=4 )
       CALL Info( 'CaffeSolveFabric', 'Starting assembly...',Level=4 )
!------------------------------------------------------------------------------
       DO t=1,Solver % NumberOFActiveElements
!------------------------------------------------------------------------------

         IF ( RealTime() - at0 > 1.0 ) THEN
           WRITE(Message,'(a,i3,a)' ) '   Assembly: ', INT(100.0 - 100.0 * &
            (Solver % NumberOfActiveElements-t) / &
               (1.0*Solver % NumberOfActiveElements)), ' % done'
                       
           CALL Info( 'CaffeSolveFabric', Message, Level=5 )
           at0 = RealTime()
         END IF

         CurrentElement => GetActiveElement(t)
         CALL GetElementNodes( ElementNodes )
         n = GetElementDOFs( Indexes )
         n = GetElementNOFNodes()
         NodeIndexes => CurrentElement % NodeIndexes

         
         Material => GetMaterial()
         IF (.NOT. ASSOCIATED(Material)) THEN
         WRITE(Message,'(A,I5,A)') &
         'No material for bulk element no. ',t,' found.'
         CALL FATAL('CAFFE (Fabric solver)',Message)
         END IF

         body_id = CurrentElement % BodyId
!------------------------------------------------------------------------------
!        Read in material constants from Material section
!------------------------------------------------------------------------------
		IF (body_id /= old_body) Then 
           old_body = body_id
           CALL GetMaterialDefs()

		 Equation => GetEquation()
		 IF (.NOT.ASSOCIATED(Equation)) THEN
		 WRITE (Message, *) 'No Equation  found in model'
		 CALL FATAL('CAFFE (Fabric solver)',Message)
		 END IF
		END IF
	  	  
!------------------------------------------------------------------------------
!        Get element local stiffness & mass matrices
!------------------------------------------------------------------------------
        
         K1(1:n) = CurrFabric( 5*(Solver % Variable % Perm(Indexes(1:n))-1)+1 )
         K2(1:n) = CurrFabric( 5*(Solver % Variable % Perm(Indexes(1:n))-1)+2 )
         E1(1:n) = CurrFabric( 5*(Solver % Variable % Perm(Indexes(1:n))-1)+3 )
         E2(1:n) = CurrFabric( 5*(Solver % Variable % Perm(Indexes(1:n))-1)+4 )
         E3(1:n) = CurrFabric( 5*(Solver % Variable % Perm(Indexes(1:n))-1)+5 )

         k = FlowVariable % DOFs
         Velocity = 0.0d0
         DO i=1,k-1
            Velocity(i,1:n) = FlowValues(k*(FlowPerm(NodeIndexes)-1)+i)
         END DO

!------------meshvelocity
         MeshVelocity=0._dp
         IF (ASSOCIATED(MeshVeloVariable)) Then
           k = MeshVeloVariable % DOFs
           DO i=1,k
              MeshVelocity(i,1:n) = MeshVeloValues(k*(MeshVeloPerm(NodeIndexes)-1)+i)
           END DO
         EndIF
!----------------------------------

         CALL LocalMatrix( COMP, MASS, STIFF, FORCE, LOAD, K1, K2, E1, &
           E2, E3, Velocity, MeshVelocity, CurrentElement, n, ElementNodes, iota)

!------------------------------------------------------------------------------
!        Update global matrices from local matrices 
!------------------------------------------------------------------------------
         IF ( TransientSimulation )  CALL Default1stOrderTime(MASS,STIFF,FORCE)
         CALL DefaultUpdateEquations( STIFF, FORCE )
!------------------------------------------------------------------------------
      END DO
      CALL Info( 'CaffeSolveFabric', 'Assembly done', Level=4 )
!------------------------------------------------------------------------------
   
!------------------------------------------------------------------------------
!     Assembly of the edge terms
!------------------------------------------------------------------------------
!
!!Used for the discountinus galerkin method!
!3D => Edges => Faces
      If (dim.eq.3) then 
      DO t=1,Solver % Mesh % NumberOfFaces
         Edge => Solver % Mesh % Faces(t)
         IF ( .NOT. ActiveBoundaryElement(Edge) ) CYCLE
       
         LeftParent  => Edge % BoundaryInfo % Left
         RightParent => Edge % BoundaryInfo % Right
         IF ( ASSOCIATED(RightParent) .AND. ASSOCIATED(LeftParent) ) THEN
            n  = GetElementNOFNodes( Edge )
            n1 = GetElementNOFNodes( LeftParent )
            n2 = GetElementNOFNodes( RightParent )

            k = FlowVariable % DOFs
            Velocity = 0.0d0
            DO i=1,k-1
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
            CALL LocalJumps( STIFF,Edge,n,LeftParent,n1,RightParent,n2,Velocity,MeshVelocity )
            IF ( TransientSimulation )  CALL Default1stOrderTime(MASS, STIFF, FORCE)
            CALL DefaultUpdateEquations( STIFF, FORCE, Edge )
         END IF
      END DO
!
!  2D
      Else

      DO t=1,Solver % Mesh % NumberOfEdges
         Edge => Solver % Mesh % Edges(t)
         IF ( .NOT. ActiveBoundaryElement(Edge) ) CYCLE
       
         LeftParent  => Edge % BoundaryInfo % Left
         RightParent => Edge % BoundaryInfo % Right
         IF ( ASSOCIATED(RightParent) .AND. ASSOCIATED(LeftParent) ) THEN
            n  = GetElementNOFNodes( Edge )
            n1 = GetElementNOFNodes( LeftParent )
            n2 = GetElementNOFNodes( RightParent )

            k = FlowVariable % DOFs
            Velocity = 0.0d0
            DO i=1,k-1
               Velocity(i,1:n) = FlowValues(k*(FlowPerm(Edge % NodeIndexes(1:n))-1)+i)
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
            CALL LocalJumps( STIFF,Edge,n,LeftParent,n1,RightParent,n2,Velocity,MeshVelocity )
            IF ( TransientSimulation )  CALL Default1stOrderTime(MASS, STIFF, FORCE)
            CALL DefaultUpdateEquations( STIFF, FORCE, Edge )
         END IF
       END DO

      END IF

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
         DO i=1,k-1
            Velocity(i,1:n) = FlowValues(k*(FlowPerm(Element % NodeIndexes(1:n))-1)+i)
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
            LOAD(1,1:n) = GetReal( BC, ComponentName('Fabric', Comp) , GotIt )
         END IF

         MASS = 0.0d0
         CALL LocalMatrixBoundary(  STIFF, FORCE, LOAD(1,1:n), &
                              Element, n, ParentElement, n1, Velocity,MeshVelocity, GotIt )

         IF ( TransientSimulation )  CALL Default1stOrderTime(MASS, STIFF, FORCE)
         CALL DefaultUpdateEquations( STIFF, FORCE )
      END DO

      CALL DefaultFinishAssembly()
!------------------------------------------------------------------------------
      CALL Info( 'CaffeSolveFabric', 'Set boundaries done', Level=4 )

!------------------------------------------------------------------------------
!     Solve the system and check for convergence
!------------------------------------------------------------------------------
      Unorm = DefaultSolve()
!      CurrFabric( COMP::5 ) = Solver % Variable % Values
      WRITE(Message,*) 'solve done', minval( solver % variable % values), maxval( Solver % variable % values)
      CALL Info( 'CaffeSolveFabric', Message, Level=4 )
      
      n1 = Solver % Mesh % NumberOfNodes
      ALLOCATE( Ref(n1) )
      Ref = 0
      !
      ! FabricValues( COMP::5 ) = 0 !fab sinon remet toutes les
      ! fabriques ?? 0 dans le cas ou on a 2 domaines dont l'un a la
      ! fabrique fixe
      TempFabVal(COMP::5 ) = 0. !fab
      
      DO t=1,Solver % NumberOfActiveElements
         Element => GetActiveElement(t) 
         n = GetElementDOFs( Indexes )
         n = GetElementNOFNodes()
         
         DO i=1,n
            k = Element % NodeIndexes(i)
            TempFabVal( 5*(FabricPerm(k)-1) + COMP ) =    & 
            TempFabVal( 5*(FabricPerm(k)-1) + COMP ) + &
            Solver % Variable % Values( Solver % Variable % Perm(Indexes(i)) )
            FabricValues( 5*(FabricPerm(k)-1) + COMP ) = &
                          TempFabVal(5*(FabricPerm(k)-1) + COMP ) 
            Ref(k) = Ref(k) + 1
         END DO
      END DO

      DO i=1,n1
         j=FabricPerm(i)
         IF (j < 1) CYCLE
         IF ( Ref(i) > 0 ) THEN
            FabricValues( 5*(j-1)+COMP ) = &
                   FabricValues( 5*(j-1)+COMP ) / Ref(i)
         END IF
      END DO

      DEALLOCATE( Ref )

      SELECT CASE( Comp ) 
      CASE(1)
      FabricValues( COMP:SIZE(FabricValues):5 ) = &
          MIN(MAX( FabricValues( COMP:SIZE(FabricValues):5 ) , 0._dp),1._dp)
          
      CASE(2)
       FabricValues( COMP:SIZE(FabricValues):5 ) = &
       MIN(MAX( FabricValues( COMP:SIZE(FabricValues):5 ) , 0._dp),1._dp)
       
       !DO i=1,SIZE(FabricValues),5 
       !  IF((FabricValues(i)+FabricValues(i+1)).GT.1._dp) THEN
       !      A1plusA2=FabricValues(i)+FabricValues(i+1)
       !      FabricValues(i)= &
       !        FabricValues(i)/A1plusA2
       !      FabricValues(i+1)= &
       !        FabricValues(i+1)/A1plusA2
       !   END IF
       ! END DO
        

      CASE(3:5)
       DO i=COMP,SIZE(FabricValues),5
         IF(FabricValues(i).GT.0._dp) THEN
               FabricValues(i) = MIN( FabricValues(i) , 0.5_dp)
         ELSE
               FabricValues(i) = MAX( FabricValues(i) , -0.5_dp)
         END IF
       END DO
      END SELECT

      END DO ! End DO Comp

       DO i=1,Solver % NumberOFActiveElements
          CurrentElement => GetActiveElement(i)   
          n = GetElementDOFs( Indexes )
          n = GetElementNOFNodes()
          NodeIndexes => CurrentElement % NodeIndexes
          Indexes(1:n) = Solver % Variable % Perm( Indexes(1:n) )
          DO COMP=1,5
            CurrFabric(5*(Indexes(1:n)-1)+COMP) = &
                        FabricValues(5*(FabricPerm(NodeIndexes(1:n))-1)+COMP)
          END DO
       END DO


      Unorm = SQRT( SUM( FabricValues**2 ) / SIZE(FabricValues) )
      Solver % Variable % Norm = Unorm  

      IF ( PrevUNorm + UNorm /= 0.0d0 ) THEN
         RelativeChange = 2.0d0 * ABS( PrevUNorm - UNorm) / ( PrevUnorm + UNorm)
      ELSE
         RelativeChange = 0.0d0
      END IF

      WRITE( Message, * ) 'Result Norm   : ',UNorm
      CALL Info( 'CaffeSolveFabric', Message, Level=4 )
      WRITE( Message, * ) 'Relative Change : ',RelativeChange
      CALL Info( 'CaffeSolveFabric', Message, Level=4 )

      
!------------------------------------------------------------------------------
      IF ( RelativeChange < NewtonTol .OR. &
            iter > NewtonIter ) NewtonLinearization = .TRUE.

      IF ( RelativeChange < NonLinearTol ) EXIT

!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
    END DO ! of nonlinear iter
!------------------------------------------------------------------------------
CONTAINS

!------------------------------------------------------------------------------  
   SUBROUTINE GetMaterialDefs()
!------------------------------------------------------------------------------

   iota = ListGetConstReal( Material,'Grain Rotation parameter', GotIt )
   IF (.NOT. GotIt) THEN
       WRITE(Message,'(A)') 'No grain rotation parameter (iota) found in Material '
       CALL Fatal('CAFFE (Fabric solver)', Message)
   ELSE 
       WRITE(Message,'(A,F10.4)') 'Iota parameter = ', iota
       CALL INFO('CAFFE (Fabric solver)', Message, Level = 9)
   END IF 
!------------------------------------------------------------------------------
   END SUBROUTINE GetMaterialDefs
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
   SUBROUTINE LocalMatrix( Comp, MASS, STIFF, FORCE, LOAD, &
          NodalK1, NodalK2, NodalEuler1, NodalEuler2, NodalEuler3, & 
          NodalVelo, NodMeshVel, &
          Element, n, Nodes, iota)
!------------------------------------------------------------------------------

     REAL(KIND=dp) :: STIFF(:,:),MASS(:,:)
     REAL(KIND=dp) :: LOAD(:,:), NodalVelo(:,:),NodMeshVel(:,:)
     REAL(KIND=dp), DIMENSION(:) :: FORCE, NodalK1, NodalK2, NodalEuler1, &
               NodalEuler2, NodalEuler3

     TYPE(Nodes_t) :: Nodes
     TYPE(Element_t) :: Element
     INTEGER :: n, Comp
!------------------------------------------------------------------------------
!
     REAL(KIND=dp) :: Basis(2*n),ddBasisddx(1,1,1)
     REAL(KIND=dp) :: dBasisdx(2*n,3),SqrtElementMetric

     REAL(KIND=dp) :: A1, A2, A3, E1, E2, E3, Theta

     REAL(KIND=dp) :: A, M, hK,unorm,C0
     REAL(KIND=dp) :: LoadAtIp
     REAL(KIND=dp) :: iota, Deq, ai(6),a4(9),hmax

     INTEGER :: i,j,k,p,q,t,dim,NBasis,ind(3), DOFs = 1

     REAL(KIND=dp) :: s,u,v,w,Radius
     REAL(KIND=dp) :: Velo(3),Spin(3)

     REAL(KIND=dp) :: LGrad(3,3),StrainRate(3,3), angle(3),epsi
     REAL(KIND=dp) :: ap(3),C(6,6),Spin1(3,3)
     Integer :: INDi(6),INDj(6)
     LOGICAL :: CSymmetry
     
     LOGICAL :: Fond
              
!              
     INTEGER :: N_Integ
     REAL(KIND=dp), DIMENSION(:), POINTER :: U_Integ,V_Integ,W_Integ,S_Integ

     LOGICAL :: stat
     TYPE(GaussIntegrationPoints_t), TARGET :: IntegStuff

     INTERFACE
        
        SUBROUTINE IBOF(ai,a4)
           USE Types
           REAL(KIND=dp),intent(in) :: ai(6)
           REAL(KIND=dp),intent(out) :: a4(9)
        END SUBROUTINE
 
     END INTERFACE
!------------------------------------------------------------------------------
      Fond=.False.
      hmax = maxval (Nodes % y(1:n))
        
     dim = CoordinateSystemDimension()

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
!     Orientation parameters at the integration point:
!     ------------------------------------------------
      A1 = SUM( NodalK1(1:n) * Basis(1:n) ) 
      A2 = SUM( NodalK2(1:n) * Basis(1:n) )
      A3 = 1._dp - A1 - A2

      E1 = SUM( NodalEuler1(1:n) * Basis(1:n) )
      E2 = SUM( NodalEuler2(1:n) * Basis(1:n) )
      E3 = SUM( NodalEuler3(1:n) * Basis(1:n) )

!      Strain-Rate, Stresses and Spin

      CSymmetry = CurrentCoordinateSystem() == AxisSymmetric
      
      StrainRate = 0.0
      Spin1 = 0.0
!
!    Material parameters at that point
!    ---------------------------------
!
      ai(1)=A1
      ai(2)=A2
      ai(3)=A3
      ai(4)=E1
      ai(5)=E2
      ai(6)=E3

!    fourth order orientation tensor
      call IBOF(ai,a4)

!     A2 expressed in the orthotropic frame
!
!      call R2Ro(ai,dim,ap,angle)

!    Compute strainRate and Spin :
!    -----------------------------

      LGrad = MATMUL( NodalVelo(:,1:n), dBasisdx(1:n,:) )

      StrainRate = 0.5 * ( LGrad + TRANSPOSE(LGrad) )

!      Spin1 = 0.5 * ( LGrad - TRANSPOSE(LGrad) )

      IF ( CSymmetry ) THEN
        StrainRate(1,3) = 0.0
        StrainRate(2,3) = 0.0
        StrainRate(3,1) = 0.0
        StrainRate(3,2) = 0.0
        StrainRate(3,3) = 0.0

        Radius = SUM( Nodes % x(1:n) * Basis(1:n) )

        IF ( Radius > 10*AEPS ) THEN
         StrainRate(3,3) = SUM( Nodalvelo(1,1:n) * Basis(1:n) ) / Radius
        END IF

        epsi = StrainRate(1,1)+StrainRate(2,2)+StrainRate(3,3)
        DO i=1,3   
          StrainRate(i,i) = StrainRate(i,i) - epsi/3.0
        END DO

      ELSE
        epsi = StrainRate(1,1)+StrainRate(2,2)+StrainRate(3,3)
        DO i=1,dim 
          StrainRate(i,i) = StrainRate(i,i) - epsi/dim
        END DO

      END IF
      

      INDi(1:6) = (/ 1, 2, 3, 1, 2, 3 /)
      INDj(1:6) = (/ 1, 2, 3, 2, 3, 1 /)

   
     Do i=1,2*dim-3
       Spin(i)=Spin1(INDi(i+3),INDj(i+3))
     End do

!     Velocity :
!     ----------
      Velo = 0.0d0
      DO i=1,dim
         Velo(i) = SUM( Basis(1:n) * (NodalVelo(i,1:n) - NodMeshVel(i,1:n)) )
      END DO
      Unorm = SQRT( SUM( Velo**2._dp ) )

!
!     Reaction coefficient:
!     ---------------------
      SELECT CASE(comp)
      CASE(1)
        C0 = 2._dp*iota*(StrainRate(1,1)-StrainRate(3,3)) 

      CASE(2)
        C0 = 2._dp*iota*(StrainRate(2,2)-StrainRate(3,3))
       
      CASE(3)
        C0 = iota*(StrainRate(1,1)+StrainRate(2,2)-2._dp*StrainRate(3,3))
        
      CASE(4)
        C0 = iota*(StrainRate(2,2)-StrainRate(3,3))

      CASE(5)
        C0 = iota*(StrainRate(1,1)-StrainRate(3,3)) 
        

        
      END SELECT

      If (Fond) C0=0._dp

!     Loop over basis functions (of both unknowns and weights):
!     ---------------------------------------------------------
      DO p=1,NBasis
         DO q=1,NBasis
            A = 0.0d0
            M = Basis(p) * Basis(q)
!
!           Reaction terms:
!           ---------------
            A = A + C0 * Basis(q) * Basis(p)

            !
            ! Advection terms:
            ! ----------------
            DO j=1,dim
               A = A - Velo(j) * Basis(q) * dBasisdx(p,j)
            END DO

!           Add nodal matrix to element matrix:
!           -----------------------------------
            MASS( p,q )  = MASS( p,q )  + s * M
            STIFF( p,q ) = STIFF( p,q ) + s * A
         END DO


!
!        The righthand side...:
!        ----------------------
         SELECT CASE(comp)
         
         CASE(1)
          IF(dim == 3) THEN 
            LoadAtIp = 2._dp*(Spin(1)*E1-Spin(3)*E3)+ &
                 2._dp*iota*(StrainRate(1,2)*(2._dp*a4(7)-E1)+ &
                 StrainRate(3,1)*(2._dp*a4(6)-E3)+a4(1)*StrainRate(1,1)+ &
                 a4(3)*StrainRate(2,2)-(a4(1)+a4(3))*StrainRate(3,3)+&
                 2._dp*a4(4)*StrainRate(2,3))
          ELSE
              LoadAtIp = 2._dp*Spin(1)*E1+ &
                 2._dp*iota*(StrainRate(1,2)*(2._dp*a4(7)-E1)+ &
                 a4(1)*StrainRate(1,1)+a4(3)*StrainRate(2,2)- &
                 (a4(1)+a4(3))*StrainRate(3,3)) 
          END IF
         
         CASE(2)
          IF(dim == 3) THEN
            LoadAtIp = 2._dp*(-Spin(1)*E1+Spin(2)*E2)+ &
                 2._dp*iota*(StrainRate(1,2)*(2._dp*a4(9)-E1)+ &
                 StrainRate(2,3)*(2._dp*a4(8)-E2)+a4(3)*StrainRate(1,1)+ &
                 a4(2)*StrainRate(2,2)-(a4(2)+a4(3))*StrainRate(3,3)+&
                 2._dp*a4(5)*StrainRate(3,1)) 
          ELSE
            LoadAtIp = -2._dp*Spin(1)*E1+ &
                 2._dp*iota*(StrainRate(1,2)*(2._dp*a4(9)-E1)+ &
                 a4(3)*StrainRate(1,1)+a4(2)*StrainRate(2,2)- &
                 (a4(2)+a4(3))*StrainRate(3,3)) 
          END IF

         CASE(3)
          IF(dim == 3) THEN
            LoadAtIp = Spin(1)*(A2-A1)-Spin(3)*E2+Spin(2)*E3+ &
                 StrainRate(1,2)*iota*(4._dp*a4(3)-(A1+A2))+ &
                 StrainRate(3,1)*iota*(4._dp*a4(4)-E2)+ &
                 StrainRate(2,3)*iota*(4._dp*a4(5)-E3)+ &
                 2._dp*iota*(a4(7)*StrainRate(1,1)+a4(9)*StrainRate(2,2)- &
                 (a4(7)+a4(9))*StrainRate(3,3))  
          ELSE
              LoadAtIp = Spin(1)*(A2-A1)+ &
                 StrainRate(1,2)*iota*(4._dp*a4(3)-(A1+A2))+ &
                 2._dp*iota*(a4(7)*StrainRate(1,1)+a4(9)*StrainRate(2,2)- &
                 (a4(7)+a4(9))*StrainRate(3,3))
          END IF

         CASE(4)
          LoadAtIp = Spin(2)*(A3-A2)+Spin(3)*E1-Spin(1)*E3+ &
                 StrainRate(1,2)*iota*(4._dp*a4(5)-E3)+ &
                 StrainRate(2,3)*iota*(3._dp*A2-A3-4._dp*(a4(2)+a4(3)))+ &
                 StrainRate(3,1)*iota*(3._dp*E1-4._dp*(a4(7)+a4(9)))+ &
                 2._dp*iota*(a4(4)*StrainRate(1,1)+a4(8)*StrainRate(2,2)- &
                 (a4(4)+a4(8))*StrainRate(3,3))
         
         CASE(5)
          LoadAtIp = Spin(3)*(A1-A3)+Spin(1)*E2-Spin(2)*E1+ &
                 StrainRate(1,2)*iota*(4._dp*a4(4)-E2) + &
                 StrainRate(3,1)*iota*(3._dp*A1-A3-4._dp*(a4(1)+a4(3)))+ &
                 StrainRate(2,3)*iota*(3._dp*E1-4._dp*(a4(7)+a4(9)))+ &
                 2._dp*iota*(a4(6)*StrainRate(1,1)+a4(5)*StrainRate(2,2)- &
                 (a4(6)+a4(5))*StrainRate(3,3))

         END SELECT

        If (Fond) LoadAtIp=0._dp
        LoadAtIp= LoadAtIp * Basis(p)
        FORCE(p) = FORCE(p) + s*LoadAtIp
      END DO
   END DO 

 1000 FORMAT((a),x,i2,x,6(e13.5,x)) 
 1001 FORMAT(6(e13.5,x))

!------------------------------------------------------------------------------
    END SUBROUTINE LocalMatrix
!------------------------------------------------------------------------------


!!! Compute fourth order tensor a4 from a2 with closure function IBOF (Chung, 2002)
!!! a2 enters in the order : 11, 22, 33, 12, 23 ,13
!!! Output for a4 is in the order : 1111, 2222, 1122, 1123, 2231, 1131, 1112, 2223, 2212
!!! Code modified from Gillet-Chaullet source
!------------------------------------------------------------------------------
    SUBROUTINE IBOF(a2,a4)
!------------------------------------------------------------------------------

      USE Types
       
       implicit none
       Real(dp),dimension(6),intent(in):: a2  
       Real(dp),dimension(9),intent(out):: a4  
       Real(dp):: a_11,a_22,a_33,a_12,a_13,a_23
       Real(dp):: b_11,b_22,b_12,b_13,b_23
       Real(dp):: aPlusa

       Real(dp),dimension(21) :: vec
       Real(dp),dimension(3,21) :: Mat
       Real(dp),dimension(6) :: beta
       Real(dp) :: Inv2,Inv3
       integer :: i,j      

       
       a_11=a2(1)
       a_22=a2(2)
       a_33=a2(3)
       a_12=a2(4)
       a_23=a2(5)
       a_13=a2(6)


      !Coefficiants 

      Mat(1,1)=0.217774509809788e+02_dp
      Mat(1,2)=-.297570854171128e+03_dp
      Mat(1,3)=0.188686077307885e+04_dp
      Mat(1,4)=-.272941724578513e+03_dp
      Mat(1,5)=0.417148493642195e+03_dp
      Mat(1,6)=0.152038182241196e+04_dp
      Mat(1,7)=-.137643852992708e+04_dp
      Mat(1,8)=-.628895857556395e+03_dp
      Mat(1,9)=-.526081007711996e+04_dp
      Mat(1,10)=-.266096234984017e+03_dp
      Mat(1,11)=-.196278098216953e+04_dp
      Mat(1,12)=-.505266963449819e+03_dp
      Mat(1,13)=-.110483041928547e+03_dp
      Mat(1,14)=0.430488193758786e+04_dp
      Mat(1,15)=-.139197970442470e+02_dp
      Mat(1,16)=-.144351781922013e+04_dp
      Mat(1,17)=-.265701301773249e+03_dp
      Mat(1,18)=-.428821699139210e+02_dp
      Mat(1,19)=-.443236656693991e+01_dp
      Mat(1,20)=0.309742340203200e+04_dp
      Mat(1,21)=0.386473912295113e+00_dp
      Mat(2,1)=-.514850598717222e+00_dp
      Mat(2,2)=0.213316362570669e+02_dp
      Mat(2,3)=-.302865564916568e+03_dp
      Mat(2,4)=-.198569416607029e+02_dp
      Mat(2,5)=-.460306750911640e+02_dp
      Mat(2,6)=0.270825710321281e+01_dp
      Mat(2,7)=0.184510695601404e+03_dp
      Mat(2,8)=0.156537424620061e+03_dp
      Mat(2,9)=0.190613131168980e+04_dp
      Mat(2,10)=0.277006550460850e+03_dp
      Mat(2,11)=-.568117055198608e+02_dp
      Mat(2,12)=0.428921546783467e+03_dp
      Mat(2,13)=0.142494945404341e+03_dp
      Mat(2,14)=-.541945228489881e+04_dp
      Mat(2,15)=0.233351898912768e+02_dp
      Mat(2,16)=0.104183218654671e+04_dp
      Mat(2,17)=0.331489412844667e+03_dp
      Mat(2,18)=0.660002154209991e+02_dp
      Mat(2,19)=0.997500770521877e+01_dp
      Mat(2,20)=0.560508628472486e+04_dp
      Mat(2,21)=0.209909225990756e+01_dp
      Mat(3,1)=0.203814051719994e+02_dp
      Mat(3,2)=-.283958093739548e+03_dp
      Mat(3,3)=0.173908241235198e+04_dp
      Mat(3,4)=-.195566197110461e+03_dp
      Mat(3,5)=-.138012943339611e+03_dp
      Mat(3,6)=0.523629892715050e+03_dp
      Mat(3,7)=0.859266451736379e+03_dp
      Mat(3,8)=-.805606471979730e+02_dp
      Mat(3,9)=-.468711180560599e+04_dp
      Mat(3,10)=0.889580760829066e+01_dp
      Mat(3,11)=-.782994158054881e+02_dp
      Mat(3,12)=-.437214580089117e+02_dp
      Mat(3,13)=0.112996386047623e+01_dp
      Mat(3,14)=0.401746416262936e+04_dp
      Mat(3,15)=0.104927789918320e+01_dp
      Mat(3,16)=-.139340154288711e+03_dp
      Mat(3,17)=-.170995948015951e+02_dp
      Mat(3,18)=0.545784716783902e+00_dp
      Mat(3,19)=0.971126767581517e+00_dp
      Mat(3,20)=0.141909512967882e+04_dp
      Mat(3,21)=0.994142892628410e+00_dp

       
      ! Compute the invariants
      Inv2=0.5_dp*(1._dp-(a_11*a_11+a_22*a_22+a_33*a_33+ &
            2._dp*(a_12*a_12+a_13*a_13+a_23*a_23)))
            
       Inv3=a_11*(a_22*a_33-a_23*a_23)+a_12*(a_23*a_13-a_12*a_33)+ &
             a_13*(a_12*a_23-a_22*a_13)
       
     ! complete polynome of degree 5 for the 2 invariants.
         vec(1)=1._dp
         vec(2)=Inv2
         vec(3)=vec(2)*vec(2)
         vec(4)=Inv3
         vec(5)=vec(4)*vec(4)
         vec(6)=vec(2)*vec(4)
         vec(7)=vec(3)*vec(4)
         vec(8)=vec(2)*vec(5)
         vec(9)=vec(2)*vec(3)
         vec(10)=vec(5)*vec(4)
         vec(11)=vec(9)*vec(4)
         vec(12)=vec(3)*vec(5)
         vec(13)=vec(2)*vec(10)
         vec(14)=vec(3)*vec(3)
         vec(15)=vec(5)*vec(5)
         vec(16)=vec(14)*vec(4)
         vec(17)=vec(12)*vec(2)
         vec(18)=vec(12)*vec(4)
         vec(19)=vec(2)*vec(15)
         vec(20)=vec(14)*vec(2)
         vec(21)=vec(15)*vec(4)

       ! Compites beta_bar (cf annexe C Chung)
       ! Warning: beta(1)=beta_bar_3 (Chung); beta(2)=beta_bar_4; beta(3)=beta_bar_6
       !           beta(4)=beta_bar_1        ; beta(5)=beta_bar_2; beta(6)=beta_bar_5

       ! calcul the three betas in terms of the polynomes
         beta(:)=0._dp
         Do i=1,3
          Do j=1,21
            beta(i)=beta(i)+Mat(i,j)*vec(j)
          End do
         End do
          
       ! calcul the other 3 to get the normalisation
         beta(4)=3._dp*(-1._dp/7._dp+beta(1)*(1._dp/7._dp+4._dp*Inv2/7._dp+8._dp*Inv3/3._dp)/5._dp- &
                  beta(2)*(0.2_dp-8._dp*Inv2/15._dp-14._dp*Inv3/15._dp)- &
                  beta(3)*(1._dp/35._dp-24._dp*Inv3/105._dp-4._dp*Inv2/35._dp+ &
                  16._dp*Inv2*Inv3/15._dp+8._dp*Inv2*Inv2/35._dp))/5._dp

         beta(5)=6._dp*(1._dp-0.2_dp*beta(1)*(1._dp+4._dp*Inv2)+ &
                  7._dp*beta(2)*(1._dp/6._dp-Inv2)/5._dp- &
                  beta(3)*(-0.2_dp+2._dp*Inv3/3._dp+4._dp*Inv2/5._dp- &
                  8._dp*Inv2*Inv2/5._dp))/7._dp

         beta(6)=-4._dp*beta(1)/5._dp-7._dp*beta(2)/5._dp- &
                   6._dp*beta(3)*(1._dp-4._dp*Inv2/3._dp)/5._dp

        !beta_bar
        Do i=1,6
         beta(i)=beta(i)/3._dp
        End do
         beta(2)=beta(2)/2._dp
         beta(5)=beta(5)/2._dp
         beta(6)=beta(6)/2._dp

        !! Compute 5 b=a.a
        b_11=a_11*a_11+a_12*a_12+a_13*a_13
        b_22=a_22*a_22+a_12*a_12+a_23*a_23
        b_12=a_11*a_12+a_12*a_22+a_13*a_23
        b_13=a_11*a_13+a_12*a_23+a_13*a_33
        b_23=a_12*a_13+a_22*a_23+a_23*a_33

        !Compute the  9 terms of a4

        a4(1)=3._dp*beta(4)+6._dp*beta(5)*a_11+3._dp*beta(1)*a_11*a_11+&
         6._dp*beta(2)*b_11+6._dp*beta(6)*a_11*b_11+3._dp*beta(3)*b_11*b_11
        a4(2)=3._dp*beta(4)+6._dp*beta(5)*a_22+3._dp*beta(1)*a_22*a_22+&
         6._dp*beta(2)*b_22+6._dp*beta(6)*a_22*b_22+3._dp*beta(3)*b_22*b_22

        a4(3)=beta(4)+beta(5)*(a_22+a_11)+beta(1)*(a_11*a_22+2._dp*a_12*a_12)+&
         beta(2)*(b_22+b_11)+beta(6)*(a_11*b_22+a_22*b_11+4._dp*a_12*b_12)+&
         beta(3)*(b_11*b_22+2._dp*b_12*b_12)


         a4(4)=beta(5)*a_23+beta(1)*(a_11*a_23+2._dp*a_12*a_13)+beta(2)*b_23+&
          beta(6)*(a_11*b_23+a_23*b_11+2._dp*(a_12*b_13+a_13*b_12))+beta(3)*&
          (b_11*b_23+2._dp*b_12*b_13)
         a4(5)=beta(5)*a_13+beta(1)*(a_22*a_13+2._dp*a_12*a_23)+beta(2)*b_13+&
          beta(6)*(a_22*b_13+a_13*b_22+2._dp*(a_12*b_23+a_23*b_12))+beta(3)*&
          (b_22*b_13+2._dp*b_12*b_23)


         a4(6)=3._dp*beta(5)*a_13+3._dp*beta(1)*a_11*a_13+3._dp*beta(2)*b_13+&
          3._dp*beta(6)*(a_11*b_13+a_13*b_11)+3._dp*beta(3)*b_11*b_13
         a4(7)=3._dp*beta(5)*a_12+3._dp*beta(1)*a_11*a_12+3._dp*beta(2)*b_12+&
          3._dp*beta(6)*(a_11*b_12+a_12*b_11)+3._dp*beta(3)*b_11*b_12
         a4(8)=3._dp*beta(5)*a_23+3._dp*beta(1)*a_22*a_23+3._dp*beta(2)*b_23+&
          3._dp*beta(6)*(a_22*b_23+a_23*b_22)+3._dp*beta(3)*b_22*b_23
         a4(9)=3._dp*beta(5)*a_12+3._dp*beta(1)*a_22*a_12+3._dp*beta(2)*b_12+&
          3._dp*beta(6)*(a_22*b_12+a_12*b_22)+3._dp*beta(3)*b_22*b_12
!------------------------------------------------------------------------------
    END SUBROUTINE IBOF
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

      TYPE(Nodes_t) :: EdgeNodes, LeftParentNodes, RightParentNodes

      Save EdgeNodes, LeftParentNodes, RightParentNodes
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
            STIFF(p,q) = STIFF(p,q) + s * Udotn * Average(q) * Jump(p)
            STIFF(p,q) = STIFF(p,q) + s * ABS(Udotn)/2 * Jump(q) * Jump(p)
          END DO
        END DO
     END DO
!------------------------------------------------------------------------------
   END SUBROUTINE LocalJumps
!------------------------------------------------------------------------------



!------------------------------------------------------------------------------
   SUBROUTINE LocalMatrixBoundary( STIFF, FORCE, LOAD, &
        Element, n, ParentElement, np, Velo,MeshVelo, InFlowBC )
!------------------------------------------------------------------------------
     REAL(KIND=dp) :: STIFF(:,:),  FORCE(:), LOAD(:), Velo(:,:),MeshVelo(:,:)
     INTEGER :: n, np
     LOGICAL :: InFlowBC
     TYPE(Element_t), POINTER :: Element, ParentElement
!------------------------------------------------------------------------------
     REAL(KIND=dp) :: Basis(n), dBasisdx(n,3), ddBasisddx(n,3,3)
     REAL(KIND=dp) :: ParentBasis(np), ParentdBasisdx(np,3), ParentddBasisddx(np,3,3)
     INTEGER :: i,j,p,q,t,dim

     REAL(KIND=dp) :: Normal(3), g, L, Udotn, UdotnA, cu(3), detJ,U,V,W,S
     LOGICAL :: Stat
     TYPE(GaussIntegrationPoints_t) :: IntegStuff

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

       L = SUM( LOAD(1:n) * Basis(1:n) )
       cu = 0.0d0
       DO i=1,dim
          cu(i) = SUM( (Velo(i,1:n)-MeshVelo(i,1:n)) * Basis(1:n) )
       END DO
       Udotn = SUM( Normal * cu )

       DO p = 1,np
        IF (InFlowBC .And. (UdotnA < 0.) ) THEN
            FORCE(p) = FORCE(p) - s * Udotn*L*ParentBasis(p)
         ELSE
           DO q=1,np
             STIFF(p,q) = STIFF(p,q) + s*Udotn*ParentBasis(q)*ParentBasis(p)
           END DO
         END IF
       END DO
     END DO
!------------------------------------------------------------------------------
   END SUBROUTINE LocalMatrixBoundary
!------------------------------------------------------------------------------



!------------------------------------------------------------------------------
 END SUBROUTINE CaffeSolveFabric
!------------------------------------------------------------------------------
