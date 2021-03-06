! *****************************************************************************/
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
! * Module containing a solver for computing deviatoric or Cauchy 
! *          stress from flow solution (Restricted to NS solver) 
! * 2D SDOFs = 4 (S11, S22, S33, S12)
! * 3D SDOFs = 6 (S11, S22, S33, S12, S23, S31)
! * Keywords : Cauchy (Logical), 
! *
! ******************************************************************************
! *
! *                    Author:       Juha Ruokolainen 
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
! *                Modified by:  Fabien / OG 
! *                              Hakime Seddik
! *
! *       Date of modification: 13/10/05 from version 1.5
! *
! *****************************************************************************/
!------------------------------------------------------------------------------
   RECURSIVE SUBROUTINE ComputeDevStress( Model,Solver,dt,TransientSimulation )
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
     TYPE(Solver_t), POINTER :: PSolver

     TYPE(Matrix_t),POINTER :: StiffMatrix

     INTEGER :: i, j, k, l, n, t, iter, NDeg
     INTEGER :: dim, STDOFs, StressDOFs, LocalNodes, istat

     TYPE(ValueList_t),POINTER :: Material, BC
     TYPE(Nodes_t) :: ElementNodes
     TYPE(Element_t),POINTER :: CurrentElement

     REAL(KIND=dp) :: RelativeChange, UNorm, PrevUNorm

     REAL(KIND=dp), ALLOCATABLE :: Basis(:),ddBasisddx(:,:,:)
     REAL(KIND=dp), ALLOCATABLE :: dBasisdx(:,:)
     REAL(KIND=dp) :: u, v, w, detJ

     LOGICAL :: stat, CSymmetry 

     INTEGER :: NewtonIter, NonlinearIter

     TYPE(Variable_t), POINTER :: StressSol, FlowVariable 

     REAL(KIND=dp), POINTER ::  Stress(:), Solution(:), &
           ForceVector(:), FlowValues(:) 

     INTEGER, POINTER :: StressPerm(:), NodeIndexes(:), &
                         FlowPerm(:)

     INTEGER :: body_id
     INTEGER :: old_body = -1

     LOGICAL :: Isotropic, AllocationsDone = .FALSE.,  &
                Requal0
     LOGICAL :: GotIt,  Cauchy = .FALSE.

     REAL(KIND=dp), ALLOCATABLE:: LocalMassMatrix(:,:), &
       LocalStiffMatrix(:,:), LocalForce(:), &
       LocalP(:),  &
       LocalVelo(:,:), LocalViscosity(:)

     INTEGER :: NumberOfBoundaryNodes, COMP
     INTEGER, POINTER :: BoundaryReorder(:)

     REAL(KIND=dp), POINTER :: BoundaryNormals(:,:), &
         BoundaryTangent1(:,:), BoundaryTangent2(:,:)
     CHARACTER(LEN=MAX_NAME_LEN) :: FlowSolverName

     REAL(KIND=dp) :: at, at0, CPUTime, RealTime


!------------------------------------------------------------------------------
     SAVE NumberOfBoundaryNodes, BoundaryReorder, BoundaryNormals, &
              BoundaryTangent1, BoundaryTangent2

     SAVE Basis, dBasisdx, ddBasisddx
     SAVE LocalMassMatrix, LocalStiffMatrix, LocalForce, &
          ElementNodes,  &
          AllocationsDone,  &
          old_body, &
          LocalViscosity, Cauchy
     SAVE LocalVelo, LocalP, dim

!------------------------------------------------------------------------------
!  Read the name of the Flow Solver (NS)
!------------------------------------------------------------------------------

     FlowSolverName = GetString( Solver % Values, 'Flow Solver Name', GotIt )    
     IF (.NOT.Gotit) FlowSolverName = 'Flow Solution'
     FlowVariable => VariableGet( Solver % Mesh % Variables, FlowSolverName )
     IF ( ASSOCIATED( FlowVariable ) ) THEN
       FlowPerm    => FlowVariable % Perm
       FlowValues  => FlowVariable % Values
     ELSE
       CALL Info('ComputeDevStress', &
                      & 'No variable for velocity associated.', Level=4)
     END IF
!              
!------------------------------------------------------------------------------
!  Read constants from constants section of SIF file
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!    Get variables needed for solution
!------------------------------------------------------------------------------
      IF ( .NOT. ASSOCIATED( Solver % Matrix ) ) RETURN

      Solution => Solver % Variable % Values
      STDOFs   =  Solver % Variable % DOFs

      IF ( STDOFs /=1 ) THEN
              CALL Fatal( 'ComputeDevStress', 'DOF must be equal to 1' )
      END IF

      StressSol => VariableGet( Solver % Mesh % Variables, 'Stress' )
      StressPerm => StressSol % Perm
      StressDOFs = StressSol % DOFs
      Stress => StressSol % Values

      StiffMatrix => Solver % Matrix
      ForceVector => StiffMatrix % RHS
      Unorm = SQRT( SUM( Stress**2 ) / SIZE(Stress) )

!------------------------------------------------------------------------------
!     Allocate some permanent storage, this is done first time only
!------------------------------------------------------------------------------
      IF ( .NOT. AllocationsDone .OR. Solver % Mesh % Changed) THEN
        N = Model % MaxElementNodes
        dim = CoordinateSystemDimension()

       IF ( AllocationsDone ) THEN
         DEALLOCATE( ElementNodes % x,     &
                     ElementNodes % y,     &
                     ElementNodes % z,     &
                     LocalVelo, LocalP,    &                      
                     Basis, ddBasisddx,    &
                     dBasisdx,             &
                     LocalMassMatrix,      &
                     LocalStiffMatrix,     &
                     LocalForce,           &
                     LocalViscosity )
       END IF

       ALLOCATE( ElementNodes % x( N ), &
                 ElementNodes % y( N ), &
                 ElementNodes % z( N ), &
                 LocalVelo( 3,N ), LocalP( N ), &                                     
                 Basis( 2*N ),ddBasisddx(1,1,1), dBasisdx( 2*N,3 ), &
                 LocalMassMatrix( 2*STDOFs*N,2*STDOFs*N ),  &
                 LocalStiffMatrix( 2*STDOFs*N,2*STDOFs*N ),  &
                 LocalForce( 2*STDOFs*N ),  &
                 LocalViscosity(N), STAT=istat )

       IF ( istat /= 0 ) THEN
          CALL Fatal( 'ComputeDevStress', 'Memory allocation error.' )
       END IF
!------------------------------------------------------------------------------

       AllocationsDone = .TRUE.
      END IF


!------------------------------------------------------------------------------
      NonlinearIter = 1
      DO iter=1,NonlinearIter

       at  = CPUTime()
       at0 = RealTime()

       CALL Info( 'ComputeDevStress', ' ', Level=4 )
       CALL Info( 'ComputeDevStress', ' ', Level=4 )
       CALL Info( 'ComputeDevStress', ' ', Level=4 )
       CALL Info( 'ComputeDevStress', ' ', Level=4 )
       CALL Info( 'ComputeDevStress', 'Starting assembly...',Level=4 )

! Loop over the Stress components [Sxx, Syy, Szz, Sxy, Syz, Szx] 

       PrevUNorm = UNorm

       DO COMP = 1, 2*dim

        WRITE(Message,'(a,i3)' ) ' Component : ', COMP  
        CALL Info( 'ComputeDevStress', Message, Level=5 )


!------------------------------------------------------------------------------
       CALL DefaultInitialize()
!------------------------------------------------------------------------------
       DO t=1,Solver % NumberOFActiveElements

         IF ( RealTime() - at0 > 1.0 ) THEN
           WRITE(Message,'(a,i3,a)' ) '   Assembly: ',  &
             INT(100.0 - 100.0 * (Solver % NumberOfActiveElements-t) / &
             (1.0*Solver % NumberOfActiveElements)), ' % done'
           CALL Info( 'ComputeDevStress', Message, Level=5 )
           at0 = RealTime()
         END IF

         CurrentElement => GetActiveElement(t)
         n = GetElementNOFNodes()
         NodeIndexes => CurrentElement % NodeIndexes

         ElementNodes % x(1:n) = Model % Nodes % x(NodeIndexes(1:n))
         ElementNodes % y(1:n) = Model % Nodes % y(NodeIndexes(1:n))
         ElementNodes % z(1:n) = Model % Nodes % z(NodeIndexes(1:n))

         Material => GetMaterial()

!------------------------------------------------------------------------------
!    Read in material constants from Material section
!------------------------------------------------------------------------------


!!!! Restricted to the Power Law case

          Cauchy = ListGetLogical( Material , 'Cauchy', Gotit )
          IF (.NOT.Gotit) THEN
              Cauchy = .FALSE.
              WRITE(Message,'(A)') 'Cauchy set to False'
              CALL INFO('ComputeDevStress', Message, Level = 20)
          END IF

          LocalViscosity(1:n) = ListGetReal( Material, &
                         'Viscosity', n, NodeIndexes, GotIt )
          IF (.NOT.GotIt) THEN
            WRITE(Message,'(A)') 'Variable Viscosity not found. &
                            &Setting to 1.0'
            CALL INFO('ComputeStress', Message, Level = 20)
            LocalViscosity(1:n) = 1.0_dp
          END IF

          LocalVelo = 0.0_dp
          DO i=1, dim
             LocalVelo(i,1:n) = FlowValues((dim+1)*(FlowPerm(NodeIndexes(1:n))-1) + i)
          END DO
          LocalP(1:n) = FlowValues((dim+1)*(FlowPerm(NodeIndexes(1:n))-1)+(dim+1))

          CALL LocalNSMatrix(COMP, LocalMassMatrix, LocalStiffMatrix, &
              LocalForce,  LocalVelo, LocalP, &
              LocalViscosity, CurrentElement, n, &
              ElementNodes, Cauchy)

!------------------------------------------------------------------------------
!        Update global matrices from local matrices 
!------------------------------------------------------------------------------
         CALL DefaultUpdateEquations( LocalStiffMatrix, LocalForce )

      END DO

      CALL Info( 'ComputeDevStress', 'Assembly done', Level=4 )


      CALL DefaultFinishAssembly()

!------------------------------------------------------------------------------
!     Dirichlet boundary conditions
!------------------------------------------------------------------------------
      CALL DefaultDirichletBCs()

!------------------------------------------------------------------------------

      CALL Info( 'ComputeDevStress', 'Set boundaries done', Level=4 )

!------------------------------------------------------------------------------
!     Solve the system and check for convergence
!------------------------------------------------------------------------------
      PrevUNorm = UNorm

      UNorm = DefaultSolve()

      DO t=1,Solver % NumberOfActiveElements
         CurrentElement => GetActiveElement(t) 
         n = GetElementNOFNodes()
         DO i=1,n
            k = CurrentElement % NodeIndexes(i)
            Stress( StressDOFs*(StressPerm(k)-1) + COMP ) =    & 
            Solver % Variable % Values( Solver % Variable % Perm(k) )
         END DO
      END DO

      END DO ! End DO Comp


      Unorm = SQRT( SUM( Stress**2 ) / SIZE(Stress) )
      Solver % Variable % Norm = Unorm  

      IF ( PrevUNorm + UNorm /= 0.0d0 ) THEN
         RelativeChange = 2.0d0 * ABS( PrevUNorm - UNorm) / ( PrevUnorm + UNorm)
      ELSE
         RelativeChange = 0.0d0
      END IF

      WRITE( Message, * ) 'Result Norm   : ',UNorm, PrevUNorm
      CALL Info( 'ComputeDevStress', Message, Level=4 )
      WRITE( Message, * ) 'Relative Change : ',RelativeChange
      CALL Info( 'ComputeDevStress', Message, Level=4 )


!------------------------------------------------------------------------------
    END DO ! of nonlinear iter
!------------------------------------------------------------------------------


CONTAINS


!------------------------------------------------------------------------------
      SUBROUTINE LocalNSMatrix(COMP, MassMatrix, StiffMatrix, ForceVector, &
               NodalVelo, NodalP, NodalViscosity, &
              Element, n, Nodes, Cauchy )
!------------------------------------------------------------------------------

     USE MaterialModels

     REAL(KIND=dp) :: StiffMatrix(:,:), MassMatrix(:,:)
     REAL(KIND=dp) ::  NodalVelo(:,:)
     REAL(KIND=dp), DIMENSION(:) :: ForceVector,  &
                              NodalViscosity, NodalP
     TYPE(Nodes_t) :: Nodes
     TYPE(Element_t), POINTER :: Element
     LOGICAL ::  Cauchy
     INTEGER :: n, COMP
!------------------------------------------------------------------------------
!
     REAL(KIND=dp) :: Basis(2*n),ddBasisddx(1,1,1)
     REAL(KIND=dp) :: dBasisdx(2*n,3),detJ, pBasis(n)

     REAL(KIND=dp) :: Stress, epsi

     REAL(KIND=dp) :: Pressure
     REAL(KIND=dp) :: LGrad(3,3), SR(3,3)

     INTEGER :: i, j, k, p, q, t, dim, cc, NBasis,  LinearBasis

     REAL(KIND=dp) :: s, u, v, w, Radius, eta

     REAL(KIND=dp) :: dDispldx(3,3), Viscosity
     TYPE(GaussIntegrationPoints_t), TARGET :: IntegStuff

     INTEGER, POINTER :: EdgeMap(:,:)
     INTEGER :: N_Integ, nd
     INTEGER, DIMENSION(6), PARAMETER :: indx = (/1, 2, 3, 1, 2, 3/), &
                                         indy = (/1, 2, 3, 2, 3, 1/)

     REAL(KIND=dp), DIMENSION(:), POINTER :: U_Integ,V_Integ,W_Integ,S_Integ

     LOGICAL :: stat, CSymmetry

!------------------------------------------------------------------------------
      dim = CoordinateSystemDimension()
      cc=2*dim


      ForceVector = 0.0_dp
      StiffMatrix = 0.0_dp
      MassMatrix  = 0.0_dp

      IntegStuff = GaussPoints( Element )

      U_Integ => IntegStuff % u
      V_Integ => IntegStuff % v
      W_Integ => IntegStuff % w
      S_Integ => IntegStuff % s
      N_Integ =  IntegStuff % n
!
!   Now we start integrating
!
    DO t=1,N_Integ

      u = U_Integ(t)
      v = V_Integ(t)
      w = W_Integ(t)

!------------------------------------------------------------------------------
!     Basis function values & derivatives at the integration point
!------------------------------------------------------------------------------
       stat = ElementInfo(Element,Nodes,u,v,w,detJ, &
                 Basis,dBasisdx,ddBasisddx,.FALSE.,.FALSE.)

       s = detJ * S_Integ(t)

       Radius = SUM( Nodes % x(1:n) * Basis(1:n) )
       CSymmetry = CurrentCoordinateSystem() == AxisSymmetric
       IF ( CSymmetry ) s = s * Radius
!
!
        Viscosity = SUM( NodalViscosity(1:n)*Basis(1:n) )
        Viscosity = EffectiveViscosity( Viscosity, 1.0_dp, NodalVelo(1,1:n), NodalVelo(2,1:n), NodalVelo(3,1:n), &
                    Element, Nodes, n, n, u, v, w )

!
! Strain-Rate
!
        LGrad = MATMUL( NodalVelo(:,1:n), dBasisdx(1:n,:) )
        SR = 0.5 * ( LGrad + TRANSPOSE(LGrad) )
        IF ( CSymmetry ) THEN
          SR(1,3) = 0.0_dp
          SR(2,3) = 0.0_dp
          SR(3,1) = 0.0_dp
          SR(3,2) = 0.0_dp
          SR(3,3) = 0.0_dp
          IF ( Radius > 10*AEPS ) THEN
            SR(3,3) = SUM( Nodalvelo(1,1:n) * Basis(1:n) ) /Radius

          END IF
          epsi = SR(1,1)+SR(2,2)+SR(3,3)
          DO i=1,3   
            SR(i,i) = SR(i,i) - epsi/3.0_dp
          END DO
        ELSE
          epsi = SR(1,1)+SR(2,2)+SR(3,3)
          DO i=1,dim 
            SR(i,i) = SR(i,i) - epsi/dim
          END DO
        END IF
       ! write(*,*)'SR'
       ! write(*,*)SR(1,1), SR(1,2)
       ! write(*,*)SR(2,1), SR(2,2)
       ! write(*,*)SR(3,3)

!
!    Compute deviatoric stresses or Cauchy stresses: 
!    ----------------------------

      Stress = 2.0 * Viscosity * SR(indx(COMP),indy(COMP))
!      Write(*,*)Stress, Viscosity, SR(indx(COMP),indy(COMP))

      IF ((Cauchy).AND.(COMP.LE.3)) THEN
          Pressure = SUM( NodalP(1:n)*Basis(1:n) )
          Stress = Stress - Pressure
      END IF

      DO p=1,n         
        DO q=1,n        
          StiffMatrix(p,q) =  &
               StiffMatrix(p,q) + s*Basis(q)*Basis(p)
        END DO
          ForceVector(p) =  &
                     ForceVector(p) + s*Stress*Basis(p) 
      END DO

    END DO 

!------------------------------------------------------------------------------
      END SUBROUTINE LocalNSMatrix
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
      END SUBROUTINE ComputeDevStress
!------------------------------------------------------------------------------
