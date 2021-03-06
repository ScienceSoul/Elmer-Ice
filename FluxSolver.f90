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
! *****************************************************************************/
!
!/******************************************************************************
! *
! *  Module for computing fluxes of Poisson type of equations
! *
! ******************************************************************************
! *
! *  Authors: Peter R�back, Juha Ruokolainen, Hakime Seddik
! *  Email:   Peter.Raback@csc.fi
! *  Web:     http://www.csc.fi/elmer
! *  Address: CSC - IT Center for Science Ltd.
! *           Keilaranta 14
! *           02101 Espoo, Finland 
! *
! *  Original Date: 20.06.2007
! *
! *****************************************************************************/

!------------------------------------------------------------------------------
  SUBROUTINE FluxSolver_Init( Model,Solver,dt,Transient )
!------------------------------------------------------------------------------
    USE DefUtils

    TYPE(Model_t)  :: Model
    TYPE(Solver_t) :: Solver
    REAL(KIND=dp) :: DT
    LOGICAL :: Transient
!------------------------------------------------------------------------------
    TYPE(ValueList_t), POINTER :: SolverParams
    INTEGER :: dim
!------------------------------------------------------------------------------
    SolverParams => GetSolverParams()
    dim = CoordinateSystemDimension()

    IF ( .NOT. ListCheckPresent( SolverParams,'Variable') ) THEN
      CALL ListAddInteger( SolverParams, 'Variable DOFs', 1 )
      IF(ListCheckPresent(SolverParams,'Flux Component')) THEN
        CALL ListAddString( SolverParams,'Variable','Flux' )
      ELSE
        CALL ListAddString( SolverParams, 'Variable','-nooutput flux_temp' )
        CALL ListAddString(SolverParams, 'Flux Result Variable','F')
        IF(dim == 2) THEN
          CALL ListAddString( SolverParams,&
               NextFreeKeyword('Exported Variable',SolverParams),'F[Flux:2]')
        ELSE IF(dim == 3) THEN
          CALL ListAddString( SolverParams,&
               NextFreeKeyword('Exported Variable',SolverParams),'F[Flux:3]')
        ELSE
          CALL Fatal('VortictySolver_init','Flux computation makes sense only in 2D and 3D')
        END IF
      END IF
    END IF
    CALL ListAddInteger( SolverParams, 'Time derivative order', 0 )

    ! Add linear system defaults: cg+diagonal
    IF(.NOT. ListCheckPresent(SolverParams,'Linear System Solver')) &
      CALL ListAddString(SolverParams,'Linear System Solver','Iterative')
    IF(.NOT. ListCheckPresent(SolverParams,'Linear System Iterative Method')) &
      CALL ListAddString(SolverParams,'Linear System Iterative Method','cg')
    IF(.NOT. ListCheckPresent(SolverParams,'Linear System Preconditioning')) &
      CALL ListAddString(SolverParams,'Linear System Preconditioning','diagonal')
    IF(.NOT. ListCheckPresent(SolverParams,'Linear System Max Iterations')) &
      CALL ListAddInteger(SolverParams,'Linear System Max Iterations',500)
    IF(.NOT. ListCheckPresent(SolverParams,'Linear System Residual Output')) &
      CALL ListAddInteger(SolverParams,'Linear System Residual Output',10)
    IF(.NOT. ListCheckPresent(SolverParams,'Linear System Convergence Tolerance')) &
      CALL ListAddConstReal(SolverParams,'Linear System Convergence Tolerance',1.0e-10_dp)

!------------------------------------------------------------------------------
  END SUBROUTINE FluxSolver_Init
!------------------------------------------------------------------------------




!------------------------------------------------------------------------------
SUBROUTINE FluxSolver( Model,Solver,dt,Transient )
!------------------------------------------------------------------------------

  USE CoordinateSystems
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
  LOGICAL ::  Transient
  REAL(KIND=dp) :: dt
!------------------------------------------------------------------------------
!    Local variables
!------------------------------------------------------------------------------
  TYPE(ValueList_t),POINTER :: SolverParams
  CHARACTER(LEN=MAX_NAME_LEN) :: VarName, CondName
  INTEGER :: i,j,dim,DOFs,ActiveDir
  LOGICAL :: ConstantBulkMatrix, ConstantBulkMatrixInUse, CSymmetry
  LOGICAL :: GotIt, Visited = .FALSE., DirMask(3)
  REAL(KIND=dp) :: Unorm, Totnorm
  REAL(KIND=dp), POINTER :: ForceVector(:,:), SaveRHS(:)
  REAL(KIND=dp) :: at0,at1,at2,CPUTime,RealTime
  TYPE(Variable_t), POINTER :: FluxSol
  
  SAVE Visited

 
  CALL Info( 'FluxSolver', '-------------------------------------',Level=4 )
  CALL Info( 'FluxSolver','Computing the flux',Level=4 )
  CALL Info( 'FluxSolver', '-------------------------------------',Level=4 )

  dim = CoordinateSystemDimension()
!------------------------------------------------------------------------------
!    Get variables needed for solution
!------------------------------------------------------------------------------
  IF ( .NOT. ASSOCIATED( Solver % Matrix ) ) RETURN
  IF ( COUNT( Solver % Variable % Perm > 0 ) <= 0 ) RETURN
  
  SolverParams => GetSolverParams()

!-------------------------------------------------------------------------------
! If only one component is used use the scalar equation, otherwise use an
! auxiliary variable to store all the dimensions
!-------------------------------------------------------------------------------
  ActiveDir = GetInteger( SolverParams, 'Flux Component', GotIt) 
  IF(GotIt) THEN
    Dofs = 1
  ELSE
    VarName = GetString(SolverParams,'Flux Result Variable',GotIt )
    IF(.NOT. gotIt) VarName = 'Flux'
    FluxSol => VariableGet( Solver % Mesh % Variables, VarName )
    IF(ASSOCIATED(FluxSol)) THEN
      Dofs = FluxSol % DOFs
      IF(Dofs /= DIM) THEN
        CALL Fatal('FluxSolver','The flux should have DOFs equal to DIM')
      END IF
    ELSE
      CALL Fatal('FluxSolver','Flux Result Variable is missing: '//TRIM(VarName))      
    END IF
    DirMask = .TRUE.
  END IF
  
  CSymmetry = CurrentCoordinateSystem() == AxisSymmetric .OR. &
      CurrentCoordinateSystem() == CylindricSymmetric
  
  VarName = GetString(SolverParams,'Flux Variable',GotIt )
  IF(.NOT. gotIt) VarName = TRIM('Temperature')

  CondName = ListGetString(SolverParams,'Flux Coefficient',GotIt )
  IF(.NOT. gotIt) CondName = TRIM('Heat Conductivity')
  
  at0 = RealTime()
  
  ConstantBulkMatrix = GetLogical( SolverParams, 'Constant Bulk Matrix', GotIt )
  ConstantBulkMatrixInUse = ConstantBulkMatrix .AND. &
      ASSOCIATED(Solver % Matrix % BulkValues)
  
  IF ( ConstantBulkMatrixInUse ) THEN
    Solver % Matrix % Values = Solver % Matrix % BulkValues        
    Solver % Matrix % rhs = 0.0_dp
  ELSE
    CALL DefaultInitialize()
  END IF
  
  IF(Dofs > 1) THEN
    ALLOCATE(ForceVector(Dofs-1,SIZE(Solver % Matrix % RHS)))  
    ForceVector = 0.0_dp
    SaveRHS => Solver % Matrix % RHS
  END IF
  
  CALL BulkAssembly()
  CALL DefaultFinishAssembly()
  
  at1 = RealTime()
  WRITE(Message,* ) 'Assembly Time: ',at1-at0
  CALL Info( 'FluxSolver', Message, Level=5 )
!        
!------------------------------------------------------------------------------     

  IF(Dofs > 1) THEN
    TotNorm = 0.0_dp
    DO i=1,Dofs
      IF(i==1) THEN
        Solver % Matrix % RHS => SaveRHS
      ELSE
        Solver % Matrix % RHS => ForceVector(i-1,:)
      END IF
      UNorm = DefaultSolve()
      TotNorm = TotNorm + Unorm ** 2
      DO j=1,Solver % Matrix % NumberOfRows
        FluxSol % Values(DOFs*(j-1)+i) = Solver % Variable % Values(j)
      END DO
    END DO
    DEALLOCATE( ForceVector )  
    Solver % Matrix % RHS => SaveRHS
    TotNorm = SQRT(TotNorm)
    Solver % Variable % Norm = Totnorm
  ELSE
    TotNorm = DefaultSolve()
  END IF
!------------------------------------------------------------------------------     

  at2 = RealTime()
  WRITE(Message,* ) 'Solution Time: ',at2-at1
  CALL Info( 'FluxSolver', Message, Level=5 )
  
  WRITE( Message, * ) 'Result Norm: ',TotNorm
  CALL Info( 'FluxSolver', Message, Level=4 )
  
CONTAINS


!------------------------------------------------------------------------------
  SUBROUTINE BulkAssembly()
!------------------------------------------------------------------------------
       
    INTEGER :: elem,t,i,j,p,q,n,nd, Rank
    REAL(KIND=dp), ALLOCATABLE :: STIFF(:,:), FORCE(:,:)
    TYPE(GaussIntegrationPoints_t), TARGET :: IntegStuff
    TYPE(Nodes_t) :: Nodes
    TYPE(Element_t), POINTER :: Element
    REAL(KIND=dp) :: weight,grad(3),C(3,3),coeff,detJ
    REAL(KIND=dp), ALLOCATABLE :: Basis(:), dBasisdx(:,:)
    REAL(KIND=dp), ALLOCATABLE :: LocalPotential(:)
    LOGICAL :: Found
    TYPE(ValueList_t), POINTER :: Material
    REAL(KIND=dp), POINTER :: Conductivity(:,:,:)=>NULL()
    
    SAVE Conductivity, Nodes
    
    n = MAX( Solver % Mesh % MaxElementDOFs, Solver % Mesh % MaxElementNodes )
    ALLOCATE( STIFF(n,n), FORCE(dim,n) )
    ALLOCATE( LocalPotential(n), Basis(n), dBasisdx(n,3) )

    DO elem = 1,Solver % NumberOFActiveElements
         
      ! Element information
      ! ---------------------
      Element => GetActiveElement(elem)
      CALL GetElementNodes( Nodes )
      nd = GetElementNOFDOFs()
      n  = GetElementNOFNodes()
      
      CALL GetRealArray( GetMaterial(), Conductivity, CondName, Found )
      Rank = 0
      IF ( Found ) Rank = GetTensorRank(Conductivity)
      CALL GetScalarLocalSolution( LocalPotential, VarName )
      
      ! Integrate local stresses:
      ! -------------------------
      IntegStuff = GaussPoints( Element )
      STIFF  = 0.0_dp
      FORCE  = 0.0_dp

      C = 0.0_dp
      DO i=1,dim
        C(i,i) = 1.0_dp
      END DO
      
      DO t=1,IntegStuff % n
        Found = ElementInfo( Element, Nodes, IntegStuff % u(t), &
            IntegStuff % v(t), IntegStuff % w(t), detJ, Basis, dBasisdx )
        
        Weight = IntegStuff % s(t) * detJ
        IF ( CSymmetry ) Weight = Weight * SUM( Basis(1:n) * Nodes % x(1:n) )
        
        IF ( .NOT. ConstantBulkMatrixInUse ) THEN
          DO p=1,nd
            DO q=1,nd
              STIFF(p,q) = STIFF(p,q) + Weight * Basis(q) * Basis(p)
            END DO
          END DO
        END IF
        
        Grad(1:dim) = MATMUL( LocalPotential(1:nd), dBasisdx(1:nd,1:dim) )
        
        SELECT CASE(Rank)
        CASE(0)
        CASE(1)
          DO i=1,dim
            C(i,i) = SUM( Basis(1:n) * Conductivity(1,1,1:n) )
          END DO
        CASE(2)
          DO i=1,dim
            C(i,i) = SUM( Basis(1:n) * Conductivity(i,1,1:n) )
          END DO
        CASE DEFAULT
          DO i=1,dim
            DO j=1,dim
              C(i,j) = SUM( Basis(1:n) * Conductivity(i,j,1:n) )
            END DO
          END DO
        END SELECT
        
        IF(Dofs == 1) THEN
          Coeff = Weight * SUM( C(ActiveDir,1:dim) * Grad(1:dim) )
          FORCE(1,1:nd) = FORCE(1,1:nd) - Coeff*Basis(1:nd)
        ELSE
          DO i=1,dim
            Coeff = Weight * SUM( C(i,1:dim) * Grad(1:dim) )
            FORCE(i,1:nd) = FORCE(i,1:nd) - Coeff*Basis(1:nd)
          END DO
        END IF
      END DO

!------------------------------------------------------------------------------
!      Update global matrices from local matrices 
!------------------------------------------------------------------------------
      IF ( .NOT. ConstantBulkMatrixInUse ) THEN
        IF(Dofs > 1) Solver % Matrix % RHS => SaveRHS
        CALL DefaultUpdateEquations( STIFF, FORCE(1,1:nd), &
            BulkUpdate=ConstantBulkMatrix )
      ELSE
        CALL DefaultUpdateForce( FORCE(1,1:nd) )       
      END IF

      DO i=2,Dofs
        Solver % Matrix % RHS => ForceVector(i-1,:)
        CALL DefaultUpdateForce( FORCE(i,1:nd) )
      END DO
    END DO
    
    DEALLOCATE( LocalPotential, STIFF, FORCE, Basis, dBasisdx )
!------------------------------------------------------------------------------
  END SUBROUTINE BulkAssembly
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
  FUNCTION GetTensorRank( Tensor ) RESULT ( Rank )
!------------------------------------------------------------------------------
    REAL(KIND=dp), POINTER :: Tensor(:,:,:)
    INTEGER :: Rank
    
    IF ( SIZE(Tensor,1) == 1 ) THEN
      Rank = 1
    ELSE IF ( SIZE(Tensor,2) == 1 ) THEN
      Rank = 2
    ELSE
      Rank = 3
    END IF
!-----------------------------------------------------------------------------
  END FUNCTION GetTensorRank
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
END SUBROUTINE FluxSolver
!------------------------------------------------------------------------------

