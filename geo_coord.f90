!-------------------------------------------------------------------------------
!> Computation of longitude lambda and latitude phi for position (x,y) in the
!! numerical domain.  Returns only phi.
!<------------------------------------------------------------------------------
FUNCTION geo_coord_phi(Model, Node, var) RESULT(phii)

USE Types  
USE CoordinateSystems  
USE SolverUtils  
USE ElementDescription

IMPLICIT NONE

TYPE(Model_t) :: Model
INTEGER :: Node

REAL(KIND=dp):: x_val, y_val

REAL(KIND=dp):: lambda_val, phi_val

REAL(KIND=dp) :: K, phii, var

REAL(KIND=dp), PARAMETER :: phi0 = 1.239184d+00       ! Standard parallel +71 deg (71 deg N) in rad
REAL(KIND=dp), PARAMETER :: lambda0 = -0.6806784d+00  ! Reference longitude -39 deg (39 deg W) in rad

REAL(KIND=dp), PARAMETER :: A = 6.378137d+06          ! Semi-major axis of the Earth = 6378137 m
REAL(KIND=dp), PARAMETER :: B = 6.3567523142d+06      ! Semi-minor axis of the Earth = 6356752.3142 m  

INTEGER, PARAMETER :: GRID = 1

x_val = Model % Nodes % x (Node)
y_val = Model % Nodes % y (Node)


IF (GRID == 0 .OR.  GRID == 1) THEN

!-------- Inverse stereographic projection --------

call stereo_inv_ellipsoid(x_val, y_val, A, B, &
                          lambda0, phi0, lambda_val, phi_val)

phii = phi_val

ELSE IF (GRID == 2) THEN

!-------- Identity --------

lambda_val = x_val
phi_val    = y_val

phii = phi_val


END IF

CONTAINS

!-------------------------------------------------------------------------------
!> Inverse stereographic projection for an ellipsoidal planet.
!<------------------------------------------------------------------------------
SUBROUTINE stereo_inv_ellipsoid(x_val, y_val, A, B, &
                                lambda0, phi0, lambda_val, phi_val)

IMPLICIT NONE

REAL(kind=dp), INTENT(IN)  :: x_val, y_val, A, B, lambda0, phi0
REAL(KIND=dp), INTENT(OUT) :: lambda_val, phi_val

INTEGER :: l
INTEGER :: sign_phi0
REAL(KIND=dp) :: phi_aux, phi0_aux
REAL(KIND=dp) :: e, mc, t, tc, kp, rho, phi_p, residual
REAL(KIND=dp) :: sinphi0, sinlambda0, cosphi0, coslambda0

REAL(KIND=dp), PARAMETER :: eps = 1.0e-05_dp, &
                         eps_residual = 1.0e-09_dp

IF (phi0 > eps) THEN   ! for northern hemisphere
   sign_phi0 =  1
ELSE IF (phi0 < (-eps)) THEN  ! for southern hemisphere
   sign_phi0 = -1
ELSE
   stop ' geo_coord: phi0 must be different from zero!'
END IF

phi0_aux = phi0 * sign_phi0

e=SQRT((A**2-B**2)/(A**2))

sinphi0    = SIN(phi0_aux)
sinlambda0 = SIN(lambda0)
cosphi0    = COS(phi0_aux)
coslambda0 = COS(lambda0)
  
tc=SQRT(((1.0_dp-sinphi0)/(1.0_dp+sinphi0))* &
       ((1.0_dp+e*sinphi0)/(1.0_dp-e*sinphi0))**e)
mc=cosphi0/SQRT(1.0_dp-e*e*sinphi0*sinphi0)
rho=SQRT(x_val*x_val+y_val*y_val)
t=rho*tc/(A*mc)

lambda_val = lambda0 + sign_phi0*ATAN2(y_val,x_val) + 0.5_dp*pi

!  fix point iteration

phi_p=0.5_dp*pi-2.0_dp*ATAN(t)
l=0
residual=3600.0_dp
DO WHILE(residual >= eps_residual)
   l=l+1
   phi_aux=0.5_dp*pi-2.0_dp*ATAN(t*((1.0_dp-e*SIN(phi_p))/ &
           (1.0_dp+e*SIN(phi_p)))**(0.5_dp*e))
   residual=ABS(phi_aux-phi_p)
   phi_p=phi_aux
END DO

phi_val = phi_aux * sign_phi0

IF (lambda_val < 0.0_dp) THEN
   lambda_val = lambda_val + 2.0_dp*pi
ELSE IF (lambda_val >= (2.0_dp*pi)) THEN
   lambda_val = lambda_val - 2.0_dp*pi
END IF

END SUBROUTINE stereo_inv_ellipsoid


END FUNCTION geo_coord_phi






!-------------------------------------------------------------------------------
!> Computation of longitude lambda and latitude phi for position (x,y) in the
!! numerical domain. Returns only lambda.
!<------------------------------------------------------------------------------
FUNCTION geo_coord_lambda(Model, Node, var) RESULT(lambda)

USE Types  
USE CoordinateSystems  
USE SolverUtils  
USE ElementDescription

IMPLICIT NONE

TYPE(Model_t) :: Model
INTEGER :: Node

REAL(KIND=dp):: x_val, y_val

REAL(KIND=dp):: lambda_val, phi_val

REAL(KIND=dp) :: K, lambda, var

REAL(KIND=dp), PARAMETER :: phi0 = 1.239184d+00       ! Standard parallel +71 deg (71 deg N) in rad
REAL(KIND=dp), PARAMETER :: lambda0 = -0.6806784d+00  ! Reference longitude -39 deg (39 deg W) in rad

REAL(KIND=dp), PARAMETER :: A = 6.378137d+06          ! Semi-major axis of the Earth = 6378137 m
REAL(KIND=dp), PARAMETER :: B = 6.3567523142d+06      ! Semi-minor axis of the Earth = 6356752.3142 m  

INTEGER, PARAMETER :: GRID = 1

x_val = Model % Nodes % x (Node)
y_val = Model % Nodes % y (Node)

IF (GRID == 0 .OR.  GRID == 1) THEN

!-------- Inverse stereographic projection --------

call stereo_inv_ellipsoid(x_val, y_val, A, B, &
                          lambda0, phi0, lambda_val, phi_val)

lambda = lambda_val

ELSE IF (GRID == 2) THEN

!-------- Identity --------

lambda_val = x_val
phi_val    = y_val

lambda = lambda_val


END IF

CONTAINS

!-------------------------------------------------------------------------------
!> Inverse stereographic projection for an ellipsoidal planet.
!<------------------------------------------------------------------------------
SUBROUTINE stereo_inv_ellipsoid(x_val, y_val, A, B, &
                                lambda0, phi0, lambda_val, phi_val)

IMPLICIT NONE

REAL(kind=dp), INTENT(IN)  :: x_val, y_val, A, B, lambda0, phi0
REAL(KIND=dp), INTENT(OUT) :: lambda_val, phi_val

INTEGER :: l
INTEGER :: sign_phi0
REAL(KIND=dp) :: phi_aux, phi0_aux
REAL(KIND=dp) :: e, mc, t, tc, kp, rho, phi_p, residual
REAL(KIND=dp) :: sinphi0, sinlambda0, cosphi0, coslambda0

REAL(KIND=dp), PARAMETER :: eps = 1.0e-05_dp, &
                         eps_residual = 1.0e-09_dp

IF (phi0 > eps) THEN   ! for northern hemisphere
   sign_phi0 =  1
ELSE IF (phi0 < (-eps)) THEN  ! for southern hemisphere
   sign_phi0 = -1
ELSE
   stop ' geo_coord: phi0 must be different from zero!'
END IF

phi0_aux = phi0 * sign_phi0

e=SQRT((A**2-B**2)/(A**2))

sinphi0    = SIN(phi0_aux)
sinlambda0 = SIN(lambda0)
cosphi0    = COS(phi0_aux)
coslambda0 = COS(lambda0)
  
tc=SQRT(((1.0_dp-sinphi0)/(1.0_dp+sinphi0))* &
       ((1.0_dp+e*sinphi0)/(1.0_dp-e*sinphi0))**e)
mc=cosphi0/SQRT(1.0_dp-e*e*sinphi0*sinphi0)
rho=SQRT(x_val*x_val+y_val*y_val)
t=rho*tc/(A*mc)

lambda_val = lambda0 + sign_phi0*ATAN2(y_val,x_val) + 0.5_dp*pi

!  fix point iteration

phi_p=0.5_dp*pi-2.0_dp*ATAN(t)
l=0
residual=3600.0_dp
DO WHILE(residual >= eps_residual)
   l=l+1
   phi_aux=0.5_dp*pi-2.0_dp*ATAN(t*((1.0_dp-e*SIN(phi_p))/ &
           (1.0_dp+e*SIN(phi_p)))**(0.5_dp*e))
   residual=ABS(phi_aux-phi_p)
   phi_p=phi_aux
END DO

phi_val = phi_aux * sign_phi0

IF (lambda_val < 0.0_dp) THEN
   lambda_val = lambda_val + 2.0_dp*pi
ELSE IF (lambda_val >= (2.0_dp*pi)) THEN
   lambda_val = lambda_val - 2.0_dp*pi
END IF

END SUBROUTINE stereo_inv_ellipsoid


END FUNCTION geo_coord_lambda
