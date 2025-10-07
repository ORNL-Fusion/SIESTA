!*******************************************************************************
!>  @file fourier.f90
!>  @brief Contains module @ref fourier
!
!  Note separating the Doxygen comment block here so detailed decription is
!  found in the Module not the file.
!
!>  Module contains subroutines for computing FFT on parallel radial intervals.
!>  Converts quantities between real and fourier space. Note @ref fixarray must
!>  be called once before call any of the Fourier subroutines to calculate all
!>  necessary sine and cosine factors on a fixed mesh of collocation angles.
!*******************************************************************************
      MODULE fourier
      USE stel_kinds
      USE stel_constants
      USE timer_mod
      USE descriptor_mod, ONLY: INHESSIAN, DIAGONALDONE
      
      IMPLICIT NONE

!*******************************************************************************
!  fourier module parameters
!*******************************************************************************
!>  Cosine parity.
      INTEGER, PARAMETER   :: f_cos = 0
!>  Sine parity.
      INTEGER, PARAMETER   :: f_sin = 1

!  Bit locations for fourier control flags.
!>  Do not sum fouier real space transformed quantity. This is used when a
!>  stellarator symmetric parity is being converted to real space.
      INTEGER, PARAMETER   :: f_none = 0               !000000
!>  Radial derivative flag.
      INTEGER, PARAMETER   :: f_ds   = 1               !000001
!>  Poloidal derivative flag.
      INTEGER, PARAMETER   :: f_du   = 2               !000010
!>  Toroidal derivative flag.
      INTEGER, PARAMETER   :: f_dv   = 4               !000100
!>  Combined toroidal and poloidal derivatives.
      INTEGER, PARAMETER   :: f_dudv = IOR(f_du, f_dv) !000110
!>  Sum fouier real space transformed quantity. This is used when a
!>  non-stellarator symmetric parity is being converted to real space.
      INTEGER, PARAMETER   :: f_sum  = 8               !001000
!>  Quantity, contains jacobian.
      INTEGER, PARAMETER   :: f_jac  = 16              !010000
!>  Covariant basis.
      INTEGER, PARAMETER   :: f_con  = 32              !100000

!  Test bit positions for Derivative flags.
!>  Bit position of the @ref f_ds flag.
      INTEGER, PARAMETER   :: b_ds   = 0
!>  Bit position of the @ref f_du flag.
      INTEGER, PARAMETER   :: b_du   = 1
!>  Bit position of the @ref f_dv flag.
      INTEGER, PARAMETER   :: b_dv   = 2
!>  Bit position of the @ref f_sum flag.
      INTEGER, PARAMETER   :: b_sum  = 3
!>  Bit position of the @ref f_jac flag.
      INTEGER, PARAMETER   :: b_jac  = 4
!>  Bit position of the @ref b_con flag.
      INTEGER, PARAMETER   :: b_con  = 5

!  Special Fourier modes.
!  Poloidal modes.
!>  m = 0 mode.
      INTEGER, PARAMETER   :: m0 = 0
!>  m = 1 mode.
      INTEGER, PARAMETER   :: m1 = 1
!>  m = 2 mode.
      INTEGER, PARAMETER   :: m2 = 2
!  Toroidal modes.
!>  n = 0 mode.
      INTEGER, PARAMETER   :: n0 = 0
!>  n = 1 mode.
      INTEGER, PARAMETER   :: n1 = 1

!*******************************************************************************
!  DERIVED-TYPE DECLARATIONS
!  1) fourier base class
!
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  Base class containing fourier memory.
!-------------------------------------------------------------------------------
      TYPE :: fourier_class
!>  Orthonorm factors for normalization.
         REAL (dp), DIMENSION(:,:), POINTER :: orthonorm => null()
!>  Cosine poloidal m and u values.
         REAL (dp), DIMENSION(:,:), POINTER :: cosmu => null()
!>  Cosine derivative with respect to u for poloidal m and u values.
         REAL (dp), DIMENSION(:,:), POINTER :: cosmum => null()
!>  Normalized Cosine poloidal m and u values.
         REAL (dp), DIMENSION(:,:), POINTER :: cosmui => null()
!>  Cosine of nv.
         REAL (dp), DIMENSION(:,:), POINTER :: cosnv => null()
!>  Cosine derivative with respect to v.
         REAL (dp), DIMENSION(:,:), POINTER :: cosnvn => null()
!>  Sine poloidal m and u values.
         REAL (dp), DIMENSION(:,:), POINTER :: sinmu => null()
!>  Sine derivative with respect to u for poloidal m and u values.
         REAL (dp), DIMENSION(:,:), POINTER :: sinmum => null()
!>  Normalized Sine poloidal m and u values.
         REAL (dp), DIMENSION(:,:), POINTER :: sinmui => null()
!>  Sine of nv.
         REAL (dp), DIMENSION(:,:), POINTER :: sinnv => null()
!>  Sine derivative with respect to v.
         REAL (dp), DIMENSION(:,:), POINTER :: sinnvn => null()
!>  Working buffer for first mj terms.
         REAL (dp), DIMENSION(:,:), POINTER :: workmj1 => null()
!>  Working buffer for second mj terms.
         REAL (dp), DIMENSION(:,:), POINTER :: workmj2 => null()
!>  Working buffer for first in terms.
         REAL (dp), DIMENSION(:,:), POINTER :: workin1 => null()
!>  Working buffer for second in terms.
         REAL (dp), DIMENSION(:,:), POINTER :: workin2 => null()
!>  Toroidal mode numbers.
         INTEGER, DIMENSION(:), POINTER     :: tor_modes => null()
!>  Found modes resonances.
         LOGICAL, DIMENSION(:,:), POINTER   :: found_modes => null()
      CONTAINS
         PROCEDURE, PASS :: toijsp_2d => fourier_toijsp_2d
         PROCEDURE, PASS :: toijsp_3d => fourier_toijsp_3d
         GENERIC         :: toijsp => toijsp_2d, toijsp_3d
         PROCEDURE, PASS :: tomnsp_3d => fourier_tomnsp_3d
         PROCEDURE, PASS :: tomnsp_2d => fourier_tomnsp_2d
         PROCEDURE, PASS :: tomnsp_2d_u => fourier_tomnsp_2d_u
         PROCEDURE, PASS :: tomnsp_2d_v => fourier_tomnsp_2d_v
         PROCEDURE, PASS :: tomnsp_3d_pest => fourier_tomnsp_3d_pest
         PROCEDURE, PASS :: tomnsp_2d_pest => fourier_tomnsp_2d_pest
         PROCEDURE, PASS :: tomnsp_2d_u_pest => fourier_tomnsp_2d_u_pest
         GENERIC         :: tomnsp => tomnsp_2d, tomnsp_3d, tomnsp_2d_u,       &
                                      tomnsp_2d_v, tomnsp_3d_pest,             &
                                      tomnsp_2d_pest, tomnsp_2d_u_pest
         PROCEDURE       :: get_origin => fourier_get_origin
         PROCEDURE       :: get_index => fourier_get_index
         FINAL           :: fourier_destruct
      END TYPE

!*******************************************************************************
!  INTERFACE BLOCKS
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  Interface for the fourier constructor.
!-------------------------------------------------------------------------------
      INTERFACE fourier_class
         MODULE PROCEDURE fourier_construct
      END INTERFACE

!-------------------------------------------------------------------------------
!>  Interface for checking results of the unit tests.
!-------------------------------------------------------------------------------
      INTERFACE check
         MODULE PROCEDURE check_all,                                           &
     &                    check_mn
      END INTERFACE

      PRIVATE :: check, check_all, check_mn

      CONTAINS

!*******************************************************************************
!  CONSTRUCTION SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  @brief Construct a @ref fourier_class object.
!>
!>  Allocates memory and initializes a @ref fourier_class object. Computes the
!>  othronorm and cosine and sine buffers.
!>
!>  This subroutine computes the cosine-sine factors that will be needed when
!>  moving between Fourier and real space. All normalizations are contained in
!>  the poloidal quantities used in the Fourier to real space transformation:
!>  * sinmui
!>  * cosmui
!>
!>  Fourier representations are assumed to have stellarator symmetry:
!>  *# cosine (iparity=0):      C(u, v) = sum_u sum_v (C_mn COS(mu + n*nfp*v))
!>  *# sine   (iparity=1):      S(u, v) = sum_u sum_v (S_mn SIN(mu + n*nfp*v))
!>
!>  The number of collocation points have been set initalially equal to the
!>  number of modes:
!>    theta_j = j*pi/M,         j = 0,...,M   Where M = mpol - 1
!>    zeta_k  = k*2*pi/(2N + 1) k = 0,...,2N  Where N = ntor.
!>
!>  @param[in] mpol      Number of poloidal modes.
!>  @param[in] ntor      Number of toroidal modes.
!>  @param[in] ntheta    Number of poloidal real grid points.
!>  @param[in] nzeta     Number of toroidal real grid points.
!>  @param[in] nfp       Number of field periods.
!>  @param[in] sym       Symmetry flag.
!>  @param[in] tor_modes Toroidal mode numbers.
!-------------------------------------------------------------------------------
      FUNCTION fourier_construct(mpol, ntor, ntheta, nzeta, nfp, sym,          &
     &                           tor_modes)

      IMPLICIT NONE

!  Declare Arguments
      CLASS (fourier_class), POINTER             :: fourier_construct
      INTEGER, INTENT(in)                        :: mpol
      INTEGER, INTENT(in)                        :: ntor
      INTEGER, INTENT(in)                        :: ntheta
      INTEGER, INTENT(in)                        :: nzeta
      INTEGER, INTENT(in)                        :: nfp
      LOGICAL, INTENT(in)                        :: sym
      INTEGER, DIMENSION(-ntor:ntor), INTENT(in) :: tor_modes

!  Local Variables
      REAL (dp)                                  :: pi_norm
      INTEGER                                    :: i
      INTEGER                                    :: j
      REAL (dp)                                  :: arg
      REAL (dp)                                  :: dnorm

!  Local Parameters
      REAL (dp), PARAMETER                       :: nfactor = 1.41421356237310_dp

!  Start of executable code.
      ALLOCATE(fourier_construct)

      IF (sym) THEN
         dnorm = one/(ntheta*nzeta)
         pi_norm = twopi/ntheta
      ELSE
         dnorm = one/((ntheta - 1)*nzeta)
         pi_norm = pi/(ntheta - 1)
      END IF

      ALLOCATE(fourier_construct%cosmu(ntheta,0:mpol))
      ALLOCATE(fourier_construct%cosmum(ntheta,0:mpol))
      ALLOCATE(fourier_construct%cosmui(ntheta,0:mpol))
      ALLOCATE(fourier_construct%sinmu(ntheta,0:mpol))
      ALLOCATE(fourier_construct%sinmum(ntheta,0:mpol))
      ALLOCATE(fourier_construct%sinmui(ntheta,0:mpol))

      DO i = 1, ntheta
!  Poloidal collocation points*m.
         arg = pi_norm*(i - 1)

         DO j = 0, mpol
            fourier_construct%cosmu(i, j) = COS(j*arg)
            fourier_construct%sinmu(i, j) = SIN(j*arg)

            CALL fourier_round_trig(fourier_construct%cosmu(i, j),             &
     &                              fourier_construct%sinmu(i, j))

            fourier_construct%cosmum(i, j) = j*fourier_construct%cosmu(i, j)
            fourier_construct%sinmum(i, j) = j*fourier_construct%sinmu(i, j)

            fourier_construct%cosmui(i, j) = dnorm*fourier_construct%cosmu(i, j)
            fourier_construct%sinmui(i, j) = dnorm*fourier_construct%sinmu(i, j)

            IF (.not.sym .and. (i .eq. 1 .or. i .eq. ntheta)) THEN
               fourier_construct%cosmui(i, j) = 0.5                            &
     &                                        * fourier_construct%cosmui(i, j)
               fourier_construct%sinmui(i, j) = 0.5                            &
     &                                        * fourier_construct%sinmui(i, j)
            END IF
         END DO
      END DO

      ALLOCATE(fourier_construct%cosnv(nzeta,0:ntor))
      ALLOCATE(fourier_construct%cosnvn(nzeta,0:ntor))
      ALLOCATE(fourier_construct%sinnv(nzeta,0:ntor))
      ALLOCATE(fourier_construct%sinnvn(nzeta,0:ntor))

!  For now we are in one field period.
      DO j = 0, ntor
         DO i = 1, nzeta
            arg = (twopi*(i - 1))/nzeta

            fourier_construct%cosnv(i,j) = COS(tor_modes(j)*arg)
            fourier_construct%sinnv(i,j) = SIN(tor_modes(j)*arg)

            CALL fourier_round_trig(fourier_construct%cosnv(i,j),              &
     &                              fourier_construct%sinnv(i,j))
         END DO

         fourier_construct%cosnvn(:,j) = tor_modes(j)*nfp*fourier_construct%cosnv(:,j)
         fourier_construct%sinnvn(:,j) = tor_modes(j)*nfp*fourier_construct%sinnv(:,j)
      END DO

!  Compute orthonorm factor for cos(mu + nv), sin(mu + nv) basis
!  (Note cos(mu)cos(nv) basis!) so that
!
!     <cos(mu + nv)*cos(m'u + n'v)>*orthonorm(m,n) = 1 (for m=m',n=n')
!
!  Orthonormal factor is Sqrt(2) for all modes except (0,0) mode. The (0,0) mode
!  is 1.
      ALLOCATE(fourier_construct%orthonorm(0:mpol,-ntor:ntor))
      fourier_construct%orthonorm = nfactor
      fourier_construct%orthonorm(m0,n0) = 1
      IF (ntheta .eq. mpol + 1) THEN
         fourier_construct%orthonorm(mpol,:) = 1
      END IF

!  Allocate working buffers.
      ALLOCATE(fourier_construct%workmj1(0:mpol,nzeta))
      ALLOCATE(fourier_construct%workmj2(0:mpol,nzeta))
      ALLOCATE(fourier_construct%workin1(ntheta,-ntor:ntor))
      ALLOCATE(fourier_construct%workin2(ntheta,-ntor:ntor))
      ALLOCATE(fourier_construct%tor_modes(-ntor:ntor))
      ALLOCATE(fourier_construct%found_modes(0:mpol,-ntor:ntor))

      fourier_construct%tor_modes(-ntor:ntor) = tor_modes(-ntor:ntor)
      fourier_construct%found_modes = .false.

      END FUNCTION

!*******************************************************************************
!  DESTRUCTION SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  @brief Deconstruct a @ref fourier_class object.
!>
!>  Deallocates memory and uninitializes a @ref fourier_class object.
!>
!>  @param[inout] this A @ref fourier_class instance.
!-------------------------------------------------------------------------------
      SUBROUTINE fourier_destruct(this)

      IMPLICIT NONE

!  Declare Arguments
      TYPE (fourier_class), INTENT(inout) :: this

!  Start of executable code.
      IF (ASSOCIATED(this%orthonorm)) THEN
         DEALLOCATE(this%orthonorm)
         this%orthonorm => null()
      END IF

      IF (ASSOCIATED(this%cosmu)) THEN
         DEALLOCATE(this%cosmu)
         this%cosmu => null()
      END IF

      IF (ASSOCIATED(this%cosmum)) THEN
         DEALLOCATE(this%cosmum)
         this%cosmum => null()
      END IF

      IF (ASSOCIATED(this%cosmui)) THEN
         DEALLOCATE(this%cosmui)
         this%cosmui => null()
      END IF

      IF (ASSOCIATED(this%cosnv)) THEN
         DEALLOCATE(this%cosnv)
         this%cosnv => null()
      END IF

      IF (ASSOCIATED(this%cosnvn)) THEN
          DEALLOCATE(this%cosnvn)
         this%cosnvn => null()
      END IF

      IF (ASSOCIATED(this%sinmu)) THEN
         DEALLOCATE(this%sinmu)
         this%sinmu => null()
      END IF

      IF (ASSOCIATED(this%sinmum)) THEN
         DEALLOCATE(this%sinmum)
         this%sinmum => null()
      END IF

      IF (ASSOCIATED(this%sinmui)) THEN
         DEALLOCATE(this%sinmui)
         this%sinmui => null()
      END IF

      IF (ASSOCIATED(this%sinnv)) THEN
         DEALLOCATE(this%sinnv)
         this%sinnv => null()
      END IF

      IF (ASSOCIATED(this%sinnvn)) THEN
         DEALLOCATE(this%sinnvn)
         this%sinnvn => null()
      END IF

      IF (ASSOCIATED(this%workmj1)) THEN
         DEALLOCATE(this%workmj1)
         this%workmj1 => null()
      END IF

      IF (ASSOCIATED(this%workmj2)) THEN
         DEALLOCATE(this%workmj2)
         this%workmj2 => null()
      END IF

      IF (ASSOCIATED(this%workin1)) THEN
         DEALLOCATE(this%workin1)
         this%workin1 => null()
      END IF

      IF (ASSOCIATED(this%workin2)) THEN
         DEALLOCATE(this%workin2)
         this%workin2 => null()
      END IF

      IF (ASSOCIATED(this%tor_modes)) THEN
         DEALLOCATE(this%tor_modes)
         this%tor_modes => null()
      END IF

      IF (ASSOCIATED(this%found_modes)) THEN
         DEALLOCATE(this%found_modes)
         this%found_modes => null()
      END IF

      END SUBROUTINE

!*******************************************************************************
!  UTILITY SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  @brief Convert a quantity from real to fourier space in parallel.
!>
!>  This subroutine moves a 3D quantity to fourier space by summing over its
!>  real represenetation for each surface. @see fourier_tomnsp_2d
!>
!>  @param[inout] this   A @ref fourier_class instance.
!>  @param[in]    xuv    Real space quantity.
!>  @param[out]   xmn    Fourier space quantity.
!>  @param[in]    parity Fourier parity flag.
!-------------------------------------------------------------------------------
      SUBROUTINE fourier_tomnsp_3d(this, xuv, xmn, parity)

      IMPLICIT NONE

!  Declare Arguments
      CLASS (fourier_class), INTENT(inout)     :: this
      REAL (dp), DIMENSION(:,:,:), INTENT(in)  :: xuv
      REAL (dp), DIMENSION(:,:,:), INTENT(out) :: xmn
      INTEGER, INTENT(in)                      :: parity

!  Local Variables
      INTEGER                                  :: js
      REAL (dp)                                :: skston
      REAL (dp)                                :: skstoff

!  Start of executable code.
      CALL second0(skston)

      DO js = 1, MIN(SIZE(xmn,3), SIZE(xuv,3))
         CALL this%tomnsp(xuv(:,:,js), xmn(:,:,js), parity)
      END DO

      CALL second0(skstoff)
      IF (DIAGONALDONE .and. INHESSIAN) THEN
         time_tomnsp = time_tomnsp + (skstoff - skston)
      END IF

      time_tomnsp = time_tomnsp + (skstoff - skston)

      END SUBROUTINE

!-------------------------------------------------------------------------------
!>  @brief Convert a quantity from real to fourier space in pest coordinates.
!>
!>  This subroutine moves a 3D quantity to fourier space by summing over its
!>  real represenetation for each surface. @see fourier_tomnsp_2d_pest.
!>
!>  @param[inout] this   A @ref fourier_class instance.
!>  @param[in]    xuv    Real space quantity.
!>  @param[out]   xmn    Fourier space quantity.
!>  @param[in]    lmns   Stellarator symmetric lambda.
!>  @param[in]    lmnc   Stellarator asymmetric lambda.
!>  @param[in]    parity Fourier parity flag.
!>  @param[in]    asym   Add stellarator asymmetric terms.
!-------------------------------------------------------------------------------
      SUBROUTINE fourier_tomnsp_3d_pest(this, xuv, xmn, lmns, lmnc,            &
                                        parity, asym)

      IMPLICIT NONE

!  Declare Arguments
      CLASS (fourier_class), INTENT(inout)     :: this
      REAL (dp), DIMENSION(:,:,:), INTENT(in)  :: xuv
      REAL (dp), DIMENSION(:,:,:), INTENT(out) :: xmn
      REAL (dp), DIMENSION(:,:,:), INTENT(in)  :: lmns
      REAL (dp), DIMENSION(:,:,:), INTENT(in)  :: lmnc
      INTEGER, INTENT(in)                      :: parity
      LOGICAL, INTENT(in)                      :: asym

!  Local Variables
      INTEGER                                  :: js
      REAL (dp), ALLOCATABLE, DIMENSION(:,:,:) :: lambda
      REAL (dp), ALLOCATABLE, DIMENSION(:,:,:) :: dupdu
      INTEGER                                  :: ns
      REAL (dp)                                :: ton
      REAL (dp)                                :: toff

!  Start of executable code
      CALL second0(ton)

      ns = MIN(SIZE(xmn,3), SIZE(xuv,3), SIZE(lmns, 3))

!  Compute lambda.
      ALLOCATE(lambda(SIZE(xuv,1),SIZE(xuv,2),ns))
      ALLOCATE(dupdu(SIZE(xuv,1),SIZE(xuv,2),ns))

      CALL this%toijsp(lmns, lambda, f_none, f_sin)
      CALL this%toijsp(lmns, dupdu, f_du, f_sin)
      IF (asym) THEN
         CALL this%toijsp(lmnc, lambda, f_sum, f_cos)
         CALL this%toijsp(lmnc, dupdu, IOR(f_sum,f_du), f_cos)
      END IF

!  Transform from up to u
      dupdu = 1 + dupdu

!  Interpolate to the full grid.
      lambda(:,:,1) = lambda(:,:,2)
      dupdu(:,:,1) = dupdu(:,:,2)
      DO js = 2, ns - 1
         lambda(:,:,js) = 0.5*(lambda(:,:,js) + lambda(:,:,js + 1))
         lambda(:,:,js) = 0.5*(lambda(:,:,js) + lambda(:,:,js + 1))
      END DO
      lambda(:,:,ns) = 2*lambda(:,:,ns) - lambda(:,:,ns - 1)
      dupdu(:,:,ns) = 2*dupdu(:,:,ns) - dupdu (:,:,ns - 1)

      DO js = 1, ns
         CALL this%tomnsp(xuv(:,:,js), xmn(:,:,js), lambda(:,:,js),            &
                          dupdu(:,:,js), parity, asym)
      END DO

      DEALLOCATE(lambda)
      DEALLOCATE(dupdu)

      CALL second0(toff)
      IF (DIAGONALDONE.AND.INHESSIAN) THEN
         tomnsp_time = tomnsp_time + (toff - ton)
      END IF

      time_tomnsp = time_tomnsp + (toff - ton)

      END SUBROUTINE

!-------------------------------------------------------------------------------
!>  @brief Fourier transform a 2D quantity to fourier space.
!>
!>  Works by reshaping to a 3D array then passing to the normal 3D version.
!>  Have the 3D function call the 2D function for each slice.
!>
!>  @param[inout] this   A @ref fourier_class instance.
!>  @param[in]    xuv    Fourier space quantity.
!>  @param[out]   xmn    Real space quantity.
!>  @param[in]    parity Fourier parity flag.
!-------------------------------------------------------------------------------
      PURE SUBROUTINE fourier_tomnsp_2d(this, xuv, xmn, parity)

      IMPLICIT NONE

!  Declare Arguments
      CLASS (fourier_class), INTENT(inout)   :: this
      REAL (dp), DIMENSION(:,:), INTENT(in)  :: xuv
      REAL (dp), DIMENSION(:,:), INTENT(out) :: xmn
      INTEGER, INTENT(in)                    :: parity

!  Start of executable code.
!  First, add over poloidal collocation angles. Calls fourier_tomnsp_2d_u.
      CALL this%tomnsp(xuv)
!  Then add over toroidal collocation angles. Calls fourier_tomnsp_2d_v.
      CALL this%tomnsp(xmn, parity)

      END SUBROUTINE

!-------------------------------------------------------------------------------
!>  @brief Fourier transform a 2D quantity to fourier space.
!>
!>  Works by reshaping to a 3D array then passing to the normal 3D version.
!>  Have the 3D function call the 2D function for each slice.
!>
!>  @param[inout] this   A @ref fourier_class instance.
!>  @param[in]    xuv    Fourier space quantity.
!>  @param[out]   xmn    Real space quantity.
!>  @param[in]    lambda Stellarator symmetric lambda.
!>  @param[in]    dupdu  Stellarator asymmetric lambda.
!>  @param[in]    parity Fourier parity flag.
!>  @param[in]    asym   Add stellarator asymmetric terms.
!-------------------------------------------------------------------------------
      PURE SUBROUTINE fourier_tomnsp_2d_pest(this, xuv, xmn, lambda, dupdu,    &
                                             parity, asym)

      IMPLICIT NONE

!  Declare Arguments
      CLASS (fourier_class), INTENT(inout)   :: this
      REAL (dp), DIMENSION(:,:), INTENT(in)  :: xuv
      REAL (dp), DIMENSION(:,:), INTENT(out) :: xmn
      REAL (dp), DIMENSION(:,:), INTENT(in)  :: lambda
      REAL (dp), DIMENSION(:,:), INTENT(in)  :: dupdu
      INTEGER, INTENT(in)                    :: parity
      LOGICAL, INTENT(in)                    :: asym

!  Start of executable code.
!  First, add over poloidal collocation angles. Calls fourier_tomnsp_2d_u.
      CALL this%tomnsp(xuv, lambda, dupdu, asym)
!  Then add over toroidal collocation angles. Calls fourier_tomnsp_2d_v.
      CALL this%tomnsp(xmn, parity)

      END SUBROUTINE

!-------------------------------------------------------------------------------
!>  @brief Sum over the poloidal part.
!>
!>  Pest transforms use the same sum over v but different sums over u. To
!>  encourage code reuse break the u and v sums into separate methods. This
!>  version uses the regular poloidal sum.
!>
!>  @param[inout] this A @ref fourier_class instance.
!>  @param[in]    xuv  Fourier space quantity.
!-------------------------------------------------------------------------------
      PURE SUBROUTINE fourier_tomnsp_2d_u(this, xuv)

      IMPLICIT NONE

!  Declare Arguments
      CLASS (fourier_class), INTENT(inout)  :: this
      REAL (dp), DIMENSION(:,:), INTENT(in) :: xuv

!  Local Variables
      INTEGER                               :: j
      INTEGER                               :: m

!  Start of executable code.

!  First, add over poloidal collocation angles. Note use of cosmui/sinmui
!  include normalization factors.
      DO j = 1, SIZE(xuv, 2) !  nzeta
         DO m = m0, UBOUND(this%cosmu, 2) !  mpol
            this%workmj1(m,j) = DOT_PRODUCT(xuv(:,j), this%cosmui(:,m))
            this%workmj2(m,j) = DOT_PRODUCT(xuv(:,j), this%sinmui(:,m))
         END DO
      END DO

      END SUBROUTINE

!-------------------------------------------------------------------------------
!>  @brief Sum over the poloidal part.
!>
!>  Pest transforms use the same sum over v but different sums over u. To
!>  encourage code reuse break the u and v sums into separate methods. This
!>  version uses the lambda poloidal sum.
!>
!>  @param[inout] this A @ref fourier_class instance.
!>  @param[in]    xuv  Fourier space quantity.
!>  @param[in]    lmns Stellarator symmetric lambda.
!>  @param[in]    lmnc Stellarator asymmetric lambda.
!>  @param[in]    asym Add stellarator asymmetric terms.
!-------------------------------------------------------------------------------
      PURE SUBROUTINE fourier_tomnsp_2d_u_pest(this, xuv, lambda, dupdu, asym)

      IMPLICIT NONE

!  Declare Arguments
      CLASS (fourier_class), INTENT(inout)   :: this
      REAL (dp), DIMENSION(:,:), INTENT(in)  :: xuv
      REAL (dp), DIMENSION(:,:), INTENT(in)  :: lambda
      REAL (dp), DIMENSION(:,:), INTENT(in)  :: dupdu
      LOGICAL, INTENT(in)                    :: asym

!  Local Variables
      INTEGER                                :: i
      INTEGER                                :: j
      INTEGER                                :: m
      REAL (dp)                              :: cosml
      REAL (dp)                              :: sinml
      REAL (dp)                              :: templ
      REAL (dp)                              :: cosmuip
      REAL (dp)                              :: sinmuip

!  Start of executable code.

!  Use recurrence relations below for efficiency.
!
!  cosmuip => cos(m*u') = cosmiu(u,m)*cos(m*lambda) - sinmui(u,m)*sin(m*lambda)
!  sinmuip => sin(m*u') = sinmui(u,m)*cos(m*lambda) + cosmui(u,m)*sin(m*lambda)
!
!  Let cos(m*lambda) == cosml(u,v,m)
!    sin(m*lambda) == sinml(u,v,m)
!  use recurrence formulae to compute these efficiently
!    cosml(m=0) = 1
!    sinml(m=0) = 0
!  compute
!    cosml(m=1)
!    sinml(m=1)
!  then
!    cosml(m+1) = cosml(m)*cosml(1) - sinml(m)*sinml(1)
!    sinml(m+1) = sinml(m)*cosml(1) + cosml(m)*sinml(1)

!  First, add over poloidal collocation angles. Note use of cosmui/sinmui
!  include normalization factors.
      DO j = 1, SIZE(xuv, 2) !  nzeta
         DO m = m0, UBOUND(this%cosmu, 2) !  mpol
            this%workmj1(m,j) = 0
            this%workmj2(m,j) = 0
            DO i = 1, SIZE(xuv, 1) !  ntheta
               IF (m .eq. m0) THEN
                  cosml = dupdu(i,j)
                  sinml = 0
               ELSE IF (m .eq. m1) THEN
                  cosml = COS(lambda(i,j))*dupdu(i,j)
                  sinml = SIN(lambda(i,j))*dupdu(i,j)
               ELSE
                  templ = cosml
                  cosml = cosml*COS(lambda(i,j)) - sinml*SIN(lambda(i,j))
                  sinml = sinml*COS(lambda(i,j)) + templ*SIN(lambda(i,j))
               END IF

               cosmuip = this%cosmui(i,m)*cosml - this%sinmui(i,m)*sinml
               sinmuip = this%sinmui(i,m)*cosml + this%cosmui(i,m)*sinml

               this%workmj1(m,j) = this%workmj1(m,j) + xuv(i,j)*cosmuip
               this%workmj2(m,j) = this%workmj2(m,j) + xuv(i,j)*sinmuip
            END DO
         END DO
      END DO

      END SUBROUTINE

!-------------------------------------------------------------------------------
!>  @brief Sum over toroidal angle.
!>
!>  Pest transforms use the same sum over v but different sums over u. To
!>  encourage code reuse break the u and v sums into separate methods. This
!>  impliments the common toroidal sum.
!>
!>  @param[in]  this   A @ref fourier_class instance.
!>  @param[in]  xuv    Fourier space quantity.
!>  @param[out] xmn    Real space quantity.
!>  @param[in]  parity Fourier parity flag.
!-------------------------------------------------------------------------------
      PURE SUBROUTINE fourier_tomnsp_2d_v(this, xmn, parity)

      IMPLICIT NONE

!  Declare Arguments
      CLASS (fourier_class), INTENT(inout)   :: this
      REAL (dp), DIMENSION(:,:), INTENT(out) :: xmn
      INTEGER, INTENT(in)                    :: parity

!  Local Variables
      REAL (dp)                               :: workcos
      REAL (dp)                               :: worksin
      INTEGER                                 :: j
      INTEGER                                 :: m
      INTEGER                                 :: n
      INTEGER                                 :: moff
      INTEGER                                 :: noff

!  Start of executable code.
      moff = LBOUND(xmn,1)
      noff = LBOUND(xmn,2) + UBOUND(this%cosnv,2)

!  Then add over toroidal collocation angles.
      xmn(:,:) = 0
      IF (parity .eq. f_cos) THEN
         DO j = 1, SIZE(this%cosnv, 1) !  nzeta
            DO m = m0, UBOUND(this%cosmu, 2) ! mpol
               workcos = this%workmj1(m,j)*this%cosnv(j,n0)
               worksin = this%workmj2(m,j)*this%sinnv(j,n0)

               xmn(m + moff,n0 + noff) = xmn(m + moff,n0 + noff)               &
                                       + workcos - worksin
            END DO
         END DO

         DO n = n1, UBOUND(this%cosnv, 2) !  ntor
            DO j = 1, SIZE(this%cosnv, 1) !  nzeta
               workcos = this%workmj1(m0,j)*this%cosnv(j,n)
               worksin = this%workmj2(m0,j)*this%sinnv(j,n)

               xmn(m0 + moff,n + noff) = xmn(m0 + moff,n + noff)               &
                                       + workcos - worksin

               DO m = m1, UBOUND(this%cosmu, 2) ! mpol
                  workcos = this%workmj1(m,j)*this%cosnv(j,n)
                  worksin = this%workmj2(m,j)*this%sinnv(j,n)

                  xmn(m + moff,n + noff) = xmn(m + moff,n + noff)              &
                                         + workcos - worksin
                  xmn(m + moff,-n + noff) = xmn(m + moff,-n + noff)            &
                                          + workcos + worksin
               END DO
            END DO
         END DO
         xmn(m0 + moff,:-n1 + noff) = zero
      ELSE
         DO j = 1, SIZE(this%cosnv, 1) !  nzeta
            DO m = m0, UBOUND(this%cosmu, 2) ! mpol
               workcos = this%workmj2(m,j)*this%cosnv(j,n0)
               worksin = this%workmj1(m,j)*this%sinnv(j,n0)

               xmn(m + moff,n0 + noff) = xmn(m + moff,n0 + noff)               &
                                       + workcos + worksin
            END DO
         END DO

         DO n = n1, UBOUND(this%cosnv, 2) !  ntor
            DO j = 1, SIZE(this%cosnv, 1) !  nzeta
               workcos = this%workmj2(m0,j)*this%cosnv(j,n)
               worksin = this%workmj1(m0,j)*this%sinnv(j,n)

               xmn(m0 + moff,n + noff) = xmn(m0 + moff,n + noff)               &
                                       + workcos + worksin

               DO m = m1, UBOUND(this%cosmu, 2) ! mpol
                  workcos = this%workmj2(m,j)*this%cosnv(j,n)
                  worksin = this%workmj1(m,j)*this%sinnv(j,n)

                  xmn(m + moff,n + noff) = xmn(m + moff,n + noff)              &
                                         + workcos + worksin
                  xmn(m + moff,-n + noff) = xmn(m + moff,-n + noff)            &
                                          + workcos - worksin
               END DO
            END DO
         END DO
         xmn(m0 + moff,:n0 + noff) = zero
      END IF

      xmn = this%orthonorm*xmn

      END SUBROUTINE

!-------------------------------------------------------------------------------
!>  @brief Convert a quantity from fourier to real space in parallel.
!>
!>  This subroutine moves a quantity to real space by summing over its Fourier
!>  harmonics. If the @ref f_du flag is set take the poloidal derivative. If the
!>  @ref f_dv flag is set take the toroidal derivative. If the f_sum flag is set
!>  add the real quanity to the previous value.
!>
!>  @param[in]    this    A @ref fourier_class instance.
!>  @param[in]    xmn     Fourier space quantity.
!>  @param[inout] xuv     Real space quantity.
!>  @param[in]    dflag   Derivative and sum control flags.
!>  @param[in]    parity Fourier parity flag.
!-------------------------------------------------------------------------------
      SUBROUTINE fourier_toijsp_3d(this, xmn, xuv, dflag, parity)

      IMPLICIT NONE

!  Declare Arguments
      CLASS (fourier_class), INTENT(inout)       :: this
      REAL (dp), DIMENSION(:,:,:), INTENT(in)    :: xmn
      REAL (dp), DIMENSION(:,:,:), INTENT(inout) :: xuv
      INTEGER, INTENT(in)                        :: dflag
      INTEGER, INTENT(in)                        :: parity

!  Local Variables
      INTEGER                                    :: js
      REAL (dp)                                  :: skston
      REAL (dp)                                  :: skstoff

!  Start of executable code.
      CALL second0(skston)

      DO js = 1, MIN(SIZE(xmn,3), SIZE(xuv,3))
         CALL this%toijsp(xmn(:,:,js), xuv(:,:,js), dflag, parity)
      END DO

      CALL second0(skstoff)
      IF (DIAGONALDONE .AND. INHESSIAN) THEN
         toijsp_time = toijsp_time + (skstoff - skston)
      END IF

      time_toijsp = time_toijsp + (skstoff - skston)

      END SUBROUTINE

!-------------------------------------------------------------------------------
!>  @brief Fourier transform a 2D quantity to real space.
!>
!>  Works by reshaping to a 3D array then passing to the normal 3D version.
!>
!>  @param[in]    this   A @ref fourier_class instance.
!>  @param[in]    xmn    Fourier space quantity.
!>  @param[inout] xuv    Real space quantity.
!>  @param[in]    dflag  Derivative and sum control flags.
!>  @param[in]    parity Fourier parity flag.
!-------------------------------------------------------------------------------
      PURE SUBROUTINE fourier_toijsp_2d(this, xmn, xuv, dflag, parity)

      IMPLICIT NONE

!  Declare Arguments
      CLASS (fourier_class), INTENT(inout)     :: this
      REAL (dp), DIMENSION(:,:), INTENT(in)    :: xmn
      REAL (dp), DIMENSION(:,:), INTENT(inout) :: xuv
      INTEGER, INTENT(in)                      :: parity
      INTEGER, INTENT(in)                      :: dflag

!  Local Variables
      INTEGER                                  :: j
      INTEGER                                  :: m
      INTEGER                                  :: n
      INTEGER                                  :: moff
      INTEGER                                  :: noff
      REAL (dp)                                :: workmn
      REAL (dp), DIMENSION(:,:), POINTER       :: anglem1
      REAL (dp), DIMENSION(:,:), POINTER       :: anglem2
      INTEGER                                  :: msign
      REAL (dp), DIMENSION(:,:), POINTER       :: anglen1
      REAL (dp), DIMENSION(:,:), POINTER       :: anglen2
      INTEGER                                  :: nsign1
      INTEGER                                  :: nsign2
      INTEGER                                  :: nsign3

!  Start of executable code.
      moff = LBOUND(xmn,1)
      noff = LBOUND(xmn,2) + UBOUND(this%cosnv,2)

!  Poloidal derivative requested use cosmum, sinmum.
      IF (BTEST(dflag, b_du)) THEN
         anglem1 => this%sinmum
         anglem2 => this%cosmum
         msign = -1
      ELSE
         anglem1 => this%cosmu
         anglem2 => this%sinmu
         msign = 1
      END IF

      this%workin1 = 0
      this%workin2 = 0

!  Sum over poloidal angles.
      DO n = LBOUND(this%workin1, 2), UBOUND(this%workin1, 2) !  -ntor:ntor
         DO m = 0, UBOUND(this%cosmu, 2) !  mpol
            workmn = this%orthonorm(m,n)*xmn(m + moff,n + noff)

            this%workin1(:,n) = this%workin1(:,n) + msign*workmn*anglem1(:,m)
            this%workin2(:,n) = this%workin2(:,n) + workmn*anglem2(:,m)
         END DO
      END DO

      IF (parity .eq. f_cos) THEN
         nsign3 = 1
         IF (BTEST(dflag, b_dv)) THEN
            anglen1 => this%sinnvn
            anglen2 => this%cosnvn
            nsign1 = -1
            nsign2 = -1
         ELSE
            anglen1 => this%cosnv
            anglen2 => this%sinnv
            nsign1 = 1
            nsign2 = -1
         END IF
      ELSE
         nsign3 = -1
         IF (BTEST(dflag, b_dv)) THEN
            anglen1 => this%cosnvn
            anglen2 => this%sinnvn
            nsign1 = 1
            nsign2 = -1
         ELSE
            anglen1 => this%sinnv
            anglen2 => this%cosnv
            nsign1 = 1
            nsign2 = 1
         END IF
      END IF

      IF (.not.BTEST(dflag, b_sum)) THEN
         xuv(:,:) = 0
      END IF

!  First do the n = 0 mode.
      IF (.not.BTEST(dflag, b_dv)) THEN
         IF (parity .eq. f_cos) THEN
            DO j = 1, SIZE(xuv, 2) !  nzeta
               xuv(:,j) = this%workin1(:,n0) + xuv(:,j)
            END DO
         ELSE
            DO j = 1, SIZE(xuv, 2) !  nzeta
               xuv(:,j) = this%workin2(:,n0) + xuv(:,j)
            END DO
         END IF
      END IF

!  Sum over remaining toroidal modes.
      DO n = n1, UBOUND(this%cosnv, 2) !  ntor
         DO j = 1, SIZE(xuv, 2) !  nzeta
            xuv(:,j) = xuv(:,j)                                                &
                     + nsign1*(this%workin1(:,n) +                             &
                               nsign3*this%workin1(:,-n))*anglen1(j,n)         &
                     + nsign2*(this%workin2(:,n) -                             &
                               nsign3*this%workin2(:,-n))*anglen2(j,n)
         END DO
      END DO

      END SUBROUTINE

!-------------------------------------------------------------------------------
!>  @brief Computes the origin value of a half-mesh quantity.
!>
!>  The origin quantity is found from the m mode of the js = 2 value.
!>
!>  @param[inout] this A @ref fourier_class instance.
!>  @param[inout] xuv  Real space quantity.
!>  @param[in]    mode Poloidal mode number to be retained.
!>  @param[in]    asym Add stellarator asymmetric terms.
!-------------------------------------------------------------------------------
      SUBROUTINE fourier_get_origin(this, xuv, mode, asym)

      IMPLICIT NONE

!  Declare Arguments
      CLASS (fourier_class), INTENT(inout)                   :: this
      REAL(dp), DIMENSION(:,:,:), ALLOCATABLE, INTENT(inout) :: xuv
      INTEGER, INTENT(in)                                    :: mode
      LOGICAL, INTENT(in)                                    :: asym

!  Local Variables
      REAL(dp), DIMENSION(:,:), ALLOCATABLE                  :: xmn_cos
      REAL(dp), DIMENSION(:,:), ALLOCATABLE                  :: xmn_sin

!  Local Parameters
      INTEGER, PARAMETER                                     :: js1 = 1
      INTEGER, PARAMETER                                     :: js2 = 2

!  Start of executable code.

!  This routine is only called once from init_quantities so just allocate this
!  buffer in the call and not the fourier_class.
      ALLOCATE(xmn_cos(0:UBOUND(this%cosmu,2),                                 &
                       -UBOUND(this%cosnv,2):UBOUND(this%cosnv,2)))
      CALL this%tomnsp(xuv(:,:,js2), xmn_cos(:,:), f_cos)

      IF (mode .gt. 0) THEN
         xmn_cos(0:mode - 1,:) = 0
      END IF
      IF (mode .lt. UBOUND(this%cosmu,2)) THEN
         xmn_cos(mode + 1,:) = 0
      END IF

      IF (asym) THEN
         ALLOCATE(xmn_sin(0:UBOUND(this%cosmu,2),                              &
                          -UBOUND(this%cosnv,2):UBOUND(this%cosnv,2)))
         CALL this%tomnsp(xuv(:,:,js2), xmn_sin(:,:), f_sin)

         IF (mode .gt. 0) THEN
            xmn_sin(0:mode - 1,:) = 0
         END IF
         IF (mode .lt. UBOUND(this%cosmu,2)) THEN
            xmn_sin(mode + 1,:) = 0
         END IF
      END IF

      CALL this%toijsp(xmn_cos(:,:), xuv(:,:,js1), f_none, f_cos)
      DEALLOCATE(xmn_cos)

      IF (asym) THEN
         CALL this%toijsp(xmn_sin(:,:), xuv(:,:,js1), f_sum, f_sin)
         DEALLOCATE(xmn_sin)
      END IF

      END SUBROUTINE

!-------------------------------------------------------------------------------
!>  @brief Round trig values to whole number for near numbers.
!>
!>  @param[inout] cosi Cosine value to round.
!>  @param[inout] sini Sine value to round.
!-------------------------------------------------------------------------------
      SUBROUTINE fourier_round_trig(cosi, sini)

      IMPLICIT NONE

!  Declare Arguments
      REAL (dp), INTENT(inout) :: cosi
      REAL (dp), INTENT(inout) :: sini

!  Local variables
      REAL(dp)                 :: eps

!  Start of executable code
      eps = 10*EPSILON(eps)

      IF (ABS(cosi - 1)      .le. eps) THEN
         cosi = 1
         sini = 0
      ELSE IF (ABS(cosi + 1) .le. eps) THEN
         cosi = -1
         sini = 0
      ELSE IF (ABS(sini - 1) .le. eps) THEN
         cosi = 0
         sini = 1
      ELSE IF (ABS(sini + 1) .le. eps) THEN
         cosi = 0
         sini = -1
      END IF

      END SUBROUTINE

!-------------------------------------------------------------------------------
!>  @brief The the index position of the toroidal mode if it exists.
!>
!>  @param[inout] this A @ref fourier_class instance.
!>  @param[inout] n    Toroidal mode number to check the index on exit.
!>  @returns True if the mode number was found.
!-------------------------------------------------------------------------------
      FUNCTION fourier_get_index(this, n)

      IMPLICIT NONE

!  Declare Arguments
      LOGICAL                              :: fourier_get_index
      CLASS (fourier_class), INTENT(inout) :: this
      INTEGER, INTENT(inout)               :: n

!  Local variables
      INTEGER                              :: i

!  Start of executable code
      fourier_get_index = .false.

      DO i = 0, SIGN(UBOUND(this%tor_modes, 1), n), SIGN(1, n)
         IF (this%tor_modes(i) .eq. n) THEN
            n = i
            fourier_get_index = .true.
            EXIT
         END IF
      END DO

      END FUNCTION

!*******************************************************************************
!  UNIT TESTS
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  @brief Test fourier routines.
!>
!>  This tests the routines by fourier transforming form flux to real space then
!>  back again.
!-------------------------------------------------------------------------------
      FUNCTION test_fourier()
      USE timer_mod, ONLY: test_fourier_time

      IMPLICIT NONE

!  Declare Arguments
      LOGICAL                               :: test_fourier

!  Local variables
      INTEGER                               :: m
      INTEGER                               :: n
      REAL(dp), ALLOCATABLE, DIMENSION(:,:) :: testij
      REAL(dp), ALLOCATABLE, DIMENSION(:,:) :: testmn
      REAL(dp), ALLOCATABLE, DIMENSION(:,:) :: result
      CLASS (fourier_class), POINTER        :: four => null()
      INTEGER, ALLOCATABLE, DIMENSION(:)    :: tor_modes

!  Local PARAMETERS
      INTEGER, PARAMETER                    :: pol = 10
      INTEGER, PARAMETER                    :: tor = 10
      INTEGER, PARAMETER                    :: theta = 33
      INTEGER, PARAMETER                    :: zeta = 34
      INTEGER, PARAMETER                    :: nfp = 5

!  Start of executable code.
      ALLOCATE(testij(theta,zeta))
      ALLOCATE(testmn(0:pol,-tor:tor))
      ALLOCATE(result(0:pol,-tor:tor))
      ALLOCATE(tor_modes(-tor:tor))

      DO n = -tor, tor
         tor_modes(n) = n
      END DO

      four => fourier_class(pol, tor, theta, zeta, nfp, .false., tor_modes)
      test_fourier = .true.

      DO n = -tor, tor
         DO m = 0, pol
            IF (m .eq. m0 .and. n .lt. n0) THEN
               CYCLE
            END IF

!  Test cosine parity transform.
            testmn = 0
            testmn(m,n) = 1
            CALL four%toijsp(testmn, testij, f_none, f_cos)
            CALL four%tomnsp(testij, result, f_cos)
            test_fourier = test_fourier .and.                                  &
                           check(testmn, result, pol, tor,                     &
     &                           tor_modes, 'cosine_parity_test')

!  Test sine parity transform.
            IF (m .eq. m0 .and. n .eq. n0) THEN
               CYCLE
            END IF
            testmn = 0
            testmn(m,n) = 1
            CALL four%toijsp(testmn, testij, f_none, f_sin)
            CALL four%tomnsp(testij, result, f_sin)
            test_fourier = test_fourier .and.                                  &
                           check(testmn, result, pol, tor,                     &
     &                           tor_modes, 'sine_parity_test')

!  Test u derivatives of cosine parity.
            testmn = 0
            testmn(m,n) = 1
            CALL four%toijsp(testmn, testij, f_du, f_cos)
            CALL four%tomnsp(testij, result, f_sin)
            testmn(m,n) = -m
            test_fourier = test_fourier .and.                                  &
                           check(testmn, result, pol, tor,                     &
     &                           tor_modes, 'du_cosine_parity_test')

!  Test v derivatives of cosine parity.
            testmn = 0
            testmn(m,n) = 1
            CALL four%toijsp(testmn, testij, f_dv, f_cos)
            CALL four%tomnsp(testij, result, f_sin)
            testmn(m,n) = -n*nfp
            test_fourier = test_fourier .and.                                  &
                           check(testmn, result, pol, tor,                     &
     &                           tor_modes, 'dv_cosine_parity_test')

!  Test u derivatives of sine parity.
            testmn = 0
            testmn(m,n) = 1
            CALL four%toijsp(testmn, testij, f_du, f_sin)
            CALL four%tomnsp(testij, result, f_cos)
            testmn(m,n) = m
            test_fourier = test_fourier .and.                                  &
                           check(testmn, result, pol, tor,                     &
     &                           tor_modes, 'du_sine_parity_test')

!  Test v derivatives of sine parity.
            testmn = 0
            testmn(m,n) = 1
            CALL four%toijsp(testmn, testij, f_dv, f_sin)
            CALL four%tomnsp(testij, result, f_cos)
            testmn(m,n) = n*nfp
            test_fourier = test_fourier .and.                                  &
     &                     check(testmn, result, pol, tor,                     &
     &                           tor_modes, 'dv_sine_parity_test')
         END DO
      END DO

      DEALLOCATE(four)

      DEALLOCATE(testij)
      DEALLOCATE(testmn)
      DEALLOCATE(result)

      DEALLOCATE(tor_modes)

      END FUNCTION

!*******************************************************************************
!  CHECK FUNCTIONS
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  @brief Check test values.
!>
!>  @param[in] expected Expected value for the test.
!>  @param[in] received Recieved value for the test.
!>  @param[in] m        Poloidal mode to check.
!>  @param[in] n        Toroidal mode to check.
!>  @param[in] name     Name of the test.
!-------------------------------------------------------------------------------
      FUNCTION check_mn(expected, received, m, n, name)

      IMPLICIT NONE

!  Declare Arguments
      LOGICAL                       :: check_mn
      REAL (rprec), INTENT(in)      :: expected
      REAL (rprec), INTENT(in)      :: received
      INTEGER, INTENT(in)           :: m
      INTEGER, INTENT(in)           :: n
      CHARACTER (len=*), INTENT(in) :: name

!  Local Parameters
      REAL (rprec), PARAMETER       :: eps = 1.E-12_dp

!  Start of executable code.
      check_mn = ABS(expected - received) .lt. eps
      IF (.not.check_mn) THEN
         WRITE (*,1000) TRIM(name), m, n, expected, received
      END IF

1000  FORMAT('fourier.f90: 'a,' m = ',i3,' n = ',i3,' test failed.',/,         &
             'Expected ',e12.3,' Recieved ',e12.3)

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Check all test values.
!>
!>  @param[in] expected  Expected value for the test.
!>  @param[in] received  Recieved value for the test.
!>  @param[in] mpol      Number of poloidal modes.
!>  @param[in] ntor      Number of toroidal modes.
!>  @param[in] tor_modes Toroidal modes.
!>  @param[in] name      Name of the test.
!-------------------------------------------------------------------------------
      FUNCTION check_all(expected, received, mpol, ntor, tor_modes, name)

      IMPLICIT NONE

!  Declare Arguments
      LOGICAL                                    :: check_all
      REAL (rprec), DIMENSION(:,:), INTENT(in)   :: expected
      REAL (rprec), DIMENSION(:,:), INTENT(in)   :: received
      INTEGER, INTENT(in)                        :: mpol
      INTEGER, INTENT(in)                        :: ntor
      INTEGER, DIMENSION(-ntor:ntor), INTENT(in) :: tor_modes
      CHARACTER (len=*), INTENT(in)              :: name

!  Local Variables
      INTEGER                                    :: m
      INTEGER                                    :: n
      INTEGER                                    :: n_mode
      INTEGER                                    :: moff
      INTEGER                                    :: noff

!  Start of executable code.
      moff = LBOUND(expected, 1)
      noff = LBOUND(expected, 2) + ntor

      check_all = .true.

      
      DO n = -ntor, ntor
         DO m = 0, mpol
            check_all = check_all .and. check(expected(m + moff, n + noff),    &
                                              received(m + moff, n + noff),    &
                                              m, tor_modes(n), name)
         END DO
      END DO

      END FUNCTION

      END MODULE
