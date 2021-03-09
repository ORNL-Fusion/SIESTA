!*******************************************************************************
!>  @file siesta_currents.f90
!>  @brief Contains routines to compute the current density from Curl(B).
!
!  Note separating the Doxygen comment block here so detailed decription is
!  found in the Module not the file.
!
!>  Takes the Curl(B) to return the contravariant current density. Contravariant
!>  components of the magnetic field must be converted to covariant components
!>  first.
!*******************************************************************************
      MODULE siesta_currents
      USE v3_utilities, ONLY: assert
      USE stel_kinds
      USE stel_constants
      USE descriptor_mod, ONLY: SIESTA_COMM
      USE shared_data, ONLY: l_update_state, lasym, ste,                       &
                             jsju_ratio_s, jsju_ratio_a 
      USE fourier
      USE nscalingtools, ONLY: nranks, startglobrow, endglobrow, MPI_ERR
      USE metrics, ONLY: tolowerf, tolowerh
      USE quantities

      IMPLICIT NONE

      CONTAINS

!*******************************************************************************
!  UTILITY SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  @brief Compute current density from Curl(B).
!>
!>  Current density K is defined as
!>
!>     K = Curl(B)/mu0                                                       (1)
!>
!>  In curvilinear coordinates this becomes
!>
!>     K^s = 1/sqrt(g)(dB_v/du - dB_u/dv)                                   (2a)
!>
!>     K^u = 1/sqrt(g)(dB_s/dv - dB_v/ds)                                   (2b)
!>
!>     K^v = 1/sqrt(g)(dB_u/ds - dB_s/du)                                   (2c)
!>
!>  where sqrt(g) is the coordinate system Jacobian. Since there are derivatives
!>  with respect to s, these routines take the half grid values as input. The
!>  contravariant current computed will be on the full grid. These routines also
!>  output currents in real space. Note the computed Kmn are [sqrt(g)*K]mn with
!>  the embedded jacobian factor.
!>
!>  @param[in]    bsupsijh B^s on the half grid.
!>  @param[in]    bsupuijh B^u on the half grid.
!>  @param[in]    bsupvijh B^v on the half grid.
!>  @param[inout] ksupsijf K^s on the full grid.
!>  @param[inout] ksupuijf K^u on the full grid.
!>  @param[inout] ksupvijf K^v on the full grid.
!>  @param[in]    lmagen   Logical controlling calc of magnetic energy wb.
!>  @param[in]    lcurr    UNKNOWN
!-------------------------------------------------------------------------------
      SUBROUTINE cv_currents(bsupsijh, bsupuijh, bsupvijh,                     &
                             ksupsijf, ksupuijf, ksupvijf,                     &
                             lmagen, lcurr)
      USE island_params, ONLY: nu_i, rmajor_i, fourier_context
      USE Hessian, ONLY: mupar_norm
      USE descriptor_mod, ONLY: iam
      USE quantities, ONLY: bsq
      USE timer_mod
      USE mpi_inc
      USE descriptor_mod, ONLY: SIESTA_COMM
      USE shared_data, ONLY: siesta_curtor, unit_out

      IMPLICIT NONE

!  Declare Arguments
      REAL (dp), DIMENSION(:,:,:), ALLOCATABLE, INTENT(in)    :: bsupsijh
      REAL (dp), DIMENSION(:,:,:), ALLOCATABLE, INTENT(in)    :: bsupuijh
      REAL (dp), DIMENSION(:,:,:), ALLOCATABLE, INTENT(in)    :: bsupvijh
      REAL (dp), DIMENSION(:,:,:), ALLOCATABLE, INTENT(inout) :: ksupsijf
      REAL (dp), DIMENSION(:,:,:), ALLOCATABLE, INTENT(inout) :: ksupuijf
      REAL (dp), DIMENSION(:,:,:), ALLOCATABLE, INTENT(inout) :: ksupvijf
      LOGICAL, INTENT(in)                                     :: lmagen
      LOGICAL, INTENT(in)                                     :: lcurr

!  Local Variables
      INTEGER                                                 :: istat
      REAL (dp)                                               :: beta
      REAL (dp)                                               :: ton
      REAL (dp)                                               :: toff
      REAL (dp), DIMENSION(6)                                 :: temp(6)
      REAL (dp)                                               :: tmp
      INTEGER                                                 :: nsmin
      INTEGER                                                 :: nsmax
      INTEGER                                                 :: nmin
      INTEGER                                                 :: nmax
      REAL (dp), ALLOCATABLE, DIMENSION(:,:,:)                :: bsubsijh
      REAL (dp), ALLOCATABLE, DIMENSION(:,:,:)                :: bsubuijh
      REAL (dp), ALLOCATABLE, DIMENSION(:,:,:)                :: bsubvijh
      REAL (dp), ALLOCATABLE, DIMENSION(:,:,:)                :: bs_filter
      REAL (dp), ALLOCATABLE, DIMENSION(:,:,:)                :: bu_filter
      REAL (dp), ALLOCATABLE, DIMENSION(:,:,:)                :: bv_filter
      REAL (dp)                                               :: edge_cur

!  Start of executable code
      CALL second0(ton)

      nmin = MAX(1, startglobrow)
      nmax = MIN(ns, endglobrow)
      nsmin = LBOUND(bsupsijh,3)
      nsmax = UBOUND(bsupsijh,3)

!  Calculate covariant BsubXh on half mesh
      ALLOCATE(bsubsijh(ntheta,nzeta,nsmin:nsmax),                             &
               bsubuijh(ntheta,nzeta,nsmin:nsmax),                             &
               bsubvijh(ntheta,nzeta,nsmin:nsmax), stat=istat)
      CALL assert_eq(0, istat, 'Allocation1 failed in cv_currents')

      CALL tolowerh(bsupsijh, bsupuijh, bsupvijh,                              &
                    bsubsijh, bsubuijh, bsubvijh, nsmin, nsmax)

!  Compute parallel damping scaling cofficient <bsupv^2/B.B> ~ k^2
      IF (.not.ALLOCATED(mupar_norm)) THEN
         ALLOCATE(mupar_norm(nsmin:nsmax), stat=istat)
         CALL assert_eq(0, istat,'mupar_norm allocation failed!')
         beta = MAX(1.E-5_dp, wp_i/wb_i)
         mupar_norm = signjac*beta/(rmajor_i*rmajor_i)  !/ vp_f(nsmin:nsmax)     !add 1/jac correction factor in AddDamping
         IF (iam .eq. 0 .and. lverbose) THEN
            WRITE (*,1000) rmajor_i
         END IF
      END IF

      IF (lmagen) THEN
         bsq(:,:,nmin:nmax) = bsupsijh(:,:,nmin:nmax)*bsubsijh(:,:,nmin:nmax)  &
                            + bsupuijh(:,:,nmin:nmax)*bsubuijh(:,:,nmin:nmax)  &
                            + bsupvijh(:,:,nmin:nmax)*bsubvijh(:,:,nmin:nmax)
         IF (nmin .EQ. 1) THEN
            bsq(:,:,1) = 0
         END IF

         CALL assert(nmin                   .eq. 1 .or.                        &
                     ALL(bsq(:,:,nmin:nmax) .gt. zero),                        &
                     'BSQ <= 0 in cv_currents!')

         wb = SUM(bsq(:,:,nmin:nmax)*jacobh(:,:,nmin:nmax)*wint(:,:,nmin:nmax))

         IF (PARSOLVER) THEN
#if defined(MPI_OPT)
!  FIXME: All reduce is not deterministic. This causes a divergent run sequence.
            CALL MPI_ALLREDUCE(MPI_IN_PLACE, wb, 1, MPI_REAL8, MPI_SUM,        &
                               SIESTA_COMM, MPI_ERR)
#endif
         END IF
         wb = signjac*(twopi*pi*hs_i)*wb
      END IF

!  Calculate BsubXmnh and KsupXmnf harmonics. Deallocated in calling routine.
      IF (l_update_state) THEN
         ALLOCATE(bs_filter(ntheta,nzeta,nsmin:nsmax),                         &
                  bu_filter(ntheta,nzeta,nsmin:nsmax),                         &
                  bv_filter(ntheta,nzeta,nsmin:nsmax), stat=istat)
         CALL assert_eq(0, istat, 'Allocation2 failed in CURLB')
      END IF

      edge_cur = 0.0
      CALL curlb(ksupsmnsf, ksupumncf, ksupvmncf,                              &
                 bsubsijh, bsubuijh, bsubvijh,                                 &
                 bs_filter, bu_filter, bv_filter, edge_cur,                    &
                 jsju_ratio_s, f_sin, lcurr, nsmax, nsmin, siesta_curtor)
      IF (lasym) THEN
         CALL curlb(ksupsmncf, ksupumnsf, ksupvmnsf,                           &
                    bsubsijh, bsubuijh, bsubvijh,                              &
                    bs_filter, bu_filter, bv_filter, edge_cur,                 &
                    jsju_ratio_a, f_cos, lcurr, nsmax, nsmin,                  &
                    siesta_curtor)
      END IF

      IF (nsmax .eq. ns) THEN
         siesta_curtor = -fourier_context%orthonorm(m0,n0)*siesta_curtor       &
                       / b_factor
      END IF

      IF (lcurr) THEN
         nsmax = MIN(ns, endglobrow + 1)
      ELSE
         nsmax = MIN(ns, endglobrow)
      END IF

!  Calculate full mesh, contravariant real-space current components
!  KsupXF = sqrt(g)*JsupXF.
      CALL getcv(ksupsmnsf, ksupumncf, ksupvmncf,                              &
                 ksupsijf, ksupuijf, ksupvijf, f_sin, nsmin, nsmax)
      IF (lasym) THEN
         CALL getcv(ksupsmncf, ksupumnsf, ksupvmnsf,                           &
                    ksupsijf, ksupuijf, ksupvijf, f_cos, nsmin, nsmax)
      END IF

!  Covariant K components (with jacob factor still there) in real space needed
!  for current diffusion
      IF (lcurr) THEN
         CALL tolowerf(ksupsijf, ksupuijf, ksupvijf,                           &
                       ksubsijf, ksubuijf, ksubvijf, nsmin, nsmax)
         IF (nsmin .eq. 1) THEN
            ksubuijf(:,:,1) = 0
         END IF
      END IF

!  Diagnostic output (only compute when state is updated)
      DIAGNO: IF (l_update_state) THEN

!  Compute spectral truncation error measure (Use ksupsijf for temporary
!  filtered bsubXijh values)
         IF (nsmax .eq. ns) THEN
            tmp = SUM(wint(:,:,ns)*(bsubvijh(:,:,ns)*bsubvijh(:,:,ns) +        &
                                    bsubuijh(:,:,ns)*bsubuijh(:,:,ns)))
            IF (tmp .ne. zero) THEN
               tmp = SQRT(edge_cur/tmp)
            END IF
         END IF
#if defined(MPI_OPT)
         IF (PARSOLVER) THEN
            CALL MPI_BCAST(tmp, 1, MPI_REAL8, nranks - 1, SIESTA_COMM, MPI_ERR)
         END IF
#endif
         IF (iam .EQ. 0) THEN
            IF (lverbose) THEN
               WRITE (*, 1001) tmp
            END IF
            WRITE (unit_out, 1001) tmp
            CALL FLUSH(unit_out)
         END IF

         temp(1) = SUM(wint(:,:,nmin:nmax) *                                   &
                       (bs_filter - bsubsijh(:,:,nmin:nmax))**2)
         temp(2) = SUM(wint(:,:,nmin:nmax) *                                   &
                       (bs_filter + bsubsijh(:,:,nmin:nmax))**2)
         temp(3) = SUM(wint(:,:,nmin:nmax) *                                   &
                       (bu_filter - bsubuijh(:,:,nmin:nmax))**2)
         temp(4) = SUM(wint(:,:,nmin:nmax) *                                   &
                       (bu_filter + bsubuijh(:,:,nmin:nmax))**2)
         temp(5) = SUM(wint(:,:,nmin:nmax) *                                   &
                       (bv_filter - bsubvijh(:,:,nmin:nmax))**2)
         temp(6) = SUM(wint(:,:,nmin:nmax) *                                   &
                       (bv_filter + bsubvijh(:,:,nmin:nmax))**2)
#if defined(MPI_OPT)
         IF (PARSOLVER) THEN
            CALL MPI_ALLREDUCE(MPI_IN_PLACE, temp, 6, MPI_REAL8, MPI_SUM,      &
                               SIESTA_COMM, MPI_ERR)
         END IF
#endif
         IF (temp(2) .ne. zero) THEN
            ste(2) = SQRT(temp(1)/temp(2))
         END IF
         IF (temp(4) .ne. zero) THEN
            ste(3) = SQRT(temp(3)/temp(4))
         END IF
         IF (temp(6) .ne. zero) THEN
            ste(4) = SQRT(temp(5)/temp(6))
         END IF
         DEALLOCATE (bs_filter, bu_filter, bv_filter)
      END IF DIAGNO

      DEALLOCATE (bsubsijh, bsubuijh, bsubvijh) 

      CALL second0(toff)
      cv_current_time = cv_current_time + (toff - ton)
      time_current = time_current + (toff - ton)

1000  FORMAT(' Rmajor_i = ',1p,e12.3)
1001  FORMAT(' RMS EDGE CURRENT ||KSUPS(NS)|| = ',1p,e10.3)

      END SUBROUTINE

!-------------------------------------------------------------------------------
!>  @brief Compute fourier components of contravariant currents on full radial
!>  grid.
!>
!>  @param[inout] ksupsmnf   Contravarant s component of the current with
!>                           jacobian.
!>  @param[inout] ksupumnf   Contravarant u component of the current with
!>                           jacobian.
!>  @param[inout] ksupvmnf   Contravarant v component of the current with
!>                           jacobian.
!>  @param[in]    bsubsijh   Covariant s component of the magnetic field.
!>  @param[in]    bsubuijh   Covariant u component of the magnetic field.
!>  @param[inout] bsubvijh   Covariant v component of the magnetic field.
!>  @param[inout] bs_filter
!>  @param[inout] bu_filter
!>  @param[inout] bv_filter
!>  @param[inout] edge_cur   Edge current.
!>  @param[out]   jsju_ratio Ratio of <K^s + r12*K^u>/<K^s - r12*K^u>.
!>  @param[in]    iparity    Parity flag for the fouier components.
!>  @param[in]    lcurr      UNKNOWN
!>  @param[in]    nsmax      Maximum radial index.
!>  @param[in]    nsmin      Minimum radial index.
!>  @param[inout] curtor     Current enclosed at the boundary.
!-------------------------------------------------------------------------------
      SUBROUTINE curlb(ksupsmnf, ksupumnf, ksupvmnf,                           &
                       bsubsijh, bsubuijh, bsubvijh,                           &
                       bs_filter, bu_filter, bv_filter, edge_cur,              &
                       jsju_ratio, iparity, lcurr, nsmax, nsmin, curtor)
      USE utilities, ONLY: curl_htof
      USE island_params, ONLY: fourier_context

      IMPLICIT NONE

!  Declare Arguments
      REAL (dp), ALLOCATABLE, DIMENSION(:,:,:), INTENT(inout) :: ksupsmnf
      REAL (dp), ALLOCATABLE, DIMENSION(:,:,:), INTENT(inout) :: ksupumnf
      REAL (dp), ALLOCATABLE, DIMENSION(:,:,:), INTENT(inout) :: ksupvmnf
      REAL (dp), ALLOCATABLE, DIMENSION(:,:,:), INTENT(in)    :: bsubsijh
      REAL (dp), ALLOCATABLE, DIMENSION(:,:,:), INTENT(in)    :: bsubuijh
      REAL (dp), ALLOCATABLE, DIMENSION(:,:,:), INTENT(inout) :: bsubvijh
      REAL (dp), DIMENSION(:,:,:), INTENT(inout)              :: bs_filter
      REAL (dp), DIMENSION(:,:,:), INTENT(inout)              :: bu_filter
      REAL (dp), DIMENSION(:,:,:), INTENT(inout)              :: bv_filter
      REAL (dp), INTENT(inout)                                :: edge_cur
      REAL (dp), INTENT(out)                                  :: jsju_ratio
      INTEGER, INTENT(in)                                     :: iparity
      LOGICAL, INTENT(in)                                     :: lcurr
      INTEGER, INTENT(in)                                     :: nsmax
      INTEGER, INTENT(in)                                     :: nsmin
      REAL (dp), INTENT(inout)                                :: curtor

!  local variables
      INTEGER                                                 :: fours
      INTEGER                                                 :: fouruv
      INTEGER                                                 :: fcomb
      INTEGER                                                 :: nmin
      INTEGER                                                 :: nmax
      INTEGER                                                 :: istat
      REAL (dp)                                               :: r12
      REAL (dp), DIMENSION(-ntor:ntor)                        :: temp
      REAL (dp), ALLOCATABLE, DIMENSION(:,:,:)                :: bsubsmnh
      REAL (dp), ALLOCATABLE, DIMENSION(:,:,:)                :: bsubumnh
      REAL (dp), ALLOCATABLE, DIMENSION(:,:,:)                :: bsubvmnh

!  Start of executable code
      nmin = MAX(1, startglobrow)
      nmax = MIN(ns, endglobrow)

      IF (iparity .EQ. f_sin) THEN
         fcomb = f_none
         fours = f_sin
         fouruv = f_cos
         r12 = 0.5*hs_i
      ELSE
         fcomb = f_sum
         fours = f_cos
         fouruv = f_sin
         r12 = -0.5*hs_i
      END IF

      ALLOCATE(bsubsmnh(0:mpol,-ntor:ntor,nsmin:nsmax),                        &
               bsubumnh(0:mpol,-ntor:ntor,nsmin:nsmax),                        &
               bsubvmnh(0:mpol,-ntor:ntor,nsmin:nsmax), stat=istat)   
      CALL assert_eq(0, istat, 'Allocation1 failed in CURLB')

!  Compute bsubxmnh
      CALL fourier_context%tomnsp(bsubsijh, bsubsmnh, fours)
      CALL fourier_context%tomnsp(bsubuijh, bsubumnh, fouruv)
      CALL fourier_context%tomnsp(bsubvijh, bsubvmnh, fouruv)

      IF (lcurr) THEN
         nmax = MIN(endglobrow + 1, ns)
      ELSE
         nmax = MIN(endglobrow, ns)
      END IF

      CALL CURL_HtoF(bsubsmnh, bsubumnh, bsubvmnh,                             &
                     ksupsmnf, ksupumnf, ksupvmnf,                             &
                     iparity, nsmin, nsmax, nmax, curtor)

!  Diagnostics
      DIAGNO: IF (l_update_state) THEN
         IF (nmax .eq. ns) THEN
            edge_cur = edge_cur + SUM(ksupsmnf(:,:,ns)*ksupsmnf(:,:,ns))
         END IF

         CALL fourier_context%toijsp(bsubsmnh(:,:,nmin:nmax), bs_filter,       &
                                     fcomb, fours)
         CALL fourier_context%toijsp(bsubumnh(:,:,nmin:nmax), bu_filter,       &
                                     fcomb, fouruv)
         CALL fourier_context%toijsp(bsubvmnh(:,:,nmin:nmax), bv_filter,       &
                                     fcomb, fouruv)
	     IF (nmin .eq. 1) THEN
            bv_filter(:,:,1) = bv_filter(:,:,2)
            bsubvijh(:,:,1) = bsubvijh(:,:,2)
         END IF

!  K^s-r12*K^u at origin BEFORE bc applied
         IF (nsmin .eq. 1) THEN
            temp = ksupsmnf(m1,:,1) + r12*ksupumnf(m1,:,1)
            jsju_ratio = SUM(temp*temp)
            IF (jsju_ratio .gt. zero) THEN
               temp = ksupsmnf(m1,:,1) - r12*ksupumnf(m1,:,1)
               jsju_ratio = SUM(temp*temp)/jsju_ratio
               jsju_ratio = SQRT(jsju_ratio)
            END IF
         END IF
      END IF DIAGNO

      DEALLOCATE (bsubsmnh, bsubumnh, bsubvmnh, stat=istat)   
      CALL assert_eq(0, istat, 'Deallocation failed in CURLB')

      END SUBROUTINE

!-------------------------------------------------------------------------------
!>  @brief Transform full grid contravariant currents to real space.
!>
!>  @param[in]    ksupsmnf Contravariant fourier current in the s direction.
!>  @param[in]    ksupumnf Contravariant fourier current in the u direction.
!>  @param[in]    ksupvmnf Contravariant fourier current in the v direction.
!>  @param[inout] ksupsmnf Contravariant real space current in the s direction.
!>  @param[inout] ksupumnf Contravariant real space current in the u direction.
!>  @param[inout] ksupvmnf Contravariant real space current in the v direction.
!>  @param[in]    parity   Partity flag fourier components.
!>  @param[in]    nsmin    Minimum radial index.
!>  @param[in]    nsmax    Maximum radial index.
!-------------------------------------------------------------------------------
      SUBROUTINE getcv(ksupsmnf, ksupumnf, ksupvmnf,                           &
                       ksupsijf, ksupuijf, ksupvijf, parity, nsmin, nsmax)
      USE island_params, ONLY: fourier_context

      IMPLICIT NONE

!  Declare Arguments
      REAL(dp), DIMENSION(:,:,:), ALLOCATABLE, INTENT(in)    :: ksupsmnf
      REAL(dp), DIMENSION(:,:,:), ALLOCATABLE, INTENT(in)    :: ksupumnf
      REAL(dp), DIMENSION(:,:,:), ALLOCATABLE, INTENT(in)    :: ksupvmnf
      REAL(dp), DIMENSION(:,:,:), ALLOCATABLE, INTENT(inout) :: ksupsijf
      REAL(dp), DIMENSION(:,:,:), ALLOCATABLE, INTENT(inout) :: ksupuijf
      REAL(dp), DIMENSION(:,:,:), ALLOCATABLE, INTENT(inout) :: ksupvijf
      INTEGER, INTENT(in)                                    :: parity
      INTEGER, INTENT(in)                                    :: nsmin
      INTEGER, INTENT(in)                                    :: nsmax

!  local variables
      INTEGER                                                :: fours
      INTEGER                                                :: fouruv
      INTEGER                                                :: fcomb

!  Start of executable code
      IF (parity .eq. f_sin) THEN
         fcomb = f_none
         fours = f_sin
         fouruv = f_cos
      ELSE
         fcomb = f_sum
         fours = f_cos
         fouruv = f_sin
      END IF

      CALL fourier_context%toijsp(ksupsmnf(:,:,nsmin:nsmax), ksupsijf, fcomb,  &
                                  fours)
      CALL fourier_context%toijsp(ksupumnf(:,:,nsmin:nsmax), ksupuijf, fcomb,  &
                                  fouruv)
      CALL fourier_context%toijsp(ksupvmnf(:,:,nsmin:nsmax), ksupvijf, fcomb,  &
                                  fouruv)

      END SUBROUTINE

      END MODULE
