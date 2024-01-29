!*******************************************************************************
!>  @file metrics.f90
!>  @brief Contains module @ref metrics.
!
!  Note separating the Doxygen comment block here so detailed decription is
!  found in the Module not the file.
!
!>  Module for reading in VMEC input and computing metric elements on the SIESTA
!>  (sqrt-flux) radial grid. Computes and stores the real-space metric elements,
!>  Jacobian based on square root of fluc splined vmec coordinate system.
!*******************************************************************************

      MODULE metrics
      USE stel_kinds
      USE stel_constants
      USE island_params
      USE shared_data, ONLY: lasym, r1_i, z1_i, ru_i, zu_i, rv_i, zv_i,  &
                             lverbose, jsupvdotA, nsp
      USE read_wout_mod, ns_vmec=>ns, mpol_vmec=>mpol, ntor_vmec=>ntor,  &
          rmnc_vmec=>rmnc, zmns_vmec=>zmns, lmns_vmec=>lmns,             &
          xm_vmec=>xm, xn_vmec=>xn, chipf_vmec=>chipf,                   &  ! MRC 4/1/2016
          rmns_vmec=>rmns, zmnc_vmec=>zmnc, lmnc_vmec=>lmnc,             &  ! MRC 12/1/2016
          phipf_vmec=>phipf, presf_vmec=>presf, nfp_vmec=>nfp,           &
          wb_vmec=>wb, wp_vmec=>wp, gamma_vmec=>gamma,                   &
          volume_vmec=>volume, raxis_vmec=>raxis, lasym_vmec=>lasym,     &
          iasym_vmec=>iasym
      USE descriptor_mod, ONLY: iam
      USE timer_mod
      USE utilities, ONLY: to_full_mesh
      USE siesta_namelist, ONLY: l_vessel

      IMPLICIT NONE

!*******************************************************************************
!  metrics module variables
!*******************************************************************************
!>  Jacobian on half grid
      REAL(dp), DIMENSION(:,:), ALLOCATABLE :: sqrtg
!>  Lower metric tensor half grid. e_s . e_s
      REAL(dp), DIMENSION(:,:), ALLOCATABLE :: gss
!>  Lower metric tensor half grid. e_s . e_u
      REAL(dp), DIMENSION(:,:), ALLOCATABLE :: gsu
!>  Lower metric tensor half grid. e_s . e_v
      REAL(dp), DIMENSION(:,:), ALLOCATABLE :: gsv
!>  Lower metric tensor half grid. e_u . e_u
      REAL(dp), DIMENSION(:,:), ALLOCATABLE :: guu
!>  Lower metric tensor half grid. e_u . e_v
      REAL(dp), DIMENSION(:,:), ALLOCATABLE :: guv
!>  Lower metric tensor half grid. e_v . e_v
      REAL(dp), DIMENSION(:,:), ALLOCATABLE :: gvv
!>  Upper metric tensor half grid. e^s . e^s
      REAL(dp), DIMENSION(:,:), ALLOCATABLE :: hss
!>  Upper metric tensor. e^s . e^u
      REAL(dp), DIMENSION(:,:), ALLOCATABLE :: hsu
!>  Upper metric tensor. e^s . e^v
      REAL(dp), DIMENSION(:,:), ALLOCATABLE :: hsv
!>  Upper metric tensor. e^u . e^u
      REAL(dp), DIMENSION(:,:), ALLOCATABLE :: huu
!>  Upper metric tensor. e^u . e^v
      REAL(dp), DIMENSION(:,:), ALLOCATABLE :: huv
!>  Upper metric tensor. e^v . e^v
      REAL(dp), DIMENSION(:,:), ALLOCATABLE :: hvv
!>  Lower metric tensor full grid. e_s . e_s
      REAL(dp), DIMENSION(:,:), ALLOCATABLE :: gssf
!>  Lower metric tensor full grid. e_s . e_u
      REAL(dp), DIMENSION(:,:), ALLOCATABLE :: gsuf
!>  Lower metric tensor full grid. e_s . e_v
      REAL(dp), DIMENSION(:,:), ALLOCATABLE :: gsvf
!>  Lower metric tensor full grid. e_u . e_u
      REAL(dp), DIMENSION(:,:), ALLOCATABLE :: guuf
!>  Lower metric tensor full grid. e_u . e_v
      REAL(dp), DIMENSION(:,:), ALLOCATABLE :: guvf
!>  Lower metric tensor full grid. e_v . e_v
      REAL(dp), DIMENSION(:,:), ALLOCATABLE :: gvvf

!>  Maximum of the grid R inside the last closed flux surface.
      REAL(dp)                              :: rmax
!>  Minimum of the grid R inside the last closed flux surface.
      REAL(dp)                              :: rmin
!>  Maximum of the grid Z inside the last closed flux surface.
      REAL(dp)                              :: zmax
!>  Minimum of the grid Z inside the last closed flux surface.
      REAL(dp)                              :: zmin

!>  Start timer
      REAL(dp)                              :: skston
!>  Stop timer
      REAL(dp)                              :: skstoff

      CONTAINS

!*******************************************************************************
!  UTILITY SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  @brief Initialize the metric elements.
!>
!>  Loads values on to the SIESTA meshes. This involves splining and
!>  interpolating the VMEC quantities from the VMEC radial grid to the SIESTA
!>  radial grid.
!-------------------------------------------------------------------------------
      SUBROUTINE init_metric_elements()
      USE timer_mod, ONLY: init_timers
      USE shared_data, ONLY: torflux, polflux
      USE island_params, hs=>hs_i, ns=>ns_i

      IMPLICIT NONE

!  Start of executable code
!  Compute half-mesh lower and upper metric elements and jacobian.
      CALL half_mesh_metrics(r1_i, ru_i, rv_i, z1_i, zu_i, zv_i)
!  Convert half grid elements to the full grid.
      CALL full_mesh_metrics

      END SUBROUTINE init_metric_elements

!-------------------------------------------------------------------------------
!>  @brief Set grid sizes.
!>
!>  The real space grid is determined from the number of toroidal and poloidal
!>  modes.
!>
!>  @param[in] mpol_in   Number of SIESTA poloidal modes.
!>  @param[in] ntor_in   Number of SIESTA toroidal modes.
!>  @param[in] nfp_in    Number of field periods.
!>  @param[in] tor_modes Toroidal mode numbers.
!-------------------------------------------------------------------------------
      SUBROUTINE set_grid_sizes(mpol_in, ntor_in, nfp_in, tor_modes)
      USE island_params
      USE shared_data, ONLY: mblk_size, ndims, lasym

      IMPLICIT NONE

!  Declare Arguments
      INTEGER, INTENT(IN)                              :: mpol_in
      INTEGER, INTENT(IN)                              :: ntor_in
      INTEGER, INTENT(IN)                              :: nfp_in
      INTEGER, DIMENSION(-ntor_in:ntor_in), INTENT(in) :: tor_modes

!  Start of executable code
      ndims = 3
      IF (lasym) THEN
         ndims = 2*ndims
      END IF

      mpol_i = mpol_in
      ntor_i = ntor_in

!  Set number of points == number of modes for now! (May want mid-points for
!  flux conservation)
      nu_i = mpol_i + 2
      nv_i = 2*tor_modes(ntor_i) + 2

!USE 3/2 (ORSZAG) RULE FOR ANTI-ALIASING OF EVOLUTION EQNS
!Suppresses RADIAL grid separation in pressure
      nu_i = 3*nu_i/2
!SPH051617      nu_i = nu_i+MOD(nu_i,2)
      IF (lasym) THEN
         nu_i = 2*nu_i
      ELSE IF (MOD(nu_i, 2) .ne. 1) THEN
         nu_i = nu_i + 1
      END IF

      nv_i = 3*nv_i/2
      nv_i = nv_i + MOD(nv_i, 2)

      IF (ntor_i .EQ. 0) THEN
         nv_i = 1
      END IF
      nuv_i = nu_i*nv_i
      mnmax_i = (mpol_i + 1)*(2*ntor_i + 1) ! Added RS. Contains total number of modes.

      mblk_size = ndims*mnmax_i
      CALL ASSERT(mblk_size .NE. 0, 'mblk_size = 0 in set_grid_sizes')

      END SUBROUTINE

!-------------------------------------------------------------------------------
!>  @brief Surface average a quantity.
!>
!>  @param[out] average Surface average
!>  @param[in]  q3d     3D quantity is real space.
!>  @param[in]  nsmin   Minimum radial index.
!>  @param[in]  nsmax   Maximum radial index.
!-------------------------------------------------------------------------------
      SUBROUTINE SurfAvg(average, q3d, nsmin, nsmax)
      USE island_params, ONLY: fourier_context

      IMPLICIT NONE

!  Declare Arguments
      REAL(dp), DIMENSION(nsmin:nsmax), INTENT(out)          :: average
      REAL(dp), DIMENSION(nu_i,nv_i,nsmin:nsmax), INTENT(IN) :: q3d
      INTEGER, INTENT(IN)                                    :: nsmin
      INTEGER, INTENT(IN)                                    :: nsmax

!  Local variables
      INTEGER                                                :: l
      INTEGER                                                :: js
      REAL(dp), ALLOCATABLE, DIMENSION(:,:,:)                :: wint

!  Start of executable code
      ALLOCATE(wint(nu_i,nv_i,nsmin:nsmax))
      DO l = 1, nu_i
         wint(l,:,:) = fourier_context%cosmui(l,0)
      END DO
      DO js = nsmin, nsmax
         average(js) = SUM(q3d(:,:,js)*wint(:,:,js))
      END DO
      DEALLOCATE(wint)

      END SUBROUTINE

!-------------------------------------------------------------------------------
!>  @brief Load the R, Z and lambda arrays from VMEC.
!>
!>  R, Z and Lambda are recast onto the SIESTA grid.
!>
!>  @param[out] istat Error status.
!-------------------------------------------------------------------------------
      SUBROUTINE LoadGrid(istat)
      USE vmec_info

      IMPLICIT NONE

!  Declare Arguments
      INTEGER, INTENT(out) :: istat

!  Start executable code.
      ALLOCATE (r1_i(nu_i, nv_i, ns_i), z1_i(nu_i, nv_i, ns_i),         &
                ru_i(nu_i, nv_i, ns_i), zu_i(nu_i, nv_i, ns_i),         &
                rv_i(nu_i, nv_i, ns_i), zv_i(nu_i, nv_i, ns_i),         &
                stat = istat)

      CALL rz_to_ijsp(rmnc_i, zmns_i, .false.)
      IF (lasym) THEN
         CALL rz_to_ijsp(rmns_i, zmnc_i, .true.)
      END IF

      rmax = MAXVAL(r1_i)
      rmin = MINVAL(r1_i)
      zmax = MAXVAL(z1_i)
      zmin = MINVAL(z1_i)
      IF (.NOT.lasym) THEN
         zmin = MIN(zmin, -zmax) !  ONLY TOP HALF
      END IF

      END SUBROUTINE

!-------------------------------------------------------------------------------
!>  @brief Transform from fourier R Z to real space qauntities.
!>
!>  Computes R, Z, dr/du, dr/dv, dz/du and dz/dv from the fourier
!>  representation.
!>
!>  @param[in] rmn   Radial fourier amplitudes.
!>  @param[in] zmn   Vertical fourier amplitudes.
!>  @param[in] lasym A symmetric flag.
!-------------------------------------------------------------------------------
      SUBROUTINE rz_to_ijsp(rmn, zmn, lasym)
      USE island_params, ONLY: fourier_context
      USE fourier, ONLY: f_cos, f_sin, f_none, f_sum, f_du, f_dv

      IMPLICIT NONE

!  Declare Arguments
      REAL(dp), DIMENSION(:,:,:), INTENT(in) :: rmn
      REAL(dp), DIMENSION(:,:,:), INTENT(in) :: zmn
      LOGICAL, INTENT(in)                    :: lasym

!  Local variables
      INTEGER                                :: sum
      INTEGER                                :: m
      INTEGER                                :: n
      INTEGER                                :: rmode
      INTEGER                                :: zmode

!  Start executable code.
      IF (lasym) THEN
         sum = f_sum
         rmode = f_sin
         zmode = f_cos
      ELSE
         sum = f_none
         rmode = f_cos
         zmode = f_sin
      END IF
      m = IOR(sum, f_du)
      n = IOR(sum, f_dv)

      CALL fourier_context%toijsp(rmn, r1_i, sum, rmode)
      CALL fourier_context%toijsp(rmn, ru_i, m,   rmode)
      CALL fourier_context%toijsp(rmn, rv_i, n,   rmode)

      CALL fourier_context%toijsp(zmn, z1_i, sum, zmode)
      CALL fourier_context%toijsp(zmn, zu_i, m,   zmode)
      CALL fourier_context%toijsp(zmn, zv_i, n,   zmode)

      END SUBROUTINE

!-------------------------------------------------------------------------------
!>  @brief Compute the metric elements on the half mesh.
!>
!>  @param[in] r1_i R on the full mesh.
!>  @param[in] ru_i dRdu on the full mesh.
!>  @param[in] rv_i dRdv on the full mesh.
!>  @param[in] z1_i Z on the full mesh.
!>  @param[in] zu_i dZdu on the full mesh.
!>  @param[in] zv_i dZdv on the full mesh.
!-------------------------------------------------------------------------------
      SUBROUTINE half_mesh_metrics(r1_i, ru_i, rv_i, z1_i, zu_i, zv_i)
      USE island_params, nuv=>nuv_i, nu=>nu_i, nv=>nv_i!, ns=>ns_i

      IMPLICIT NONE  

!  Declare Arguments
      REAL(dp), DIMENSION(nuv,ns_i), INTENT(in) :: r1_i
      REAL(dp), DIMENSION(nuv,ns_i), INTENT(in) :: ru_i
      REAL(dp), DIMENSION(nuv,ns_i), INTENT(in) :: rv_i
      REAL(dp), DIMENSION(nuv,ns_i), INTENT(in) :: z1_i
      REAL(dp), DIMENSION(nuv,ns_i), INTENT(in) :: zu_i
      REAL(dp), DIMENSION(nuv,ns_i), INTENT(in) :: zv_i

!  Local variables
      INTEGER                                   :: js
      INTEGER                                   :: istat
      INTEGER                                   :: imin(2)
      INTEGER                                   :: imax(2)
      REAL(dp)                                  :: mintest
      REAL(dp)                                  :: maxtest
      REAL(dp), DIMENSION(nuv)                  :: r12
      REAL(dp), DIMENSION(nuv)                  :: rs12
      REAL(dp), DIMENSION(nuv)                  :: ru12
      REAL(dp), DIMENSION(nuv)                  :: rv12
      REAL(dp), DIMENSION(nuv)                  :: zs12
      REAL(dp), DIMENSION(nuv)                  :: zu12
      REAL(dp), DIMENSION(nuv)                  :: zv12

!  Local parameters.
      REAL(dp), PARAMETER                     :: p5 = 1.d0/2.d0
      REAL(dp), PARAMETER                     :: zero = 0

!  Start executable code.

!  ALLOCATE METRIC ELEMENT ARRAYS
      ALLOCATE(sqrtg(nuv,ns_i),                                                &
               gss(nuv,ns_i), gsu(nuv,ns_i), gsv(nuv,ns_i),                    &
               guu(nuv,ns_i), guv(nuv,ns_i), gvv(nuv,ns_i), stat=istat)

!  COMPUTE ALL ON THE HALF MESH
      sqrtg(:,1) = 0

      DO js = 2, ns_i
         r12  = p5*(r1_i(:,js) + r1_i(:,js - 1))
         rs12 =    (r1_i(:,js) - r1_i(:,js - 1))*ohs_i
         ru12 = p5*(ru_i(:,js) + ru_i(:,js - 1))
         rv12 = p5*(rv_i(:,js) + rv_i(:,js - 1))
         zs12 =    (z1_i(:,js) - z1_i(:,js - 1))*ohs_i
         zu12 = p5*(zu_i(:,js) + zu_i(:,js - 1))
         zv12 = p5*(zv_i(:,js) + zv_i(:,js - 1))
         gss(:,js) = rs12*rs12 + zs12*zs12
         gsu(:,js) = rs12*ru12 + zs12*zu12
         gsv(:,js) = rs12*rv12 + zs12*zv12
         guu(:,js) = ru12*ru12 + zu12*zu12
         guv(:,js) = ru12*rv12 + zu12*zv12
         gvv(:,js) = rv12*rv12 + r12*r12 + zv12*zv12
         sqrtg(:,js) = r12*(ru12*zs12 - rs12*zu12)
      END DO

      sqrtg(:,1) = sqrtg(:,2)

!NEED THESE FOR CURRENT CALCULATION AT ORIGIN
      gss(:,1) = gss(:,2)
      gsu(:,1) = gsu(:,2)
      gsv(:,1) = gsv(:,2)
      guu(:,1) = guu(:,2)
      guv(:,1) = guv(:,2)
      gvv(:,1) = gvv(:,2)

      mintest = MINVAL(sqrtg(:,2:))
      maxtest = MAXVAL(sqrtg(:,2:))
      
      IF (mintest*maxtest .LE. zero) THEN
         imin = MINLOC(sqrtg(:,2:))
         imax = MAXLOC(sqrtg(:,2:))
         IF (iam .EQ. 0) THEN
            WRITE(*,1000) mintest, imin(2), maxtest, imax(2)
         END IF
         CALL ASSERT(mintest*maxtest.GT.zero,                                  &
                     ' Jacobian changed sign in half_mesh_metrics')
      END IF

1000  FORMAT(' MIN-G: ',1p,e12.4,' AT JS: ',i5,                                &
             ' MAX-G: ',1p,e12.4,' AT JS: ',i5)

      END SUBROUTINE

!-------------------------------------------------------------------------------
!>  @brief Gets the full grid metric elements.
!>
!>  This subroutine gets the lower metric elements on the full mesh. Preserves
!>  positive definiteness of metric tensor.
!-------------------------------------------------------------------------------
      SUBROUTINE full_mesh_metrics
      USE island_params, nuv=>nuv_i

      IMPLICIT NONE  

!  Local variables
      INTEGER                                :: istat
      INTEGER                                :: js
      REAL (dp), DIMENSION(:), ALLOCATABLE   :: det
      REAL (dp), DIMENSION(:,:), ALLOCATABLE :: temph

!  Start executable code.
      ALLOCATE(gssf(nuv,ns_i), guuf(nuv,ns_i), gvvf(nuv,ns_i),                 &
               gsuf(nuv,ns_i), gsvf(nuv,ns_i), guvf(nuv,ns_i),                 &
               hss(nuv,ns_i),  huu(nuv,ns_i),  hvv(nuv,ns_i),                  &
               hsu(nuv,ns_i),  hsv(nuv,ns_i),  huv(nuv,ns_i), stat=istat)

      CALL to_full_mesh(gss, gssf)
      CALL to_full_mesh(guu, guuf)
      CALL to_full_mesh(gvv, gvvf)
      CALL to_full_mesh(gsu, gsuf)
      CALL to_full_mesh(gsv, gsvf)
      CALL to_full_mesh(guv, guvf)

      ALLOCATE(det(nuv))

!  Compute upper metric elements (inverse of lower matrix) on full mesh.
!  NOTE: DET == 1/det|Gij| == 1/sqrtg(full)**2
      DO js = 1, ns_i
         det(:) = gsvf(:,js)*gsvf(:,js)*guuf(:,js)                             &
                - 2.0*gsuf(:,js)*gsvf(:,js)*guvf(:,js)                         &
                + gssf(:,js)*guvf(:,js)*guvf(:,js)                             &
                + gsuf(:,js)*gsuf(:,js)*gvvf(:,js)                             &
                - gssf(:,js)*guuf(:,js)*gvvf(:,js)

          hss(:,js) = (guvf(:,js)*guvf(:,js) - guuf(:,js)*gvvf(:,js))/det
          hsu(:,js) = (gsuf(:,js)*gvvf(:,js) - gsvf(:,js)*guvf(:,js))/det
          hsv(:,js) = (gsvf(:,js)*guuf(:,js) - gsuf(:,js)*guvf(:,js))/det
          huu(:,js) = (gsvf(:,js)*gsvf(:,js) - gssf(:,js)*gvvf(:,js))/det
          huv(:,js) = (gssf(:,js)*guvf(:,js) - gsuf(:,js)*gsvf(:,js))/det
          hvv(:,js) = (gsuf(:,js)*gsuf(:,js) - gssf(:,js)*guuf(:,js))/det
      END DO

      guuf(:,1) = 0
      gsuf(:,1) = 0
      guvf(:,1) = 0

      hss(:,1) = 0
      hsu(:,1) = 0
      hsv(:,1) = 0
      huv(:,1) = 0
      hvv(:,1) = 0

      DEALLOCATE(det)

      CALL check_metrics

1000  FORMAT('Determinant |gijf| <= 0 at js = ',i4)

      END SUBROUTINE full_mesh_metrics

!-------------------------------------------------------------------------------
!>  @ brief Test to check that we computed the upper metric elements correctly.
!>
!>  Check that the metric elements conform to the following relation to machine
!>  precission.
!>
!>    hij*gjk = delta_ik                                                     (1)
!>
!>  Note the repeated index implys a sum.
!-------------------------------------------------------------------------------
      SUBROUTINE check_metrics

      IMPLICIT NONE

!  Local Parameters
      REAL (dp), PARAMETER :: tolarance = 1.0E15

!  Start executable code.
      CALL ASSERT(ALL(ABS(hss*gssf + hsu*gsuf + hsv*gsvf - 1) .lt. tolarance), &
                      's Diagonal metric element failed.')
      CALL ASSERT(ALL(ABS(hsu*gsuf + huu*guuf + huv*guvf - 1) .lt. tolarance), &
                      'u Diagonal metric element failed.')
      CALL ASSERT(ALL(ABS(hsv*gsvf + huv*guvf + hvv*gvvf - 1) .lt. tolarance), &
                      'v Diagonal metric element failed.')
      CALL ASSERT(ALL(ABS(hss*gsuf + hsu*guuf + hsv*guvf) .lt. tolarance),     &
                      'su Off diagonal metric element failed.')
      CALL ASSERT(ALL(ABS(hss*gsvf + hsu*guvf + hsv*gvvf) .lt. tolarance),     &
                      'sv Off diagonal metric element failed.')
      CALL ASSERT(ALL(ABS(hsu*gsvf + huu*guvf + huv*gvvf) .lt. tolarance),     &
                      'uv Off diagonal metric element failed.')

      END SUBROUTINE

!-------------------------------------------------------------------------------
!>  @brief Converts to contravariant component full grid.
!>
!>  Computes contravariant (upper) from covariant (lower)  components on the
!>  full mesh.
!>
!>    xsupsij = hss*xsubsij + hsu*xsupuij + hsv*xsubvij
!>    xsupuij = hsu*xsubsij + huu*xsupuij + huv*xsubvij
!>    xsupvij = hsv*xsubsij + huv*xsupuij + hvv*xsubvij
!>
!>  @param[in]  xsubsij Covariant s component.
!>  @param[in]  xsubuij Covariant u component.
!>  @param[in]  xsubvij Covariant v component.
!>  @param[out] xsupsij Contravariant s component.
!>  @param[out] xsupuij Contravariant u component.
!>  @param[out] xsupvij Contravariant v component.
!>  @param[in]  nsmin   Min radial index.
!>  @param[in]  nsmax   Max radial index.
!-------------------------------------------------------------------------------
      SUBROUTINE toupper(xsubsij, xsubuij, xsubvij,                            &
                         xsupsij, xsupuij, xsupvij, nsmin, nsmax)
      USE stel_kinds
      USE island_params, nuv=>nuv_i

      IMPLICIT NONE

!  Declare Arguments
      REAL(dp), DIMENSION(nuv,nsmin:nsmax), INTENT(in)  :: xsubsij
      REAL(dp), DIMENSION(nuv,nsmin:nsmax), INTENT(in)  :: xsubuij
      REAL(dp), DIMENSION(nuv,nsmin:nsmax), INTENT(in)  :: xsubvij
      REAL(dp), DIMENSION(nuv,nsmin:nsmax), INTENT(out) :: xsupsij
      REAL(dp), DIMENSION(nuv,nsmin:nsmax), INTENT(out) :: xsupuij
      REAL(dp), DIMENSION(nuv,nsmin:nsmax), INTENT(out) :: xsupvij
      INTEGER, INTENT(in)                               :: nsmin
      INTEGER, INTENT(in)                               :: nsmax

!  Local variables
        INTEGER                                         :: js

!  Start executable code.
      DO js = nsmin, nsmax
         xsupsij(:,js) = hss(:,js)*xsubsij(:,js)                               &
                       + hsu(:,js)*xsubuij(:,js)                               &
                       + hsv(:,js)*xsubvij(:,js)
         xsupuij(:,js) = hsu(:,js)*xsubsij(:,js)                               &
                       + huu(:,js)*xsubuij(:,js)                               &
                       + huv(:,js)*xsubvij(:,js)
         xsupvij(:,js) = hsv(:,js)*xsubsij(:,js)                               &
                       + huv(:,js)*xsubuij(:,js)                               &
                       + hvv(:,js)*xsubvij(:,js)
      END DO

      END SUBROUTINE

!-------------------------------------------------------------------------------
!>  @brief Converts to covariant component half grid.
!>
!>  Computes covariant (lower) from contravariant (upper) magnetic field
!>  components on the half mesh.
!>
!>    xsubsij = gss*xsupsij + gsu*xsupuij + gsv*xsupvij
!>    xsubuij = gsu*xsupsij + guu*xsupuij + guv*xsupvij
!>    xsubvij = gsv*xsupsij + guv*xsupuij + gvv*xsupvij
!>
!>  @param[in]  xsupsij Contravariant s component.
!>  @param[in]  xsupuij Contravariant u component.
!>  @param[in]  xsupvij Contravariant v component.
!>  @param[out] xsubsij Covariant s component.
!>  @param[out] xsubuij Covariant u component.
!>  @param[out] xsubvij Covariant v component.
!>  @param[in]  nsmin   Min radial index.
!>  @param[in]  nsmax   Max radial index.
!-------------------------------------------------------------------------------
      SUBROUTINE tolowerh(xsupsij, xsupuij, xsupvij,                      &
                          xsubsij, xsubuij, xsubvij, nsmin, nsmax)
      USE stel_kinds
      USE island_params, ONLY: nuv=>nuv_i

      IMPLICIT NONE

!  Declare Arguments
      REAL(dp), DIMENSION(nuv,nsmin:nsmax), INTENT(in)  :: xsupsij
      REAL(dp), DIMENSION(nuv,nsmin:nsmax), INTENT(in)  :: xsupuij
      REAL(dp), DIMENSION(nuv,nsmin:nsmax), INTENT(in)  :: xsupvij
      REAL(dp), DIMENSION(nuv,nsmin:nsmax), INTENT(out) :: xsubsij
      REAL(dp), DIMENSION(nuv,nsmin:nsmax), INTENT(out) :: xsubuij
      REAL(dp), DIMENSION(nuv,nsmin:nsmax), INTENT(out) :: xsubvij
      INTEGER, INTENT(IN)                               :: nsmin
      INTEGER, INTENT(IN)                               :: nsmax

!  Start executable code.
      xsubsij = gss(:,nsmin:nsmax)*xsupsij                               &
              + gsu(:,nsmin:nsmax)*xsupuij                               &
              + gsv(:,nsmin:nsmax)*xsupvij

      xsubuij = gsu(:,nsmin:nsmax)*xsupsij                               &
              + guu(:,nsmin:nsmax)*xsupuij                               &
              + guv(:,nsmin:nsmax)*xsupvij
      
      xsubvij = gsv(:,nsmin:nsmax)*xsupsij                               &
              + guv(:,nsmin:nsmax)*xsupuij                               &
              + gvv(:,nsmin:nsmax)*xsupvij

      END SUBROUTINE tolowerh

!-------------------------------------------------------------------------------
!>  @brief Converts to covariant component full grid.
!>
!>  Computes covariant (lower) from contravariant (upper) magnetic field
!>  components on the full mesh.
!>
!>    xsubsij = gssf*xsupsij + gsuf*xsupuij + gsvf*xsupvij
!>    xsubuij = gsuf*xsupsij + guuf*xsupuij + guvf*xsupvij
!>    xsubvij = gsvf*xsupsij + guvf*xsupuij + gvvf*xsupvij
!>
!>  @param[in]  xsupsij Contravariant s component.
!>  @param[in]  xsupuij Contravariant u component.
!>  @param[in]  xsupvij Contravariant v component.
!>  @param[out] xsubsij Covariant s component.
!>  @param[out] xsubuij Covariant u component.
!>  @param[out] xsubvij Covariant v component.
!>  @param[in]  nsmin   Min radial index.
!>  @param[in]  nsmax   Max radial index.
!-------------------------------------------------------------------------------
      SUBROUTINE tolowerf(xsupsij, xsupuij, xsupvij,                     &
                          xsubsij, xsubuij, xsubvij, nsmin, nsmax)
      USE stel_kinds
      USE island_params, ONLY: ns=>ns_i, nuv=>nuv_i

      IMPLICIT NONE

!  Declare Arguments
      REAL(dp), DIMENSION(nuv,ns), INTENT(in)  :: xsupsij  !  MRC: Why ns and not nsmin:nsmax here?
      REAL(dp), DIMENSION(nuv,ns), INTENT(in)  :: xsupuij
      REAL(dp), DIMENSION(nuv,ns), INTENT(in)  :: xsupvij
      REAL(dp), DIMENSION(nuv,ns), INTENT(out) :: xsubsij
      REAL(dp), DIMENSION(nuv,ns), INTENT(out) :: xsubuij
      REAL(dp), DIMENSION(nuv,ns), INTENT(out) :: xsubvij
      INTEGER, INTENT(in)                        :: nsmin
      INTEGER, INTENT(in)                        :: nsmax

!  Local variables
      INTEGER                                    :: js
      INTEGER                                    :: js1

!  Start executable code.
      DO js = nsmin, nsmax
         js1 = js - nsmin + 1
         xsubsij(:,js1) = gssf(:,js)*xsupsij(:,js1)                      &
                        + gsuf(:,js)*xsupuij(:,js1)                      &
                        + gsvf(:,js)*xsupvij(:,js1)

         xsubuij(:,js1) = gsuf(:,js)*xsupsij(:,js1)                      &
                        + guuf(:,js)*xsupuij(:,js1)                      &
                        + guvf(:,js)*xsupvij(:,js1)

         xsubvij(:,js1) = gsvf(:,js)*xsupsij(:,js1)                      &
                        + guvf(:,js)*xsupuij(:,js1)                      &
                        + gvvf(:,js)*xsupvij(:,js1)
      END DO

!  FIXME: CHECK THIS
      IF (nsmin .eq. 1) THEN
         xsubuij(:,1) = 0
      END IF

      END SUBROUTINE tolowerf

!-------------------------------------------------------------------------------
!>  @brief Deallocate memory containing metric elements on the half mesh.
!>
!>  Also deallocates the grid variables.
!>
!>  @note
!>  sqrtg is deallocated in @ref init_bcovar and stored in @ref jacob variable
!-------------------------------------------------------------------------------
      SUBROUTINE cleanup_metric_elements
      USE island_params, ONLY: fourier_context
      USE vmec_info

      IMPLICIT NONE

!  Local variables
      INTEGER :: istat

!  Start executable code.
      DEALLOCATE(gss, gsu, gsv, guu, guv, gvv,                           &
                 hss, hsu, hsv, huu, huv, hvv,                           &
                 phipf_i, chipf_i, presf_i, stat=istat)
      
      DEALLOCATE (r1_i, z1_i, ru_i, zu_i, rv_i, zv_i)

      DEALLOCATE(fourier_context)
      CALL vmec_info_destruct_island

      END SUBROUTINE cleanup_metric_elements

!-------------------------------------------------------------------------------
!>  @brief Deallocate memory containing metric elements on the full mesh.
!-------------------------------------------------------------------------------
      SUBROUTINE dealloc_full_lower_metrics
      USE stel_kinds

      IMPLICIT NONE

!  Local variables
      INTEGER :: istat

!  Start executable code.
      DEALLOCATE(gssf, guuf, gvvf, gsuf, gsvf, guvf, stat = istat)
!      IF (istat .NE. 0) CALL ErrorAbort('Problem dealloc. in LOWER_METRIC_FULL')

      END SUBROUTINE dealloc_full_lower_metrics

      END MODULE
