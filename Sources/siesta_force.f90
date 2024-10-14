!*******************************************************************************
!>  @file siesta_force.f90
!>  @brief Contains the @ref siesta_force module.
!
!  Note separating the Doxygen comment block here so detailed decription is
!  found in the Module not the file.
!
!>  Compute the JxB - Grad(p) covariant force components. The plasma is in
!>  equilibrium when the force of the magnetic pressure balances the plasma
!> pressure.
!>
!>     F = J X B - Grad(p) = 0                                               (1)
!>
!>  where J is the current density. Updates for B (dB) and p (dp) were computed
!>  in @ref siesta_bfield and @ref siesta_pressure, respectively.
!>
!>     B_i+1 = B_i + dB                                                     (2a)
!>
!>     p_i+1 = p_i + dp                                                     (2b)
!>
!>  Note that the magnetic field and pressure terms contain a Jacobian factor. 
!>  Updated current density can be found from the updated magnetic field.
!>
!>     J = Curl(B)/mu0                                                       (3)
!>
!>  Actual computation of the contravariant current components is implimented in
!>  @ref cv_currents. The curvilinear coordinates, the J X B term becomes
!>
!>     (J X B)_s = (K^u*B^v - K^v*B^u)                                     (4a)
!>
!>     (J X B)_u = (K^v*B^s - K^s*B^v)                                     (4b)
!>
!>     (J X B)_v = (K^s*B^u - K^u*B^s)                                     (4c)
!>
!>  where K^i = sqrt(g) J^i.
!*******************************************************************************
      MODULE siesta_force
      USE stel_kinds
      USE stel_constants
      USE v3_utilities, ONLY: assert
      USE descriptor_mod, ONLY: SIESTA_COMM, iam, nprocs
      USE fourier 
      USE island_params, ONLY: nu_i, hs_i
      USE timer_mod
      USE shared_data, ONLY: l_linearize, l_getwmhd, lverbose, l_pedge
      USE nscalingtools, ONLY: startglobrow, endglobrow, MPI_ERR
      USE utilities, ONLY: GradientFull, to_full_mesh
      USE quantities
      USE mpi_inc

      IMPLICIT NONE

      CONTAINS

!*******************************************************************************
!  UTILITY SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  @brief Update MHD forces on full radial mesh.
!>
!>  Updates force using the advanced values of B and p obtained from calls to
!>  @ref siesta_pressure::update_pres and @ref siesta_bfield::update_bfield
!>  Update values of fsub*mn*f stored in @ref quantities. Linearized force is
!>
!>    Flin ~ delta_v
!>
!>  has only the linear part. Not DC part and not non-linear.
!-------------------------------------------------------------------------------
      SUBROUTINE update_force
      USE shared_data, ONLY: gc, col_scale
      USE siesta_currents, ONLY: cv_currents
      USE hessian, ONLY: Apply_ColScale
      USE bhtobf

      IMPLICIT NONE

!  local variables
      REAL (dp), DIMENSION(:,:,:), ALLOCATABLE :: bsupsijh
      REAL (dp), DIMENSION(:,:,:), ALLOCATABLE :: bsupuijh
      REAL (dp), DIMENSION(:,:,:), ALLOCATABLE :: bsupvijh
      REAL (dp), DIMENSION(:,:,:), ALLOCATABLE :: pijh
      REAL (dp), DIMENSION(:,:,:), ALLOCATABLE :: pijf1
      REAL (dp), DIMENSION(:,:,:), ALLOCATABLE :: KxBsij
      REAL (dp), DIMENSION(:,:,:), ALLOCATABLE :: KxBuij
      REAL (dp), DIMENSION(:,:,:), ALLOCATABLE :: KxBvij
      REAL (dp), DIMENSION(:,:), ALLOCATABLE   :: pardamp
      INTEGER                                  :: istat            ,m,n
      INTEGER                                  :: n1
      INTEGER                                  :: n2
      REAL (dp)                                :: ton
      REAL (dp)                                :: toff
      REAL (dp)                                :: skston
      REAL (dp)                                :: skstoff
      INTEGER                                  :: nsmin
      INTEGER                                  :: nsmax

!  Start of executable code
      CALL second0(skston)
      nsmin = MAX(1, startglobrow)
      nsmax = MIN(ns, endglobrow + 1)

      ALLOCATE(bsupsijh(ntheta,nzeta,nsmin:nsmax),                             &
               bsupuijh(ntheta,nzeta,nsmin:nsmax),                             &
               bsupvijh(ntheta,nzeta,nsmin:nsmax),                             &
               pijh(ntheta,nzeta,nsmin:nsmax), stat=istat)
      CALL assert_eq(0, istat, 'Allocation1 failed in UPDATE_FORCE')

!  Compute real-space contravariant components of B and p. Fields are updated
!  in update state, not here. Unperturbed fields are unchanged until
!  update_state call. Perturbed fields were calculated in update_bfield and
!  update_press.
      CALL IncFields(jbsupsmnsh,  jbsupumnch,  jbsupvmnch,  jpmnch,            &
                     djbsupsmnsh, djbsupumnch, djbsupvmnch, djpmnch,           &
                     bsupsijh, bsupuijh, bsupvijh, pijh, f_sin)
      IF (lasym) THEN
         CALL IncFields(jbsupsmnch,  jbsupumnsh,  jbsupvmnsh,  jpmnsh,         &
                        djbsupsmnch, djbsupumnsh, djbsupvmnsh, djpmnsh,        &
                        bsupsijh, bsupuijh, bsupvijh, pijh, f_cos)
      END IF

!  Update thermal energy (pressure based). pijh contains the jacobian term at
!  this point.
      WPRES: IF (l_getwmhd) THEN
         CALL assert(.not.INHESSIAN,                                          &
                     'l_getwmhd must be set to FALSE in Hessian')
         n2 = MIN(ns, endglobrow)
         n1 = MAX(2, nsmin)
         wp = signjac*(twopi*twopi*hs_i)*SUM(pijh(:,:,n1:n2)*wint(:,:,n1:n2))
#if defined(MPI_OPT)
         IF (PARSOLVER) THEN
!  FIXME: All reduce is not deterministic. This causes a divergent run sequence.
            CALL MPI_ALLREDUCE(MPI_IN_PLACE, wp, 1, MPI_REAL8, MPI_SUM,        &
                               SIESTA_COMM, MPI_ERR)
         END IF
#endif
      END IF WPRES

      ALLOCATE(pijf1(ntheta,nzeta,nsmin:nsmax), stat=istat)
      CALL assert_eq(0, istat, 'Allocation2 failed in UPDATE_FORCE')

!  Update full mesh components of bsup* and pressure in real-space by averaging.
!  This removes the jacobian term.
      CALL bhalftobfull(bsupsijh, bsupuijh, bsupvijh,                          &
                        bsupsijf, bsupuijf, bsupvijf, pijh, pijf1)

!  Calculate contravariant components of J = curl B (Ksup = gsqrt*Jsup)
!  and update the magnetic energy if nonlinear force is being computed.
      CALL cv_currents(bsupsijh, bsupuijh, bsupvijh,                           &
                       ksupsijf, ksupuijf, ksupvijf,                           &
                       .not.l_linearize .OR. l_getwmhd, .false.)

      DEALLOCATE(bsupsijh, bsupuijh, bsupvijh, stat=istat)

      nsmax = MIN(ns, endglobrow)

      ALLOCATE(KxBsij(ntheta,nzeta,nsmin:nsmax),                               &
               KxBuij(ntheta,nzeta,nsmin:nsmax),                               &
               KxBvij(ntheta,nzeta,nsmin:nsmax), stat=istat)
      CALL assert_eq(0, istat, 'Allocation3 failed in UPDATE_FORCE')

!  KxBsij used as a scratch array.
      CALL InitParDamping(KxBsij, pardamp)

!  Compute curvilinear Lorentz force components
      IF (l_linearize) THEN
!  dJ X B0
         CALL Lorentz(bsupsijf0, bsupuijf0, bsupvijf0,                         &
                      ksupsijf, ksupuijf, ksupvijf,                            &
                      KxBsij, KxBuij, KxBvij, f_none)
!  J0 X dB
         CALL Lorentz(bsupsijf, bsupuijf, bsupvijf,                            &
                      ksupsijf0, ksupuijf0, ksupvijf0,                         &
                      KxBsij, KxBuij, KxBvij, f_sum)
      ELSE
!  (J0+dJ) X (B0+dB)
         CALL Lorentz(bsupsijf, bsupuijf, bsupvijf,                            &
                      ksupsijf, ksupuijf, ksupvijf,                            &
                      KxBsij, KxBuij, KxBvij, f_none)
      END IF

      IF (nsmax .eq. ns .and. .not.l_vessel) THEN
         CALL assert(ALL(bsupsijf(:,:,ns) .eq. zero),                          &
                     'bsupsijf(ns) != 0 in UPDATE_FORCE')
      END IF

      CALL GetMHDForce(fsubsmncf, fsubumnsf, fsubvmnsf,                        &
                       pijh, KxBsij, KxBuij, KxBvij, pardamp, f_cos)
      IF (lasym) THEN
         CALL GetMHDForce(fsubsmnsf, fsubumncf, fsubvmncf,                     &
                          pijh, KxBsij, KxBuij, KxBvij, pardamp, f_sin)
      END IF

      CALL Apply_ColScale(gc, col_scale, nsmin, nsmax)

      DEALLOCATE(pijf1, KxBsij, KxBuij, KxBvij, pijh, stat=istat)
      IF (ALLOCATED(pardamp)) THEN
         DEALLOCATE(pardamp)
      END IF

      CALL second0(skstoff)
      time_update_force = time_update_force + (skstoff - skston)

      END SUBROUTINE

!-------------------------------------------------------------------------------
!>  @brief Compute nonlinear or linearized contravariant magnetic field
!>         components and pressure in real space.
!>
!>  @param[in]    jbsupsmnh   Contravariant B field in the s direction.
!>  @param[in]    jbsupumnh   Contravariant B field in the u direction.
!>  @param[in]    jbsupvmnh   Contravariant B field in the v direction.
!>  @param[in]    jpmnh       Pressure on the half grid.
!>  @param[in]    djbsupsmnh  Contravariant B field perturbation in the s
!>                            direction.
!>  @param[in]    djbsupumnh  Contravariant B field perturbation in the u
!>                            direction.
!>  @param[in]    djbsupvmnh  Contravariant B field perturbation in the v
!>                            direction.
!>  @param[in]    djpmnh      Pressure perturbation.
!>  @param[out]   jbsupsijh    Real space contravariant B field in the s
!>                            direction.
!>  @param[out]   jbsupuijh    Real space contravariant B field in the u
!>                            !direction.
!>  @param[out]   jbsupvijh    Real space contravariant B field in the v
!>                            direction.
!>  @param[in]    iparity     Fourier parity.
!-------------------------------------------------------------------------------
      SUBROUTINE IncFields(jbsupsmnh, jbsupumnh, jbsupvmnh, jpmnh,             &
                           djbsupsmnh, djbsupumnh, djbsupvmnh, djpmnh,         &
                           jbsupsijh, jbsupuijh, jbsupvijh, jpijh,             &
                           parity)
      USE island_params, ONLY: fourier_context

      IMPLICIT NONE

!  Declare Arguments
      REAL (dp), DIMENSION(:,:,:), ALLOCATABLE, INTENT(in) :: jbsupsmnh
      REAL (dp), DIMENSION(:,:,:), ALLOCATABLE, INTENT(in) :: jbsupumnh
      REAL (dp), DIMENSION(:,:,:), ALLOCATABLE, INTENT(in) :: jbsupvmnh
      REAL (dp), DIMENSION(:,:,:), ALLOCATABLE, INTENT(in) :: jpmnh
      REAL (dp), DIMENSION(:,:,:), ALLOCATABLE, INTENT(in) :: djbsupsmnh
      REAL (dp), DIMENSION(:,:,:), ALLOCATABLE, INTENT(in) :: djbsupumnh
      REAL (dp), DIMENSION(:,:,:), ALLOCATABLE, INTENT(in) :: djbsupvmnh
      REAL (dp), DIMENSION(:,:,:), ALLOCATABLE, INTENT(in) :: djpmnh
      REAL (dp), DIMENSION(:,:,:), INTENT(out)             :: jbsupsijh
      REAL (dp), DIMENSION(:,:,:), INTENT(out)             :: jbsupuijh
      REAL (dp), DIMENSION(:,:,:), INTENT(out)             :: jbsupvijh
      REAL (dp), DIMENSION(:,:,:), INTENT(out)             :: jpijh
      INTEGER, INTENT(in)                                  :: parity

!  local variables
      INTEGER                                              :: fours
      INTEGER                                              :: fouruv
      INTEGER                                              :: fcomb
      INTEGER                                              :: istat
      INTEGER                                              :: nsmin
      INTEGER                                              :: nsmax
      REAL (dp), DIMENSION(:,:,:), ALLOCATABLE             :: jbsupsmnh_i
      REAL (dp), DIMENSION(:,:,:), ALLOCATABLE             :: jbsupumnh_i
      REAL (dp), DIMENSION(:,:,:), ALLOCATABLE             :: jbsupvmnh_i
      REAL (dp), DIMENSION(:,:,:), ALLOCATABLE             :: jpmnh_i

!  Start of executable code
      nsmin = MAX(1, startglobrow)
      nsmax = MIN(ns, endglobrow + 1)

      IF (parity .EQ. f_sin) THEN
         fcomb = f_none
         fours = f_sin
         fouruv = f_cos
      ELSE
         fcomb = f_sum
         fours = f_cos
         fouruv = f_sin
      END IF

!  Allocate space for total values of bsup*mnh_i and pmnh_i.
      ALLOCATE(jbsupsmnh_i(0:mpol,-ntor:ntor,nsmin:nsmax),                     &
               jbsupumnh_i(0:mpol,-ntor:ntor,nsmin:nsmax),                     &
               jbsupvmnh_i(0:mpol,-ntor:ntor,nsmin:nsmax),                     &
               jpmnh_i(0:mpol,-ntor:ntor,nsmin:nsmax), stat=istat)
      CALL assert_eq(0, istat, 'Allocation failed in IncFields')

!  For nonlinear force, add perturbation to unperturbed value. For linear force,
!  compute JUST perturbed fields here.
      jbsupsmnh_i = djbsupsmnh
      jbsupumnh_i = djbsupumnh
      jbsupvmnh_i = djbsupvmnh
      jpmnh_i     = djpmnh
      IF (.not.l_linearize) THEN
         jbsupsmnh_i = jbsupsmnh_i + jbsupsmnh(:,:,nsmin:nsmax)
         jbsupumnh_i = jbsupumnh_i + jbsupumnh(:,:,nsmin:nsmax)
         jbsupvmnh_i = jbsupvmnh_i + jbsupvmnh(:,:,nsmin:nsmax)
         jpmnh_i     = jpmnh_i     + jpmnh(:,:,nsmin:nsmax)
      END IF

!  js=1 origin values. Note: This assumed that the magnetic axis and the grid
!  axis align. Maybe the cause of poor convergence.
      IF (nsmin .EQ. 1) THEN
         jbsupsmnh_i(:,:,1) = 0
!         jbsupsmnh_i(m1,:,1) = jbsupsmnh_i(m1,:,2)
         jbsupumnh_i(:,:,1) = 0
!         jbsupumnh_i(m0:m1,:,1) = jbsupumnh_i(m0:m1,:,2) MRC
         jbsupvmnh_i(:,:,1) = 0
!         jbsupvmnh_i(m0,:,1) = jbsupvmnh_i(m0,:,2)
         jpmnh_i(:,:,1) = 0
!         jpmnh_i(m0,:,1) = jpmnh_i(m0,:,2)
      END IF

!  Calculate real-space jbsup*ihh and jpijh half mesh.
      CALL fourier_context%toijsp(jbsupsmnh_i, jbsupsijh, fcomb, fours)
      CALL fourier_context%toijsp(jbsupumnh_i, jbsupuijh, fcomb, fouruv)
      CALL fourier_context%toijsp(jbsupvmnh_i, jbsupvijh, fcomb, fouruv)
      CALL fourier_context%toijsp(jpmnh_i, jpijh, fcomb, fouruv)

      DEALLOCATE(jbsupsmnh_i, jbsupumnh_i, jbsupvmnh_i, jpmnh_i)

      END SUBROUTINE

!-------------------------------------------------------------------------------
!>  @brief Compute covariant (sub) real-space components of the Lorentz K X B
!>         force.
!>
!>  @param[in]    bsupsijf Contravariant component of the B field in the s
!>                         direction.
!>  @param[in]    bsupuijf Contravariant component of the B field in the u
!>                         direction.
!>  @param[in]    bsupvijf Contravariant component of the B field in the v
!>                         direction.
!>  @param[in]    ksupsijf Contravariant component of the current in the s
!>                         direction.
!>  @param[in]    ksupuijf Contravariant component of the current in the u
!>                         direction.
!>  @param[in]    ksupvijf Contravariant component of the current in the v
!>                         direction.
!>  @param[inout] KxBsij   Covariant component of the lorentz force in the s
!>                         direction.
!>  @param[inout] KxBuij   Covariant component of the lorentz force in the u
!>                         direction.
!>  @param[inout] KxBvij   Covariant component of the lorentz force in the v
!>                         direction.
!>  @param[in]    fcomb    Control flag to sum the components with previous
!>                         values.
!-------------------------------------------------------------------------------
      SUBROUTINE Lorentz(bsupsijf, bsupuijf, bsupvijf,                         &
                         ksupsijf, ksupuijf, ksupvijf,                         &
                         KxBsij, KxBuij, KxBvij, fcomb)

      IMPLICIT NONE

!  Declare Arguments
      REAL (dp), DIMENSION(:,:,:), ALLOCATABLE, INTENT(in)    :: bsupsijf
      REAL (dp), DIMENSION(:,:,:), ALLOCATABLE, INTENT(in)    :: bsupuijf
      REAL (dp), DIMENSION(:,:,:), ALLOCATABLE, INTENT(in)    :: bsupvijf
      REAL (dp), DIMENSION(:,:,:), ALLOCATABLE, INTENT(in)    :: ksupsijf
      REAL (dp), DIMENSION(:,:,:), ALLOCATABLE, INTENT(in)    :: ksupuijf
      REAL (dp), DIMENSION(:,:,:), ALLOCATABLE, INTENT(in)    :: ksupvijf
      REAL (dp), DIMENSION(:,:,:), ALLOCATABLE, INTENT(inout) :: KxBsij
      REAL (dp), DIMENSION(:,:,:), ALLOCATABLE, INTENT(inout) :: KxBuij
      REAL (dp), DIMENSION(:,:,:), ALLOCATABLE, INTENT(inout) :: KxBvij
      INTEGER, INTENT(in)                                     :: fcomb

!  local variables
      INTEGER                                                 :: nsmin
      INTEGER                                                 :: nsmax

!  Start of executable code
      nsmin = MAX(1, startglobrow)
      nsmax = MIN(ns, endglobrow)
      
!  Initialize to avoid NaN errors.
      IF (fcomb .eq. f_none) THEN
         KxBsij = 0
         KxBuij = 0
         KxBvij = 0
      END IF


      KxBsij = KxBsij                                                          &
             + ksupuijf(:,:,nsmin:nsmax)*bsupvijf(:,:,nsmin:nsmax)             &
             - ksupvijf(:,:,nsmin:nsmax)*bsupuijf(:,:,nsmin:nsmax)
      KxBuij = KxBuij                                                          &
             + ksupvijf(:,:,nsmin:nsmax)*bsupsijf(:,:,nsmin:nsmax)             &
             - ksupsijf(:,:,nsmin:nsmax)*bsupvijf(:,:,nsmin:nsmax)
      KxBvij = KxBvij                                                          &
             + ksupsijf(:,:,nsmin:nsmax)*bsupuijf(:,:,nsmin:nsmax)             &
             - ksupuijf(:,:,nsmin:nsmax)*bsupsijf(:,:,nsmin:nsmax)
      
      END SUBROUTINE

!-------------------------------------------------------------------------------
!>  @brief Compute covariant (sub) components of the MHD force J X B - grad(p).
!>
!>  The MHD forces are sum of Lorentz and pressure gradient forces
!>
!>     F = J X B - Grad(p)                                                 (1)
!>
!>  where J is the plasma current density. In curvilinear coodinates, the forces
!>  components becomes
!>
!>     F_s = (K^u*B^v - K^v*B^u) - dp/ds                                   (2a)
!>
!>     F_u = (K^v*B^s - K^s*B^v) - dp/du                                   (2b)
!>
!>     F_v = (K^s*B^u - K^u*B^s) - dp/dv                                   (2c)
!>
!>  where K^i = sqrt(g)*J^i . 
!>  Note this routine works on one parity at a time. When l_linearize is true,
!>  calculate the perturbed Force only (Ignores F(xc=0) part). When l_linearize
!>  is false, calculate full non-linear force(F(xc=0) part only)
!>
!>
!>  @param[inout] fsubsmnf Covariant force in the s direction.
!>  @param[inout] fsubumnf Covariant force in the s direction.
!>  @param[inout] fsubvmnf Covariant force in the s direction.
!>  @param[in]    pijh     Real space pressure.
!>  @param[inout] KxBsij   Covariant component of the lorentz force in the s
!>                         direction.
!>  @param[inout] KxBuij   Covariant component of the lorentz force in the u
!>                         direction.
!>  @param[inout] KxBvij   Covariant component of the lorentz force in the v
!>                         direction.
!>  @param[inout] pardamp  Parallel damping.
!>  @param[in]    parity   Parity of the fourier components.
!-------------------------------------------------------------------------------
      SUBROUTINE GetMHDForce(fsubsmnf, fsubumnf, fsubvmnf, pijh,               &
                             KxBsij, KxBuij, KxBvij, pardamp, parity)
      USE island_params, ONLY: ohs=>ohs_i, fourier_context
      USE utilities, ONLY: set_bndy_half_to_full, set_bndy_half_to_full_ds,    &
                           set_bndy_full_origin, set_bndy_fouier_m0
      USE fourier, ONLY: f_none

      IMPLICIT NONE

!  Declare Arguments
      REAL (dp), DIMENSION(:,:,:), INTENT(inout)            :: fsubsmnf
      REAL (dp), DIMENSION(:,:,:), INTENT(inout)            :: fsubumnf
      REAL (dp), DIMENSION(:,:,:), INTENT(inout)            :: fsubvmnf
      REAL (dp), DIMENSION(:,:,:), ALLOCATABLE, INTENT(in)  :: pijh
      REAL (dp), DIMENSION(:,:,:), ALLOCATABLE, INTENT(in)  :: KxBsij
      REAL (dp), DIMENSION(:,:,:), ALLOCATABLE, INTENT(in)  :: KxBuij
      REAL (dp), DIMENSION(:,:,:), ALLOCATABLE, INTENT(in)  :: KxBvij
      REAL (dp), DIMENSION(:,:), ALLOCATABLE, INTENT(inout) :: pardamp
      INTEGER, INTENT(in)                                   :: parity

!  local variables
      INTEGER                                               :: fours
      INTEGER                                               :: fouruv
      INTEGER                                               :: sparity
      INTEGER                                               :: m
      INTEGER                                               :: n
      INTEGER                                               :: n_mode
      INTEGER                                               :: moff
      INTEGER                                               :: noff
      REAL (dp), DIMENSION(:,:,:), ALLOCATABLE              :: pmnh
      REAL (dp), DIMENSION(:,:,:), ALLOCATABLE              :: pmnf
      REAL (dp), DIMENSION(:,:,:), ALLOCATABLE              :: pmnf_ds
      INTEGER                                               :: nsmin
      INTEGER                                               :: nsmax
      INTEGER                                               :: istat

!  Start of executable code
      nsmax = MIN(ns, endglobrow + 1)
      nsmin = MAX(1, startglobrow)

      ALLOCATE(pmnf_ds(0:mpol,-ntor:ntor,nsmin:nsmax),                         &
               pmnf(0:mpol,-ntor:ntor,nsmin:nsmax),                            &
               pmnh(0:mpol,-ntor:ntor,nsmin:nsmax), stat=istat)
      CALL assert_eq(0, istat, 'Allocation failed before GetMHDForce')

      nsmax = MIN(ns, endglobrow)

      IF (parity .EQ. f_cos) THEN
         fours = f_cos
         fouruv = f_sin
         sparity = -1
      ELSE
         fours = f_sin
         fouruv = f_cos
         sparity = 1
      END IF

!  Harmonics of covariant Lorentz force components on full grid.
      CALL fourier_context%tomnsp(KxBsij, fsubsmnf(:,:,nsmin:nsmax), fours)
      CALL fourier_context%tomnsp(KxBuij, fsubumnf(:,:,nsmin:nsmax), fouruv)
      CALL fourier_context%tomnsp(KxBvij, fsubvmnf(:,:,nsmin:nsmax), fouruv)

      IF (nsmin .eq. 1) THEN
         CALL set_bndy_full_origin(fsubsmnf, fsubumnf, fsubvmnf, f_none)
      END IF

!  Harmonics of true pressure. The jacobian factor has been divided out of pijh.
      CALL fourier_context%tomnsp(pijh, pmnh, fours)
      CALL to_full_mesh(pmnh, pmnf)
      CALL set_bndy_half_to_full(pmnh, pmnf, nsmin, fours, f_none)

!  Radial derivative of pressure.
      CALL GradientFull(pmnf_ds, pmnh)
      CALL set_bndy_half_to_full_ds(pmnh, pmnf_ds, nsmin, fours, f_none)

      IF (nsmax .eq. ns .and. l_pedge .and. .not.l_vessel) THEN
         pmnf(:,:,ns) = 0
         pmnf(m0,n0,ns) = pmnh(m0,n0,ns)
      END IF

!  Add pressure gradient to the lorentz force.
      fsubsmnf(:,:,nsmin:nsmax) = fsubsmnf(:,:,nsmin:nsmax) - pmnf_ds(:,:,nsmin:nsmax)
      DO m = 0, mpol
         moff = m + LBOUND(fsubumnf, 1)
         fsubumnf(moff,:,nsmin:nsmax) = fsubumnf(moff,:,nsmin:nsmax)           &
                                      - m*sparity*pmnf(m,:,nsmin:nsmax)
      END DO
      DO n = -ntor, ntor
         n_mode = fourier_context%tor_modes(n)
         noff = n + ntor + LBOUND(fsubvmnf, 2)
         fsubvmnf(:,noff,nsmin:nsmax) = fsubvmnf(:,noff,nsmin:nsmax)           &
     &                                - (n_mode*nfp*sparity *                  &
     &                                   pmnf(:,n,nsmin:nsmax))
      END DO

      DEALLOCATE(pmnf, pmnf_ds, pmnh)

      CALL set_bndy_fouier_m0(fsubsmnf, fsubumnf, fsubvmnf, parity)
      IF (nsmin .eq. 1) THEN
         CALL set_bndy_full_origin(fsubsmnf, fsubumnf, fsubvmnf, f_none)
      END IF

      CALL get_force_harmonics(pardamp, fsubsmnf, fsubumnf, fsubvmnf, parity)

      END SUBROUTINE

!-------------------------------------------------------------------------------
!>  @brief Compute scaling factors for parallel flow damping.
!>
!>  @param[inout] parscale Scratch array.
!>  @param[inout] pardamp  Parallel damping.
!-------------------------------------------------------------------------------
      SUBROUTINE InitParDamping(parscale, pardamp)
      USE hessian, ONLY:  l_Compute_Hessian, mupar, mupar_norm

!  Declare Arguments
      REAL (dp), DIMENSION(:,:,:), ALLOCATABLE, INTENT(inout) :: parscale
      REAL (dp), DIMENSION(:,:), ALLOCATABLE, INTENT(inout)   :: pardamp

!  local variables
      INTEGER                                                 :: nsmin
      INTEGER                                                 :: nsmax

!  Start of executable code
      IF (muPar .eq. zero .or. .not.l_Compute_Hessian) THEN
         RETURN
      ELSE
         nsmin = MAX(1, startglobrow)
         nsmax = MIN(ns, endglobrow)

         CALL assert_eq(nsmin, LBOUND(parscale, 3), 'PARSCALE LBOUND WRONG!')
         CALL assert_eq(nsmax, UBOUND(parscale, 3), 'PARSCALE UBOUND WRONG!')
         ALLOCATE(pardamp(nsmin:nsmax,3))

!  Update mupar_norm in cv_currents and put in 1/jacobh factor. Simpler and
!  faster approximation, keep only diagonal displacement terms from v dot D.
         muPar = ABS(muPar)

         parscale = bsubsijf(:,:,nsmin:nsmax)*bsubsijf(:,:,nsmin:nsmax)        &
                  / jacobf(:,:,nsmin:nsmax)
         CALL SurfAverage(pardamp(nsmin:nsmax,1), parscale, nsmin, nsmax)
         pardamp(:,1) = pardamp(:,1)*muPar*mupar_norm(nsmin:nsmax)

         parscale = bsubuijf(:,:,nsmin:nsmax)*bsubuijf(:,:,nsmin:nsmax)        &
                  / jacobf(:,:,nsmin:nsmax)
         CALL SurfAverage(pardamp(nsmin:nsmax,2), parscale, nsmin, nsmax)
         pardamp(:,2) = pardamp(:,2)*muPar*mupar_norm(nsmin:nsmax)
      
         parscale = bsubvijf(:,:,nsmin:nsmax)*bsubvijf(:,:,nsmin:nsmax)        &
                  / jacobf(:,:,nsmin:nsmax)
         CALL SurfAverage(pardamp(nsmin:nsmax,3), parscale, nsmin, nsmax)
         pardamp(:,3) = pardamp(:,3)*muPar*mupar_norm(nsmin:nsmax)
      END IF

      END SUBROUTINE

!-------------------------------------------------------------------------------
!>  @brief Final computation of the MHD covariant Fourier force components.
!>
!>  Add pressure gradients and parallel damping (to Hessian), and scale 
!>  and print boundary values. See @ref siesta_init::GetFields for boundary
!>  conditions.
!>
!>  @param[inout] pardamp   Parallel damping.
!>  @param[inout] f_smnf  Covariant component of the mhd force in the s
!>                        direction.
!>  @param[inout] f_umnf  Covariant component of the mhd force in the u
!>                        direction.
!>  @param[inout] f_vmnf  Covariant component of the mhd force in the v
!>                        direction.
!>  @param[in]    parity  Parity of the fourier components.
!-------------------------------------------------------------------------------
      SUBROUTINE get_force_harmonics(pardamp, f_smnf, f_umnf, f_vmnf, parity)
      USE hessian, ONLY: muPar, l_Compute_Hessian
      USE shared_data, ONLY: nprecon, l_PrintOriginForces,                     &
                             l_push_s, l_push_u, l_push_v,                     &
                             l_push_edge, col_scale
      USE nscalingtools, ONLY: FIRSTPARTITION_GE_2, LASTPARTITION_GE_2
!  These were computed in siesta_displacment and are used for || damping
      USE quantities, ONLY: gvsupsmncf => fsupsmncf, gvsupumnsf => fsupumnsf,  &
                            gvsupvmnsf => fsupvmnsf, gvsupsmnsf => fsupsmnsf,  &
                            gvsupumncf => fsupumncf, gvsupvmncf => fsupvmncf

      IMPLICIT NONE

!  Declare Arguments
      REAL (dp), DIMENSION(:,:), ALLOCATABLE, INTENT(inout) :: pardamp
      REAL (dp), DIMENSION(:,:,:), INTENT(inout)            :: f_smnf
      REAL (dp), DIMENSION(:,:,:), INTENT(inout)            :: f_umnf
      REAL (dp), DIMENSION(:,:,:), INTENT(inout)            :: f_vmnf
      INTEGER, INTENT(in)                                   :: parity

!  local variables
      REAL (dp)                                             :: skston
      REAL (dp)                                             :: skstoff
      REAL (dp), DIMENSION(10)                              :: tmps
      INTEGER                                               :: nsmin
      INTEGER                                               :: nsmax
      INTEGER                                               :: moff
      INTEGER                                               :: noff

!  local parameters
      REAL (dp), PARAMETER                                  :: p5 = 0.5_dp

!  Start of executable code
      CALL second0(skston)

      nsmin = MAX(1, startglobrow)
      nsmax = MIN(ns,endglobrow)

      moff = 1
      noff = 1 + ntor

      CALL assert(ALL(f_smnf(m0+moff,-ntor+noff:noff-1,nsmin:nsmax) .eq. zero),&
                  'f_smnf != 0, m=0, n<0')
      CALL assert(ALL(f_umnf(m0+moff,-ntor+noff:noff-1,nsmin:nsmax) .eq. zero),&
                  'f_umnf != 0, m=0, n<0')
      CALL assert(ALL(f_vmnf(m0+moff,-ntor+noff:noff-1,nsmin:nsmax) .eq. zero),&
                  'f_vmnf != 0, m=0, n<0')
      IF (ALLOCATED(pardamp)) THEN
         IF (parity .eq. f_cos) THEN
            CALL AddParDamping(pardamp,                                        &
                               f_smnf(:,:,nsmin:nsmax),                        &
                               f_umnf(:,:,nsmin:nsmax),                        &
                               f_vmnf(:,:,nsmin:nsmax),                        &
                               gvsupsmncf(:,:,nsmin:nsmax),                    &
                               gvsupumnsf(:,:,nsmin:nsmax),                    &
                               gvsupvmnsf(:,:,nsmin:nsmax))
            CALL assert(ALL(f_umnf(m0+moff,noff,:) .eq. zero),                 &
                        'f_umnf(m0,0) != 0')
            CALL assert(ALL(f_vmnf(m0+moff,noff,:) .eq. zero),                 &
                        'f_vmnf(m0,0) != 0')
         ELSE
            CALL AddParDamping(pardamp,                                        &
                               f_smnf(:,:,nsmin:nsmax),                        &
                               f_umnf(:,:,nsmin:nsmax),                        &
                               f_vmnf(:,:,nsmin:nsmax),                        &
                               gvsupsmnsf(:,:,nsmin:nsmax),                    &
                               gvsupumncf(:,:,nsmin:nsmax),                    &
                               gvsupvmncf(:,:,nsmin:nsmax))
            CALL ASSERT(ALL(f_smnf(m0+moff,noff,:) .eq. zero),                 &
                        'f_smnf(m0,0) != 0')
         END IF
      END IF

!  GATHER BDY FORCES AND PRINT OUT ON PROC=0
      PRINT_O: IF (l_PrintOriginForces) THEN
         tmps = 0
!  Assume s=0 is on the 0th processor.
         IF (nsmin .le. 2 .and. nsmax .ge. 2) THEN
            tmps(1) = SQRT(SUM(f_smnf(m1+moff,:,2)**2))
            tmps(2) = SQRT(SUM(f_smnf(m0+moff,:,2)**2) +                       &
                           SUM(f_smnf(m2+moff:,:,2)**2))
            tmps(3) = SQRT(SUM(f_umnf(m1+moff,:,2)**2))
            tmps(4) = SQRT(SUM(f_umnf(m0+moff,:,2)**2) +                       &
                           SUM(f_umnf(m2+moff:,:,2)**2))
            tmps(5) = SQRT(SUM(f_vmnf(m0+moff,:,2)**2))
            tmps(6) = SQRT(SUM(f_vmnf(m1+moff:,:,2)**2))
         END IF

         IF (nsmin .le. ns - 1 .and. nsmax .ge. ns - 1) THEN
            tmps(7) = SQRT(SUM(f_umnf(:,:,ns-1)**2))
            tmps(8) = SQRT(SUM(f_vmnf(:,:,ns-1)**2))
         END IF

         IF (nsmax .eq. ns) THEN
            tmps(9) = SQRT(SUM(f_umnf(:,:,ns)**2))
            tmps(10) = SQRT(SUM(f_vmnf(:,:,ns)**2))
         END IF

#if defined(MPI_OPT)
         IF (PARSOLVER) THEN
            CALL MPI_ALLREDUCE(MPI_IN_PLACE, tmps, 10, MPI_REAL8, MPI_SUM,     &
                               SIESTA_COMM, MPI_ERR)
         END IF
#endif

         IF (iam .eq. 0 .and. lverbose) THEN
            IF (parity .eq. f_cos) THEN
               WRITE (*,1000)
            ELSE
               WRITE (*,1001)
            END IF
            WRITE (*,1002) SQRT(SUM(f_smnf(m1+moff,:,1)**2)), tmps(1), tmps(2)
            WRITE (*,1003) SQRT(SUM(f_umnf(m1+moff,:,1)**2)), tmps(3), tmps(4)
            WRITE (*,1004) SQRT(SUM(f_vmnf(m0+moff,:,1)**2)), tmps(5), tmps(6)
            WRITE (*,1005) tmps(9), tmps(7)
            WRITE (*,1006) tmps(10), tmps(8)
         END IF
      END IF PRINT_O

!  Origin boundary conditions. 1/2 factors below are needed to preserve hessian
!  symmetry.
       IF (nsmin .eq. 1) THEN

          IF (.not.l_push_s) THEN
             f_smnf(m1 + moff,:,1) = 0
          ELSE
             f_smnf(m1 + moff,:,1) = p5*f_smnf(m1 + moff,:,1)
          END IF

          f_umnf(:,:,1) = 0

          IF (.not.l_push_v) THEN
             f_vmnf(m0 + moff,:,1) = 0
          ELSE
             f_vmnf(m0 + moff,:,1) = p5*f_vmnf(m0 + moff,:,1)
          END IF
       END IF

!  Edge boundary conditions.
       IF (nsmax .eq. ns) THEN
          f_smnf(:,:,ns) = 0
          IF (l_pedge .or. .not.l_push_edge .or. l_vessel) THEN
             f_umnf(:,:,ns) = 0
          ELSE
             f_umnf(:,:,ns) = p5*f_umnf(:,:,ns)
          END IF
          IF (.not.l_push_edge .or. l_vessel) THEN
             f_vmnf(:,:,ns) = 0
          ELSE
             f_vmnf(:,:,ns) = p5*f_vmnf(:,:,ns)
          END IF
       END IF

       CALL second0(skstoff)
       get_force_harmonics_time = get_force_harmonics_time + (skstoff - skston)

1000  FORMAT(/,1x,'SYMMETRIC FORCES: ')
1001  FORMAT(/,1x,'ASYMMETRIC FORCES: ')
1002  FORMAT(' fs(1,m=1): ',1p,e12.3,                                          &
             ' fs(2,m=1): ',1pe12.3,' fs(2,m!=1):',1pe12.3)
1003  FORMAT(' fu(1,m=1): ',1p,e12.3,                                          &
             ' fu(2,m=1): ',1pe12.3,' fu(2,m!=1):',1pe12.3)
1004  FORMAT(' fv(1,m=0): ',1p,e12.3,                                          &
             ' fv(2,m=0): ',1pe12.3,' fv(2,m>0): ',1pe12.3)
1005  FORMAT(' fu(ns-1) : ',1p,e12.3,' fu(ns)   : ',1pe12.3)
1006  FORMAT(' fv(ns-1) : ',1p,e12.3,' fv(ns)   : ',1pe12.3,/)

      END SUBROUTINE

!-------------------------------------------------------------------------------
!>  @brief Add parallel flow damping terms to forces (for Hessian calculation)
!>         to suppress null space at resonances.
!>
!>  Compute parallel (to B) damping force and add it to (linearized) forces to
!>  suppress null space in preconditioner.
!>
!>  @param[inout] f_smnf    Covariant component of the mhd force in the s
!>                          direction.
!>  @param[inout] f_umnf    Covariant component of the mhd force in the u
!>                          direction.
!>  @param[inout] f_vmnf    Covariant component of the mhd force in the v
!>                          direction.
!>  @param[in]    jvsupsmnf Contravariant displacement in the s direction.
!>  @param[in]    jvsupsmnf Contravariant displacement in the u direction.
!>  @param[in]    jvsupsmnf Contravariant displacement in the v direction.
!-------------------------------------------------------------------------------
      SUBROUTINE addpardamping(pardamp, f_smnf, f_umnf, f_vmnf,                &
                               jvsupsmnf, jvsupumnf, jvsupvmnf)
      USE shared_data, ONLY:  col_scale

      IMPLICIT NONE

!  Declare Arguments
      REAL (dp), DIMENSION(:,:), INTENT(in)      :: pardamp
      REAL (dp), DIMENSION(:,:,:), INTENT(inout) :: f_smnf
      REAL (dp), DIMENSION(:,:,:), INTENT(inout) :: f_umnf
      REAL (dp), DIMENSION(:,:,:), INTENT(inout) :: f_vmnf
      REAL (dp), DIMENSION(:,:,:), INTENT(in)    :: jvsupsmnf
      REAL (dp), DIMENSION(:,:,:), INTENT(in)    :: jvsupumnf
      REAL (dp), DIMENSION(:,:,:), INTENT(in)    :: jvsupvmnf

!  local variables
      INTEGER                                    :: js

!  Start of executable code
!  Update mupar_norm in cv_currents and put in 1/jacobh factor. Simpler, faster
!  approcimately diagonal form.
      DO js = 1, SIZE(pardamp,1)
         f_smnf(:,:,js) = f_smnf(:,:,js) + jvsupsmnf(:,:,js)*pardamp(js,1)
         f_umnf(:,:,js) = f_umnf(:,:,js) + jvsupumnf(:,:,js)*pardamp(js,2)
         f_vmnf(:,:,js) = f_vmnf(:,:,js) + jvsupvmnf(:,:,js)*pardamp(js,3)
      END DO

      END SUBROUTINE

      END MODULE
