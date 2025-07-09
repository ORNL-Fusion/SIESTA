!*******************************************************************************
!>  @file siesta_bfield.f90
!>  @brief Contains module @ref siesta_bfield.
!
!  Note separating the Doxygen comment block here so detailed decription is
!  found in the Module not the file.
!
!>  This file contains subroutines for updating half-grid magnetic fields.
!*******************************************************************************
      MODULE siesta_bfield
      USE v3_utilities, ONLY: assert
      USE stel_kinds
      USE stel_constants
      USE nscalingtools, ONLY: startglobrow, endglobrow
      USE timer_mod
      USE fourier

      IMPLICIT NONE

      CONTAINS

!*******************************************************************************
!  UTILITY SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  @brief Update contravariant componets of the magnetic field.
!>
!>  Advances half-mesh magnetic field components in nm space from time t to
!>  t + dt. Used to update values of bsup*mnh and bsupijh stored in the
!>  @ref quantities module. Discretized version of Faraday's law
!>
!>    dB/dt = -Curl(E)                                                       (1)
!>
!>  @param[in] l_add_res Control flag to add resonant field.
!-------------------------------------------------------------------------------
      SUBROUTINE update_bfield(l_add_res)
      USE quantities
      USE descriptor_mod, ONLY: INHESSIAN
      USE siesta_namelist, ONLY: lresistive, eta_factor
      USE shared_data, ONLY: lasym, fsq_res, fsq_total, buv_res
USE utilities

      IMPLICIT NONE

!  Declare Arguments
      LOGICAL, INTENT(in)                      :: l_add_res

!  local variables
      INTEGER                                  :: istat
      INTEGER                                  :: js
      INTEGER                                  :: nsmin
      INTEGER                                  :: nsmax
      REAL (dp), DIMENSION(:,:,:), ALLOCATABLE :: resistivity
      REAL (dp)                                :: ton
      REAL (dp)                                :: toff
      REAL (dp)                                :: rho
      REAL (dp)                                :: delt_cfl
      REAL (dp)                                :: eta_prof
      REAL (dp), DIMENSION(:,:,:), ALLOCATABLE :: esubsijf
      REAL (dp), DIMENSION(:,:,:), ALLOCATABLE :: esubuijf
      REAL (dp), DIMENSION(:,:,:), ALLOCATABLE :: esubvijf
      REAL (dp), DIMENSION(:,:,:), ALLOCATABLE :: TEMP
      REAL (dp), DIMENSION(:,:,:), POINTER     :: esubsmnf
      REAL (dp), DIMENSION(:,:,:), POINTER     :: esubumnf
      REAL (dp), DIMENSION(:,:,:), POINTER     :: esubvmnf

!  Start of executable code
      CALL second0(ton)

      nsmin = MAX(1, startglobrow - 1)
      nsmax = MIN(ns, endglobrow + 1)

!  Calculate covariant components of (non-resistive) electric field E = -v X B
!  on full-space mesh.
      ALLOCATE(esubsijf(ntheta,nzeta,nsmin:nsmax),                             &
               esubuijf(ntheta,nzeta,nsmin:nsmax),                             &
               esubvijf(ntheta,nzeta,nsmin:nsmax), stat=istat)
      CALL assert_eq(0, istat, 'Allocation1 failed in UPDATE_BFIELD')

!  Compute electric field.
!
!   E = V X B
!
!  Note V contains a jacobian factor so E has none since it cancels out with the
!  cross product.
!  -E_s = V^uB^v - V^vB^u : jv -> jac*v
      esubsijf = -(jvsupuijf(:,:,nsmin:nsmax)*bsupvijf0(:,:,nsmin:nsmax) -     &
                   jvsupvijf(:,:,nsmin:nsmax)*bsupuijf0(:,:,nsmin:nsmax))
!  -E_u = V^vB^s - V^sB^v
      esubuijf = -(jvsupvijf(:,:,nsmin:nsmax)*bsupsijf0(:,:,nsmin:nsmax) -     &
                   jvsupsijf(:,:,nsmin:nsmax)*bsupvijf0(:,:,nsmin:nsmax))
!  -E_v = V^sB^u - V^uB^s
      esubvijf = -(jvsupsijf(:,:,nsmin:nsmax)*bsupuijf0(:,:,nsmin:nsmax) -     &
                   jvsupuijf(:,:,nsmin:nsmax)*bsupsijf0(:,:,nsmin:nsmax))

!  Verify boundary condition: esubu,v(s=1) = 0 (tangential E vanishes at bdy =>
!  dB^s = 0). FIXME: Does this still hold in free boundary?
      IF (nsmax .eq. ns) THEN
         CALL assert(ALL(esubuijf(:,:,ns).EQ.zero),                            &
                     'esubuijf(ns) != 0 in UPDATE_BFIELD')
      END IF
      IF (nsmax .eq. ns) THEN
         CALL assert(ALL(esubvijf(:,:,ns).EQ.zero),                            &
                     'esubvijf(ns) != 0 in UPDATE_BFIELD')
      END IF

!  Note this will LOWER the energy due to eta*|J|||**2 heating
!
!    esubX(resistive) = eta(JdotB/B^2)*BsubX
!
!  so the magnetic energy decreased due to this term. Note ksubX=JsubX are the
!  covariant components of the current.
      IF (lresistive .and. (l_add_res .or. ALLOCATED(buv_res))) THEN
         delt_cfl = hs_i*hs_i*ABS(eta_factor)
         IF (fsq_total .lt. fsq_res) THEN
            delt_cfl = delt_cfl*SQRT(fsq_total/fsq_res)
         END IF

         ALLOCATE(resistivity(ntheta,nzeta,nsmin:nsmax), stat=istat)
         CALL assert_eq(0, istat, 'Allocation2 failed in update_bfield')

         DO js = nsmin, nsmax
            rho = hs_i*(js - 1)
            eta_prof = rho*rho*(1 - rho)
            resistivity(:,:,js) = delt_cfl*eta_prof
         END DO

         IF (ALLOCATED(buv_res)) THEN
            resistivity = resistivity*buv_res(:,:,nsmin:nsmax)
         ELSE
!  Divide out jacobf factor from cv_currents.
            resistivity = resistivity/jacobf(:,:,nsmin:nsmax)
         END IF

!  Isotropic resistivity, E ~ eta*J. When adding the perturbation, K ~ B in
!  init_state.
         esubsijf = esubsijf + resistivity*ksubsijf(:,:,nsmin:nsmax)
         esubuijf = esubuijf + resistivity*ksubuijf(:,:,nsmin:nsmax)
         esubvijf = esubvijf + resistivity*ksubvijf(:,:,nsmin:nsmax)

         DEALLOCATE (resistivity)
      END IF

!  Update Bfield using Faraday's Law
      CALL Faraday(djbsupsmnsh, djbsupumnch, djbsupvmnch,                      &
                   esubsijf, esubuijf, esubvijf, f_sin, nsmin, nsmax)
      IF (lasym) THEN
         CALL Faraday(djbsupsmnch, djbsupumnsh, djbsupvmnsh,                   &
                      esubsijf, esubuijf, esubvijf, f_cos, nsmin, nsmax)
      END IF

      DEALLOCATE(esubsijf, esubuijf, esubvijf, stat=istat)
      CALL assert_eq(0, istat, 'Deallocation failed in update_bfield')

      CALL second0(toff)
      time_update_bfield = time_update_bfield + (toff - ton)

      END SUBROUTINE

!-------------------------------------------------------------------------------
!>  @brief Use Faraday's law dB = -curl(E) to compute magnitic field
!>         perturbation.
!>
!>  The arrays must be passed in as ALLOCATABLE to preserve the array bounds.
!>
!>  @param[inout] djbsupsmnh Contravariant displacement in the s direction.
!>  @param[inout] djbsupsmnh Contravariant displacement in the u direction.
!>  @param[inout] djbsupsmnh Contravariant displacement in the v direction.
!>  @param[in]    esubsijf   Covariant real space electric field in the s
!>                           direction.
!>  @param[in]    esubuijf   Covariant real space electric field in the s
!>                           direction.
!>  @param[in]    esubvijf   Covariant real space electric field in the s
!>                           direction.
!>  @param[in]    parity     Fourier parity control flag.
!>  @param[in]    nsmin      Minimum radial index.
!>  @param[in]    nsmax      Maximum radial index.
!-------------------------------------------------------------------------------
      SUBROUTINE Faraday(djbsupsmnh, djbsupumnh, djbsupvmnh,                   &
                         esubsijf, esubuijf, esubvijf, parity, nsmin, nsmax)
      USE shared_data, ONLY: l_natural, delta_t
      USE hessian, ONLY: l_compute_hessian
      USE utilities, ONLY: Curl_FtoH, set_bndy_full_origin, set_bndy_fouier_m0
      USE island_params, ONLY: mpol => mpol_i, ntor => ntor_i, fourier_context

      IMPLICIT NONE

!  Declare Arguments
      REAL (dp), DIMENSION(:,:,:), ALLOCATABLE, INTENT(inout) :: djbsupsmnh
      REAL (dp), DIMENSION(:,:,:), ALLOCATABLE, INTENT(inout) :: djbsupumnh
      REAL (dp), DIMENSION(:,:,:), ALLOCATABLE, INTENT(inout) :: djbsupvmnh
      REAL (dp), DIMENSION(:,:,:), INTENT(in)                 :: esubsijf
      REAL (dp), DIMENSION(:,:,:), INTENT(in)                 :: esubuijf
      REAL (dp), DIMENSION(:,:,:), INTENT(in)                 :: esubvijf
      INTEGER, INTENT(in)                                     :: parity
      INTEGER, INTENT(in)                                     :: nsmin
      INTEGER, INTENT(in)                                     :: nsmax

!  local variables
      INTEGER                                                 :: fours
      INTEGER                                                 :: fouruv
      REAL (dp)                                               :: sparity
      REAL (dp)                                               :: r12
      REAL (dp), DIMENSION(:,:,:), POINTER                    :: esubsmnf
      REAL (dp), DIMENSION(:,:,:), POINTER                    :: esubumnf
      REAL (dp), DIMENSION(:,:,:), POINTER                    :: esubvmnf

!  Start of executable code
      IF (parity .EQ. f_sin) THEN
!  e_s (sin), e_u (cos), e_v (cos)
         fours = f_sin
         fouruv = f_cos
         sparity = 1
      ELSE
!  e_s (cos), e_u (sin), e_v (sin)
         fours = f_cos
         fouruv = f_sin
         sparity = -1
      END IF

!  Allocate working arrays used in Faraday. FIXME: Move to inside faraday.
      ALLOCATE(esubsmnf(0:mpol,-ntor:ntor,nsmin:nsmax),                        &
               esubumnf(0:mpol,-ntor:ntor,nsmin:nsmax),                        &
               esubvmnf(0:mpol,-ntor:ntor,nsmin:nsmax))

!  Calculate harmonics of electric field
      CALL fourier_context%tomnsp(esubsijf, esubsmnf, fours)
      CALL fourier_context%tomnsp(esubuijf, esubumnf, fouruv)
      CALL fourier_context%tomnsp(esubvijf, esubvmnf, fouruv)

!  Impose boundary consition at first half-grid point.
      IF (startglobrow .EQ. 1) THEN
!  These conditions are needed so that delta-W jogs (CheckForces) agree with
!  forces at origin.
         CALL set_bndy_full_origin(esubsmnf, esubumnf, esubvmnf, f_none)

!         IF (.not.l_natural) THEN
!            r12 = sparity*hs_i/2
!            djbsupsmnh(m1,:,2) = (esubsmnf(m1,:,1) + esubsmnf(m1,:,2))/2
!            djbsupumnh(m1,:,2) = (esubumnf(m1,:,1) + esubumnf(m1,:,2))/2
!            djbsupsmnh(m1,:,1) = r12*djbsupsmnh(m1,:,2) - djbsupumnh(m1,:,2) !-> 0
!  This constrains [esubs(2)*r12 - esubu(2)] = 0 at s=r12 (first half grid pt)
!  needed to make djb^s ~ r12*djb^u there
!            esubsmnf(m1,:,1) = esubsmnf(m1,:,1) - djbsupsmnh(m1,:,1)/r12
!            esubumnf(m1,:,1) = esubumnf(m1,:,1) + djbsupsmnh(m1,:,1)
!         END IF
      END IF

!  dB = -Curl(E)
      CALL curl_ftoh(esubsmnf, esubumnf, esubvmnf,                             &
                     djbsupsmnh, djbsupumnh, djbsupvmnh,                       &
                     parity, nsmin, nsmax)

      DEALLOCATE(esubsmnf, esubumnf, esubvmnf)

!  Boundary conditions at origin (Used to compute bfields at origin in
!  siesta_init and siesta_forces). See siesta_init::GetFields for boundary
!  conditions.
      IF (nsmin .eq. 1) THEN
         djbsupsmnh(:,:,1) = djbsupsmnh(:,:,2)
         djbsupumnh(:,:,1) = djbsupumnh(:,:,2)
         djbsupvmnh(:,:,1) = djbsupvmnh(:,:,2)

         CALL set_bndy_full_origin(djbsupsmnh, djbsupumnh, djbsupvmnh, f_jac)
      END IF
      CALL set_bndy_fouier_m0(djbsupsmnh, djbsupumnh, djbsupvmnh, parity)

!  Calculate increment of the magnetic field harmonics: use Euler scheme with dt
!  given by vsub*mn advance equations. Results computed here for half-grid
!  perturbations in mn space are used in the calculation of the force and to
!  advance the B's in UPDATE_STATE (in siesta_state module). Minus sign since
!  Faraday is dB = -Curl(E).
      djbsupsmnh = -delta_t*djbsupsmnh
      djbsupumnh = -delta_t*djbsupumnh
      djbsupvmnh = -delta_t*djbsupvmnh

      END SUBROUTINE

      END MODULE
