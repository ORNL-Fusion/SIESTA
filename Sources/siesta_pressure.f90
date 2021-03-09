!*******************************************************************************
!>  @file siesta_pressure.f90
!>  @brief Contains module @ref siesta_pressure
!
!  Note separating the Doxygen comment block here so detailed decription is
!  found in the Module not the file.
!
!>  Updates pressure.
!*******************************************************************************

      MODULE siesta_pressure
      USE v3_utilities, ONLY: assert
      USE stel_kinds
      USE stel_constants
      USE quantities
      USE nscalingtools, ONLY: startglobrow, endglobrow
      USE utilities, ONLY: GradientHalf, to_half_mesh

      IMPLICIT NONE

      CONTAINS

!*******************************************************************************
!  UTILITY SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  @brief advance pressure on the half radial grid from t to t + delta_t.
!>
!>  Changes in pressure are governed by the continuity equation.
!>
!>     dn/dt + Div(nv) = dn/dt + n*Div(v) + v.Grad(n) = 0                    (1)
!>
!>  where n is the density and v is the velocity. From the adiabatic principle,
!>
!>     pV^gamma = Constant                                                   (2)
!>
!>  where p is the pressure and V is the volume. From the ideal gas law V scales
!>  as V ~ 1/n. Using this equation (2) can be rewritten in terms of
!>
!>     p/n^gamma = Constant                                                  (3)
!>
!>  Putting equation (3) into equation (1), the continuty equation can be
!>  written in terms of the pressure.
!>
!>     dp/dt + gamma*p*Div(v) + v.Grad(p) = 0                                (4)
!>
!>  Using the idenity
!>
!>     Div(fA) = f*Div(A) + A.Grad(f)                                        (5)
!>
!>  equation (4) can be written as.
!>
!>     dp/dt = (gamma - 1)v.Grad(p) - gamma*Div(pv)                          (6)
!>
!>  This leads to equation (2.3) in Hirshman et. al. doi:10.1063/1.3597155.
!>  Expanding out equation (6) gives
!>
!>     dp/dt = (gamma - 1)*(v^s*dp/ds + v^u*dp/du + v^v*dp/dv)               (7)
!>           - gamma/J*(dJpv^s/ds + dJpv^u/du + dJpv^v/dv)
!-------------------------------------------------------------------------------
      SUBROUTINE update_pres
      USE timer_mod, ONLY: time_update_pres
      USE fourier, ONLY: f_cos, f_sin, f_sum, f_du, f_dv
      USE quantities, ONLY: gvsupumnsf => fsupumnsf, gvsupvmnsf => fsupvmnsf,  &
                            gvsupumncf => fsupumncf, gvsupvmncf => fsupvmncf
      USE utilities, ONLY: set_bndy_fouier_m0, m0
      USE shared_data, ONLY: delta_t
      USE island_params, ONLY: fourier_context

      IMPLICIT NONE

!  local variables
      INTEGER                                  :: istat
      REAL (dp)                                :: ton
      REAL (dp)                                :: toff
      REAL (dp), DIMENSION(:,:,:), ALLOCATABLE :: jvsupuijh
      REAL (dp), DIMENSION(:,:,:), ALLOCATABLE :: jvsupvijh
      REAL (dp), DIMENSION(:,:,:), ALLOCATABLE :: workij1
      REAL (dp), DIMENSION(:,:,:), ALLOCATABLE :: workij2
      REAL (dp), DIMENSION(:,:,:), ALLOCATABLE :: workij3
      REAL (dp), DIMENSION(:,:,:), ALLOCATABLE :: workmn4
      REAL (dp), DIMENSION(:,:,:), ALLOCATABLE :: workmn5
      REAL (dp), DIMENSION(:,:,:), ALLOCATABLE :: vgradph
      INTEGER                                  :: nsmin
      INTEGER                                  :: nsmax

!  local parameters
!  Does not work well, keep it off.
      LOGICAL, PARAMETER                       :: l_upwind = .FALSE.

!  Start of executable code
      CALL second0(ton)
      nsmin = MAX(1, startglobrow - 1)
      nsmax = MIN(ns, endglobrow + 1)

!  Compute first term ~(gamma-1)V * grad(p) on r.h.s. of pressure equation the
!  u,v derivatives should be consistent with div(pv) term and combine to cancel
!  if div(V) = 0
      ALLOCATE(workij1(ntheta,nzeta,nsmin:nsmax),                              &
               workij2(ntheta,nzeta,nsmin:nsmax),                              &
               workij3(ntheta,nzeta,nsmin:nsmax),                              &
               vgradph(ntheta,nzeta,nsmin:nsmax),  stat=istat)
      CALL assert_eq(0, istat, 'Allocation1 failed in update_pres')

      ALLOCATE(workmn4(0:mpol,-ntor:ntor,nsmin:nsmax),                         &
               workmn5(0:mpol,-ntor:ntor,nsmin:nsmax))
      CALL assert_eq(0, istat, 'Allocation2 failed in update_pres')

      ALLOCATE (jvsupuijh(ntheta,nzeta,nsmin:nsmax),                           &
                jvsupvijh(ntheta,nzeta,nsmin:nsmax), stat=istat)
      CALL assert_eq(0, istat, 'Allocation3 failed in update_pres')

!  On exit, jsupXh valid on [nsmin+1:nsmax]
      CALL to_half_mesh(jvsupuijf, jvsupuijh)
      CALL to_half_mesh(jvsupvijf, jvsupvijh)

!  On exit, vgradph is valid on [nsmin+1:nsmax]. Mixed half/full form: Diagonal
!  in angular derivatives. Exactly cancel Div(pV) term when div(V) = 0.
      workij1 = pijf0_ds(:,:,nsmin:nsmax)*jvsupsijf(:,:,nsmin:nsmax)
      CALL to_half_mesh(workij1, vgradph)
      vgradph = vgradph                                                        &
              + jvsupuijh(:,:,nsmin:nsmax)*pijh0_du(:,:,nsmin:nsmax)           &
              + jvsupvijh(:,:,nsmin:nsmax)*pijh0_dv(:,:,nsmin:nsmax)

!  d/ds Jv^s
      workij2(:,:,nsmin:nsmax) = jvsupsijf(:,:,nsmin:nsmax)

!      IF (nsmin .EQ. 1) THEN
!  Enforce bc at origin in Div(pv) term.
!         workij2(:,:,1) = 0
!      END IF

      CALL GradientHalf(workij1, workij2)

!  Get djvsupudu and djvsupvdv. d/du Jv^u and d/dv Jv^v
      CALL to_half_mesh(gvsupumnsf, workmn4)
      CALL to_half_mesh(gvsupvmnsf, workmn5)

      CALL set_bndy_fouier_m0(workmn4, f_sin)
      CALL set_bndy_fouier_m0(workmn5, f_sin)

      CALL fourier_context%toijsp(workmn4(:,:,nsmin:nsmax),                    &
                                  workij2(:,:,nsmin:nsmax), f_du, f_sin)
      CALL fourier_context%toijsp(workmn5(:,:,nsmin:nsmax),                    &
                                  workij3(:,:,nsmin:nsmax), f_dv, f_sin)
      IF (lasym) THEN
         CALL to_half_mesh(gvsupumncf, workmn4)
         CALL to_half_mesh(gvsupvmncf, workmn5)

         CALL set_bndy_fouier_m0(workmn4, f_cos)
         CALL set_bndy_fouier_m0(workmn5, f_cos)

         CALL fourier_context%toijsp(workmn4(:,:,nsmin:nsmax),                 &
                                     workij2(:,:,nsmin:nsmax),                 &
                                     IOR(f_du, f_sum), f_cos)
         CALL fourier_context%toijsp(workmn5(:,:,nsmin:nsmax),                 &
                                     workij3(:,:,nsmin:nsmax),                 &
                                     IOR(f_dv, f_sum), f_cos)
      END IF

      nsmin = MAX(2, startglobrow)

!  First half of the pressure variation.
!  gamma*p*d/ds Jv^s + gamma*p*d/du Jv^u + gamma*p*d/dv Jv^v
      workij1(:,:,nsmin:nsmax) = gamma*pijh0(:,:,nsmin:nsmax)                  &
                               * (workij1(:,:,nsmin:nsmax) +                   &
                                  workij2(:,:,nsmin:nsmax) +                   &
                                  workij3(:,:,nsmin:nsmax))

!  Second half the the pressure variation.
!  J dp/dt = -Jv^s d/ds p - Jv^u d/du -Jv^v d/dv
!          - gamma*p*d/ds Jv^s + gamma*p*d/du Jv^u + gamma*p*d/dv Jv^v
      workij1(:,:,nsmin:nsmax) = -vgradph(:,:,nsmin:nsmax)                     &
                               - workij1(:,:,nsmin:nsmax)

!  Convert to Fourier space.
      CALL fourier_context%tomnsp(workij1(:,:,nsmin:nsmax),                    &
                                  djpmnch(:,:,nsmin:nsmax), f_cos)
      djpmnch(:,:,nsmin:nsmax) = delta_t*djpmnch(:,:,nsmin:nsmax)
      IF (lasym) THEN
         CALL fourier_context%tomnsp(workij1(:,:,nsmin:nsmax),                 &
                                     djpmnsh(:,:,nsmin:nsmax), f_sin)
         djpmnsh(:,:,nsmin:nsmax) = delta_t*djpmnsh(:,:,nsmin:nsmax)
      END IF

      IF (startglobrow .eq. 1) THEN
         djpmnch(:,:,1) = 0

         IF (lasym) THEN
            djpmnsh(:,:,1) = 0
         END IF
      END IF

      DEALLOCATE(workij1, workij2, workij3, vgradph)
      DEALLOCATE(workmn4, workmn5)

      CALL second0(toff)
      time_update_pres = time_update_pres + (toff - ton)

      END SUBROUTINE update_pres

      END MODULE
