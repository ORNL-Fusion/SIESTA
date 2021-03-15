!*******************************************************************************
!>  @file siesta_displacement.f90
!>  @brief Contains module @ref siesta_displacement
!
!  Note separating the Doxygen comment block here so detailed decription is
!  found in the Module not the file.
!
!>  Module contains subroutines for computing the real-space contravariant (full
!>  mesh) displacements as part of the SIESTA project. Displacements are defined
!>  by:
!>
!>     J*v*dt                                                                (1)
!>
!>  where J is the coordinate system Jacobian and v is the dispalcement
!>  velocity. The velocity harmonics are obtained from solving the Fourier
!>  components of the MHD force. Fourier harmonics of the displacement velocity
!>  are converted to real space.
!*******************************************************************************
      MODULE siesta_displacement
      USE island_params, ONLY: hs_i, ns=>ns_i
      USE stel_kinds
      USE stel_constants        
      USE nscalingtools, ONLY: startglobrow, endglobrow
      USE timer_mod, ONLY: time_update_upperv

      IMPLICIT NONE

      CONTAINS

!*******************************************************************************
!  UTILITY SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  @brief Initialize updated displacements in real space.
!>
!>  Although sqrt(g) = 0 at origin, jvsups and jvsupv are defined with jacobf
!>  which is finite (=jacobh(2)) at the origin to allow asymptotically correct
!>  variations in the s and v components.
!>
!>  Pertubations to the displacement are stored in the xc array. But the gc_sup
!>  Arrays is used as a temp buffer to scale the displacements. This is why the
!>  gvsup*mn*f components point to fsup*mn*f and not jvsup*mn*f.
!-------------------------------------------------------------------------------
      SUBROUTINE update_upperv
      USE shared_data, ONLY: lasym, gc_sup, xc, col_scale
      USE siesta_state, ONLY: clear_field_perts
      USE quantities, ONLY: gvsupsmncf => fsupsmncf, gvsupumnsf => fsupumnsf,  &
                            gvsupvmnsf => fsupvmnsf, gvsupsmnsf => fsupsmnsf,  &
                            gvsupumncf => fsupumncf, gvsupvmncf => fsupvmncf
      USE fourier, ONLY: f_sin, f_cos

      IMPLICIT NONE

!  local variables
      INTEGER   :: istat
      REAL (dp) :: ton
      REAL (dp) :: toff

!  Start of executable code
      CALL second0(ton)

      CALL Clear_Field_Perts

!  Use components of gc_sup as scratch array. gvsup*mn*f componets point to
!  fsup*mn*f values which in turn point to parts of the gc_sup array.
!  @ref ScaleDisplacement sets gc_sup to the xc array then applies the column
!  factors so the gc_sup array is treated as the xc here.
      CALL ScaleDisplacement(gc_sup, xc, col_scale)

      CALL InitDisplacement(gvsupsmncf, gvsupumnsf, gvsupvmnsf, f_cos)
      IF (lasym) THEN
         CALL InitDisplacement(gvsupsmnsf, gvsupumncf, gvsupvmncf, f_sin)
      END IF

      CALL second0(toff)
      time_update_upperv = time_update_upperv + (toff - ton)

      END SUBROUTINE update_upperv

!-------------------------------------------------------------------------------
!>  @param Initalized variational displacements for a single parity.
!>
!>  Applied boundary controls. See @ref siesta_init::GetFields for boundary
!>  conditions.
!>
!>  @param[inout] gvsupsmnf Contravariant dispalcement in s direction.
!>  @param[inout] gvsupumnf Contravariant dispalcement in u direction.
!>  @param[inout] gvsupvmnf Contravariant dispalcement in v direction.
!>  @param[in]    Fourier iparity   Parity of the displacements.
!-------------------------------------------------------------------------------
      SUBROUTINE InitDisplacement(gvsupsmnf, gvsupumnf, gvsupvmnf, iparity)
      USE fourier, ONLY: f_cos, f_sin, f_none, f_sum, m0, m1, m2, n1, f_jac
      USE quantities, ONLY: jvsupsijf, jvsupuijf, jvsupvijf, mpol, ntor
      USE shared_data, ONLY: l_push_s, l_push_u, l_push_v, l_pedge,     &
                             l_push_edge, col_scale
      USE v3_utilities, ONLY: assert_eq
      USE island_params, ONLY: fourier_context
      USE utilities, ONLY: set_bndy_fouier_m0, set_bndy_full_origin
      USE siesta_namelist, ONLY: l_vessel

      IMPLICIT NONE

!  Declare Arguments
      REAL (dp), DIMENSION(0:mpol,-ntor:ntor,ns), INTENT(inout) :: gvsupsmnf
      REAL (dp), DIMENSION(0:mpol,-ntor:ntor,ns), INTENT(inout) :: gvsupumnf
      REAL (dp), DIMENSION(0:mpol,-ntor:ntor,ns), INTENT(inout) :: gvsupvmnf
      INTEGER, INTENT(in)                                       :: iparity

!  local variables
      INTEGER                                                   :: fours
      INTEGER                                                   :: fouruv
      INTEGER                                                   :: fcomb
      INTEGER                                                   :: nsmin
      INTEGER                                                   :: nsmax

!  Start of executable code
      nsmin = MAX(1, startglobrow - 1)
      nsmax = MIN(endglobrow + 1, ns)

      IF (iparity .eq. f_cos) THEN
         fours = f_cos
         fouruv = f_sin
         fcomb = f_none
      ELSE
         fours = f_sin
         fouruv = f_cos
         fcomb = f_sum
      END IF

      CALL assert_eq(1, LBOUND(gvsupsmnf,3), 'LBOUND WRONG IN InitDisplacement')
      CALL assert_eq(ns, UBOUND(gvsupsmnf,3),                                  &
                     'UBOUND WRONG IN InitDisplacement')

      CALL set_bndy_fouier_m0(gvsupsmnf, gvsupumnf, gvsupvmnf, fours)

!  Origin boundary conditions (evolve m = 1 F_u, F_s, m = 0 F_v), See
!  siesta_init::GetFields for boundary conditions.
      IF (nsmin .eq. 1) THEN
         CALL set_bndy_full_origin(gvsupsmnf, gvsupumnf, gvsupvmnf, f_jac)

!         IF (.not.l_push_s) THEN
!            gvsupsmnf(m1,:,1) = 0
!         END IF
!         IF (.not.l_push_u) THEN
!            gvsupumnf(m0,:,1) = 0
!            gvsupumnf(m1,:,1) = 0
!         END IF
!         IF (.not.l_push_v) THEN
!            gvsupvmnf(m0,:,1) = 0
!         END IF
      END IF

!  Edge boundary conditions (B_tangential = 0 => vsups = 0; J^s = 0 => F_v = 0)
!  Note that F_v and F_u are both ~J^s at edge, so only one can be varied to
!  prevent a dependent row in the Hessian matrix.
!  FIXME: These boundary conditions will only hold for non free boundary cases.
      IF (nsmax .eq. ns) THEN
         gvsupsmnf(:,:,ns) = 0
         IF (l_pedge .or. .not.l_push_edge .or. l_vessel) THEN
            gvsupumnf(:,:,ns) = 0
         END IF
         IF (.not.l_push_edge .or. l_vessel) THEN
            gvsupvmnf(:,:,ns) = 0
         END IF
      END IF

!  Calculate contravariant (SUP) velocities in real space on full mesh
      CALL fourier_context%toijsp(gvsupsmnf(:,:,nsmin:nsmax),                  &
                                  jvsupsijf(:,:,nsmin:nsmax), fcomb, fours)
      CALL fourier_context%toijsp(gvsupumnf(:,:,nsmin:nsmax),                  &
                                  jvsupuijf(:,:,nsmin:nsmax), fcomb, fouruv)
      CALL fourier_context%toijsp(gvsupvmnf(:,:,nsmin:nsmax),                  &
                                  jvsupvijf(:,:,nsmin:nsmax), fcomb, fouruv)

      END SUBROUTINE
      
!-------------------------------------------------------------------------------
!>  @brief Remove columns scaling to get the correct displacements.
!>
!>  @ref update_upperv has the gvsup*mn*f pointing to the quantities pointed to
!>  by the xc_scratch array.
!>
!>  @param[out] xc_scratch Scratch buffer to scale columns in.
!>  @param[in]  xc         Current displacements.
!>  @param[in]  colscale   Column scaling factors.
!-------------------------------------------------------------------------------
      SUBROUTINE ScaleDisplacement(xc_scratch, xc, colscale)
      USE hessian, ONLY: apply_colscale
      USE shared_data, ONLY: mblk_size

      IMPLICIT NONE

!  Declare Arguments
      REAL (dp), DIMENSION(mblk_size,ns), INTENT(out) :: xc_scratch
      REAL (dp), DIMENSION(mblk_size,ns), INTENT(in)  :: xc
      REAL (dp), DIMENSION(mblk_size,ns), INTENT(in)  :: colscale

!  local variables
      INTEGER                                         :: nsmin
      INTEGER                                         :: nsmax

!  Start of executable code
      nsmin = MAX(1, startglobrow - 1)
      nsmax = MIN(endglobrow + 1, ns)

      xc_scratch(:,nsmin:nsmax) = xc(:,nsmin:nsmax)
      CALL Apply_ColScale(xc_scratch, colscale, nsmin, nsmax)

      END SUBROUTINE

      END MODULE
