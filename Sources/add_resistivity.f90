!>  @brief updates resistive contribution to B with resistive currents
      SUBROUTINE add_resistive_E
      USE stel_kinds
      USE shared_data, ONLY: l_getfsq, l_init_state, lverbose,                 &
                             l_update_state, unit_out, nprecon,                &
                             fsq_total1, xc, gc
      USE shared_functions, ONLY: funct_island
      USE siesta_namelist, ONLY: ftol, eta_factor
      USE descriptor_mod, ONLY: iam
      USE siesta_bfield, ONLY: update_bfield
      USE siesta_displacement, ONLY: update_upperv
      USE siesta_init, ONLY: init_state
      USE siesta_state, ONLY: update_state
      IMPLICIT NONE
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      INTEGER, PARAMETER   :: nits = 10
      REAL(dp), PARAMETER  :: zero = 0
      REAL(dp)             :: save_eta, fmhd, eps
      INTEGER              :: itime, j
      LOGICAL, PARAMETER   :: lcurrent_only=.TRUE.
      LOGICAL              :: lverb_save
      
      INTEGER, SAVE        :: nCheck = 0
!-----------------------------------------------
      IF (fsq_total1 .LE. ftol) RETURN

      lverb_save = lverbose
      lverbose = .FALSE.
      
      fmhd = fsq_total1

      l_update_state = .FALSE.  
      
!     No perturbation, vsupX = 0
      xc = 0
      CALL init_state(.FALSE.)
      CALL update_upperv

      save_eta = eta_factor
      eta_factor = eta_factor/MAX(1,nits)

!     Diffuse B-field (eta*j) from last update_state call
      DO itime=1,nits

!  Use new B-field to get KsubXij currents
         CALL init_state (lcurrent_only)
!  Resistive update to B-field
         CALL update_bfield (.TRUE.)
!  Updates B-field with resistive piece
         CALL update_state (.FALSE., zero, zero)

      END DO

      nCheck = nCheck+1
!set iteration value > 0 at which force check is done
      IF (nCheck .EQ. 0) THEN
         CALL CheckForces(xc, gc)
      END IF

      l_getfsq = .TRUE.
      l_init_state = .TRUE.
      CALL funct_island
      
      
      eta_factor = save_eta
      lverbose = lverb_save
      
      IF (iam .EQ. 0) THEN
         DO j = 6, unit_out, unit_out-6
         IF (.NOT.lverbose .AND. j.EQ.6) CYCLE
            WRITE (j,'(/,a,i3)') ' UPDATING RESISTIVE E-FIELD: ITERATIONS=',nits
            WRITE (j, '(a,1p2e12.3)')                                          &
            ' JUMP IN FSQ DUE TO RESISTIVE DIFFUSION: ', fmhd, fsq_total1
         END DO
      END IF

!SPH101116: limit jump next time
      IF (fsq_total1 .GT. 4*fmhd) THEN
         eta_factor = eta_factor*(4*fmhd/fsq_total1)
      ELSE IF (fsq_total1 .LT. 1.1_dp*fmhd) THEN
         eta_factor = eta_factor*(1.1_dp*fmhd/fsq_total1)
      END IF
      
      END SUBROUTINE add_resistive_E
