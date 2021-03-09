!*******************************************************************************
!>  @file utilities.f90
!>  @brief Contains module @ref utilities.
!
!  Note separating the Doxygen comment block here so detailed decription is
!  found in the Module not the file.
!
!>  This file contains subroutines for aupdating from t to t + delta_t the
!>  magnetic field and pressure as part of the SIESTA project. Stores updated
!>  values of JPMN*H, JBSUBXMN*H in Quantities Module.
!*******************************************************************************
      MODULE siesta_state
#define DUMP_STATE
      USE stel_kinds
      USE quantities
      USE descriptor_mod, ONLY: iam, SIESTA_COMM
      USE diagnostics_mod
      USE shared_data, ONLY: ste, bs0, bu0, bsbu_ratio_s, bsbu_ratio_a, &
                             unit_out
      USE timer_mod, ONLY: time_update_state
      USE siesta_init, ONLY: init_state
      USE nscalingtools, ONLY: startglobrow, endglobrow, MPI_ERR,       &
                               PARSOLVER
      USE mpi_inc
      USE fourier, ONLY: f_sin, f_cos

      IMPLICIT NONE

!*******************************************************************************
!  Module variables
!*******************************************************************************
!>  Flag to indicate that this was the first run.
      LOGICAL, PUBLIC :: lfirst = .true.

      CONTAINS

!*******************************************************************************
!  UTILITY SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  @brief Update the SIESTA state.
!>
!>  Adds the perturbation to the field. This is called before taking the next
!>  step.
!>
!>  @param[in] lprint    Controls screen output.
!>  @param[in] fsq_total Total force residual.
!>  @param[in] ftol      Force residual tolarance.
!-------------------------------------------------------------------------------
      SUBROUTINE update_state(lprint, fsq_total, ftol)

      IMPLICIT NONE

!  Declare Arguments
      LOGICAL, INTENT(in)   :: lprint
      REAL (dp), INTENT(in) :: fsq_total
      REAL (dp), INTENT(in) :: ftol

!  local variables
      REAL (dp)             :: ton
      REAL (dp)             :: toff
      INTEGER               :: nsmin
      INTEGER               :: nsmax

!  Start of executable code
      CALL second0(ton)

      CALL UpdateFields(jbsupsmnsh,  jbsupumnch,  jbsupvmnch,  jpmnch,         &
                        djbsupsmnsh, djbsupumnch, djbsupvmnch, djpmnch)
      IF (lasym) THEN
         CALL UpdateFields(jbsupsmnch,  jbsupumnsh,  jbsupvmnsh,  jpmnsh,      &
                           djbsupsmnch, djbsupumnsh, djbsupvmnsh, djpmnsh)
      END IF

!  Gather new total fields onto all processors
      IF (PARSOLVER) THEN
         CALL GatherFields
      END IF

      CALL second0(toff)
      time_update_state = time_update_state + (toff - ton)

!  Reset perturbations
      CALL Clear_Field_Perts

      CALL second0(toff)
      time_update_state = time_update_state + (toff - ton)

      IF (lprint) THEN
         CALL second0(ton)

         CALL Update_Diagnostics(jbsupsmnsh, jbsupumnch, jbsupumnch, jpmnch,   &
                                 bs0(1:6), bu0(1:6), bsbu_ratio_s, pwr_spec_s, &
                                 f_sin)
         IF (lasym) THEN
            CALL Update_Diagnostics(jbsupsmnch, jbsupumnsh, jbsupvmnsh,        &
                                    jpmnsh, bs0(7:12), bu0(7:12),              &
                                    bsbu_ratio_a, pwr_spec_a, f_cos)
         END IF
         lfirst = .false.

         nsmin = MAX(1, startglobrow)
         nsmax = MIN(endglobrow, ns)

         CALL init_state(.true.)

!  Compute co-variant and contravariant components of current (times jacobian)
!  need for divJ, BdotJ
         lcurr_init = .true.
         CALL divb(nsmin, nsmax)
         CALL divj(nsmin, nsmax)
         CALL bgradp(nsmin, nsmax)
         CALL tflux
         CALL bdotj_par
         lcurr_init = .false.

         toroidal_flux = toroidal_flux - toroidal_flux0
         IF (iam .eq. 0) THEN
            IF (lverbose) THEN
               WRITE (*,1000) ste(1), ste(2), ste(3), ste(4), divb_rms,        &
                              toroidal_flux, wp/wb, bgradp_rms, max_bgradp,    &
                              min_bgradp, bdotj_rms, bdotj2_rms, divj_rms
            END IF
            WRITE (unit_out,1000) ste(1), ste(2), ste(3), ste(4), divb_rms,    &
                                  toroidal_flux, wp/wb, bgradp_rms,            &
                                  max_bgradp, min_bgradp, bdotj_rms,           &
                                  bdotj2_rms, divj_rms
         END IF

         CALL second0(toff)
         time_update_state = time_update_state + (toff - ton)
      END IF

1000  FORMAT(' SPECTRAL TRUNC ERROR - p: ',1pe11.3,' B_s: ',1pe11.3,           &
             ' B_u: ',1pe11.3,' B_v: ',1pe11.3,/,                              &
             ' DIV-B (rms): ',1pe11.3, ' DEL_TFLUX: ',1pe11.3,/,               &
             ' <BETA>: ', 1pe11.3,' B.GRAD-P (rms): ',1pe11.3,                 &
             ' B.GRAD-P (max): ',1pe11.3,' B.GRAD-P (min): ',1pe11.3,/,        &
             ' (J*B)/|JxB| (rms): ', 1pe11.3,' (J_par)/|J_tot| (rms): ',       &
             1pe11.3,'   DIV-J (rms): ',1pe11.3)

      END SUBROUTINE

!-------------------------------------------------------------------------------
!>  @brief Update the magentic and pressure for the current SIESTA state.
!>
!>  @param[inout] jbsupsmnh  Contravariant magnetic field in the s direction.
!>  @param[inout] jbsupumnh  Contravariant magnetic field in the u direction.
!>  @param[inout] jbsupvmnh  Contravariant magnetic field in the v direction.
!>  @param[inout] jpmnh      Pressure Fourier coeffients.
!>  @param[inout] djbsupsmnh Contravariant magnetic field perturbation in the s
!>                           direction.
!>  @param[inout] djbsupumnh Contravariant magnetic field perturbation in the u
!>                           direction.
!>  @param[inout] djbsupvmnh Contravariant magnetic field perturbation in the v
!>                           direction.
!>  @param[inout] djpmnh     Pressure perturbation.
!-------------------------------------------------------------------------------
      SUBROUTINE UpdateFields(jbsupsmnh, jbsupumnh, jbsupvmnh, jpmnh,          &
                              djbsupsmnh, djbsupumnh, djbsupvmnh, djpmnh)

      IMPLICIT NONE

!  Declare Arguments
      REAL (dp), ALLOCATABLE, DIMENSION(:,:,:), INTENT(inout) :: jbsupsmnh
      REAL (dp), ALLOCATABLE, DIMENSION(:,:,:), INTENT(inout) :: jbsupumnh
      REAL (dp), ALLOCATABLE, DIMENSION(:,:,:), INTENT(inout) :: jbsupvmnh
      REAL (dp), ALLOCATABLE, DIMENSION(:,:,:), INTENT(inout) :: jpmnh
      REAL (dp), ALLOCATABLE, DIMENSION(:,:,:), INTENT(inout) :: djbsupsmnh
      REAL (dp), ALLOCATABLE, DIMENSION(:,:,:), INTENT(inout) :: djbsupumnh
      REAL (dp), ALLOCATABLE, DIMENSION(:,:,:), INTENT(inout) :: djbsupvmnh
      REAL (dp), ALLOCATABLE, DIMENSION(:,:,:), INTENT(inout) :: djpmnh

!  Local Variables
      INTEGER                                                 :: n1
      INTEGER                                                 :: n2

!  Start of executable code.
      n1 = startglobrow
      n2 = endglobrow

      jbsupsmnh(:,:,n1:n2) = jbsupsmnh(:,:,n1:n2) + djbsupsmnh(:,:,n1:n2)
      jbsupumnh(:,:,n1:n2) = jbsupumnh(:,:,n1:n2) + djbsupumnh(:,:,n1:n2)
      jbsupvmnh(:,:,n1:n2) = jbsupvmnh(:,:,n1:n2) + djbsupvmnh(:,:,n1:n2)
      jpmnh(:,:,n1:n2)     = jpmnh(:,:,n1:n2)     + djpmnh(:,:,n1:n2)

      END SUBROUTINE

!-------------------------------------------------------------------------------
!>  @brief Write SIESTA diagnostic (screen) output per iteration.
!>
!>  @param[in]    jbsupsmnh  Contravariant magnetic field in the s direction.
!>  @param[in]    jbsupumnh  Contravariant magnetic field in the u direction.
!>  @param[in]    jbsupvmnh  Contravariant magnetic field in the v direction.
!>  @param[in]    jpmnh      Pressure Fourier coeffients.
!>  @param[out]   bs         FIXME: UNKNOWN
!>  @param[out]   bu         FIXME: UNKNOWN
!>  @param[out]   bsbu_ratio FIXME: UNKNOWN
!>  @param[inout] pwr_spec   FIXME: UNKNOWN
!>  @param[in]    iparity    Fourier parity.
!-------------------------------------------------------------------------------
      SUBROUTINE Update_Diagnostics(jbsupsmnh, jbsupumnh, jbsupvmnh,    &
                                    jpmnh, bs, bu, bsbu_ratio,          &
                                    pwr_spec, iparity)
      USE island_params, ns=>ns_i, hs=>hs_i
      USE fourier, ONLY: m1

      IMPLICIT NONE

!  Declare Arguments
      REAL (dp), DIMENSION(0:mpol,-ntor:ntor,ns), INTENT(IN)      :: jbsupsmnh
      REAL (dp), DIMENSION(0:mpol,-ntor:ntor,ns), INTENT(IN)      :: jbsupumnh
      REAL (dp), DIMENSION(0:mpol,-ntor:ntor,ns), INTENT(IN)      :: jbsupvmnh
      REAL (dp), DIMENSION(0:mpol,-ntor:ntor,ns), INTENT(IN)      :: jpmnh
      REAL (dp), DIMENSION(6), INTENT(out)                        :: bs
      REAL (dp), DIMENSION(6), INTENT(out)                        :: bu
      REAL (dp), INTENT(out)                                      :: bsbu_ratio
      REAL (dp), DIMENSION(0:mpol,-ntor:ntor,ns,4), INTENT(inout) :: pwr_spec
      INTEGER, INTENT(in)                                         :: iparity

!  Local Variables
      REAL (dp)                                                   :: r12
      REAL (dp)                                                   :: diag_b
      REAL (dp)                                                   :: diag_p
      REAL (dp), DIMENSION(12)                                    :: d1
      INTEGER                                                     :: js
      INTEGER                                                     :: itype
      INTEGER                                                     :: n1
      INTEGER                                                     :: n2

!  Local Parameters
      CHARACTER (len=4), DIMENSION(2), PARAMETER                  ::           &
         def = (/'SYM ', 'ASYM'/)

!  Start of executable code.
      n1 = startglobrow
      n2 = endglobrow

      IF (iparity .eq. f_sin) THEN
         itype = 1
         r12 = hs/2
      ELSE
         itype = 2
         r12 = -hs/2
      END IF

#undef _TEST_STATE
!#define _TEST_STATE
#if defined(_TEST_STATE)
      IF (iam .eq. 0) THEN
         WRITE (*,1000) def(itype)
         DO js = -ntor,ntor
            WRITE (*,1001) js, jbsupsmnh(m1,js,2), r12*jbsupumnh(m1,js,2)
         END DO
      END IF
#endif

      bs(1) = SQRT(SUM((jbsupsmnh(m1,:,2) - r12*jbsupumnh(m1,:,2))**2))/r12
      bu(1) = SQRT(SUM((jbsupsmnh(m1,:,2) + r12*jbsupumnh(m1,:,2))**2))/r12
      
#if defined(MPI_OPT)
      IF (PARSOLVER) THEN
         d1(1) = bs(1)
         d1(2) = bu(1)
         CALL MPI_BCAST(d1, 2, MPI_REAL8, 0, SIESTA_COMM, MPI_ERR)
         bs(1) = d1(1)
         bu(1) = d1(2)
      END IF
#endif
      IF (ABS(bu(1)) .GT. 1.E-10_dp) THEN
         bsbu_ratio = bs(1)/bu(1)
         d1 = 0
         DO js=n1, MIN(6, n2)
            d1(js)   = SQRT(SUM(jbsupsmnh(m1,:,js)**2))/ABS(r12)
            d1(js+6) = SQRT(SUM(jbsupumnh(m1,:,js)**2))
         END DO
#if defined(MPI_OPT)
         IF (PARSOLVER) THEN
            CALL MPI_ALLREDUCE(MPI_IN_PLACE, d1, 12, MPI_REAL8, MPI_SUM,       &
                               SIESTA_COMM, MPI_ERR)
         END IF
#endif
         bs(1:6) = d1(1:6)
         bu(1:6) = d1(7:12)
      ELSE 
         bs(1:6) = 0
         bu(1:6) = 0
         bsbu_ratio_s = 1
      END IF

#if defined(DUMP_STATE)
      IF (lfirst) THEN
         pwr_spec(:,:,n1:n2,1) = jbsupsmnh(:,:,n1:n2)
         pwr_spec(:,:,n1:n2,2) = jbsupumnh(:,:,n1:n2)
         pwr_spec(:,:,n1:n2,3) = jbsupvmnh(:,:,n1:n2)
         pwr_spec(:,:,n1:n2,4) = jpmnh(:,:,n1:n2)
      ELSE
         d1(1) = SUM((jbsupsmnh(:,:,n1:n2) - pwr_spec(:,:,n1:n2,1))**2 +       &
                     (jbsupumnh(:,:,n1:n2) - pwr_spec(:,:,n1:n2,2))**2 +       &
                     (jbsupvmnh(:,:,n1:n2) - pwr_spec(:,:,n1:n2,3))**2)
         d1(2) = SUM((jpmnh(:,:,n1:n2) - pwr_spec(:,:,n1:n2,4))**2)
#if defined(MPI_OPT)
         IF (PARSOLVER) THEN
            CALL MPI_ALLREDUCE(MPI_IN_PLACE, d1, 2, MPI_REAL8, MPI_SUM,        &
                               SIESTA_COMM, MPI_ERR)
		 END IF
#endif		 
         diag_b = d1(1)
         diag_p = d1(2)
         IF (iam .EQ. 0) THEN
            DO js = 6, unit_out, unit_out - 6
               IF (.NOT.lverbose .AND. js.EQ.6) CYCLE
               WRITE(js, 1002) def(itype), diag_b, diag_p
            END DO
         END IF
      END IF

1000  FORMAT('ipar: ',a,/,'  n       B^s        r12*B^u')
1001  FORMAT(i4,1p,2e14.4)
1002  FORMAT(1x,'POWER SPECTRA(',a,') -- dB: ',1p,e10.3,' dP: ',1p,e10.3)
#endif
      END SUBROUTINE

!-------------------------------------------------------------------------------
!>  @brief Reset the perturbations of the SIESTA state.
!-------------------------------------------------------------------------------
      SUBROUTINE Clear_Field_Perts

      IMPLICIT NONE

!  Start of executable code.
      djbsupsmnsh = 0
      djbsupumnch = 0
      djbsupvmnch = 0
      djpmnch     = 0

      IF (lasym) THEN
         djbsupsmnch = 0
         djbsupumnsh = 0
         djbsupvmnsh = 0
         djpmnsh     = 0
      END IF

      END SUBROUTINE

      END MODULE
