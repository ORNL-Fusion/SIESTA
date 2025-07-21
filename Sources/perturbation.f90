!*******************************************************************************
!>  @file utilities.f90
!>  @brief Contains module @ref perturbation
!
!  Note separating the Doxygen comment block here so detailed decription is
!  found in the Module not the file.
!
!>  Subroutines and functionsto apply the helical perturbation to open up
!>  islands.
!*******************************************************************************

      MODULE perturbation
      USE v3_utilities, ONLY: assert, assert_eq
      USE quantities
      USE descriptor_mod, ONLY: iam
      USE siesta_state, ONLY: update_state, update_state
      USE shared_data, ONLY: ngmres_type, iortho, lcolscale, lasym,            &
                             hesspass_test, buv_res, unit_out
      USE nscalingtools, ONLY: startglobrow, endglobrow

      IMPLICIT NONE
      
      PRIVATE

!*******************************************************************************
!  MODULE VARIABLES
!*******************************************************************************
!>  Minimum radial index.
      INTEGER                                 :: nsmin
!>  Maximum radial index.
      INTEGER                                 :: nsmax
!>  Width of radial indicies the helial perturbation is applied to.
      INTEGER                                 :: irad_scale = 11
!>  Maximum number of iterations.
      INTEGER, PUBLIC                         :: niter_max
#if defined(JDOTB_PERT)
!>  Stellarator symmetric resonant components of jdotb.
      REAL(dp), ALLOCATABLE, DIMENSION(:,:,:) :: bmnc_res
!>  Stellarator asymmetric resonant components of jdotb.
      REAL(dp), ALLOCATABLE, DIMENSION(:,:,:) :: bmns_res
#endif

      PUBLIC :: init_data, add_perturb

      CONTAINS

!*******************************************************************************
!  UTILITY SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  @brief Load the namelist
!>
!>  READS siesta.jcf (JOB CONTROL FILE) for data (generalized to
!>  siesta_<name>.jcf, where <name> is the command line argument (ignored if not
!>  present!)
!>
!>  @param[in] namelist_file Path to the siesta namelist input file.
!>  @returns The preconditioner control value.
!>  FIXME: This doesn't really have much todo with the applying the helical
!>         perturbation. This should be somewhere else.
!-------------------------------------------------------------------------------
      FUNCTION init_data(namelist_file)
      USE date_and_computer
      USE siesta_namelist

      IMPLICIT NONE

!  Declare Arguments
      INTEGER                       :: init_data
      CHARACTER (len=*), INTENT(in) :: namelist_file

!  Local variables
      INTEGER                       :: istat
      INTEGER                       :: index1
      INTEGER                       :: imon
      CHARACTER(LEN=10)             :: date0
      CHARACTER(LEN=10)             :: time0
      CHARACTER(LEN=10)             :: zone0
      CHARACTER(LEN=256)            :: temp
      CHARACTER(LEN=256)            :: filename

!  Start of executable code
      CALL siesta_namelist_read(namelist_file)
      init_data = nprecon

      irad_scale = MIN(irad_scale, nsin)
      ngmres_type = MIN(2, MAX(1, ngmres_type))
      iortho = MIN(3, MAX(0, iortho))
      niter_max = niter

!  Separate out the name form the wout file name. Wout files sould have the form
!  of wout_name.nc. Remove the leading wout_ and trailing .nc.
      IF ((wout_file(1:5) .ne. 'WOUT_' .and.                                   &
           wout_file(1:5) .ne. 'wout_') .or.                                   &
          wout_file(LEN_TRIM(wout_file) - 2:LEN_TRIM(wout_file)) .ne.          &
          '.nc') THEN
         CALL assert(.false., 'ERROR: ' // TRIM(wout_file) //                  &
     &                        ' is not a valid wout file name.')
      ELSE
         temp = wout_file(6:LEN_TRIM(wout_file) - 3)
      END IF

      istat = 0
      IF (iam .eq. 0) THEN
	     WRITE (filename,1000) TRIM(temp), nsin, mpolin, ntorin
	     OPEN (unit=unit_out, file=filename, iostat=istat)
      END IF
      CALL assert_eq(istat, 0, 'ERROR WRITING SIESTA OUTPUT FILE')
      IF (iam .ne. 0) THEN
         RETURN
      END IF

      IF (lverbose) THEN
         WRITE (*,1001) '4.0', 100917
      END IF
      WRITE (unit_out,1001) '4.0', 100917

      WRITE (unit_out, 1002) TRIM(temp)
      IF (lverbose) THEN
         WRITE (*,1003) TRIM(temp)
      END IF

!  Format date and time.
      CALL DATE_AND_TIME(date0, time0, zone0)
      READ (date0(5:6),2000) imon
      WRITE (unit_out,1004) months(imon), date0(7:8), date0(1:4),              &
                            time0(1:2), time0(3:4), time0(5:6)

      DO istat = 1, SIZE(HelPert)
         IF (HelPert(istat)  .ne. zero) THEN
            WRITE (unit_out, 1005) istat, mres(istat), HelPert(istat),         &
                                   HelPhase(istat)
         END IF
      END DO

      WRITE (unit_out, 1006) ngmres_type, iortho, lColScale
      CALL FLUSH(unit_out)

1000  FORMAT('output_',a,'_',i4.4,2('X',i3.3),'.txt')
1001  FORMAT(62('-'),/,                                                        &
             1x,'SIESTA MHD EQUILIBRIUM CODE v',a,' (',i6,')', /,1x,           &
             'Scalable Island Equilibrium Solver for Toroidal Applications',/, &
             62('-'),/)
1002  FORMAT('CASE: ',a)
1003  FORMAT('OUTPUT FILE (SCREEN DUMP): output_',a,'.txt')
1004  FORMAT(' DATE = ',a3,' ',a2,',',a4,' ',' TIME = ',2(a2,':'),a2,/)
1005  FORMAT(i2,' mres: ',i4,' HelPert: ',1pe9.2,' HelPhase: ',1pe9.2)
1006  FORMAT(/,' ngmres_type: ', i4,' iOrtho: ', i4, ' lColScale: ', l2)

2000  FORMAT(i2)

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Add helical perturbation at the resonance.
!>
!>  @brief[out] xc      Displacement vector.
!>  @brief[in]  getwmhd Callback function to compute the total stored energy.
!-------------------------------------------------------------------------------
      SUBROUTINE add_perturb(xc, getwmhd)
      USE siesta_error
      USE fourier, ONLY: f_cos, f_sin, f_none
      USE siesta_namelist, ONLY: HelPert, HelPhase, eta_factor, lresistive, mres
      USE island_params, ONLY: fourier_context
      USE stel_constants

      IMPLICIT NONE

!  Declare Arguments
      REAL (dp), DIMENSION(:), INTENT(out)    :: xc
! FIXME: Not sure why were using a call back function here.
!        There is only one getwmhd function and it never changes.
      INTERFACE
         FUNCTION getwmhd(p)
         USE stel_kinds
         REAL (dp)                           :: getwmhd
         REAL (dp), DIMENSION(:), INTENT(in) :: p
         END FUNCTION
      END INTERFACE

!  Local Variables
      INTEGER                                 :: js
      INTEGER                                 :: l
      INTEGER                                 :: iprint
      INTEGER                                 :: mres0
      INTEGER                                 :: nres0
      INTEGER                                 :: n
      INTEGER                                 :: icount
      INTEGER                                 :: irscale
      INTEGER                                 :: imin
      INTEGER                                 :: jstart
      INTEGER                                 :: istat
      INTEGER                                 :: jsave
      REAL (dp)                               :: w0
      REAL (dp)                               :: w1
      REAL (dp)                               :: wmhd
      REAL (dp)                               :: eta_factor_save
      REAL (dp)                               :: p_width_min
      REAL (dp)                               :: normal
      REAL (dp)                               :: HelPert0
      REAL (dp)                               :: HelPert0A
      REAL (dp)                               :: rad
      REAL (dp)                               :: chip0
      REAL (dp)                               :: phip0
      LOGICAL                                 :: lresist_save
      REAL (dp)                               :: p_width
#if defined(JDOTB_PERT)
      REAL(dp), ALLOCATABLE, DIMENSION(:,:,:) :: jdotb
#endif

!  Start of executable code.

!  Add helical RESISTIVE flux perturbation on full mesh. Approximate nonideal
!  E_sub ~ f(mu+nv) B_sub in update_bfield (E ~ B).
      IF (ntor .eq. 0 .or. (ALL(HelPert  .eq. zero))) THEN
         RETURN
      END IF

!  Initial state. Allocate buv_res so bsubXijf will be computed in init_state.
      nsmin = MAX(1, startglobrow - 1)
      nsmax = MIN(ns, endglobrow + 1)

      ALLOCATE(buv_res(ntheta,nzeta,nsmin:nsmax), stat=istat)
      CALL ASSERT(istat .eq. 0,'ALLOCATION ERROR IN add_perturbation')

#if defined(JDOTB_PERT)
      ALLOCATE(jdotb(ntheta,nzeta,nsmin:nsmax), stat=istat)
      CALL ASSERT(istat .eq. 0,'ALLOCATION ERROR IN add_perturbation')
#endif

      lresist_save = lresistive
      lresistive = .false.
      eta_factor_save = eta_factor
      eta_factor = 1
      xc = 0

!  Get latest current and bsubXijf.
      w0 = getwmhd(xc)

      IF (iam .eq. 0) THEN
         DO iprint = 6, unit_out, unit_out - 6
            IF (lverbose .or. iprint .ne. 6) THEN
               WRITE (iprint, 1000)
            END IF
         END DO
      ENDIF

#if defined(JDOTB_PERT)
!  Not used at the moment.
!  Compute jdotb to extract resonant components of sqrt(g)*J dot B.
      jdotb = ksupsijf0(:,:,nsmin:nsmax)*bsubsijf(:,:,nsmin:nsmax)             &
            + ksupuijf0(:,:,nsmin:nsmax)*bsubuijf(:,:,nsmin:nsmax)             &
            + ksupvijf0(:,:,nsmin:nsmax)*bsubvijf(:,:,nsmin:nsmax)

      ALLOCATE (bmnc_res(0:mpol,-ntor:ntor, nsmin:nsmax),                      &
                bmns_res(0:mpol,-ntor:ntor, nsmin:nsmax), stat=istat)
      CALL ASSERT(istat .eq. 0, 'ISTAT != 0 IN ADD_PERT')
      
      bmns_res = 0
      CALL fourier_context%tomnsp(jdotb, bmnc_res, f_cos)
      IF (lasym) THEN
         CALL fourier_context%tomnsp(jdotb, bmns_res, f_sin)
      END IF
#endif

!  Compute responses to resonant resistitive perturbations
      lresistive = .true.
      normal = hs_i*hs_i

      DO icount = 1, SIZE(HelPert)
         mres0 = ABS(mres(icount))
         IF (HelPert(icount) .eq. zero .or. mres0 .gt. mpol) THEN
            CYCLE
         END IF
         IF (lasym) THEN
            HelPert0 =  ABS(HelPert(icount)/normal)*COS(HelPhase(icount)*pi/180.0)
            HelPert0A= -ABS(HelPert(icount)/normal)*SIN(HelPhase(icount)*pi/180.0)
         ELSE
            HelPert0 = ABS(HelPert(icount)/normal)
            HelPert0A = 0
            IF (HelPhase(icount) .eq. 180.0) THEN
               HelPert0 = -HelPert0
            END IF
         END IF

!  Scan in radius to determine primary resonance.
         DO n = -ntor, ntor
!  Avoid 0/0 chip/phip at origin
            jstart = 3

            nres0 = fourier_context%tor_modes(n)

!  Look for multiple resonances
            DO WHILE (jstart .lt. ns - 1)

               IF (.not.FindResonance(mres0, nres0, rad, jstart)) THEN
                  EXIT
               END IF

               w1 = 2*w0
               imin = 1

!  NOTE: in update_bfield, where perturbed b-field is computed, this is
!  multiplied by eta_prof = rho*(1 - rho) so it vanishes at both ends. Multiply
!  again here by that factor to assure radial derivatives entering in dB also
!  vanish at endpoints s=0,1
 
               DO irscale = 1, irad_scale
                  CALL GetResPert(irscale, mres0, n, rad, HelPert0, HelPert0A, &
                                  chip0, phip0, p_width)
                  xc = 0
!  FIXME: Why call this again with the same inputs?
                  wmhd = getwmhd(xc)
                  IF (wmhd .lt. w1) THEN
                     imin = irscale
                     w1 = wmhd
                     p_width_min = p_width
                  END IF
               END DO

!  Make sure energy decreases.
               wmhd = w1
!  Recompute perturbation buv_res here.
               irscale = imin
         
               IF (iam .eq. 0) THEN
                  DO iprint = 6, unit_out, unit_out - 6
                     IF (lverbose .or. iprint .ne. 6) THEN
                        WRITE (iprint, 1001) 10**6*(wmhd - w0), mres0, nres0,  &
                                             HelPert0*normal, rad,             &
                                             ABS(mres0*chip0 +                 &
                                                 nres0*nfp*phip0),             &
                                             chip0/phip0, p_width_min
                     END IF
                  END DO

                  fourier_context%found_modes(mres0, n) = .true.
                  CALL FLUSH(unit_out)
               END IF
         
               CALL GetResPert(irscale, mres0, n, rad, HelPert0, HelPert0A,    &
                               chip0, phip0, p_width)

               wmhd = getwmhd(xc)
               IF (ABS(wmhd - w1) .gt. 1.E-12_dp) THEN
                  CALL siesta_error_set_error(siesta_error_general,            &
                                              'Error1 in Perturb')
               END IF

               xc = 0
               CALL update_state(.false., zero, zero)
!  Compute dW relative to perturbed state.
               w0 = wmhd
   
            END DO ! Look for multiple resonances
         END DO
      END DO

      IF (iam .eq. 0) THEN
         DO iprint = 6, unit_out, unit_out - 6
            IF (lverbose .or. iprint .ne. 6) THEN
               WRITE (iprint, *)
            END IF
         END DO
      END IF

#if defined(JDOTB_PERT)
      DEALLOCATE(jdotb, bmnc_res, bmns_res)
#endif

      DEALLOCATE(buv_res)
      lresistive = lresist_save
      eta_factor = eta_factor_save

1000  FORMAT(' Adding helical magnetic field perturbations'/' 10^6 X Del-W',   &
             '    mres    nres     HelPert     rad  |m*chip+n*phip|',          &
             '     iota   radial width')
1001  FORMAT(1p,e12.3,2(i8),3x,1pe10.2,0p,f8.2,f11.2,4x,2(f11.2))

      END SUBROUTINE

!-------------------------------------------------------------------------------
!>  @brief Find radial index of the resonance.
!>
!>  Starting a given index this searched for the point where
!>
!>     m*chi' + n*nfp*phi' = 0                                               (1)
!>
!>  @param[in]    mres0  Poloidal mode number resonance.
!>  @param[in]    nres0  Toroidal mode number resonance.
!>  @param[out]   resrad Radial position of the resonance.
!>  @param[inout] jstart Starting serach index and final radial index of the
!>                       Resonance.
!>  @returns Ture if the resonace is located.
!-------------------------------------------------------------------------------
      FUNCTION FindResonance(mres0, nres0, resrad, jstart)

      IMPLICIT NONE

!  Declare Arguments
      LOGICAL                :: FindResonance
      INTEGER, INTENT(in)    :: mres0
      INTEGER, INTENT(in)    :: nres0
      INTEGER, INTENT(inout) :: jstart
      REAL(dp), INTENT(out)  :: resrad

!  Local Variables
      INTEGER                :: js
      INTEGER                :: jmin
      REAL (dp)              :: del1
      REAL (dp)              :: del2
      REAL (dp)              :: rho1
      REAL (dp)              :: rho2
      REAL (dp)              :: delmin

!  Start of executable code.
      del2 = 0
      resrad = 0
      delmin = -1
      IF (mres0*nres0 .eq. 0) THEN
          FindResonance = .false.
          RETURN
      END IF

!  Check for zero crossing
      jmin = jstart
      DO js = jstart, ns - 1
         del1 = mres0*chipf_i(js) + nres0*nfp*phipf_i(js)
         del1 = del1/MAX(ABS(mres0*chipf_i(js)),                               &
                         ABS(nres0*nfp*phipf_i(js)),                           &
                         1.E-20_dp)
         IF (delmin .eq. -1 .or. ABS(del1) .lt. delmin) THEN
            IF (delmin .eq. -1) THEN
               del2 = del1
            END IF
            delmin = ABS(del1)
            jmin = js
         END IF
         
         IF (del1*del2 .lt. zero) THEN
            jstart = js + 1
            rho2 = hs_i*(js - 2)
            rho1 = hs_i*(js - 1)
            resrad = (rho2*ABS(del1) + rho1*ABS(del2))/(ABS(del1) + ABS(del2))
            delmin = 0
            EXIT
         END IF

         del2 = del1
      END DO

! If no zero-crossing, resonance might be very close.
      IF (delmin .lt. 1.E-3_dp) THEN
         jstart = jmin + 1
         resrad = hs_i*(jmin - 1)
         FindResonance = .true.
      ELSE
         jstart = ns - 1
         resrad = 0
         FindResonance = .false.
      END IF

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Get the resonant perturbation.
!>
!>  @param[in]  irscale
!>  @param[in]  mres0     Poloidal mode number for resonace.
!>  @param[in]  nres0     Toroidal mode number for resonace.
!>  @param[in]  HelPert0  Helical pertubation amplitude.
!>  @param[in]  HelPert0A Asymmetric Helical pertubation amplitude.
!>  @param[out] chip0     Poloidal flux derivative at resonance.
!>  @param[out] phip0     Toroidal flux derivative at resonance.
!>  @param[out] p_width   With of perturbation.
!-------------------------------------------------------------------------------
      SUBROUTINE GetResPert(irscale, mres0, nres0, rad, HelPert0, HelPert0A,   &
     &                      chip0, phip0, p_width)
      USE fourier, ONLY: f_cos, f_sin, f_sum, f_none
      USE siesta_namelist, ONLY: mpolin, ntorin
      USE island_params, ONLY: fourier_context

      IMPLICIT NONE

!  Declare Arguments
      INTEGER, INTENT(in)                      :: irscale
      INTEGER, INTENT(in)                      :: mres0
      INTEGER, INTENT(in)                      :: nres0
      REAL(dp), INTENT(in)                     :: rad
      REAL(dp), INTENT(in)                     :: HelPert0
      REAL(dp), INTENT(in)                     :: HelPert0A
      REAL(dp), INTENT(out)                    :: chip0
      REAL(dp), INTENT(out)                    :: phip0
      REAL(dp), INTENT(out)                    :: p_width

!  Local Variables
      INTEGER                                  :: js
      INTEGER                                  :: istat
      INTEGER                                  :: n
      REAL (dp)                                :: rho
      REAL (dp)                                :: rhores
      REAL (dp), DIMENSION(:), ALLOCATABLE     :: pert_prof
      REAL (dp)                                :: locrad
      REAL (dp), DIMENSION(:,:,:), ALLOCATABLE :: bmn_res
#if defined(JDOTB_PERT)
      REAL (dp), DIMENSION(:,:,:), ALLOCATABLE :: temp_res
#endif

!  Start of executable code.
      IF (mres0 .gt. mpolin .or. ABS(nres0) .gt. ntorin) THEN
         RETURN
      END IF

!  Compute radial form factor (window-function)
      js = INT(rad/hs_i) + 1
      IF (js .lt. 1 .or. js .gt. ns) THEN
         RETURN
      END IF
      chip0 = chipf_i(js)
      phip0 = phipf_i(js)
      locrad = REAL(irscale - 1, dp)/(irad_scale - 1)
!  Decay length, ~delta-fcn for irscale=1 to broad (full radius)
      p_width = 0.1_dp + 0.9_dp*locrad
      locrad = one/p_width
      ALLOCATE(pert_prof(nsmin:nsmax))
      DO js = nsmin, nsmax
         rho = hs_i*(js - 1)
         rhores = rho - rad
         rhores = locrad*rhores
!  Tearing parity, multiply by alpha for equal areas.
         pert_prof(js) = one/(one + rhores*rhores)
      END DO

#if defined(JDOTB_PERT)
!  This doesn't seem to work as well as the simple form factor do keep
!  jdotb_pert undefined until we can improve it.
      ALLOCATE(temp_res(0:mpol,-ntor:ntor,nsmin:nsmax))

      rho = MAX(MAXVAL(ABS(bmnc_res(mres0,nres0,:))),                   &
                MAXVAL(ABS(bmns_res(mres0,nres0,:))))

      temp_res = 0
      temp_res(mres0,nres0,:) = HelPert0*bmnc_res(mres0,nres0,:)
      CALL fourier_context%toijsp(temp_res, buv_res, f_none, f_cos)
      
      IF (lasym) THEN
         temp_res(mres0,nres0,:) = HelPert0A*bmns_res(mres0,nres0,:)
         CALL fourier_context%toijsp(temp_res, buv_res, f_sum, f_sin)
      END IF
      
      IF (rho .ne. zero) THEN
         buv_res = buv_res/rho
      END IF

!	Apply radial form-factor
      DO js = nsmin, nsmax
         buv_res(:,:,js) = buv_res(:,:,js)*pert_prof(js)
      END DO

      DEALLOCATE(temp_res, stat=istat)
#else
      ALLOCATE(bmn_res(0:mpol,-ntor:ntor,nsmin:nsmax), stat=istat)
      CALL ASSERT(istat .eq. 0, 'ISTAT != 0 IN GETRESPERT')
      bmn_res = 0
      buv_res = 0

      IF (HelPert0 .ne. 0.0) THEN
         bmn_res(mres0,nres0,:) = HelPert0*pert_prof(:)
         CALL fourier_context%toijsp(bmn_res, buv_res, f_none, f_cos)
      END IF

      IF (lasym .and. HelPert0A .ne. 0.0) THEN
         bmn_res(mres0,nres0,:) = HelPert0A*pert_prof(:)
         CALL fourier_context%toijsp(bmn_res, buv_res, f_sum, f_sin)
      END IF
      
      DEALLOCATE(bmn_res, stat=istat)
#endif
      DEALLOCATE(pert_prof)

      END SUBROUTINE

      END MODULE
