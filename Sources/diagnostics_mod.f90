      MODULE diagnostics_mod
      USE v3_utilities, ONLY: assert
      USE fourier, ONLY: f_sum, f_cos, f_sin, f_none, f_du, f_dv
      USE descriptor_mod, ONLY: iam, nprocs, SIESTA_COMM
      USE nscalingtools, ONLY: startglobrow, endglobrow, MPI_ERR
      USE timer_mod
      USE quantities
      USE metrics, ONLY: tolowerh, nsh
      USE utilities, ONLY: to_half_mesh
#if defined(MPI_OPT)
      USE mpi_inc
#endif
      USE island_params, ONLY: fourier_context

      IMPLICIT NONE

      INTEGER, PRIVATE  :: unit35=35, unit36=36, nsmin, nsmax
      REAL(dp) :: divb_rms, avbeta, bdotj_rms, divj_rms, bdotj2_rms
      REAL(dp) :: toroidal_flux, toroidal_flux0=0, bgradp_rms, wbgradp
      REAL(dp) :: max_bgradp, min_bgradp
      REAL(dp) :: tnorm_bgradp
!
!     The next three allocatable variables are kept in memory and used to
!     produce output whenever L_DUMP_DIAGNO = .TRUE. 
!      
      REAL(dp), ALLOCATABLE, DIMENSION(:,:,:) :: bdotjmnch, bdotjijh
      REAL(dp), ALLOCATABLE, DIMENSION(:,:,:) :: bdotjmnsh
      REAL(dp), ALLOCATABLE, DIMENSION(:,:,:) :: divjmnsh, divbmnsf
      REAL(dp), ALLOCATABLE, DIMENSION(:,:,:) :: divjmnch, divbmncf
      REAL(dp), ALLOCATABLE, DIMENSION(:,:,:) :: bgradpf, jpmnch_loc
      REAL(dp), ALLOCATABLE, DIMENSION(:,:,:), TARGET :: fgradpmnch
      REAL(dp), DIMENSION(:), POINTER     :: xc, gc
      LOGICAL :: lextend = .FALSE.                   !Set in call to EVOLVE_BGRADP
      LOGICAL :: lcurr_init = .FALSE.
!
!     CALL THIS AFTER THE B-FIELD COMPONENTS bi = Jac*B-sup-i HAVE BEEN UPDATED
!     IN UPDATE_STATE (CALL FROM UPDATE_STATE IS SAFEST WAY!) 

      CONTAINS

!*******************************************************************************
!>  UTILITY SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  @brief Compute divergence of B on the full mesh points.
!>
!>  @param[in] ns_min Lower radial index.
!>  @param[in] ns_max Upper radial index.
!-------------------------------------------------------------------------------
      SUBROUTINE divb(ns_min, ns_max)
      USE island_params, ONLY: ohs=>ohs_i
      USE hessian, ONLY: gather_array
      USE utilities, ONLY: div_h

!  Declare Arguments
      INTEGER, INTENT(in)     :: ns_min
      INTEGER, INTENT(in)     :: ns_max

!  Local Variables
      INTEGER                 :: istat
      INTEGER                 :: nl
      INTEGER                 :: nh
INTEGER                 :: s,m,n
      REAL (dp)               :: tnorm
      REAL (dp), DIMENSION(2) :: temp
      REAL (dp)               :: ton
      REAL (dp)               :: toff

!  Start of executable code
      CALL second0(ton)

      nsmin = ns_min
      nsmax = ns_max
      nl = MAX(1, nsmin)
      nh = MIN(ns - 1, nsmax)

      IF (ALLOCATED(divbmnsf)) THEN
         DEALLOCATE(divbmnsf)
      END IF
      IF (ALLOCATED(divbmncf)) THEN
         DEALLOCATE(divbmncf)
      END IF
      ALLOCATE(divbmnsf(0:mpol,-ntor:ntor,nsmin:nsmax),                      &
               divbmncf(0:mpol,-ntor:ntor,nsmin:nsmax), stat=istat)
      CALL ASSERT(istat.EQ.0,'Allocation failed in divb')

      divbmnsf = 0
      divbmncf = 0

      CALL div_h(jbsupsmnsh, jbsupumnch, jbsupvmnch, divbmnsf,                 &
                 f_sin, nsmin, nsmax)
      IF (lasym) THEN
         CALL div_h(jbsupsmnch, jbsupumnsh, jbsupvmnsh, divbmncf,              &
                    f_cos, nsmin, nsmax)
      END IF

      tnorm = hs_i*SUM(jbsupumnch(:,:,nsmin:nsmax)**2                          &
            +          jbsupvmnch(:,:,nsmin:nsmax)**2                          &
            +          jbsupsmnsh(:,:,nsmin:nsmax)**2)

      divb_rms = SUM(divbmnsf(:,:,nl:nh)**2 + divbmncf(:,:,nl:nh)**2)
#if defined(MPI_OPT)
      IF (PARSOLVER) THEN
         temp(1) = divb_rms
         temp(2) = tnorm
         CALL MPI_ALLREDUCE(MPI_IN_PLACE,temp,2,MPI_REAL8, MPI_SUM,            &
                            SIESTA_COMM,MPI_ERR)
         divb_rms = temp(1)
         tnorm = temp(2)
      END IF
#endif
      IF (tnorm .NE. zero) THEN
         divb_rms = SQRT(divb_rms/tnorm)
         divbmnsf = divbmnsf/SQRT(tnorm)
         divbmncf = divbmncf/SQRT(tnorm)
      END IF

      CALL second0(toff)
      time_divb = time_divb + (toff-ton)

      END SUBROUTINE

!-------------------------------------------------------------------------------
!>  @brief Compute divergence of dB on the full mesh points.
!>
!>  @param[in] ns_min Lower radial index.
!>  @param[in] ns_max Upper radial index.
!-------------------------------------------------------------------------------
      FUNCTION divdb(ns_min, ns_max)
      USE island_params, ONLY: ohs=>ohs_i, fourier_context
      USE hessian, ONLY: gather_array

!  Declare Arguments
      REAL (rprec)            :: divdb
      INTEGER, INTENT(in)     :: ns_min
      INTEGER, INTENT(in)     :: ns_max

!  Local Variables
      INTEGER                 :: istat
      INTEGER                 :: m
      INTEGER                 :: n
      INTEGER                 :: n_mode
      INTEGER                 :: nl
      INTEGER                 :: nh
      REAL (dp)               :: tnorm
      REAL (dp), DIMENSION(2) :: temp
      REAL (dp)               :: ton
      REAL (dp)               :: toff

!  Start of executable code
      CALL second0(ton)

      nsmin = ns_min
      nsmax = ns_max
      nl = MAX(2, nsmin)
      nh = MIN(ns - 1, nsmax)

      IF (ALLOCATED(divbmnsf)) THEN
         DEALLOCATE(divbmnsf)
      END IF
      IF (ALLOCATED(divbmncf)) THEN
         DEALLOCATE(divbmncf)
      END IF
      ALLOCATE(divbmnsf(0:mpol,-ntor:ntor,nsmin:nsmax),                        &
               divbmncf(0:mpol,-ntor:ntor,nsmin:nsmax), stat=istat)
      CALL ASSERT(istat.EQ.0,'Allocation failed in divb')

      divbmnsf = 0
      divbmncf = 0

      DO n = -ntor, ntor
         n_mode = fourier_context%tor_modes(n)*nfp
         DO m = 0, mpol
            divbmnsf(m,n,nl:nh) =     -m*(djbsupumnch(m,n,nl+1:nh+1) +         &
     &                                    djbsupumnch(m,n,nl:nh))*0.5          &
     &                          - n_mode*(djbsupvmnch(m,n,nl+1:nh+1) +         &
     &                                    djbsupvmnch(m,n,nl:nh))*0.5          &
     &                          +   ohs*(djbsupsmnsh(m,n,nl+1:nh+1) -          &
     &                                   djbsupsmnsh(m,n,nl:nh))
            IF (lasym) THEN
               divbmncf(m,n,nl:nh) =      m*(djbsupumnsh(m,n,nl+1:nh+1) +      &
     &                                       djbsupumnsh(m,n,nl:nh))*0.5       &
     &                             + n_mode*(djbsupvmnsh(m,n,nl+1:nh+1) +      &
     &                                       djbsupvmnsh(m,n,nl:nh))*0.5       &
     &                             +    ohs*(djbsupsmnch(m,n,nl+1:nh+1) -      &
     &                                       djbsupsmnch(m,n,nl:nh))
            END IF
         END DO
      END DO

      tnorm = hs_i*SUM(djbsupumnch(:,:,nsmin:nsmax)**2                         &
            +          djbsupvmnch(:,:,nsmin:nsmax)**2                         &
            +          djbsupsmnsh(:,:,nsmin:nsmax)**2)

      divdb = SUM(divbmnsf**2 + divbmncf**2)
#if defined(MPI_OPT)
      temp(1) = divdb
      temp(2) = tnorm
      CALL MPI_ALLREDUCE(MPI_IN_PLACE,temp,2,MPI_REAL8, MPI_SUM,               &
                         SIESTA_COMM,MPI_ERR)
      divdb = temp(1)
      tnorm = temp(2)
#endif

      IF (tnorm .NE. zero) THEN
         divdb = SQRT(divdb/tnorm)
         divbmnsf = divbmnsf/SQRT(tnorm)
         divbmncf = divbmncf/SQRT(tnorm)
      END IF

      CALL second0(toff)
      time_divb = time_divb + (toff-ton)

      END FUNCTION

        SUBROUTINE WRITE_PROFILES(fsq_total)
        USE safe_open_mod
        USE siesta_init, ONLY: init_state
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
        REAL(dp) :: fsq_total
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
        INTEGER, PARAMETER       :: ifull=0, ihalf=1
        INTEGER                  :: istat, js
        REAL(dp), ALLOCATABLE, DIMENSION(:,:,:) :: pijh, pmnh
        REAL(dp)                 :: average(ns)
!-----------------------------------------------
!
!       IT IS ASSUMED COMING IN HERE THAT THE PROFILES HAVE BEEN PREVIOUSLY
!       GATHERED (IN UPDATE_STATE) AND PREPARE_QUANTITIES HAS BEEN CALLED
!
        IF (iam .NE. 0) RETURN
        
        ALLOCATE (pijh(ntheta,nzeta,ns),                                &
                  pmnh(0:mpol,-ntor:ntor,ns), stat=istat)      
        CALL ASSERT(istat.EQ.0,'Allocation failed in WRITE_PROFILES')

        CALL safe_open(unit35, istat, "siesta_profiles.txt", 'replace', &
                       'formatted')
        CALL safe_open(unit36, istat, "siesta_profiles_pest.txt",       &
                       'replace','formatted')

        CALL ASSERT(istat.EQ.0,'Error opening siesta profiles data file')

!Check that init_state was initialized so start/end globrow = (1,ns)
        CALL ASSERT(startglobrow.EQ.1 .AND. endglobrow.EQ.ns,           &
          ' INIT_STATE NOT INITIALIZED IN WRITE_PROFILES!')
          
!USE pijh, pmnh scratch arrays to compute pressure and bsupX components on half-mesh

!pressure
        CALL WRITE_SPECTRA(pijf0, pmnh, f_cos, f_sin, ifull, 'p', fsq_total)

!BSUPS
        CALL WRITE_SPECTRA(bsupsijf0, pmnh, f_sin, f_cos, ifull, 'B^s', fsq_total)

!BSUPU
        CALL WRITE_SPECTRA(bsupuijf0, pmnh, f_cos, f_sin, ifull, 'B^u', fsq_total)

!BSUPV
        CALL WRITE_SPECTRA(bsupvijf0, pmnh, f_cos, f_sin, ifull, 'B^v', fsq_total)

!KSUPS CRCook
        pijh = ksupsijf0/jacobf
        CALL WRITE_SPECTRA(pijh, pmnh, f_sin, f_cos, ifull, 'J^s', fsq_total)

!KSUPU CRCook
        pijh = ksupuijf0/jacobf
        CALL WRITE_SPECTRA(pijh, pmnh, f_cos, f_sin, ifull, 'J^u', fsq_total)

!KSUPV CRCook
        pijh = ksupvijf0/jacobf
        CALL WRITE_SPECTRA(pijh, pmnh, f_cos, f_sin, ifull, 'J^v', fsq_total)

!K|| = (J*B)/B**2 (SPH:040915)
        CALL bdotj(bdotjijh)
        pijh = bdotjijh*jacobh
        DEALLOCATE(bdotjijh)
        CALL WRITE_SPECTRA(pijh, pmnh, f_cos, f_sin, ihalf, 'K||', fsq_total)


!<K||-nonaxi> : RWILCOX metric in PEST coordinates
!need to extend this to lasym=T
        CALL fourier_context%tomnsp(pijh, pmnh, f_cos)
        pmnh(:,0,:) = 0
        CALL fourier_context%toijsp(pmnh, pijh, f_none, f_cos)
        CALL SurfAverage(average, ABS(pijh), 1, ns)
        WRITE (unit36, *) 
        WRITE (unit36, *) '<ABS(K||) - nonaxisymmetric>'
        DO js = 2,ns
           WRITE (unit36, '(i4,1p,e10.2)') js, average(js)
        END DO
        
!Displacements X jac (SPH:052114)
        CALL siesta_profiles(unit35, jvsupsmncf, 'vel^s_cosmn;  fsq_total: ', fsq_total)
        IF (lasym) CALL siesta_profiles(unit35, jvsupsmnsf, 'vel^s_sinmn;  fsq_total: ', fsq_total)
        CALL siesta_profiles(unit35, jvsupumnsf, 'vel^u_sinmn;  fsq_total: ', fsq_total)
        IF (lasym) CALL siesta_profiles(unit35, jvsupumncf, 'vel^u_cosmn;  fsq_total: ', fsq_total)
        CALL siesta_profiles(unit35, jvsupvmnsf, 'vel^v_sinmn;  fsq_total: ', fsq_total)
        IF (lasym) CALL siesta_profiles (unit35, jvsupvmncf, 'vel^v_cosmn;  fsq_total: ', fsq_total)

        DEALLOCATE(pijh, pmnh)

        CLOSE (unit35)
        CLOSE (unit36)

      END SUBROUTINE WRITE_PROFILES
      
!-------------------------------------------------------------------------------
!>  @brief  Write fourier spectra for both regular and pest coordinates.
!>
!>  There is no lmn* quantity when using free boundary SIESTA so pest output is
!>  disabled. Arrays must be ALLOCATABLE to preserve array bounds.
!>
!>  @param[inout] xij       Real space quantity.
!>  @param[inout] xmn       Fourier space buffer.
!>  @param[in]    fsym      Symmetric parity.
!>  @param[in]    fasym     Asymmetric parity.
!>  @param[in]    igrid     Radial grid type. This should either be 0 for full
!>                          or 1 for half grid.
!>  @param[in]    label     Name of the quantity in xij.
!>  @param[in]    fsq_total
!>------------------------------------------------------------------------------
      SUBROUTINE write_spectra(xij, xmn, fsym, fasym, igrid, label, fsq_total)
      USE vmec_info, ONLY:  lmns => lmns_i, lmnc => lmnc_i
      USE siesta_namelist, ONLY: l_vessel
      USE utilities, ONLY: to_full_mesh

      IMPLICIT NONE

!  Declare Arguments
      REAL (dp), DIMENSION(:,:,:), ALLOCATABLE, INTENT(inout) :: xij
      REAL (dp), DIMENSION(:,:,:), ALLOCATABLE, INTENT(inout) :: xmn
      INTEGER, INTENT(IN)                                     :: fsym
      INTEGER, INTENT(IN)                                     :: fasym
      INTEGER, INTENT(IN)                                     :: igrid
      CHARACTER (len=*), INTENT(in)                           :: label
      REAL (dp), INTENT(IN)                                   :: fsq_total

!  Local variables
      REAL (dp), DIMENSION(:,:,:), ALLOCATABLE                :: xijf
      CHARACTER (LEN=128)                                     :: mylabel

!  Start of executable code.
      IF (igrid .EQ. 1) THEN
         ALLOCATE(xijf(SIZE(xij,1), SIZE(xij,2), SIZE(xij,3)))
         CALL to_full_mesh(xij, xijf)
         xij = xijf/jacobf
         DEALLOCATE(xijf)
      END IF

      IF (fsym .EQ. f_cos) THEN
         WRITE (mylabel,1000) TRIM(LABEL)
      ELSE
         WRITE (mylabel,1001) TRIM(LABEL)
      END IF
      
      CALL fourier_context%tomnsp(xij, xmn, fsym)
      CALL siesta_profiles(unit35, xmn, TRIM(mylabel), fsq_total)
      IF (.not.l_vessel) THEN
         CALL fourier_context%tomnsp(xij, xmn, lmns, lmnc, fsym, lasym)
         CALL siesta_profiles(unit36, xmn, TRIM(mylabel), fsq_total)
      END IF
      IF (lasym) THEN

         IF (fasym .EQ. f_cos) THEN
            WRITE (mylabel,1000) TRIM(LABEL)
         ELSE
            WRITE (mylabel,1001) TRIM(LABEL)
         END IF

         CALL fourier_context%tomnsp(xij, xmn, fasym)
         CALL siesta_profiles(unit35, xmn, TRIM(mylabel), fsq_total)
         IF (.not.l_vessel) THEN
            CALL fourier_context%tomnsp(xij, xmn, lmns, lmnc, fasym, lasym)
            CALL siesta_profiles(unit36, xmn, TRIM(mylabel), fsq_total)
         END IF
      END IF

1000  FORMAT(a,'_cosmn; fsq_total: ')
1001  FORMAT(a,'_sinmn; fsq_total: ')

      END SUBROUTINE

      SUBROUTINE SIESTA_PROFILES (iunit, arr_value, label, fsq_total)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER, INTENT(IN)       :: iunit
      REAL(dp), INTENT(IN)      :: arr_value(0:mpol,-ntor:ntor,ns)
      REAL(dp), INTENT(IN)      :: fsq_total
      CHARACTER*(*), INTENT(IN) :: label
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER  :: n, m, js
      INTEGER  :: mnmax
      REAL(dp) :: sj, rms_value
      CHARACTER (LEN=36) :: p_format
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: m_array
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: n_array
!-----------------------------------------------

      WRITE (iunit, *)
      WRITE (iunit, '(1x,a,1p,e12.4)') TRIM(label), fsq_total

      ALLOCATE(m_array(0:mpol,-ntor:ntor))
      ALLOCATE(n_array(0:mpol,-ntor:ntor))

      DO n = -ntor, ntor
         DO m = 0, mpol
            m_array(m,n) = m
         END DO
      END DO

      mnmax = (mpol + 1)*(2*ntor + 1)

      WRITE (p_format,1000) mnmax
      WRITE (iunit,p_format) '      ','   ','MPOL--->', m_array
      WRITE (p_format,1000) mnmax
      WRITE (iunit,p_format) 'RADIUS','RMS','NTOR--->',                        &
     &                       fourier_context%tor_modes

      WRITE (p_format,1002) mnmax
      DO js = 2, ns
         sj = hs_i*(js - 1)
         rms_value = SQRT(SUM(arr_value(:,:,js)**2)/mnmax)
         WRITE(iunit,p_format) sj, rms_value, arr_value(:,:,js)
      END DO

      DEALLOCATE(m_array, n_array)

1000  FORMAT('(a,8x,a,6x,a,',i3,'(2x,i12))')
1002  FORMAT('(f6.3,2x,es12.5,14x,',i4,'(2x,es12.5))')

      END SUBROUTINE SIESTA_PROFILES


      SUBROUTINE bgradp(ns_min, ns_max)
!     
!     WRITTEN 03-21-13 BY S. HIRSHMAN AS PART OF THE ORNL SIESTA PROJECT (c)
!     
!     PURPOSE: Computes B*GRAD P (normalized to B*P) on the FULL mesh
!
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER, INTENT(IN) :: ns_min, ns_max
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER     :: istat
      REAL(dp)    :: s1(2), r1(2), ton, toff
!-----------------------------------------------
!
!     COMPUTE B DOT GRAD P AT FULL-GRID POINTS (sans endpts)
!
      nsmin=ns_min; nsmax=ns_max
      
      ALLOCATE(bgradpf(ntheta,nzeta,nsmin:nsmax), stat=istat)
      CALL ASSERT(istat.EQ.0,'Allocation error in BGRADP')
      CALL get_bgradp(nsmin, nsmax)

      DEALLOCATE (bgradpf, stat=istat)

      END SUBROUTINE bgradp

!-------------------------------------------------------------------------------
!>  @brief Compute B.Grad(p) on the full mesh sans end points.
!>
!>  @param[in] ns_min Minimum radial index.
!>  @param[in] ns_max Maximum radial index.
!-------------------------------------------------------------------------------
      SUBROUTINE get_bgradp(ns_min, ns_max)

      IMPLICIT NONE

!  Declare Arguments
      INTEGER, INTENT(in)                      :: ns_min
      INTEGER, INTENT(in)                      :: ns_max

!  local variables
      INTEGER                                  :: istat
      INTEGER                                  :: n1
      INTEGER                                  :: n2
      REAL (dp), DIMENSION(:,:,:), ALLOCATABLE :: dpduijf
      REAL (dp), DIMENSION(:,:,:), ALLOCATABLE :: dpdvijf
      REAL (dp), DIMENSION(2)                  :: temp(2)
      REAL (dp)                                :: ton
      REAL (dp)                                :: toff

!  Start of executable code
      CALL second0(ton)

      CALL assert(ALLOCATED(jbsupsmnsh),'jbsupXmn unallocated in get_bgradp')

      nsmin = ns_min
      nsmax = MIN(ns_max + 1, ns)

      bgradpf = 0 

!  Compute pressure and angle derivatives
      ALLOCATE(dpduijf(ntheta,nzeta,nsmin:nsmax),                              &
               dpdvijf(ntheta,nzeta,nsmin:nsmax), stat=istat)
      CALL assert_eq(0, istat, 'Allocation failed in get_bgradp')

      CALL to_full_mesh(pijh0_du, dpduijf)
      CALL to_full_mesh(pijh0_dv, dpdvijf)

      nsmax = MIN(endglobrow, ns)
      n1 = MAX(2, nsmin)
      n2 = MIN(nsh, nsmax)
!  Get all values except first and
      bgradpf(:,:,n1:n2) = bsupsijf0(:,:,n1:n2)*pijf0_ds(:,:,n1:n2)            &
                         + bsupuijf0(:,:,n1:n2)*dpduijf(:,:,n1:n2)             &
                         + bsupvijf0(:,:,n1:n2)*dpdvijf(:,:,n1:n2)

      tnorm_bgradp = SUM((bsupsijf0(:,:,n1:nsmax)**2 +                         &
                          bsupuijf0(:,:,n1:nsmax)**2 +                         &
                          bsupvijf0(:,:,n1:nsmax)**2)*                         &
                         pijf0(:,:,n1:nsmax)**2*wint(:,:,n1:nsmax))

!  Ignore inner s<0.1 for this diagnostic
      n1 = MAX(ns/20, nsmin)
      wbgradp = SUM(bgradpf(:,:,n1:n2)**2*wint(:,:,n1:n2))
      temp(1) = tnorm_bgradp
      temp(2) = wbgradp
#if defined(MPI_OPT)
      IF (PARSOLVER) THEN
         CALL MPI_ALLREDUCE(MPI_IN_PLACE, temp, 2, MPI_REAL8, MPI_SUM,         &
                            SIESTA_COMM, MPI_ERR)
      END IF
#endif
      tnorm_bgradp = temp(1)
      wbgradp = temp(2)

      IF (tnorm_bgradp .GT. zero) THEN
         bgradp_rms = SQRT(wbgradp/tnorm_bgradp)
         tnorm_bgradp = SQRT(tnorm_bgradp/ns)

!  Ignore first point
         max_bgradp = MAXVAL(bgradpf(:,:,n1:n2))/tnorm_bgradp
         min_bgradp = MINVAL(bgradpf(:,:,n1:n2))/tnorm_bgradp
         temp(1) = max_bgradp
         temp(2) = -min_bgradp
#if defined(MPI_OPT)
         IF (PARSOLVER) THEN
            CALL MPI_ALLREDUCE(MPI_IN_PLACE, temp, 2, MPI_REAL8, MPI_MAX,      &
                               SIESTA_COMM, MPI_ERR)
         END IF
#endif
         max_bgradp = temp(1)
         min_bgradp = -temp(2)
         IF (.not.lasym) THEN
!  Reflect onto u = pi - 2*pi: bgradp has sin parity.
            max_bgradp = MAX(max_bgradp, -min_bgradp)
            min_bgradp = MIN(min_bgradp, -max_bgradp)
         END IF
      END IF

      DEALLOCATE (dpdvijf, dpduijf, stat=istat)

      CALL second0(toff)
      time_bgradp = time_bgradp + (toff - ton)

      END SUBROUTINE

!-------------------------------------------------------------------------------
!>  @brief Plasma beta.
!-------------------------------------------------------------------------------
      SUBROUTINE BETA

      avbeta = wp/wb
 
      END SUBROUTINE
 
      
      SUBROUTINE TFLUX
      USE fourier, ONLY: m0, n0

!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      REAL(dp)           :: part_sum
!-----------------------------------------------
      nsmin=MAX(2,startglobrow); nsmax=MIN(endglobrow,ns)

!Averages over all toroidal cross sections (which should be the same)
      toroidal_flux = SUM(jbsupvmnch(m0,n0,nsmin:nsmax))
#if defined(MPI_OPT)
      IF (PARSOLVER) THEN
         CALL MPI_ALLREDUCE(MPI_IN_PLACE,toroidal_flux,1,MPI_REAL8, MPI_SUM,   &
                            SIESTA_COMM,MPI_ERR)
      END IF
#endif
      toroidal_flux = signjac*twopi*toroidal_flux*hs_i/b_factor  
      
      IF (toroidal_flux0 .EQ. zero) toroidal_flux0 = toroidal_flux

      END SUBROUTINE TFLUX


      SUBROUTINE bdotj (bdotjijhA)
!
!     WRITTEN 08-20-07 BY R. SANCHEZ AS PART OF THE ORNL SIESTA PROJECT (c)
!     
!     PURPOSE: Computes parallel current normalized to the total current (J*B/|JxB|) on the half mesh
!-----------------------------------------------
!   D u m m y   A r g u m e n t
!-----------------------------------------------
        REAL(dp), ALLOCATABLE, DIMENSION(:,:,:) :: bdotjijhA
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
        INTEGER      :: istat
        REAL(dp), ALLOCATABLE, DIMENSION(:,:,:) ::                      &
!          bsupsijh, bsupuijh, bsupvijh,                                 &
          bsubsijh, bsubuijh, bsubvijh, ksubsijh, ksubuijh,             &
          ksubvijh, ksupsijh, ksupuijh, ksupvijh, work1, work2,         &
          work3, work4, work5, work6
        REAL(dp)  :: tnorm, tnorm2, ton, toff
!-----------------------------------------------
        CALL second0(ton)
        CALL ASSERT(SIZE(bsupsijh0,3).EQ.ns, 'bsupsijh0 wrong size in bdotj')
! 
!   Compute covariant field. Used for normalization.
!        
        ALLOCATE(bsubsijh(ntheta,nzeta,ns),                             &
                 bsubuijh(ntheta,nzeta,ns),                             &
                 bsubvijh(ntheta,nzeta,ns), stat=istat)
        CALL ASSERT(istat.EQ.0, 'Allocation #2 failed in BDOTJ')
        bsubsijh=0; bsubvijh=0; bsubuijh=0
      
        CALL tolowerh(bsupsijh0, bsupuijh0, bsupvijh0,                  &
                      bsubsijh, bsubuijh, bsubvijh, 1, ns)

!
!    In initial state, calculate currents (remember they are multiplied by a jacobian!)
!
!
!    Move contravariant currents (*jacob) to the half mesh
!       
        ALLOCATE(ksupsijh(ntheta,nzeta,ns),                             &
                 ksupuijh(ntheta,nzeta,ns),                             &
                 ksupvijh(ntheta,nzeta,ns), stat=istat)
        CALL ASSERT(istat.EQ.0,'Allocation #3 failed in BDOTJ')

        CALL to_half_mesh(ksupsijf0, ksupsijh)
        CALL to_half_mesh(ksupuijf0, ksupuijh)
        CALL to_half_mesh(ksupvijf0, ksupvijh)
!
!   Remove jacobian 
!        
        ksupsijh = ksupsijh/jacobh      
        ksupuijh = ksupuijh/jacobh  
        ksupvijh = ksupvijh/jacobh
!          
!   Get covariant currents
!  
        ALLOCATE(ksubsijh(ntheta,nzeta,ns),                             &
                 ksubuijh(ntheta,nzeta,ns),                             &
                 ksubvijh(ntheta,nzeta,ns), stat=istat)
        CALL ASSERT(istat.EQ.0,'Allocation #4 failed in BDOTJ')
        ksubsijh=0; ksubuijh=0; ksubvijh=0

        CALL tolowerh(ksupsijh, ksupuijh, ksupvijh,                     &
                      ksubsijh, ksubuijh, ksubvijh, 1, ns)
!
!   Compute B*J in the half mesh.
!        
        IF (ALLOCATED(bdotjmnch)) THEN
           IF (SIZE(bdotjmnch,3).NE. ns) THEN
              DEALLOCATE(bdotjmnch)
           END IF
        END IF
        IF (.NOT. ALLOCATED(bdotjmnch)) THEN                             ! Should already be allocated in QUANTITIES       
          ALLOCATE(bdotjmnch(0:mpol,-ntor:ntor,ns), stat=istat)          ! Otherwise, allocate first time called
          CALL ASSERT(istat.EQ.0,'Allocation #5 failed in BDOTJ')
        ENDIF
        IF (lasym) THEN
        IF (ALLOCATED(bdotjmnsh) .AND. SIZE(bdotjmnsh,3).NE. ns)        &
            DEALLOCATE(bdotjmnsh)
        IF (.NOT. ALLOCATED(bdotjmnsh)) THEN                             ! Should already be allocated in QUANTITIES       
          ALLOCATE(bdotjmnsh(0:mpol,-ntor:ntor,ns), stat=istat)          ! Otherwise, allocate first time called
          CALL ASSERT(istat.EQ.0,'Allocation #5 failed in BDOTJ')
        ENDIF
        END IF
        ALLOCATE(bdotjijhA(ntheta,nzeta,ns), stat=istat)
        CALL ASSERT(istat.EQ.0, 'Allocation #6 failed in BDOTJ')
        bdotjijhA=0
        
        bdotjijhA(:,:,3:nsh-1) =                                        &  ! RS: If I include 2 and nsh, it gets huge
          bsupsijh0(:,:,3:nsh-1)*ksubsijh(:,:,3:nsh-1) +                &  ! something funny happens at the boundaries.
          bsupuijh0(:,:,3:nsh-1)*ksubuijh(:,:,3:nsh-1) +                &
          bsupvijh0(:,:,3:nsh-1)*ksubvijh(:,:,3:nsh-1)                       ! B*J = (B^s*J_s + B^u*J_u + B^v*J_v)
! 
!   Compute tnorm = |JXB|^2 for normalization; thus we print out J_parallel/|J_perp|.
!
        ALLOCATE(work1(ntheta,nzeta,ns), work2(ntheta,nzeta,ns),        &
                 work3(ntheta,nzeta,ns), work4(ntheta,nzeta,ns),        &
                 work5(ntheta,nzeta,ns), work6(ntheta,nzeta,ns),        &
                 stat=istat)
        CALL ASSERT(istat.EQ.0,'Allocation #7 failed in BDOTJ')
        work1=0; work2=0; work3=0
        work4=0; work5=0; work6=0
        
        work1 = (bsubuijh*ksubvijh - bsubvijh*ksubuijh)/jacobh              ! (JxB)^s  = (1/sqrt(g))*(K_uB_b-K_vB_u)
        work2 = (bsubvijh*ksubsijh - bsubsijh*ksubvijh)/jacobh              ! (JxB)^u  = .....
        work3 = (bsubsijh*ksubuijh - bsubuijh*ksubsijh)/jacobh              ! (JxB)^v  = .....
        CALL tolowerh(work1, work2, work3, work4, work5, work6, 1, ns)
        
        work1 = work1*work4 + work2*work5 + work3*work6                     ! |JxB|**2
        tnorm = SUM(work1(:,:,3:nsh-1))
        IF (iam.EQ.0 .AND. ANY(work1 .LT. zero))                        &
           PRINT *,'ERROR1 IN BDOTJ, JXB^2 < 0!'
        tnorm = SQRT(ABS(tnorm))
        DEALLOCATE(work3, work4, work5, work6)
        
        work1 = bsupsijh0*bsubsijh + bsupuijh0*bsubuijh +               &  ! |B|**2
                bsupvijh0*bsubvijh
        work2 = ksupsijh*ksubsijh + ksupuijh*ksubuijh +                 &  ! |J|**2
                ksupvijh*ksubvijh
       
        tnorm2 = SUM(work1(:,:,3:nsh-1)*work2(:,:,3:nsh-1))
        work1(:,:,1) = 0;  work2(:,:,1) = 0
        IF (iam.EQ.0 .AND. ANY(work1.LT.zero))                          &
           PRINT *, 'ERROR2 IN BDOTJ: B^2 < 0!'               
        IF (iam.EQ.0 .AND. ANY(work2.LT.zero))                          &
           PRINT *, 'ERROR3 IN BDOTJ: J^2 < 0!'               
        tnorm2 = SQRT(ABS(tnorm2))                                           ! |J|*|B|
        bdotj_rms = SUM(bdotjijhA*bdotjijhA)

        DEALLOCATE (bsubsijh, bsubuijh, bsubvijh, ksupsijh, ksupuijh,   &
                    ksupvijh, ksubsijh, ksubuijh, ksubvijh,             &      
                    work1, work2)

        tnorm = MAX(tnorm, EPSILON(tnorm))
        tnorm2= MAX(tnorm2,EPSILON(tnorm2))
        bdotj2_rms = SQRT(ABS(bdotj_rms))/tnorm2                           ! RMS of bdotj/|J|*|B|  
        bdotj_rms  = SQRT(ABS(bdotj_rms))/tnorm                            ! RMS of bdotj/|JxB|  
        CALL fourier_context%tomnsp(bdotjijhA, bdotjmnch, f_cos)                  ! Keep harmonics of BDOTJ for output
        bdotjmnch = bdotjmnch/tnorm
        IF (lasym) THEN
        CALL fourier_context%tomnsp(bdotjijhA, bdotjmnsh, f_cos)                  ! Keep harmonics of BDOTJ for output
        bdotjmnsh = bdotjmnsh/tnorm
        END IF
        CALL second0(toff)
        time_bdotj = time_bdotj+(toff-ton)

        bdotjijh = bdotjijhA

        END SUBROUTINE bdotj

      SUBROUTINE bdotj_par
      USE island_params, ONLY: fourier_context

      IMPLICIT NONE

!     
!     WRITTEN 08-20-07 BY R. SANCHEZ AS PART OF THE ORNL SIESTA PROJECT (c)
!     UPDATED 03-21-13 BY S. HIRSHMAN FOR PARALLEL EXECUTION
!     
!     PURPOSE: Computes parallel current normalized to the total current (J*B/|JxB|) on the half mesh
!
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER      :: istat, n1, n2, ns_span
      REAL(dp), ALLOCATABLE, DIMENSION(:,:,:) :: bsupsijh,              &
          bsupuijh, bsupvijh,                                           &
          bsubsijh, bsubuijh, bsubvijh, ksubsijh, ksubuijh,             &
          ksubvijh, ksupsijh, ksupuijh, ksupvijh, work1, work2,         &
          work3, work4, work5, work6
      REAL(dp)  :: tnorm, tnorm2, temp(3), ton, toff
!-----------------------------------------------
      CALL second0(ton)

      CALL ASSERT(lcurr_init, 'MUST CALL init_state(.TRUE.) BEFORE bdotj_par')
      nsmin=MAX(1,startglobrow); nsmax=MIN(ns,endglobrow)
      ns_span=nsmax-nsmin+1
! 
!     Compute contravariant components of magnetic field
!       
        ALLOCATE(bsupsijh(ntheta,nzeta,nsmin:nsmax),                    &
                 bsupuijh(ntheta,nzeta,nsmin:nsmax),                    &
                 bsupvijh(ntheta,nzeta,nsmin:nsmax), stat=istat)
        CALL ASSERT(istat.EQ.0,'Allocation #1 failed in BDOTJ_PAR')

        CALL fourier_context%toijsp(jbsupsmnsh(:,:,nsmin:nsmax), bsupsijh,     &
       &                            f_none, f_sin)
        CALL fourier_context%toijsp(jbsupumnch(:,:,nsmin:nsmax), bsupuijh,     &
                                    f_none, f_cos)
        CALL fourier_context%toijsp(jbsupvmnch(:,:,nsmin:nsmax), bsupvijh,     &
                                    f_none, f_cos)
        IF (lasym) THEN
           CALL fourier_context%toijsp(jbsupsmnch(:,:,nsmin:nsmax), bsupsijh,  &
                                       f_sum, f_cos)
           CALL fourier_context%toijsp(jbsupumnsh(:,:,nsmin:nsmax), bsupuijh,  &
                                       f_sum, f_sin)
           CALL fourier_context%toijsp(jbsupvmnsh(:,:,nsmin:nsmax), bsupvijh,  &
                                       f_sum, f_sin)
        END IF

        bsupsijh = bsupsijh/jacobh(:,:,nsmin:nsmax)
        bsupuijh = bsupuijh/jacobh(:,:,nsmin:nsmax)
        bsupvijh = bsupvijh/jacobh(:,:,nsmin:nsmax)
!
!   Compute covariant field. Used for normalization.
!        
        ALLOCATE(bsubsijh(ntheta,nzeta,nsmin:nsmax),                    &
                 bsubuijh(ntheta,nzeta,nsmin:nsmax),                    &
                 bsubvijh(ntheta,nzeta,nsmin:nsmax), stat=istat)
        CALL ASSERT(istat.EQ.0,'Allocation #2 failed in BDOTJ_PAR')
      
        CALL tolowerh(bsupsijh, bsupuijh, bsupvijh,                     &
                      bsubsijh, bsubuijh, bsubvijh, nsmin, nsmax)
!
!    In initial state, calculate currents (remember they are multiplied by a jacobian!)
!
!
!    Move contravariant currents (*jacob) to the half mesh
!       
        nsmin=MAX(1,startglobrow-1)
        ALLOCATE(ksupsijh(ntheta,nzeta,nsmin:nsmax),                    &
                 ksupuijh(ntheta,nzeta,nsmin:nsmax),                    &
                 ksupvijh(ntheta,nzeta,nsmin:nsmax), stat=istat)
        CALL ASSERT(istat.EQ.0,'Allocation #3 failed in BDOTJ_PAR')

        CALL to_half_mesh(ksupsijf0, ksupsijh)
        CALL to_half_mesh(ksupuijf0, ksupuijh)
        CALL to_half_mesh(ksupvijf0, ksupvijh)

!
!   Remove jacobian 
!        
        ksupsijh = ksupsijh/jacobh(:,:,nsmin:nsmax)
        ksupuijh = ksupuijh/jacobh(:,:,nsmin:nsmax)
        ksupvijh = ksupvijh/jacobh(:,:,nsmin:nsmax)

!          
!   Get covariant currents
!  
        ALLOCATE(ksubsijh(ntheta,nzeta,nsmin:nsmax),                    &
                 ksubuijh(ntheta,nzeta,nsmin:nsmax),                    &
                 ksubvijh(ntheta,nzeta,nsmin:nsmax), stat=istat)
        CALL ASSERT(istat.EQ.0,'Allocation #4 failed in BDOTJ_PAR')

        CALL tolowerh(ksupsijh, ksupuijh, ksupvijh,                     &
                      ksubsijh, ksubuijh, ksubvijh, nsmin, nsmax)
        nsmin=MAX(1,startglobrow)

!
!   Compute B*J in the half mesh.
!      
        IF (ALLOCATED(bdotjmnch)) THEN
           IF (SIZE(bdotjmnch,3).NE.ns_span) THEN
              DEALLOCATE(bdotjmnch)
           END IF
        END IF
        IF (.NOT. ALLOCATED(bdotjmnch)) THEN                            
          ALLOCATE(bdotjmnch(0:mpol,-ntor:ntor,nsmin:nsmax), stat=istat)
          CALL ASSERT(istat.EQ.0,'Allocation #5A failed in BDOTJ_PAR')
        ENDIF
        IF (lasym) THEN
        IF (ALLOCATED(bdotjmnsh) .AND. SIZE(bdotjmnsh,3).NE.ns_span)    &
            DEALLOCATE(bdotjmnsh)
        IF (.NOT. ALLOCATED(bdotjmnsh)) THEN                            
          ALLOCATE(bdotjmnsh(0:mpol,-ntor:ntor,nsmin:nsmax), stat=istat)
          CALL ASSERT(istat.EQ.0,'Allocation #5B failed in BDOTJ_PAR')
        ENDIF
        ENDIF

        ALLOCATE(bdotjijh(ntheta,nzeta,nsmin:nsmax), stat=istat)
        CALL ASSERT(istat.EQ.0,'Allocation #6 failed in BDOTJ_PAR')
        bdotjijh=0
        
        n1 = MAX(3,nsmin); n2 = MIN(nsh-1, nsmax)
        bdotjijh(:,:,n1:n2) =                                           &  ! RS: If include 2 and nsh, it gets huge
          bsupsijh(:,:,n1:n2)*ksubsijh(:,:,n1:n2) +                     &  ! something funny happens at the boundaries.
          bsupuijh(:,:,n1:n2)*ksubuijh(:,:,n1:n2) +                     &
          bsupvijh(:,:,n1:n2)*ksubvijh(:,:,n1:n2)                           ! B*J = (B^s*J_s + B^u*J_u + B^v*J_v)
        
! 
!   Compute tnorm = |JXB|^2 for normalization; thus we print out J_parallel/|J_perp|.
!
        ALLOCATE(work1(ntheta,nzeta,nsmin:nsmax),                       &
                 work2(ntheta,nzeta,nsmin:nsmax),                       &
                 work3(ntheta,nzeta,nsmin:nsmax),                       &
                 work4(ntheta,nzeta,nsmin:nsmax),                       &
                 work5(ntheta,nzeta,nsmin:nsmax),                       &
                 work6(ntheta,nzeta,nsmin:nsmax),                       &
                 stat=istat)
        CALL ASSERT(istat.EQ.0,'Allocation #7 failed in BDOTJ_PAR')
        
        work1 = bsubuijh*ksubvijh(:,:,nsmin:nsmax)                      &
              - bsubvijh*ksubuijh(:,:,nsmin:nsmax)                        ! (JxB)^s  = (1/sqrt(g))*(K_uB_b-K_vB_u)
        work2 = bsubvijh*ksubsijh(:,:,nsmin:nsmax)                      & 
              - bsubsijh*ksubvijh(:,:,nsmin:nsmax)                        ! (JxB)^u  = .....
        work3 = bsubsijh*ksubuijh(:,:,nsmin:nsmax)                      &
              - bsubuijh*ksubsijh(:,:,nsmin:nsmax)                        ! (JxB)^v  = .....
        work1 = work1/jacobh(:,:,nsmin:nsmax)
        work2 = work2/jacobh(:,:,nsmin:nsmax)
        work3 = work3/jacobh(:,:,nsmin:nsmax)
        CALL tolowerh(work1, work2, work3,                              &
                      work4, work5, work6, nsmin, nsmax)
        
        work1 = work1*work4 + work2*work5 + work3*work6                     ! |JxB|**2
        tnorm = SUM(work1(:,:,n1:n2))
        IF (nsmin .EQ. 1) work1(:,:,1) = 0
        IF (iam.EQ.0 .AND. ANY(work1 .LT. zero))                        &
           PRINT *,'ERROR1 IN BDOTJ_PAR, JXB^2 < 0!'
        DEALLOCATE(work3, work4, work5, work6)
        
        work1 = bsupsijh*bsubsijh + bsupuijh*bsubuijh +                 & ! |B|**2
                bsupvijh*bsubvijh
        work2 = ksupsijh(:,:,nsmin:nsmax)*ksubsijh(:,:,nsmin:nsmax)     &
              + ksupuijh(:,:,nsmin:nsmax)*ksubuijh(:,:,nsmin:nsmax)     &
              + ksupvijh(:,:,nsmin:nsmax)*ksubvijh(:,:,nsmin:nsmax)       ! |J|**2
       
        IF (nsmin .EQ. 1) THEN
           work1(:,:,1) = 0;  work2(:,:,1) = 0
        END IF
        IF (iam.EQ.0 .AND. ANY(work1.LT.zero))                          &
           PRINT *, 'ERROR2 IN BDOTJ_PAR: B^2 < 0!'               
        IF (iam.EQ.0 .AND. ANY(work2.LT.zero))                          &
           PRINT *, 'ERROR3 IN BDOTJ_PAR: J^2 < 0!'               
        tnorm2 = SUM(work1(:,:,n1:n2)*work2(:,:,n1:n2))
        bdotj_rms = SUM(bdotjijh(:,:,n1:n2)*bdotjijh(:,:,n1:n2))
#if defined(MPI_OPT)
        temp(1)=tnorm; temp(2)=tnorm2; temp(3)=bdotj_rms
        CALL MPI_ALLREDUCE(MPI_IN_PLACE,temp,3,MPI_REAL8, MPI_SUM,                  &
                           SIESTA_COMM,MPI_ERR)
        tnorm  = SQRT(ABS(temp(1)))
        tnorm2 = SQRT(ABS(temp(2)))                                         ! |J|*|B|
        bdotj_rms = temp(3)
#else
        tnorm = SQRT(ABS(tnorm)); tnorm2 = SQRT(ABS(tnorm2))
#endif
        DEALLOCATE (bsupsijh, bsupuijh, bsupvijh, bsubsijh,             &
                    bsubuijh, bsubvijh, ksupsijh, ksupuijh,             &
                    ksupvijh, ksubsijh, ksubuijh, ksubvijh,             &      
                    work1, work2)


        tnorm2= MAX(tnorm2,EPSILON(tnorm2))
        bdotj2_rms = SQRT(ABS(bdotj_rms))/tnorm2                          ! RMS of bdotj/|J|*|B|  
        bdotj_rms  = SQRT(ABS(bdotj_rms))/tnorm                           ! RMS of bdotj/|JxB|  

        bdotjijh = bdotjijh/tnorm
        
        CALL fourier_context%tomnsp(bdotjijh, bdotjmnch, f_cos)                       ! Keep harmonics of BDOTJ for output
        IF (lasym) THEN
           CALL fourier_context%tomnsp(bdotjijh, bdotjmnsh, f_sin)                       ! Keep harmonics of BDOTJ for output
        END IF

        DEALLOCATE (bdotjijh)

        CALL second0(toff)
        time_bdotj = time_bdotj+(toff-ton)

        END SUBROUTINE bdotj_par
      
!-------------------------------------------------------------------------------
!>  @brief Compute the divergence of J.
!>
!>  Computes divergence of the current normalized to the total current on the
!>  half mesh.
!>
!>  @param[in] ns_min Lower radial index.
!>  @param[in] ns_max Upper radial index.
!-------------------------------------------------------------------------------
      SUBROUTINE divj(ns_min, ns_max)
      USE island_params, ONLY: ohs=>ohs_i

      IMPLICIT NONE

!  Declare Arguments
      INTEGER, INTENT(in)                     :: ns_min
      INTEGER, INTENT(in)                     :: ns_max

!  Local variables
      INTEGER                                 :: m
      INTEGER                                 :: n
      INTEGER                                 :: istat
      INTEGER                                 :: n1
      INTEGER                                 :: n2
      REAL(dp), DIMENSION(:,:,:), ALLOCATABLE :: temp1
      REAL(dp), DIMENSION(:,:,:), ALLOCATABLE :: temp2
      REAL(dp)                                :: tnorm
      REAL(dp), DIMENSION(2)                  :: temp
      REAL(dp)                                :: ton
      REAL(dp)                                :: toff

!  Start of executable code.
      CALL second0(ton)

      CALL ASSERT(lcurr_init,'MUST CALL init_state(.TRUE.) BEFORE divJ_par')

      nsmin = ns_min
      nsmax = ns_max

      IF (ALLOCATED(divjmnsh)) THEN
         DEALLOCATE(divjmnsh)
      END IF
      IF (ALLOCATED(divjmnch)) THEN
         DEALLOCATE(divjmnch)
      END IF
      ALLOCATE(divjmnsh(0:mpol,-ntor:ntor,nsmin:nsmax),                        &
               divjmnch(0:mpol,-ntor:ntor,nsmin:nsmax), stat=istat)
      CALL assert_eq(istat, 0, 'Allocation failed in DIVJ')
      divjmnsh = 0
      divjmnch = 0

!  Compute divergence of J. It is a SINE parity quantity.
      n1 = MAX(3, nsmin)
      n2 = MIN(ns-1, nsmax)
      DO m = 0, mpol
         DO n = -ntor,ntor
            divjmnsh(m,n,n1:n2) =                                              &
                    -m*(ksupumncf(m,n,n1:n2) + ksupumncf(m,n,n1-1:n2-1))/2     &
               - n*nfp*(ksupvmncf(m,n,n1:n2) + ksupvmncf(m,n,n1-1:n2-1))/2     &
               +   ohs*(ksupsmnsf(m,n,n1:n2) - ksupsmnsf(m,n,n1-1:n2-1))
            IF (lasym) THEN
                divjmnch(m,n,n1:n2) =                                          &
                         m*(ksupumnsf(m,n,n1:n2) + ksupumnsf(m,n,n1-1:n2-1))/2 &
                   + n*nfp*(ksupvmnsf(m,n,n1:n2) + ksupvmnsf(m,n,n1-1:n2-1))/2 &
                   +   ohs*(ksupsmncf(m,n,n1:n2) - ksupsmncf(m,n,n1-1:n2-1))
            END IF
         END DO
      END DO
!
!    Compute covariant components of current (times jacobian). Used for normalization. 
!  
      nsmin = MAX(1, startglobrow - 1)
      nsmax = MIN(endglobrow + 1, ns)
      ALLOCATE(temp1(ntheta,nzeta,nsmin:nsmax),                                &
               temp2(ntheta,nzeta,nsmin:nsmax), stat=istat)
      CALL assert_eq(istat, 0, 'Allocation failed in DIVJ')
      temp1 = 0
      temp2 = 0

!  Norm: |J|^2
      temp1 = ksupsijf0(:,:,nsmin:nsmax)*ksubsijf(:,:,nsmin:nsmax)             &
            + ksupuijf0(:,:,nsmin:nsmax)*ksubuijf(:,:,nsmin:nsmax)             &
            + ksupvijf0(:,:,nsmin:nsmax)*ksubvijf(:,:,nsmin:nsmax)

      CALL to_half_mesh(temp1, temp2)

      temp2 = temp2/(jacobh(:,:,nsmin:nsmax)*jacobh(:,:,nsmin:nsmax))

      nsmin = MAX(1, startglobrow)
      nsmax = MIN(endglobrow, ns)
      tnorm = SUM(temp2(:,:,nsmin:nsmax))

      DEALLOCATE(temp1, temp2, stat=istat)

      divj_rms = SUM(divjmnsh(:,:,n1:n2)**2                                    &
               +     divjmnch(:,:,n1:n2)**2)
#if defined(MPI_OPT)
      IF (PARSOLVER) THEN
         temp(1) = divj_rms
         temp(2) = tnorm
         CALL MPI_ALLREDUCE(MPI_IN_PLACE,temp,2,MPI_REAL8, MPI_SUM,            &
                            SIESTA_COMM,MPI_ERR)
         divj_rms = temp(1)
         tnorm = temp(2)
      END IF
#endif

!  Compute rms of divergence of J
      IF (tnorm .ge. zero) THEN
         divjmnsh = divjmnsh/SQRT(tnorm)
         divjmnch = divjmnch/SQRT(tnorm)
         divj_rms = SQRT(divj_rms/tnorm)
      END IF

      CALL second0(toff)
      time_divj = time_divj + (toff - ton)

      END SUBROUTINE

!-------------------------------------------------------------------------------
!>  @brief Deallocate diagnostic buffers.
!-------------------------------------------------------------------------------
      SUBROUTINE dealloc_diagnostics

      IF (ALLOCATED(divbmnsf)) THEN
         DEALLOCATE(divbmnsf)
      END IF
      IF (ALLOCATED(divbmncf)) THEN
         DEALLOCATE(divbmncf)
      END IF
      IF (ALLOCATED(divjmnsh)) THEN
         DEALLOCATE(divjmnsh)
      END IF
      IF (ALLOCATED(divjmnch)) THEN
         DEALLOCATE(divjmnch)
      END IF
      IF (ALLOCATED(bdotjmnch)) THEN
         DEALLOCATE(bdotjmnch)
      END IF
      IF (ALLOCATED(bdotjmnsh)) THEN
         DEALLOCATE(bdotjmnsh)
      END IF

      END SUBROUTINE

      END MODULE
