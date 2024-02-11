!*******************************************************************************
!>  @file vmec_info.f90
!>  @brief Contains module @ref vmec_info
!
!  Note separating the Doxygen comment block here so detailed decription is
!  found in the Module not the file.
!
!>  Holds memory references for splined and island VMEC quantities.
!*******************************************************************************
      MODULE vmec_info
      USE stel_kinds
      USE read_wout_mod, ns_vmec=>ns, mpol_vmec=>mpol, ntor_vmec=>ntor,        &
                         rmnc_vmec=>rmnc, zmns_vmec=>zmns, lmns_vmec=>lmns,    &
                         xm_vmec=>xm, xn_vmec=>xn, chipf_vmec=>chipf,          &  ! MRC 4/1/2016
                         rmns_vmec=>rmns, zmnc_vmec=>zmnc, lmnc_vmec=>lmnc,    &  ! MRC 12/1/2016
                         phipf_vmec=>phipf, presf_vmec=>presf, nfp_vmec=>nfp,  &
                         wb_vmec=>wb, wp_vmec=>wp, gamma_vmec=>gamma,          &
                         volume_vmec=>volume, raxis_vmec=>raxis,               &
                         lasym_vmec=>lasym, iasym_vmec=>iasym,                 &
                         vmec_curtor=>Itor
      USE fourier, ONLY: f_cos, f_sin

      IMPLICIT NONE

!*******************************************************************************
!  vmec_info module variables
!*******************************************************************************
!>  Splined Fourier coefficents for VMEC R cosine parity.
      REAL (dp), DIMENSION(:,:), ALLOCATABLE   :: rmnc_spline
!>  Splined Fourier coefficents for VMEC Z sine parity.
      REAL (dp), DIMENSION(:,:), ALLOCATABLE   :: zmns_spline
!>  Splined Fourier coefficents for VMEC lambda sine parity.
      REAL (dp), DIMENSION(:,:), ALLOCATABLE   :: lmns_spline
!>  Splined Fourier coefficents for VMEC R sine parity.
      REAL (dp), DIMENSION(:,:), ALLOCATABLE   :: rmns_spline
!>  Splined Fourier coefficents for VMEC Z cosine parity.
      REAL (dp), DIMENSION(:,:), ALLOCATABLE   :: zmnc_spline
!>  Splined Fourier coefficents for VMEC lambda cosine parity.
      REAL (dp), DIMENSION(:,:), ALLOCATABLE   :: lmnc_spline
!>  Splined Fourier coefficents for VMEC J^u cosine parity.
      REAL (dp), DIMENSION(:,:), ALLOCATABLE   :: currumnc_spline
!>  Splined Fourier coefficents for VMEC J^v cosine parity.
      REAL (dp), DIMENSION(:,:), ALLOCATABLE   :: currvmnc_spline
!>  Splined Fourier coefficents for VMEC J^u sine parity.
      REAL (dp), DIMENSION(:,:), ALLOCATABLE   :: currumns_spline
!>  Splined Fourier coefficents for VMEC J^v sine parity.
      REAL (dp), DIMENSION(:,:), ALLOCATABLE   :: currvmns_spline

!>  Fourier coefficents for VMEC R cosine parity.
      REAL (dp), DIMENSION(:,:,:), ALLOCATABLE :: rmnc_i
!>  Fourier coefficents for VMEC Z sine parity.
      REAL (dp), DIMENSION(:,:,:), ALLOCATABLE :: zmns_i
!>  Fourier coefficents for VMEC lambda sine parity.
      REAL (dp), DIMENSION(:,:,:), ALLOCATABLE :: lmns_i
!>  Fourier coefficents for VMEC R sine parity.
      REAL (dp), DIMENSION(:,:,:), ALLOCATABLE :: rmns_i
!>  Fourier coefficents for VMEC Z cosine parity.
      REAL (dp), DIMENSION(:,:,:), ALLOCATABLE :: zmnc_i
!>  Fourier coefficents for VMEC lambda cosine parity.
      REAL (dp), DIMENSION(:,:,:), ALLOCATABLE :: lmnc_i
!>  Fourier coefficents for VMEC J^u cosine parity.
      REAL (dp), DIMENSION(:,:,:), ALLOCATABLE :: jcurrumnc
!>  Fourier coefficents for VMEC J^v cosine parity.
      REAL (dp), DIMENSION(:,:,:), ALLOCATABLE :: jcurrvmnc
!>  Fourier coefficents for VMEC J^u sine parity.
      REAL (dp), DIMENSION(:,:,:), ALLOCATABLE :: jcurrumns
!>  Fourier coefficents for VMEC J^v sine parity.
      REAL (dp), DIMENSION(:,:,:), ALLOCATABLE :: jcurrvmns

!*******************************************************************************
!  vmec_info module parameters
!*******************************************************************************
!>     SET iflipj -> -1 FOR POSITIVE JACOBIAN SYSTEM
      INTEGER, PARAMETER                       :: iflipj = 1
!>  Spline lower boundary condition. Uses last four points.
      REAL (dp), PARAMETER                     :: lower_b = -1.e30_dp
!>  Spline upper boundary condition. Uses last four points.
      REAL (dp), PARAMETER                     :: upper_b = -1.e30_dp

      CONTAINS
!*******************************************************************************
!  CONSTRUCTION SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  @brief Allocate spline buffers.
!>
!>  Allocated memory of VMEC quantities spline buffers. Spline arrays take the
!>  dimensions of the total number of Fourier modes and the size of the radial
!>  grid. Lambda is allocated with an extra radial grid point.
!>
!>  @param[in] mnmax   Total number of Fourier modes.
!>  @param[in] ns      Total number of radial surfaces.
!>  @param[in] is_asym Flag to mark if not using stellarator symmetry.
!-------------------------------------------------------------------------------
      SUBROUTINE vmec_info_construct_spline(mnmax, ns, is_asym)

      IMPLICIT NONE

!  Declare Arguments
      INTEGER, INTENT(in) :: mnmax
      INTEGER, INTENT(in) :: ns
      LOGICAL, INTENT(in) :: is_asym

!  Start of executable code

!  Stellarator symmetric quantities.
      ALLOCATE(rmnc_spline(mnmax,ns))
      ALLOCATE(zmns_spline(mnmax,ns))
      ALLOCATE(lmns_spline(mnmax,ns + 1))
      ALLOCATE(currumnc_spline(mnmax,ns))
      ALLOCATE(currvmnc_spline(mnmax,ns))

!  Asymmetric quantities.
      IF (is_asym) THEN
         ALLOCATE(rmns_spline(mnmax,ns))
         ALLOCATE(zmnc_spline(mnmax,ns))
         ALLOCATE(lmnc_spline(mnmax,ns + 1))
         ALLOCATE(currumns_spline(mnmax,ns))
         ALLOCATE(currvmns_spline(mnmax,ns))
      END IF

      END SUBROUTINE

!-------------------------------------------------------------------------------
!>  @brief Allocate island buffers.
!>
!>  Allocated memory of VMEC quantities island buffers. Island arrays take the
!>  dimensions of 0:mpol, -ntor:ntor and the size of the radial grid. Lambda is
!>  allocated with an extra radial grid point.
!>
!>  @param[in] mpol    Number of poloidal modes.
!>  @param[in] ntor    Number of Toroidal modes.
!>  @param[in] ns      Number of radial surfaces.
!>  @param[in] is_asym Flag to mark if not using stellarator symmetry.
!-------------------------------------------------------------------------------
      SUBROUTINE vmec_info_construct_island(mpol, ntor, ns, is_asym)

      IMPLICIT NONE

!  Declare Arguments
      INTEGER, INTENT(in) :: mpol
      INTEGER, INTENT(in) :: ntor
      INTEGER, INTENT(in) :: ns
      LOGICAL, INTENT(in) :: is_asym

!  Start of executable code

!  Stellarator symmetric quantities.
      ALLOCATE(rmnc_i(0:mpol,-ntor:ntor,ns))
      ALLOCATE(zmns_i(0:mpol,-ntor:ntor,ns))
      ALLOCATE(lmns_i(0:mpol,-ntor:ntor,ns + 1))
      ALLOCATE(jcurrumnc(0:mpol,-ntor:ntor,ns))
      ALLOCATE(jcurrvmnc(0:mpol,-ntor:ntor,ns))

      rmnc_i = 0
      zmns_i = 0
      lmns_i = 0
      jcurrumnc = 0
      jcurrvmnc = 0

!  Asymmetric quantities.
      IF (is_asym) THEN
         ALLOCATE(rmns_i(0:mpol,-ntor:ntor,ns))
         ALLOCATE(zmnc_i(0:mpol,-ntor:ntor,ns))
         ALLOCATE(lmnc_i(0:mpol,-ntor:ntor,ns + 1))
         ALLOCATE(jcurrumns(0:mpol,-ntor:ntor,ns))
         ALLOCATE(jcurrvmns(0:mpol,-ntor:ntor,ns))

         rmns_i = 0
         zmnc_i = 0
         lmnc_i = 0
         jcurrumns = 0
         jcurrvmns = 0
      END IF

      END SUBROUTINE

!*******************************************************************************
!  DESTRUCTION SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  @brief Deallocate spline buffers.
!>
!>  Deallocated memory of VMEC quantities spline buffers.
!-------------------------------------------------------------------------------
      SUBROUTINE vmec_info_destruct_spline

      IMPLICIT NONE

!  Start of executable code

!  Stellarator symmetric quantities.
      DEALLOCATE(rmnc_spline)
      DEALLOCATE(zmns_spline)
      DEALLOCATE(lmns_spline)
      DEALLOCATE(currumnc_spline)
      DEALLOCATE(currvmnc_spline)

!  Asymmetric quantities.
      IF (ALLOCATED(rmns_spline)) THEN
         DEALLOCATE(rmns_spline)
         DEALLOCATE(zmnc_spline)
         DEALLOCATE(lmnc_spline)
         DEALLOCATE(currumns_spline)
         DEALLOCATE(currvmns_spline)
      END IF

      END SUBROUTINE

!-------------------------------------------------------------------------------
!>  @brief Deallocate island buffers.
!>
!>  Deallocated memory of VMEC quantities island buffers.
!-------------------------------------------------------------------------------
      SUBROUTINE vmec_info_destruct_island

      IMPLICIT NONE

!  Start of executable code

!  Stellarator symmetric quantities.
      DEALLOCATE(rmnc_i)
      DEALLOCATE(zmns_i)
      DEALLOCATE(lmns_i)
      DEALLOCATE(jcurrumnc)
      DEALLOCATE(jcurrvmnc)

!  Asymmetric quantities.
      IF (ALLOCATED(rmns_i)) THEN
         DEALLOCATE(rmns_i)
         DEALLOCATE(zmnc_i)
         DEALLOCATE(lmnc_i)
         DEALLOCATE(jcurrumns)
         DEALLOCATE(jcurrvmns)
      END IF

      END SUBROUTINE

!*******************************************************************************
!  SETTER SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  @brief Set the vmec quantities from a wout file.
!>
!>  @param[in] wout_file     Filename to load vmec quantities from.
!>  @param[in] ns_in         Number of SIESTA radial surfaces.
!>  @param[in] mpol_in       Number of SIESTA poloidal modes.
!>  @param[in] ntor_in       Number of SIESTA toroidal modes.
!>  @param[in] nfp_in        Number of SIESTA field periods.
!>  @param[in] ntor_modes_in SIESTA Toroidal mode numbers.
!-------------------------------------------------------------------------------
      SUBROUTINE vmec_info_set_wout(wout_file, ns_in, mpol_in, ntor_in,        &
     &                              nfp_in, ntor_modes_in)
      USE descriptor_mod, ONLY: iam
      USE v3_utilities, ONLY: assert_eq, assert
      USE island_params
      USE stel_constants, ONLY: twopi, zero, mu0
      USE island_params, hs=>hs_i, ns=>ns_i
      USE shared_data, ONLY: lasym

      IMPLICIT NONE

!  Declare Arguments
      CHARACTER(len=*), INTENT(in)                     :: wout_file
      INTEGER, INTENT(IN)                              :: ns_in
      INTEGER, INTENT(IN)                              :: mpol_in
      INTEGER, INTENT(IN)                              :: ntor_in
      INTEGER, INTENT(IN)                              :: nfp_in
      INTEGER, DIMENSION(-ntor_in:ntor_in), INTENT(in) :: ntor_modes_in

!  Local variables
      INTEGER                                          :: istat
      INTEGER                                          :: js
      REAL (dp)                                        :: temp

!  Start of executable code

!  Load wout file.
      CALL read_wout_file(wout_file, istat)
      CALL assert_eq(0, istat, 'Read-wout error in vmec_info_set_wout')

      IF (nfp_in .lt. 1) THEN
         nfp_i = nfp_vmec
      ELSE
         nfp_i = nfp_in
      END IF

      wb_i  = (twopi*twopi)*wb_vmec
      wp_i  = (twopi*twopi)*wp_vmec
      volume_i = volume_vmec
      rmajor_i = raxis_vmec(0,1)

      CALL assert(wb_i .gt. zero, 'wb_vmec = 0!')

      ns = ns_in
      nsh = ns_in - 1

      fourier_context => fourier_class(mpol_in, ntor_in, nu_i, nv_i, nfp_i,    &
                                       lasym, ntor_modes_in)

!  Spline r, z, and l Fourier components in s from original vmec mesh (s ~ phi,
!  ns_vmec points) to a "polar" mesh s ~ sqrt(phi) for better axis resolution.
      CALL vmec_info_construct_spline(mnmax, ns, lasym)
      CALL vmec_info_set_RZL_splines(rmnc_vmec, zmns_vmec, lmns_vmec,          &
                                     currumnc, currvmnc,                       &
                                     rmnc_spline, zmns_spline, lmns_spline,    &
                                     currumnc_spline, currvmnc_spline, f_cos)
      IF (lasym) THEN
         CALL vmec_info_set_RZL_splines(rmns_vmec, zmnc_vmec, lmnc_vmec,       &
                                        currumns, currvmns,                    &
                                        rmns_spline, zmnc_spline, lmnc_spline, &
                                        currumns_spline, currvmns_spline,      &
                                        f_sin)
      END IF

!  Spline 1-D arrays: careful -> convert phipf VMEC and multiply by
!  ds-vmec/ds-island, since phipf_i = d(PHI)/ds-island
      ALLOCATE(phipf_i(ns), chipf_i(ns), presf_i(ns), stat=istat)
      CALL vmec_info_spline_oned_array(chipf_vmec, chipf_i, istat)
      CALL vmec_info_spline_oned_array(phipf_vmec, phipf_i, istat)
      presf_vmec = mu0*presf_vmec
      CALL vmec_info_spline_oned_array(presf_vmec, presf_i, istat)

!  Pessure should never be negative.
      WHERE (presf_i .lt. 0)
         presf_i = 0
      END WHERE

!  Scale phipf_i and convert to sqrt(flux) mesh by multiplying by 2*s
      phipf_i = phipf_i/twopi
      chipf_i = chipf_i/twopi

!  Mapping s_vmec = (s_siesta)^2, so d(s_vmec)/d(s_siesta) = 2 s_siesta, where
!  s_siesta = hs_i*(js-1)
      DO js = 1, ns
         phipf_i(js) = 2.0*hs*(js - 1)*phipf_i(js)
         chipf_i(js) = 2.0*hs*(js - 1)*chipf_i(js)
      END DO

      IF (iam .EQ. 0) THEN
         temp = twopi*hs*(SUM(phipf_i(2:ns)) + SUM(phipf_i(1:ns-1)))/2
         WRITE (33, 1000) volume_i, 1.E-6_dp*(wb_i+wp_i/(gamma-1))/mu0,        &
                          1.E-6_dp*wp_i/mu0, temp, wp_i/wb_i,                  &
                          chipf_i(2)/phipf_i(2), chipf_i(ns)/phipf_i(ns),      &
                          gamma
      END IF

      CALL vmec_info_set_RZL

!  Construct R
1000  FORMAT(/,' INITIAL PARAMETERS (FROM VMEC)',/,                            &
               ' PHYSICS QUANTITIES',/,                                        &
               ' PLASMA VOLUME (M^3): ',1pE12.4,                               &
               ' TOTAL MHD ENERGY (MJ):   ', 1pE16.8,/,                        &
               ' THERMAL ENERGY (MJ): ', 1pE12.4,                              &
               ' EDGE TOROIDAL FLUX (Wb): ', 1pE12.4,/,                        &
               ' DIMENSIONLESS QUANTITIES',/,                                  &
               ' <BETA>: ',1pE12.4, ' IOTA(0): ',1pE12.4,                      &
               ' IOTA(1) ',1pE12.4, ' GAMMA: ',1pE12.4,/, 21('-'),/)

      END SUBROUTINE

!-------------------------------------------------------------------------------
!>  @brief Splines the VMEC quantities.
!>
!>  Recasts the VMEC quantities to the SIESTA mesh.
!>
!>  @param[in]    rmn_vmec       VMEC R Fourier amplitudes.
!>  @param[in]    zmn_vmec       VMEC Z Fourier amplitudes.
!>  @param[in]    lmn_vmec       VMEC lambda Fourier amplitudes.
!>  @param[in]    currumn_vmec   VMEC J^u Fourier amplitudes.
!>  @param[in]    currvmn_vmec   VMEC J^v Fourier amplitudes.
!>  @param[inout] rmn_spline     Spline R Fourier amplitudes.
!>  @param[inout] zmn_spline     Spline Z Fourier amplitudes.
!>  @param[inout] lmn_spline     Spline lambda Fourier amplitudes.
!>  @param[inout] currumn_spline Spline J^u Fourier amplitudes.
!>  @param[inout] currvmn_spline Spline J^v Fourier amplitudes.
!>  @param[in]    parity         Parity of the Fourier amplitudes.
!-------------------------------------------------------------------------------
      SUBROUTINE vmec_info_set_RZL_splines(rmn_vmec, zmn_vmec, lmn_vmec,       &
                                           currumn_vmec, currvmn_vmec,         &
                                           rmn_spline, zmn_spline, lmn_spline, &
                                           currumn_spline, currvmn_spline,     &
                                           parity)
      USE shared_data, ONLY: jsupvdotA
      USE stel_constants, ONLY: mu0
      USE v3_utilities, ONLY: assert, assert_eq

      IMPLICIT NONE

!  Declare Arguments
      REAL (dp), DIMENSION(:,:), INTENT(in)                 :: rmn_vmec
      REAL (dp), DIMENSION(:,:), INTENT(in)                 :: zmn_vmec
      REAL (dp), DIMENSION(:,:), INTENT(in)                 :: lmn_vmec
      REAL (dp), DIMENSION(:,:), INTENT(in)                 :: currumn_vmec
      REAL (dp), DIMENSION(:,:), INTENT(in)                 :: currvmn_vmec
      REAL (dp), ALLOCATABLE, DIMENSION(:,:), INTENT(inout) :: rmn_spline
      REAL (dp), ALLOCATABLE, DIMENSION(:,:), INTENT(inout) :: zmn_spline
      REAL (dp), ALLOCATABLE, DIMENSION(:,:), INTENT(inout) :: lmn_spline
      REAL (dp), ALLOCATABLE, DIMENSION(:,:), INTENT(inout) :: currumn_spline
      REAL (dp), ALLOCATABLE, DIMENSION(:,:), INTENT(inout) :: currvmn_spline
      INTEGER, INTENT(in)                                   :: parity

!  Local variables
      INTEGER                                               :: istat
      INTEGER                                               :: modes
      INTEGER                                               :: modes_nyq
      INTEGER                                               :: n_off
      INTEGER                                               :: n_max
      REAL (dp), DIMENSION(:,:), ALLOCATABLE                :: tempmn1
      REAL (dp), DIMENSION(:,:), ALLOCATABLE                :: tempmn2

!  Start of executable code
      IF (parity .eq. f_sin .and. .not.lasym_vmec) THEN
         rmn_spline = 0
         zmn_spline = 0
         lmn_spline = 0
         currumn_spline = 0
         currvmn_spline = 0
         RETURN
      END IF

      istat = 0

      CALL vmec_info_spline_modes(rmn_vmec, rmn_spline, xm_vmec, istat)
      CALL assert_eq(0, istat, 'Error splining rmn')

      CALL vmec_info_spline_modes(zmn_vmec, zmn_spline, xm_vmec, istat)
      CALL assert_eq(0, istat, 'Error splining zmn')

      ALLOCATE(tempmn1(mnmax,ns_vmec + 1), stat=istat)
      CALL vmec_info_set_ghost_points(lmn_vmec, tempmn1)
      CALL vmec_info_spline_modes(tempmn1, lmn_spline, xm_vmec, istat)
      CALL assert_eq(0, istat, 'Error splining lmn')
      DEALLOCATE(tempmn1)

!  The currents use the _nyq spectrum so we need to truncae down to just the
!  non-nyquest modes.
      ALLOCATE (tempmn1(mnmax,ns_vmec),                                        &
                tempmn2(mnmax,ns_vmec), stat=istat)
      tempmn1 = 0
      tempmn2 = 0

!  Assuming here that _nyq modes will be larger than the regular modes.
      CALL assert(mnmax_nyq .ge. mnmax, 'ERROR: More modes than _nyq modes')

      n_max = INT(MAXVAL(xn_vmec))
      n_off = 2*INT((MAXVAL(xn_nyq) - MAXVAL(xn_vmec))/nfp_vmec) + 1

      modes_nyq = 1
      DO modes = 1, mnmax
         CALL assert(xm_vmec(modes) .eq. xm_nyq(modes_nyq), 'm mode mismatch')
         CALL assert(xn_vmec(modes) .eq. xn_nyq(modes_nyq), 'n mode mismatch')

         tempmn1(modes,:) = mu0*currumn_vmec(modes_nyq,:)
         tempmn2(modes,:) = mu0*currvmn_vmec(modes_nyq,:)

         IF (INT(xn_vmec(modes)) .eq. n_max) THEN
            modes_nyq = modes_nyq + n_off
         ELSE
            modes_nyq = modes_nyq + 1
         END IF
      END DO
      IF (parity .eq. f_cos) THEN
         jsupvdotA = mu0*(SUM(tempmn2(1,2:ns_vmec - 1))                        &
                   + (tempmn2(1, 1) + tempmn2(1, ns_vmec))/2)
         jsupvdotA = jsupvdotA/(ns_vmec - 1)
      END IF

      CALL vmec_info_spline_modes(tempmn1, currumn_spline, xm_vmec, istat)
      CALL vmec_info_spline_modes(tempmn2, currvmn_spline, xm_vmec, istat)

      DEALLOCATE(tempmn1, tempmn2)

      END SUBROUTINE

!-------------------------------------------------------------------------------
!>  @brief Add ghost points to js=1 (s=0) and js = ns + 1 (s=1).
!>
!>  Extends the VMEC lmn to the ns+1 grid.
!>
!>  @param[in]    VMEC Fourier quantity.
!>  @param[inout] Ghost Fourier quantity.
!-------------------------------------------------------------------------------
     SUBROUTINE vmec_info_set_ghost_points(ymn_vmec, ymn_ghost)

     IMPLICIT NONE

!  Declare Arguments
     REAL (dp), DIMENSION(:,:), INTENT(in)  :: ymn_vmec
     REAL (dp), DIMENSION(:,:), INTENT(out) :: ymn_ghost

!  Start of executable code
     ymn_ghost(:,2:ns_vmec) = ymn_vmec(:,2:ns_vmec)
     ymn_ghost(:,ns_vmec + 1) = 2*ymn_vmec(:,ns_vmec)                          &
    &                         - ymn_vmec(:,ns_vmec - 1)

     WHERE (xm_vmec .eq. 0)
        ymn_ghost(:,1) = 2*ymn_vmec(:,2) - ymn_vmec(:,3)
     ELSE WHERE
        ymn_ghost(:,1) = 0
     END WHERE

     END SUBROUTINE

!-------------------------------------------------------------------------------
!>  @brief Load the R, Z and lambda arrays from VMEC.
!>
!>  R, Z and Lambda are recast onto the SIESTA grid. Repack splined arrays from
!>  vmec-ordering to siesta-ordering and swap n_vmec -> -n_siesta modes to be
!>  consistent with mu + nv siesta representation.
!-------------------------------------------------------------------------------
      SUBROUTINE vmec_info_set_RZL()
      USE shared_data, ONLY: lasym
      USE island_params, ONLY: mpol=>mpol_i, ntor=>ntor_i, ns=>ns_i

      IMPLICIT NONE

!  Start executable code.
      CALL vmec_info_construct_island(mpol, ntor, ns, lasym)
      CALL vmec_info_repack(rmnc_i, zmns_i, lmns_i, jcurrumnc, jcurrvmnc,      &
                            rmnc_spline, zmns_spline, lmns_spline,             &
                            currumnc_spline, currvmnc_spline, f_cos)
      IF (lasym) THEN
         CALL vmec_info_repack(rmns_i, zmnc_i, lmnc_i, jcurrumns, jcurrvmns,   &
                               rmns_spline, zmnc_spline, lmnc_spline,          &
                               currumns_spline, currvmns_spline, f_sin)
      END IF
      CALL vmec_info_destruct_spline

      END SUBROUTINE

!*******************************************************************************
!  Utilitiy SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  @brief Resets VMEC Fourier quanties to use SIESTA convension.
!>
!>  The splined arrays (rmnc_spl, etc) are all in VMEC ordering.
!>
!>    rmnc_spl(1:mnmax,ns_i)
!>
!>  Now pack them so they can be used by island Fourier routines.
!>
!>    rmn(0:mpol, -ntor:ntor, ns_i)
!>
!>  Arrays passed into the spline arguments get deallocated.
!>
!>  @param[inout] rmn         Repacked R Fourier coefficents.
!>  @param[inout] zmn         Repacked Z Fourier coefficents.
!>  @param[inout] lmn         Repacked Lambda Fourier coefficents.
!>  @param[inout] currumn     Repacked jsupu Fourier coefficents.
!>  @param[inout] currvmn     Repacked jsupv Fourier coefficents.
!>  @param[in]    rmn_spl     Splined R Fourier coefficents.
!>  @param[in]    zmn_spl     Splined Z Fourier coefficents.
!>  @param[in]    lmn_spl     Splined Lambda Fourier coefficents.
!>  @param[in]    currumn_spl Splined jsupu Fourier coefficents.
!>  @param[in]    currvmn_spl Splined jsupv Fourier coefficents.
!>  @param[in]    parity      True if using stellarator symmetry.
!>
!>  @note
!>  mpol_i == mpol_vmec, ntor_i == ntor_vmec here. Also, must flip the sign of
!>  'n' terms since VMEC uses mu-nv and SIESTA uses mu+nv for args of cos and
!>  sin.
!>
!>  Keep lmn on half mesh and add extrapolation points at s=1 (rho=0)(ns+1 -> 1)
!>
!>  Rountines here assume that *mn have 0:mpol_i and -ntor_i:ntor_i array
!>  bounds.
!-------------------------------------------------------------------------------
      SUBROUTINE vmec_info_repack(rmn, zmn, lmn, currumn, currvmn,             &
                                  rmn_spl, zmn_spl, lmn_spl,                   &
                                  currumn_spl, currvmn_spl, parity)
      USE fourier, ONLY: m1
      USE island_params, ONLY: mpol=>mpol_i, ntor=>ntor_i, ns=>ns_i,           &
     &                         nfp=>nfp_i, hs=>hs_i, fourier_context
      USE descriptor_mod, ONLY: iam
      USE shared_data, ONLY: lverbose
      USE utilities

      IMPLICIT NONE

!  Declare Arguments
      REAL (dp), ALLOCATABLE, DIMENSION(:,:,:), INTENT(inout) :: rmn
      REAL (dp), ALLOCATABLE, DIMENSION(:,:,:), INTENT(inout) :: zmn
      REAL (dp), ALLOCATABLE, DIMENSION(:,:,:), INTENT(inout) :: lmn
      REAL (dp), ALLOCATABLE, DIMENSION(:,:,:), INTENT(inout) :: currumn
      REAL (dp), ALLOCATABLE, DIMENSION(:,:,:), INTENT(inout) :: currvmn
      REAL (dp), ALLOCATABLE, DIMENSION(:,:), INTENT(in)      :: rmn_spl
      REAL (dp), ALLOCATABLE, DIMENSION(:,:), INTENT(in)      :: zmn_spl
      REAL (dp), ALLOCATABLE, DIMENSION(:,:), INTENT(in)      :: lmn_spl
      REAL (dp), ALLOCATABLE, DIMENSION(:,:), INTENT(in)      :: currumn_spl
      REAL (dp), ALLOCATABLE, DIMENSION(:,:), INTENT(in)      :: currvmn_spl
      INTEGER, INTENT(in)                                     :: parity

!  Local variables
      INTEGER                                                 :: mn_mode
      INTEGER                                                 :: js
      INTEGER                                                 :: m
      INTEGER                                                 :: n
      INTEGER                                                 :: n_abs
      INTEGER                                                 :: s1
      INTEGER                                                 :: s2
      INTEGER                                                 :: istat

!  Start executable code.
      CALL assert_eq(mnmax, SIZE(rmn_spl,1),'rmn_spl wrong size1 in REPACK')
      CALL assert_eq(ns, SIZE(rmn_spl,2),'rmn_spl wrong size2 in REPACK')
      CALL assert_eq(ns + 1, SIZE(lmn_spl,2),'lmn_spl wrong size2 in REPACK')

      IF (parity .eq. f_cos) THEN
         s1 = 1
      ELSE
         s2 = 1
      END IF

      IF (iam .eq. 0 .and. lverbose .and. parity .eq. f_cos) THEN
         m = INT(MAXVAL(xm_vmec))
         n = INT(MAXVAL(ABS(xn_vmec)))/nfp
         IF (m .gt. mpol .or. .not.fourier_context%get_index(n)) THEN
            WRITE(*,1000) m, n
         END IF
         CALL assert_eq(0, MOD(nfp_vmec, nfp), 'nfpin should be an even ' //   &
                                               'divisor of the VMEC nfp.')
      END IF

!  Flip vmec n to siesta -n so trig arg is in siesta notation. Note that the
!  siesta n modes do not need to match the vmec n modes.
      DO mn_mode = 1,mnmax
         m = xm_vmec(mn_mode)
         n = xn_vmec(mn_mode)/nfp
         IF (.not.fourier_context%get_index(n)) THEN
            CYCLE
         END IF
         n_abs = ABS(n)
         IF (m .le. mpol .and. n_abs .le. ntor) THEN

!  For m=0 load only the positive n values.
            IF (m .eq. 0) THEN
               IF (parity .eq. f_cos) THEN
                  s2 = -SIGN(1,n)
               END IF
               IF (parity .eq. f_sin) THEN
                  s1 = -SIGN(1,n)
               END IF

               rmn(m,n_abs,:) = s1*rmn_spl(mn_mode,:)
               zmn(m,n_abs,:) = s2*zmn_spl(mn_mode,:)
               lmn(m,n_abs,:) = s2*lmn_spl(mn_mode,:)*iflipj
               currumn(m,n_abs,:) = s1*currumn_spl(mn_mode,:)
               currvmn(m,n_abs,:) = s1*currvmn_spl(mn_mode,:)
            ELSE
               rmn(m,-n*iflipj,:) = rmn_spl(mn_mode,:)
               zmn(m,-n*iflipj,:) = zmn_spl(mn_mode,:)*iflipj
               lmn(m,-n*iflipj,:) = lmn_spl(mn_mode,:)
               currumn(m,-n*iflipj,:) = currumn_spl(mn_mode,:)
               currvmn(m,-n*iflipj,:) = currvmn_spl(mn_mode,:)
            END IF
         END IF
      END DO

!  Divide out the orthonorm factor used by the Fourier routines. mu0 is already
!  multipled into the current in vmec_info_set_RZL_splines.
      DO js = 1, ns
         lmn(:,:,js) = lmn(:,:,js)/fourier_context%orthonorm(0:mpol,-ntor:ntor)
         rmn(:,:,js) = rmn(:,:,js)/fourier_context%orthonorm(0:mpol,-ntor:ntor)
         zmn(:,:,js) = zmn(:,:,js)/fourier_context%orthonorm(0:mpol,-ntor:ntor)
         currumn(:,:,js) = 2*hs*(js - 1)*currumn(:,:,js)                       &
                         / fourier_context%orthonorm(0:mpol,-ntor:ntor)
         currvmn(:,:,js) = 2*hs*(js - 1)*currvmn(:,:,js)                       &
                         / fourier_context%orthonorm(0:mpol,-ntor:ntor)
      END DO

      lmn(:,:,ns + 1) = lmn(:,:,ns + 1)                                        &
                      / fourier_context%orthonorm(0:mpol,-ntor:ntor)

      rmn(m1:,:,1) = 0
      zmn(m1:,:,1) = 0

      IF (parity .eq. f_cos) THEN
         CALL set_bndy_fouier_m0(lmn, f_sin)
         CALL set_bndy_fouier_m0(rmn, f_cos)
         CALL set_bndy_fouier_m0(zmn, f_sin)
         CALL set_bndy_fouier_m0(currumn, f_cos)
         CALL set_bndy_fouier_m0(currvmn, f_cos)
      ELSE
         CALL set_bndy_fouier_m0(lmn, f_cos)
         CALL set_bndy_fouier_m0(rmn, f_sin)
         CALL set_bndy_fouier_m0(zmn, f_cos)
         CALL set_bndy_fouier_m0(currumn, f_sin)
         CALL set_bndy_fouier_m0(currvmn, f_sin)
      END IF

      currumn(:,:,ns) = 0
      currvmn(:,:,ns) = 0

1000  FORMAT(' It is recommended to increase the number of modes to at' &
             ' least VMEC values mpol=',i4,' ntor=',i4)

      END SUBROUTINE

!-------------------------------------------------------------------------------
!>  @brief Spline reinterpolate a VMEC Fourier quantity.
!>
!>  Redistributes the VMEC quantity from the VMEC radial grid to he SIESTA
!>  radial grid.
!>
!>  @param[in]    y_vmec   Orginal VMEC Fourier quantity.
!>  @param[out]   y_spline Splined VMEC Fourier quantity.
!>  @param[in]    xm_in    Poloidal mode numbers.
!>  @param[inout] istat    Error status.
!>
!>  @note
!>  When called using the lambda quantities, this is one element larger in the
!>  radial.
!-------------------------------------------------------------------------------
      SUBROUTINE vmec_info_spline_modes(ymn_vmec, ymn_spline, xm_in, istat)
      USE island_params, hs=>hs_i, ohs=>ohs_i, ns=>ns_i
      USE descriptor_mod, ONLY: iam

      IMPLICIT NONE

!  Declare Arguments
      REAL (dp), DIMENSION(:,:), INTENT(in)  :: ymn_vmec
      REAL (dp), DIMENSION(:,:), INTENT(out) :: ymn_spline
      REAL (dp), DIMENSION(:), INTENT(in)    :: xm_in
      INTEGER, INTENT(inout)                 :: istat

!  Local variables
      INTEGER                                :: js
      INTEGER                                :: modes
      INTEGER                                :: mp
      REAL (dp), ALLOCATABLE, DIMENSION(:)   :: snodes_vmec
      REAL (dp), ALLOCATABLE, DIMENSION(:)   :: y2_vmec
      REAL (dp), ALLOCATABLE, DIMENSION(:)   :: fac1
      REAL (dp), ALLOCATABLE, DIMENSION(:)   :: y_vmec
      REAL (dp), ALLOCATABLE, DIMENSION(:)   :: snodes
      REAL (dp), ALLOCATABLE, DIMENSION(:)   :: fac2
      REAL (dp), ALLOCATABLE, DIMENSION(:)   :: y_spline
      REAL (dp)                              :: expm
      REAL (dp)                              :: offset

!  Local parameters
      REAL (dp), PARAMETER                   :: one = 1

!  Start of executable code

!  Check array sizes.
      IF (SIZE(ymn_vmec, 2) .le. 1) THEN
         istat = 1
         RETURN
      END IF

      IF (ns .LE. 1) THEN
         istat = 2
         RETURN
      END IF

      ALLOCATE(snodes_vmec(SIZE(ymn_vmec, 2)))
      ALLOCATE(y2_vmec(SIZE(ymn_vmec, 2)))
      ALLOCATE(fac1(SIZE(ymn_vmec, 2)))
      ALLOCATE(y_vmec(SIZE(ymn_vmec, 2)))

      ALLOCATE(snodes(SIZE(ymn_spline, 2)))
      ALLOCATE(fac2(SIZE(ymn_spline, 2)))
      ALLOCATE(y_spline(SIZE(ymn_spline, 2)))

      IF (SIZE(snodes_vmec) .eq. ns_vmec) THEN
         offset = 1.0_dp
      ELSE
         offset = 1.5_dp
         snodes_vmec(SIZE(snodes_vmec)) = (SIZE(snodes_vmec) - offset)         &
                                        / (ns_vmec - 1.0_dp)
         fac1(SIZE(snodes_vmec)) = 1;
         fac2(SIZE(snodes)) = 1
         snodes(SIZE(snodes)) = (SIZE(snodes) - offset)/(ns_vmec - 1.0_dp)
      END IF

!  Set up "knots" on initial (vmec, svmec ~ phi) mesh and factor for taking out
!  (putting back) [sqrt(s)]**m factor for odd m.
      snodes_vmec(1) = 0
      fac1(1) = 1
      DO js = 2, ns_vmec
         snodes_vmec(js) = (js - offset)/(ns_vmec - 1.0_dp)
      END DO
!  Regularization factor: sqrt(FLUX) for odd modes
      fac1(2:ns_vmec) = one/SQRT(snodes_vmec(2:ns_vmec))

!  Set up s-nodes on final (snodes, splined) mesh [s_siesta(sfinal) ~ sfinal**2]
!  if s_siesta ~ sfinal, this is then the original vmec mesh.
      ohs = ns - 1
      hs = one/ohs
      snodes(1) = 0
      DO js = 2, ns
         snodes(js) = hs*(js - offset)
      END DO

      snodes = snodes*snodes !SIESTA s==fac2 ~ sqrt(s-vmec) mesh
      fac2 = SQRT(snodes)

      DO modes = 1, SIZE(xm_in)
         expm = 0  ! MRC: Check this. Orginal code looked like this would always
                   !      be non zero if the first loop iteration reset it.

         y_vmec = ymn_vmec(modes,:)
         mp = xm_in(modes)

         IF (istat .eq. 0 .and. mp .gt. 0) THEN
            IF (MOD(mp,2) .eq. 1) THEN
               expm = 1
            ELSE
               expm = 2
            END IF
            y_vmec = y_vmec*(fac1**expm)
            IF (mp .le. 2) THEN
               y_vmec(1) = 2*y_vmec(2) - y_vmec(3)
            END IF
         END IF

!  Initialize spline for each mode amplitude (factor out sqrt(s) factor for
!  odd-m) (spline routines in LIBSTELL Miscel folder).
         CALL spline(snodes_vmec, y_vmec, SIZE(snodes_vmec),                   &
                     lower_b, upper_b, y2_vmec)

!  Interpolate onto snodes mesh.
         DO js = 1, SIZE(snodes)
            CALL splint(snodes_vmec, y_vmec, y2_vmec,                          &
                        SIZE(snodes_vmec), snodes(js), y_spline(js))
         END DO

         IF (istat .eq. 0 .and. mp .gt. 0 .and. expm .ne. 0) THEN
            y_spline = y_spline*(fac2**expm)
         END IF
         ymn_spline(modes,:) = y_spline(:)
      END DO

      istat = 0

      DEALLOCATE(snodes_vmec)
      DEALLOCATE(y2_vmec)
      DEALLOCATE(fac1)
      DEALLOCATE(y_vmec)

      DEALLOCATE(snodes)
      DEALLOCATE(fac2)
      DEALLOCATE(y_spline)

      END SUBROUTINE

!-------------------------------------------------------------------------------
!>  @brief Spline reinterpolate a VMEC quantity.
!>
!>  Redistributes the VMEC quantity from the VMEC radial grid to he SIESTA
!>  radial grid.
!>
!>  @param[in]  y_vmec   Orginal VMEC quantity.
!>  @param[out] y_spline Splined VMEC quantity.
!>  @param[out] istat    Error status.
!-------------------------------------------------------------------------------
      SUBROUTINE vmec_info_spline_oned_array(y_vmec, y_spline, istat)
      USE island_params, ONLY: hs=>hs_i, ns=>ns_i

      IMPLICIT NONE

!  Declare Arguments
      REAL (dp), DIMENSION(ns_vmec), INTENT(in) :: y_vmec
      REAL (dp), DIMENSION(ns), INTENT(out)     :: y_spline
      INTEGER, INTENT(OUT)                      :: istat

!  Local variables
      INTEGER                                   :: js
      REAL (dp), DIMENSION(ns_vmec)             :: snodes_vmec
      REAL (dp), DIMENSION(ns_vmec)             :: y2_vmec
      REAL (dp), DIMENSION(ns)                  :: snodes
      REAL (dp)                                 :: hs_vmec

!  Local parameters
      REAL (dp), PARAMETER                      :: one = 1

!  Start of executable code

!  Check grid dimensions.
      IF (ns_vmec .le. 1) THEN
         istat = 1
         RETURN
      END IF

      IF (ns .LE. 1) THEN
         istat = 2
         RETURN
      END IF

!  Set up "knots" on initial (vmec, svmec ~ phi) mesh.
      hs_vmec = one/(ns_vmec - 1)
      DO js = 1, ns_vmec
         snodes_vmec(js) = hs_vmec*(js - 1)
      END DO

!  Set up s-nodes on final (snodes, splined) mesh [s_siesta(sfinal) ~ sfinal**2]
!  if s_siesta ~ sfinal, this is then the original vmec mesh.
!      ohs = ns - 1
!      hs = one/ohs
      DO js = 1, ns
         snodes(js) = hs*(js - 1) !vmec mesh   s-siesta~s-vmec
         snodes(js) = snodes(js)*snodes(js) !polar mesh, s-siesta~s-vmec**2
      END DO

!  Initialize spline for each mode amplitude (factor out sqrt(s) factor for
!  odd-m)
      CALL spline(snodes_vmec, y_vmec, ns_vmec, lower_b, upper_b, y2_vmec)

!  Interpolate onto snodes mesh
      DO js = 1, ns
         CALL splint(snodes_vmec, y_vmec, y2_vmec, ns_vmec,                    &
                     snodes(js), y_spline(js))
      END DO

      END SUBROUTINE

      END MODULE
