!*******************************************************************************
!>  @file utilities.f90
!>  @brief Contains module @ref utilities.
!
!  Note separating the Doxygen comment block here so detailed decription is
!  found in the Module not the file.
!
!>  This file contains subroutines for allocating and initializing curvilinear
!>  magnetic covariant and pressure Fourier harmonics on the half-radial mesh.
!*******************************************************************************
      MODULE quantities
      USE stel_kinds
      USE stel_constants      
      USE v3_utilities, ONLY: assert, assert_eq
      USE island_params, ns=>ns_i, ntheta=>nu_i, nzeta=>nv_i,           &
          mpol=>mpol_i, ntor=>ntor_i, nuv=>nuv_i, mnmax=>mnmax_i,       &
          nfp=>nfp_i
      USE descriptor_mod, ONLY: iam, nprocs, SIESTA_COMM
      USE shared_data, ONLY: lasym, lverbose
      USE nscalingtools, ONLY: startglobrow, endglobrow, PARSOLVER, MPI_ERR
      USE utilities, ONLY: to_half_mesh, to_full_mesh
      USE siesta_namelist, ONLY: l_vessel, l_lambda

      IMPLICIT NONE

!*******************************************************************************
!  Module variables
!*******************************************************************************
!>  Magnetic scaling factor to scale internal energy to 1.
      REAL (dp) :: b_factor
!>  Pressure scaling factor to scale internal energy to 1.
      REAL (dp) :: p_factor
!>  Sign of the jacobian.
      REAL (dp) :: signjac
!>  Energy due to the pressure.
      REAL (dp) :: wp
!>  Energy due to the magnetic field.
      REAL (dp) :: wb

!>  Fouier amplitudes of the contravariant magnetic field and jacobian for
!>  stellarator symmetric parity in the s direction.
      REAL (dp), ALLOCATABLE, DIMENSION(:,:,:) :: jbsupsmnsh
!>  Fouier amplitudes of the contravariant magnetic field and jacobian for
!>  stellarator symmetric parity in the u direction.
      REAL (dp), ALLOCATABLE, DIMENSION(:,:,:) :: jbsupumnch
!>  Fouier amplitudes of the contravariant magnetic field and jacobian for
!>  stellarator symmetric parity in the v direction.
      REAL (dp), ALLOCATABLE, DIMENSION(:,:,:) :: jbsupvmnch
!>  Fouier amplitudes of the pressure and jacobian for stellarator symmetric
!>  parity.
      REAL (dp), ALLOCATABLE, DIMENSION(:,:,:) :: jpmnch
!>  Fouier amplitudes of the contravariant current for stellarator symmetric
!>  parity in the s direction.
      REAL (dp), ALLOCATABLE, DIMENSION(:,:,:) :: ksupsmnsf
!>  Fouier amplitudes of the contravariant current for stellarator symmetric
!>  parity in the u direction.
      REAL (dp), ALLOCATABLE, DIMENSION(:,:,:) :: ksupumncf
!>  Fouier amplitudes of the contravariant current for stellarator symmetric
!>  parity in the v direction.
      REAL (dp), ALLOCATABLE, DIMENSION(:,:,:) :: ksupvmncf
!>  Fouier amplitudes of the pressure perturbation and jacobian for stellarator
!>  symmetric parity.
      REAL (dp), ALLOCATABLE, DIMENSION(:,:,:) :: djpmnch
!>  Fouier amplitudes of the contravariant magnetic field perturbation and
!>  jacobian for stellarator symmetric parity in the s direction.
      REAL (dp), ALLOCATABLE, DIMENSION(:,:,:) :: djbsupsmnsh
!>  Fouier amplitudes of the contravariant magnetic field perturbation and
!>  jacobian for stellarator symmetric parity in the u direction.
      REAL (dp), ALLOCATABLE, DIMENSION(:,:,:) :: djbsupumnch
!>  Fouier amplitudes of the contravariant magnetic field perturbation and
!>  jacobian for stellarator symmetric parity in the v direction.
      REAL (dp), ALLOCATABLE, DIMENSION(:,:,:) :: djbsupvmnch

!>  Fouier amplitudes of the contravariant displacement vector and jacobian for
!>  stellarator symmetric parity in the s direction.
      REAL (dp), POINTER, DIMENSION(:,:,:) :: jvsupsmncf
!>  Fouier amplitudes of the contravariant displacement vector and jacobian for
!>  stellarator symmetric parity in the u direction.
      REAL (dp), POINTER, DIMENSION(:,:,:) :: jvsupumnsf
!>  Fouier amplitudes of the contravariant displacement vector and jacobian for
!>  stellarator symmetric parity in the v direction.
      REAL (dp), POINTER, DIMENSION(:,:,:) :: jvsupvmnsf
!>  Fouier amplitudes of the covariant force for stellarator symmetric parity
!>  in the s direction.
      REAL (dp), POINTER, DIMENSION(:,:,:) :: fsubsmncf
!>  Fouier amplitudes of the covariant force for stellarator symmetric parity
!>  in the u direction.
      REAL (dp), POINTER, DIMENSION(:,:,:) :: fsubumnsf
!>  Fouier amplitudes of the covariant force for stellarator symmetric parity
!>  in the v direction.
      REAL (dp), POINTER, DIMENSION(:,:,:) :: fsubvmnsf
!>  Fouier amplitudes of the contravariant force for stellarator symmetric
!>  parity in the s direction.
      REAL (dp), POINTER, DIMENSION(:,:,:) :: fsupsmncf
!>  Fouier amplitudes of the contravariant force for stellarator symmetric
!>  parity in the u direction.
      REAL (dp), POINTER, DIMENSION(:,:,:) :: fsupumnsf
!>  Fouier amplitudes of the contravariant force for stellarator symmetric
!>  parity in the v direction.
      REAL (dp), POINTER, DIMENSION(:,:,:) :: fsupvmnsf

!>  Fouier amplitudes of the contravariant magnetic field and jacobian for
!>  stellarator asymmetric parity in the s direction.
      REAL (dp), ALLOCATABLE, DIMENSION(:,:,:) :: jbsupsmnch
!>  Fouier amplitudes of the contravariant magnetic field and jacobian for
!>  stellarator asymmetric parity in the u direction.
      REAL (dp), ALLOCATABLE, DIMENSION(:,:,:) :: jbsupumnsh
!>  Fouier amplitudes of the contravariant magnetic field and jacobian for
!>  stellarator asymmetric parity in the v direction.
      REAL (dp), ALLOCATABLE, DIMENSION(:,:,:) :: jbsupvmnsh
!>  Fouier amplitudes of the pressure and jacobian for stellarator asymmetric
!>  parity.
      REAL (dp), ALLOCATABLE, DIMENSION(:,:,:) :: jpmnsh
!>  Fouier amplitudes of the contravariant current for stellarator asymmetric
!>  parity in the s direction.
      REAL (dp), ALLOCATABLE, DIMENSION(:,:,:) :: ksupsmncf
!>  Fouier amplitudes of the contravariant current for stellarator asymmetric
!>  parity in the u direction.
      REAL (dp), ALLOCATABLE, DIMENSION(:,:,:) :: ksupumnsf
!>  Fouier amplitudes of the contravariant current for stellarator asymmetric
!>  parity in the v direction.
      REAL (dp), ALLOCATABLE, DIMENSION(:,:,:) :: ksupvmnsf
!>  Fouier amplitudes of the pressure perturbation and jacobian for stellarator
!>  asymmetric parity.
      REAL (dp), ALLOCATABLE, DIMENSION(:,:,:) :: djpmnsh
!>  Fouier amplitudes of the contravariant magnetic field perturbation and
!>  jacobian for astellarator symmetric parity in the s direction.
      REAL (dp), ALLOCATABLE, DIMENSION(:,:,:) :: djbsupsmnch
!>  Fouier amplitudes of the contravariant magnetic field perturbation and
!>  jacobian for astellarator symmetric parity in the u direction.
      REAL (dp), ALLOCATABLE, DIMENSION(:,:,:) :: djbsupumnsh
!>  Fouier amplitudes of the contravariant magnetic field perturbation and
!>  jacobian for astellarator symmetric parity in the v direction.
      REAL (dp), ALLOCATABLE, DIMENSION(:,:,:) :: djbsupvmnsh

!>  Fouier amplitudes of the contravariant displacement vector and jacobian for
!>  stellarator asymmetric parity in the s direction.
      REAL (dp), POINTER, DIMENSION(:,:,:)     :: jvsupsmnsf
!>  Fouier amplitudes of the contravariant displacement vector and jacobian for
!>  stellarator asymmetric parity in the u direction.
      REAL (dp), POINTER, DIMENSION(:,:,:)     :: jvsupumncf
!>  Fouier amplitudes of the contravariant displacement vector and jacobian for
!>  stellarator asymmetric parity in the v direction.
      REAL (dp), POINTER, DIMENSION(:,:,:)     :: jvsupvmncf
!>  Fouier amplitudes of the covariant force for stellarator asymmetric parity
!>  in the s direction.
      REAL (dp), POINTER, DIMENSION(:,:,:)     :: fsubsmnsf
!>  Fouier amplitudes of the covariant force for stellarator asymmetric parity
!>  in the u direction.
      REAL (dp), POINTER, DIMENSION(:,:,:)     :: fsubumncf
!>  Fouier amplitudes of the covariant force for stellarator asymmetric parity
!>  in the v direction.
      REAL (dp), POINTER, DIMENSION(:,:,:)     :: fsubvmncf
!>  Fouier amplitudes of the contravariant force for stellarator asymmetric
!>  parity in the s direction.
      REAL (dp), POINTER, DIMENSION(:,:,:)     :: fsupsmnsf
!>  Fouier amplitudes of the contravariant force for stellarator asymmetric
!>  parity in the u direction.
      REAL (dp), POINTER, DIMENSION(:,:,:)     :: fsupumncf
!>  Fouier amplitudes of the contravariant force for stellarator asymmetric
!>  parity in the v direction.
      REAL (dp), POINTER, DIMENSION(:,:,:)     :: fsupvmncf

!>  Initial power spectrum storage for stellarator symmetric parity.
      REAL(dp), ALLOCATABLE, DIMENSION(:,:,:,:) :: pwr_spec_s
!>  Initial power spectrum storage for stellarator asymmetric parity.
      REAL(dp), ALLOCATABLE, DIMENSION(:,:,:,:) :: pwr_spec_a

!>  Real space contravariant displacement vector in the s direction.
      REAL (dp), ALLOCATABLE, DIMENSION(:,:,:) :: jvsupsijf
!>  Real space contravariant displacement vector in the u direction.
      REAL (dp), ALLOCATABLE, DIMENSION(:,:,:) :: jvsupuijf
!>  Real space contravariant displacement vector in the v direction.
      REAL (dp), ALLOCATABLE, DIMENSION(:,:,:) :: jvsupvijf
!>  Real space jacobian on the half grid.
      REAL (dp), ALLOCATABLE, DIMENSION(:,:,:) :: jacobh
!>  Real space jacobian on the full grid.
      REAL (dp), ALLOCATABLE, DIMENSION(:,:,:) :: jacobf
!>  Volumn integration element.
      REAL (dp), ALLOCATABLE, DIMENSION(:,:,:) :: wint

!>  Unperturbed realspace contravariant magnetic field in the s direction.
      REAL (dp), ALLOCATABLE, DIMENSION(:,:,:) :: bsupsijf0
!>  Unperturbed realspace contravariant magnetic field in the v direction.
      REAL (dp), ALLOCATABLE, DIMENSION(:,:,:) :: bsupuijf0
!>  Unperturbed realspace contravariant magnetic field in the u direction.
      REAL (dp), ALLOCATABLE, DIMENSION(:,:,:) :: bsupvijf0
!>  Unperturbed half grid realspace contravariant magnetic field in the s
!>  direction.
      REAL (dp), ALLOCATABLE, DIMENSION(:,:,:) :: bsupsijh0
!>  Unperturbed half grid realspace contravariant magnetic field in the u
!>  direction.
      REAL (dp), ALLOCATABLE, DIMENSION(:,:,:) :: bsupuijh0
!>  Unperturbed half grid realspace contravariant magnetic field in the v
!>  direction.
      REAL (dp), ALLOCATABLE, DIMENSION(:,:,:) :: bsupvijh0
!>  Perturbed full grid realspace contravariant magnetic field in the s
!>  direction.
      REAL (dp), ALLOCATABLE, DIMENSION(:,:,:) :: bsupsijf
!>  Perturbed full grid realspace contravariant magnetic field in the u
!>  direction.
      REAL (dp), ALLOCATABLE, DIMENSION(:,:,:) :: bsupuijf
!>  Perturbed full grid realspace contravariant magnetic field in the v
!>  direction.
      REAL (dp), ALLOCATABLE, DIMENSION(:,:,:) :: bsupvijf
!>  Perturbed full grid realspace covariant magnetic field in the s direction.
      REAL (dp), ALLOCATABLE, DIMENSION(:,:,:) :: bsubsijf
!>  Perturbed full grid realspace covariant magnetic field in the u direction.
      REAL (dp), ALLOCATABLE, DIMENSION(:,:,:) :: bsubuijf
!>  Perturbed full grid realspace covariant magnetic field in the v direction.
      REAL (dp), ALLOCATABLE, DIMENSION(:,:,:) :: bsubvijf
!>  Realspace magnetude of |B|^2.
      REAL (dp), ALLOCATABLE, DIMENSION(:,:,:) :: bsq
!>  Perturbed full grid realspace covariant current in the s direction.
      REAL (dp), ALLOCATABLE, DIMENSION(:,:,:) :: ksubsijf
!>  Perturbed full grid realspace covariant current in the u direction.
      REAL (dp), ALLOCATABLE, DIMENSION(:,:,:) :: ksubuijf
!>  Perturbed full grid realspace covariant current in the v direction.
      REAL (dp), ALLOCATABLE, DIMENSION(:,:,:) :: ksubvijf
!>  Unperturbed full grid realspace contravariant current in the s direction.
      REAL (dp), ALLOCATABLE, DIMENSION(:,:,:) :: ksupsijf0
!>  Unperturbed full grid realspace contravariant current in the u direction.
      REAL (dp), ALLOCATABLE, DIMENSION(:,:,:) :: ksupuijf0
!>  Unperturbed full grid realspace contravariant current in the v direction.
      REAL (dp), ALLOCATABLE, DIMENSION(:,:,:) :: ksupvijf0
!>  Perturbed full grid realspace contravariant current in the s direction.
      REAL (dp), ALLOCATABLE, DIMENSION(:,:,:) :: ksupsijf
!>  Perturbed full grid realspace contravariant current in the u direction.
      REAL (dp), ALLOCATABLE, DIMENSION(:,:,:) :: ksupuijf
!>  Perturbed full grid realspace contravariant current in the v direction.
      REAL (dp), ALLOCATABLE, DIMENSION(:,:,:) :: ksupvijf

!>  Unperturbed half grid realspace pressure.
      REAL (dp), ALLOCATABLE, DIMENSION(:,:,:) :: pijh0
!>  Unperturbed half grid realspace u gradient of the pressure.
      REAL (dp), ALLOCATABLE, DIMENSION(:,:,:) :: pijh0_du
!>  Unperturbed half grid realspace v gradient of the pressure.
      REAL (dp), ALLOCATABLE, DIMENSION(:,:,:) :: pijh0_dv
!>  Unperturbed full grid realspace pressure.
      REAL (dp), ALLOCATABLE, DIMENSION(:,:,:) :: pijf0
!>  Unperturbed full grid realspace s gradient of the pressure.
      REAL (dp), ALLOCATABLE, DIMENSION(:,:,:) :: pijf0_ds

!>  Saved boundary harmonics.
      REAL (dp) :: fbdy(13)

      CONTAINS

!*******************************************************************************
!  UTILITY SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  @brief Intialize quantities.
!>
!>  @param[in] restarted Flag to signal that we are in restart mode.
!-------------------------------------------------------------------------------
      SUBROUTINE init_quantities(restarted)
      USE descriptor_mod, ONLY: LSCALAPACK
      USE island_params, ONLY: fourier_context
      USE fourier, ONLY: f_none, f_cos, f_sin, f_sum, m0
      USE metrics, ONLY: sqrtg
      USE shared_data, ONLY: mblk_size, ndims, lasym
#if defined(MPI_OPT)        
      USE prof_mod, ONLY: SetUpScalingAllGather
#endif
      USE blocktridiagonalsolver_s, ONLY: Initialize

      IMPLICIT NONE

!  Declare Arguments
      LOGICAL                                  :: restarted

!  Local variables.
      INTEGER                                  :: istat
      INTEGER                                  :: l
      REAL (dp)                                :: sum1
      REAL (dp), DIMENSION(:,:,:), ALLOCATABLE :: tempmn_sym
      REAL (dp), DIMENSION(:,:,:), ALLOCATABLE :: tempmn_asym

!  Start of executable code.
      IF (PARSOLVER) THEN
         CALL Initialize(ns, mblk_size)
      ELSE IF (LSCALAPACK) THEN
#if defined(MPI_OPT)
         CALL SetUpScalingAllGather(mblk_size)
#endif
      ELSE
         startglobrow = 1
         endglobrow = ns
      END IF

      IF (.not.restarted) THEN
!  Allocate and initialize dependent variables
         CALL alloc_quantities

         p_factor = 1.0_dp/ABS(wb_i)
         b_factor = SQRT(p_factor)
         gnorm_i = ABS(wb_i)/(wb_i + wp_i/(gamma - 1.0))
      END IF

!  Reshapes sqrtg in local form
      jacobh(:,:,:) = RESHAPE(sqrtg, SHAPE(jacobh))

!  Avoid divide by zero and store averaged value over theta at first half grid
!  point at origin in jacobh.
      CALL fourier_context%get_origin(jacobh, m0, lasym)
      signjac = jacobh(1,1,2)/ABS(jacobh(1,1,2))

!  Filter jacobian ro preserve the unperturbed (constant in s) part of the
!  pressure.
      ALLOCATE(tempmn_sym(0:mpol,-ntor:ntor,SIZE(jacobh,3)))
      CALL fourier_context%tomnsp(jacobh, tempmn_sym, f_cos)
      IF (lasym) THEN
         ALLOCATE(tempmn_asym(0:mpol,-ntor:ntor,SIZE(jacobh,3)))
         CALL fourier_context%tomnsp(jacobh, tempmn_asym, f_sin)
      END IF

      CALL fourier_context%toijsp(tempmn_sym, jacobh, f_none, f_cos)
      DEALLOCATE(tempmn_sym)
      IF (lasym) THEN
         CALL fourier_context%toijsp(tempmn_asym, jacobh, f_sum, f_sin)
         DEALLOCATE(tempmn_asym)
      END IF

      CALL assert(ALL(jacobh*signjac .gt. 0), 'FILTERED JACOBIAN CHANGED SIGN!')

      IF (iam .eq. 0 .and. lverbose) THEN
         sum1 = SUM((jacobh(:,:,2:) - jacobf(:,:,2:))**2 /                     &
                    (jacobh(:,:,2:) + jacobf(:,:,2:))**2)
         WRITE (*,1000) SQRT(sum1/SIZE(jacobh))
      END IF

      CALL to_full_mesh(jacobh, jacobf)
!  Must be consistent with variational form of pressure evolution equation. Do
!  not use 2pt extrapoltaion.
      jacobf(:,:,1) = jacobh(:,:,1)
      jacobf(:,:,ns)= jacobh(:,:,ns)

      DO l = 1, ntheta
         wint(l,:,:) = fourier_context%cosmui(l,m0)
      END DO

      IF (.not.ALLOCATED(vp_f)) THEN
         ALLOCATE (vp_f(ns), stat=istat)
         CALL assert_eq(0, istat, 'Allocation error in init_quantities')
      END IF
      CALL SurfAverage(vp_f, jacobf, 1, ns)

!  Compar volumes for a sanity check that the siesta grid and orginal VMEC
!  volume match.
      sum1 = 0.5*hs_i*(SUM(vp_f(2:ns)) + SUM(vp_f(1:ns - 1)))*signjac
      sum1 = twopi*twopi*sum1
      IF (ABS((sum1 - volume_i)/volume_i) .GT. 1.E-3_dp) THEN
         IF (IAM .EQ. 0 .and. lverbose) THEN
            WRITE (*,1001) volume_i, sum1
         END IF
      END IF

1000  FORMAT(' JACOBIAN SPECTRAL TRUNCATION ERROR: ',1pe8.2)
1001  FORMAT(' VMEC VOLUME: ', 1pe8.2, ' SIESTA VOLUME: ',1pe8.2)

      END SUBROUTINE

!-------------------------------------------------------------------------------
!>  @brief Initalize half-mesh covariant components of JBsup* and JP.
!>
!>  Fields are gathered to all processors at the end. Only called if lrestart is
!>  false otherwise data is read directly from the restart file.
!-------------------------------------------------------------------------------
        SUBROUTINE Init_Fields
        USE island_params, ONLY: fourier_context
        USE fourier, ONLY: f_sin, f_cos, m0, m1
        USE vmec_info, ONLY: lmns_i, lmnc_i, iflipj
        USE metrics, ONLY: sqrtg

        IMPLICIT NONE

!  Local variables.
        INTEGER                                 :: istat
        INTEGER                                 :: nsmin
        INTEGER                                 :: nsmax
        INTEGER                                 :: nloc
        INTEGER                                 :: js
        REAL(dp), DIMENSION(:), ALLOCATABLE     :: phiph
        REAL(dp), DIMENSION(:), ALLOCATABLE     :: chiph
        REAL(dp), DIMENSION(:,:,:), ALLOCATABLE :: presif
        REAL(dp), DIMENSION(:,:,:), ALLOCATABLE :: presih

!  Start of executable code.
        nsmin = MAX(1, startglobrow)
        nsmax = MIN(endglobrow, ns)

        ALLOCATE(phiph(ns), chiph(ns))

        phipf_i = signjac*iflipj*phipf_i
        chipf_i = signjac*iflipj*chipf_i

        phiph(2:ns) = (phipf_i(2:ns) + phipf_i(1:nsh))/2
        chiph(2:ns) = (chipf_i(2:ns) + chipf_i(1:nsh))/2
        phiph(1) = 0
        chiph(1) = 0

        IF (.not.l_vessel .and. l_lambda) THEN
           CALL ReCompute_Lambda(lmns_i, lmnc_i, jacobh,                       &
                                 fourier_context%orthonorm, phiph, chiph,      &
                                 nsmin, nsmax)
        END IF
         
  		CALL Init_Bfield(jbsupsmnsh, jbsupumnch, jbsupvmnch,                   &
                         lmns_i, phiph, chiph, nsmin, nsmax, f_sin)
		IF (lasym) THEN
           CALL Init_Bfield(jbsupsmnch, jbsupumnsh, jbsupvmnsh,                &
                            lmnc_i, phiph, chiph, nsmin, nsmax, f_cos)
        END IF

        DEALLOCATE (phiph, chiph, stat=istat)
        CALL ASSERT(istat.EQ.0,'Deallocate error #1 in init_fields')

!  Initialize half mesh pressure
        nloc = MAX(1, startglobrow - 1)
        ALLOCATE(presih(ntheta,nzeta,nloc:nsmax),                              &
                 presif(ntheta,nzeta,nloc:nsmax), stat=istat)
        CALL ASSERT(istat.EQ.0, 'Allocate error #1 in init_fields')

        DO js = nloc, nsmax
           presif(:,:,js) = p_factor*presf_i(js)
        END DO

        CALL to_half_mesh(presif, presih)

!  Initialize (internally) jacob*(real pressure), which is evolved
        presih(:,:,nsmin:nsmax) =                                              &
       &   presih(:,:,nsmin:nsmax)*jacobh(:,:,nsmin:nsmax)

!  Initialize internal pressure harmonics (half mesh)
        CALL fourier_context%tomnsp(presih(:,:,nsmin:nsmax),                   &
                                    jpmnch(:,:,nsmin:nsmax), f_cos)
        jpmnch(:,:,1) = 0
        jpmnch(m0,:,1) = jpmnch(m0,:,2)
        IF (lasym) THEN
           CALL fourier_context%tomnsp(presih(:,:,nsmin:nsmax),                &
                                       jpmnsh(:,:,nsmin:nsmax), f_sin)
           jpmnsh(:,:,1) = 0
           jpmnsh(m0,:,1) = jpmnsh(m0,:,2)
        END IF

        DEALLOCATE(presif, presih, stat=istat)
        CALL assert_eq(istat, 0, 'Deallocate error #2 in init_fields')

!  Gather jbsupmnX's, jpmn onto all processors. This same as in update_state but
!  for initial values.
        CALL GatherFields

        END SUBROUTINE

!-------------------------------------------------------------------------------
!>  @brief Initialized magnetic field.
!>
!>  @param[inout] jbsupsmnh Contravariant magnetic field in the s direction.
!>  @param[inout] jbsupumnh Contravariant magnetic field in the u direction.
!>  @param[inout] jbsupvmnh Contravariant magnetic field in the v direction.
!>  @param[in]    lmn       Lambda Fourier coeffients.
!>  @param[in]    phiph     Radial toroidal flux derivative.
!>  @param[in]    chiph     Radial poloidal flux derivative.
!>  @param[in]    nsmin     Minimum radial index.
!>  @param[in]    nsmax     Maximum radial index.
!>  @param[in]    parity    Fourier parity.
!-------------------------------------------------------------------------------
      SUBROUTINE Init_Bfield(jbsupsmnh, jbsupumnh, jbsupvmnh,                  &
                             lmn, phiph, chiph, nsmin, nsmax, parity)
      USE fourier, ONLY: m0, m1, n0, f_sin, f_jac
      USE siesta_namelist, ONLY: l_vessel
      USE utilities, ONLY: set_bndy_fouier_m0, set_bndy_full_origin
      USE island_params, ONLY: fourier_context

      IMPLICIT NONE

!  Declare Arguments
      REAL (dp), ALLOCATABLE, DIMENSION(:,:,:), INTENT(inout) :: jbsupsmnh
      REAL (dp), ALLOCATABLE, DIMENSION(:,:,:), INTENT(inout) :: jbsupumnh
      REAL (dp), ALLOCATABLE, DIMENSION(:,:,:), INTENT(inout) :: jbsupvmnh
      REAL (dp), ALLOCATABLE, DIMENSION(:,:,:), INTENT(in)    :: lmn
      REAL (dp), DIMENSION(:), INTENT(in)                     :: phiph
      REAL (dp), DIMENSION(:), INTENT(in)                     :: chiph
	  INTEGER, INTENT(in)                                     :: nsmin
      INTEGER, INTENT(in)                                     :: nsmax
      INTEGER, INTENT(in)                                     :: parity

!  local variables
      INTEGER                                                 :: m
      INTEGER                                                 :: n
      INTEGER                                                 :: n_mode
      INTEGER                                                 :: i
      INTEGER                                                 :: sparity
      INTEGER                                                 :: mp
      INTEGER                                                 :: np

!  Start of executable code.
      IF (parity .eq. f_sin) THEN
         sparity = 1
      ELSE
         sparity = -1
      END IF

      IF (.not.l_vessel) THEN

         DO i = nsmin, nsmax
            DO n = -ntor, ntor
               n_mode = fourier_context%tor_modes(n)
               np = n_mode*nfp*sparity
               DO m = 0, mpol
                  mp = m*sparity
                  jbsupumnh(m,n,i) = phiph(i)*(-np*lmn(m,n,i))
                  jbsupvmnh(m,n,i) = phiph(i)*( mp*lmn(m,n,i))
               END DO
            END DO

            IF (parity .eq. f_sin) THEN
               jbsupumnh(m0,n0,i) = chiph(i)                                   &
                                  / fourier_context%orthonorm(m0,n0)
               jbsupvmnh(m0,n0,i) = phiph(i)                                   &
                                  / fourier_context%orthonorm(m0,n0)
            END IF
         END DO

         jbsupsmnh = 0
      END IF
        
      jbsupsmnh(:,:,nsmin:nsmax) = b_factor*jbsupsmnh(:,:,nsmin:nsmax)
      jbsupumnh(:,:,nsmin:nsmax) = b_factor*jbsupumnh(:,:,nsmin:nsmax)
      jbsupvmnh(:,:,nsmin:nsmax) = b_factor*jbsupvmnh(:,:,nsmin:nsmax)

      CALL set_bndy_fouier_m0(jbsupsmnh, jbsupumnh, jbsupvmnh, parity)

!  Store m=0 component of jbsupvmnh (and m=0,1 of jbsupumnh) at origin of
!  half-grid js=1.
      IF (nsmin .eq. 1) THEN
         jbsupsmnh(:,:,1) = jbsupsmnh(:,:,2)
         jbsupumnh(:,:,1) = jbsupumnh(:,:,2)
         jbsupvmnh(:,:,1) = jbsupvmnh(:,:,2)

         CALL set_bndy_full_origin(jbsupsmnh, jbsupumnh, jbsupvmnh, f_jac)
      END IF

      END SUBROUTINE

!-------------------------------------------------------------------------------
!>  @brief Recompute lambda on the SIESTA mesh.
!>
!>  Solves J^s = 0 to find a new lambda consistent on the SIESTA mesh. The J^s
!>  component is equal to zero when assuming a nested flux surface solution.
!>
!>    sqrt(g)J^s = Curl(B) = dB_v/du - dB_u/dv = 0                           (1)
!>
!>  From the fouier representation, equation 1 becomes.
!>
!>    0 = m*B^mn_v - n*nfp*B^mn_u                                            (2)
!>
!>  The magnetic field can be defined in terms of lambda by
!>
!>    B_ij^u = (chi' - phi'*dlambda_ij/dv)/sqrt(g)                           (3)
!>
!>    B_ij^v = (phi' + phi'*dlambda_ij/du)/sqrt(g)                           (4)
!>
!>  Using the Fourier definition equations 3 and 4 become
!>
!>    B_ij^u = (chi' - phi'*SUM(n*nfp*lambda_mn*cos(m*u + n*nfp*v)))/sqrt(g) (5)
!>
!>    B_ij^v = (phi' + phi'*SUM(m*lambda_mn*cos(m*u + n*nfp*v))/sqrt(g)      (6)
!>
!>  Where chi' is the radial derivative of the poloidal flux and phi' is the
!>  radial derivative of the toroidal flux. The covariant components are then
!>
!>    B^ij_u = g_uu*B_ij^u + g_uv*B_ij^v                                     (7)
!>
!>    B^ij_v = g_uv*B_ij^u + g_vv*B_ij^v                                     (8)
!>
!>  B_u and B_v need to be transformed into Fourier space for use in equation
!>  2.
!>
!>  @param[out] lmns      Lambda for stellarator symmetric parity.
!>  @param[out] lmnc      Lambda for stellarator asymmetric parity.
!>  @param[in]  jacobh    Jacobian on the half mesh.
!>  @param[in]  orthonorm Fouier normalization factors.
!>  @param[in]  phiph     Radial dervative of the toroidal flux.
!>  @param[in]  chiph     Radial dervative of the poloidal flux.
!>  @param[in]  nsmin     Minimum radial index.
!>  @param[in]  nsmax     Maximum radial index.
!-------------------------------------------------------------------------------
      SUBROUTINE ReCompute_Lambda(lmns, lmnc, jacobh, orthonorm,               &
                                  phiph, chiph, nsmin, nsmax)
      USE metrics, ONLY: guu, guv, gvv

      IMPLICIT NONE

!  Declare Arguments
      REAL (dp), DIMENSION(mnmax,ns), INTENT(out) :: lmns
      REAL (dp), DIMENSION(mnmax,ns), INTENT(out) :: lmnc
      REAL (dp), DIMENSION(nuv,ns), INTENT(in)    :: jacobh
      REAL (dp), DIMENSION(mnmax), INTENT(in)     :: orthonorm
      REAL (dp), DIMENSION(ns), INTENT(in)        :: phiph
      REAL (dp), DIMENSION(ns), INTENT(in)        :: chiph
      INTEGER, INTENT(in)                         :: nsmin
      INTEGER, INTENT(in)                         :: nsmax

!  Local Variables
      INTEGER                                     :: js            , i
      INTEGER                                     :: m
      INTEGER                                     :: n
      INTEGER                                     :: n_mode
      INTEGER                                     :: mp
      INTEGER                                     :: np
      INTEGER                                     :: np_mode
      INTEGER                                     :: mn
      INTEGER                                     :: mnp
      INTEGER                                     :: lk
      INTEGER                                     :: lu
      INTEGER                                     :: lv
      INTEGER                                     :: info
      INTEGER                                     :: isign
      INTEGER                                     :: n2
      REAL (dp), DIMENSION(:,:), ALLOCATABLE      :: amat
      REAL (dp), DIMENSION(:,:), ALLOCATABLE      :: brhs
      REAL (dp), DIMENSION(:), ALLOCATABLE        :: atemp
      REAL (dp), DIMENSION(:), ALLOCATABLE        :: btemp
      REAL (dp), DIMENSION(:,:), ALLOCATABLE      :: cosmni(:,:)
      REAL (dp), DIMENSION(:,:), ALLOCATABLE      :: cosmnp(:,:)
      REAL (dp), DIMENSION(:,:), ALLOCATABLE      :: sinmni(:,:)
      REAL (dp), DIMENSION(:,:), ALLOCATABLE      :: sinmnp(:,:)
      REAL (dp)                                   :: ton
      REAL (dp)                                   :: toff

!  Start of executable code
      CALL second0(ton)

      lmns(:,1) = 0
      n2 = 1
      IF (lasym) THEN
         lmnc(:,1) = 0
         n2 = 2
      END IF

      ALLOCATE(atemp(nuv), btemp(nuv), cosmni(nuv,mnmax), cosmnp(nuv,mnmax),   &
               brhs(mnmax,n2), amat(mnmax*n2, mnmax*n2), stat=info)
      CALL ASSERT(info.EQ.0,'Allocation error in Recompute_Lambda')

      IF (lasym) THEN
         ALLOCATE (sinmni(nuv,mnmax), sinmnp(nuv,mnmax), stat=info)
         CALL ASSERT(info.EQ.0,'Allocation error in Recompute_Lambda')
      END IF

!  Compute matrix elements surface by surface.
      mn = 0
      DO n = -ntor, ntor
         isign = 1
         IF (n .lt. 0) THEN
            isign = -1
         END IF

         np = ABS(n)
         DO m = 0, mpol
            mn = mn + 1
            lk = 0

            DO lv = 1, nzeta
               DO lu = 1, ntheta
                  lk = lk + 1

                  cosmni(lk,mn) =       fourier_context%cosmui(lu,m) *         &
                                        fourier_context%cosnv(lv,np)           &
                                - isign*fourier_context%sinmui(lu,m) *         &
                                        fourier_context%sinnv(lv,np)
                  cosmnp(lk,mn) =       fourier_context%cosmu(lu,m) *          &
                                        fourier_context%cosnv(lv,np)           &
                                - isign*fourier_context%sinmu(lu,m) *          &
                                        fourier_context%sinnv(lv,np)

                  IF (lasym) THEN
                     sinmni(lk,mn) =       fourier_context%sinmui(lu,m) *      &
                                           fourier_context%cosnv(lv,np)        &
                                   + isign*fourier_context%cosmui(lu,m) *      &
                                           fourier_context%sinnv(lv,np)
                     sinmnp(lk,mn) =       fourier_context%sinmu(lu,m) *       &
                                           fourier_context%cosnv(lv,np)        &
                                   + isign*fourier_context%cosmu(lu,m) *       &
                                           fourier_context%sinnv(lv,np)
                  END IF
               END DO
            END DO

            IF (m .eq. 0 .and. n .lt. 0) THEN
               cosmni(:,mn) = 0
               cosmnp(:,mn) = 0
               IF (lasym) THEN
                  sinmni(:,mn) = 0
                  sinmnp(:,mn) = 0
               END IF
            END IF
         END DO
      END DO

!  Compute lambda from sqrt(g)*J^s == mB_v - nB_u = 0
      DO js = MAX(nsmin, 2), nsmax
         mn = 0

         DO n = -ntor, ntor
            n_mode = fourier_context%tor_modes(n)
            DO m = 0, mpol

!  Right hand side terms sinmi and -cosmi parts.
               mn = mn + 1
               btemp = (-m*         (chiph(js)*guv(:,js) +                     &
     &                               phiph(js)*gvv(:,js)) +                    &
     &                   n_mode*nfp*(chiph(js)*guu(:,js) +                     &
     &                               phiph(js)*guv(:,js)))                     &
                     / jacobh(:,js)

               brhs(mn,1) = SUM(cosmni(:,mn)*btemp)
               IF (lasym) THEN
                  brhs(mn,2) = -SUM(sinmni(:,mn)*btemp)
               END IF

!  Matrix elements. Coefficients of lambda(ms, js)
               mnp = 0

               DO np = -ntor, ntor
                  np_mode = fourier_context%tor_modes(np)
                  DO mp = 0, mpol
                     mnp = mnp + 1
!  Node the phip term is folded into the solve for lambda.
                     atemp = (m*mp*gvv(:,js) +                                 &
     &                        n_mode*nfp*np_mode*nfp*guu(:,js) -               &
     &                        nfp*(m*np_mode + mp*n_mode)*guv(:,js))           &
     &                     / jacobh(:,js)

                     amat(mn,mnp) = SUM(cosmni(:,mn)*cosmnp(:,mnp)*atemp)
	                 IF (mnp .eq. mn .and. amat(mn,mnp) .eq. zero) THEN
                        amat(mn,mnp) = signjac
                     END IF

                     IF (lasym) THEN
                        amat(mn,mnp+mnmax) =                                   &
                           -SUM(cosmni(:,mn)*sinmnp(:,mnp)*atemp)
                        amat(mn+mnmax,mnp) =                                   &
                           -SUM(sinmni(:,mn)*cosmnp(:,mnp)*atemp)
                        amat(mn+mnmax,mnp+mnmax) =                             &
                           SUM(sinmni(:,mn)*sinmnp(:,mnp)*atemp)
	                    IF (mnp .eq. mn .and.                                  &
                            amat(mn+mnmax,mnp+mnmax).eq. zero) THEN
                           amat(mn+mnmax,mnp+mnmax) = signjac
                        END IF
                     END IF
                  END DO
	           END DO
	        END DO
         END DO

         CALL solver(amat, brhs, mnmax*n2, 1, info)
         CALL ASSERT(info.EQ.0,'INFO != 0 IN RECOMPUTE_LAMBDA')

!  The we really solved for phip*lambda so divide off that component.
         lmns(:,js) = brhs(:,1)/(orthonorm*phiph(js))
         IF (lasym) THEN
            lmnc(:,js) = brhs(:,2)/(orthonorm*phiph(js))
         END IF
      END DO

      CALL assert(mn .eq. mnmax .and. mnp .eq. mnmax,                          &
                  'mn or mnp != mnmax in RECOMPUTE_LAMBDA')

      DEALLOCATE(atemp, btemp, cosmni, cosmnp, brhs, amat)
      IF (lasym) THEN
         DEALLOCATE(sinmni, sinmnp)
      END IF

      CALL second0(toff)
      IF (iam .EQ. 0 .and. lverbose) THEN
         WRITE (*,1000) toff - ton
      END IF

1000  FORMAT(' RECOMPUTE_LAMBDA - TIME: ',f10.2,' S')
      END SUBROUTINE

!-------------------------------------------------------------------------------
!>  @brief Allocates dependant variables.
!>
!>  All variables allocated here are then deallocated in dealloc_quantities,
!>  which should be called at the end of the run.
!-------------------------------------------------------------------------------
      SUBROUTINE alloc_quantities

      IMPLICIT NONE

!  local variables
      INTEGER :: istat

!  Start of executable code
!  Half mesh quantities (harmonics), B^s (sine), B^u (cosine), B^v (cosine)
      IF (.not.ALLOCATED(jbsupsmnsh)) THEN
         ALLOCATE(jbsupsmnsh(0:mpol,-ntor:ntor,ns),                            &
                  jbsupumnch(0:mpol,-ntor:ntor,ns),                            &
                  jbsupvmnch(0:mpol,-ntor:ntor,ns),                            &
                  jpmnch    (0:mpol,-ntor:ntor,ns), stat=istat)
         CALL assert_eq(0, istat, 'Allocation #1 failed in alloc_quantities')
      END IF
      jbsupsmnsh = zero
      jbsupumnch = zero
      jbsupvmnch = zero
      jpmnch = zero

      IF (.not.ALLOCATED(pwr_spec_s)) THEN
         ALLOCATE(pwr_spec_s(0:mpol,-ntor:ntor,ns,4),                          &
                  pwr_spec_a(0:mpol,-ntor:ntor,ns,4), stat=istat)
         CALL assert_eq(0, istat, 'Allocation #2 failed in alloc_quantities')
      END IF

      IF (lasym) THEN
!  Half mesh quantities (harmonics), B^s (sine), B^u (cosine), B^v (cosine)
         IF (.not.ALLOCATED(jbsupsmnch)) THEN
            ALLOCATE(jbsupsmnch(0:mpol,-ntor:ntor,ns),                         &
                     jbsupumnsh(0:mpol,-ntor:ntor,ns),                         &
                     jbsupvmnsh(0:mpol,-ntor:ntor,ns),                         &
                     jpmnsh    (0:mpol,-ntor:ntor,ns), stat=istat)
            CALL assert_eq(0, istat,                                           &
     &                     'Allocation #3 failed in alloc_quantities')
         END IF
         jbsupsmnch = zero
         jbsupumnsh = zero
         jbsupvmnsh = zero
         jpmnsh = zero
      END IF

      IF (.not.ALLOCATED(bsq)) THEN
         ALLOCATE(bsq(ntheta,nzeta,ns), jacobh(ntheta,nzeta,ns),               &
                  jacobf(ntheta,nzeta,ns), wint(ntheta,nzeta,ns), stat=istat)
         CALL assert_eq(0, istat, 'Allocation #4 failed in alloc_quantities')
      END IF
      bsq = 0.0

      END SUBROUTINE

!-------------------------------------------------------------------------------
!>  @brief Deallocates dependant variables.
!>
!>  Deallocates all variables in @ref alloc_quantities and some additional
!>  variables.
!-------------------------------------------------------------------------------
      SUBROUTINE dealloc_quantities

      IMPLICIT NONE

!  Start of executable code
!  Quantities allocated in alloc_quantities
      IF (ALLOCATED(jpmnch)) THEN
         DEALLOCATE(jpmnch, jbsupsmnsh, jbsupumnch, jbsupvmnch)
      END IF
      IF (ALLOCATED(pwr_spec_s)) THEN
         DEALLOCATE(pwr_spec_s, pwr_spec_a)
      END IF

!  Other quantities
      IF (ALLOCATED(jpmnch)) THEN
         DEALLOCATE(djpmnch, djbsupsmnsh, djbsupumnch, djbsupvmnch,            &
                    ksupsmnsf, ksupumncf, ksupvmncf)
      END IF
      IF (ALLOCATED(jpmnch)) THEN
         DEALLOCATE(jacobh, jacobf, jvsupsijf, jvsupuijf, jvsupvijf,           &
                    ksubsijf, ksubuijf, ksubvijf, pijh0, pijf0,                &
                    pijf0_ds, bsupuijf0, vp_f, bsq, wint,                      &
                    bsupvijf0, bsupsijf0, pijh0_du, pijh0_dv)
      END IF
      IF (ALLOCATED(bsupsijf0)) THEN
         DEALLOCATE(bsupsijf0, bsupuijf0, bsupvijf0)
      END IF
      IF (ALLOCATED(bsupsijh0)) THEN
         DEALLOCATE(bsupsijh0, bsupuijh0, bsupvijh0)
      END IF
      IF (ALLOCATED(ksupsijf0)) THEN
         DEALLOCATE(ksupsijf0, ksupuijf0, ksupvijf0)
      END IF
      IF (ALLOCATED(ksupsijf)) THEN
         DEALLOCATE(ksupsijf, ksupuijf, ksupvijf)
      END IF
      IF (ALLOCATED(bsubsijf)) THEN
         DEALLOCATE(bsubsijf, bsubuijf, bsubvijf)
      END IF
      IF (ALLOCATED(bsupsijf)) THEN
         DEALLOCATE(bsupsijf, bsupuijf, bsupvijf)
      END IF

!  Asymmetric quantities.
      IF (lasym) THEN
!  Quantities allocated in alloc_quantities
         IF (ALLOCATED(jpmnsh)) THEN
            DEALLOCATE(jpmnsh, jbsupsmnch, jbsupumnsh, jbsupvmnsh)
         END IF
!  Other quantities
         IF (ALLOCATED(djpmnsh)) THEN
            DEALLOCATE(djpmnsh, djbsupsmnch, djbsupumnsh, djbsupvmnsh,         &
                       ksupsmncf, ksupumnsf, ksupvmnsf)
         END IF
      END IF

      END SUBROUTINE

!-------------------------------------------------------------------------------
!>  @brief Routine to compute the contravariant force components.
!-------------------------------------------------------------------------------
      SUBROUTINE toupper_forces
      USE island_params, ONLY: ntheta=>nu_i, nzeta=>nv_i,                      &
                               mpol=>mpol_i, ntor=>ntor_i,                     &
                               fourier_context
      USE fourier, ONLY: f_none, f_cos, f_sin, f_sum, f_con
      USE metrics, ONLY: toupper
      USE utilities, ONLY: set_bndy_full_origin
      USE hessian, ONLY: gather_array

!  Local Variables
      REAL(dp), ALLOCATABLE, DIMENSION(:,:,:) :: fsubsijf
      REAL(dp), ALLOCATABLE, DIMENSION(:,:,:) :: fsubuijf
      REAL(dp), ALLOCATABLE, DIMENSION(:,:,:) :: fsubvijf
      REAL(dp), ALLOCATABLE, DIMENSION(:,:,:) :: fsupsijf
      REAL(dp), ALLOCATABLE, DIMENSION(:,:,:) :: fsupuijf
      REAL(dp), ALLOCATABLE, DIMENSION(:,:,:) :: fsupvijf
      INTEGER                                 :: nsmin
      INTEGER                                 :: nsmax
      INTEGER                                 :: istat

!  Start of executable code.
      nsmin = MAX(1, startglobrow)
      nsmax = MIN(ns, endglobrow)

      ALLOCATE(fsubsijf(ntheta,nzeta,nsmin:nsmax),                             &
               fsubuijf(ntheta,nzeta,nsmin:nsmax),                             &
               fsubvijf(ntheta,nzeta,nsmin:nsmax),                             &
               fsupsijf(ntheta,nzeta,nsmin:nsmax),                             &
               fsupuijf(ntheta,nzeta,nsmin:nsmax),                             &
               fsupvijf(ntheta,nzeta,nsmin:nsmax), stat=istat)
      CALL ASSERT(istat.eq.0,' Allocation error in toupper_forces')

      CALL fourier_context%toijsp(fsubsmncf(:,:,nsmin:nsmax), fsubsijf,        &
                                  f_none, f_cos)
      CALL fourier_context%toijsp(fsubumnsf(:,:,nsmin:nsmax), fsubuijf,        &
                                  f_none, f_sin)
      CALL fourier_context%toijsp(fsubvmnsf(:,:,nsmin:nsmax), fsubvijf,        &
                                  f_none, f_sin)

      IF (lasym) THEN
         CALL fourier_context%toijsp(fsubsmnsf(:,:,nsmin:nsmax), fsubsijf,     &
                                     f_sum, f_sin)
         CALL fourier_context%toijsp(fsubumncf(:,:,nsmin:nsmax), fsubuijf,     &
                                     f_sum, f_cos)
         CALL fourier_context%toijsp(fsubvmncf(:,:,nsmin:nsmax), fsubvijf,     &
                                     f_sum, f_cos)
      END IF

      IF (nsmin .eq. 1) THEN
         fsubuijf(:,:,1) = 0
      END IF

!  Compute contravariant forces needed to compute <||F||^2>.
      CALL toupper(fsubsijf, fsubuijf, fsubvijf,                               &
                   fsupsijf, fsupuijf, fsupvijf, nsmin, nsmax) 
      DEALLOCATE(fsubsijf, fsubuijf, fsubvijf)

      CALL fourier_context%tomnsp(fsupsijf, fsupsmncf(:,:,nsmin:nsmax), f_cos)
      CALL fourier_context%tomnsp(fsupuijf, fsupumnsf(:,:,nsmin:nsmax), f_sin)
      CALL fourier_context%tomnsp(fsupvijf, fsupvmnsf(:,:,nsmin:nsmax), f_sin)

      IF (nsmin .eq. 1) THEN
         CALL set_bndy_full_origin(fsupsmncf, fsupumnsf, fsupvmnsf, f_con)
      END IF

      CALL gather_array(fsupsmncf)
      CALL gather_array(fsupumnsf)
      CALL gather_array(fsupvmnsf)

      IF (lasym) THEN
         CALL fourier_context%tomnsp(fsupsijf, fsupsmnsf(:,:,nsmin:nsmax),     &
                                     f_sin)
         CALL fourier_context%tomnsp(fsupuijf, fsupumncf(:,:,nsmin:nsmax),     &
                                     f_cos)
         CALL fourier_context%tomnsp(fsupvijf, fsupvmncf(:,:,nsmin:nsmax),     &
                                     f_cos)

         IF (nsmin .eq. 1) THEN
            CALL set_bndy_full_origin(fsupsmnsf, fsupumncf, fsupvmncf, f_con)
         END IF

         CALL gather_array(fsupsmnsf)
         CALL gather_array(fsupumncf)
         CALL gather_array(fsupvmncf)
      END IF

      DEALLOCATE(fsupsijf, fsupuijf, fsupvijf, stat=istat)

      END SUBROUTINE

!-------------------------------------------------------------------------------
!>  @brief Gather fields from all processors if running in parallel.
!-------------------------------------------------------------------------------
      SUBROUTINE GATHERFIELDS

      IMPLICIT NONE

!  Start of executable code
      IF (PARSOLVER) THEN
         CALL GATHER_FIELDS(jbsupsmnsh, jbsupumnch, jbsupvmnch, jpmnch)
         IF (lasym) THEN
            CALL GATHER_FIELDS(jbsupsmnch, jbsupumnsh, jbsupvmnsh, jpmnsh)
         END IF
      END IF
      
      END SUBROUTINE

!-------------------------------------------------------------------------------
!>  @brief Gather fields for a partity.
!>
!>  @param[inout] jbsupsmnh Fouier amplitudes of the contravariant magnetic
!>                          field and jacobian in the s direction.
!>  @param[inout] jbsupumnh Fouier amplitudes of the contravariant magnetic
!>                          field and jacobian in the udirection.
!>  @param[inout] jbsupvmnh Fouier amplitudes of the contravariant magnetic
!>                          field and jacobian in the vdirection.
!>  @param[inout] jpmnh     Fouier amplitudes of the pressure and jacobian.
!-------------------------------------------------------------------------------
      SUBROUTINE GATHER_FIELDS(jbsupsmnh, jbsupumnh, jbsupvmnh, jpmnh)
      USE hessian, ONLY: gather_array

      IMPLICIT NONE

!  Start of executable code
      REAL (dp), DIMENSION(:,:,:), INTENT(inout) :: jbsupsmnh
      REAL (dp), DIMENSION(:,:,:), INTENT(inout) :: jbsupumnh
      REAL (dp), DIMENSION(:,:,:), INTENT(inout) :: jbsupvmnh
      REAL (dp), DIMENSION(:,:,:), INTENT(inout) :: jpmnh

!  Start of executable code.
      CALL assert_eq(ns, SIZE(jbsupsmnh,3),                                    &
                     'Radial array size needs to be ns in GATHER_FIELDS')

      CALL gather_array(jbsupsmnh)
      CALL gather_array(jbsupumnh)
      CALL gather_array(jbsupvmnh)
      CALL gather_array(jpmnh)

      END SUBROUTINE

!-------------------------------------------------------------------------------
!>  @brief Get the surface average of quantity.
!>
!>  @param[inout] jbsupsmnh Fouier amplitudes of the contravariant magnetic
!>                          field and jacobian in the s direction.
!>  @param[inout] jbsupumnh Fouier amplitudes of the contravariant magnetic
!>                          field and jacobian in the udirection.
!>  @param[inout] jbsupvmnh Fouier amplitudes of the contravariant magnetic
!>                          field and jacobian in the vdirection.
!>  @param[inout] jpmnh     Fouier amplitudes of the pressure and jacobian.
!-------------------------------------------------------------------------------
      SUBROUTINE SurfAverage(average, q3d, nsmin, nsmax)
      USE stel_kinds, ONLY: dp
      USE island_params, ONLY: ns=>ns_i, ntheta=>nu_i, nzeta=>nv_i

      IMPLICIT NONE

!  Declare Arguments.
      REAL(dp), DIMENSION(nsmin:nsmax), INTENT(out)             :: average
      REAL(dp), DIMENSION(ntheta,nzeta,nsmin:nsmax), INTENT(in) :: q3d
      INTEGER, INTENT(in)                                       :: nsmin
      INTEGER, INTENT(in)                                       :: nsmax

!  Local variables.
      INTEGER                                                   :: js

!  Start of executable code.
      DO js = nsmin, nsmax
		 average(js) = SUM(q3d(:,:,js)*wint(:,:,js))
      END DO

      END SUBROUTINE
     
      END MODULE
