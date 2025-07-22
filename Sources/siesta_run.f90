!*******************************************************************************
!>  @file siesta_run.f90
!>  @brief Contains module @ref siesta_run.
!
!  Note separating the Doxygen comment block here so detailed decription is
!  found in the Module not the file.
!
!>  Defines the base class of the type @ref siesta_run_class. This module
!>  contains all the code necessary to interface setup and run SIESTA.
!*******************************************************************************

      MODULE siesta_run
      USE stel_kinds
      USE shared_data, ONLY: unit_out
      USE timer_mod

      IMPLICIT NONE

!*******************************************************************************
!  siesta run module parameters
!*******************************************************************************
!>  Clear all control bits.
      INTEGER, PARAMETER :: siesta_run_control_clear = 0

!>  Bit position to initalize mpi.
      INTEGER, PARAMETER :: siesta_run_control_mpi = 1
!>  Bit position to close the wout file.
      INTEGER, PARAMETER :: siesta_run_control_wout = 2
!>  Bit position to sync the namelist inputs.
      INTEGER, PARAMETER :: siesta_run_sync_namelist = 3
!>  Bit position to write diagnostics outputs.
      INTEGER, PARAMETER :: siesta_run_write_profiles = 4

!*******************************************************************************
!  DERIVED-TYPE DECLARATIONSf
!  1) siesta_run base class
!
!*******************************************************************************
      TYPE :: siesta_run_class
!>  Timer On
         REAL (dp) :: time_on
!>  Control state.
         INTEGER   :: control_state
      CONTAINS
!  FIXME: Remove _temp after merge of V3FIT fixes.
         PROCEDURE :: set_vmec_ => siesta_run_set_vmec
         PROCEDURE :: set_vmec_temp => siesta_run_set_vmec_temp
         GENERIC   :: set_vmec => set_vmec_, set_vmec_temp
         PROCEDURE :: set_restart => siesta_run_set_restart
         PROCEDURE :: set_1d => siesta_run_set_1d
         GENERIC   :: set => set_1d
         PROCEDURE :: set_control_state => siesta_run_set_control_state
         PROCEDURE :: is_state_set => siesta_run_is_state_set
         PROCEDURE :: clear_state => siesta_run_clear_state
         PROCEDURE :: get_1d => siesta_run_get_1d
         GENERIC   :: get => get_1d
         PROCEDURE :: converge => siesta_run_converge
         PROCEDURE :: sync => siesta_run_sync
         FINAL     :: siesta_run_destruct
      END TYPE

!*******************************************************************************
!  INTERFACE BLOCKS
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  Interface for the siesta_run constructor.
!-------------------------------------------------------------------------------
      INTERFACE siesta_run_class
         MODULE PROCEDURE siesta_run_construct
      END INTERFACE

      CONTAINS

!*******************************************************************************
!  CONSTRUCTION SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  @brief Construct new @ref siesta_run_class object.
!>
!>  Allocates memory and initializes a @ref siesta_run_class object. Performs
!>  all the initialization needed to operate SIESTA.
!>
!>  @param[in] run_comm      MPI Communicator to use.
!>  @param[in] verbose       Control the code screen output.
!>  @param[in] init_mpi      Instructs if MPI should be initialized.
!>  @param[in] close_wout    Instructs if wout file should be closed on
!>                           deallocate.
!>  @param[in] namelist_file Path to the siesta namelist input file.
!>  @returns A pointer to a constructed @ref siesta_run_class object.
!-------------------------------------------------------------------------------
      FUNCTION siesta_run_construct(run_comm, verbose, init_mpi, close_wout,   &
                                    namelist_file)
      USE nscalingtools, ONLY: MyEnvVariables, PARSOLVER, SIESTA_COMM
      USE descriptor_mod, ONLY: LSCALAPACK, INHESSIAN, iam, nprocs,            &
                                nprow, npcol, icontxt, icontxt_1xp,            &
                                icontxt_px1, icontxt_global, isroot,           &
                                myrow, mycol
#if defined(MPI_OPT)
      USE prof_mod, ONLY: profinit
#endif      
      USE shared_data, ONLY: nprecon, lverbose
      USE perturbation, ONLY: init_data
      USE hessian, ONLY: HESSPASS
      USE siesta_state, ONLY: lfirst
      USE siesta_error
      USE diagnostics_mod, ONLY: toroidal_flux0
      USE blocktridiagonalsolver_s, ONLY: Initialize, GetRanks
      USE siesta_namelist, ONLY: mpolin, ntorin, nfpin, ntor_modes
      USE metrics, ONLY: set_grid_sizes

      IMPLICIT NONE

!  Declare Arguments
      CLASS (siesta_run_class), POINTER :: siesta_run_construct
      INTEGER, INTENT(in)               :: run_comm
      LOGICAL, INTENT(in)               :: verbose
      LOGICAL, INTENT(in)               :: init_mpi
      LOGICAL, INTENT(in)               :: close_wout
      CHARACTER (len=*), INTENT(in)     :: namelist_file

!  local variables.
      REAL (dp)                         :: ton
      REAL (dp)                         :: toff

!  Start of executable code
      ALLOCATE(siesta_run_construct)

      siesta_run_construct%control_state = siesta_run_control_clear

      CALL siesta_run_construct%set_control_state(siesta_run_write_profiles)
      IF (init_mpi) THEN
         CALL siesta_run_construct%set_control_state(siesta_run_control_mpi)
      END IF

      IF (close_wout) THEN
         CALL siesta_run_construct%set_control_state(siesta_run_control_wout)
      END IF

      CALL siesta_error_clear_all

      lverbose = verbose
      lfirst = .TRUE.
      toroidal_flux0=0

      SIESTA_COMM = run_comm
      HESSPASS = 0
  	  PARSOLVER=.FALSE.

#if defined(MPI_OPT)
      LSCALAPACK = .TRUE.
      INHESSIAN = .FALSE.

      CALL MyEnvVariables
      IF (PARSOLVER) THEN
         LSCALAPACK=.FALSE.
         CALL Initialize(0, 0)        !Must call this to initialize mpi timer
         CALL GetRanks(iam, nprocs)
         CALL sks_timers
      END IF

      IF (LSCALAPACK) THEN
         CALL profinit()

! setup blacs
         CALL blacs_pinfo(iam, nprocs)
         CALL blacs_setup(iam, nprocs)
         DO nprow = INT(SQRT(DBLE(nprocs))) + 1, 1, -1
            npcol = INT( nprocs/nprow )
            IF (nprow*npcol .eq. nprocs) EXIT
         END DO

         CALL blacs_get(0, 0, icontxt)
         CALL blacs_gridinit(icontxt, 'C', nprow, npcol)

         CALL blacs_get(0, 0, icontxt_1xp)
         CALL blacs_gridinit(icontxt_1xp, 'R', 1, nprow*npcol)

         CALL blacs_get(0, 0, icontxt_px1)
         CALL blacs_gridinit(icontxt_px1, 'C', nprow*npcol, 1)

         icontxt_global = icontxt_1xp

         CALL blacs_gridinfo(icontxt, nprow, npcol, myrow, mycol)
         isroot = (myrow .eq. 0) .and. (mycol .eq. 0)
         IF (isroot .and. lverbose) THEN
            WRITE (*,1000) nprow, npcol
            WRITE (*,1001) icontxt_global, icontxt, icontxt_1xp, icontxt_px1
         END IF
      END IF
#else
      LSCALAPACK = .FALSE.
      iam = 0
      nprocs = 1
#endif

      CALL second0(ton)
      siesta_run_construct%time_on = ton

      nprecon = init_data(namelist_file)
      CALL set_grid_sizes(mpolin, ntorin, nfpin, ntor_modes(-ntorin:ntorin))

      CALL second0(toff)

      init_data_time = toff - ton

1000  FORMAT('nprow,npcol ',2(x,i5))
1001  FORMAT('icontxt_global,icontxt,icontxt_1xp,icontxt_px1 ',4(x,i5))

      END FUNCTION

!*******************************************************************************
!  DESTRUCTION SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  @brief Deconstruct a @ref siesta_run_class object.
!>
!>  Deallocates memory and uninitializes a @ref siesta_run_class object. This 
!>  also performs the shut down procedures for SIESTA.
!>
!>  @param[inout] this         A @ref siesta_run_class instance.
!-------------------------------------------------------------------------------
      SUBROUTINE siesta_run_destruct(this)
      USE nscalingtools, ONLY: OUTPUT_TIMINGS, GetTimes, PARSOLVER,     &
                               PARFUNCTISL, PARGMRES, SKSDBG, TOFU,     &
                               finalizeremap
      USE blocktridiagonalsolver_s, ONLY: finalize
      USE descriptor_mod, ONLY: iam, icontxt, LSCALAPACK, myrow, mycol, nprocs
      USE shared_data, ONLY: nprecon, xc, gc, fsq_total1, mblk_size,    &
                             niter, col_scale, lverbose, ngmres_type,   &
                             lColScale, gc_sup
      USE hessian, ONLY: levmarq_param, mupar, asym_index, dealloc_hessian
      USE quantities, ONLY: fbdy, dealloc_quantities
      USE dumping_mod, ONLY: write_output
      USE diagnostics_mod, ONLY: write_profiles, dealloc_diagnostics
      USE siesta_namelist, ONLY: wout_file
      USE island_params, ONLY: ns_i
      USE metrics, ONLY: cleanup_metric_elements, dealloc_full_lower_metrics
#if defined(MPI_OPT)
      USE prof_mod, ONLY: profstat
#endif
      IMPLICIT NONE

!  Declare Arguments
      TYPE (siesta_run_class), INTENT(inout) :: this

!  local variables.
      REAL (dp)                              :: ton
      REAL (dp)                              :: toff
      INTEGER                                :: i

!  Start of executable code
      CALL second0(toff)
      time_total = toff - this%time_on

      IF (OUTPUT_TIMINGS) THEN
         CALL GetTimes
      END IF

      IF (iam .EQ. 0) THEN
         DO i = 6, unit_out, unit_out-6
            IF (.NOT.lverbose .AND. i.EQ.6) CYCLE
            WRITE (i, 1000) nprecon, levmarq_param, mupar, asym_index
            WRITE (i,'(a,i5)') ' Number processors: ', nprocs
            WRITE (i, 1001)
            WRITE (i, 1002) time_total,               fbdy(1),          &
                            time_init,                fbdy(2),          &
                            time_diag_prec,           fbdy(3),          &
                            time_block_prec,          fbdy(4),          &
                            time_factor_blocks,       fbdy(5),          &
                            time_toijsp,              fbdy(6),          &
                            time_tomnsp,                                &
                            gmres_time,                                 &
                            conj_grad_time,                             &
                            time_update_pres,         fbdy(7),          &
                            time_update_bfield,       fbdy(8),          &
                            time_current,             fbdy(9),          &
                            get_force_harmonics_time, fbdy(10),         &
                            time_update_force,        fbdy(11),         &
                            time_update_upperv,       fbdy(12),         &
                            time_update_state,        fbdy(13),         &
                            time_funci,                                 &
                            time_apply_precon,                          &
                            (diag_add_pert_time + block_add_pert_time), &
                            time_init_state
            WRITE (i,*)
            IF (PARSOLVER) THEN
               WRITE (i,*)        'PARSOLVER=T    : NSCALED'
            ELSE IF (LSCALAPACK) THEN
               WRITE (i,*)        'PARSOLVER=F    : SCALAPACK'
            ELSE 
               WRITE (i,*)        'PARSOLVER=F    : SERIAL'
            END IF
            WRITE (i,*)           'PARFUNCTISL    :', PARFUNCTISL
            WRITE (i,*)           'COLUMN SCALING :', lColScale
            WRITE (i,'(1x,a,L2,a,i2)') 'PARGMRES       :', PARGMRES,    &
                                  ' GMRES_TYPE: ', ngmres_type

            WRITE (i,*)           'OUTPUT_TIMINGS :', OUTPUT_TIMINGS

            WRITE (i,*)
            WRITE (i,101)' TIME DIVB      : ', time_divb
            WRITE (i,101)' TIME DIVJ      : ', time_divj
            WRITE (i,101)' TIME BGRADP    : ', time_bgradp
            WRITE (i,101)' TIME BDOTJ     : ', time_bdotj
            WRITE (i,*)
            WRITE (i,102)' M (block size) :', mblk_size
            WRITE (i,102)' N (block rows) :', ns_i
            WRITE (i,102)' P (processors) :', nprocs
         END DO
      ENDIF
 101  FORMAT(a,1p,e10.3)
 102  FORMAT(a,i5)
 
      CALL write_output(wout_file, niter,                                      &
                        this%is_state_set(siesta_run_control_wout))
      IF (this%is_state_set(siesta_run_write_profiles)) THEN
         CALL write_profiles(fsq_total1)  ! SPH: write pmn, bsupXmn, ksubXmn, jvsupXmn profiles
      END IF
      IF (iam .EQ. 0) THEN
         IF (lverbose) PRINT *,' Writing output to "siesta_profiles.txt" is finished!'
         CLOSE (UNIT=unit_out)
      ENDIF

!CLEAN UP ARRAYS
      CALL cleanup_metric_elements

      CALL dealloc_quantities
      CALL dealloc_hessian
      CALL dealloc_diagnostics
      CALL dealloc_full_lower_metrics
      IF (ALLOCATED(xc)) DEALLOCATE(xc, gc, gc_sup)
      IF (ALLOCATED(col_scale)) DEALLOCATE(col_scale)

#if defined(MPI_OPT)
      IF (LSCALAPACK) THEN
         CALL blacs_barrier(icontxt, 'All')
         CALL blacs_gridexit(icontxt)

         CALL blacs_exit(.not.this%is_state_set(siesta_run_control_mpi))
         IF ((myrow.EQ.0) .AND. (mycol.EQ.0)) THEN
            CALL profstat
         END IF

      ELSE
         IF (SKSDBG) THEN
            WRITE(TOFU,*) 'Called finalizeRemap and Finalize'
            FLUSH(TOFU)
         END IF
         CALL finalizeRemap
         CALL Finalize(this%is_state_set(siesta_run_control_mpi))

      END IF
#endif

1000  FORMAT(/,' nprecon: ',i3,' LM parameter: ',1pe9.2,' mu||: ',1pe9.2,      &
             ' Symmetry Index: ',1pe9.2)
1001  FORMAT(/,'==============================',14x,'=======================', &
             /,                                                                &
             /,' TIMING INFORMATION           ',14x,' RMS BOUNDARY FORCES',    &
             /,                                                                &
             /,'==============================',14x,'=======================')
1002  FORMAT(' Total runtime  : ',  f12.3,15x,'fs(1,m=1)  :',1pe10.2/,         &
             ' Initialization : ',0pf12.3,15x,'fs(2,m=1)  :',1pe10.2/,         &
             ' Diagonal prec  : ',0pf12.3,15x,'fs(2,m!=1) :',1pe10.2/,         &
             ' Compute blocks : ',0pf12.3,15x,'fu(1,m=1)  :',1pe10.2/,         &
             ' Factor blocks  : ',0pf12.3,15x,'fu(2,m=1)  :',1pe10.2/,         &
             ' Toijsp         : ',0pf12.3,15x,'fu(2,m!=1) :',1pe10.2/,         &
             ' Tomnsp         : ',0pf12.3,/,                                   &
             ' GMRES          : ',0pf12.3,/,                                   &
             ' Conj Gradient  : ',0pf12.3,//,                                  &
             ' SUBROUTINES     ',/,                                            &
             ' Update Pressure: ',0pf12.3,15x,'fv(1,m=0)  :',1pe10.2/,         &
             ' Update Bfield  : ',0pf12.3,15x,'fv(2,m=0)  :',1pe10.2/,         &
             ' CV Currents    : ',0pf12.3,15x,'fv(2,m!=0) :',1pe10.2/,         &
             ' Force Harmonics: ',0pf12.3,15x,'fu(ns)     :',1pe10.2/,         &
             ' Update Force   : ',0pf12.3,15x,'fu(ns-1)   :',1pe10.2/,         &
             ' Update UpperV  : ',0pf12.3,15x,'fv(ns)     :',1pe10.2/,         &
             ' Update State   : ',0pf12.3,15x,'fv(ns-1)   :',1pe10.2/,         &
             ' Funct Island   : ',0pf12.3/,                                    &
             ' Apply Precon   : ',0pf12.3/,                                    &
             ' Add Perturb    : ',0pf12.3/,                                    &
             ' Init State     : ',0pf12.3,                                     &
             /,'==============================',14x,'=======================')

      END SUBROUTINE

!*******************************************************************************
!  SETTER SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  @brief Sets quanties based on VMEC values.
!>
!>  This method loads and sets variables based on the VMEC equilibrium. VMEC 
!>  controls the metric elements and coordinate system jacobian.
!>
!>  @param[inout] this      A @ref siesta_run_class instance.
!>  @param[in]    load_wout Flag to load the wout file.
!-------------------------------------------------------------------------------
!  FIXME: Temp routine to make sure CI tests still. Remove once V3FIT changes
!         finished.
      SUBROUTINE siesta_run_set_vmec(this)
      IMPLICIT NONE
      CLASS (siesta_run_class), INTENT(inout) :: this
      CALL siesta_run_set_vmec_temp(this, .true.)
      END SUBROUTINE

      SUBROUTINE siesta_run_set_vmec_temp(this, load_wout)
      USE siesta_namelist, ONLY: nsin, mpolin, ntorin, nfpin, wout_file,       &
                                 l_vessel, ntor_modes
      USE metrics, ONLY: init_metric_elements, LoadGrid, sqrtg
      USE quantities, ONLY: init_quantities, init_fields
      USE island_params, ns=>ns_i, mpol=>mpol_i, ntor=>ntor_i
      USE shared_data, ONLY: nsp
      USE vmec_info
      USE grid_extension, ONLY: grid_extender, read_vessel_file
      USE pchelms
      USE evolution, ONLY: init_evolution
      USE nscalingtools, ONLY: initRemap

      IMPLICIT NONE

!  Declare Arguments
      CLASS (siesta_run_class), INTENT(inout) :: this
      LOGICAL, INTENT(in)                     :: load_wout

!  local variables
      INTEGER                                 :: istat

!  Start of executable code
      CALL this%sync

      CALL vmec_info_set_wout(wout_file, nsin, mpolin, ntorin, nfpin,          &
     &                        ntor_modes(-ntorin:ntorin), load_wout)

!  CONSTRUCT R, Z, L REAL-SPACE ARRAYS ON SQRT(FLUX) - "POLAR" - MESH AND
!  COMPUTE METRIC ELEMENTS AND JACOBIAN
      CALL LoadGrid(istat)
      CALL assert_eq(0, istat, 'LoadRZL error in siesta_run_set_vmec')

!  If the l_vessel is true try to read a vessel file and extend the grid.
!  Otherwise proceed with fixed boundary equilibrium.
      IF (l_vessel) THEN
         IF (read_vessel_file() .eq. 0) THEN
            CALL grid_extender(wout_file, 'quad', istat)
         END IF
      ELSE
         l_vessel = .false.
      END IF

      CALL init_metric_elements()
      CALL init_quantities(.false.) !Initializes BCYCLIC
      CALL init_evolution !neqs is set here
      CALL initRemap(mpol, ntor, ns, nprocs, iam)
      CALL InitHess

      IF (l_vessel) THEN
         CALL run_pchelms
      END IF
      CALL Init_Fields

      DEALLOCATE (sqrtg)

      END SUBROUTINE

!-------------------------------------------------------------------------------
!>  @brief Sets quantities for restarting.
!>
!>  This loads and sets variables from a restarted SIESTA equilibrium. The
!>  restart file constains the metric elements, coordinate system jacobian and
!>  the
!>
!>  @param[inout] this A @ref siesta_run_class instance.
!-------------------------------------------------------------------------------
      SUBROUTINE siesta_run_set_restart(this)
      USE restart_mod, ONLY: restart_read
      USE siesta_namelist, ONLY: lrestart, mpolin, ntorin, nsin, restart_ext,  &
     &                           wout_file, nfpin, nsin_ext, l_vessel,         &
     &                           ntor_modes
      USE shared_data, ONLY: nprecon
      USE metrics, ONLY: init_metric_elements
      USE quantities, ONLY: init_quantities
      USE descriptor_mod, ONLY: iam, nprocs
      USE evolution, ONLY: init_evolution
      USE nscalingtools, ONLY: initRemap
      USE island_params, ONLY: ns=>ns_i, hs=>hs_i, ohs=>ohs_i, nsh
      USE hessian, ONLY: InitHess

      IMPLICIT NONE

!  Declare Arguments
      CLASS (siesta_run_class), INTENT(inout) :: this

!  Start of executable code
      CALL this%sync

      IF (l_vessel) THEN
         ns = nsin + nsin_ext
      ELSE
         ns = nsin
      END IF
      nsh = ns - 1
      ohs = nsin - 1
      hs = 1.0/ohs

      nprecon = restart_read(restart_ext, wout_file, mpolin, ntorin, nfpin,    &
                             ns, ntor_modes(-ntorin:ntorin))

      CALL init_metric_elements()
      CALL init_quantities(.true.) !Initializes BCYCLIC
      CALL init_evolution !neqs is set here
      CALL initRemap(mpolin, ntorin, ns, nprocs, iam)
      CALL InitHess

      END SUBROUTINE

!-------------------------------------------------------------------------------
!>  @brief Set 1D parameter value.
!>
!>  @param[inout] this       @ref siesta_run_class instance.
!>  @param[in]    param_name Name of the parameter to set.
!>  @param[in]    value      Value to set the parameter to.
!>  @param[in]    index      Array index of to set the value to.
!-------------------------------------------------------------------------------
      SUBROUTINE siesta_run_set_1d(this, param_name, value, index)
      USE siesta_namelist, ONLY: helpert, helphase
      USE siesta_error

      IMPLICIT NONE

!  Declare Arguments
      CLASS (siesta_run_class), INTENT(inout) :: this
      CHARACTER (len=*), INTENT(in)           :: param_name
      REAL (dp)                               :: value
      INTEGER, INTENT(in)                     :: index

!  Start of executable code.
      SELECT CASE (TRIM(param_name))

         CASE ('helpert')
            helpert(index) = value
            CALL this%set_control_state(siesta_run_sync_namelist)

         CASE ('helphase')
            helphase(index) = value
            CALL this%set_control_state(siesta_run_sync_namelist)

         CASE DEFAULT
            CALL siesta_error_set_error(siesta_error_param,                    &
                                        'Warning unknown parameter name ' //   &
                                        TRIM(param_name))

      END SELECT

      END SUBROUTINE

!-------------------------------------------------------------------------------
!>  @brief Set a control state.
!>
!>  @param[inout] this  @ref siesta_run_class instance.
!>  @param[in]    state State bit to set.
!-------------------------------------------------------------------------------
      SUBROUTINE siesta_run_set_control_state(this, state)

      IMPLICIT NONE

!  Declare Arguments
      CLASS (siesta_run_class), INTENT(inout) :: this
      INTEGER, INTENT(in)                     :: state

!  Start of executable code.
      this%control_state = IBSET(this%control_state, state)

      END SUBROUTINE

!*******************************************************************************
!  GETTER SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  @brief Get 1D parameter value.
!>
!>  @param[inout] this       @ref siesta_run_class instance.
!>  @param[in]    param_name Name of the parameter to set.
!>  @param[in]    index      Array index of to set the value to.
!>  @returns Value of the parameter at the index.
!-------------------------------------------------------------------------------
      FUNCTION siesta_run_get_1d(this, param_name, index)
      USE siesta_namelist, ONLY: helpert, helphase

      IMPLICIT NONE

!  Declare Arguments
      REAL (dp)                            :: siesta_run_get_1d
      CLASS (siesta_run_class), INTENT(in) :: this
      CHARACTER (len=*), INTENT(in)        :: param_name
      INTEGER, INTENT(in)                  :: index

!  Start of executable code.
      SELECT CASE (TRIM(param_name))

         CASE ('helpert')
            siesta_run_get_1d = helpert(index)

         CASE ('helphase')
            siesta_run_get_1d = helphase(index)

         CASE DEFAULT
            siesta_run_get_1d = 0.0

      END SELECT

      END FUNCTION

!*******************************************************************************
!  QUERY SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  @brief Query if state is set.
!>
!>  @param[in] this  @ref siesta_run_class instance.
!>  @param[in] state State bit to query.
!>  @returns True if the state is set.
!-------------------------------------------------------------------------------
      FUNCTION siesta_run_is_state_set(this, state)

      IMPLICIT NONE

!  Declare Arguments
      LOGICAL                              :: siesta_run_is_state_set
      CLASS (siesta_run_class), INTENT(in) :: this
      INTEGER, INTENT(in)                  :: state

!  Start of executable code.
      siesta_run_is_state_set = BTEST(this%control_state, state)

      END FUNCTION

!*******************************************************************************
!  UTILITY SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  @brief Clear state bit.
!>
!>  @param[in] this  @ref siesta_run_class instance.
!>  @param[in] state State bit to query.
!-------------------------------------------------------------------------------
      SUBROUTINE siesta_run_clear_state(this, state)

      IMPLICIT NONE

!  Declare Arguments
      CLASS (siesta_run_class), INTENT(inout) :: this
      INTEGER, INTENT(in)                     :: state

!  Start of executable code.
      this%control_state = IBCLR(this%control_state, state)

      END SUBROUTINE

!-------------------------------------------------------------------------------
!>  @brief Solves the SIESTA equilibrium.
!>
!>  This method converges the siesta_run_class equilibrium. Initializes the 
!>  evolution module variables and converges the force residual using a diagonal 
!>  preconditioner. Applies the external island perturbation if possible. Loads 
!>  restart file before convergence.
!>
!>  @param[inout] this A @ref siesta_run_class instance.
!-------------------------------------------------------------------------------
      SUBROUTINE siesta_run_converge(this)
      USE evolution, ONLY: converge_diagonal, converge_blocks, pert_added
      USE descriptor_mod, ONLY: DIAGONALDONE
      USE siesta_namelist, ONLY: ftol, wout_file, ladd_pert
      USE utilities, ONLY: test_utilities
      USE shared_data, ONLY: xc

      IMPLICIT NONE

!  Declare Arguments
      CLASS (siesta_run_class), INTENT(inout) :: this

!  Local variables
      LOGICAL                                 :: passed

!  Start of executable code
      pert_added = .not.ladd_pert

      IF (test_utilities()) THEN
         WRITE (*,*) 'Failed Diverence Test.'
!         STOP
      END IF

!  Converge initial residues with diagonal preconditioner
      xc = 0
      DIAGONALDONE = .false.
      CALL converge_diagonal(wout_file, ftol)
      DIAGONALDONE = .true.

!  Converge using block preconditioner
      CALL converge_blocks(wout_file, ftol)

      END SUBROUTINE

!*******************************************************************************
!  MPI SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  @brief Sync child equilibria.
!>
!>  @param[inout] this A @ref siesta_run_class instance.
!-------------------------------------------------------------------------------
      SUBROUTINE siesta_run_sync(this)
      USE siesta_namelist, ONLY: helpert, helphase
      USE nscalingtools, ONLY: SIESTA_COMM, MPI_ERR
#if defined(MPI_OPT)
      USE mpi_inc
#endif

      IMPLICIT NONE

!  Declate Arguments
      CLASS (siesta_run_class), INTENT(inout) :: this

!  Start of executable code
#if defined(MPI_OPT)
      CALL MPI_BCAST(this%control_state, 1, MPI_INTEGER, 0, SIESTA_COMM,       &
                     MPI_ERR)
      IF (this%is_state_set(siesta_run_sync_namelist)) THEN
         CALL MPI_BCAST(helpert, SIZE(helpert), MPI_REAL8, 0, SIESTA_COMM,     &
                        MPI_ERR)
         CALL MPI_BCAST(helphase, SIZE(helphase), MPI_REAL8, 0, SIESTA_COMM,   &
                        MPI_ERR)
         CALL this%clear_state(siesta_run_sync_namelist)
      END IF
#endif

      END SUBROUTINE

      END MODULE
