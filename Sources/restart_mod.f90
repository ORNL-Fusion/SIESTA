!-------------------------------------------------------------------------------
!  The @header, @table_section, @table_subsection, @item and @end_table commands
!  are custom defined commands in Doxygen.in. They are defined under ALIASES.
!  For the page created here, the 80 column limit is exceeded. Arguments of
!  aliases are separated by ','. If you intended ',' to be a string you must use
!  an escaped comma '\,'.
!
!>  @page siesta_restart_sec Discription of the SIESTA restart file.
!>
!>  @tableofcontents
!>  @section siesta_restart_intro_sec Introduction
!>  This page documents the contents of a SIESTA restart file.
!>
!>  @section siesta_restart_var_sec Namelist Variables
!>  @header{Input variable, Description, Code Reference}
!>
!>  @table_section{siesta_restart_flag_sec, Control Flags}
!>     @item{state_flags, State flags of the model. Aditionally stores the
!>                        version number in the first 31 bits.
!>                        * restart_version
!>                        * restart_lasym,               }
!>     @item{wout_file,   Name of the inital wout_file., siesta_namelist::wout_file}
!>  @end_table
!>
!>  @table_section{siesta_restart_size_sec, Dimension sizes.}
!>     @item{mpol, Number of torodial modes.,  siesta_namelist::mpolin}
!>     @item{nfp,  Number of field periods.,   siesta_namelist::nfpin}
!>     @item{nrad, Number of radial surfaces., island_params::ns_i}
!>     @item{ntor, Maximum torodal mode.,      siesta_namelist::ntorin}
!>  @end_table
!>
!>  @table_section{siesta_restart_norm_fact_sec, Normalization factors.}
!>     @item{b_factor, Normalization factor for the magnetic field., quantities::b_factor}
!>     @item{p_factor, Normalization factor for the pressure.,       quantities::p_factor}
!>  @end_table
!>
!>  @table_section{siesta_restart_scalar_sec, Scalar quantities.}
!>     @item{curtor, The toroidal current.,                shared_data::siesta_curtor}
!>     @item{p_max,  Maximum pressure.,                    }
!>     @item{p_min,  Minimum pressure.,                    }
!>     @item{rmajor, Major radius.,                        island_params::rmajor_i}
!>     @item{wb,     Energy stored in the magnetid field., quantities::wb}
!>     @item{wp,     Energy stored in the pressure.,       quantities::wb}
!>  @end_table
!>
!>  @table_section{siesta_restart_1D_arrays_sec, 1D profiles.}
!>     @item{chipf_r_, Radial derivative of the poloidal flux., island_params::chipf_i}
!>     @item{phipf_r_, Radial derivative of the toroidal flux., island_params::phipf_i}
!>  @end_table
!>
!>  @table_sections{siesta_restart_2D_arrays_sec, 3D arrays.}
!>     @table_subsection{siesta_restart_internal_arrays_sec, Internal arrays}.
!>        @item{JBsupssh_m_n_r_, Normalized JB^s component sine parity.,   quantities::jbsupsmnsh}
!>        @item{JBsupuch_m_n_r_, Normalized JB^u component cosine parity., quantities::jbsupumnch}
!>        @item{JBsupvch_m_n_r_, Normalized JB^v component cosine parity., quantities::jbsupvmnch}
!>        @item{jpresch_m_n_r_,  Normalized JP   component cosine parity., quantities::jpmnch}
!>     @table_subsection{siesta_restart_internal_arrays_asym_sec, Internal arrays asym.}
!>        @item{JBsupsch_m_n_r_, Normalized JB^s component cosine parity., quantities::jbsupsmnch}
!>        @item{JBsupush_m_n_r_, Normalized JB^u component sine parity.,   quantities::jbsupumnsh}
!>        @item{JBsupvsh_m_n_r_, Normalized JB^v component sine parity.,   quantities::jbsupvmnsh}
!>        @item{jpressh_m_n_r_,  Normalized JB^v component sine parity.,   quantities::jpmnsh}
!>     @table_subsection{siesta_restart_grid_arrays_sec, Grid arrays}.
!>        @item{rmnc_m_n_r_, R cosine parity., vmec_info::rmnc_i}
!>        @item{zmns_m_n_r_, Z sine parity.,   vmec_info::zmns_i}
!>     @table_subsection{siesta_restart_grid_arrays_asym_sec, Gird arrays asym.}
!>        @item{rmns_m_n_r_, R sine parity.,   vmec_info::rmns_i}
!>        @item{zmnc_m_n_r_, Z cosine parity., vmec_info::zmnc_i}
!>     @table_subsection{siesta_restart_vmec_arrays_sec, VMEC arrays.}
!>        @item{lmns_m_n_r_, Lambda sine parity.,   vmec_info::lmns_i}
!>     @table_subsection{siesta_restart_vmec_arrays_asym_sec, VMEC arrays asym.}
!>        @item{lmnc_m_n_r_, Lambda cosine parity., vmec_info::lmnc_i}
!>     @table_subsection{siesta_restart_mag_arrays_sec, Magnetic fields.}
!>        @item{bsubsmnsh_m_n_r_, B_s component sine parity.,   }
!>        @item{bsubumnch_m_n_r_, B_u component cosine parity., }
!>        @item{bsubvmnch_m_n_r_, B_v component cosine parity., }
!>        @item{bsupsmnsh_m_n_r_, B^s component sine parity.,   }
!>        @item{bsupumnch_m_n_r_, B^u component cosine parity., }
!>        @item{bsupvmnch_m_n_r_, B^v component cosine parity., }
!>     @table_subsection{siesta_restart_mag_arrays_asym_sec, Magnetic fields.}
!>        @item{bsubsmnch_m_n_r_, B_s component cosine parity., }
!>        @item{bsubumnsh_m_n_r_, B_u component sine parity.,   }
!>        @item{bsubvmnsh_m_n_r_, B_v component sine parity.,   }
!>        @item{bsupsmnch_m_n_r_, B^s component cosine parity., }
!>        @item{bsupumnsh_m_n_r_, B^u component sine parity.,   }
!>        @item{bsupvmnsh_m_n_r_, B^v component sine parity.,   }
!>     @table_subsection{siesta_restart_pres_arrays_sec, Pressure.}
!>        @item{pmnch_m_n_r_, Pressure cosine parity., }
!>     @table_subsection{siesta_restart_mag_arrays_asym_sec, Magnetic fields.}
!>        @item{pmnsh_m_n_r_, Pressure sine parity., }
!>     @table_subsection{siesta_restart_mag_arrays_sec, Magnetic fields.}
!>        @item{fsubsmnsh_m_n_r_, F_s component sine parity.,   quantities::fsubsmnsf}
!>        @item{fsubumnch_m_n_r_, F_u component cosine parity., quantities::fsubumncf}
!>        @item{fsubvmnch_m_n_r_, F_v component cosine parity., quantities::fsubvmncf}
!>        @item{fsupsmnsh_m_n_r_, F^s component sine parity.,   quantities::fsupsmncf}
!>        @item{fsupumnch_m_n_r_, F^u component cosine parity., quantities::fsupumnsf}
!>        @item{fsupvmnch_m_n_r_, F^v component cosine parity., quantities::fsupvmnsf}
!>     @table_subsection{siesta_restart_mag_arrays_asym_sec, Magnetic fields.}
!>        @item{fsubsmnch_m_n_r_, F_s component cosine parity., quantities::fsubsmncf}
!>        @item{fsubumnsh_m_n_r_, F_u component sine parity.,   quantities::fsubumnsf}
!>        @item{fsubvmnsh_m_n_r_, F_v component sine parity.,   quantities::fsubvmnsf}
!>        @item{fsupsmnch_m_n_r_, F^s component cosine parity., quantities::fsupsmncf}
!>        @item{bsupumnsh_m_n_r_, F^u component sine parity.,   quantities::fsupumnsf}
!>        @item{bsupvmnsh_m_n_r_, F^v component sine parity.,   quantities::fsupvmnsf}
!>  @end_table
!>
!>  @section siesta_restart_prog_ref_sec Programmers Reference
!>  Reference material for the coding to implement this restart file is found in
!>  the @ref restart_mod module.
!-------------------------------------------------------------------------------
!*******************************************************************************
!>  @file restart_mod.f90
!>  @brief Contains module @ref restart_mod.
!
!  Note separating the Doxygen comment block here so detailed decription is
!  found in the Module not the file.
!
!>  Contains routines for writting the restart file.
!*******************************************************************************
      MODULE restart_mod
      USE ezcdf
      USE stel_kinds
      USE metrics, ONLY: tolowerh
      USE descriptor_mod, ONLY: iam
      USE shared_data, ONLY: lasym, unit_out, wtotal0
      USE v3_utilities, ONLY: assert_eq
      USE stel_constants, ONLY: mu0

      IMPLICIT NONE

!*******************************************************************************
!  restart module parameters
!*******************************************************************************
!>  Version number.
      INTEGER, PARAMETER           :: restart_version = 0

!  Flag will be placed in the last bits. This way the version number and flags 
!  can fit in the same memory.
!>  Bit position for the lasym flag.
      INTEGER, PARAMETER           :: restart_lasym = 31

!>  Radial Dimension names.
      CHARACTER (LEN=*), DIMENSION(1), PARAMETER ::                            &
     &   radial_dim = (/ 'radius' /)
!>  Fourier Dimension names.
      CHARACTER (LEN=*), DIMENSION(3), PARAMETER ::                            &
     &   restart_dims = (/ 'm-mode', 'n-mode', radial_dim(1) /)

!>  Name for the restart file number of radial points.
      CHARACTER (len=*), PARAMETER :: vn_nsin = 'nrad'
!>  Name for the restart file number of poloidal modes.
      CHARACTER (len=*), PARAMETER :: vn_mpolin = 'mpol'
!>  Name for the restart file number of toroidal modes.
      CHARACTER (len=*), PARAMETER :: vn_ntorin = 'ntor'
!>  Name for the restart file number of field periods.
      CHARACTER (len=*), PARAMETER :: vn_nfpin = 'nfp'
!>  Name for the restart file number of wout file modes.
      CHARACTER (len=*), PARAMETER :: vn_wout = 'wout_file'
!>  Name for the restart file number of state flags modes.
      CHARACTER (len=*), PARAMETER :: vn_flags = 'state_flags'
!>  Name for the restart file toroidal modes array.
      CHARACTER (len=*), PARAMETER :: vn_tor_modes = 'tor_modes'

!>  Name for the restart file jbsupss.
      CHARACTER (len=*), PARAMETER :: vn_jbsupss = 'JBsupssh_m_n_r_'
!>  Name for the restart file jbsupuc.
      CHARACTER (len=*), PARAMETER :: vn_jbsupuc = 'JBsupuch_m_n_r_'
!>  Name for the restart file jbsupvc.
      CHARACTER (len=*), PARAMETER :: vn_jbsupvc = 'JBsupvch_m_n_r_'
!>  Name for the restart file jbsupsc.
      CHARACTER (len=*), PARAMETER :: vn_jbsupsc = 'JBsupsch_m_n_r_'
!>  Name for the restart file jbsupus.
      CHARACTER (len=*), PARAMETER :: vn_jbsupus = 'JBsupush_m_n_r_'
!>  Name for the restart file jbsupus.
      CHARACTER (len=*), PARAMETER :: vn_jbsupvs = 'JBsupvsh_m_n_r_'
!>  Name for the restart file jpresc.
      CHARACTER (len=*), PARAMETER :: vn_jpresc = 'jpresch_m_n_r_'
!>  Name for the restart file jpress.
      CHARACTER (len=*), PARAMETER :: vn_jpress = 'jpressh_m_n_r_'

!>  Name for the restart file rmnc.
      CHARACTER (len=*), PARAMETER :: vn_rmnc = 'rmnc_m_n_r_'
!>  Name for the restart file rmns.
      CHARACTER (len=*), PARAMETER :: vn_rmns = 'rmns_m_n_r_'
!>  Name for the restart file zmnc.
      CHARACTER (len=*), PARAMETER :: vn_zmnc = 'zmnc_m_n_r_'
!>  Name for the restart file zmns.
      CHARACTER (len=*), PARAMETER :: vn_zmns = 'zmns_m_n_r_'

!>  Name for the restart file chipf.
      CHARACTER (len=*), PARAMETER :: vn_chipf = 'chipf(r)'
!>  Name for the restart file phipf.
      CHARACTER (len=*), PARAMETER :: vn_phipf = 'phipf(r)'

!>  Name for the restart file bsupsmns.
      CHARACTER (len=*), PARAMETER :: vn_bsupsmns = 'bsupsmnsh_m_n_r_'
!>  Name for the restart file bsupsmnc.
      CHARACTER (len=*), PARAMETER :: vn_bsupsmnc = 'bsupsmnch_m_n_r_'
!>  Name for the restart file bsupumns.
      CHARACTER (len=*), PARAMETER :: vn_bsupumns = 'bsupumnsh_m_n_r_'
!>  Name for the restart file bsupumnc.
      CHARACTER (len=*), PARAMETER :: vn_bsupumnc = 'bsupumnch_m_n_r_'
!>  Name for the restart file bsupvmns.
      CHARACTER (len=*), PARAMETER :: vn_bsupvmns = 'bsupvmnsh_m_n_r_'
!>  Name for the restart file bsupvmnc.
      CHARACTER (len=*), PARAMETER :: vn_bsupvmnc = 'bsupvmnch_m_n_r_'
!>  Name for the restart file bsubsmns.
      CHARACTER (len=*), PARAMETER :: vn_bsubsmns = 'bsubsmnsh_m_n_r_'
!>  Name for the restart file bsubsmnc.
      CHARACTER (len=*), PARAMETER :: vn_bsubsmnc = 'bsubsmnch_m_n_r_'
!>  Name for the restart file bsubumns.
      CHARACTER (len=*), PARAMETER :: vn_bsubumns = 'bsubumnsh_m_n_r_'
!>  Name for the restart file bsubumnc.
      CHARACTER (len=*), PARAMETER :: vn_bsubumnc = 'bsubumnch_m_n_r_'
!>  Name for the restart file bsubvmns.
      CHARACTER (len=*), PARAMETER :: vn_bsubvmns = 'bsubvmnsh_m_n_r_'
!>  Name for the restart file bsubvmnc.
      CHARACTER (len=*), PARAMETER :: vn_bsubvmnc = 'bsubvmnch_m_n_r_'
!>  Name for the restart file pmns.
      CHARACTER (len=*), PARAMETER :: vn_pmns = 'pmnsh_m_n_r_'
!>  Name for the restart file pmnc.
      CHARACTER (len=*), PARAMETER :: vn_pmnc = 'pmnch_m_n_r_'
!>  Name for the restart file jksupsmns.
      CHARACTER (len=*), PARAMETER :: vn_jksupsmns = 'jksupsmnsf_m_n_r_'
!>  Name for the restart file jksupsmnc.
      CHARACTER (len=*), PARAMETER :: vn_jksupsmnc = 'jksupsmncf_m_n_r_'
!>  Name for the restart file jksupumns.
      CHARACTER (len=*), PARAMETER :: vn_jksupumns = 'jksupumnsf_m_n_r_'
!>  Name for the restart file jksupumnc.
      CHARACTER (len=*), PARAMETER :: vn_jksupumnc = 'jksupumncf_m_n_r_'
!>  Name for the restart file jksupvmns.
      CHARACTER (len=*), PARAMETER :: vn_jksupvmns = 'jksupvmnsf_m_n_r_'
!>  Name for the restart file jksupvmnc.
      CHARACTER (len=*), PARAMETER :: vn_jksupvmnc = 'jksupvmncf_m_n_r_'
!>  Name for the restart file fsupsmns.
      CHARACTER (len=*), PARAMETER :: vn_fsupsmns = 'fsupsmnsf_m_n_r_'
!>  Name for the restart file fsupsmnc.
      CHARACTER (len=*), PARAMETER :: vn_fsupsmnc = 'fsupsmncf_m_n_r_'
!>  Name for the restart file fsupumns.
      CHARACTER (len=*), PARAMETER :: vn_fsupumns = 'fsupumnsf_m_n_r_'
!>  Name for the restart file fsupumnc.
      CHARACTER (len=*), PARAMETER :: vn_fsupumnc = 'fsupumncf_m_n_r_'
!>  Name for the restart file fsupvmns.
      CHARACTER (len=*), PARAMETER :: vn_fsupvmns = 'fsupvmnsf_m_n_r_'
!>  Name for the restart file fsupvmnc.
      CHARACTER (len=*), PARAMETER :: vn_fsupvmnc = 'fsupvmncf_m_n_r_'
!>  Name for the restart file fsubsmns.
      CHARACTER (len=*), PARAMETER :: vn_fsubsmns = 'fsubsmnsf_m_n_r_'
!>  Name for the restart file fsubsmnc.
      CHARACTER (len=*), PARAMETER :: vn_fsubsmnc = 'fsubsmncf_m_n_r_'
!>  Name for the restart file fsubumns.
      CHARACTER (len=*), PARAMETER :: vn_fsubumns = 'fsubumnsf_m_n_r_'
!>  Name for the restart file fsubumnc.
      CHARACTER (len=*), PARAMETER :: vn_fsubumnc = 'fsubumncf_m_n_r_'
!>  Name for the restart file fsubvmns.
      CHARACTER (len=*), PARAMETER :: vn_fsubvmns = 'fsubvmnsf_m_n_r_'
!>  Name for the restart file fsubvmnc.
      CHARACTER (len=*), PARAMETER :: vn_fsubvmnc = 'fsubvmncf_m_n_r_'

!>  Name for the restart file lmns.
      CHARACTER (len=*), PARAMETER :: vn_lmns = 'lmns_m_n_r_'
!>  Name for the restart file lmnc.
      CHARACTER (len=*), PARAMETER :: vn_lmnc = 'lmnc_m_n_r_'

!>  Name for the restart file p_factor.
      CHARACTER (len=*), PARAMETER :: vn_p_factor = 'p_factor'
!>  Name for the restart file b_factor.
      CHARACTER (len=*), PARAMETER :: vn_b_factor = 'b_factor'
!>  Name for the restart file wb.
      CHARACTER (len=*), PARAMETER :: vn_wb = 'wb'
!>  Name for the restart file wb.
      CHARACTER (len=*), PARAMETER :: vn_wp = 'wp'
!>  Name for the restart file rmajor.
      CHARACTER (len=*), PARAMETER :: vn_rmajor = 'rmajor'
!>  Name for the restart file total toroidal current.
      CHARACTER (len=*), PARAMETER :: vn_curtor = 'curtor'

!>  Name for the restart file p_max.
      CHARACTER (len=*), PARAMETER :: vn_p_max = 'p_max'
!>  Name for the restart file p_min.
      CHARACTER (len=*), PARAMETER :: vn_p_min = 'p_min'

!>  Inital stored energy.
      CHARACTER (len=*), PARAMETER :: vn_wtotal0 = 'wtotal0'

      CONTAINS

!*******************************************************************************
!  UTILITY SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  @brief Reads the restart file.
!>
!>  Reads the restart information and initalizes SIESTA quantities.
!>
!>  @param[in]    restart_ext Restart extention.
!>  @param[inout] wout_file   Name of the wout file.
!>  @param[in]    mpolin      Namelist number of polodal modes.
!>  @param[in]    ntorin      Namelist number of toroidal modes.
!>  @param[in]    nfpin       Namelist number of field periods.
!>  @param[in]    nsin        Namelist number of radial grid points.
!>  @param[in]    tor_modesin Namelist torodial modes.
!>  @returns Preconditioner control flag.
!-------------------------------------------------------------------------------
      FUNCTION restart_read(restart_ext, wout_file, mpolin, ntorin, nfpin,     &
                            nsin, tor_modesin)
      USE quantities, ONLY: jbsupsmnsh, jbsupsmnch,                            &
                            jbsupumnsh, jbsupumnch,                            &
                            jbsupvmnsh, jbsupvmnch,                            &
                            jpmnsh,     jpmnch,                                &
                            b_factor, p_factor, alloc_quantities,              &
                            dealloc_quantities
      USE island_params, ONLY: chipf => chipf_i, phipf => phipf_i,             &
     &                         wb => wb_i, wp => wp_i, nfp_i, gamma,           &
     &                         gnorm => gnorm_i, rmajor => rmajor_i,           &
     &                         fourier_context, nu_i, nv_i
      USE vmec_info, ONLY: rmnc => rmnc_i, zmns => zmns_i,                     &
     &                     rmns => rmns_i, zmnc => zmnc_i,                     &
     &                     vmec_info_construct_island,                         &
     &                     vmec_info_destruct_island
      USE metrics, ONLY: LoadGrid
      USE fourier, ONLY: m0, m1, fourier_class
      USE stel_constants, ONLY: one
      USE siesta_namelist, ONLY: lrestart

      IMPLICIT NONE

!  Declare Arguments
      INTEGER                                        :: restart_read
      CHARACTER (len=*), INTENT(in)                  :: restart_ext
      CHARACTER (len=*), INTENT(inout)               :: wout_file
      INTEGER, INTENT(in)                            :: mpolin
      INTEGER, INTENT(in)                            :: ntorin
      INTEGER, INTENT(in)                            :: nsin
      INTEGER, INTENT(in)                            :: nfpin
      INTEGER, DIMENSION(-ntorin:ntorin), INTENT(in) :: tor_modesin

!  local variables
      INTEGER                                        :: flags
      INTEGER                                        :: ncid
      INTEGER                                        :: ns
      INTEGER                                        :: mpol
      INTEGER                                        :: ntor
      INTEGER                                        :: nfp
      REAL (dp), DIMENSION(:,:,:), ALLOCATABLE       :: tempmn_r
      REAL (dp), DIMENSION(:), ALLOCATABLE           :: temp_r
      INTEGER, DIMENSION(:), ALLOCATABLE             :: temp_modes
      INTEGER                                        :: nmin
      INTEGER                                        :: status
      CHARACTER (LEN=256)                            :: filename

!  Start of executable code
      restart_read = 1

      CALL dealloc_quantities
      CALL alloc_quantities

      filename = 'siesta_' // TRIM(restart_ext) // '.nc'
      CALL cdf_open(ncid, TRIM(filename), 'r', status)
      CALL assert_eq(0, status, 'Failed to read restart file: ' //             &
     &                          TRIM(filename))

      CALL cdf_read(ncid, vn_flags, flags)

!  The namelist inut file can change the size of the radial grid and the 
!  poloidal and toroidal modes.
      CALL cdf_read(ncid, vn_nsin, ns)
      CALL cdf_read(ncid, vn_mpolin, mpol)
      CALL cdf_read(ncid, vn_ntorin, ntor)
      CALL cdf_read(ncid, vn_nfpin, nfp)

      CALL cdf_read(ncid, vn_wtotal0, wtotal0)

      ALLOCATE(temp_modes(-ntor:ntor))
      CALL cdf_read(ncid, vn_tor_modes, temp_modes)

      IF (nfpin .lt. 1) THEN
         nfp_i = nfp
      ELSE
         nfp_i = nfpin
      END IF

      fourier_context => fourier_class(mpolin, ntorin, nu_i, nv_i, nfp_i,      &
     &                                 lasym, tor_modesin)

      nmin = MIN(ntor, ntorin)

      CALL cdf_read(ncid, vn_wout, wout_file)

!  The namelist input file may turn the asymmetric terms on and off.
      ALLOCATE(tempmn_r(0:mpol,-ntor:ntor,ns))
      ALLOCATE(temp_r(ns))
      CALL vmec_info_destruct_island
      CALL vmec_info_construct_island(mpolin, ntorin, nsin, lasym)

      IF (.not.ALLOCATED(chipf)) THEN
         ALLOCATE(chipf(nsin))
      END IF
      CALL cdf_read(ncid, vn_chipf, temp_r)
      CALL interpit_1d(temp_r, chipf, ns, nsin, .false., 1)

      IF (.not.ALLOCATED(phipf)) THEN
         ALLOCATE(phipf(nsin))
      END IF
      CALL cdf_read(ncid, vn_phipf, temp_r)
      CALL interpit_1d(temp_r, phipf, ns, nsin, .false., 1)

      CALL cdf_read(ncid, vn_jbsupss, tempmn_r)
      jbsupsmnsh(:,:,1) = 0
      CALL interpit(tempmn_r, jbsupsmnsh, ns, nsin, mpol, mpolin,              &
     &              ntor, ntorin, nfp, nfp_i, temp_modes, tor_modesin,         &
     &              .true.)

      CALL cdf_read(ncid, vn_jbsupuc, tempmn_r)
      jbsupumnch(:,:,1) = 0
      CALL interpit(tempmn_r, jbsupumnch, ns, nsin, mpol, mpolin,              &
     &              ntor, ntorin, nfp, nfp_i, temp_modes, tor_modesin,         &
     &              .true.)

      CALL cdf_read(ncid, vn_jbsupvc, tempmn_r)
      jbsupvmnch(:,:,1) = 0
      CALL interpit(tempmn_r, jbsupvmnch, ns, nsin, mpol, mpolin,              &
     &              ntor, ntorin, nfp, nfp_i, temp_modes, tor_modesin,         &
     &              .true.)

      CALL cdf_read(ncid, vn_jpresc,  tempmn_r)
      jpmnch(:,:,1) = 0
      CALL interpit(tempmn_r, jpmnch,     ns, nsin, mpol, mpolin,              &
     &              ntor, ntorin, nfp, nfp_i, temp_modes, tor_modesin,         &
     &              .true.)

      CALL cdf_read(ncid, vn_rmnc, tempmn_r)
      CALL interpit(tempmn_r, rmnc, ns, nsin, mpol, mpolin,                    &
     &              ntor, ntorin, nfp, nfp_i, temp_modes, tor_modesin,         &
     &              .false.)

      CALL cdf_read(ncid, vn_zmns, tempmn_r)
      CALL interpit(tempmn_r, zmns, ns, nsin, mpol, mpolin,                    &
     &              ntor, ntorin, nfp, nfp_i, temp_modes, tor_modesin,         &
     &              .false.)

      IF (BTEST(flags, restart_lasym) .and. lasym) THEN
         jbsupsmnch(:,:,1) = 0
         CALL cdf_read(ncid, vn_jbsupsc, tempmn_r)
         CALL interpit(tempmn_r, jbsupsmnch, ns, nsin, mpol, mpolin,           &
     &                 ntor, ntorin, nfp, nfp_i, temp_modes, tor_modesin,      &
     &                 .true.)

         CALL cdf_read(ncid, vn_jbsupus, tempmn_r)
         CALL interpit(tempmn_r, jbsupumnsh, ns, nsin, mpol, mpolin,           &
     &                 ntor, ntorin, nfp, nfp_i, temp_modes, tor_modesin,      &
     &                 .true.)

         CALL cdf_read(ncid, vn_jbsupvs, tempmn_r)
         CALL interpit(tempmn_r, jbsupvmnsh, ns, nsin, mpol, mpolin,           &
     &                 ntor, ntorin, nfp, nfp_i, temp_modes, tor_modesin,      &
     &                 .true.)

         CALL cdf_read(ncid, vn_jpress,  tempmn_r)
         CALL interpit(tempmn_r, jpmnsh,     ns, nsin, mpol, mpolin,           &
     &                 ntor, ntorin, nfp, nfp_i, temp_modes, tor_modesin,      &
     &                 .true.)

         CALL cdf_read(ncid, vn_rmns, tempmn_r)
         CALL interpit(tempmn_r, rmns, ns, nsin, mpol, mpolin,                 &
     &                 ntor, ntorin, nfp, nfp_i, temp_modes, tor_modesin,      &
     &                 .false.)

         CALL cdf_read(ncid, vn_zmns, tempmn_r)
         CALL interpit(tempmn_r, zmnc, ns, nsin, mpol, mpolin,                 &
     &                 ntor, ntorin, nfp, nfp_i, temp_modes, tor_modesin,      &
     &                 .false.)

      ELSE IF (lasym) THEN
         jbsupsmnch = 0
         jbsupumnsh = 0
         jbsupvmnsh = 0
         jpmnsh = 0
         rmns = 0
         rmnc = 0
      END IF

!  Read normalization factors.
      CALL cdf_read(ncid, vn_wb, wb)
      CALL cdf_read(ncid, vn_wp, wp)
      CALL cdf_read(ncid, vn_rmajor, rmajor)

      CALL cdf_close(ncid)

      CALL restart_normalize(rmnc, one)
      CALL restart_normalize(zmns, one)
      IF (lasym) THEN
         CALL restart_normalize(rmns, one)
         CALL restart_normalize(zmnc, one)
      END IF

      CALL LoadGrid(status)

      DEALLOCATE(tempmn_r)
      DEALLOCATE(temp_r)
      DEALLOCATE(temp_modes)

!  Init quantities.
      p_factor = 1.0/ABS(wb)
      b_factor = SQRT(p_factor)
      gnorm = ABS(wb)/(wb + wp/(gamma - 1.0))

      WRITE (unit_out, 1000) mpol, ntor, ns

      restart_read = 2

1000  FORMAT(/,' RESTARTED FROM RUN PARAMETERS M: ',i3,' N: ',i3,' NS: ', i3,/)

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Write the restart file.
!>
!>  Writes the restart information.
!>
!>  @param[in] restart_ext Restart file extension.
!>  @param[in] wout_file   Name of the wout file.
!-------------------------------------------------------------------------------
      SUBROUTINE restart_write(restart_ext, wout_file)
      USE quantities, ONLY: jbsupsmnsh, jbsupsmnch,                            &
                            jbsupumnsh, jbsupumnch,                            &
                            jbsupvmnsh, jbsupvmnch,                            &
                            jpmnsh,     jpmnch,                                &
                            fsubsmncf, fsubumnsf, fsubvmnsf,                   &
                            fsubsmnsf, fsubumncf, fsubvmncf,                   &
                            fsupsmncf, fsupumnsf, fsupvmnsf,                   &
                            fsupsmnsf, fsupumncf, fsupvmncf,                   &
                            b_factor,   p_factor, jacobh
      USE fourier, ONLY: f_cos, f_sin, f_sum, f_none, n0, m0
      USE island_params, ONLY: nfp => nfp_i, chipf => chipf_i,                 &
                               phipf => phipf_i, wb => wb_i, wp => wp_i,       &
                               rmajor => rmajor_i, fourier_context,            &
                               mpol=>mpol_i, ntor=>ntor_i, ns=>ns_i,           &
                               ntheta=>nu_i, nzeta=>nv_i
      USE utilities, ONLY: curl_htof
      USE stel_constants, ONLY: one
      USE vmec_info, ONLY: rmnc => rmnc_i, zmns => zmns_i,                     &
     &                     rmns => rmns_i, zmnc => zmnc_i,                     &
     &                     lmns => lmns_i, lmnc => lmnc_i

      IMPLICIT NONE

!  Declare Arguments
      CHARACTER (len=*), INTENT(in)            :: restart_ext
      CHARACTER (len=*), INTENT(in)            :: wout_file

!  Local Variables
      INTEGER                                  :: flags
      INTEGER                                  :: ncid
      INTEGER                                  :: nblock
      INTEGER                                  :: s
      INTEGER                                  :: m
      INTEGER                                  :: n
      INTEGER                                  :: status

      REAL (dp)                                :: r0
      REAL (dp)                                :: tempscalar
      REAL (dp), DIMENSION(:,:,:), ALLOCATABLE :: bsupsijh
      REAL (dp), DIMENSION(:,:,:), ALLOCATABLE :: bsupuijh
      REAL (dp), DIMENSION(:,:,:), ALLOCATABLE :: bsupvijh
      REAL (dp), DIMENSION(:,:,:), ALLOCATABLE :: bsubsijh
      REAL (dp), DIMENSION(:,:,:), ALLOCATABLE :: bsubuijh
      REAL (dp), DIMENSION(:,:,:), ALLOCATABLE :: bsubvijh
      REAL (dp), DIMENSION(:,:,:), ALLOCATABLE :: bsubsmnh
      REAL (dp), DIMENSION(:,:,:), ALLOCATABLE :: bsubumnh
      REAL (dp), DIMENSION(:,:,:), ALLOCATABLE :: bsubvmnh
      REAL (dp), DIMENSION(:,:,:), ALLOCATABLE :: jksupsmnf
      REAL (dp), DIMENSION(:,:,:), ALLOCATABLE :: jksupumnf
      REAL (dp), DIMENSION(:,:,:), ALLOCATABLE :: jksupvmnf
      REAL (dp), DIMENSION(:,:,:), ALLOCATABLE :: tempmn_w
      REAL (dp), DIMENSION(:,:,:), ALLOCATABLE :: pijh

!  Start of executable code.
      CALL cdf_open(ncid, 'siesta_' // TRIM(restart_ext) // '.nc', 'w', status)
      CALL assert_eq(0, status, 'Failed to write restart file.')

      flags = restart_version
      IF (lasym) THEN
         flags = IBSET(flags, restart_lasym)
      END IF

      ALLOCATE(tempmn_w(0:mpol,-ntor:ntor,ns))

      CALL cdf_define(ncid, vn_flags, flags)
      CALL cdf_define(ncid, vn_nsin, ns)
      CALL cdf_define(ncid, vn_mpolin, mpol)
      CALL cdf_define(ncid, vn_ntorin, ntor)
      CALL cdf_define(ncid, vn_nfpin, nfp)
      CALL cdf_define(ncid, vn_wout, wout_file)

      CALL cdf_define(ncid, vn_wtotal0, wtotal0)

      CALL cdf_define(ncid, vn_tor_modes, fourier_context%tor_modes)

      CALL cdf_define(ncid, vn_jbsupss, jbsupsmnsh, dimname=restart_dims)
      CALL cdf_define(ncid, vn_jbsupuc, jbsupumnch, dimname=restart_dims)
      CALL cdf_define(ncid, vn_jbsupvc, jbsupvmnch, dimname=restart_dims)
      CALL cdf_define(ncid, vn_jpresc,  jpmnch,     dimname=restart_dims)

      CALL cdf_define(ncid, vn_rmnc, tempmn_w, dimname=restart_dims)
      CALL cdf_define(ncid, vn_zmns, tempmn_w, dimname=restart_dims)

      CALL cdf_define(ncid, vn_chipf, chipf, dimname=radial_dim)
      CALL cdf_define(ncid, vn_phipf, phipf, dimname=radial_dim)

      CALL cdf_define(ncid, vn_wb, wb)
      CALL cdf_define(ncid, vn_wp, wp)

!  Using variable with jacobian only for the dimension sizes and types. Before
!  writting, the jacobian normalization will be removed.
      CALL cdf_define(ncid, vn_bsupsmns,  jbsupsmnsh, dimname=restart_dims)
      CALL cdf_define(ncid, vn_bsupumnc,  jbsupumnch, dimname=restart_dims)
      CALL cdf_define(ncid, vn_bsupvmnc,  jbsupvmnch, dimname=restart_dims)
      CALL cdf_define(ncid, vn_bsubsmns,  jbsupsmnsh, dimname=restart_dims)
      CALL cdf_define(ncid, vn_bsubumnc,  jbsupumnch, dimname=restart_dims)
      CALL cdf_define(ncid, vn_bsubvmnc,  jbsupvmnch, dimname=restart_dims)
      CALL cdf_define(ncid, vn_jksupsmns, jbsupsmnsh, dimname=restart_dims)
      CALL cdf_define(ncid, vn_jksupumnc, jbsupumnch, dimname=restart_dims)
      CALL cdf_define(ncid, vn_jksupvmnc, jbsupvmnch, dimname=restart_dims)
      CALL cdf_define(ncid, vn_fsupsmnc,  jbsupsmnsh, dimname=restart_dims)
      CALL cdf_define(ncid, vn_fsupumns,  jbsupumnch, dimname=restart_dims)
      CALL cdf_define(ncid, vn_fsupvmns,  jbsupvmnch, dimname=restart_dims)
      CALL cdf_define(ncid, vn_fsubsmnc,  jbsupsmnsh, dimname=restart_dims)
      CALL cdf_define(ncid, vn_fsubumns,  jbsupumnch, dimname=restart_dims)
      CALL cdf_define(ncid, vn_fsubvmns,  jbsupvmnch, dimname=restart_dims)
      CALL cdf_define(ncid, vn_lmns,      lmns,       dimname=restart_dims)
      CALL cdf_define(ncid, vn_pmnc,      jpmnch,     dimname=restart_dims)
      CALL cdf_define(ncid, vn_b_factor,  b_factor)
      CALL cdf_define(ncid, vn_p_factor,  p_factor)
      CALL cdf_define(ncid, vn_rmajor,    rmajor)
      CALL cdf_define(ncid, vn_p_max,     b_factor)
      CALL cdf_define(ncid, vn_p_min,     p_factor)
      CALL cdf_define(ncid, vn_curtor,    tempscalar)

      IF (lasym) THEN
         CALL cdf_define(ncid, vn_jbsupsc, jbsupsmnch, dimname=restart_dims)
         CALL cdf_define(ncid, vn_jbsupus, jbsupumnsh, dimname=restart_dims)
         CALL cdf_define(ncid, vn_jbsupvs, jbsupvmnsh, dimname=restart_dims)
         CALL cdf_define(ncid, vn_jpress,  jpmnsh,     dimname=restart_dims)

         CALL cdf_define(ncid, vn_rmns, tempmn_w, dimname=restart_dims)
         CALL cdf_define(ncid, vn_zmnc, tempmn_w, dimname=restart_dims)

         CALL cdf_define(ncid, vn_bsupsmnc,  jbsupsmnch, dimname=restart_dims)
         CALL cdf_define(ncid, vn_bsupumns,  jbsupumnsh, dimname=restart_dims)
         CALL cdf_define(ncid, vn_bsupvmns,  jbsupvmnsh, dimname=restart_dims)
         CALL cdf_define(ncid, vn_bsubsmnc,  jbsupsmnch, dimname=restart_dims)
         CALL cdf_define(ncid, vn_bsubumns,  jbsupumnsh, dimname=restart_dims)
         CALL cdf_define(ncid, vn_bsubvmns,  jbsupvmnsh, dimname=restart_dims)
         CALL cdf_define(ncid, vn_jksupsmnc, jbsupsmnch, dimname=restart_dims)
         CALL cdf_define(ncid, vn_jksupumns, jbsupumnsh, dimname=restart_dims)
         CALL cdf_define(ncid, vn_jksupvmns, jbsupvmnsh, dimname=restart_dims)
         CALL cdf_define(ncid, vn_fsupsmns,  jbsupsmnsh, dimname=restart_dims)
         CALL cdf_define(ncid, vn_fsupumnc,  jbsupumnch, dimname=restart_dims)
         CALL cdf_define(ncid, vn_fsupvmnc,  jbsupvmnch, dimname=restart_dims)
         CALL cdf_define(ncid, vn_fsubsmns,  jbsupsmnsh, dimname=restart_dims)
         CALL cdf_define(ncid, vn_fsubumnc,  jbsupumnch, dimname=restart_dims)
         CALL cdf_define(ncid, vn_fsubvmnc,  jbsupvmnch, dimname=restart_dims)
         CALL cdf_define(ncid, vn_lmnc,      lmnc,       dimname=restart_dims)
         CALL cdf_define(ncid, vn_pmns, jpmnsh, dimname=restart_dims)
      END IF

      CALL cdf_write(ncid, vn_flags, flags)
      CALL cdf_write(ncid, vn_nsin, ns)
      CALL cdf_write(ncid, vn_mpolin, mpol)
      CALL cdf_write(ncid, vn_ntorin, ntor)
      CALL cdf_write(ncid, vn_nfpin, nfp)
      CALL cdf_write(ncid, vn_wout, wout_file)

      CALL cdf_write(ncid, vn_wtotal0, wtotal0)

      CALL cdf_write(ncid, vn_tor_modes, fourier_context%tor_modes)

      CALL cdf_write(ncid, vn_jbsupss, jbsupsmnsh)
      CALL cdf_write(ncid, vn_jbsupuc, jbsupumnch)
      CALL cdf_write(ncid, vn_jbsupvc, jbsupvmnch)
      CALL cdf_write(ncid, vn_jpresc,  jpmnch)

      CALL cdf_write(ncid, vn_chipf, chipf)
      CALL cdf_write(ncid, vn_phipf, phipf)

      tempmn_w = rmnc
      CALL restart_denormalize(tempmn_w, one)
      CALL cdf_write(ncid, vn_rmnc, tempmn_w)
      tempmn_w = zmns
      CALL restart_denormalize(tempmn_w, one)
      CALL cdf_write(ncid, vn_zmns, tempmn_w)

      IF (lasym) THEN
         CALL cdf_write(ncid, vn_jbsupsc, jbsupsmnch)
         CALL cdf_write(ncid, vn_jbsupus, jbsupumnsh)
         CALL cdf_write(ncid, vn_jbsupvs, jbsupvmnsh)
         CALL cdf_write(ncid, vn_jpress,  jpmnsh)

         tempmn_w = rmns
         CALL restart_denormalize(tempmn_w, one)
         CALL cdf_write(ncid, vn_rmns, tempmn_w)
         tempmn_w = zmnc
         CALL restart_denormalize(tempmn_w, one)
         CALL cdf_write(ncid, vn_zmnc, tempmn_w)
      END IF

      CALL cdf_write(ncid, vn_wb, wb)
      CALL cdf_write(ncid, vn_wp, wp)
      CALL cdf_write(ncid, vn_rmajor, rmajor)

      CALL cdf_write(ncid, vn_p_factor, p_factor)
      CALL cdf_write(ncid, vn_b_factor, b_factor)

      ALLOCATE(bsupsijh(ntheta,nzeta,ns))
      ALLOCATE(bsupuijh(ntheta,nzeta,ns))
      ALLOCATE(bsupvijh(ntheta,nzeta,ns))
      ALLOCATE(pijh(ntheta,nzeta,ns))

      bsupsijh = 0
      bsupuijh = 0
      bsupvijh = 0
      pijh = 0

      CALL fourier_context%toijsp(jbsupsmnsh, bsupsijh, f_none, f_sin)
      CALL fourier_context%toijsp(jbsupumnch, bsupuijh, f_none, f_cos)
      CALL fourier_context%toijsp(jbsupvmnch, bsupvijh, f_none, f_cos)
      CALL fourier_context%toijsp(jpmnch, pijh, f_none, f_cos)

      IF (lasym) THEN
         CALL fourier_context%toijsp(jbsupsmnch, bsupsijh, f_sum, f_cos)
         CALL fourier_context%toijsp(jbsupumnsh, bsupuijh, f_sum, f_sin)
         CALL fourier_context%toijsp(jbsupvmnsh, bsupvijh, f_sum, f_sin)
         CALL fourier_context%toijsp(jpmnsh, pijh, f_sum, f_sin)
      END IF

      bsupsijh(:,:,1) = 0
      bsupuijh(:,:,1) = 0
      bsupvijh(:,:,1) = 0
      pijh(:,:,1) = 0

!  Remove the jacobian term.
      bsupsijh = bsupsijh/jacobh
      bsupuijh = bsupuijh/jacobh
      bsupvijh = bsupvijh/jacobh
      pijh = pijh/jacobh

      tempscalar = MAXVAL(pijh(:,:,2:))
      CALL cdf_write(ncid, vn_p_max, tempscalar/(p_factor*mu0))
      tempscalar = MINVAL(pijh(:,:,2:))
      CALL cdf_write(ncid, vn_p_min, tempscalar/(p_factor*mu0))

      tempmn_w = 0

!  Remove the orthonorm and b_factor so these quantities can be directly summed.
      CALL fourier_context%tomnsp(bsupsijh, tempmn_w, f_sin)
      CALL restart_denormalize(tempmn_w, b_factor)
      CALL cdf_write(ncid, vn_bsupsmns, tempmn_w)

      CALL fourier_context%tomnsp(bsupuijh, tempmn_w, f_cos)
      CALL restart_denormalize(tempmn_w, b_factor)
      CALL cdf_write(ncid, vn_bsupumnc, tempmn_w)

      CALL fourier_context%tomnsp(bsupvijh, tempmn_w, f_cos)
      CALL restart_denormalize(tempmn_w, b_factor)
      CALL cdf_write(ncid, vn_bsupvmnc, tempmn_w)

      tempmn_w = fsubsmncf
      CALL restart_denormalize(tempmn_w, 1.0_dp)
      CALL cdf_write(ncid, vn_fsubsmnc, tempmn_w)

      tempmn_w = fsubumnsf
      CALL restart_denormalize(tempmn_w, 1.0_dp)
      CALL cdf_write(ncid, vn_fsubumns, tempmn_w)

      tempmn_w = fsubvmnsf
      CALL restart_denormalize(tempmn_w, 1.0_dp)
      CALL cdf_write(ncid, vn_fsubvmns, tempmn_w)

      tempmn_w = fsupsmncf
      CALL restart_denormalize(tempmn_w, 1.0_dp)
      CALL cdf_write(ncid, vn_fsupsmnc, tempmn_w)

      tempmn_w = fsupumnsf
      CALL restart_denormalize(tempmn_w, 1.0_dp)
      CALL cdf_write(ncid, vn_fsupumns, tempmn_w)

      tempmn_w = fsupvmnsf
      CALL restart_denormalize(tempmn_w, 1.0_dp)
      CALL cdf_write(ncid, vn_fsupvmns, tempmn_w)

      tempmn_w = lmns
      CALL restart_denormalize(tempmn_w, 1.0_dp)
      CALL cdf_write(ncid, vn_lmns, tempmn_w)

!  Remove the orthonorm and p_factor so this quantity can be directly summed.
      CALL fourier_context%tomnsp(pijh, tempmn_w, f_cos)
      CALL restart_denormalize(tempmn_w, p_factor*mu0)
      CALL cdf_write(ncid, vn_pmnc, tempmn_w)

      ALLOCATE(bsubsijh(ntheta,nzeta,ns))
      ALLOCATE(bsubuijh(ntheta,nzeta,ns))
      ALLOCATE(bsubvijh(ntheta,nzeta,ns))

      CALL tolowerh(bsupsijh, bsupuijh, bsupvijh,                              &
                    bsubsijh, bsubuijh, bsubvijh, 1, ns)

!  Need these to compute the currents later.
      ALLOCATE(bsubsmnh(0:mpol,-ntor:ntor,ns))
      ALLOCATE(bsubumnh(0:mpol,-ntor:ntor,ns))
      ALLOCATE(bsubvmnh(0:mpol,-ntor:ntor,ns))

!  Need to compute current before denormalizing quantities.
      CALL fourier_context%tomnsp(bsubsijh, bsubsmnh, f_sin)
      CALL fourier_context%tomnsp(bsubuijh, bsubumnh, f_cos)
      CALL fourier_context%tomnsp(bsubvijh, bsubvmnh, f_cos)

!  Compute currents on the full mesh.
      ALLOCATE(jksupsmnf(0:mpol,-ntor:ntor,ns))
      ALLOCATE(jksupumnf(0:mpol,-ntor:ntor,ns))
      ALLOCATE(jksupvmnf(0:mpol,-ntor:ntor,ns))

      CALL curl_htof(bsubsmnh, bsubumnh, bsubvmnh,                             &
     &               jksupsmnf, jksupumnf, jksupvmnf,                          &
     &               f_sin, 1, ns, ns, tempscalar)

!  Remove the orthonorm and b_factor so these quantities can be directly summed.
      CALL restart_denormalize(bsubsmnh, b_factor)
      CALL restart_denormalize(bsubvmnh, b_factor)
      CALL restart_denormalize(bsubumnh, b_factor)

      CALL cdf_write(ncid, vn_bsubsmns, bsubsmnh)
      CALL cdf_write(ncid, vn_bsubumnc, bsubumnh)
      CALL cdf_write(ncid, vn_bsubvmnc, bsubvmnh)

      CALL restart_denormalize(jksupsmnf, b_factor)
      CALL restart_denormalize(jksupumnf, b_factor)
      CALL restart_denormalize(jksupvmnf, b_factor)

      CALL cdf_write(ncid, vn_jksupsmns, jksupsmnf)
      CALL cdf_write(ncid, vn_jksupumnc, jksupumnf)
      CALL cdf_write(ncid, vn_jksupvmnc, jksupvmnf)

      tempscalar = -fourier_context%orthonorm(m0,n0)*tempscalar/b_factor
      CALL cdf_write(ncid, vn_curtor, tempscalar)

      IF (lasym) THEN
!  Remove the orthonorm and b_factor so these quantities can be directly summed.
         CALL fourier_context%tomnsp(bsupsijh, tempmn_w, f_cos)
         CALL restart_denormalize(tempmn_w, b_factor)
         CALL cdf_write(ncid, vn_bsupsmnc, tempmn_w)

         CALL fourier_context%tomnsp(bsupuijh, tempmn_w, f_sin)
         CALL restart_denormalize(tempmn_w, b_factor)
         CALL cdf_write(ncid, vn_bsupumns, tempmn_w)

         CALL fourier_context%tomnsp(bsupvijh, tempmn_w, f_sin)
         CALL restart_denormalize(tempmn_w, b_factor)
         CALL cdf_write(ncid, vn_bsupvmns, tempmn_w)

         tempmn_w = fsubsmnsf
         CALL restart_denormalize(tempmn_w, 1.0_dp)
         CALL cdf_write(ncid, vn_fsubsmns, tempmn_w)

         tempmn_w = fsubumncf
         CALL restart_denormalize(tempmn_w, 1.0_dp)
         CALL cdf_write(ncid, vn_fsubumnc, tempmn_w)

         tempmn_w = fsubvmncf
         CALL restart_denormalize(tempmn_w, 1.0_dp)
         CALL cdf_write(ncid, vn_fsubvmnc, tempmn_w)

         tempmn_w = fsupsmnsf
         CALL restart_denormalize(tempmn_w, 1.0_dp)
         CALL cdf_write(ncid, vn_fsupsmns, tempmn_w)

         tempmn_w = fsupumncf
         CALL restart_denormalize(tempmn_w, 1.0_dp)
         CALL cdf_write(ncid, vn_fsupumnc, tempmn_w)

         tempmn_w = fsupvmncf
         CALL restart_denormalize(tempmn_w, 1.0_dp)
         CALL cdf_write(ncid, vn_fsupvmnc, tempmn_w)

         tempmn_w = lmnc
         CALL restart_denormalize(tempmn_w, 1.0_dp)
         CALL cdf_write(ncid, vn_lmnc, tempmn_w)

!  Remove the orthonorm and p_factor so this quantity can be directly summed.
         CALL fourier_context%tomnsp(pijh, tempmn_w, f_sin)
         CALL restart_denormalize(tempmn_w, p_factor*mu0)
         CALL cdf_write(ncid, vn_pmns, tempmn_w)

!  Need to compute current before denormalizing quantities.
         CALL fourier_context%tomnsp(bsubsijh, bsubsmnh, f_cos)
         CALL fourier_context%tomnsp(bsubuijh, bsubumnh, f_sin)
         CALL fourier_context%tomnsp(bsubvijh, bsubvmnh, f_sin)

         CALL curl_htof(bsubsmnh, bsubumnh, bsubvmnh,                          &
     &                  jksupsmnf, jksupumnf, jksupvmnf,                       &
     &                  f_cos, 1, ns, ns, tempscalar)

!  Remove the orthonorm and b_factor so these quantities can be directly summed.
         CALL restart_denormalize(bsubsmnh, b_factor)
         CALL restart_denormalize(bsubumnh, b_factor)
         CALL restart_denormalize(bsubvmnh, b_factor)

         CALL cdf_write(ncid, vn_bsubsmnc, bsubsmnh)
         CALL cdf_write(ncid, vn_bsubumns, bsubumnh)
         CALL cdf_write(ncid, vn_bsubvmns, bsubvmnh)

         CALL restart_denormalize(jksupsmnf, b_factor)
         CALL restart_denormalize(jksupumnf, b_factor)
         CALL restart_denormalize(jksupvmnf, b_factor)

         CALL cdf_write(ncid, vn_jksupsmnc, jksupsmnf)
         CALL cdf_write(ncid, vn_jksupumns, jksupumnf)
         CALL cdf_write(ncid, vn_jksupvmns, jksupvmnf)
      END IF

      DEALLOCATE(bsupsijh)
      DEALLOCATE(bsupuijh)
      DEALLOCATE(bsupvijh)
      DEALLOCATE(pijh)

      DEALLOCATE(bsubsmnh)
      DEALLOCATE(bsubumnh)
      DEALLOCATE(bsubvmnh)

      DEALLOCATE(jksupsmnf)
      DEALLOCATE(jksupumnf)
      DEALLOCATE(jksupvmnf)

      DEALLOCATE(tempmn_w)

      CALL cdf_close(ncid)

      END SUBROUTINE

!-------------------------------------------------------------------------------
!>  @brief Interpolate fourier quantites from the restart file.
!>
!>  Restart files can change the number of surfaces and modes.
!>
!>  @param[in]  aold          Value from the restart file.
!>  @param[out] anew          New interpolated value.
!>  @param[in]  ns_old        Radial grid size in the restart file.
!>  @param[in]  ns_new        New radial grid size.
!>  @param[in]  mpol_old      Number of poloidal modes in the restart file.
!>  @param[in]  mpol_new      New number of poloidal modes.
!>  @param[in]  ntor_old      Number of totoidal modes in the restart file.
!>  @param[in]  ntor_new      New number of totoidal modes.
!>  @param[in]  nfp_old       Number of field periods in the restart file.
!>  @param[in]  nfp_new       New number of field periods.
!>  @param[in]  tor_modes_old Toroidal modes in the restart file.
!>  @param[in]  tor_modes_new New toroidal modes.
!>  @param[in]  lhalf         Grid type.
!-------------------------------------------------------------------------------
      SUBROUTINE interpit(aold,     anew,                                      &
     &                    ns_old,   ns_new,                                    &
     &                    mpol_old, mpol_new,                                  &
     &                    ntor_old, ntor_new,                                  &
     &                    nfp_old,  nfp_new,                                   &
     &                    tor_modes_old, tor_modes_new,                        &
     &                    lhalf)
      USE stel_constants, ONLY: one, zero
      USE island_params, ONLY: fourier_context

      IMPLICIT NONE

!  Declare Arguments
      REAL (dp), DIMENSION(:,:,:), INTENT(in)            :: aold
      REAL (dp), DIMENSION(:,:,:), INTENT(out)           :: anew
      INTEGER, INTENT(in)                                :: ns_old
      INTEGER, INTENT(in)                                :: ns_new
      INTEGER, INTENT(in)                                :: mpol_old
      INTEGER, INTENT(in)                                :: mpol_new
      INTEGER, INTENT(in)                                :: ntor_old
      INTEGER, INTENT(in)                                :: ntor_new
      INTEGER, INTENT(in)                                :: nfp_old
      INTEGER, INTENT(in)                                :: nfp_new
      INTEGER, DIMENSION(-ntor_old:ntor_old), INTENT(in) :: tor_modes_old
      INTEGER, DIMENSION(-ntor_new:ntor_new), INTENT(in) :: tor_modes_new
      LOGICAL, INTENT(in)                                :: lhalf

!  local variables
      INTEGER                                            :: i_n
      INTEGER                                            :: n
      INTEGER                                            :: i_m
      INTEGER                                            :: parity
      INTEGER                                            :: mpol_min
      INTEGER                                            :: mold_off
      INTEGER                                            :: mnew_off
      INTEGER                                            :: nold_off
      INTEGER                                            :: nnew_off

!  Start of executable code
      mold_off = LBOUND(aold,1)
      mnew_off = LBOUND(anew,1)

      nold_off = LBOUND(aold,2) + ntor_old
      nnew_off = LBOUND(anew,2) + ntor_new

      mpol_min = MIN(mpol_new, mpol_old)

!  If the Fourier dimensions of the both arrays are the same just copy them.
!  Otherwise copy only the matcing modes.
      IF (ns_old   .eq. ns_new                  .and.                          &
     &    mpol_old .eq. mpol_new                .and.                          &
     &    ALL(tor_modes_old .eq. tor_modes_new) .and.                          &
     &    nfp_old  .eq. nfp_new) THEN
         anew = aold
      ELSE
         anew = 0
!  Only some of the toroidal modes will match if nfp_new and nfp_old are
!  different.
         DO i_n = -ntor_old, ntor_old
            n = tor_modes_old(i_n)/nfp_old
            IF (.not.fourier_context%get_index(n)) THEN
               CYCLE
            END IF

            DO i_m = 0, mpol_min
               IF (MOD(i_m, 2) .EQ. 0) THEN
                  parity = -1
               ELSE
                  parity = 1
               END IF

               CALL interpit_1d(aold(i_m + mold_off, i_n + nold_off, :),       &
     &                          anew(i_m + mnew_off, n   + nnew_off, :),       &
     &                          ns_old, ns_new, lhalf, parity)
            END DO
         END DO
      END IF

      END SUBROUTINE

!-------------------------------------------------------------------------------
!>  @brief Interpolate radial quantites from the restart file.
!>
!>  Restart files can change the number of surfaces.
!>
!>  @param[in]  aold     Value from the restart file.
!>  @param[out] anew     New interpolated value.
!>  @param[in]  ns_old   Radial grid size in the restart file.
!>  @param[in]  ns_new   New radial grid size.
!-------------------------------------------------------------------------------
      SUBROUTINE interpit_1d(aold, anew, ns_old, ns_new, lhalf, parity)
      USE stel_constants, ONLY: one, zero

      IMPLICIT NONE

!  Declare Arguments
      REAL (dp), DIMENSION(:), INTENT(in)  :: aold
      REAL (dp), DIMENSION(:), INTENT(out) :: anew
      INTEGER, INTENT(in)                  :: ns_old
      INTEGER, INTENT(in)                  :: ns_new
      LOGICAL, INTENT(in)                  :: lhalf
      INTEGER, INTENT(in)                  :: parity

!  local variables
      REAL (dp), ALLOCATABLE, DIMENSION(:) :: temp_y
      REAL (dp), ALLOCATABLE, DIMENSION(:) :: temp_x
      REAL (dp), ALLOCATABLE, DIMENSION(:) :: temp_new
      REAL (dp), ALLOCATABLE, DIMENSION(:) :: temp
      REAL (dp)                            :: ds
      INTEGER                              :: i_s

!  local parameters
      REAL (dp), PARAMETER                 :: lower_b = -1.e30_dp
      REAL (dp), PARAMETER                 :: upper_b = -1.e30_dp

!  Start of executable code
      anew = 0
      IF (ns_old .EQ. ns_new) THEN
         anew = aold
         RETURN
      END IF

!  Need to radially interpolate. Extend the arrays from -ns to ns to get the
!  correct behavior at s=0.
      ds = 1.0/(ns_old - 1)
      IF (lhalf) THEN
         ALLOCATE(temp_y(2*(ns_old - 1)))
         ALLOCATE(temp_x(2*(ns_old - 1)))
         ALLOCATE(temp(2*(ns_old - 1)))

         DO i_s = 2, ns_old
            temp_x(ns_old - i_s + 1) = -ds*(i_s - 1.5)
            temp_y(ns_old + i_s - 2) = ds*(i_s - 1.5)
            temp_y(ns_old - i_s + 1) = parity*aold(i_s)
            temp_y(ns_old + i_s - 2) = aold(i_s)
         END DO
      ELSE
         ALLOCATE(temp_y(2*ns_old - 1))
         ALLOCATE(temp_x(2*ns_old - 1))
         ALLOCATE(temp(2*ns_old - 1))

         DO i_s = 1, ns_old
            temp_x(ns_old - i_s + 1) = -ds*(i_s - 1)
            temp_x(ns_old + i_s - 1) = ds*(i_s - 1)
            temp_y(ns_old - i_s + 1) = parity*aold(i_s)
            temp_y(ns_old + i_s - 1) = aold(i_s)
         END DO
      END IF

      CALL spline(temp_x, temp_y, SIZE(temp_x), lower_b, upper_b, temp)

      ALLOCATE(temp_new(ns_new))
      ds = 1.0/(ns_new - 1)
      IF (lhalf) THEN
         DO i_s = 2, ns_new
            temp_new(i_s) = ds*(i_s - 1.5)
         END DO
      ELSE
         DO i_s = 1, ns_new
            temp_new(i_s) = ds*(i_s - 1)
         END DO
      END IF

      CALL splINT(temp_x, temp_y, temp, SIZE(temp_x), temp_new, anew)

      DEALLOCATE(temp_y)
      DEALLOCATE(temp_x)
      DEALLOCATE(temp)

      END SUBROUTINE

!-------------------------------------------------------------------------------
!>  @brief Denormalize a quantity so the value written to the restart file can
!>         be summed directly.
!>
!>  This removes the orthonorm and energy scale factors for a fourier quantity.
!>
!>  @param[inout] xmn    Fourier quantity to denormalize.
!>  @param[in]    factor Energy scale factor.
!-------------------------------------------------------------------------------
      SUBROUTINE restart_normalize(xmn, factor)
      USE island_params, ONLY: fourier_context

      IMPLICIT NONE

!  Declare Arguments
      REAL (dp), DIMENSION(:,:,:) :: xmn
      REAL (dp)                   :: factor

!  Local Variables
      INTEGER                     :: s

!  Start of executable code.
      DO s = 1, SIZE(xmn,3)
         xmn(:,:,s) = xmn(:,:,s)*factor/fourier_context%orthonorm
      END DO

      END SUBROUTINE

!-------------------------------------------------------------------------------
!>  @brief Denormalize a quantity so the value written to the restart file can
!>         be summed directly.
!>
!>  This removes the orthonorm and energy scale factors for a fourier quantity.
!>
!>  @param[inout] xmn    Fourier quantity to denormalize.
!>  @param[in]    factor Energy scale factor.
!-------------------------------------------------------------------------------
      SUBROUTINE restart_denormalize(xmn, factor)
      USE island_params, ONLY: fourier_context

      IMPLICIT NONE

!  Declare Arguments
      REAL (dp), DIMENSION(:,:,:) :: xmn
      REAL (dp)                   :: factor

!  Local Variables
      INTEGER                     :: s

!  Start of executable code.
      DO s = 1, SIZE(xmn,3)
         xmn(:,:,s) = fourier_context%orthonorm*xmn(:,:,s)/factor
      END DO

      END SUBROUTINE

      END MODULE
