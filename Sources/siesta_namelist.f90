!-------------------------------------------------------------------------------
!  The @header, @table_section, @table_subsection, @item and @end_table commands
!  are custom defined commands in Doxygen.in. They are defined under ALIASES.
!  For the page created here, the 80 column limit is exceeded. Arguments of
!  aliases are separated by ','. If you intended ',' to be a string you must use
!  an escaped comma '\,'.
!
!>  @page siesta_namelist_sec Namelist siesta_main_nli definition
!>  
!>  @tableofcontents
!>  @section siesta_namelist_intro_sec Introduction
!>  This page documents the contents of a SIESTA namelist input file. SIESTA 
!>  namelist variables are defined in the @fixed_width{siesta_info} common 
!>  block.
!>
!>  @section siesta_namelist_var_sec Namelist Variables
!>  @header{Input variable, Description, Code Reference}
!>
!>  @table_section{siesta_flag_sec, Control Flags}
!>     @item{ladd_pert,        Use helical perturbation.,                        siesta_namelist::ladd_pert}
!>     @item{lresistive,       Use resistive perturbation.,                      siesta_namelist::lresistive}
!>     @item{lrestart,         Restart SIESTA from pervious run.,                siesta_namelist::lrestart}
!>     @item{l_tracing,        Produce output file for fieldline tracing.,       siesta_namelist::l_tracing}
!>     @item{lcolscale,        Use column scaling.,                              shared_data::lcolscale}
!>     @item{l_silo_output,    Produce silo output.,                             siesta_namelist::l_silo_output}
!>     @item{l_silo3D,         Produce silo 3D output.,                          siesta_namelist::l_silo3d}
!>     @item{l_output_alliter, Write output files on all iterations instead of
!>                             only iterations that lower the MHD energy and
!>                             force residuals.,                                 siesta_namelist::l_output_alliter}
!>     @item{l_VMEC_Uniform,   UNKNOWN,                                          siesta_namelist::l_vmec_uniform}
!>     @item{lasym,            Use non stellarator symmetric terms.,             shared_data::lasym}
!>     @item{l_vessel,         Use free-boundary grid extension to vessel wall., siesta_namelist::l_vessel}
!>     @item{lrecon,           Add additional output to the restart file when
!>                             used in a reconstruction context.
!>                             DEPRICATED,                                       shared_data::lrecon}
!>  @end_table
!>
!>  @table_section{siesta_algrothim_sec, Algrothim Control Variables}
!>     @item{niter,         Maximum number of iterations after diagonal prec.,  siesta_namelist::niter}
!>     @item{ftol,          Minimum force tolarance for converged solution.,    siesta_namelist::ftol}
!>     @item{mupar,         Resistivity factor.,                                hessian::mupar}
!>     @item{levmarq_param, Inital value of Levenberg-Marquardt parameter.,     hessian::levmarq_param}
!>     @item{eta_factor,    Resistivity value.,                                 siesta_namelist::eta_factor}
!>     @item{nprecon,       Skip diagonal preconditioner if greater than zero., siesta_namelist::nprecon}
!>     @item{ngmres_type,   Preconditioner type.,                               shared_data::ngmres_type}
!>     @item{iortho,        UNKNOWN,                                            shared_data::iortho}
!>  @end_table
!>
!>  @table_section{siesta_island_sec, Island Parameters}
!>     @item{mres,    M numbers of island resonances.,    siesta_namelist::mres}
!>     @item{helpert, Sizes of the helical perturbation., siesta_namelist::helpert}
!>  @end_table
!>
!>  @table_section{siesta_grid_size_sec, Grid Sizes}
!>     @item{nsin,       Size of plasma radial grid.,                            siesta_namelist::nsin}
!>     @item{nsin_ext,   Size of extended radial grid.,                          siesta_namelist::nsin_ext}
!>     @item{ntor_type,  Type of toroidal modes. Available types are
!>                       -# @fixed_width{'dense'}  Use all toroidal modes.
!>                       -# @fixed_width{'sparse'} Use selected toroidal modes., siesta_namelist::ntor_type}
!>     @item{mpolin,     Number of poloidal modes.,                              siesta_namelist::mpolin}
!>     @item{ntorin,     Number of toroidal modes.,                              siesta_namelist::ntorin}
!>     @item{nfpin,      Number of field periods to use. Setting this to
!>                       anything less than one will use the value form the wout
!>                       file.,                                                  siesta_namelist::nfpin}
!>     @item{ntor_modes, Sparse toroidal mode numbers.,                          siesta_namelist::ntor_modes}
!>     @table_subsection{siesta_grid_size_out_sec, Output Grid Sizes}
!>        @item{nphis, Number of cylindrical phi planes.,     siesta_namelist::nphis}
!>        @item{nrs,   Number of radial grid points.,         siesta_namelist::nrs}
!>        @item{nzs,   Number of vertical grid points.,       siesta_namelist::nzs}
!>        @item{nvs,   Number of flux space toroidal points., siesta_namelist::nvs}
!>        @item{nus,   Number of flux space poloidal points., siesta_namelist::nus}
!>        @item{nss,   Number of flux space radial points.,   siesta_namelist::nss}
!>  @end_table
!>
!>  @table_section{siesta_file_name_sec, File Names}
!>      @item{wout_file,   Filename of the VMEC woutfile.,      siesta_namelist::wout_file}
!>      @item{restart_ext, Name of the restart file extension., siesta_namelist::restart_ext}
!>      @item{mgrid_file,  Filename of the MGRID file.,         siesta_namelist::mgrid_file}
!>      @item{vessel_file, Filename of the extended surfaces.,  siesta_namelist::vessel_file}
!>  @end_table
!>
!>  @table_section{siesta_test_sec, Test Controls}
!>      @item{hesspass_test, UNKNOWN, shared_data::hesspass_test}
!>      @item{mupar_test,    UNKNOWN, shared_data::mupar_test}
!>  @end_table
!>
!>  @section siesta_namelist_prog_ref_sec Programmers Reference
!>  Reference material for the coding to implement this namelist is found in the
!>  @ref siesta_namelist module.
!-------------------------------------------------------------------------------
!*******************************************************************************
!>  @file siesta_namelist.f90
!>  @brief Contains module @ref siesta_namelist.
!
!  Note separating the Doxygen comment block here so detailed decription is
!  found in the Module not the file.
!
!>  This file contains all the variables and maximum sizes of the inputs for a
!>  SIESTA namelist input file. The module contained within does not represent 
!>  an object instance. Instead all variables are contained in a global context.
!>  This is required due to limitations of FORTRAN 95 and namelist inputs.
!>
!>  @ref siesta_namelist_sec "Namelist siesta_info definition"
!>
!>  @note Some of the references are missing here. This is due to a bug in
!>  Doxygen when variable decalarations span multiple lines.
!*******************************************************************************
      MODULE siesta_namelist
      USE shared_data, ONLY: ngmres_type, iortho, lcolscale, lasym,            &
                             hesspass_test, mupar_test, lrecon
      USE Hessian, ONLY: levmarq_param, mupar
      USE stel_kinds

      IMPLICIT NONE

!*******************************************************************************
!  siesta_namelist input module parameters
!*******************************************************************************
!>  Input string length.
      INTEGER, PARAMETER ::  siesta_namelist_name_length = 256

!*******************************************************************************
!  DERIVED-TYPE DECLARATIONS
!  1) siesta_namelist_class
!
!*******************************************************************************
!  Control Flags
!>  Use helical perturbation.
      LOGICAL :: ladd_pert = .TRUE.
!>  Use resistive perturbaton.
      LOGICAL :: lresistive = .TRUE.
!>  Restart SIESTA from pervious run.
      LOGICAL :: lrestart = .FALSE.
!>  Produce output file for fieldline tracing.
      LOGICAL :: l_tracing = .FALSE.
!>  Produce silo output.
      LOGICAL :: l_silo_output = .FALSE.
!>  Produce silo 3D output.
      LOGICAL :: l_silo3D = .FALSE.
!>  Write output files on all iterations.
      LOGICAL :: l_output_alliter = .FALSE.
!>  FIXME: Unknown
      LOGICAL :: l_VMEC_Uniform
!>  If extended grid is to be used using an available vessel file
      LOGICAL :: l_vessel = .FALSE.
!>  Recompute lambda on the SIESTA grid.
      LOGICAL :: l_lambda = .TRUE.

!  Algrothim Control Variables}
!>  Maximum number of iterations after diagonal prec.
      INTEGER  :: niter = 10
!>  Force tolarance.
      REAL(dp) :: ftol = 1.E-20_dp
!>  Resistivity value.
      REAL(dp) :: eta_factor = 1.E-2_dp
!>  Skip diagonal preconditioner if greater than zero.
      INTEGER  :: nprecon = 0

!  Island parameters
!>  Sizes of the helical perturbation.
      INTEGER, DIMENSION(20)  :: mres = 0
!>  Sizes of the helical perturbation.
      REAL(dp), DIMENSION(20) :: HelPert = 0.0
!>  Sizes of the helical perturbation.
      REAL(dp), DIMENSION(20) :: HelPertA = 0.0

!  Grid Sizes
!>  Radial size of the plasma grid.
      INTEGER :: nsin = 101
!>  Radial size of the extended grid.
      INTEGER :: nsin_ext = 0
!>  Type of toroidal modes.
      CHARACTER(LEN=siesta_namelist_name_length) :: ntor_type = 'dense'
!>  Number of poloidal modes.
      INTEGER :: mpolin = 12
!>  Number of toroidal modes.
      INTEGER :: ntorin = 3
!>  Number of field periods to use. -1 means set this to the value in the wout
!>  file
      INTEGER :: nfpin = 0
!>  Sparse toroidal modes.
      INTEGER, DIMENSION(-50:50) :: ntor_modes

!  Output Grid Sizes
!>  Number of cylindrical phi planes.
      INTEGER :: nphis = 2
!>  Number of radial grid points.
      INTEGER :: nrs = 200
!>  Number of vertical grid points.
      INTEGER :: nzs = 200
!>  Number of flux space toroidal points.
      INTEGER :: nvs = 150
!>  Number of flux space poloidal points.
      INTEGER :: nus = 150
!>  Number of flux space radial points.
      INTEGER :: nss = 100

!  File Names
!>  Filename of the VMEC woutfile.
      CHARACTER(LEN=siesta_namelist_name_length) :: wout_file = ''
!>  Name of the restart file extension.
      CHARACTER(LEN=siesta_namelist_name_length) :: restart_ext = ''
!>  Filename of the VMEC woutfile.
      CHARACTER(LEN=siesta_namelist_name_length) :: mgrid_file = ''
!>  Name of the restart file extension.
      CHARACTER(LEN=siesta_namelist_name_length) :: vessel_file = ''

!  Declare Namelist
      NAMELIST/siesta_info/                                                    &
!  Control flags
        ladd_pert, lresistive, lrestart, l_tracing, lcolscale,                 &
        l_silo_output, l_silo_output, l_silo3D, l_output_alliter,              &
        l_VMEC_Uniform, lasym, lrecon, l_vessel, l_lambda,                     &
!  Algrothim Control Variables
        niter, ftol, mupar, levmarq_param, eta_factor, nprecon,                &
        ngmres_type, iortho,                                                   &
!  Island parameters Island Parameters
        mres, HelPert, HelPertA,                                               &
!  Input grid sizes
        nsin, nsin_ext, ntor_type, mpolin, ntorin, nfpin, ntor_modes,          &
!  Output grid sizes
        nphis, nrs, nzs, nvs, nus, nss,                                        &
!  File names
        wout_file, restart_ext, mgrid_file, vessel_file,                       &
!  Test controls
        hesspass_test, mupar_test

      CONTAINS

!*******************************************************************************
!  UTILITY SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  @brief Reads the namelist input file.
!>
!>  Reads the namelist input file.
!>
!>  @param[in] namelist_file The file name of the namelist input file.
!-------------------------------------------------------------------------------
      SUBROUTINE siesta_namelist_read(namelist_file)
      USE safe_open_mod
      USE v3_utilities
      USE Hessian, ONLY: levmarq_param, mupar, levmarq_param0, mupar0

      IMPLICIT NONE

!  Declare Arguments
      CHARACTER (len=*), INTENT(in) :: namelist_file

!  local variables
      INTEGER                       :: iou_mnli
      INTEGER                       :: status
      INTEGER                       :: i

!  Start of executable code
      levmarq_param = 1.E-3_dp
      mupar         = 0
      niter         = 10
      mres          = 0
      HelPert       = 0
      HelPertA      = 0
      lcolscale     = .TRUE.
      mupar_test    = 0
      ntor_type     = 'dense'
      ntor_modes    = 0

!  Initalize a default value of the I\O unit. SIESTA increments from there.
      iou_mnli = 0
      CALL safe_open(iou_mnli, status, TRIM(namelist_file),                    &
                    'old', 'formatted')
      CALL assert_eq(0, status, 'siesta_namelist_read' //                      &
        ': Safe_open of ' // TRIM(namelist_file) // ' failed')

!  Read the namelist input file.
      READ (iou_mnli, nml=siesta_info)
      CLOSE (iou_mnli, iostat=status)
      CALL assert_eq(0, status, 'siesta_namelist_read' //                      &
        ': Error closing ' // TRIM(namelist_file) // ' failed')

      levmarq_param0 = levmarq_param
      mupar0         = mupar

      IF (TRIM(ntor_type) .eq. 'dense') THEN
         DO i = -ntorin, ntorin
            ntor_modes(i) = i
         END DO
      END IF

      END SUBROUTINE

!-------------------------------------------------------------------------------
!>  @brief Writes the namelist input file.
!>
!>  Writes the namelist input file.
!>
!>  @param[in] namelist_file The file name of the namelist input file.
!-------------------------------------------------------------------------------
      SUBROUTINE siesta_namelist_write(namelist_file)
      USE safe_open_mod
      USE v3_utilities

      IMPLICIT NONE

!  Declare Arguments
      CHARACTER (len=*), INTENT(in) :: namelist_file

!  local variables
      INTEGER                       :: iou_mnli
      INTEGER                       :: status

!  Start of executable code

!  Initalize a default value of the I\O unit. SIESTA increments from there.
      iou_mnli = 0
      CALL safe_open(iou_mnli, status, TRIM(namelist_file),                    &
     &               'replace', 'formatted', delim_in='quote')
      CALL assert_eq(0, status, 'siesta_namelist_write' //                     &
     &   ': Safe_open of ' // TRIM(namelist_file) // ' failed')

!  Write the namelist input file.
      WRITE (iou_mnli, nml=siesta_info)
      CLOSE (iou_mnli, iostat=status)
      CALL assert_eq(0, status, 'siesta_namelist_read' //                      &
     &   ': Error closing ' // TRIM(namelist_file) // ' failed')

      END SUBROUTINE

      END MODULE
