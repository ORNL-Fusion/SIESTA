!*******************************************************************************
!>  @file shared_data.f90
!>  @brief Contains module @ref shared_data.
!
!  Note separating the Doxygen comment block here so detailed decription is
!  found in the Module not the file.
!
!>  This file contains variables and parameters used by many modules in SIESTA.
!*******************************************************************************
      MODULE shared_data
      USE stel_kinds
      USE stel_constants

      IMPLICIT NONE

!*******************************************************************************
!  shared_data module parameters
!*******************************************************************************
!  Solver parameters
!>  Max number of gmres steps (10-100) should scale with ns
      INTEGER, PARAMETER   :: ngmres_steps = 100
!>  Diagonal preconditioning flag.
      INTEGER, PARAMETER   :: PREDIAG = 1
!>  Block preconditioning flag.
      INTEGER, PARAMETER   :: PREBLOCK = 2
!>  Threshold force for turning off resistive perturbations.
      REAL (dp), PARAMETER :: fsq_res = 1.E-16_dp
!>  FIXME: UNKNOWN
      REAL (dp), PARAMETER :: levm_ped = 1.E-10_dp
!>  Pedestal value of levenberg/mu. Should be between 10^-5 and 10^-10.
      REAL (dp), PARAMETER :: mu_ped = 1.E-8_dp
!>  GMRES peak improvement
      INTEGER, PARAMETER   :: gmres_peak = 2
!>  GMRES no_peak improvement.
      INTEGER, PARAMETER   :: gmres_no_peak = 1

!  IO Parameters
!>  File output io unit.
      INTEGER, PARAMETER   :: unit_out = 336

!  Force calulation flags.
!>  Preserve s=1 as iso-pressure surface.
      LOGICAL, PARAMETER   :: l_pedge = .true.
!>  Natural boundry condition flag.
!>    * True Use natural bc at s=0 (preferred: maintains symmetry)
!>    * False Evolves all forces at origin.
      LOGICAL, PARAMETER   :: l_natural = .TRUE.

!*******************************************************************************
!  shared_data module variables
!*******************************************************************************
!  Size variables.
!>  Number of elements in the xc array.
      INTEGER                                     :: neqs
!>  Number of independent variables.
!>    * 3 Stellarator symmetry.
!>    * 6 Stellarator symmetry.
      INTEGER                                     :: ndims
!>  Total number of iteration to run.
      INTEGER                                     :: niter
!>  Block size. (mpol + 1)*(2*ntor + 1)*ndims
      INTEGER                                     :: mblk_size
!>  Total radial grid size in the VMEC region.
      INTEGER                                     :: nsp

!  Solver control variables.
!>  Preconditioner flag.
      INTEGER                                     :: nprecon
!>  Preconditioner type.
      INTEGER                                     :: nprecon_type
!>  GMRES control flag.
!>    * @ref gmres_peak
!>    * @ref gmres_no_peak
      INTEGER                                     :: ngmres_type = gmres_peak
!>  Orthogonalization in GMRES.
!>    * 3 ICGS (Recommended)
!>    * 2 CGS
!>    * 1 IMGS
!>    * 0 MGS
!>  FIXME: Make this a parameter.
      INTEGER                                     :: iortho = 3
!>  Dump block and data files for testing.
      INTEGER                                     :: hessPass_Test = -1
!>  FIXME UNKNOWN
      INTEGER                                     :: in_hess_nfunct
!>  FIXME UNKNOWN
      INTEGER                                     :: out_hess_nfunct
!>  FIXME UNKNOWN
      REAL (dp)                                   :: mupar_test
!>  |F|^2 WITH column scaling.
      REAL (dp)                                   :: fsq_total
!>  |F|^2 WITHOUT column scaling.
      REAL (dp)                                   :: fsq_total1

!  Solver work variables.
!>  1D array of Fourier mode displacement components.
      REAL(dp), DIMENSION(:), ALLOCATABLE         :: xc
!>  1D Array of Fourier mode MHD force components
      REAL(dp), DIMENSION(:), ALLOCATABLE, TARGET :: gc
!>  1D Array of Fourier mode MHD force components, FIXME Check if this is really
!>  needed.
      REAL(dp), DIMENSION(:), ALLOCATABLE, TARGET :: gc_sup
!>  Resonant magnetic field perturbation.
      REAL(dp), DIMENSION(:,:,:), ALLOCATABLE     :: buv_res
!>  Saved fouier displacements.
      REAL(dp), DIMENSION(:), ALLOCATABLE         :: xc0
!>  Saved fouier MHD forces.
      REAL(dp), DIMENSION(:), ALLOCATABLE         :: gc0
!>  Column scaling factors.
      REAL(dp), DIMENSION(:,:,:,:), ALLOCATABLE   :: col_scale  !<column scaling factor

!>  |F|^2 for GMRES iterations.
      REAL (dp)                                   :: fsq_gmres
!>  Linear |F|^2.
      REAL (dp)                                   :: fsq_lin
!>  FIXME: UNKNOWN
      REAL (dp)                                   :: etak_tol
!>  FIXME: UNKNOWN
      REAL (dp)                                   :: levm_scale = 1
!>  MHD energy sum of magnetic and thermal.
      REAL (dp)                                   :: wtotal
!>  Saved MHD energy sum of magnetic and thermal.
      REAL (dp)                                   :: wtotal0
!>  Time step.
      REAL (dp)                                   :: delta_t
!>  |F|^2 for s components.
      REAL (dp)                                   :: fsqvs
!>  |F|^2 for u components.
      REAL (dp)                                   :: fsqvu
!>  |F|^2 for v components.
      REAL (dp)                                   :: fsqvv
!>  Spectral Truncation RMS error,
      REAL (dp)                                   :: ste(4)
!>  FIXME: UNKNOWN
      REAL (dp)                                   :: bs0(12)
!>  FIXME: UNKNOWN
      REAL (dp)                                   :: bu0(12)
!>  FIXME: UNKNOWN
      REAL (dp)                                   :: bsbu_ratio_s
!>  FIXME: UNKNOWN
      REAL (dp)                                   :: jsju_ratio_s
!>  FIXME: UNKNOWN
      REAL (dp)                                   :: bsbu_ratio_a
!>  FIXME: UNKNOWN
      REAL (dp)                                   :: jsju_ratio_a
!>  FIXME: UNKNOWN
      REAL (dp)                                   :: scale_s
!>  FIXME: UNKNOWN
      REAL (dp)                                   :: scale_u

!  Cordinate basis.
!>  R coordinates of the computational grid.
      REAL(dp), DIMENSION(:,:,:), ALLOCATABLE     :: r1_i
!>  Z coordinates of the computational grid.
      REAL(dp), DIMENSION(:,:,:), ALLOCATABLE     :: z1_i
!>  dR/du coordinates of the computational grid.
      REAL(dp), DIMENSION(:,:,:), ALLOCATABLE     :: ru_i
!>  dZ/du coordinates of the computational grid.
      REAL(dp), DIMENSION(:,:,:), ALLOCATABLE     :: zu_i
!>  dR/dv coordinates of the computational grid.
      REAL(dp), DIMENSION(:,:,:), ALLOCATABLE     :: rv_i
!>  dZ/dv coordinates of the computational grid.
      REAL(dp), DIMENSION(:,:,:), ALLOCATABLE     :: zv_i

!  Shared quantities.
!>  FIXME: UNKNOWN
      REAL (dp)                                   :: jsupvdotA
!>  Toroidal flux profile.
      REAL (dp), DIMENSION(:), ALLOCATABLE        :: torflux
!>  Poloidal flux profile.
      REAL (dp), DIMENSION(:), ALLOCATABLE        :: polflux
!>  Covariant vector potential for stellator symmetric s component on full grid.
      REAL (dp), DIMENSION(:,:,:), POINTER        :: asubsmnsf
!>  Covariant vector potential for non-stellator symmetric s component on full
!>  grid.
      REAL (dp), DIMENSION(:,:,:), POINTER        :: asubsmncf
!>  Covariant vector potential for stellator symmetric u component on full grid.
      REAL (dp), DIMENSION(:,:,:), POINTER        :: asubumncf
!>  Covariant vector potential for non-stellator symmetric u component on full
!>  grid.
      REAL (dp), DIMENSION(:,:,:), POINTER        :: asubumnsf
!>  Covariant vector potential for stellator symmetric v component on full grid.
      REAL (dp), DIMENSION(:,:,:), POINTER        :: asubvmncf
!>  Covariant vector potential for non-stellator symmetric v component on full
!>  grid.
      REAL (dp), DIMENSION(:,:,:), POINTER        :: asubvmnsf
      
!  Flags to control evolving origin and edge.
!>  Solve u,v components at s=1.
      LOGICAL :: l_push_edge
!>  Solve for s component at origin.
      LOGICAL :: l_push_s
!>  Solve for u component at origin.
      LOGICAL :: l_push_u
!>  Solve for v component at origin.
      LOGICAL :: l_push_v
!>  Use linearized forces.
      LOGICAL :: l_linearize
!>  FIXME: UNKNOWN
      LOGICAL :: l_conjgrad
!>  Compute MHD energy.
      LOGICAL :: l_getwmhd
!>  Compute |F|^2.
      LOGICAL :: l_getfsq
!>  Apply preconditioner.
      LOGICAL :: l_ApplyPrecon
!>  Print forces at the origin.
      LOGICAL :: l_PrintOriginForces = .false.
!>  Store initial field/pressure state.
      LOGICAL :: l_init_state
!>  Update the @ref ste array.
      LOGICAL :: l_update_state = .false.
!>  Parallel allocated quantities? FIXME: check this.
      LOGICAL :: l_par_state
!>  Apply column scaling to hessian.
      LOGICAL :: lcolscale
!>  Use non-stellarator symmetry.
      LOGICAL :: lasym = .false.
!>  Output extra information to the restart file that will be used by V3FIT.
!>  @DEPRICATED
      LOGICAL :: lrecon = .false.
!>  Use verbose screen output.
      LOGICAL :: lverbose = .true.
!>  Equilibrate matrix with col 1-norm
      LOGICAL :: lequi1 = .true.

!>  Total toroidal current
      REAL (dp) :: siesta_curtor = 0.0

      END MODULE shared_data
