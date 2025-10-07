!*******************************************************************************
!>  @file island_params.f90
!>  @brief Contains module @ref island_params.
!
!  Note separating the Doxygen comment block here so detailed decription is
!  found in the Module not the file.
!
!>  This file contains fix parameters related to the computational grids.
!*******************************************************************************
      MODULE island_params
      USE stel_kinds
      USE fourier, ONLY: fourier_class

      IMPLICIT NONE

!*******************************************************************************
!  island_params module variables
!*******************************************************************************
!>  Number of radial grid points.
	  INTEGER                             :: ns_i
!>  Number of polodial grid points. FIXME: This should go in the fourier_class
      INTEGER                             :: nu_i
!>  Number of toroidal grid points. FIXME: This should go in the fourier_class
      INTEGER                             :: nv_i
!>  Number of polodial and toroidal grid points. FIXME: Is this really needed?
      INTEGER                             :: nuv_i
!>  Number of poloidal modes. FIXME: This should go in the fourier_class
      INTEGER                             :: mpol_i
!>  Number of toroidal modes. FIXME: This should go in the fourier_class
      INTEGER                             :: ntor_i
!>  Number of radial grid vmec grid points. FIXME: This should be in vmec_info
      INTEGER                             :: ns_v
!>  Number of poloidal vmec modes. FIXME: This should be in vmec_info
      INTEGER                             :: mpol_v
!>  Number of toroidal vmec modes. FIXME: This should be in vmec_info
      INTEGER                             :: ntor_v
!>  Number of field periods.  FIXME: This should go in the fourier_class
      INTEGER                             :: nfp_i
!>  Number of total toroidal and poloidal vmec modes. FIXME: This should be in vmec_info
      INTEGER                             :: mnmax_i
!>  Number of radial points on the half mesh. FIXME: is this really needed?
      INTEGER                             :: nsh
!  FIXME: These don't get used. Delete them.
!      INTEGER, ALLOCATABLE                :: irefu(:)
!      INTEGER, ALLOCATABLE                :: jrefv(:)                          ! -u, -v indices
!>  Radial grid spacing. hs = s_i+1 - s_i
      REAL (dp)                           :: hs_i
!>  Radial grid derivative factor.
!>    d X/ ds = (x_i+1 - x_i)/(s_i+1 - s_i)                                  (1)
!>  Where ohs = 1/(s_i+1 - s_i) = ns - 1
      REAL (dp)                           :: ohs_i
!>  Differential area element.
!>    dudv 1/(ntheta*nzeta)                                                  (1)
      REAL (dp)                           :: dnorm_i
!>  
      REAL (dp)                           :: gnorm_i
!>  Magnetic energy. 2Pi*wb @ref read_wout_mod::wb.
      REAL (dp)                           :: wb_i
!>  Thermal energy. 2Pi*wp @ref read_wout_mod::wp.
      REAL (dp)                           :: wp_i
!>  Equilibrium volume @ref read_wout_mod::volume.
      REAL (dp)                           :: volume_i
!>  Major radius.
      REAL (dp)                           :: rmajor_i
!>  Adiabatic constant.
      REAL(dp), PARAMETER                 :: gamma = 5._dp/3._dp
!>  Radial toroidal flux derivative.
      REAL(dp), ALLOCATABLE, DIMENSION(:) :: phipf_i
!>  Radial poloidal flux derivative.
      REAL(dp), ALLOCATABLE, DIMENSION(:) :: chipf_i
!>  Radial toroidal flux.
      REAL(dp), ALLOCATABLE, DIMENSION(:) :: phif_i
!>  Radial poloidal flux.
      REAL(dp), ALLOCATABLE, DIMENSION(:) :: chif_i
!>  Radial pressure. FIXME: Check if this is really needed.
      REAL(dp), ALLOCATABLE, DIMENSION(:) :: presf_i
!>  Volume of a radial grid surface.
      REAL(dp), ALLOCATABLE, DIMENSION(:) :: vp_f
!>  Fourier transform object.
      TYPE (fourier_class), POINTER       :: fourier_context => null()

      END MODULE
