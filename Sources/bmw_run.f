!*******************************************************************************
!>  @file bmw.f
!>  @brief Contains the main routines for Biot-Savart Magnetic Vmec
!>  Vector-Potential (BMW) code.
!
!  Note separating the Doxygen comment block here so detailed decription is
!  found in the Module not the file.
!
!>  BMW is a code for extending fields belond the VMEC domain in a manner that
!>  ensures divergence free fields. BMW does a Biot-Savart volume integration of
!>  of the equilibrium current density to obtain a continous vector potential
!>  every where on the mgrid grid.
!>
!>  @author Mark Cianciosa
!>
!>  Below is a brief discription of the major top level objects of the code. For
!>  discriptions of lower level objects consult the referenced top level
!>  objects.
!>
!>  @ref m_grid        Object containing the vaccum field information.
!>  @ref primed_grid   Object containing the plasma currents and primed grid
!>                     positions.
!>  @ref unprimed_grid Object containing the plasma vector potential response.
!>  @ref bmw_context   Defines the main pointers to all allocated memory and
!>                     objects.
!*******************************************************************************
      MODULE bmw_run
      USE bmw_context
      USE bmw_commandline_parser
      USE, INTRINSIC :: iso_fortran_env, Only : output_unit
      USE bmw_state_flags

      IMPLICIT NONE

      CONTAINS 

!-------------------------------------------------------------------------------
!>  @brief BMW main program.
!>
!>  Highest level BMW routine. This computes the vector potential on the mgrid
!>  grid.
!>
!>  @param[in]  mgrid_file_name Path to the mgrid file.
!>  @param[in]  wout_file_mame  Path to the wout file.
!>  @param[in]  r_grid          Radial grid points to compute the vector
!>                              potential on.
!>  @param[in]  z_grid          Vertial grid points to compute the vector
!>                              potential on.
!>  @param[in]  dphi            Grid size in the phi direction.
!>  @param[out] A_r             Vector potential in the cyclindical R direction.
!>  @param[out] A_p             Vector potential in the cyclindical Phi
!>                              direction.
!>  @param[out] A_z             Vector potential in the cyclindical Z direction.
!-------------------------------------------------------------------------------
      SUBROUTINE bmw_exec(mgrid_file_name, wout_file_name,                     &
     &                    r_grid, z_grid, dphi, A_r, A_p, A_z)
      USE read_wout_mod, ONLY: nfp_vmec=>nfp
      USE island_params, ONLY: nfp_i

!  Declare Arguments
      CHARACTER (len=*), INTENT(in)                 :: mgrid_file_name
      CHARACTER (len=*), INTENT(in)                 :: wout_file_name
      REAL (rprec), DIMENSION(:,:,:), INTENT(in)    :: r_grid
      REAL (rprec), DIMENSION(:,:,:), INTENT(in)    :: z_grid
      REAL (rprec), INTENT(in)                      :: dphi
      REAL (rprec), DIMENSION(:,:,:), INTENT(out)   :: A_r
      REAL (rprec), DIMENSION(:,:,:), INTENT(out)   :: A_p
      REAL (rprec), DIMENSION(:,:,:), INTENT(out)   :: A_z

!  Local Variables
      CLASS (bmw_parallel_context_class), POINTER :: parallel => null()
      CLASS (bmw_context_class), POINTER            :: context => null()
      CLASS (bmw_commandline_parser_class), POINTER ::                         &
     &   cl_parser => null()
      INTEGER                                       :: flags
      INTEGER                                       :: num_p
      INTEGER                                       :: status
      INTEGER, DIMENSION(3)                         :: dims

!  Start of executable code
      CALL profiler_construct

      parallel => bmw_parallel_context_class(MPI_COMM_WORLD)
      flags = bmw_state_flags_off

      CALL parallel%set_threads(-1)
      CALL parallel%report(output_unit)

!  BMW uses num_p to set the number of grid points for the primed grid. In doing
!  so, it mutliplies by the number of field periods in the mgrid file. If using
!  a user defined number of field periods that is less than the amount in mgrid,
!  we need to scale num_p.
      num_p = SIZE(A_r, 2)*2.0/(nfp_vmec/nfp_i)
      context => bmw_context_class(mgrid_file_name, wout_file_name,            &
     &                             '', flags, num_p, parallel,                 &
     &                             output_unit)

      dims(1) = SIZE(r_grid, 1)
      dims(2) = SIZE(r_grid, 2)
      dims(3) = SIZE(r_grid, 3)

      CALL context%set_up_grid(RESHAPE(r_grid, dims),                          &
     &                         RESHAPE(z_grid, dims),                          &
     &                         dphi, parallel, output_unit)

      A_r = context%up_grid%a_r
      A_p = context%up_grid%a_p
      A_z = context%up_grid%a_z

      DEALLOCATE(context)
      CALL profiler_destruct

      DEALLOCATE(parallel)

1000  FORMAT('BMW ',i4,' Series.')

      END SUBROUTINE

      END MODULE
