
!*******************************************************************************
!>  @file siesta.f90
!>  @brief Contains main routines for SIESTA.
!
!  Note separating the Doxygen comment block here so detailed decription is
!  found in the Module not the file.
!
!>  SIESTA is an ideal MHD equilibrium code that allows for islands and
!>  stocastic fields. SIESTA has been created with contributions from.
!>
!>  @author Steve Hirshman
!>  @author R. Sanchez
!>  @author Sudip Seal
!>  @author Mark Cianciosa
!>
!>  The theory behind SIESTA is outlined in
!>     S.P. Hirshman et. al. https://doi.org/10.1063/1.3597155
!>  The design of the parallelization method is outlined in
!>     S.K. Seal at. al. https://doi.org/10.1002/cpe.2919
!>
!>  Please report any bugs to ciancisamr@ornl.gov
!*******************************************************************************
!*******************************************************************************
!  MAIN PROGRAM
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  @brief SIESTA Main program.
!>
!>  Highest level SIESTA routine.
!-------------------------------------------------------------------------------
      PROGRAM SIESTA
      USE siesta_run
      USE mpi_inc
      USE siesta_namelist, ONLY: lrestart

      IMPLICIT NONE

!  Local variables
      CLASS (siesta_run_class), POINTER :: context => NULL()
      CHARACTER(LEN=256)                :: temp
      CHARACTER(LEN=100)                :: command_arg
      INTEGER                           :: numargs
      INTEGER                           :: length

!  Local parameters
#if defined(MPI_OPT)
      INTEGER, PARAMETER                :: run_comm = MPI_COMM_WORLD
#else
      INTEGER, PARAMETER                :: run_comm = 0
#endif

!  Read Command line arguments if present.
      CALL getcarg(1, command_arg, numargs)

      IF (numargs .gt. 0) THEN
         IF (LEN_TRIM(temp) .ne. 0) THEN
            IF (command_arg(1:7) .ne. 'siesta_' .AND.                          &
                command_arg(1:7) .ne. 'SIESTA_') THEN
               temp = 'siesta_' // TRIM(command_arg)
            END IF
            length = LEN_TRIM(temp)
            IF (temp(length - 3:length) .ne. '.jcf' .AND.                        &
                temp(length - 3:length) .ne. '.JCF') THEN
               temp = TRIM(temp) // '.jcf'
            END IF
         END IF
      ELSE
         temp = 'siesta.jcf'
      END IF

!  Start of executable code.
!  Create a context and read the SIESTA namelist input file.
      context => siesta_run_class(run_comm, .true., .true., .true., temp)

!  Check the restart flag in from the namelist input. If true, restart from an
!  existing restart file. Otherwise initalize from a VMEC file.
      IF (lrestart) THEN
         CALL context%set_restart
      ELSE
         CALL context%set_vmec(.true.)
      END IF

      CALL context%converge
      DEALLOCATE(context)

      END PROGRAM
