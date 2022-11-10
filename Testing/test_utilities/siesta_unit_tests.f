!*******************************************************************************
!>  @file siesta_unit_tests.f90
!>  @brief Contains main routines for SIESTA unit tests.
!
!  Note separating the Doxygen comment block here so detailed decription is
!  found in the Module not the file.
!
!>  These run the embeded units tests in SIESTA.
!*******************************************************************************
!*******************************************************************************
!  MAIN PROGRAM
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  @brief SIESTA Unit tests main program.
!>
!>  Highest level SIESTA unit test routine.
!-------------------------------------------------------------------------------
      PROGRAM siesta_unit_tests
      USE fourier

      IMPLICIT NONE

!  Local variables
      LOGICAL :: test_passed

!  Start of executable code.
      test_passed = test_fourier()
      IF (.not.test_passed) THEN
         WRITE (*,*) 'test_fourier failed.'
      END IF

      IF (.not.test_passed) THEN
         CALL exit(1)
      END IF

      END PROGRAM
