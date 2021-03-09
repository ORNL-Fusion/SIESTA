!*******************************************************************************
!>  @file siesta_run.f90
!>  @brief Contains module @ref siesta_run.
!
!  Note separating the Doxygen comment block here so detailed decription is
!  found in the Module not the file.
!
!>  This module contains all the code needed to define an error.
!*******************************************************************************
      MODULE siesta_error
      USE v3_utilities, ONLY: assert
      USE shared_data, Only: lverbose

      IMPLICIT NONE

!*******************************************************************************
!  module parameters
!*******************************************************************************
!>  All error flags off.
      INTEGER, PARAMETER :: siesta_error_no_error   = 0

!  Bit positions for the various error conditions.
!>  Assertion Error.
      INTEGER, PARAMETER :: siesta_error_assert     = 0
!>  Assertion Error.
      INTEGER, PARAMETER :: siesta_error_block_tri  = 1
!>  General Error.
      INTEGER, PARAMETER :: siesta_error_general    = 2
!>  Allocation Error.
      INTEGER, PARAMETER :: siesta_error_allocation = 3
!>  IO Error.
      INTEGER, PARAMETER :: siesta_error_io         = 4
!>  Parameter error.
      INTEGER, PARAMETER :: siesta_error_param      = 5

!*******************************************************************************
!  module variables.
!*******************************************************************************
      INTEGER :: siesta_error_state = siesta_error_no_error

      CONTAINS
!*******************************************************************************
!  SETTER SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  @brief Sets the bit flag for the error code.
!>
!>  Sets the bits for an error. Errors are defined by the module parameters. If
!>  this is reconstruction context, only report the error. Otherwise exit.
!>
!>  @param[in] error_code Error code to set.
!>  @param[in] message    Message to report.
!-------------------------------------------------------------------------------
      SUBROUTINE siesta_error_set_error(error_code, message)
      USE descriptor_mod, ONLY: iam

      IMPLICIT NONE

!  Declare Arguments
      INTEGER, INTENT(in)           :: error_code
      CHARACTER (len=*), INTENT(in) :: message

!  Start of executable code.
      siesta_error_state = IBSET(siesta_error_state, error_code)

      IF (iam .eq. 0 .and. lverbose) THEN
         WRITE (*,*) TRIM(message)
      END IF

      END SUBROUTINE

!-------------------------------------------------------------------------------
!>  @brief Clears the bit flag for the error code.
!>
!>  Clears the bits for an error. Errors are defined by the module parameters.
!>
!>  @param[in] error_code Error code to set.
!-------------------------------------------------------------------------------
      SUBROUTINE siesta_error_clear_error(error_code)

      IMPLICIT NONE

!  Declare Arguments
      INTEGER, INTENT(in) :: error_code

!  Start of executable code.
      siesta_error_state = IBCLR(siesta_error_state, error_code)

      END SUBROUTINE

!-------------------------------------------------------------------------------
!>  @brief Clears all error code.
!>
!>  Sets all bit positions to zero.
!-------------------------------------------------------------------------------
      SUBROUTINE siesta_error_clear_all

      IMPLICIT NONE

!  Start of executable code.
      siesta_error_state = siesta_error_no_error

      END SUBROUTINE

      END MODULE
