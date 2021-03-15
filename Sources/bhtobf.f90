!*******************************************************************************
!>  @file siesta_currents.f90
!>  @brief Contains routines to convert B fields from half grids to full grids.
!
!  Note separating the Doxygen comment block here so detailed decription is
!  found in the Module not the file.
!
!>  Converts half grid B field to full grid and removes the jacobian factor.
!*******************************************************************************
      MODULE bhtobf

      IMPLICIT NONE

      CONTAINS

!*******************************************************************************
!  UTILITY SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  @brief Convert half mesh fields to the full mesh and remove jacobian factor.
!>
!>  On entry the full mesh quantities contain an extra jacobian factor that will
!>  be removed after to coversion to the full mesh.
!>
!>  @param[inout] bsupsijh Contravariant s component of B field on half mesh.
!>  @param[inout] bsupsijh Contravariant u component of B field on half mesh.
!>  @param[inout] bsupsijh Contravariant v component of B field on half mesh.
!>  @param[inout] bsupsijh Contravariant s component of B field on full mesh.
!>  @param[inout] bsupsijh Contravariant u component of B field on full mesh.
!>  @param[inout] bsupsijh Contravariant v component of B field on full mesh.
!>  @param[inout] pijh     Pressure on half mesh.
!>  @param[inout] pijh     Pressure on full mesh.
!-------------------------------------------------------------------------------
      SUBROUTINE bhalftobfull(bsupsijh, bsupuijh, bsupvijh,                    &
                              bsupsijf, bsupuijf, bsupvijf,                    &
                              pijh, pijf)
      USE stel_kinds
      USE v3_utilities, ONLY: assert_eq
      USE quantities, ONLY: ns, jacobf, jacobh
      USE timer_mod
      USE utilities, ONLY: to_full_mesh, set_bndy_half_to_full
      USE siesta_namelist, ONLY: l_vessel
      USE island_params, ONLY: fourier_context, mpol=>mpol_i, ntor=>ntor_i
      USE fourier
      USE shared_data, ONLY: lasym

      IMPLICIT NONE

!  Declare Arguments
      REAL (dp), DIMENSION(:,:,:), ALLOCATABLE, INTENT(inout) :: bsupsijh
      REAL (dp), DIMENSION(:,:,:), ALLOCATABLE, INTENT(inout) :: bsupuijh
      REAL (dp), DIMENSION(:,:,:), ALLOCATABLE, INTENT(inout) :: bsupvijh
      REAL (dp), DIMENSION(:,:,:), ALLOCATABLE, INTENT(inout) :: bsupsijf
      REAL (dp), DIMENSION(:,:,:), ALLOCATABLE, INTENT(inout) :: bsupuijf
      REAL (dp), DIMENSION(:,:,:), ALLOCATABLE, INTENT(inout) :: bsupvijf
      REAL (dp), DIMENSION(:,:,:), ALLOCATABLE, INTENT(inout) :: pijh
      REAL (dp), DIMENSION(:,:,:), ALLOCATABLE, INTENT(inout) :: pijf

!  local variables
      INTEGER                                                 :: nsmin
      INTEGER                                                 :: nsmax

!  Start of executable code
      nsmin = LBOUND(bsupsijh,3)
      nsmax = UBOUND(bsupsijh,3)

      bsupsijh = bsupsijh/jacobh(:,:,nsmin:nsmax)
      CALL to_full_mesh(bsupsijh, bsupsijf)
      IF (nsmin .EQ. 1) THEN
         bsupsijf(:,:,1) = bsupsijh(:,:,1)
      END IF

      bsupuijh = bsupuijh/jacobh(:,:,nsmin:nsmax)
      CALL to_full_mesh(bsupuijh, bsupuijf)
      IF (nsmin .EQ. 1) THEN
         bsupuijf(:,:,1) = bsupuijh(:,:,1)
      END IF

      bsupvijh = bsupvijh/jacobh(:,:,nsmin:nsmax)
      CALL to_full_mesh(bsupvijh, bsupvijf)
      IF (nsmin .EQ. 1) THEN
         bsupvijf(:,:,1) = bsupvijh(:,:,1)
      END IF

      CALL assert_eq(SIZE(bsupsijh,3), SIZE(pijh,3), 'bhtobf pijh SIZE WRONG')
      pijh = pijh/jacobh(:,:,nsmin:nsmax)
      CALL to_full_mesh(pijh, pijf)
      IF (nsmin .EQ. 1) THEN
         pijf(:,:,1) = pijh(:,:,1)
      END IF

      IF (nsmax .eq. ns .and. .not.l_vessel) THEN
         bsupsijf(:,:,ns) = 0
      END IF

      END SUBROUTINE

      END MODULE
