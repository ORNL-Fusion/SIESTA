!*******************************************************************************
!>  @file utilities.f90
!>  @brief Contains module @ref utilities.
!
!  Note separating the Doxygen comment block here so detailed decription is
!  found in the Module not the file.
!
!>  This file contains utilities for converting the siesta grids.
!*******************************************************************************
      MODULE utilities
      USE stel_kinds
      USE island_params, ONLY: ohs=>ohs_i, ns=>ns_i, mpol=>mpol_i,             &
                               ntor=>ntor_i, nfp=>nfp_i, fourier_context
      USE v3_utilities, ONLY: assert, assert_eq
      USE fourier, ONLY: m0, m1, m2, n0, n1, f_sin, f_cos, f_ds, f_none,       &
                         b_ds, b_jac, b_con, f_jac

      IMPLICIT NONE

!*******************************************************************************
!  INTERFACE BLOCKS
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  Interface to convert between half and full mesh using either
!>  @ref to_full_mesh2 or @ref to_full_mesh3
!-------------------------------------------------------------------------------
      INTERFACE to_full_mesh
         MODULE PROCEDURE to_full_mesh2, to_full_mesh3
      END INTERFACE

!-------------------------------------------------------------------------------
!>  Interface to convert between full and half mesh using either
!>  @ref to_half_mesh or @ref to_half_mesh_p
!-------------------------------------------------------------------------------
      INTERFACE to_half_mesh
         MODULE PROCEDURE to_half_mesh, to_half_mesh_p
      END INTERFACE

!-------------------------------------------------------------------------------
!>  Interface to take the gradianet from full to half mesh using either
!>  @ref GradientHalf or @ref GradientHalf_p
!-------------------------------------------------------------------------------
      INTERFACE GradientHalf
         MODULE PROCEDURE GradientHalf, GradientHalf_p
      END INTERFACE

!-------------------------------------------------------------------------------
!>  Interface to take the curl on a half mesh quantity and return a full mesh
!>  using either @ref curl_htof or @ref curl_htof_targ
!-------------------------------------------------------------------------------
      INTERFACE curl_htof
         MODULE PROCEDURE curl_htof, curl_htof_targ
      END INTERFACE

!  FIXME: Move to fourier
!-------------------------------------------------------------------------------
!>  Interface to set the fouier conditions on a quantity for either scalar
!>  @ref set_bndy_fouier_m0 or vector @ref set_bndy_fouier_m0_vec quantities.
!-------------------------------------------------------------------------------
      INTERFACE set_bndy_fouier_m0
         MODULE PROCEDURE set_bndy_fouier_m0, set_bndy_fouier_m0_vec
      END INTERFACE

      CONTAINS

!*******************************************************************************
!  UTILITY SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  @brief Calculate the gradient to the half mesh from a full mesh quantity.
!>
!>  The arrays must be passed in as ALLOCATABLE to preserve the array bounds.
!>
!>  @param[inout] gradienth Gradient on the half grid.
!>  @param[in]    vecf      Full grid quantity.
!-------------------------------------------------------------------------------
      SUBROUTINE GradientHalf(gradienth, vecf)

      IMPLICIT NONE

!  Declare Arguments
      REAL (dp), DIMENSION(:,:,:), ALLOCATABLE, INTENT(inout) :: gradienth
      REAL (dp), DIMENSION(:,:,:), ALLOCATABLE, INTENT(in)    :: vecf

!  local variables
      INTEGER                                                 :: nsmin
      INTEGER                                                 :: nsmax

!  Start of executable code
      nsmin = LBOUND(gradienth,3)
      nsmax = UBOUND(gradienth,3)

      gradienth(:,:,nsmin + 1:nsmax) = ohs*(vecf(:,:,nsmin + 1:nsmax)          &
                                     -      vecf(:,:,nsmin:nsmax - 1))
      gradienth(:,:,nsmin) = 0

      END SUBROUTINE

!-------------------------------------------------------------------------------
!>  @brief Calculate the gradient to the half mesh from a full mesh quantity.
!>
!>  The arrays must be passed in as ALLOCATABLE to preserve the array bounds.
!>  Full vector is passed in as a pointer.
!>
!>  @param[inout] gradienth Gradient on the half grid.
!>  @param[in]    vecf      Full grid quantity.
!-------------------------------------------------------------------------------
      SUBROUTINE GradientHalf_p(gradienth, vecf)

      IMPLICIT NONE

!  Declare Arguments
      REAL (dp), DIMENSION(:,:,:), ALLOCATABLE, INTENT(inout) :: gradienth
      REAL (dp), DIMENSION(:,:,:), POINTER                    :: vecf

!  local variables
      INTEGER                                                 :: nsmin
      INTEGER                                                 :: nsmax

!  Start of executable code
      nsmin = LBOUND(gradienth,3)
      nsmax = UBOUND(gradienth,3)

      gradienth(:,:,nsmin + 1:nsmax) = ohs*(vecf(:,:,nsmin + 1:nsmax)          &
                                     -      vecf(:,:,nsmin:nsmax - 1))
      gradienth(:,:,nsmin) = 0

      END SUBROUTINE

!-------------------------------------------------------------------------------
!>  @brief Calculate the gradient to the full mesh from a half mesh quantity.
!>
!>  The arrays must be passed in as ALLOCATABLE to preserve the array bounds.
!>
!>  @param[inout] gradientf Gradient on the full grid.
!>  @param[in]    vech      Half grid quantity.
!-------------------------------------------------------------------------------
      SUBROUTINE GradientFull(gradientf, vech)

      IMPLICIT NONE

!  Declare Arguments
      REAL (dp), DIMENSION(:,:,:), ALLOCATABLE, INTENT(inout) :: gradientf
      REAL (dp), DIMENSION(:,:,:), ALLOCATABLE, INTENT(in)    :: vech

!  local variables
      INTEGER                                                 :: nsmin
      INTEGER                                                 :: nsmax

!  Start of executable code
      nsmin = LBOUND(gradientf,3)
      nsmax = UBOUND(gradientf,3)

      gradientf(:,:,nsmin:nsmax - 1) = ohs*(vech(:,:,nsmin + 1:nsmax)          &
                                     -      vech(:,:,nsmin:nsmax - 1))

      IF (UBOUND(gradientf,3) .ge. nsmax) THEN
         gradientf(:,:,nsmax:) = 0
      END IF

      IF (nsmin .eq. 1) THEN
         gradientf(:,:,1) = 0
      END IF

      IF (nsmax .eq. ns) THEN
         gradientf(:,:,ns) = gradientf(:,:,ns - 1)
      END IF

      END SUBROUTINE

!-------------------------------------------------------------------------------
!>  @brief Migrates a full mesh quantity to the half mesh.
!>
!>  The arrays must be passed in as ALLOCATABLE to preserve the array bounds.
!>
!>  @param[in]    vecf Full grid quantity.
!>  @param[inout] vech Half grid quantity.
!-------------------------------------------------------------------------------
      SUBROUTINE to_half_mesh(vecf, vech)

      IMPLICIT NONE

!  Declare Arguments
      REAL (dp), DIMENSION(:,:,:), ALLOCATABLE, INTENT(in)    :: vecf
      REAL (dp), DIMENSION(:,:,:), ALLOCATABLE, INTENT(inout) :: vech

!  local variables
      INTEGER                                                 :: nsmin
      INTEGER                                                 :: nsmax

!  Start of executable code
      nsmin = LBOUND(vech,3)
      nsmax = UBOUND(vech,3)

      vech(:,:,nsmin + 1:nsmax) = (vecf(:,:,nsmin:nsmax - 1)                   &
                                +  vecf(:,:,nsmin + 1:nsmax))*0.5
      vech(:,:,nsmin) = 0

      END SUBROUTINE

!-------------------------------------------------------------------------------
!>  @brief Migrates a full mesh quantity to the half mesh.
!>
!>  The arrays must be passed in as ALLOCATABLE to preserve the array bounds.
!>  Full vector grid is passed in as pointer.
!>
!>  @param[in]    vecf Full grid quantity.
!>  @param[inout] vech Half grid quantity.
!-------------------------------------------------------------------------------
      SUBROUTINE to_half_mesh_p(vecf, vech)

      IMPLICIT NONE

!  Declare Arguments
      REAL (dp), DIMENSION(:,:,:), POINTER                    :: vecf
      REAL (dp), DIMENSION(:,:,:), ALLOCATABLE, INTENT(inout) :: vech

!  local variables
      INTEGER                                                 :: nsmin
      INTEGER                                                 :: nsmax

!  Start of executable code
      nsmin = LBOUND(vech,3)
      nsmax = UBOUND(vech,3)
      vech(:,:,nsmin + 1:nsmax) = (vecf(:,:,nsmin:nsmax - 1)                   &
                                +  vecf(:,:,nsmin + 1:nsmax))*0.5
      vech(:,:,nsmin) = 0

      END SUBROUTINE

!-------------------------------------------------------------------------------
!>  @brief Migrate a three dimensional quantity form the half mesh to the full
!>         mesh.
!>
!>  First and last point are extrapolated linearly. The arrays must be passed in
!>  as ALLOCATABLE to preserve the array bounds.
!>
!>  @param[in]    vech Half mesh qantity.
!>  @param[inout] vecf Full mesh quantity.
!-------------------------------------------------------------------------------
      SUBROUTINE to_full_mesh3(vech, vecf)

      IMPLICIT NONE

!  Declare Arguments
      REAL (dp), DIMENSION(:,:,:), ALLOCATABLE, INTENT(in)    :: vech
      REAL (dp), DIMENSION(:,:,:), ALLOCATABLE, INTENT(inout) :: vecf

!  local variables
      INTEGER                                                 :: nsmin
      INTEGER                                                 :: nsmax

!  Start of executable code
      nsmin = LBOUND(vecf, 3)
      CALL assert(nsmin .ge. LBOUND(vech,3), 'LBOUND WRONG IN to_full_mesh3')
      nsmax = UBOUND(vecf, 3)
      CALL assert(nsmax .le. UBOUND(vech,3), 'UBOUND WRONG IN to_full_mesh3')

      vecf(:,:,nsmin:nsmax - 1) = 0.5*(vech(:,:,nsmin + 1:nsmax)               &
                                +      vech(:,:,nsmin:nsmax - 1))
      vecf(:,:,nsmax) = 0                              
      IF (nsmin .eq. 1)  THEN
         vecf(:,:,1) = vech(:,:,2)
      END IF
      IF (nsmax .eq. ns) THEN
         vecf(:,:,ns) = vech(:,:,ns)
      END IF

      END SUBROUTINE

!-------------------------------------------------------------------------------
!>  @brief Migrate a two dimensional quantity form the half mesh to the full
!>         mesh.
!>
!>  First and last point are extrapolated linearly. The arrays must be passed in
!>  as ALLOCATABLE to preserve the array bounds.
!>
!>  @param[in]    vech Half mesh qantity.
!>  @param[inout] vecf Full mesh quantity.
!-------------------------------------------------------------------------------
      SUBROUTINE to_full_mesh2(vech, vecf)

      IMPLICIT NONE

!  Declare Arguments
      REAL (dp), DIMENSION(:,:), ALLOCATABLE, INTENT(in)    :: vech
      REAL (dp), DIMENSION(:,:), ALLOCATABLE, INTENT(inout) :: vecf

!  local variables
      INTEGER                                               :: nsmin
      INTEGER                                               :: nsmax

!  Start of executable code
      nsmin = LBOUND(vecf,2)
      CALL ASSERT(nsmin .ge. LBOUND(vech,2), 'LBOUND WRONG IN to_full_mesh2')
      nsmax = UBOUND(vecf,2)
      CALL ASSERT(nsmax .le. UBOUND(vech,2), 'UBOUND WRONG IN to_full_mesh2')

      vecf(:,nsmin:nsmax - 1) = 0.5*(vech(:,nsmin:nsmax - 1)                   &
                              +      vech(:,nsmin + 1:nsmax))

      vecf(:,nsmax) = 0                              
      IF (nsmin .EQ. 1)  THEN
         vecf(:,1) = vech(:,2)
      END IF
      IF (nsmax .EQ. ns) THEN
         vecf(:,ns)= vech(:,ns)
      END IF

      END SUBROUTINE

!  FIXME: Move to a new curl module.
!-------------------------------------------------------------------------------
!>  @brief Compute curl of the full grid quantity a on the half grid.
!>
!>  The contravariant components of jac*Curl A are equal to
!>
!>    sqrt(g)*B^k = dA_j/du^i - dA_i/du^j                                    (j)
!>
!>  for cyclic permutations of i,j,k. This routine does not remove the jacobian
!>  factor. The arrays must be passed in as ALLOCATABLE or POINTER to preserve
!>  the array bounds.
!>
!>  @param[in]    asubsmnf  Covariant vector in the s direction full grid.
!>  @param[in]    asubsmnf  Covariant vector in the u direction full grid.
!>  @param[in]    asubsmnf  Covariant vector in the v direction full grid.
!>  @param[inout] jbsubsmnh Contravariant vector in the s direction with
!>                          jacobian factor half grid.
!>  @param[inout] jbsubsmnh Contravariant vector in the u direction with
!>                          jacobian factor half grid.
!>  @param[inout] jbsubsmnh Contravariant vector in the v direction with
!>                          jacobian factor half grid.
!>  @param[in]    parity    Parity flag.
!>  @param[in]    nsmin     Minimum radial index.
!>  @param[in]    nsmax     Maximum radial index.
!-------------------------------------------------------------------------------
      SUBROUTINE curl_ftoh(asubsmnf, asubumnf, asubvmnf,                       &
                           jbsupsmnh, jbsupumnh, jbsupvmnh, parity,            &
                           nsmin, nsmax)
      USE nscalingtools, ONLY: startglobrow

      IMPLICIT NONE

!  Declare Arguments
      REAL (dp), DIMENSION(:,:,:), POINTER                    :: asubsmnf
      REAL (dp), DIMENSION(:,:,:), POINTER                    :: asubumnf
      REAL (dp), DIMENSION(:,:,:), POINTER                    :: asubvmnf
      REAL (dp), DIMENSION(:,:,:), ALLOCATABLE, INTENT(inout) :: jbsupsmnh
      REAL (dp), DIMENSION(:,:,:), ALLOCATABLE, INTENT(inout) :: jbsupumnh
      REAL (dp), DIMENSION(:,:,:), ALLOCATABLE, INTENT(inout) :: jbsupvmnh
      INTEGER, INTENT(IN)                                     :: parity
      INTEGER, INTENT(IN)                                     :: nsmin
      INTEGER, INTENT(IN)                                     :: nsmax

!  Local variables
      INTEGER                                                 :: s
      INTEGER                                                 :: m
      INTEGER                                                 :: n
      INTEGER                                                 :: mp
      INTEGER                                                 :: np
      INTEGER                                                 :: sparity
      INTEGER                                                 :: fouruv
      INTEGER                                                 :: istat
      INTEGER                                                 :: nmin
      REAL (dp), ALLOCATABLE, DIMENSION(:,:,:)                :: asubsmnh
      REAL (dp), ALLOCATABLE, DIMENSION(:,:,:)                :: asubumnh
      REAL (dp), ALLOCATABLE, DIMENSION(:,:,:)                :: asubvmnh
      REAL (dp), ALLOCATABLE, DIMENSION(:,:,:)                :: da_uds
      REAL (dp), ALLOCATABLE, DIMENSION(:,:,:)                :: da_vds
LOGICAL :: test
!  Start of executable code
      CALL assert(nsmax .le. UBOUND(jbsupsmnh,3), 'UBOUND wrong in curl_ftoh')

      IF (parity .EQ. f_sin) THEN
         sparity = 1
         fouruv = f_cos
      ELSE
         sparity = -1
         fouruv = f_sin
      END IF

      ALLOCATE(asubsmnh(0:mpol,-ntor:ntor,nsmin:nsmax),                        &
               asubumnh(0:mpol,-ntor:ntor,nsmin:nsmax),                        &
               asubvmnh(0:mpol,-ntor:ntor,nsmin:nsmax), stat=istat)
      CALL assert_eq(0, istat, 'Allocation1 failed in curl_ftoh')

      CALL to_half_mesh(asubsmnf, asubsmnh)                       
      CALL to_half_mesh(asubumnf, asubumnh)
      CALL to_half_mesh(asubvmnf, asubvmnh)

      CALL set_bndy_fouier_m0(asubsmnh, asubumnh, asubvmnh, parity)

!  Compute sqrt(g)B^X contravariant B-field components in Fourier space
      ALLOCATE(da_uds(0:mpol,-ntor:ntor,nsmin:nsmax),                          &
               da_vds(0:mpol,-ntor:ntor,nsmin:nsmax), stat=istat)
      CALL assert_eq(0, istat, 'Allocation2 failed in curl_ftoh')

!  Radial gradients of asubmnu,v on half grid excluding origin and edge. Origin
!  and edge will be treated separately.
      CALL GradientHalf(da_uds, asubumnf)
      CALL GradientHalf(da_vds, asubvmnf)
      CALL set_bndy_fouier_m0(da_uds, fouruv)
      CALL set_bndy_fouier_m0(da_vds, fouruv)

      nmin = MAX(2, startglobrow)
      CALL assert(nmin .ge. LBOUND(jbsupsmnh,3), 'LBOUND wrong in curl_ftoh')

!  (sqrt(g)*B^X)mn
      DO s = nmin, nsmax
         DO n = -ntor, ntor
            np = sparity*fourier_context%tor_modes(n)*nfp
            DO m = 0, mpol
               mp = sparity*m
               jbsupsmnh(m,n,s) = np*asubumnh(m,n,s) - mp*asubvmnh(m,n,s)
               jbsupumnh(m,n,s) = np*asubsmnh(m,n,s) - da_vds(m,n,s)
               jbsupvmnh(m,n,s) = da_uds(m,n,s) - mp*asubsmnh(m,n,s)
            END DO
         END DO
      END DO

      DEALLOCATE (asubsmnh, asubumnh, asubvmnh, da_vds, da_uds, stat=istat)   
      CALL assert_eq(0, istat, 'Deallocation failed in curl_ftoh')

      IF (startglobrow .eq. 1) THEN
         jbsupsmnh(:,:,1) = 0
         jbsupumnh(:,:,1) = 0
         jbsupvmnh(:,:,1) = 0
      END IF

      CALL set_bndy_fouier_m0(jbsupsmnh, jbsupumnh, jbsupvmnh, parity)

!test = test_div(jbsupsmnh, jbsupumnh, jbsupvmnh, parity, nsmin, nsmax)

      END SUBROUTINE

!  FIXME: Move to a new curl module.
!-------------------------------------------------------------------------------
!>  @brief Compute curl of the half grid quantity a on the full grid.
!>
!>  The contravariant components of jac*Curl A are equal to
!>
!>    sqrt(g)*K^k = dB_j/du^i - dB_i/du^j                                    (j)
!>
!>  for cyclic permutations of i,j,k. This routine does not remove the jacobian
!>  factor. The arrays must be passed in as ALLOCATABLE or POINTER to preserve
!>  the array bounds.
!>
!>  @param[in]    bsubsmnh  Covariant vector in the s direction half grid.
!>  @param[in]    bsubsmnh  Covariant vector in the u direction half grid.
!>  @param[in]    bsubsmnh  Covariant vector in the v direction half grid.
!>  @param[inout] jKsubsmnf Contravariant vector in the s direction with
!>                          jacobian factor full grid.
!>  @param[inout] jKsubsmnf Contravariant vector in the u direction with
!>                          jacobian factor full grid.
!>  @param[inout] jKsubsmnf Contravariant vector in the v direction with
!>                          jacobian factor full grid.
!>  @param[in]    parity    Parity flag.
!>  @param[in]    nsmin     Minimum radial index.
!>  @param[in]    nsmax     Maximum radial index.
!>  @param[in]    nsend     Ending row index.
!>  @param[inout] curtor    Current enclosed at the boundary.
!-------------------------------------------------------------------------------
      SUBROUTINE curl_htof(bsubsmnh, bsubumnh, bsubvmnh,                       &
                           jksupsmnf, jksupumnf, jksupvmnf, parity,            &
                           nsmin, nsmax, nsend, curtor)
      USE stel_constants

      IMPLICIT NONE

!  Declare Arguments
      REAL (dp), DIMENSION(:,:,:), ALLOCATABLE, INTENT(in) :: bsubsmnh
      REAL (dp), DIMENSION(:,:,:), ALLOCATABLE, INTENT(in) :: bsubumnh
      REAL (dp), DIMENSION(:,:,:), ALLOCATABLE, INTENT(in) :: bsubvmnh
      REAL (dp), DIMENSION(:,:,:), POINTER                 :: jksupsmnf
      REAL (dp), DIMENSION(:,:,:), POINTER                 :: jksupumnf
      REAL (dp), DIMENSION(:,:,:), POINTER                 :: jksupvmnf
      INTEGER, INTENT(in)                                  :: parity
      INTEGER, INTENT(in)                                  :: nsmin
      INTEGER, INTENT(in)                                  :: nsmax
      INTEGER, INTENT(in)                                  :: nsend
      REAL (dp), INTENT(inout)                             :: curtor

!  Local variables
      INTEGER                                              :: m
      INTEGER                                              :: n
      INTEGER                                              :: mp
      INTEGER                                              :: np
      INTEGER                                              :: mo
      INTEGER                                              :: no
      INTEGER                                              :: moff
      INTEGER                                              :: noff
      INTEGER                                              :: sparity
      INTEGER                                              :: istat
      INTEGER                                              :: nmax
      REAL (dp), DIMENSION(:,:,:), ALLOCATABLE             :: bsubsmnf
      REAL (dp), DIMENSION(:,:,:), ALLOCATABLE             :: bsubumnf
      REAL (dp), DIMENSION(:,:,:), ALLOCATABLE             :: bsubvmnf
      REAL (dp), DIMENSION(:,:,:), ALLOCATABLE             :: db_uds
      REAL (dp), DIMENSION(:,:,:), ALLOCATABLE             :: db_vds
      INTEGER                                              :: fours
      INTEGER                                              :: fouruv

!  Start of executable code.
      CALL assert(nsmin .ge. LBOUND(jksupsmnf,3), 'LBOUND wrong in curl_htof')

      IF (parity .eq. f_sin) THEN
         sparity = 1
         fours = f_sin
         fouruv = f_cos
      ELSE
         sparity = -1
         fours = f_cos
         fouruv = f_sin
      END IF

      moff = LBOUND(jksupsmnf,1)
      noff = LBOUND(jksupsmnf,2) + ntor

      ALLOCATE(bsubsmnf(0:mpol,-ntor:ntor,nsmin:nsmax),                        &
               bsubumnf(0:mpol,-ntor:ntor,nsmin:nsmax),                        &
               bsubvmnf(0:mpol,-ntor:ntor,nsmin:nsmax), stat=istat)
      CALL assert_eq(0, istat, 'Allocation1 failed in curl_htof')

      CALL to_full_mesh(bsubsmnh, bsubsmnf)
      CALL to_full_mesh(bsubumnh, bsubumnf)
      CALL to_full_mesh(bsubvmnh, bsubvmnf)

      CALL set_bndy_half_to_full(bsubsmnh, bsubsmnf, nsmin, fours, f_ds)
      CALL set_bndy_half_to_full(bsubumnh, bsubumnf, nsmin, fouruv, f_none)
      CALL set_bndy_half_to_full(bsubvmnh, bsubvmnf, nsmin, fouruv, f_none)

!  Store the total enclosed current to check against the orginal VMEC current.
      IF (parity .eq. f_sin .and. nsmax .eq. ns) THEN
         curtor = twopi*bsubumnf(m0,n0,ns)/mu0
      END IF

!  Compute sqrt(g)B^X == A^X contravariant current components in Fourier space.
      ALLOCATE(db_uds(0:mpol,-ntor:ntor,nsmin:nsmax),                          &
               db_vds(0:mpol,-ntor:ntor,nsmin:nsmax), stat=istat)
      CALL assert_eq(0, istat, 'Allocation2 failed in curl_htof')

!  Radial gradients of asubmnu,v on full grid excluding origin and edge. Origin
!  and edge will be treated separately.
      CALL GradientFull(db_uds, bsubumnh)
      CALL GradientFull(db_vds, bsubvmnh)

      CALL set_bndy_half_to_full_ds(bsubumnh, db_uds, nsmin, fouruv, f_none)
      CALL set_bndy_half_to_full_ds(bsubvmnh, db_vds, nsmin, fouruv, f_none)

      nmax = MIN(ns, nsend)
      CALL assert(nmax .le. UBOUND(jksupsmnf,3), 'UBOUND wrong in curl_htof')

!  (sqrt(g)*J^X)mn
      DO n = -ntor, ntor
         np = sparity*fourier_context%tor_modes(n)*nfp
         no = n + noff
         DO m = 0, mpol
            mp = sparity*m
            mo = m + moff
            jksupsmnf(mo,no,nsmin:nmax) = np*bsubumnf(m,n,nsmin:nmax)          &
                                        - mp*bsubvmnf(m,n,nsmin:nmax)
            jksupumnf(mo,no,nsmin:nmax) = np*bsubsmnf(m,n,nsmin:nmax)          &
                                        - db_vds(m,n,nsmin:nmax)
            jksupvmnf(mo,no,nsmin:nmax) = db_uds(m,n,nsmin:nmax)               &
                                        - mp*bsubsmnf(m,n,nsmin:nmax)
         END DO
      END DO

      DEALLOCATE(bsubsmnf, bsubumnf, bsubvmnf, db_uds, db_vds)

!  Avoid fp failure in debug mode. The in comming values are one extra radial
!  index than the out going. Zero these extended terms out.
      jksupsmnf(:,:,nmax+1:UBOUND(jksupsmnf,3)) = 0
      jksupumnf(:,:,nmax+1:UBOUND(jksupsmnf,3)) = 0
      jksupvmnf(:,:,nmax+1:UBOUND(jksupsmnf,3)) = 0

      jksupsmnf(:,:,LBOUND(jksupsmnf,3):nsmin-1) = 0
      jksupumnf(:,:,LBOUND(jksupsmnf,3):nsmin-1) = 0
      jksupvmnf(:,:,LBOUND(jksupsmnf,3):nsmin-1) = 0

      CALL set_bndy_fouier_m0(jksupsmnf, jksupumnf, jksupvmnf, parity)
      IF (nsmin .eq. 1) THEN
         CALL set_bndy_full_origin(jksupsmnf, jksupumnf, jksupvmnf, f_jac)
      END IF

      END SUBROUTINE

!  FIXME: Move to a new curl module.
!-------------------------------------------------------------------------------
!>  @brief Convert allocatable to pointer for use in @ref curl_htof
!>
!>  @param[in]    Asubsmnh  Covariant vector in the s direction half grid.
!>  @param[in]    Asubsmnh  Covariant vector in the u direction half grid.
!>  @param[in]    Asubsmnh  Covariant vector in the v direction half grid.
!>  @param[inout] jBsubsmnf Contravariant vector in the s direction with
!>                          jacobian factor full grid.
!>  @param[inout] jBsubsmnf Contravariant vector in the u direction with
!>                          jacobian factor full grid.
!>  @param[inout] jBsubsmnf Contravariant vector in the v direction with
!>                          jacobian factor full grid.
!>  @param[in]    parity    Parity flag.
!>  @param[in]    nsmin     Minimum radial index.
!>  @param[in]    nsmax     Maximum radial index.
!>  @param[in]    nsend     Ending row index.
!>  @param[inout] curtor    Current enclosed at the boundary.
!-------------------------------------------------------------------------------
      SUBROUTINE curl_htof_targ(bsubsmnh, bsubumnh, bsubvmnh,                  &
                                jksupsmnf, jksupumnf, jksupvmnf,               &
                                parity, nsmin, nsmax, nsend, curtor)

!  Declare Arguments
      REAL (dp), DIMENSION(:,:,:), ALLOCATABLE, INTENT(in) :: bsubsmnh
      REAL (dp), DIMENSION(:,:,:), ALLOCATABLE, INTENT(in) :: bsubumnh
      REAL (dp), DIMENSION(:,:,:), ALLOCATABLE, INTENT(in) :: bsubvmnh
      REAL (dp), DIMENSION(:,:,:), ALLOCATABLE, TARGET, INTENT(inout) ::       &
         jksupsmnf
      REAL (dp), DIMENSION(:,:,:), ALLOCATABLE, TARGET, INTENT(inout) ::       &
         jksupumnf
      REAL (dp), DIMENSION(:,:,:), ALLOCATABLE, TARGET, INTENT(inout) ::       &
         jksupvmnf
      INTEGER, INTENT(IN)                                  :: parity
      INTEGER, INTENT(IN)                                  :: nsmin
      INTEGER, INTENT(IN)                                  :: nsmax
      INTEGER, INTENT(in)                                  :: nsend
      REAL (dp), INTENT(inout)                             :: curtor

!  Local variables
      REAL (dp), POINTER, DIMENSION(:,:,:)                 :: jksp
      REAL (dp), POINTER, DIMENSION(:,:,:)                 :: jkup
      REAL (dp), POINTER, DIMENSION(:,:,:)                 :: jkvp

!  Start of executable code
      jksp => jksupsmnf
      jkup => jksupumnf
      jkvp => jksupvmnf
      CALL curl_htof(bsubsmnh, bsubumnh, bsubvmnh, jksp, jkup, jkvp,           &
                     parity, nsmin, nsmax, nsend, curtor)

      END SUBROUTINE

!  FIXME: Move to fourier
!-------------------------------------------------------------------------------
!>  @brief Set fouier boundary conditions for full grid quantity converted from
!>         half grid.
!>
!>  At the center of the grid, the quantity should not be a function of u so all
!>  m > 0 modes should be zero. For the m = 0 terms, use the value from the half
!>  grid.
!>
!>  @param[in]    amnh   Half grid fourier quantity.
!>  @param[inout] amnf   Full grid fourier quantity.
!>  @param[in]    nsmin  Minimum radial index.
!>  @param[in]    parity Parity flag.
!>  @param[in]    flags  Control flags.
!-------------------------------------------------------------------------------
      PURE SUBROUTINE set_bndy_half_to_full(amnh, amnf, nsmin, parity, flags)

      IMPLICIT NONE

!  Declare Arguments
      REAL (dp), DIMENSION(:,:,:), INTENT(in)    :: amnh
      REAL (dp), DIMENSION(:,:,:), INTENT(inout) :: amnf
      INTEGER, INTENT(IN)                        :: nsmin
      INTEGER, INTENT(in)                        :: parity
      INTEGER, INTENT(in)                        :: flags

!  Local Variables
      INTEGER                                    :: moffh
      INTEGER                                    :: mofff

!  Start of executable code
      IF (nsmin .eq. 1) THEN
         moffh = LBOUND(amnh, 1)
         mofff = LBOUND(amnf, 1)

         amnf(:,:,1) = 0
         IF (BTEST(flags, b_ds)) THEN
            amnf(m1 + mofff,:,1) = amnh(m1 + moffh,:,2)
         ELSE
            amnf(m0 + mofff,:,1) = amnh(m0 + moffh,:,2)
         END IF
      END IF

      CALL set_bndy_fouier_m0(amnf, parity)

      END SUBROUTINE

!  FIXME: Move to fourier
!-------------------------------------------------------------------------------
!>  @brief Set fouier boundary conditions for full grid radial derivative
!>         quantity converted from half grid.
!>
!>  At the center of the grid the the m=1 quanity flips size across the s=0
!>  point. The results in a doubling of the m=1 component for the first half
!>  grid.
!>
!>  @param[in]    amnh    Half grid fourier quantity.
!>  @param[inout] amnf_ds Full grid radial derivative fourier quantity.
!>  @param[in]    nsmin   Minimum radial index.
!>  @param[in]    parity  Parity flag.
!>  @param[in]    flags   Control flags.
!-------------------------------------------------------------------------------
      PURE SUBROUTINE set_bndy_half_to_full_ds(amnh, amnf_ds, nsmin, parity,   &
                                               flags)

      IMPLICIT NONE

!  Declare Arguments
      REAL (dp), DIMENSION(:,:,:), INTENT(in)    :: amnh
      REAL (dp), DIMENSION(:,:,:), INTENT(inout) :: amnf_ds
      INTEGER, INTENT(IN)                        :: nsmin
      INTEGER, INTENT(in)                        :: parity
      INTEGER, INTENT(in)                        :: flags

!  Local Variables
      INTEGER                                    :: moffh
      INTEGER                                    :: mofff

!  Start of executable code
      IF (nsmin .eq. 1) THEN
         moffh = LBOUND(amnh, 1)
         mofff = LBOUND(amnf_ds, 1)

         amnf_ds(:,:,1) = 0.0

         IF (BTEST(flags, b_ds)) THEN
            amnf_ds(m0 + mofff,:,1) = 2.0*ohs*amnh(m0 + moffh,:,2)
         ELSE
            amnf_ds(m1 + mofff,:,1) = 2.0*ohs*amnh(m1 + moffh,:,2)
         END IF
      END IF

      CALL set_bndy_fouier_m0(amnf_ds, parity)

      END SUBROUTINE

!  FIXME: Move to fourier
!-------------------------------------------------------------------------------
!>  @brief Apply origin vector boundary conditions.
!>
!>  At the center of the grid, a vector quantity should have only m=1 modes in
!>  the s direction, no modes in the u direction and only m=0 modes in the v
!>  direction. See @ref siesta_init::GetFields.
!>
!>  @param[inout] amnsf Full grid fourier quantity for s component.
!>  @param[inout] amnuf Full grid fourier quantity for u component.
!>  @param[inout] amnvf Full grid fourier quantity for v component.
!>  @param[in]    flags Control flags.
!-------------------------------------------------------------------------------
      PURE SUBROUTINE set_bndy_full_origin(amnsf, amnuf, amnvf, flags)

      IMPLICIT NONE

!  Declare Arguments
      REAL (dp), DIMENSION(:,:,:), INTENT(inout) :: amnsf
      REAL (dp), DIMENSION(:,:,:), INTENT(inout) :: amnuf
      REAL (dp), DIMENSION(:,:,:), INTENT(inout) :: amnvf
      INTEGER, INTENT(in)                        :: flags

!  Local Variables
      INTEGER                                    :: moff

!  Start of executable code
      moff = LBOUND(amnsf, 1)

      IF (BTEST(flags, b_jac)) THEN
         amnsf(:,:,1) = 0
         amnuf(:,:,1) = 0
         amnvf(:,:,1) = 0
      ELSE
         amnsf(m0 + moff,:,1) = 0
         IF (m2 + moff .le. UBOUND(amnsf, 1)) THEN
            amnsf(m2 + moff:,:,1) = 0
         END IF

         IF (BTEST(flags, b_con)) THEN
            amnuf(m0 + moff:,:,1) = 0
            IF (m2 + moff .le. UBOUND(amnuf, 1)) THEN
               amnuf(m2 + moff:,:,1) = 0
            END IF
         ELSE
            amnuf(:,:,1) = 0
         END IF

         IF (BTEST(flags, b_con)) THEN
            amnvf(:,:,1) = 0
         ELSE
            amnvf(m1 + moff:,:,1) = 0
         END IF
      END IF

      END SUBROUTINE

!  FIXME: Move to fourier
!-------------------------------------------------------------------------------
!>  @brief Set vector fourier conditions.
!>
!>  For fouier modes for m=0 and for the negative n modes. The n=0 modes are
!>  zero for quantities with sin parity.
!>
!>  @param[inout] amn    Fourier quantity.
!>  @param[in]    parity Fouier parity flag.
!-------------------------------------------------------------------------------
      PURE SUBROUTINE set_bndy_fouier_m0(amn, parity)

      IMPLICIT NONE

!  Declare Arguments
      REAL (dp), DIMENSION(:,:,:), INTENT(inout) :: amn
      INTEGER, INTENT(in)                        :: parity

!  Local Variables
      INTEGER                                    :: moff
      INTEGER                                    :: noff

!  Start of executable code
      moff = LBOUND(amn, 1)
      noff = LBOUND(amn, 2) + ntor

      IF (parity .eq. f_cos) THEN
         amn(m0 + moff,-ntor+noff:-n1+noff,:) = 0
      ELSE
         amn(m0 + moff,-ntor+noff:n0+noff,:) = 0
      END IF

      END SUBROUTINE

!  FIXME: Move to fourier
!-------------------------------------------------------------------------------
!>  @brief Set vector fourier conditions.
!>
!>  For fouier modes for m=0 and for the negative n modes. The n=0 modes are
!>  zero for quantities with sin parity.
!>
!>  @param[inout] amns   Fourier quantity for s component.
!>  @param[inout] amnu   Fourier quantity for u component.
!>  @param[inout] amnv   Fourier quantity for v component.
!>  @param[in]    parity Fouier parity flag.
!-------------------------------------------------------------------------------
      PURE SUBROUTINE set_bndy_fouier_m0_vec(amns, amnu, amnv, parity)
      USE fourier, ONLY: f_cos, f_sin

      IMPLICIT NONE

!  Declare Arguments
      REAL (dp), DIMENSION(:,:,:), INTENT(inout) :: amns
      REAL (dp), DIMENSION(:,:,:), INTENT(inout) :: amnu
      REAL (dp), DIMENSION(:,:,:), INTENT(inout) :: amnv
      INTEGER, INTENT(in)                        :: parity

!  Start of executable code
      IF (parity .eq. f_cos) THEN
         CALL set_bndy_fouier_m0(amns, f_cos)
         CALL set_bndy_fouier_m0(amnu, f_sin)
         CALL set_bndy_fouier_m0(amnv, f_sin)
      ELSE
         CALL set_bndy_fouier_m0(amns, f_sin)
         CALL set_bndy_fouier_m0(amnu, f_cos)
         CALL set_bndy_fouier_m0(amnv, f_cos)
      END IF

      END SUBROUTINE

!-------------------------------------------------------------------------------
!>  @brief Compute the divergence of a half grid quantity.
!>
!>  This operator is used for unit testing of the curl operator. Divergence is
!>  defined as
!>
!>    DIV(B) = 1/sqrt(g)d/dx (JB^x}   for x = s,u,v                          (1)
!>
!>  Note for testing that divergence remains zero, add a 1/sqrt(g) term is not
!>  necessary to retain.
!>
!>  @param[in]    jasupsmnh Contravariant vector in the s direction half grid.
!>  @param[in]    jasupsmnh Contravariant vector in the u direction half grid.
!>  @param[in]    jasupsmnh Contravariant vector in the v direction half grid.
!>  @param[inout] divmnf    Diverence on the full grid.
!>  @param[in]    parity    Parity flag.
!>  @param[in]    nsmin     Minimum radial index.
!>  @param[in]    nsmax     Maximum radial index.
!-------------------------------------------------------------------------------
      SUBROUTINE div_h(jbsupsmnh, jbsupumnh, jbsupvmnh,                        &
                       divmnf, parity, nsmin, nsmax)
      USE descriptor_mod, ONLY: iam

      IMPLICIT NONE

!  Declare Arguments
      REAL (dp), DIMENSION(:,:,:), ALLOCATABLE, INTENT(in)    :: jbsupsmnh
      REAL (dp), DIMENSION(:,:,:), ALLOCATABLE, INTENT(in)    :: jbsupumnh
      REAL (dp), DIMENSION(:,:,:), ALLOCATABLE, INTENT(in)    :: jbsupvmnh
      REAL (dp), DIMENSION(:,:,:), ALLOCATABLE, INTENT(inout) :: divmnf
      INTEGER, INTENT(IN)                                     :: parity
      INTEGER, INTENT(IN)                                     :: nsmin
      INTEGER, INTENT(IN)                                     :: nsmax

!  Local Variables
      REAL (dp), ALLOCATABLE, DIMENSION(:,:,:)                :: jbsupumnf
      REAL (dp), ALLOCATABLE, DIMENSION(:,:,:)                :: jbsupvmnf
      INTEGER                                                 :: sparity
      INTEGER                                                 :: m
      INTEGER                                                 :: n
      INTEGER                                                 :: n_mode
      INTEGER                                                 :: s
      INTEGER                                                 :: fours
      INTEGER                                                 :: fouruv

!  Start of executable code.
      IF (parity .EQ. f_sin) THEN
         sparity = 1
         fours = f_sin
         fouruv = f_cos
      ELSE
         sparity = -1
         fours = f_cos
         fouruv = f_sin
      END IF

      CALL GradientFull(divmnf, jbsupsmnh)
      CALL set_bndy_half_to_full_ds(jbsupsmnh, divmnf, nsmin, fours, f_none)

      ALLOCATE(jbsupumnf(0:mpol,-ntor:ntor,nsmin:nsmax))
      ALLOCATE(jbsupvmnf(0:mpol,-ntor:ntor,nsmin:nsmax))

      CALL to_full_mesh(jbsupumnh, jbsupumnf)
      CALL to_full_mesh(jbsupvmnh, jbsupvmnf)

      CALL set_bndy_half_to_full(jbsupumnh, jbsupumnf, nsmin, fouruv, f_ds)
      CALL set_bndy_half_to_full(jbsupvmnh, jbsupvmnf, nsmin, fouruv, f_ds)

      DO s = nsmin, nsmax
         DO n = -ntor, ntor
            n_mode = fourier_context%tor_modes(n)*nfp
            DO m = 0, mpol
               divmnf(m,n,s) = divmnf(m,n,s)                                   &
                             - sparity*m*jbsupumnf(m,n,s)                      &
                             - sparity*n_mode*jbsupvmnf(m,n,s)
            END DO
         END DO
      END DO

      DEALLOCATE(jbsupumnf)
      DEALLOCATE(jbsupvmnf)

      END SUBROUTINE

!*******************************************************************************
!  UNIT TESTS
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  @brief Test utility routines.
!>
!>  Test the various operations of in the @ref utilities module. Current tests
!>  check:
!>
!>    * Check the div(curl(A)) is identially zero.
!-------------------------------------------------------------------------------
      FUNCTION test_utilities()
      USE descriptor_mod, ONLY: iam

      IMPLICIT NONE

!  Declare Arguments
      LOGICAL                                  :: test_utilities

!  Local variables
      REAL (dp), DIMENSION(:,:,:), POINTER     :: amnsf
      REAL (dp), DIMENSION(:,:,:), POINTER     :: amnuf
      REAL (dp), DIMENSION(:,:,:), POINTER     :: amnvf
      REAL (dp), DIMENSION(:,:,:), ALLOCATABLE :: jbmnsh
      REAL (dp), DIMENSION(:,:,:), ALLOCATABLE :: jbmnuh
      REAL (dp), DIMENSION(:,:,:), ALLOCATABLE :: jbmnvh
      INTEGER                                  :: s
      INTEGER                                  :: m
      INTEGER                                  :: n

!  Start of executable code.
      IF (iam .ne. 0) THEN
         test_utilities = .FALSE.
         RETURN
      END IF

      ALLOCATE(amnsf(0:mpol,-ntor:ntor,ns))
      ALLOCATE(amnuf(0:mpol,-ntor:ntor,ns))
      ALLOCATE(amnvf(0:mpol,-ntor:ntor,ns))

      CALL RANDOM_NUMBER(amnsf)
      CALL RANDOM_NUMBER(amnuf)
      CALL RANDOM_NUMBER(amnvf)

      ALLOCATE(jbmnsh(0:mpol,-ntor:ntor,ns))
      ALLOCATE(jbmnuh(0:mpol,-ntor:ntor,ns))
      ALLOCATE(jbmnvh(0:mpol,-ntor:ntor,ns))

      CALL set_bndy_full_origin(amnsf, amnuf, amnvf, f_none)
      CALL set_bndy_fouier_m0(amnsf, amnuf, amnvf, f_sin)

      CALL curl_ftoh(amnsf, amnuf, amnvf,                                      &
                     jbmnsh, jbmnuh, jbmnvh,                                   &
                     f_sin, 1, ns)

      test_utilities = test_div(jbmnsh, jbmnuh, jbmnvh, f_sin, 1, ns)

      DEALLOCATE(jbmnsh)
      DEALLOCATE(jbmnuh)
      DEALLOCATE(jbmnvh)

      DEALLOCATE(amnsf)
      DEALLOCATE(amnuf)
      DEALLOCATE(amnvf)

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Test divergence of a half grid quantity.
!>
!>  @param[in] jbsupsmnh Contravariant vector in the s direction half grid.
!>  @param[in] jbsupumnh Contravariant vector in the u direction half grid.
!>  @param[in] jbsupvmnh Contravariant vector in the v direction half grid.
!>  @param[in] parity    Parity flag.
!>  @param[in] nsmin     Minimum radial index.
!>  @param[in] nsmax     Maximum radial index.
!-------------------------------------------------------------------------------
      FUNCTION test_div(jbsupsmnh, jbsupumnh, jbsupvmnh, parity, nsmin, nsmax)
      USE descriptor_mod, ONLY: iam

      IMPLICIT NONE

!  Declare Arguments
      LOGICAL                                              :: test_div
      REAL (dp), DIMENSION(:,:,:), ALLOCATABLE, INTENT(in) :: jbsupsmnh
      REAL (dp), DIMENSION(:,:,:), ALLOCATABLE, INTENT(in) :: jbsupumnh
      REAL (dp), DIMENSION(:,:,:), ALLOCATABLE, INTENT(in) :: jbsupvmnh
      INTEGER, INTENT(IN)                                  :: parity
      INTEGER, INTENT(IN)                                  :: nsmin
      INTEGER, INTENT(IN)                                  :: nsmax

!  Local variables
      REAL (dp), DIMENSION(:,:,:), ALLOCATABLE             :: divbmnf
      INTEGER                                              :: s
      INTEGER                                              :: m
      INTEGER                                              :: n

!  Local parameters
      REAL (dp), PARAMETER :: tolarance = 3.0E-12

!  Start of executable code.
      ALLOCATE(divbmnf(0:mpol,-ntor:ntor,nsmin:nsmax))

      CALL div_h(jbsupsmnh, jbsupumnh, jbsupvmnh, divbmnf, f_sin, nsmin, nsmax)

      test_div = ANY(ABS(divbmnf(:,:,nsmin:nsmax-1)) .gt. tolarance)

      IF (test_div) THEN
         DO s = nsmin, nsmax-1
            DO n = -ntor, ntor
               DO m = 0, mpol
                  IF (ABS(divbmnf(m,n,s)) .gt. tolarance) THEN
                     WRITE (*,*) m, n, s, divbmnf(m,n,s)
                  END IF
               END DO
            END DO
         END DO
      END IF

      DEALLOCATE(divbmnf)

      END FUNCTION

      END MODULE
