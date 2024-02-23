!*******************************************************************************
!>  @file shared_functions.f90
!>  @brief Contains module @ref shared_functions
!
!  Note separating the Doxygen comment block here so detailed decription is
!  found in the Module not the file.
!
!>  Subroutines and functions updating MHD forces and Wmhd.
!*******************************************************************************

      MODULE shared_functions
      USE v3_utilities, ONLY: assert
      USE shared_data
      USE utilities
      USE descriptor_mod, ONLY: INHESSIAN, nprocs, iam, SIESTA_COMM
      USE island_params, ns=>ns_i, hs=>hs_i, mpol=>mpol_i, ntor=>ntor_i
      USE hessian, ONLY: apply_precond, l_Compute_Hessian, apply_colscale
      USE timer_mod, ONLY: time_init_state, time_funci
      USE nscalingtools, ONLY: startglobrow, endglobrow, MPI_ERR,       &
                               PARSOLVER
      USE mpi_inc
      USE GMRES_LIB, ONLY: Truncate
      
      IMPLICIT NONE
      
      CONTAINS

!*******************************************************************************
!  UTILITY SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  @brief Compute forces for the perturbed state.
!-------------------------------------------------------------------------------
      SUBROUTINE funct_island
      USE siesta_bfield, ONLY: update_bfield
      USE siesta_pressure, ONLY: update_pres
      USE siesta_force, ONLY: update_force
      USE siesta_displacement, ONLY: update_upperv
      USE siesta_init, ONLY: init_state
      USE quantities, ONLY: wb, wp
      USE siesta_state, ONLY: Clear_Field_Perts

      IMPLICIT NONE

!  Local Variables
      REAL (dp) :: ton
      REAL (dp) :: toff
      REAL (dp) :: skston
      REAL (dp) :: skstoff
      LOGICAL   :: ltype
      INTEGER   :: nsmin
      INTEGER   :: nsmax

!  Start if executable code.
      CALL second0(ton)
      nsmin = MAX(1, startglobrow)
      nsmax = MIN(endglobrow, ns)

      IF (INHESSIAN) THEN
         in_hess_nfunct = in_hess_nfunct + 1
      ELSE
         out_hess_nfunct = out_hess_nfunct + 1
      END IF

      l_update_state = .false.
      IF (l_init_state) THEN
         CALL second0(skston)
         CALL init_state(.false.)
         CALL second0(skstoff)
         time_init_state = time_init_state + (skstoff - skston)
         l_init_state = .FALSE.
      END IF

!  There no need to compute the pertubations if the displacement vector did not
!  change.
      IF (ANY(xc .ne. 0) .or. ALLOCATED(buv_res)) THEN
         CALL update_upperv
         CALL update_bfield(.false.)
         CALL update_pres
      ELSE
         CALL Clear_Field_Perts
      END IF
      CALL update_force

      CALL ASSERT(gamma .ne. one, 'SIESTA REQUIRES gamma != 1')
      wtotal = wb + wp/(gamma - 1)
      IF (wtotal0 .eq. -1) THEN
         wtotal0 = wtotal
      END IF

      gc = gnorm_i*gc

!  Called from gmres. Add back the inital force gc0 or store unpreconditioned
!  forces.
      IF (ALLOCATED(gc0)             .and.                                     &
          (l_getfsq .or. l_conjgrad) .and.                                     &
          l_linearize                .and.                                     &
         .not.l_Compute_Hessian) THEN
         gc = gc + gc0
      END IF

      CALL ASSERT(.not.(l_getfsq .and. inhessian),                             &
                  'l_getfsq must be set to FALSE in Hessian')

      IF (l_getfsq) THEN
!  Compute preconditioned volume-averaged force.
         fsq_total = GetFSQ(nsmin, nsmax)

!  Compute un preconfitioned volume-averaged force.
         IF (ANY(col_scale .ne. one)) THEN
            CALL Apply_ColScale(gc, one/col_scale, nsmin, nsmax)

            fsq_total1 = GetFSQ(nsmin, nsmax)

            CALL Apply_ColScale(gc, col_scale, nsmin, nsmax)
         ELSE
            fsq_total1 = fsq_total
         END IF
      END IF

!  gmres handles preconditioning itself. Do not apply when called from inside
!  gmres or when the hessian is being computed.
      IF (l_ApplyPrecon) THEN
         CALL apply_precond(gc)
      END IF

      CALL second0(toff)
      time_funci = time_funci + (toff - ton)

      END SUBROUTINE

!-------------------------------------------------------------------------------
!>  @brief Computes the force residul.
!>
!>  The force residul is found by
!>
!>    |F^2| = Fsub.Fsup = f_s*f^s + f_u*f^u + f_v*f^v                        (3)
!>
!>  @param[in] nsmin Minimum radial index.
!>  @param[in] nsmax Maximum radial index.
!>  @returns |F^2|
!-------------------------------------------------------------------------------
      FUNCTION GetFSQ(nsmin, nsmax)
      USE quantities, ONLY: fsubsmncf, fsubumnsf, fsubvmnsf,            &
                            fsupsmncf, fsupumnsf, fsupvmnsf,            &
                            fsubsmnsf, fsubumncf, fsubvmncf,            &
                            fsupsmnsf, fsupumncf, fsupvmncf,            &
                            toupper_forces

      IMPLICIT NONE

!  Declare Arguments
      REAL (dp)               :: GetFSQ
      INTEGER, INTENT(in)     :: nsmin
      INTEGER, INTENT(in)     :: nsmax

!  Local Variables
      REAL (dp), DIMENSION(4) :: temp
      INTEGER                 :: js

!  Start of executable code.
      CALL toupper_forces

      temp(1) = SUM(fsubsmncf(:,:,nsmin:nsmax)**2)
      temp(2) = SUM(fsubumnsf(:,:,nsmin:nsmax)**2)
      temp(3) = SUM(fsubvmnsf(:,:,nsmin:nsmax)**2)
      IF (lasym) THEN
         temp(1) = temp(1) + SUM(fsubsmnsf(:,:,nsmin:nsmax)**2)
         temp(2) = temp(2) + SUM(fsubumncf(:,:,nsmin:nsmax)**2)
         temp(3) = temp(3) + SUM(fsubvmncf(:,:,nsmin:nsmax)**2)
      END IF

!  Volume avereages |F|^2.
      temp(4) = 0
      DO js = nsmin, nsmax
         temp(4) = temp(4)                                                     &
                 + vp_f(js)*SUM(fsupsmncf(:,:,js)*fsubsmncf(:,:,js) +          &
                                fsupumnsf(:,:,js)*fsubumnsf(:,:,js) +          &
                                fsupvmnsf(:,:,js)*fsubvmnsf(:,:,js))
         IF (lasym) THEN
            temp(4) = temp(4)                                                  &
                    + vp_f(js)*SUM(fsupsmnsf(:,:,js)*fsubsmnsf(:,:,js) +       &
                                   fsupumncf(:,:,js)*fsubumncf(:,:,js) +       &
                                   fsupvmncf(:,:,js)*fsubvmncf(:,:,js))
         END IF
      END DO
#if defined(MPI_OPT)         
      IF (PARSOLVER) THEN
!  FIXME: All reduce is not deterministic. This causes a divergent run sequence.
         CALL MPI_ALLREDUCE(MPI_IN_PLACE, temp, 4, MPI_REAL8, MPI_SUM,         &
                            SIESTA_COMM, MPI_ERR)
      END IF
#endif

      fsqvs = hs*temp(1)
      fsqvu = hs*temp(2)
      fsqvv = hs*temp(3)

      GetFSQ = temp(4)/SUM(vp_f)

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Get the perturbed MHD and Kinetic energy.
!>
!>  @param[in] p Displacement vector.
!>  @return The total stored energy.
!-------------------------------------------------------------------------------
      FUNCTION getwmhd(p)

      IMPLICIT NONE

!  Declare Arguments
      REAL (dp)                           :: getwmhd
      REAL (dp), DIMENSION(:), INTENT(in) :: p

!  Start of executable code.
      xc = p
      l_init_state = .true.
      l_linearize  = .false.
      l_getwmhd    = .true.
      l_getfsq     = .false.

      CALL funct_island

      l_getwmhd    = .false.

      getwmhd = wtotal

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Perform line search along xc vector for minimum force residul.
!>
!>  @param[inout] xcmin   Scaled displacement vector.
!>  @param[out]   fsq_min Resdule at the minimum.
!-------------------------------------------------------------------------------
      SUBROUTINE LineSearch(xcmin, fsq_min)
      USE siesta_namelist, ONLY: ftol

      IMPLICIT NONE

!  Declare Arguments
      REAL(dp), DIMENSION(neqs), INTENT(inout) :: xcmin
      REAL(dp), INTENT(out)                    :: fsq_min

!  Local Variables
      REAL (dp)                                :: facmin
      REAL (dp)                                :: fsq_start
      INTEGER                                  :: iter
!      INTEGER                                  :: j

!  Start of executable code.
      facmin = 1 !3 !mrc
      l_init_state = .false.
      l_getfsq = .true.
      fsq_start = fsq_total1

      IF (iam .eq. 0 .and. lverbose) THEN
         WRITE (*,1000)
      END IF

      xc = facmin*xcmin
      iter = 0
      fsq_min = HUGE(fsq_min)

      CALL funct_island
      DO WHILE(fsq_total1 .lt. fsq_min)
         IF (fsq_total1 .lt. fsq_start) THEN
            xcmin = xc
         END IF
         fsq_min = fsq_total1
         iter = iter + 1
         IF (iam .eq. 0 .and. lverbose) THEN
            WRITE (*,1001) iter, fsq_total1, SQRT(SUM(xc*xc)),                 &
                           MAXVAL(ABS(xc)), facmin
            facmin = facmin/SQRT2
         END IF
         xc = xc/SQRT2
         CALL funct_island
      END DO

      fsq_min = MIN(fsq_min, fsq_start)

#if 0
      DO j = 1, 100
         CALL funct_island
         IF (fsq_total1 .lt. fsq_min) THEN
            xcmin = xc
            IF (fsq_total1 .gt. 0.98_dp*fsq_min) THEN
               iter = iter + 1
            END IF
            IF (fsq_total1 .lt. 0.85_dp*fsq_min) THEN
               iter = 0
            END IF
            fsq_min = fsq_total1
         ELSE IF (j .gt. 4) THEN
            iter = iter + 1
         END IF
         IF (iam .eq. 0 .and. lverbose) THEN
            WRITE (*,1001) j, fsq_total1, SQRT(SUM(xc*xc)), MAXVAL(ABS(xc)),   &
                           facmin
            facmin = facmin/SQRT2
         END IF
         IF (iter .gt. 2 .or. fsq_total1 .le. ftol) THEN
            EXIT
         END IF
         xc = xc/SQRT2
      END DO
#endif

1000  FORMAT(/,1x,'LINE SEARCH - SCAN ||X|| FOR MIN FSQ_NL',/,1x,15('-'),      &
             /,1x,'ITER',7x,'FSQ_NL',10x,'||X||',9x,'MAX|X|',10x,'FAC')
1001  FORMAT(i5,4(3x,1pe12.3))
 
      END SUBROUTINE

!-------------------------------------------------------------------------------
!>  @brief Initializes pointers for xc (displacement variable) and gc (force)
!>         residules.
!>
!>  @param[in]  xtarget   State array to point to.
!>  @param[out] s_ptr_sym Stellarator symmetric term for the s component.
!>  @param[out] u_ptr_sym Stellarator symmetric term for the u component.
!>  @param[out] v_ptr_sym Stellarator symmetric term for the v component.
!>  @param[out] s_ptr_asym Stellarator asymmetric term for the s component.
!>  @param[out] u_ptr_asym Stellarator asymmetric term for the u component.
!>  @param[out] v_ptr_asym Stellarator asymmetric term for the v component.
!-------------------------------------------------------------------------------
      SUBROUTINE init_ptrs(xtarget, s_ptr_sym, u_ptr_sym, v_ptr_sym,           &
                           s_ptr_asym, u_ptr_asym, v_ptr_asym)

      IMPLICIT NONE

!  Declare Arguments
      REAL (dp), DIMENSION(0:mpol,-ntor:ntor,ndims,ns), TARGET, INTENT(in) ::  &
         xtarget
      REAL (dp), DIMENSION(:,:,:), POINTER :: s_ptr_sym
      REAL (dp), DIMENSION(:,:,:), POINTER :: u_ptr_sym
      REAL (dp), DIMENSION(:,:,:), POINTER :: v_ptr_sym
      REAL (dp), DIMENSION(:,:,:), POINTER :: s_ptr_asym
      REAL (dp), DIMENSION(:,:,:), POINTER :: u_ptr_asym
      REAL (dp), DIMENSION(:,:,:), POINTER :: v_ptr_asym

!  Start of executable code.

!  Set stellarator symmetric pointers.
      s_ptr_sym => xtarget(:,:,1,:)
      u_ptr_sym => xtarget(:,:,2,:)
      v_ptr_sym => xtarget(:,:,3,:)

!  The number of dimensions if double when lasym is true. Set the stellarator
!  asymmetric pointers.
      IF (ndims .eq. 6) THEN
         s_ptr_asym => xtarget(:,:,4,:)
         u_ptr_asym => xtarget(:,:,5,:)
         v_ptr_asym => xtarget(:,:,6,:)
      END IF

      END SUBROUTINE

      END MODULE
