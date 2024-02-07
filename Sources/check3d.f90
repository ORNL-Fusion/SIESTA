      SUBROUTINE check3d_symmetry(a,b,c,ap,bp,cp,mblk,mblkc,nblocks,    &
                                  asym_index)
      USE descriptor_mod, mm_loc=>mm, csrc_loc=>csrc
!      USE prof_mod
!      USE ptrd_mod
      IMPLICIT NONE

!Written 03-12-10 by Eduardo D'Azevedo

!fake        integer, parameter :: dp = kind(1.0d0)
!fake        integer, parameter :: dlen_ = 9, m_=2,n_=3
!fake        integer, dimension(DLEN_) :: descA, descA_1xp
!fake        integer, parameter :: Locq = 10, Locp = 10
 
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER, INTENT(in) :: nblocks, mblk, mblkc
      REAL(dp), TARGET, DIMENSION(mblk,mblkc,nblocks), INTENT(inout) :: &
                            a, b, c
      REAL(dp), TARGET, DIMENSION(Locq,Locp,nblocks), INTENT(inout) ::  &
                            ap
      REAL(dp), TARGET, DIMENSION(Locq,Locp,nblocks-1), INTENT(inout) ::&
                            bp, cp

      REAL(dp), INTENT(out) :: asym_index
#if defined(MPI_OPT)
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      INTEGER :: k
     
      REAL(dp) :: w1, w2

      REAL(dp), parameter :: tol = 1.0d-6

      REAL(dp), POINTER :: amat(:,:), bmat(:,:), cmat(:,:)
      REAL(dp), POINTER :: amatp(:,:), bmatp(:,:), cmatp(:,:)
      REAL(dp), DIMENSION(:,:), POINTER :: Asrc, Bsrc, Csrc
      REAL(dp), DIMENSION(:,:), ALLOCATABLE :: Atmp
      INTEGER, DIMENSION(DLEN_) :: descAtmp

!     ---------------------------
!     check a(:,:,:) - true - or ap(:,:,:) - false
!     ---------------------------
      LOGICAL, PARAMETER :: check_matp = .FALSE.

      INTEGER :: mm,nn
      REAL(dp) :: dnorm, anorm,bnorm,cnorm,  alpha,beta
      REAL(dp), EXTERNAL :: pdlange
      REAL(dp), DIMENSION(1) :: work
!-----------------------------------------------
      IF (.NOT.LSCALAPACK) RETURN

      w1 = 0
      w2 = 0
!     -----------------------------
!     setup storage for temp matrix
!     -----------------------------
      IF (check_matp) THEN
         ALLOCATE( Atmp(Locq,Locp) )
         descAtmp = descA
      ELSE
         ALLOCATE( Atmp(mblk,mblkc) )
         descAtmp = descA_1xp
      ENDIF

!     ---------------------
!     check diagonal blocks
!     ---------------------
      DO  k=1,nblocks
         amat => a(:,:,k)
         amatp => ap(:,:,k)
         IF (check_matp) THEN
           Asrc => amatp
         ELSE
           Asrc => amat
         ENDIF


!        ---------------------------
!        copy diagonal block to Atmp
!        ---------------------------
         Atmp = 0.0d0

         alpha = 1.0d0
         beta = 0.0d0
         mm = descAtmp(M_)
         nn = descAtmp(N_)
         CALL pdgeadd('Notranspose', mm,nn,                             &
                alpha,   Asrc,1,1,descAtmp,                             &
                beta,    Atmp,1,1,descAtmp )

         mm = descAtmp(M_)
         nn = descAtmp(N_)
         anorm = pdlange('F',mm,nn,Atmp,1,1,descAtmp,work)
!        ----------------------
!        compute  Atmp <- Atmp - transpose(Asrc)
!        ----------------------

         alpha = -1.0d0
         beta = 1.0d0
         mm = descAtmp(M_)
         nn = descAtmp(N_)
         CALL pdgeadd('Transpose', mm,nn,                               &
                alpha,   Asrc,1,1,descAtmp,                             &
                beta,    Atmp,1,1,descAtmp )


         mm = descAtmp(M_)
         nn = descAtmp(N_)
         dnorm = pdlange( 'F', mm,nn,  Atmp,1,1,descAtmp, work)
         w1 = w1 + dnorm*dnorm
         w2 = w2 + anorm*anorm

!        -------------------------
!        check off-diagonal blocks
!
!        Bsrc = transpose(Csrc)
!        -----------------------------
         IF (k < nblocks) THEN
           bmat => b(:,:,k+1)
           bmatp => bp(:,:,k)

           cmat => c(:,:,k)
           cmatp => cp(:,:,k)

         IF (check_matp) THEN
             Bsrc => bmatp
             Csrc => cmatp
         ELSE
             Bsrc => bmat
             Csrc => cmat
         ENDIF

         mm = descAtmp(M_)
         nn = descAtmp(N_)
         bnorm = pdlange('F',mm,nn,Bsrc,1,1,descAtmp,work)
         cnorm = pdlange('F',mm,nn,Csrc,1,1,descAtmp,work)
    
!        ---------------------------
!        copy off-diagonal B block to Atmp
!        ---------------------------
         Atmp = 0.0d0

         alpha = 1.0d0
         beta = 0.0d0
         mm = descAtmp(M_)
         nn = descAtmp(N_)
         CALL pdgeadd('Notranspose', mm,nn,                             &
                      alpha,Bsrc,1,1,descAtmp,                          &
                      beta,Atmp,1,1,descAtmp )

!        ----------------------
!        compute  Atmp <- Atmp - transpose(Csrc)
!        ----------------------

         alpha = -1.0d0
         beta = 1.0d0
         mm = descAtmp(M_)
         nn = descAtmp(N_)
         CALL pdgeadd('Transpose', mm,nn,                               &
                      alpha,Csrc,1,1,descAtmp,                          &
                      beta,Atmp,1,1,descAtmp )


         mm = descAtmp(M_)
         nn = descAtmp(N_)
         dnorm = pdlange( 'F', mm,nn,  Atmp,1,1,descAtmp, work)

         w1 = w1 + dnorm*dnorm
         w2 = w2 + bnorm*bnorm + cnorm*cnorm

        ENDIF

       ENDDO
       
       DEALLOCATE( Atmp )
       
       IF (w2 .ne. 0) THEN
          asym_index = SQRT(w1/w2)
       ELSE
          asym_index = -1
       END IF
#else
      asym_index = -1
#endif

      END SUBROUTINE check3d_symmetry


      SUBROUTINE check3d_symmetry_serial(dblk,lblk,ublk,mpol,ntor,      &
                                         mblk,nblocks,asym_index)
      USE stel_kinds
      USE shared_data, ONLY: ndims
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER, INTENT(in)  :: mblk, nblocks, mpol, ntor
      REAL(dp), INTENT(out)  :: asym_index
      REAL(dp), INTENT(in), DIMENSION(mblk,mblk,nblocks) ::             &
         dblk, lblk, ublk
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      INTEGER              :: ntype, m, n, icol
      REAL(dp)             :: w1, w2
!-----------------------------------------------
!
!     COMPUTES asym_index = SQRT(||A - A^T||**2/||A||**2)
!
      icol = 0

      w1 = 0;  w2 = 0
      DO ntype = 1, ndims
        DO n = -ntor, ntor
          DO m = 0, mpol
            icol = icol+1
            IF (m.EQ.0 .AND. n.LT.0) CYCLE
               w1 = w1 +                                                &
                  SUM((ublk(icol,:,1:nblocks-1)                         &
                    -  lblk(:,icol,2:nblocks))**2)
               w2 = w2 + SUM(ublk(icol,:,1:nblocks-1)**2                &
                       +     lblk(:,icol,2:nblocks)**2)
               w1 = w1 + SUM((dblk(icol,:,:) -                          &
                              dblk(:,icol,:))**2)
               w2 = w2 + SUM(dblk(icol,:,:)**2)
            END DO
         END DO
      END DO

      IF (w2 .NE. 0) THEN
         asym_index = SQRT(w1/w2)
      ELSE
         asym_index = -1
      END IF

      END SUBROUTINE check3d_symmetry_serial


      SUBROUTINE CheckForces(xc, gc)
      USE stel_kinds, ONLY: dp
      USE stel_constants, ONLY: zero, twopi
      USE descriptor_mod, ONLY: iam
      USE shared_data, ONLY: wtotal, l_linearize, l_init_state,         &
                             l_getfsq, l_getwmhd, l_ApplyPrecon,        &
                             l_push_edge, l_push_s, l_push_u, l_push_v, &
                             ndims
      USE shared_functions, ONLY: funct_island
      USE quantities, ONLY: fsubsmncf, fsubumnsf, fsubvmnsf,            &
                            mpol, ntor, ns, wp0, hs_i
      USE hessian, ONLY: l_Compute_Hessian
      USE nscalingtools, ONLY: startglobrow, endglobrow
      USE v3_utilities, ONLY: assert
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      REAL(dp) :: xc(ndims*(1+mpol)*(2*ntor+1),ns),                     &
                  gc(ndims*(1+mpol)*(2*ntor+1),ns)
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      REAL(dp), ALLOCATABLE :: gc1(:,:)
      INTEGER  :: ntype, n, m, icol, js, mblk, nt1, icount
      REAL(dp) ::  wplus, wmins, force, eps, fnorm
      LOGICAL  :: lsave, lsave2, lwrite, l_edge, l_s, l_u, l_v
      CHARACTER*(10) :: force_array(3) = (/' FORCE-S: ',' FORCE-U: ',' FORCE-V: '/)
!-----------------------------------------------
!
!     COMPUTES F = -GRAD(WMHD) NUMERICALLY
!
      mblk = SIZE(xc,1)
      icol = 0
      lsave = l_linearize
      lsave2= l_Compute_Hessian
      l_edge = l_push_edge
      l_s = l_push_s
      l_u = l_push_u
      l_v = l_push_v
      l_linearize = .FALSE.
      l_Compute_Hessian = .FALSE.
      l_init_state = .FALSE.
      l_ApplyPrecon = .FALSE.
      l_getfsq = .FALSE.
      l_getwmhd = .TRUE.
      l_push_edge = .TRUE.
      l_push_s = .TRUE.; l_push_u = .TRUE.; l_push_v = .TRUE.

      xc = 0
      CALL funct_island
      ALLOCATE (gc1(mblk,ns))
      
      fnorm = (twopi*twopi)*hs_i
      gc1 = gc*fnorm

      IF (iam .EQ. 0) PRINT *,'WRITING CHECK-FORCES FILES FORT.(3000+iam)'
      
      WRITE (3000+iam, *) 'fnorm: ', fnorm,' iam: ', iam

!      eps = 10*SQRT(epsilon(eps))
      eps = 1.E-3_dp

      CALL ASSERT(wtotal.GT.zero,'WTOTAL=0!')
      
      !Need l_push_s,u,v = T for js=1; l_push_edge = T for js=ns

      RADIAL: DO icount = 1, 4
      
         IF (icount .EQ. 1) js = 1
         IF (icount .EQ. 2) js = 2
         IF (icount .EQ. 3) js = ns/2
         IF (icount .EQ. 4) js = ns
      
         lWrite = (startglobrow.LE.js .AND. endglobrow.GE.js)
         IF (lWrite)                                                    &
         WRITE (3000+iam, *)'POINT: ', js,' START: ',                   &
                             startglobrow,' END: ', endglobrow

      icol = 0

      DO ntype = 1, ndims
         DO n = -ntor, ntor
            DO m = 0, mpol
               icol = icol+1
                
               xc(icol,js) = eps
               CALL funct_island
               wplus = wtotal

               xc(icol,js) = -eps
               CALL funct_island
               wmins = wtotal
               force = -(wplus-wmins)/(2*eps)
                
               CALL ASSERT(l_getwmhd,'L_GETWMHD = FALSE!')

               xc(icol,js) = 0

               IF (lWrite) THEN
               nt1 = ntype
               IF (ntype .GT. 3) nt1 = ntype-3
                   WRITE (3000+iam,100) m, n, ntype, force, force_array(nt1),  &
                          gc1(icol,js)
               END IF

            END DO
         END DO
      END DO

      IF (lWrite) WRITE (3000+iam,*)
      CALL FLUSH(3000+iam)
      
      END DO RADIAL
      

 100  FORMAT(1x,'M: ',i4,' N: ',i4,' NTYPE: ',i4,' GRAD-W: ',1pe14.4, a, 1pe14.4)

      CALL ASSERT(icol.EQ.mblk,'icol != mblk in CheckForces')
      l_linearize = lsave
      l_Compute_Hessian = lsave2
      l_push_edge = l_edge
      l_push_s = l_s; l_push_u = l_u; l_push_v = l_v
      xc = 0
      DEALLOCATE (gc1)
      
      END SUBROUTINE CheckForces


      SUBROUTINE TriDiag_Test (xc, gc, eps)
!
!     TESTS Tridiagonal STRUCTURE OF HESSIAN at js=ns/2
!
      USE stel_kinds, ONLY: dp
      USE stel_constants, ONLY: zero
      USE quantities, ONLY: fsubsmncf, fsubumnsf, fsubvmnsf, mpol, ntor, ns
      USE hessian, ONLY: l_Compute_Hessian
      USE shared_data, ONLY: ndims, l_linearize, l_init_state,          &
                                    l_ApplyPrecon, l_getwmhd
      USE shared_functions, ONLY: funct_island
      USE descriptor_mod, ONLY: iam
      USE v3_utilities, ONLY: assert
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      REAL(dp), DIMENSION(ndims*(1+mpol)*(2*ntor+1),ns) ::              &
	            xc, gc
	  REAL(dp), INTENT(IN) :: eps
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      INTEGER :: js, icol, ntype, n, m
      REAL(dp), ALLOCATABLE :: gc0(:,:)
      LOGICAL :: l_save, l_AppSave
!-----------------------------------------------

      js = ns/2
      l_save=l_linearize
      l_linearize=.FALSE.
      l_init_state = .TRUE.
      l_AppSave = l_ApplyPrecon
      l_ApplyPrecon = .FALSE.
      l_getwmhd = .TRUE.

      ALLOCATE (gc0(SIZE(gc,1), SIZE(gc,2)))
      xc = 0
      CALL funct_island
      gc0 = gc
      l_getwmhd = .FALSE.

      icol = 0
      DO ntype = 1, ndims
         DO n = -ntor, ntor
            DO m = 0, mpol
               icol = icol+1
               IF (m .eq. 0 .and. n.lt.0) CYCLE
               xc(icol,js) = eps
               CALL funct_island
               xc(icol,js) = 0
               gc = gc - gc0
               CALL ASSERT(ALL(gc(:,1:js-2) .EQ. zero),                 &
                          'TRI_DIAGONAL FAILED FOR LOWER BANDS')
               CALL ASSERT(ALL(gc(:,js+2:ns) .EQ. zero),                &
                          'TRI_DIAGONAL FAILED FOR UPPER BANDS')
            END DO
         END DO
      END DO

      IF (iam .EQ. 0) PRINT *,'TRI_DIAGONAL TEST PASSED'
      l_linearize=l_save
      l_ApplyPrecon = l_AppSave
      DEALLOCATE (gc0)

      END SUBROUTINE TriDiag_Test
