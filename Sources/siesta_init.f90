!*******************************************************************************
!>  @file siesta_init.f90
!>  @brief Contains module @ref siesta_init
!
!  Note separating the Doxygen comment block here so detailed decription is
!  found in the Module not the file.
!
!>  Initializes unperturbed siesta fields and pressure in real space.
!*******************************************************************************
      MODULE siesta_init
      USE v3_utilities, ONLY: assert, assert_eq
      USE stel_kinds
      USE descriptor_mod, ONLY: SIESTA_COMM, iam
      USE fourier, ONLY: f_none, f_sin, f_cos, f_du, f_dv, f_sum, m0, m1, m2
      USE island_params, ONLY: fourier_context
      USE nscalingtools, ONLY: startglobrow, endglobrow, MPI_ERR
      USE quantities
      USE shared_data, ONLY: l_update_state, l_pedge
      USE mpi_inc
 
      IMPLICIT NONE

      PRIVATE
!*******************************************************************************
!  utilities module variables
!*******************************************************************************
! FIXME: Remove these if possible.
      INTEGER                                  :: nsmin
      INTEGER                                  :: nsmax
      REAL (dp), DIMENSION(:,:,:), ALLOCATABLE :: pfilter
      REAL (dp)                                :: pedge
       
      PUBLIC  :: init_state

      CONTAINS

!-------------------------------------------------------------------------------
!>  @brief Initialize equilibrium state.
!>
!>  @param[in] lcurrent_only Only update to the current.
!>  @param[in] lpar_in       Controls if this is a parallel call or a serial
!>                           call.
!-------------------------------------------------------------------------------
      SUBROUTINE init_state(lcurrent_only, lpar_in)
      USE bhtobf
      USE utilities, ONLY: GradientFull
      USE metrics, ONLY: tolowerf
      USE shared_data, ONLY: fsq_total, l_getfsq, ste, buv_res,                &
                             l_update_state
      USE hessian, ONLY: muPar
      USE siesta_currents, ONLY: cv_currents
      USE island_params, ONLY: ohs=>ohs_i

      IMPLICIT NONE

!  Declare Arguments
      LOGICAL, INTENT(IN)           :: lcurrent_only
      LOGICAL, INTENT(IN), OPTIONAL :: lpar_in

!  local variables
      INTEGER                       :: istat
      INTEGER                       :: n1
      INTEGER                       :: n2
      INTEGER                       :: nmax
      INTEGER                       :: endr
      INTEGER                       :: startr
      REAL (dp)                     :: ton
      REAL (dp)                     :: toff
      REAL (dp)                     :: temp(2)
      LOGICAL                       :: ladd_pert0
      LOGICAL                       :: lpar

!  Start of executable code
!  jbupXmnPh and jpmnPh are initially computed on full ns mesh, and are gathered
!  to that full ns mesh in update_state
      istat = 0
      IF (PRESENT(lpar_in)) THEN
         lpar = lpar_in
      ELSE
         lpar = .true.
      END IF

      IF (.not.lpar) THEN
         startr = startglobrow
         endr = endglobrow
         startglobrow = 1
         endglobrow = ns
      END IF
      nsmin = MAX(1, startglobrow - 1)
      nsmax = MIN(ns, endglobrow + 2)
       
!Unnecessary to recompute magnetic fields or currents: already stored at this point!
      ladd_pert0 = ALLOCATED(buv_res)
      IF (.not.ladd_pert0) THEN
         CALL Init_Allocate_Arrays(lpar)

!  Calculate unperturbed half mesh jac*bsupX and jac*pressure from nsmin to
!  nsmax. Sum symmetric and (for lasym) asymmetric parts. Returning from
!  update_state, jbsupXmn, jpmn have been gathered onto all processors.
         CALL GetFields(jbsupsmnsh, jbsupumnch, jbsupvmnch, jpmnch, f_sin)
         IF (lasym) THEN
            CALL GetFields(jbsupsmnch, jbsupumnsh, jbsupvmnsh, jpmnsh, f_cos)
         END IF

!  Convert jac*bsubXh to bsupXh and jac*pres to presh in real space. On exit,
!  bsupXf and pf0 ae valid on [nsmin, nsmax - 1]
         CALL bhalftobfull(bsupsijh0, bsupuijh0, bsupvijh0,                    &
                           bsupsijf0, bsupuijf0, bsupvijf0,                    &
                           pijh0, pijf0)

!  Compute and store unperturbed current components on full mesh,
!  KsupX = sqrt(g)*JsupX
         CALL cv_currents(bsupsijh0, bsupuijh0, bsupvijh0,                     &
                          ksupsijf0, ksupuijf0, ksupvijf0,                     &
                          l_getfsq, .TRUE.)

         IF (lcurrent_only) THEN
            RETURN
         END IF
      END IF

!  Need this for resistive diffusion (|| resistivity) and mupar damping
      IF (ladd_pert0 .or. muPar .ne. zero) THEN
         nmax = MIN(endglobrow + 1, ns)
!  Calculate K || B
         CALL tolowerf(bsupsijf0,bsupuijf0,bsupvijf0,                          &
                       bsubsijf, bsubuijf, bsubvijf,                           &
                       nsmin, nmax)
!  IN UPDATE-BFIELD, KSUB = jacob*JSUB: PARALLEL CURRENT
         IF (ladd_pert0) THEN
            ksubsijf(:,:,nsmin:nmax) = bsubsijf(:,:,nsmin:nmax)
            ksubuijf(:,:,nsmin:nmax) = bsubuijf(:,:,nsmin:nmax)
            ksubvijf(:,:,nsmin:nmax) = bsubvijf(:,:,nsmin:nmax)
            RETURN
         END IF
      END IF
     
!  Compute spectrally filtered dp/du and dp/dv on half radial grid [nsmin:nsmax]
      IF (l_update_state) THEN
         ALLOCATE(pfilter(ntheta, nzeta, nsmin:nsmax))
      END IF
      CALL FilterPressure(f_cos)
      IF (lasym) THEN
         CALL FilterPressure(f_sin)
      END IF

!  Compute radial pressure derivative pijf0_ds at full mesh points
!  [nsmin:nsmax-1]
      CALL GradientFull(pijf0_ds, pijh0)

!  SPH: one-sided (m=1) derivative at origin yields factor of 2 (1/(hs/2)).
!  pijh(:,:,1) should contain the value of jpmnch(m0,:,2) from Init_Fields.
      IF (nsmin .EQ. 1) THEN
         pijf0_ds(:,:,1) = 2.0*(pijh0(:,:,2) - pijh0(:,:,1))*ohs
      END IF

!  SPH11-3-16 preserve s=1 as iso-pressure contour
      IF (nsmax .eq. ns .and. l_pedge .and. .not.l_vessel) THEN
         pijf0(:,:,ns) = pedge
         pijf0_ds(:,:,ns) = 0
      END IF

!   Fourier filter (Nyquist) check
      IF (.not. l_update_state) THEN
         RETURN
      END IF

!  pijh is unfiltered
      n1 = MAX(1, startglobrow)
      n2 = MIN(ns, endglobrow)
      IF (n1 .EQ. 1) THEN
         pfilter(:,:,1) = 0
      END IF
      temp(1) = SUM((pijh0(:,:,n1:n2) - pfilter(:,:,n1:n2))**2)
      temp(2) = SUM((pijh0(:,:,n1:n2) + pfilter(:,:,n1:n2))**2)
#if defined(MPI_OPT)
      CALL MPI_ALLREDUCE(MPI_IN_PLACE, temp, 2, MPI_REAL8, MPI_SUM,            &
                         SIESTA_COMM, MPI_ERR)
#endif
      IF (temp(2) .ne. zero) THEN
         ste(1) = SQRT(temp(1)/temp(2))
      END IF
      DEALLOCATE(pfilter)

      IF (.not.lpar) THEN
         startglobrow = startr
         endglobrow = endr
      END IF

      END SUBROUTINE

!>  \brief subroutine for allocating unperturbed array
       SUBROUTINE Init_Allocate_Arrays(lpar)
       USE shared_data, ONLY: l_par_state
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
       LOGICAL, INTENT(in)     :: lpar
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
       LOGICAL                 :: lrealloc, lalloc
       INTEGER                 :: istat, n1, n2, n3, nsmin, nsmax,      &
                                  nsmin1, nsmax1
!-----------------------------------------------
       lalloc = ALLOCATED(jvsupsijf) 
       lrealloc = .TRUE.
       istat = 0
       nsmin =MAX(1,startglobrow);   nsmax =MIN(endglobrow+1,ns)
       nsmin1=MAX(1,startglobrow-1); nsmax1=MIN(endglobrow+2,ns)
       n1 = ntheta; n2 = nzeta
       
       IF (lpar) THEN
          n3 = nsmax1-nsmin1+1
          l_par_state=.TRUE.
       ELSE
          n3 = ns
          l_par_state=.FALSE.
       END IF

       IF (lalloc) THEN
          IF (SIZE(bsupsijh0,3) .EQ. n3) lrealloc=.FALSE.
          IF (lrealloc) THEN
          DEALLOCATE (jvsupsijf, jvsupuijf, jvsupvijf,                  &
                      bsupsijf0, bsupuijf0, bsupvijf0,                  &
                      bsupsijf,  bsupuijf,  bsupvijf,                     &
                      bsupsijh0, bsupuijh0, bsupvijh0,                  &
                      bsubsijf,  bsubuijf,  bsubvijf,                   &
                      ksupsijf0, ksupuijf0, ksupvijf0,                  &
                      ksupsijf,  ksupuijf,  ksupvijf,                     &
                      ksubsijf,  ksubuijf,  ksubvijf,                   &
                      pijf0, pijh0, pijf0_ds,                           &
                      pijh0_du, pijh0_dv,  stat=istat)
          CALL assert_eq(0, istat,                                             &
                         'Deallocation error #1 in INIT_ALLOCATE_ARRAYS')
          END IF
       END IF

       IF (lrealloc) THEN
          n3 = ns
          ALLOCATE (jvsupsijf(n1,n2,n3), jvsupuijf(n1,n2,n3),           & ! Full mesh quantities (real space)
                    jvsupvijf(n1,n2,n3))                                  ! jacobian*(V_s , V_u , V_v)
          IF (lpar) THEN
          ALLOCATE (bsupsijh0(n1,n2,nsmin1:nsmax1),                     & ! B^s, B^u, B^v (half)
                    bsupuijh0(n1,n2,nsmin1:nsmax1),                     &
                    bsupvijh0(n1,n2,nsmin1:nsmax1),                     &
                    bsupsijf0(n1,n2,nsmin1:nsmax1),                     & ! B^s, B^u, B^v (full)
                    bsupuijf0(n1,n2,nsmin1:nsmax1),                     &
                    bsupvijf0(n1,n2,nsmin1:nsmax1),                     &
                    ksupsijf0(n1,n2,nsmin1:nsmax),                      & ! K^s, K^u, K^v 
                    ksupuijf0(n1,n2,nsmin1:nsmax),                      &
                    ksupvijf0(n1,n2,nsmin1:nsmax),                      &
                    ksubsijf(n1,n2,nsmin1:nsmax),                       & ! Full mesh quantities (real space)
                    ksubuijf(n1,n2,nsmin1:nsmax),                       & ! K_s, K_u, K_v
                    ksubvijf(n1,n2,nsmin1:nsmax),                       &
                    bsubsijf(n1,n2,nsmin1:nsmax),                       & ! Full mesh quantities (real space)
                    bsubuijf(n1,n2,nsmin1:nsmax),                       & ! B_s, B_u, B_v
                    bsubvijf(n1,n2,nsmin1:nsmax),                       &
                    ksupsijf(n1,n2,nsmin:nsmax),                           & 
                    ksupuijf(n1,n2,nsmin:nsmax),                           & 
                    ksupvijf(n1,n2,nsmin:nsmax),                           & 
                    bsupsijf(n1,n2,nsmin:nsmax),                           &
                    bsupuijf(n1,n2,nsmin:nsmax),                           &
                    bsupvijf(n1,n2,nsmin:nsmax),                           &
                    pijh0(n1,n2,nsmin1:nsmax1),                         &
                    pijh0_du(n1,n2,nsmin1:nsmax),                       &
                    pijh0_dv(n1,n2,nsmin1:nsmax),                       &
                    pijf0(n1,n2,nsmin1:nsmax1),                         &
                    pijf0_ds(n1,n2,nsmin1:nsmax), stat=istat)
          ELSE
          ALLOCATE (ksupsijf0(n1,n2,n3), ksupuijf0(n1,n2,n3),           &
                    ksupvijf0(n1,n2,n3), ksupsijf(n1,n2,n3),               &
                    ksupuijf(n1,n2,n3),  ksupvijf(n1,n2,n3),               &
                    ksubsijf(n1,n2,n3),                                 &
                    ksubuijf(n1,n2,n3),                                 &
                    ksubvijf(n1,n2,n3),                                 &
                    bsubsijf(n1,n2,n3),                                 &
                    bsubuijf(n1,n2,n3),                                 &
                    bsubvijf(n1,n2,n3),                                 &
                    pijf0(n1,n2,n3), pijh0(n1,n2,n3),                   &
                    pijh0_du(n1,n2,n3), pijh0_dv(n1,n2,n3),             &
                    pijf0_ds(n1,n2,n3), bsupsijf0(n1,n2,n3),            &
                    bsupuijf0(n1,n2,n3),bsupvijf0(n1,n2,n3),            &
                    bsupsijf(n1,n2,n3), bsupuijf(n1,n2,n3),               &
                    bsupvijf(n1,n2,n3), bsupsijh0(n1,n2,n3),              &
                    bsupuijh0(n1,n2,n3),bsupvijh0(n1,n2,n3),            &
                    stat=istat)
          END IF
          CALL assert_eq(0, istat,                                             &
                         'Allocation error #1 in Init_Allocate_Arrays')
          jvsupsijf = 0; jvsupuijf = 0; jvsupvijf = 0                     !Need in add_resistivity loop
       END IF


       lalloc = ALLOCATED(ksupsmnsf) 
       lrealloc = .TRUE.

       IF (lpar) THEN
          n1 = mpol+1; n2 = 2*ntor+1; n3 = nsmax-nsmin1+1
       ELSE
          n1 = mpol+1; n2 = 2*ntor+1; n3 = ns
       END IF

       IF (lalloc) THEN
          IF (SIZE(ksupsmnsf,3) .EQ. n3) THEN
             lrealloc = .FALSE.
          END IF
          IF (lrealloc) THEN
             DEALLOCATE(ksupsmnsf, ksupumncf, ksupvmncf,                    &
                        djpmnch, djbsupsmnsh, djbsupumnch, djbsupvmnch,     &
                        stat=istat)
             CALL assert_eq(0, istat,                                          &
                            'Deallocation error #2 in Init_Allocate_Arrays')
             IF (lasym) THEN
                DEALLOCATE(ksupsmncf, ksupumnsf, ksupvmnsf,                 &
                           djpmnsh, djbsupsmnch, djbsupumnsh, djbsupvmnsh,  &
                           stat=istat)
                CALL assert_eq(0, istat,                                       &
                               'Deallocation error #3 in Init_Allocate_Arrays')
             END IF
          END IF
       END IF

       IF (lrealloc) THEN
          IF (lpar) THEN
                ALLOCATE(ksupsmnsf(0:mpol,-ntor:ntor,nsmin1:nsmax),            &
                         ksupumncf(0:mpol,-ntor:ntor,nsmin1:nsmax),            &
                         ksupvmncf(0:mpol,-ntor:ntor,nsmin1:nsmax),            &
                         djpmnch(0:mpol,-ntor:ntor,nsmin:nsmax),               &
                         djbsupsmnsh(0:mpol,-ntor:ntor,nsmin:nsmax),           &
                         djbsupumnch(0:mpol,-ntor:ntor,nsmin:nsmax),           &
                         djbsupvmnch(0:mpol,-ntor:ntor,nsmin:nsmax),           &
                         stat=istat)
             IF (lasym) THEN
                ALLOCATE(ksupsmncf(0:mpol,-ntor:ntor,nsmin1:nsmax),            &
                         ksupumnsf(0:mpol,-ntor:ntor,nsmin1:nsmax),            &
                         ksupvmnsf(0:mpol,-ntor:ntor,nsmin1:nsmax),            &
                         djpmnsh(0:mpol,-ntor:ntor,nsmin:nsmax),               &
                         djbsupsmnch(0:mpol,-ntor:ntor,nsmin:nsmax),           &
                         djbsupumnsh(0:mpol,-ntor:ntor,nsmin:nsmax),           &
                         djbsupvmnsh(0:mpol,-ntor:ntor,nsmin:nsmax),           &
                         stat=istat)
             END IF
          ELSE
             ALLOCATE(ksupsmnsf(0:mpol,-ntor:ntor,n3),                         &
                      ksupumncf(0:mpol,-ntor:ntor,n3),                         &
                      ksupvmncf(0:mpol,-ntor:ntor,n3),                         &
                      djpmnch(0:mpol,-ntor:ntor,n3),                           &
                      djbsupsmnsh(0:mpol,-ntor:ntor,n3),                       &
                      djbsupumnch(0:mpol,-ntor:ntor,n3),                       &
                      djbsupvmnch(0:mpol,-ntor:ntor,n3),                       &
                      stat=istat)
             IF (lasym) THEN
                ALLOCATE(ksupsmncf(0:mpol,-ntor:ntor,n3),                      &
                         ksupumnsf(0:mpol,-ntor:ntor,n3),                      &
                         ksupvmnsf(0:mpol,-ntor:ntor,n3),                      &
                         djpmnsh(0:mpol,-ntor:ntor,n3),                        &
                         djbsupsmnch(0:mpol,-ntor:ntor,n3),                    &
                         djbsupumnsh(0:mpol,-ntor:ntor,n3),                    &
                         djbsupvmnsh(0:mpol,-ntor:ntor,n3),                    &
                         stat=istat)
             END IF
          END IF
          CALL assert_eq(0, istat,                                             &
                         'Allocation error #2 in Init_Allocate_Arrays')
          djpmnch = 0
          djbsupsmnsh = 0
          djbsupumnch = 0
          djbsupvmnch = 0
          IF (lasym) THEN
             djpmnsh = 0
             djbsupsmnch = 0
             djbsupumnsh = 0
             djbsupvmnsh = 0
          END IF
       END IF

       END SUBROUTINE
       
!-------------------------------------------------------------------------------
!>  @brief Get magnetic field and pressure on the real space mesh.
!>
!>  All quantities contain a Jacobian term.
!>
!>  @param[in] jbsupsmnh Contravariant B field s component on the half grid.
!>  @param[in] jbsupumnh Contravariant B field u component on the half grid.
!>  @param[in] jbsupvmnh Contravariant B field v component on the half grid.
!>  @param[in] jpmnh     Pressure on the half grid.
!>  @param[in]    iparity   Parity of the quantities.
!-------------------------------------------------------------------------------
      SUBROUTINE GetFields(jbsupsmnh, jbsupumnh, jbsupvmnh, jpmnh, iparity)

      IMPLICIT NONE

!  Declare Arguments
      REAL (dp), DIMENSION(:,:,:), ALLOCATABLE, INTENT(in) :: jbsupsmnh
      REAL (dp), DIMENSION(:,:,:), ALLOCATABLE, INTENT(in) :: jbsupumnh
      REAL (dp), DIMENSION(:,:,:), ALLOCATABLE, INTENT(in) :: jbsupvmnh
      REAL (dp), DIMENSION(:,:,:), ALLOCATABLE, INTENT(in) :: jpmnh
      INTEGER, INTENT(in)                                  :: iparity

!  local variables
      INTEGER                                              :: fours
      INTEGER                                              :: fouruvp
      INTEGER                                              :: fcomb

!  Start of executable code
!  Check origin values. Use these boundary conditions through out the code.
!      IF (nsmin .EQ. 1) THEN
!         CALL assert(ALL(jbsupsmnh(m1,:,1)    .eq. jbsupsmnh(m1,:,2)),         &
!                     'bmnsh(1) != bmnsh(2)')
!         CALL assert(ALL(jbsupsmnh(m0,:,1)    .eq. zero), 'bmnsh(m0:,1) != 0')
!         CALL assert(ALL(jbsupsmnh(m2:,:,1)   .eq. zero), 'bmnsh(m2:,1) != 0')
!         CALL assert(ALL(jbsupumnh(m0:m1,:,1) .eq. jbsupumnh(m0:m1,:,2)),      &
!                     'bmnuh(1) != bmnuh(2)')
!         CALL assert(ALL(jbsupumnh(m2:,:,1)   .eq. zero), 'bmnuh(m2:,1) != 0')
!         CALL assert(ALL(jbsupvmnh(m0,:,1)    .eq. jbsupvmnh(m0,:,2)),         &
!                     'bmnvh(m0,1) != bmnvh(m0,2)')
!         CALL assert(ALL(jbsupvmnh(m1:,:,1)   .eq. zero), 'bmnvh(m1:,1) != 0')
!         CALL assert(ALL(jpmnh(m0,:,1)        .eq. jpmnh(m0,:,2)),             &
!                     'pmnh(m0,1) != pmnh(m0,2)')
!         CALL assert(ALL(jpmnh(m1:,:,1)       .eq. zero), 'pmnh(m1:,1) != 0')
!      END IF

      IF (iparity .EQ. f_sin) THEN
         fcomb = f_none
         fours = f_sin
         fouruvp = f_cos
      ELSE
         fcomb = f_sum
         fours = f_cos
         fouruvp = f_sin
      END IF

      CALL fourier_context%toijsp(jbsupsmnh(:,:,nsmin:nsmax), bsupsijh0,       &
                                  fcomb, fours)
      CALL fourier_context%toijsp(jbsupumnh(:,:,nsmin:nsmax), bsupuijh0,       &
                                  fcomb, fouruvp)
      CALL fourier_context%toijsp(jbsupvmnh(:,:,nsmin:nsmax), bsupvijh0,       &
                                  fcomb, fouruvp)
      CALL fourier_context%toijsp(jpmnh(:,:,nsmin:nsmax), pijh0, fcomb, fouruvp)

      END SUBROUTINE

!-------------------------------------------------------------------------------
!>  @brief Compute angular derivatives of the pressure.
!>
!>  Pressure is filtered by Fourier tranforming to mn space then back to ij
!>  space. Angles derivatives with respect to u and v are also compute here.
!>
!>  @param[in] parity Fourier parity of the pressure.
!-------------------------------------------------------------------------------
      SUBROUTINE FilterPressure(parity)
      USE island_params, ONLY: fourier_context

      IMPLICIT NONE

!  Declare Arguments
      INTEGER, INTENT(in)                      :: parity

!  local variables
      INTEGER                                  :: fcomb
      INTEGER                                  :: istat
      REAL (dp), DIMENSION(:,:,:), ALLOCATABLE :: pmnh

!  Start of executable code
      IF (parity .eq. f_cos) THEN
         fcomb = f_none
      ELSE
         fcomb = f_sum
      END IF

!  pmnh contains Nyquist-filtered (to MPOL, NTOR) Fourier harmonics of p on
!  half mesh
      ALLOCATE(pmnh(0:mpol,-ntor:ntor,nsmin:nsmax), stat=istat)
      CALL assert_eq(0, istat, 'Allocation error in FilterPressure')

      CALL fourier_context%tomnsp(pijh0, pmnh, parity)
       
      IF (parity .eq. f_cos .and. nsmax .eq. ns) THEN
         pedge = pmnh(m0,0,ns)*fourier_context%orthonorm(0,0)
      END IF

!  Spectrally filtered angle derivatives of p in real space (half mesh).
      CALL fourier_context%toijsp(pmnh, pijh0_du, IOR(f_du, fcomb), parity)
      CALL fourier_context%toijsp(pmnh, pijh0_dv, IOR(f_dv, fcomb), parity)

      IF (l_update_state) THEN
         CALL fourier_context%toijsp(pmnh, pfilter, fcomb, parity)
      END IF

      DEALLOCATE(pmnh)

      END SUBROUTINE

!-------------------------------------------------------------------------------
!>  @brief Check that p > 0.
!-------------------------------------------------------------------------------
      SUBROUTINE CheckPressure

      IMPLICIT NONE

!  Start of executable code
      IF (lverbose .and. MINVAL(pijh0) .lt. zero) THEN
         WRITE (*,1000)
      END IF

!  This may stabilize convergence (lower f2 exponent or raise lev_marq)
      WHERE (pijh0 .LT. zero)
         pijh0 = zero
      END WHERE

1000  FORMAT(' pijh < 0 in init_state. ===> Recommend lowering f2 exponent',   &
             ' in evolution')

      END SUBROUTINE

      END MODULE siesta_init
