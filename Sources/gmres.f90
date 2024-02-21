!>  \file gmres.f90
!>  \brief Contains the \ref gmres module.
!>  \authors S. P. Hirshman and S. K. Seal
!>  \date Jan, 2014
      MODULE gmres
      USE v3_utilities, ONLY: assert
      USE stel_kinds
      USE stel_constants, ONLY: zero, one
      USE descriptor_mod, ONLY: SIESTA_COMM, iam, nprocs
      USE shared_data, etak=>etak_tol
      USE shared_functions, ONLY: funct_island, LineSearch
      USE siesta_state, ONLY: update_state, update_state
      USE nscalingtools, ONLY: PARFUNCTISL, MPI_ERR, startglobrow,      &
                               endglobrow, PARSOLVER
      USE timer_mod
      USE mpi_inc

      IMPLICIT NONE

      CONTAINS

!>  \brief Uses gmres to find approximate solution of Ax=b (linearized MHD equations)
!>   that minimize the nonlinear MHD forces
      SUBROUTINE gmres_fun
      USE hessian, ONLY: levmarq_param, Dmin=>mindiag, Apply_Precond,   &
                         levmarq_param0, l_diagonal_only, gather_array, &
                         muPar
      USE siesta_namelist, ONLY: lresistive, ftol
      USE blocktridiagonalsolver_s, ONLY: RefactorHessian

!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: j, iter, istat
      REAL(dp), PARAMETER :: levscan(1) = (/16._dp/)           !!007C
!      REAL(dp), PARAMETER :: levscan(2) = (/4._dp, 16._dp/)   !!007B
      REAL(dp) :: wmhd, fmhd, skston, skstoff, eps, lm0, fsq_min,       &
                  xmax_lin, muparS, scale_fac
      REAL(dp), ALLOCATABLE :: xc1(:), xcmin(:), brhs(:)
      LOGICAL               :: bTerm
!-----------------------------------------------

!     SPH: ONLY works with 0 restart each time
!     EXTERNALLY, l_init_state = .TRUE.   

      xc = 0
      l_linearize = .FALSE.
      l_ApplyPrecon = .FALSE.    !LET GMRES APPLY PRECONDITIONER ITSELF
      muParS = muPar
      
!
!     STORE INITIAL UNPRECONDITIONED RESIDUE (BRHS=-GC AT INIT PT)
!     SO DEVIATIONS FROM THEM CAN BE COMPUTED REPEATEDLY
!
      l_getfsq = .TRUE.
      l_init_state = .TRUE.
      fmhd = fsq_total1               !!initial fsq BEFORE application of resistive diffusion

      ALLOCATE (brhs(neqs), stat=istat)
      CALL ASSERT(istat.eq.0,'ALLOCATION ERROR IN GMRES_FUN')
      l_getwmhd = .TRUE.


      CALL second0(skston)
      CALL funct_island
      brhs = -gc
      IF (PARFUNCTISL) THEN
         CALL gather_array(brhs)
      END IF
      l_getwmhd = .FALSE.
      CALL second0(skstoff)
      gmres_funct_island_time = gmres_funct_island_time + (skstoff - skston)

      wmhd = wtotal
      fsq_min = fsq_total1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! FAST CONVERGENCE - SPH 083016 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ALLOCATE(xc1(neqs), xcmin(neqs), stat=istat)
      CALL ASSERT(istat.eq.0,'Allocation error in GMRES')
      xcmin = xc

      iter = 0
      lm0 = levmarq_param
      scale_fac = 1

920   CONTINUE
!
!     FAST CONVERGENCE GUESS
!     ITERATIVELY SOLVE (A* - eps*I) x = b, WHERE A* = A + eps*I 
!     eps = lev_marq_param*Dmin
!     THEN x = x0 + x1 + ..., 
!     x0 = P*b   where P* = (A*)^-1
!     xn = eps*P*(x[n-1])
!
      LFAST: IF (.NOT.l_Diagonal_Only .AND. lcolscale) THEN

         IF (iam .EQ. 0 .AND. lverbose) THEN
            PRINT 925, Dmin
         END IF
925   FORMAT(' DMIN: ',1pe12.3,/,1x,'FAST CONVERGENCE SUMMARY',                &
             /,' -------------',/,1x,'ITER',7x,'FSQ_NL',                       &
             10x,'||X||',9x,'MAX|X|')
         CALL ASSERT(ALL(xc .EQ. zero), ' XC != 0')
         l_getfsq = .TRUE.
         l_linearize = .FALSE.
         CALL funct_island
         xc = -gc
         eps = ABS(levmarq_param)*Dmin
         IF (.not.l_ApplyPrecon) THEN
            CALL Apply_Precond(xc)
         END IF
         xc1 = xc
         DO j = 1, 100
            IF (j .GT. 1) THEN
               IF (eps .EQ. zero) THEN
                  GOTO 927
               END IF
               gc = xc1
               CALL Apply_Precond(gc)
               xc1 = eps*gc            !nth partial sum: store for next j iteration
               xc = xc + xc1
            END IF
            CALL funct_island

            IF (fsq_total1 .GT. fsq_min) THEN
               bTerm = (j .GT. 1)
            ELSE
               xcmin = xc
               bTerm = (fsq_total1 .GE. 0.95_dp*fsq_min)
               fsq_min = fsq_total1
            END IF
            IF (iam .EQ. 0 .AND. lverbose) THEN
               PRINT 912, j, fsq_total1, SQRT(SUM(xc*xc)), MAXVAL(ABS(xc))
            END IF
            IF (fsq_min .LE. ftol .OR.                                         &
                bTerm             .OR.                                         &
                eps .EQ. zero) THEN
               EXIT
            END IF
         END DO
      
!      ELSE IF (l_backslv) THEN
!         IF (.NOT.l_ApplyPrecon) CALL Apply_Precond(gc)
      END IF LFAST

 927  CONTINUE

      LGMRES: IF (fsq_min.GT.ftol .AND. fsq_min.GT.fmhd/10) THEN

!!!START GMRES CONVERGENCE
!
!     THIS IS STANDARD GMRES CALL
!
         CALL second0(skston)
         xc = xcmin
         CALL gmres_wrap (xc, brhs)
         CALL second0(skstoff)
         gmres_wrap_time=gmres_wrap_time+(skstoff-skston)

         IF (fsq_total1 .LT. fsq_min .OR. ALL(xcmin .EQ. zero)) THEN
            xcmin = xc
            fsq_min = fsq_total1
            IF (lm0 .GT. zero) THEN
               scale_fac = levmarq_param/lm0
            END IF
         END IF
!!!END GMRES CONVERGENCE
!FOR NOW, THIS IS ONLY IMPLEMENTED FOR PARSOLVER=TRUE: MUST IMPLEMENT
!RefactorHessian FOR LSCALAPACK=T
         IF (.FALSE.) THEN
!        IF (.NOT.l_Diagonal_Only .AND. PARSOLVER) THEN
            IF (fsq_min.GE.0.95_dp*fmhd .AND. iter.LT.SIZE(levscan)) THEN
               iter = iter+1
               levmarq_param = lm0*levscan(iter)
               IF (levmarq_param .LT. levm_ped) THEN
                  levmarq_param = 100*levm_ped
               END IF
               IF (levmarq_param .GE. levmarq_param0) THEN
                  levmarq_param = levmarq_param0/10._dp**iter
               END IF
               IF (iam.EQ.0 .AND. lverbose) THEN
                  PRINT 930,' Refactoring Hessian: LM = ',levmarq_param
               END IF
               CALL RefactorHessian(levmarq_param)
               GOTO 920
            END IF
         END IF
      END IF LGMRES

 912  FORMAT(i5,3(3x,1pe12.3))
 930  FORMAT (/,a,1pe12.3)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      DEALLOCATE(xc1)
      l_getfsq = .TRUE.

 1100 CONTINUE

      xmax_lin = MAXVAL(ABS(xcmin))
      IF (iam .EQ. 0) THEN
!         PRINT 1020,' ||X||-GMRES = ', SQRT(SUM(xcmin*xcmin)), ' MAX(|X|) = ',xmax_lin
         WRITE (unit_out, 1020) ' ||X||-GMRES = ', SQRT(SUM(xcmin*xcmin)), ' MAX(|X|) = ',xmax_lin
      END IF
 1020 FORMAT(2(a,1p,e11.3))

      fsq_lin = fsq_total
      l_linearize = .FALSE.

!UPDATE EVOLUTION FIRST
      IF (fsq_min .GT. ftol) THEN

!     TURN OFF NEXT 2 LINES FOR TESTING f1=levm IN EVOLUTION
         levm_scale = levm_scale*scale_fac  
         levm_scale = MIN(100._dp, levm_scale)      
!
!     SIMPLE LINE SEARCH SCAN OF |F| FOR VARIOUS |X| 
         muPar = 0
         CALL LineSearch(xcmin, fsq_min)

      END IF

      xc = xcmin
!
!     RECOMPUTE PERTURBED FIELDS, PRESSURE FOR MINIMUM STATE
!     THEN UPDATE NONLINEAR STATE BY ADDING XC (PERTURBATION)
!

      l_PrintOriginForces=.TRUE.
      l_init_state=.TRUE.
      l_getwmhd= .TRUE.
      CALL second0(skston)
      CALL funct_island
!      CALL update_state(.FALSE., zero, zero)                           !Get this in evolution, line 455
      CALL second0(skstoff)
      gmres_funct_island_time=gmres_funct_island_time+(skstoff-skston)
      l_PrintOriginForces=.FALSE.

      muPar = muParS

      DEALLOCATE (xcmin, brhs)

      END SUBROUTINE gmres_fun

      SUBROUTINE init_gmres0 (icntl, cntl, etak, m, n)
      INTEGER, INTENT(OUT) :: icntl(9), m, n
      REAL(dp),INTENT(OUT) :: cntl(5)
      REAL(dp), INTENT(IN) :: etak
      REAL(dp)             :: skston, skstoff
      INTEGER, PARAMETER   :: noPrec=0, leftPrec=1, rightPrec=2, dblePrec=3

      CALL second0(skston)
      CALL init_dgmres(icntl ,cntl)
      CALL second0(skstoff)
      gmres_init_dgmres_time=gmres_init_dgmres_time+(skstoff-skston)

!*************************
!*  Tune some parameters
!*************************
! Tolerance
       cntl(1) = etak
!       cntl(1) = 1.E-5_dp
! Write errors to fort.21
!      icntl(1) = 21   !21
! Write warnings to fort.21
      icntl(2) = 21
! Save the convergence history in file fort.20
      IF (nprecon .gt. 1) icntl(3) = 20
! Preconditioning INPUT flags (note: different than revcom flags in driver)
!      icntl(4) = leftPrec
      icntl(4) = rightPrec 
!      icntl(4) = dblePrec
!! ICGS orthogonalization
      icntl(5) = 3
! Initial guess
      icntl(6) = 0
!      icntl(6) = 1
! Maximum number of iterations at each step (~ns/5)
      icntl(7) = ngmres_steps

      icntl(8) = 1             !Default
      icntl(9) = 1             !Steps for peek at progress during rev com loop
!*********************************
! Choose the restart parameter
!*********************************
!      write(*,*) 'Restart  <', ldstrt
!      read(*,*) m
!
!     m <= n
!
      n   = neqs
      m = 200

      END SUBROUTINE init_gmres0

      SUBROUTINE gmres_wrap (x0, b)
      USE siesta_namelist, ONLY: ftol
      USE nscalingtools, ONLY: PARGMRES, rcounts, disp
      USE hessian, ONLY: apply_precond
      USE gmres_lib, ONLY: gmres_info, gmres_par, gmres_ser
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      REAL(dp), INTENT(IN)    :: b(:)
      REAL(dp), INTENT(INOUT) :: x0(:)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      TYPE(gmres_info)   :: gi
      INTEGER :: n
!-----------------------------------------------
      l_linearize = .TRUE.
      l_getfsq = .FALSE.         !So call to funct in matvec does not include gc0

!******************************************************
!*  Initialize GMRES control parameters and Load GMRES_INFO structure
!*******************************************************
      CALL init_gmres0(gi%icntl, gi%cntl, etak, gi%m, n)
      gi%ftol=ftol
      gi%lverbose = lverbose
      gi%icntl(6) = 1

!
!     EASY-TO-USE WRAPPER FOR GMRES DRIVER CALL
!
!     X0: on input, initial guess if icntl(6) == 1
!         on output, solution of Ax = b
!         NOTE: it is not overwritten UNTIL the end of this routine
      gi%iam=iam; gi%nprocs=nprocs; gi%ngmres_type=ngmres_type
      IF (PARSOLVER .AND. PARGMRES) THEN
         gi%endglobrow=endglobrow; gi%startglobrow=startglobrow
         gi%mblk_size=mblk_size; gi%rcounts=>rcounts; gi%disp=>disp

         gi%my_comm = SIESTA_COMM
         gi%my_comm_world = SIESTA_COMM
!         gi%lactive = lactive
         CALL gmres_par (n, gi, matvec_par, apply_precond, GetNLForce,  &
                         x0, b)
      ELSE
         CALL gmres_ser (n, gi, matvec, apply_precond, GetNLForce,      &
                         x0, b)
      END IF

      fsq_total1 = gi%ftol
      CALL ASSERT(ALL(xc.EQ.x0),'XC != X0 IN GMRES_WRAP')

      END SUBROUTINE gmres_wrap

      SUBROUTINE matvec_par (ploc, Ap, nloc)
      USE island_params, ONLY: ns=>ns_i
      USE hessian, ONLY: eps_factor, mupar, gather_array
      USE blocktridiagonalsolver_s, ONLY: ParMatVec, PadSides 
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER, INTENT(IN)                    :: nloc
      REAL(dp), INTENT(IN), DIMENSION(nloc)  :: ploc
      REAL(dp), INTENT(OUT), DIMENSION(nloc) :: Ap
!-----------------------------------------------
      LOGICAL, PARAMETER    :: LPARMATVEC=.FALSE.
      INTEGER               :: myrowstart, myrowend, istat
      REAL(dp), ALLOCATABLE :: p(:)
      REAL(dp)              :: delta, skston, skstoff
!-----------------------------------------------
!SPH NOTE 032713: EVENTUALLY REMOVE CALL TO ParMatVec IF THE FUNCT_ISLAND IS FASTER
!
!     NOTE: DO NOT CALL ParMatVec WHEN nprecon_type==PREDIAG
!           SINCE IN HESSIAN, WE SET ALL THE OFF-DIAGONAL COMPONENTS OF 
!           THE A-MATRIX  (USED IN ParMatVec) TO ZERO
!
!           AT PRESENT, ParMatVec WILL NOT WORK IF mupar != 0
!
      CALL ASSERT(.NOT.l_init_state,'l_init_state = T in matvec_par')
      istat = (endglobrow-startglobrow+1)*mblk_size
      CALL ASSERT(istat.EQ.nloc, 'nloc wrong in matvec_par')

      myrowstart=(startglobrow-1)*mblk_size+1
      myrowend=myrowstart+nloc-1

      CALL second0(skston)
      ALLOCATE (p(ns*mblk_size), stat=istat)
      CALL ASSERT(istat.eq.0,'Allocation error in matvec_par')
      p = 0

      p(myrowstart:myrowend) = ploc
      CALL PadSides(p, mblk_size, 1, 1)

      CALL second0(skstoff)
      gmres_wrap_allgather_time=gmres_wrap_allgather_time+(skstoff-skston)

      CALL second0(skston)

      IF (.NOT.LPARMATVEC .OR. nprecon_type.EQ.PREDIAG                  &
          .OR. mupar.NE.zero) THEN
         delta = SQRT(EPSILON(delta))*eps_factor
         xc = delta*p

         l_linearize=.TRUE.
         CALL funct_island
         Ap = gc(myrowstart:myrowend)/delta

      ELSE 
         CALL ParMatVec(p,Ap,nloc)
      END IF

      DEALLOCATE (p)

      WHERE (ABS(Ap) .LT. 1.E-10_dp) Ap = 0

      CALL second0(skstoff)

      ParyAx_time=ParyAx_time+(skstoff-skston)

      END SUBROUTINE matvec_par

      SUBROUTINE matvec (p, Ap, ndim)
      USE stel_kinds
      USE hessian, ONLY: eps_factor, gather_array
#if defined(MPI_OPT)
      USE mpi_inc
#endif
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER, INTENT(IN)   :: ndim
      REAL(dp), INTENT(IN)  :: p(ndim)
      REAL(dp), INTENT(OUT) :: Ap(ndim)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      REAL(dp), PARAMETER   :: zero=0
      REAL(dp)              :: delta
!-----------------------------------------------
      CALL ASSERT(.not.l_init_state,'l_init_state = T in matvec')

      delta = SQRT(EPSILON(delta))*eps_factor

!     Note: GMRES without Preconditioner has pnorm = 1 (checked!)
!      pnorm = SUM(p*p)
!      delta = delta/SQRT(pnorm)
!      delta = delta*(pnorm + SUM(p*xc0))/pnorm

!IF CALLED IN PAR MODE, FIRST GATHER THE xc's
      CALL ASSERT(SIZE(gc).EQ.SIZE(Ap),'gc and Ap wrong sizes')
      xc = delta*p

      l_getwmhd=.TRUE.
      CALL funct_island
      l_getwmhd=.FALSE.

      IF (l_linearize) THEN
         Ap = gc/delta
      ELSE
         Ap = (gc-gc0)/delta
      END IF

      WHERE (ABS(Ap) .LT. 1.E-10_dp) Ap = 0

      END SUBROUTINE matvec


      SUBROUTINE get_etak(fk, fkm, ftol, etak)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      REAL(dp), INTENT(in)    :: fk, fkm, ftol
      REAL(dp), INTENT(inout) :: etak
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      REAL(dp), PARAMETER :: gamma1=0.9_dp, alpha=1.5_dp, eta0 = 1.E-02_dp
      REAL(dp)            :: etakm
!-----------------------------------------------
!     ROUTINE PROVIDED L. CHACON (11/09/06)
      etakm = etak

!     Superlinear convergence
      etak = gamma1*(fk/fkm)**alpha

!     Safeguard avoid sharp decrease in etak
      etak = MIN(eta0, MAX(etak, gamma1*etakm**alpha))

!     Safeguard avoid "oversolving"
      etak = MIN(eta0, MAX(etak, gamma1*ftol/fk))

      END SUBROUTINE get_etak
 

      SUBROUTINE GetNLForce(xcstate, fsq_nl, bnorm)
#if defined(MPI_OPT)
      USE blocktridiagonalsolver_s, ONLY: PadSides 
#endif
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      REAL(dp),INTENT(IN)  :: xcstate(neqs), bnorm
      REAL(dp),INTENT(OUT) :: fsq_nl
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER              :: nloc, myrowstart, myrowend
      LOGICAL              :: l_linloc, l_Apploc, l_getloc
!-----------------------------------------------
      CALL ASSERT(.not.l_init_state,'l_init_state=T in GetNLForce!')
!Store state 
      l_linloc = l_linearize
      l_Apploc = l_ApplyPrecon
      l_getloc = l_getfsq
      
      nloc=(endglobrow-startglobrow+1)*mblk_size
      myrowstart=(startglobrow-1)*mblk_size+1
      myrowend=myrowstart+nloc-1


!Set state variables
      xc(myrowstart:myrowend) = bnorm*xcstate(myrowstart:myrowend)       !undo internal gmres normalization
#if defined(MPI_OPT)
      IF (PARSOLVER) CALL PadSides(xc, mblk_size, 1, 1)
#endif
      l_init_state = .FALSE.
      l_linearize = .FALSE.
      l_getfsq = .TRUE.
      l_ApplyPrecon=.FALSE.

      CALL funct_island
      fsq_nl = fsq_total1
 
!Restore state variables
      l_linearize = l_linloc
      l_ApplyPrecon = l_Apploc
      l_getfsq = l_getloc

      END SUBROUTINE GetNLForce

      SUBROUTINE qmr_fun
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: ndim, nlen, nlim, ierr, info(4), j
      INTEGER :: revcom, colx, colb
      REAL(dp), DIMENSION(:,:), ALLOCATABLE :: vecs
      REAL(dp) :: tol = 1.E-4_dp, gnorm
!-----------------------------------------------
!
!     STORE INITIAL POINT AND INITIAL RESIDUE (AT INIT PT)
!     SO DEVIATIONS FROM THEM CAN BE COMPUTED REPEATEDLY
!
      ALLOCATE(xc0(neqs), gc0(neqs), vecs(neqs,9), stat=j)
      CALL ASSERT(j.eq.0,'Allocation error in qmr_fun')

      ndim = SIZE(vecs,1)
      nlen = ndim
      nlim = 100

!     SPH: ONLY works with 0 restart each time??!!!     
      l_linearize = .FALSE.
      l_ApplyPrecon = .TRUE.

      xc = 0
      CALL funct_island
      xc0 = xc
      gc0 = gc

      l_linearize = .TRUE.
!
!     INITIALIZE vecs
!
      gnorm = SUM(gc*gc)
      gnorm = SQRT(gnorm)
      vecs(:ndim,2) = -gc(:ndim)/gnorm
      vecs(:ndim,3) =  gc(:ndim)/gnorm
    
      ierr = 100000
      info = 0
      info(1) = ierr

 10   CALL dutfx (ndim,nlen,nlim,vecs,tol,info)
      revcom = info(2)
      colx   = info(3)
      colb   = info(4)
      IF (revcom .eq. 1) THEN
         CALL matvec (vecs(1,colx), vecs(1,colb), ndim)
         GO TO 10
      END IF

      xc(1:ndim) = xc0(1:ndim) + gnorm*vecs(:,1)

!
!     RECOMPUTE delta current,pressure prior to update_state call
!
      CALL funct_island
      fsq_lin = fsq_total

      l_linearize = .FALSE.
      CALL funct_island

      DEALLOCATE (xc0, gc0, vecs)

      END SUBROUTINE qmr_fun

      END MODULE gmres
