      MODULE hessian
      USE island_params, mpol=>mpol_i, ntor=>ntor_i, ns=>ns_i,          &
          nfp=>nfp_i
      USE timer_mod
      USE stel_constants, ONLY: pi
      USE descriptor_mod
      USE shared_data, ONLY: col_scale, neqs, lcolscale, l_ApplyPrecon, &
                             l_linearize, l_getfsq, lverbose,           &
                             ndims, mblk_size, unit_out, nprecon
#if defined(MPI_OPT)
      USE prof_mod
      USE ptrd_mod
#endif
      USE nscalingtools, ONLY: PARSOLVER, MPI_ERR, PARFUNCTISL,         &
        rcounts, disp, startglobrow, endglobrow, nranks, rank,          &
        SKSDBG, TOFU, UPPER, LOWER, DIAG, SAVEDIAG, nrecd,              &
        SYMMETRYCHECK, totRecvs, PACKSIZE, WriteTime, send, receive,    &
        GetFullSolution
      USE blocktridiagonalsolver_s, ONLY: ApplyParallelScaling,         &
        ForwardSolve, SetMatrixRHS, BackwardSolve, GetSolutionVector,   &
        CheckSymmetry, Dmin_TRI, MAXEIGEN_TRI, FindMinMax_Tri,          &
        CheckConditionNumber
      USE mpi_inc
      USE v3_utilities, ONLY: assert, assert_eq

      IMPLICIT NONE

!-----------------------------------------------
!SPH 090116: EVENTUALLY MOVE TO nscalingtools? ASK SKS
      PUBLIC  :: gather_array
      INTERFACE gather_array
         MODULE PROCEDURE gather_array1, gather_array2, gather_array3
      END INTERFACE
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      PRIVATE
      INTEGER, PARAMETER   :: jstart(3) = (/1,2,3/)
      INTEGER, ALLOCATABLE :: ipiv_blk(:,:)

      REAL(dp), ALLOCATABLE, DIMENSION(:,:,:) :: ublk, dblk, lblk
      REAL(dp), ALLOCATABLE, DIMENSION(:,:,:) :: ublk_s, dblk_s, lblk_s
      LOGICAL, PUBLIC :: l_backslv = .FALSE.

      REAL(dp), ALLOCATABLE, DIMENSION(:) :: DataItem

      INTEGER :: myStartBlock, myEndBlock, mynumBlockRows 
      INTEGER :: mystart(3), myend(3)

      REAL(dp) :: eps, ton, toff, rcond, anorm

      REAL(dp), DIMENSION(:,:), ALLOCATABLE  :: colsum
      REAL(dp), PARAMETER, PUBLIC :: eps_factor = 1._dp 
      REAL(dp), ALLOCATABLE, PUBLIC :: mupar_norm(:)
      REAL(dp), PUBLIC  :: levmarq_param, levmarq_param0, asym_index,   &
                           mupar, mupar0, maxeigen, mindiag
      LOGICAL, PUBLIC   :: l_Compute_Hessian, l_Diagonal_Only
      LOGICAL, PUBLIC, PARAMETER :: l_Hess_sym=.FALSE.                     !<forces symmetrization of Hessian if TRUE (only implemented for Thomas, 1 proc)
      INTEGER, PUBLIC :: HESSPASS=0
      PUBLIC  :: apply_precond, dealloc_hessian, InitHess,              &
                 Compute_Hessian_Blocks, apply_colscale

      CONTAINS

!-------------------------------------------------------------------------------
!>  @brief Gather all parts of an 1D array.
!>
!>  @param[inout] buffer Buffer to gather to and from.
!-------------------------------------------------------------------------------
      SUBROUTINE gather_array1(buffer)

      IMPLICIT NONE

!  Declare Arguments
      REAL(dp), DIMENSION(:), INTENT(inout) :: buffer

!  Start of executable code
      IF (PARSOLVER) THEN
         CALL MPI_ALLGATHERV(MPI_IN_PLACE, rcounts(iam + 1), MPI_REAL8,        &
                             buffer, rcounts, disp, MPI_REAL8, SIESTA_COMM,    &
                             MPI_ERR)
      END IF

      END SUBROUTINE

!-------------------------------------------------------------------------------
!>  @brief Gather all parts of an 1D array.
!>
!>  @param[inout] buffer Buffer to gather to and from.
!-------------------------------------------------------------------------------
      SUBROUTINE gather_array2(buffer)

      IMPLICIT NONE

!  Declare Arguments
      REAL(dp), DIMENSION(:,:), INTENT(inout) :: buffer

!  Start of executable code
      IF (PARSOLVER) THEN
         CALL MPI_ALLGATHERV(MPI_IN_PLACE, rcounts(iam + 1), MPI_REAL8,        &
                             buffer, rcounts, disp, MPI_REAL8, SIESTA_COMM,    &
                             MPI_ERR)
      END IF

      END SUBROUTINE

!-------------------------------------------------------------------------------
!>  @brief Gather all parts of an 3D array.
!>
!>  @param[inout] buffer Buffer to gather to and from.
!-------------------------------------------------------------------------------
      SUBROUTINE gather_array3(buffer)
      USE nscalingtools, ONLY: mnspcounts, mnspdisps

      IMPLICIT NONE

!  Declare Arguments
      REAL (dp), DIMENSION(:,:,:), INTENT(inout) :: buffer

!  Start of executable code
      IF (PARSOLVER) THEN
         CALL MPI_ALLGATHERV(MPI_IN_PLACE, mnspcounts(iam + 1), MPI_REAL8,     &
                             buffer, mnspcounts, mnspdisps, MPI_REAL8,         &
                             SIESTA_COMM, MPI_ERR)
      END IF

      END SUBROUTINE

!-------------------------------------------------------------------------------
!>  @brief Gather column arrays.
!>
!>  @param[inout] colgath Gathered columns.
!-------------------------------------------------------------------------------
      SUBROUTINE GatherCols(colgath)

      IMPLICIT NONE

!  Declare Arguments
      REAL(dp), DIMENSION(mblk_size), INTENT(inout) :: colgath

#if defined(MPI_OPT)
!  Start of executable code
      CALL MPI_ALLREDUCE(MPI_IN_PLACE, colgath, mblk_size, MPI_REAL8,          &
                         MPI_SUM, SIESTA_COMM, MPI_ERR)
      CALL ASSERT(MPI_ERR.EQ.0, 'MPI ERR IN GatherCols')
#endif
      END SUBROUTINE

      SUBROUTINE InitHess
      INTEGER  :: icol, mesh, nsmin, nsmax
!-----------------------------------------------
      IF (.NOT.PARSOLVER) THEN
         startglobrow = 1; endglobrow = ns
      END IF
      myStartBlock=startglobrow
      myEndBlock=endglobrow
      mynumBlockRows=myEndBlock-myStartBlock+1

      nsmin = MAX(1, startglobrow-1);  nsmax = MIN(ns, endglobrow+1)
      myend = nsmax
      DO mesh = 1, 3
         icol = MOD(jstart(mesh)-nsmin, 3)
         IF (icol .LT. 0) icol = icol+3
         mystart(mesh) = nsmin+icol
      END DO

      END SUBROUTINE InitHess

!-------------------------------------------------------------------------------
!>  @brief Scale displacements.
!>
!>  @param[inout] gc       Scratch buffer to scale columns in.
!>  @param[in]    colscale Column scaling factors.
!>  @param[in]    nsmin    Minimum radial index.
!>  @param[in]    nsmax    Maximum radial index.
!-------------------------------------------------------------------------------
      SUBROUTINE Apply_ColScale(gc, colscale, nsmin, nsmax)

      IMPLICIT NONE

!  Declare Arguments
      REAL (dp), DIMENSION(mblk_size,ns), INTENT(inout) :: gc
      REAL (dp), DIMENSION(mblk_size,ns), INTENT(in)    :: colscale
      INTEGER, INTENT(in)                               :: nsmin
      INTEGER, INTENT(in)                               :: nsmax

!  Start of executable code
      gc(:,nsmin:nsmax) = gc(:,nsmin:nsmax)*colscale(:,nsmin:nsmax)
            
      END SUBROUTINE

      SUBROUTINE Compute_Hessian_Blocks (func, ldiagonal)
      USE shared_data, ONLY: xc, gc, l_linearize
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      LOGICAL, INTENT(IN)     :: ldiagonal
      EXTERNAL                :: func
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!----------------------------------------------
      INTEGER                 :: iunit
      CHARACTER(LEN=3)        :: label
      REAL(dp)                :: bsize
      LOGICAL                 :: lprint
!-----------------------------------------------
      l_Compute_Hessian = .TRUE.
      l_Diagonal_Only   = ldiagonal
      l_linearize = .TRUE.
      IF (lColScale) THEN
         col_scale = 1
      END IF

      IF (iam .EQ. 0) THEN
         DO iunit = 6, unit_out, unit_out-6
            IF (.NOT.lverbose .AND. iunit.EQ.6) CYCLE
            IF (l_Diagonal_Only) THEN
               WRITE (iunit, 90) levmarq_param, muPar
            ELSE
               WRITE (iunit, 100) levmarq_param, muPar, asym_index
            ENDIF
         END DO
         CALL FLUSH(unit_out)
      ENDIF

 90   FORMAT (/,' Computing diagonal preconditioner - ',                &
                ' LM parameter:',1pe9.2,' mu||:',1pe9.2)
 100  FORMAT (/,' Computing block preconditioner - ',                   &
                ' LM parameter:',1pe9.2,' mu||:',1pe9.2,                &
                ' Asym index:',1pe9.2)

      lprint = (iam.EQ.0 .AND. HESSPASS.EQ.0 .AND. lverbose)

      IF (PARSOLVER) THEN
        IF (PARFUNCTISL) THEN
          IF(SKSDBG) WRITE(TOFU,*) 'Executing Compute_Hessian_Blocks_With_No_Col_Redist'; IF(SKSDBG) CALL FLUSH(TOFU)
          IF (lprint) PRINT *,'HESSIAN AND INVERSE CALCULATION' //      &
                          ' USING BCYCLIC WITH NO COLUMN REDISTRIBUTION'
          CALL Compute_Hessian_Blocks_With_No_Col_Redist (xc, gc, func)
        ELSE
          IF (SKSDBG) WRITE(TOFU,*) 'Executing Compute_Hessian_Blocks_With_Col_Redist'; IF(SKSDBG) CALL FLUSH(TOFU)
          IF (lprint) PRINT *,'HESSIAN AND INVERSE CALCULATION' //      &
                             ' USING BCYCLIC WITH COLUMN REDISTRIBUTION'
          CALL Compute_Hessian_Blocks_With_Col_Redist (xc, gc, func)
        END IF 
      ELSE 
        IF (SKSDBG) WRITE(TOFU,*) 'Executing Compute_Hessian_Blocks_Thomas'; CALL FLUSH(TOFU)
        IF (lprint) PRINT *,'HESSIAN AND INVERSE CALCULATION' //        &
                            ' USING SERIAL THOMAS ALGORITHM'
        CALL Compute_Hessian_Blocks_Thomas (xc, gc, func)
      END IF

      IF (HESSPASS.EQ.0) THEN
         bsize = 3*ns*KIND(dblk)                                         !3 blocks per row
         bsize = bsize*REAL(mblk_size*mblk_size,dp)
         IF (bsize .LT. 1.E6_dp) THEN
            bsize = bsize/1.E3
            label = " Kb"
         ELSE IF (bsize .lt. 1.E9_dp) THEN
            bsize = bsize/1.E6_dp
            label = " Mb"
         ELSE
            bsize = bsize/1.E9_dp
            label = " Gb"
         END IF

         IF (iam .EQ. 0) THEN
            DO iunit = 6, unit_out, unit_out-6
               IF (.NOT.lverbose .AND. iunit.EQ.6) CYCLE
               WRITE (iunit, '(1x,a,i4,a,f6.2,a)') 'Block dim: ',       &
                  mblk_size, '^2  Preconditioner size: ', bsize,        &
                  TRIM(label)
            END DO
            CALL FLUSH(unit_out)
         ENDIF 
      END IF

      l_Compute_Hessian = .FALSE.
      l_linearize = .FALSE.
      HESSPASS = HESSPASS+1

      END SUBROUTINE Compute_Hessian_Blocks

      SUBROUTINE Compute_Hessian_Blocks_With_No_Col_Redist(xc, gc, func)
      USE blocktridiagonalsolver_s, ONLY:                                      &
     & SetMatrixRowColL, SetMatrixRowColD, SetMatrixRowColU, StoreDiagonal
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      REAL(dp), DIMENSION(mblk_size,ns) :: xc
      REAL(dp), DIMENSION(mblk_size,ns) :: gc
      EXTERNAL                          :: func
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!----------------------------------------------
      INTEGER  :: n, m, js, mesh, ntype, istat, iunit, icol
      INTEGER  :: js1

      REAL(dp) :: starttime, endtime, usedtime
      REAL(dp) :: colstarttime, colendtime, colusedtime
      INTEGER  :: nsmin, nsmax,              i
      REAL(dp) :: skston, skstoff, temp

      REAL(dp), PARAMETER    :: zero = 0
      REAL(dp), PARAMETER    :: one = 1
!----------------------------------------------

      starttime=0; endtime=0; usedtime=0; 
      colstarttime=0; colendtime=0; colusedtime=0

      !------------------------------------------------------------
      !     COMPUTE (U)pper, (D)iagonal, and (L)ower Block matrices
      !------------------------------------------------------------

      CALL second0(ton)
      eps = SQRT(EPSILON(eps))*ABS(eps_factor)
!      eps =  1.E-3_dp
      
      nsmin = MAX(1, startglobrow-1);  nsmax = MIN(ns, endglobrow+1)

      ALLOCATE (DataItem(mblk_size), stat=istat)
      CALL ASSERT(istat.EQ.0,'DataItem allocation failed:')
      DataItem = 0

      xc(:,nsmin:nsmax) = 0

      !----------------------------
      ! Column generation begins
      !----------------------------
      CALL second0(colstarttime)

      icol = 0
      VAR_TYPE_LG: DO ntype = 1, ndims
         VAR_N_LG: DO n = -ntor, ntor
            VAR_M_LG: DO m = 0, mpol

               icol = icol + 1

               MESH_3PT_LG: DO mesh = 1, 3

                  IF (.NOT.(m.EQ.0 .AND. n.LT.0)) THEN
                     DO js = mystart(mesh), myend(mesh), 3
                        xc(icol,js) = eps
                     END DO
                  END IF

                  INHESSIAN = .TRUE.
                  CALL second0(skston)

                  CALL func
                  CALL second0(skstoff)
                  hessian_funct_island_time = hessian_funct_island_time + (skstoff - skston)
                  INHESSIAN = .FALSE.

                  SKIP3_MESH_LG: DO js = mystart(mesh), myend(mesh), 3

                     xc(icol, js) = 0

                     !ublk(js-1)
                     js1 = js - 1
                     IF (startglobrow.LE.js1 .AND. js1.LE.endglobrow) THEN
                        DataItem = gc(:,js1)/eps
                        IF (l_diagonal_only) THEN
                           temp = DataItem(icol)
                           DataItem = 0
                           DataItem(icol) = temp
                        END IF
                        CALL SetMatrixRowColU(js1, DataItem, icol)
                     END IF

                     !dblk(js)
                     IF (startglobrow .LE. js .AND. js .LE. endglobrow) THEN
                        DataItem = gc(:,js)/eps
                        IF (ALL(DataItem .EQ. zero)) THEN
                           DataItem(icol) = one
                        END IF

                        IF (l_diagonal_only) THEN
                           temp = DataItem(icol)
                           DataItem = 0
                           DataItem(icol)=temp
                        END IF

                        CALL StoreDiagonal(js, icol, DataItem)
                        CALL SetMatrixRowColD(js, DataItem, icol)
                     END IF

                     !lblk(js+1)
                     js1 = js + 1
                     IF (startglobrow .LE. js1 .AND. js1 .LE. endglobrow) THEN
                        DataItem = gc(:,js1)/eps
                        IF (l_diagonal_only) THEN
                           temp = DataItem(icol)
                           DataItem = 0
                           DataItem(icol) = temp
                        END IF
                        CALL SetMatrixRowColL(js1, DataItem, icol)
                     END IF
                  END DO SKIP3_MESH_LG
               END DO MESH_3PT_LG
            END DO VAR_M_LG
         END DO VAR_N_LG
      END DO VAR_TYPE_LG

      CALL second0(colendtime)
      construct_hessian_time=construct_hessian_time+(colendtime-colstarttime)

      DEALLOCATE (DataItem, stat=istat)
      CALL ASSERT(istat.EQ.0,'DataItem deallocation error:')
!----------------------------
! Column generation ends
!----------------------------

      IF (lColScale) THEN
         CALL ApplyParallelScaling(levmarq_param, col_scale, lverbose)
      ELSE 
         CALL FindMinMax_TRI (levmarq_param)
      END IF

      mindiag = Dmin_TRI
      maxeigen = MAXEIGEN_TRI
!
! CHECK CONDITION NUMBER
!
      IF (nprocs .EQ. 1) THEN
         CALL second0(ton)
         CALL CheckConditionNumber(ns,mblk_size,anorm,rcond,info)
         CALL second0(toff)                  
         IF (INFO .EQ. 0 .and. lverbose) THEN
            PRINT '(1x,3(a,1p,e12.3))','RCOND = ', rcond,               &
            ' ||A|| = ', ANORM,' TIME: ', toff-ton
         END IF
      END IF

!
! CHECK SYMMETRY OF BLOCKS 
!
      IF (SYMMETRYCHECK) THEN
         CALL second0(skston)
         CALL CheckSymmetry(asym_index)
         CALL second0(skstoff)
         asymmetry_check_time=asymmetry_check_time+(skstoff-skston)
      ENDIF

      CALL second0(toff)
      IF (l_Diagonal_Only) THEN
         time_diag_prec = time_diag_prec+(toff-ton)
      ELSE
         time_block_prec = time_block_prec+(toff-ton)
      END IF

!      IF (l_Diagonal_Only) RETURN         !ForwardSolve will be called first time but safer this way!

      ton = toff
      skston = ton
!
!FACTORIZE (Invert) HESSIAN
!
      IF (ALLOCATED(ipiv_blk)) DEALLOCATE(ipiv_blk,stat=istat)
      CALL ForwardSolve
      CALL second0(toff)
      skstoff=toff
      block_factorization_time=block_factorization_time+(skstoff-skston)
      time_factor_blocks = time_factor_blocks + (toff-ton)

      IF (iam.EQ.0 .AND. lverbose)                                      &
         PRINT '(a,1p,e12.3)',' BLOCK FACTORIZATION TIME: ', toff-ton

      END SUBROUTINE Compute_Hessian_Blocks_With_No_Col_Redist


      SUBROUTINE Compute_Hessian_Blocks_With_Col_Redist(xc, gc, func)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      REAL(dp), DIMENSION(mblk_size,ns) :: xc, gc
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!----------------------------------------------
      REAL(dp), PARAMETER    :: zero=0, one=1
      INTEGER  :: n, m, js, js1, mesh, ntype, istat, iunit, icol, icolmpi
      EXTERNAL    :: func
      REAL(dp), ALLOCATABLE, DIMENSION(:) :: SavedDiag
      INTEGER, ALLOCATABLE, DIMENSION (:) :: sendCount, recvCount
      REAL(dp), ALLOCATABLE, DIMENSION (:) :: recvBuf
      INTEGER :: procID
      REAL(dp) :: starttime, endtime, usedtime
      REAL(dp) :: colstarttime, colendtime, colusedtime

      LOGICAL :: PROBEFLAG
      REAL(dp) :: skston, skstoff, temp

      INTEGER :: it 
!----------------------------------------------
      starttime=zero
      endtime=zero
      usedtime=zero
!------------------------------------------------------------
!     COMPUTE (U)pper, (D)iagonal, and (L)ower Block matrices
!------------------------------------------------------------

      CALL second0(ton)

      CALL ASSERT(mblk_size.EQ.SIZE(gc,1),'MBLK_SIZE IS WRONG IN HESSIAN')

      eps = SQRT(EPSILON(eps))*ABS(eps_factor)

      ALLOCATE (DataItem(mblk_size), stat=istat)
      CALL ASSERT(istat.EQ.0,'DataItem allocation failed!')
      DataItem = 0
      ALLOCATE (SavedDiag(mblk_size), stat=istat)
      CALL ASSERT(istat.EQ.0,'SavedDiag allocation failed!')
      SavedDiag = zero

      icolmpi = 0
      icol=0
      xc = 0

      IF(.NOT.ALLOCATED(sendCount)) ALLOCATE(sendCount(nprocs))
      IF(.NOT.ALLOCATED(recvCount)) ALLOCATE(recvCount(nprocs))
      IF(.NOT.ALLOCATED(recvBuf)) ALLOCATE(recvBuf(PACKSIZE))
      sendCount=0
      recvCount=0
      nrecd=0
      totRecvs=0
      PROBEFLAG=.TRUE.

      CALL second0(colstarttime)
      VAR_TYPE: DO ntype = 1, ndims
         VAR_N: DO n = -ntor, ntor
            VAR_M: DO m = 0, mpol
                  
               icol=icol+1
               IF(MOD(icol-1,nprocs)==iam) THEN
                  icolmpi = icolmpi+1

                  MESH_3PT: DO mesh = 1, 3

                  IF (.NOT.(m.EQ.0 .AND. n.LT.0)) THEN
                     DO js = jstart(mesh), ns, 3
                        xc(icol,js) = eps
                     END DO
				  END IF

                  INHESSIAN=.TRUE.
                  CALL second0(skston)
                  CALL func
                  CALL second0(skstoff)
                  hessian_funct_island_time=hessian_funct_island_time+(skstoff-skston)
                  INHESSIAN=.FALSE.

                  SKIP3_MESH: DO js = jstart(mesh), ns, 3

                     xc(icol,js) = 0

                     !ublk(js-1)
                     js1 = js-1
                     IF (js1 .gt. 0) THEN
                        DataItem = gc(:,js1)/eps
                        !m=0 constraint for n<0
                        IF (m.eq.0 .and. n.lt.0) THEN
                           DataItem = 0
                        END IF
                        IF (l_diagonal_only) THEN
                           temp = DataItem(icol)
                           DataItem = 0
                           DataItem(icol) = temp
                        END IF
                        istat = 1
                        CALL receive(PROBEFLAG)
                        CALL send(DataItem, js1, UPPER, icol, procID)
                        IF (procID-1.NE.iam) THEN
                           sendCount(procID) = sendCount(procID) + 1
                        END IF
                        CALL receive(PROBEFLAG)
                     END IF !js>1

                     !dblk(js)
                     DataItem = gc(:,js)/eps
                     IF (m.eq.0 .and. n.lt.0) THEN
                        DataItem = 0
                     END IF
                     IF (ALL(DataItem .EQ. zero)) THEN
                        DataItem(icol) = one
                     END IF
                     
                     IF (l_diagonal_only) THEN
                        temp=DataItem(icol)
                        DataItem=0
                        DataItem(icol)=temp
                     END IF
                     SavedDiag = DataItem
                     CALL receive(PROBEFLAG)
                     CALL send(SavedDiag, js, SAVEDIAG, icol, procID)
                     IF(procID-1.NE.iam) sendCount(procID) = sendCount(procID)+1
                     CALL receive(PROBEFLAG)
!                    Boundary condition at js=1 and ns (CATCH ALL OF THEM HERE)
!                    and m=0,n<0: NEED THIS to avoid (near) singular Hessian
!                    ASSUMES DIAGONALS ARE ALL NEGATIVE
                     js1 = js
                     CALL receive(PROBEFLAG)
                     CALL send(DataItem, js1, DIAG, icol, procID)
                     IF(procID-1.NE.iam) sendCount(procID) = sendCount(procID)+1
                     CALL receive(PROBEFLAG)
                   
                     !lblk(js+1)
                     js1 =js+1 
                     IF (js .lt. ns) THEN
                       DataItem = gc(:,js1)/eps
                       !m=0 constraint for n<0
                       IF (m.eq.0 .and. n.lt.0) DataItem=0
                       IF (l_diagonal_only) THEN
                          temp=DataItem(icol)
                          DataItem=0
                          DataItem(icol)=temp
                       END IF
                       CALL receive(PROBEFLAG)
                       CALL send(DataItem, js1, LOWER, icol, procID)
                       IF(procID-1.NE.iam) sendCount(procID) = sendCount(procID)+1
                       CALL receive(PROBEFLAG)
                     END IF !JS < NS

                  END DO SKIP3_MESH

               END DO MESH_3PT
               ENDIF
            END DO VAR_M
         END DO VAR_N
      END DO VAR_TYPE
      CALL second0(colendtime)
      colusedtime=colendtime-colstarttime

      construct_hessian_time=construct_hessian_time+(colendtime-colstarttime)

      CALL second0(skston)
      IF (PARSOLVER) THEN

        PROBEFLAG=.NOT.PROBEFLAG

        CALL MPI_AllTOALL(sendCount,1,MPI_INTEGER,recvCount,1,          &
		                  MPI_INTEGER,SIESTA_COMM,MPI_ERR)
        totRecvs=0
        DO js=1,nprocs,1
          totRecvs=totRecvs+recvCount(js)
        END DO

        DO WHILE (nrecd.LT.totRecvs)
          CALL receive(PROBEFLAG)
        END DO
      END IF

      DEALLOCATE (DataItem, SavedDiag, stat=istat)
      CALL ASSERT(istat.EQ.0,'DataItem deallocation error!')

      CALL second0(skstoff)
      construct_hessian_time=construct_hessian_time+(skstoff-skston)
      CALL MPI_Barrier(SIESTA_COMM, MPI_ERR)

!
!     CHECK SYMMETRY OF BLOCKS 
!
      IF (SYMMETRYCHECK) THEN
        IF (PARSOLVER) THEN 
          CALL second0(skston)
          CALL CheckSymmetry(asym_index)
          CALL second0(skstoff)
          asymmetry_check_time=asymmetry_check_time+(skstoff-skston)
        ENDIF
      ENDIF      

!
!     FACTORIZE (Invert) HESSIAN
!
      CALL second0(toff)
      IF (l_Diagonal_Only) THEN
         time_diag_prec = time_diag_prec+(toff-ton)
      ELSE
         time_block_prec = time_block_prec+(toff-ton)
      END IF

      ton=toff
      starttime=ton
      skston=ton
      CALL ForwardSolve
      CALL second0(toff)
      endtime=toff
      skstoff=toff

      block_factorization_time=block_factorization_time+(skstoff-skston)
      time_factor_blocks = time_factor_blocks + (toff-ton)
      usedtime=endtime-starttime

      IF (iam.EQ.0 .AND. lverbose)                                      &
         PRINT '(a,1p,e12.3)',' BLOCK FACTORIZATION TIME: ', toff-ton

      END SUBROUTINE Compute_Hessian_Blocks_With_Col_Redist


      SUBROUTINE Compute_Hessian_Blocks_Thomas (xc, gc, func)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      REAL(dp), DIMENSION(mblk_size,ns), INTENT(INOUT)  :: xc, gc
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!----------------------------------------------
      REAL(dp), PARAMETER    :: zero=0, one=1
      INTEGER  :: n, m, js, js1, mesh, ntype, istat, iunit, icol, icolmpi
      REAL(dp), ALLOCATABLE, DIMENSION(:,:)  :: gc1
      REAL(dp) :: skston, skstoff, temp
      EXTERNAL :: func
#if defined(MPI_OPT)
      INTEGER  :: neqtotal, numroc 
      REAL(dp) :: starttime, endtime, usedtime
      REAL(dp) :: colstarttime, colendtime, colusedtime
      EXTERNAL :: numroc
!------------------------------------------------------------

      starttime=0; endtime=0; usedtime=0
      colstarttime=0; colendtime=0; colusedtime=0
#endif

!------------------------------------------------------------
!     COMPUTE (U)pper, (D)iagonal, and (L)ower Block matrices
!------------------------------------------------------------
      CALL second0(ton)

      eps = SQRT(EPSILON(eps))*ABS(eps_factor)

      IF (ALLOCATED(ublk)) DEALLOCATE(ublk, dblk, lblk, stat=istat)
      IF (LSCALAPACK) THEN
#if defined(MPI_OPT)
        CALL blacs_gridinfo(icontxt_1xp,nprow,npcol,myrow,mycol)
        mb = mblk_size
        nb = 1

        rsrc = 0
        csrc = 0
        Locq = numroc( mblk_size, nb, mycol, csrc, npcol )
        Locp = numroc( mblk_size, mb, myrow, rsrc, nprow )
        mblk_size2=max(1,Locq)
        lld = max(1,Locp)
        call descinit(descA_1xp,mblk_size,mblk_size,mb,nb,rsrc,csrc,       &
                      icontxt_1xp,lld,info)
        if (info.NE.0) then
          write(*,*) 'myrow,mycol,nprow,npcol,desc(LLD_),lld ',            &
             myrow,mycol,nprow,npcol,descA_1xp(LLD_),lld
          write(*,*) 'Locp,m,mb ', Locp,mblk_size,mb
        endif

        call assert(info.eq.0,'descinit descA_1xp')
        ineed = max(1,lld*Locq)
        mm = mblk_size*ns
        mb=mm
        nb=nrhs1
        LocqR = numroc( nrhs1, nb, mycol,csrc,npcol)
        LocpR = numroc( mm, mb, myrow,rsrc,nprow)
        lld = max(1,LocpR)
        call descinit(descR_all,mm,nrhs1, mb,nb,rsrc,csrc,icontxt,lld,info)
        call assert(info.eq.0,'test_pdtrd:descinit return info != 0')

        call blacs_gridinfo(icontxt,nprow,npcol,myrow,mycol)

        mb = 50
        nb = mb

        rsrc = 0
        csrc = 0
        Locq = numroc( mblk_size, nb, mycol, csrc, npcol )
        Locp = numroc( mblk_size, mb, myrow, rsrc, nprow )

        lld = MAX(1,Locp)
        call descinit(descA,mblk_size,mblk_size,mb,nb,rsrc,csrc,icontxt,lld,info)
        call assert(info.eq.0,'descinit descA')
        ineed = MAX(1,lld*Locq)
        IF (ALLOCATED(ublkp)) DEALLOCATE(ublkp, dblkp, lblkp, stat=istat)
        IF (.NOT. ALLOCATED(ublkp)) THEN
           ALLOCATE(ublkp(ineed,ns-1),                                  &
                    dblkp(ineed,ns),                                    &
                    lblkp(ineed,ns-1),                                  &
                    stat=istat)
           CALL ASSERT(istat.EQ.0,'Not enough memory to store a single block!')
        END IF

        ALLOCATE(ublk(mblk_size,mblk_size2,ns),                         &
                 dblk(mblk_size,mblk_size2,ns),                         &
                 lblk(mblk_size,mblk_size2,ns),                         &
                 stat=istat)
        CALL ASSERT(istat.EQ.0,'Not enough memory to store a single block!')
#else
        CALL ASSERT(.FALSE.,'LSCALAPACK=T BUT NOT MPI!')
#endif
      ELSE   !.NOT.LSCALAPACK

        mblk_size2 = mblk_size
        ALLOCATE(ublk(mblk_size,mblk_size,ns),                          &
                 dblk(mblk_size,mblk_size,ns),                          &
                 lblk(mblk_size,mblk_size,ns),                          &
                 stat=istat)
      
      END IF       !END LSCALAPACK

      IF (ALLOCATED(ublk)) THEN
         ublk = 0; dblk = 0; lblk = 0
      END IF

      ALLOCATE (DataItem(mblk_size), stat=istat)
      CALL ASSERT(istat.EQ.0,'DataItem allocation failed!')
      DataItem = 0
      xc = 0

#if defined(MPI_OPT)
      icolmpi = 0
#endif
      icol=0

      CALL second0(skston)
      CALL func
      CALL second0(skstoff)
      hessian_funct_island_time=hessian_funct_island_time+(skstoff-skston)

#if defined(MPI_OPT)
      CALL second0(colstarttime) !Added by SKS for timing comparisons
#endif
      VAR_TYPE: DO ntype = 1, ndims
         VAR_N: DO n = -ntor, ntor
            VAR_M: DO m = 0, mpol
                  
               icol=icol+1
#if defined(MPI_OPT)
               IF(MOD(icol-1,nprocs)==iam .OR. .NOT.LSCALAPACK) THEN
                  IF (LSCALAPACK) THEN
                     icolmpi = icolmpi+1
                  ELSE 
                     icolmpi = icol
                  END IF
#else
                  icolmpi = icol
#endif
                  MESH_3PT: DO mesh = 1, 3

                  IF (.NOT.(m.EQ.0 .AND. n.LT.0)) THEN
                  DO js = jstart(mesh), ns, 3
                     xc(icol,js) = eps
                  END DO
				  END IF

                  INHESSIAN=.TRUE.
                  CALL second0(skston)
                  CALL func
                  CALL second0(skstoff)
                  hessian_funct_island_time=hessian_funct_island_time+(skstoff-skston)

                  INHESSIAN=.FALSE.
#if defined(MPI_OPT)
#endif

!OFF FOR l_linearize=T                  gc = gc-gc1

!              COMPUTE PRECONDITIONER (HESSIAN) ELEMENTS. LINEARIZED EQUATIONS
!              OF TRI-DIAGONAL (IN S) FORM 
!
!              F(j-1) = a(j-1)x(j-2) + d(j-1)x(j-1) + b(j-1)x(j)
!              F(j)   =                a(j)x(j-1)   + d(j)  x(j)  + b(j)  x(j+1)
!              F(j+1) =                               a(j+1)x(j)  + d(j+1)x(j+1) + b(j+1)x(j+2)
!
!              HESSIAN IS H(k,j) == dF(k)/dx(j); aj == lblk; dj == dblk; bj = ublk

                  SKIP3_MESH: DO js = jstart(mesh), ns, 3
                  
                     xc(icol,js) = 0

!FORCE RESPONSE AT mp,np,nptype TO m,n,js,ntype VELOCITY PERTURBATION
 
                     !ublk(js-1)
                     js1 = js-1
                     IF (js1 .GT. 0) THEN
                        DataItem = gc(:,js1)/eps

                        IF (l_Diagonal_Only) THEN
                           temp=DataItem(icol)
                           DataItem=0
                           DataItem(icol)=temp
                        END IF
                        ublk(:,icolmpi,js1) = DataItem
                     END IF

                     !dblk(js)
                     DataItem = gc(:,js)/eps
                     IF (ALL(DataItem .EQ. zero)) DataItem(icol) = one  

                     IF (l_Diagonal_Only) THEN
                        temp=DataItem(icol)
                        DataItem=0
                        DataItem(icol)=temp
                     END IF

                     dblk(:,icolmpi,js) = DataItem

                     !lblk(js+1) 
                     js1=js+1
                     IF (js .LT. ns) THEN
                        DataItem = gc(:,js1)/eps
                        IF (l_Diagonal_Only) THEN
                           temp=DataItem(icol)
                           DataItem=0
                           DataItem(icol)=temp
                        ENDIF
                        lblk(:,icolmpi,js1) = DataItem
                     END IF

                  END DO SKIP3_MESH
               END DO MESH_3PT
#if defined(MPI_OPT)
               ENDIF
#endif
            END DO VAR_M
         END DO VAR_N
      END DO VAR_TYPE
      
      CALL second0(toff)
      
      IF (l_Diagonal_Only) THEN
         time_diag_prec = time_diag_prec+(toff-ton)
      ELSE
         time_block_prec = time_block_prec+(toff-ton)
      END IF

!SPH 071416 ADDED COLUMN SCALING OPTION
      IF (lColScale) THEN
         CALL SerialScaling
      ELSE
         CALL FindMinMax_TRI_Serial
      END IF
      
      CAll second0(ton)
      colendtime = ton
      colusedtime=colendtime-colstarttime
      time_generate_blocks=time_generate_blocks+colusedtime
      construct_hessian_time=construct_hessian_time+(colendtime-colstarttime)

      DEALLOCATE (DataItem, stat=istat)
      CALL ASSERT(istat.EQ.0,'DataItem deallocation error!')

      IF (nprocs .EQ. 1) THEN

!
!SYMMETRIZE BLOCKS: REWRITE FOR nprocs>1 USING SCALAPACK (only works now for mblk_size2=mblk_size)
!
      IF (.NOT.l_Hess_sym .OR. nprocs.GT.1) GOTO 1000
      ALLOCATE(gc1(mblk_size,ns), dblk_s(mblk_size, mblk_size2,ns))

      DO icol = 1, mblk_size2      
         gc1(:,1:nsh) = (ublk(icol,:,1:nsh) + lblk(:,icol,2:ns))/2
         ublk(icol,:,1:nsh) = gc1(:,1:nsh)
         lblk(:,icol,2:ns)  = gc1(:,1:nsh)
         dblk_s(:,icol,:) = (dblk(icol,:,:) + dblk(:,icol,:))/2
      END DO
      
      dblk = dblk_s
 
      DEALLOCATE (gc1, dblk_s)
 1000 CONTINUE

!
!CHECK CONDITION NUMBER
!
         CALL second0(ton)
         CALL CheckEigenvalues_Serial(ns, mblk_size)
         CALL CheckConditionNumber_Serial(ns,mblk_size,anorm,rcond,info)
         CALL second0(toff)                  
         IF (INFO .EQ. 0 .and. lverbose) THEN
            PRINT '(1x,3(a,1p,e12.3))','RCOND = ', rcond,               &
            ' ||A|| = ', ANORM,' TIME: ', toff-ton
         END IF
      END IF

!
!CHECK SYMMETRY OF BLOCKS 
!
      CALL second0(skston)
      IF (LSCALAPACK) THEN
         CALL check3d_symmetry(dblk,lblk,ublk,dblkp,lblkp,ublkp,        &
                               mblk_size,mblk_size2,ns,asym_index)
      ELSE
         CALL check3d_symmetry_serial(dblk,lblk,ublk,mpol,ntor,         &
                               mblk_size,ns,asym_index)
      END IF
      CALL second0(skstoff)
      asymmetry_check_time=asymmetry_check_time+(skstoff-skston)

 900  CONTINUE

!      IF (l_Diagonal_Only) RETURN


!      l_backslv = (nprocs.EQ.1 .AND. nprecon.GE.2)                        !Controls when dump is made
      IF (l_backslv) THEN
         istat = 0
         IF (.NOT.ALLOCATED(dblk_s))                                    &
         ALLOCATE(ublk_s(mblk_size,mblk_size,ns),                       &
                  dblk_s(mblk_size,mblk_size,ns),                       &
                  lblk_s(mblk_size,mblk_size,ns),                       &   
                  stat=istat)
         CALL ASSERT(istat.EQ.0,'Allocation error2 in Compute_Hessian_Blocks')
         ublk_s = ublk;  dblk_s = dblk; lblk_s = lblk
      END IF
!
!     FACTORIZE (Invert) HESSIAN
!
      IF (ALLOCATED(ipiv_blk)) DEALLOCATE(ipiv_blk,stat=istat)
      CALL second0(skston)
      IF (LSCALAPACK) THEN
#if defined(MPI_OPT)
         ineed = numroc(descA(M_),descA(MB_),0,0,nprow) + descA(MB_)
         ALLOCATE( ipiv_blk(ineed, ns), stat=istat )
         CALL ASSERT(istat.EQ.0,'ipiv_blk allocation error2 in Compute_Hessian_Blocks')
         CALL blk3d_parallel(dblk,lblk,ublk,dblkp,lblkp,ublkp,mblk_size,mblk_size2,ns)
         CALL blacs_gridinfo(icontxt,nprow,npcol,myrow,mycol)
         mm = mblk_size*ns
         LocqR = numroc( nrhs1, nb, mycol,csrc,npcol)
         LocpR = numroc( mm, mb, myrow,rsrc,nprow)
         lld = MAX(1,LocpR)
         CALL descinit(descX,mm,nrhs1, mb,nb,rsrc,csrc,icontxt,lld,info)
         CALL assert(info.eq.0,'test_pdtrd:descinit return info!=0')
         descR(:) = descX

         ineedR = MAX(1,lld*LocqR)
 
         CALL blacs_barrier(icontxt,'All')
         CALL profstart('factor call:123')
         CALL pdtrdf(mblk_size,ns,dblkp,lblkp,ublkp,ipiv_blk,descA)
         CALL blacs_barrier(icontxt,'All')
         CALL profend('factor call:123')
#else
         CALL ASSERT(.FALSE.,'LSCALAPACK=T BUT NOT IN MPI!')
#endif
      ELSE    !.NOT. LSCALAPACK
         ALLOCATE (ipiv_blk(mblk_size,ns), stat=istat)
         CALL ASSERT(istat.EQ.0, 'ipiv_blk allocation error2 in Compute_Hessian_Blocks')
         CALL blk3d_factor(dblk, lblk, ublk, ipiv_blk, mblk_size, ns)
      END IF  ! LSCALAPACK

      CALL second0(skstoff)
      block_factorization_time=block_factorization_time+(skstoff-skston)
      CALL second0(toff)
      time_factor_blocks = time_factor_blocks + (toff-ton)

      IF (iam.EQ.0 .AND. lverbose)                                      &
         PRINT '(a,1p,e12.3)',' BLOCK FACTORIZATION TIME: ', toff-ton

      END SUBROUTINE Compute_Hessian_Blocks_Thomas

      
      SUBROUTINE SerialScaling
      USE shared_data, ONLY:    lequi1, niter
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!----------------------------------------------
      REAL(dp), PARAMETER   :: zero=0, one=1
      INTEGER               :: js, icol, icol2, icount
      REAL(dp)              :: eps, minScale, ton, toff, colcnd
      REAL(dp), ALLOCATABLE :: col(:)
!----------------------------------------------
      CALL ASSERT(ALL(col_scale.EQ.1), 'COL_SCALE != 1 INITIALLY IN SERIAL SCALING')
      icount = 0

      CALL second0(ton)

      IF (lequi1) THEN
         minScale = 1.E-5_dp 
      ELSE
         minScale = 1.E-8_dp
      END IF

      ALLOCATE(colsum(mblk_size, ns), col(mblk_size))

 111  CONTINUE
      colsum = 0
      BLOCK_ROW1: DO js = 1, ns
         COLS1: DO icol = iam + 1, mblk_size, nprocs
            icol2 = CEILING((1.0*icol)/nprocs)
            col = ABS(dblk(:,icol2,js))
            IF (lequi1) THEN
               colsum(icol,js) = SUM(col)
            ELSE
               colsum(icol,js) = MAXVAL(col)
            END IF
            IF (js .GT. 1) THEN
               col = ABS(ublk(:,icol2,js-1))
               IF (lequi1) THEN
                  colsum(icol,js) = colsum(icol,js) + SUM(col)
               ELSE
                  colsum(icol,js) = MAX(colsum(icol,js),MAXVAL(col))
               END IF
            END IF
            IF (js .LT. ns) THEN
               col = ABS(lblk(:,icol2,js+1))
               IF (lequi1) THEN
                  colsum(icol,js) = colsum(icol,js) + SUM(col)
               ELSE
                  colsum(icol,js) = MAX(colsum(icol,js),MAXVAL(col))
               END IF
            END IF
         END DO COLS1
         IF (LSCALAPACK) THEN
            CALL GatherCols(colsum(:,js))
         END IF
         eps = minScale*MAXVAL(colsum(:,js))                                   !need ALL(DataItem.eq.zero) check
         CALL assert(eps.GT.zero,' coltmp == 0 in SerialScaling')
         colsum(:,js) = one/SQRT(MAX(colsum(:,js), eps))
      ENDDO BLOCK_ROW1

      CALL VECTOR_COPY_SER (colsum, col_scale)

!SCALE BLOCKS
      BLOCK_ROW2: DO js = 1, ns
         COLS2: DO icol = iam + 1, mblk_size, nprocs
            icol2 = CEILING((1.0*icol)/nprocs)
            dblk(:,icol2,js) = colsum(:,js)*dblk(:,icol2,js)*colsum(icol,js)
            IF (js .GT. 1) THEN
               lblk(:,icol2,js) = colsum(:,js)*lblk(:,icol2,js)*colsum(icol,js-1)
            END IF
            IF (js .LT. ns) THEN
               ublk(:,icol2,js) = colsum(:,js)*ublk(:,icol2,js)*colsum(icol,js+1)
            END IF
         END DO COLS2
      END DO BLOCK_ROW2

!MAXIMUM COLUMN VARIATIONS (like an inverse "Condition number")
      colcnd = MAXVAL(colsum)
      IF (colcnd .NE. ZERO) THEN
         colcnd = ABS(MINVAL(colsum)/colcnd)
      ELSE
         colcnd = 1
      END IF         

      icount = icount + 1

      DEALLOCATE(colsum, col)

      CALL FindMinMax_TRI_Serial
      CALL second0(toff)
      IF (iam.EQ.0 .AND. lverbose) PRINT '(1x,3(a,1p,e10.3))',          &
         'COLUMN-SCALING TIME: ', toff-ton, ' COLCND: ', colcnd,        &
         ' DMINL: ', DMIN_TRI*ABS(levmarq_param)


      END SUBROUTINE SerialScaling
      
!>  \brief subroutine for adding levmarq parameter to diagonal matrix elements
      SUBROUTINE FindMinMax_TRI_Serial
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!----------------------------------------------
      REAL(dp), PARAMETER  :: minThreshold=1.E-6, zero=0
      INTEGER              :: js, icol, jrow
      REAL(dp)             :: eps, DminL, temp(2), sign_diag=-1
!----------------------------------------------
      maxeigen = 0
      
      Dmin_TRI = HUGE(Dmin_TRI)

      BLOCK_ROW3: DO js = 1, ns
         mindiag = 0
         jrow = 1+iam
         COLS1: DO icol = 1, mblk_size2
            eps = dblk(jrow,icol,js)
            IF (ABS(eps) .GT. minDiag) THEN
               minDiag = ABS(eps)
               sign_diag = eps/minDiag
            END IF 
            maxeigen = MAX(maxeigen, ABS(lblk(jrow,icol,js))            &
                     + ABS(dblk(jrow,icol,js)) + ABS(ublk(jrow,icol,js)))
      
            jrow = jrow + nprocs
         END DO COLS1

#if defined(MPI_OPT)      
         IF (LSCALAPACK) THEN
            temp(1) = minDiag; temp(2) = maxeigen
            CALL MPI_ALLREDUCE(MPI_IN_PLACE,temp,1,MPI_REAL8, MPI_MAX,         &
                               SIESTA_COMM, MPI_ERR)
            minDiag = temp(1);  maxeigen = temp(2)
         END IF
#endif      
         minDiag = minThreshold*minDiag
         IF (minDiag .EQ. zero) minDiag = minThreshold
         
         jrow = 1+iam
      COLS2: DO icol = 1, mblk_size2
         eps = dblk(jrow,icol,js) 
	     IF (ABS(eps) .LE. mindiag) THEN
	        dblk(jrow,icol,js) = sign_diag*mindiag
         ELSE IF (js .LT. ns) THEN
	        Dmin_TRI = MIN(MAX(1.E-12_dp,ABS(eps)), Dmin_TRI)
	     END IF
         jrow = jrow + nprocs
      END DO COLS2
      END DO BLOCK_ROW3

#if defined(MPI_OPT)      
      IF (LSCALAPACK) THEN
         Dmin_Tri = ABS(Dmin_Tri)
         CALL MPI_ALLREDUCE(MPI_IN_PLACE,Dmin_Tri,1,MPI_REAL8, MPI_MIN,        &
                            SIESTA_COMM, MPI_ERR)
      END IF
#endif
      minDiag = sign_diag*Dmin_TRI
      DminL = ABS(levmarq_param)*minDiag
      
!SPH061016: RESTORE LM PARAMETER
      DO js = 1, ns
         jrow = 1+iam
	  COLS3: DO icol = 1, mblk_size2
          eps = dblk(jrow,icol,js)
          CALL ASSERT(sign_diag*eps.GT.zero,                            &
                     ' EPS*SIGN_DIAG < 0 IN FindMinMax_TRI_Serial')
          IF (lColScale) THEN
             dblk(jrow,icol,js) = eps + SIGN(DminL,eps)
          ELSE
             dblk(jrow,icol,js) = eps*(1+ABS(levmarq_param))
          END IF
          jrow = jrow + nprocs
	  END DO COLS3
      END DO

      CALL second0(toff)
      
      END SUBROUTINE FindMinMax_TRI_Serial
      

      SUBROUTINE VECTOR_COPY_SER(colsum, colscale)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      REAL(dp), INTENT(IN)  :: COLSUM(mblk_size,ns)
      REAL(dp), INTENT(OUT) :: COLSCALE(mblk_size,ns)
!-----------------------------------------------
      
      COLSCALE = COLSUM * COLSCALE

      END SUBROUTINE VECTOR_COPY_SER
      
      
      SUBROUTINE CheckEigenvalues_Serial(nblock, bsize)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER, INTENT(IN)    :: nblock, bsize
!-----------------------------------------------
      INTEGER, SAVE          :: nCount=0
      CHARACTER*1, PARAMETER :: JOBZ='N', UPLO='U'
      INTEGER                :: KD, N, LWORK, LDZ, OFFK, OFFS, OFFSU,   &
                                ROWMIN, ROWMAX, I, J, ICOL, JS, LDAB
      REAL(dp), ALLOCATABLE, DIMENSION(:)  :: W, ACOL, WORK
      REAL(dp), ALLOCATABLE, DIMENSION(:,:) :: AB, Z, BAND
      REAL(dp)               :: ton, toff
!-----------------------------------------------
      nCount = nCount+1
      
      IF (nCount .NE. 17) RETURN
      CALL second0(ton)
      
      KD = 2*bsize-1
      N = nblock*bsize; LDAB = KD+1; LDZ = 1
      LWORK = MAX(1,3*N-2)
      offk = kd+1
      offs = 0
      offsu = bsize
      
      ALLOCATE (W(N), ACOL(N), WORK(LWORK), AB(LDAB,N), stat=I)
      CALL ASSERT(I.EQ.0,'Allocation error in CheckEigenvalues_Serial')
      IF (JOBZ .NE. 'N') ALLOCATE(Z(LDZ,N))
      AB = 0    !stores (upper) bands (A assumed symmetric)
      
!compute each col (labelled "j") for all rows ACOL(1:n1)
      DO JS = 1, nblock
         DO ICOL = 1, bsize
            j = offs + icol
            ROWMIN = offs-bsize+1; ROWMAX = offsu+bsize
            IF (JS .GT. 1)                                              &
            ACOL(rowmin:offs) = ublk(:,icol,js-1)  
            IF (JS .LT. nblock)                                         &
            ACOL(offsu+1:rowmax) = lblk(:,icol,js+1) 
            ACOL(offs+1:offsu) = dblk(:,icol,js) 
          
            DO i = MAX(1,j-kd,rowmin), MIN(j,n,rowmax)
               AB(offk+i-j,j) = ACOL(i)
            END DO
            
         END DO

         offs = offsu
         offsu = offsu + bsize
         
      END DO

!real SYMMETRIC band matrix
      CALL DSBEV( JOBZ, UPLO, N, KD, AB, LDAB, W, Z, LDZ, WORK, INFO )

      DEALLOCATE (AB, ACOL)
      
      DO J = 1, N
         WRITE (4000+ncount, '(i5,1p,2e12.4)') J, W(J), W(J)/W(1)
      END DO
      CALL FLUSH(4000+ncount)

      DEALLOCATE (W)
      IF (JOBZ .NE. 'N') DEALLOCATE(Z)
      
      CALL second0(toff)
      IF (lverbose) THEN
         PRINT 100,' Eigenvalue TIME: ', (toff-ton), ' INFO: ', info,          &
                   ' Eigenvalues written to FORT.',4000+ncount
      END IF

 100  FORMAT(a,1p,e10.2,2(a,i4))
      
      END SUBROUTINE CheckEigenvalues_Serial
      

      SUBROUTINE CheckConditionNumber_Serial(nblock,bsize,anorm,rcond,info)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER, INTENT(IN)   :: nblock, bsize
      INTEGER, INTENT(OUT)  :: info
      REAL(dp), INTENT(OUT) :: anorm, rcond
!-----------------------------------------------
      REAL(dp), PARAMETER   :: zero=0, one=1
      CHARACTER(LEN=1),PARAMETER :: CHNORM='1'                          !1 (1-norm) or 'I' for infinity norm
      INTEGER :: N1, JS, ICOL, ISTAT, OFFS, OFFSU, I, J, KU, KL,        &
                 LDAB, OFFK, ROWMIN, ROWMAX
      INTEGER, ALLOCATABLE     :: IWORK(:), IPIV(:)
      REAL(dp), ALLOCATABLE    :: WORK(:), ACOL(:), AB(:,:)
      REAL(dp)                 :: ROWCND, COLCND, AMAX
!-----------------------------------------------
      N1 = bsize*nblock

      ALLOCATE (ACOL(N1), stat=istat)     
      CALL ASSERT(ISTAT.EQ.0,'COND # ALLOCATION1 FAILED IN HESSIAN')

      ku = 2*bsize-1
      kl = 2*bsize-1
      ldab = 2*kl+ku+1
      offk = kl+ku+1
  
      ALLOCATE(AB(ldab, n1), ipiv(n1), iwork(n1), work(3*n1),         &
               stat=istat)
      CALL ASSERT(ISTAT.EQ.0,'COND # ALLOCATION2 FAILED IN HESSIAN')

      offs = 0
      offsu = bsize
      ANORM = 0
      AB = 0    !stores bands

!compute each col (labelled "j") for all rows ACOL(1:n1)
      DO JS = 1, nblock
         DO ICOL = 1, bsize
            j = offs + icol
            ROWMIN = offs-bsize+1; ROWMAX = offsu+bsize
            IF (JS .GT. 1)                                              &
            ACOL(rowmin:offs) = ublk(:,icol,js-1)  
            IF (JS .LT. nblock)                                         &
            ACOL(offsu+1:rowmax) = lblk(:,icol,js+1) 
            ACOL(offs+1:offsu) = dblk(:,icol,js) 

            IF (JS.EQ.1) THEN
            ROWMIN = offs+1
            ELSE IF (JS.EQ.nblock) THEN
            ROWMAX = offsu
            END IF
            ANORM = MAX(ANORM, SUM(ABS(ACOL(rowmin:rowmax))))      
           
            DO i = MAX(1,j-ku,rowmin), MIN(n1,j+kl,rowmax)
               AB(offk+i-j,j) = ACOL(i)
            END DO
            
         END DO

         offs = offsu
         offsu = offsu + bsize
         
      END DO
      
      DEALLOCATE (ACOL)

      CALL DGBTRF( N1, N1, KL, KU, AB, LDAB, IPIV, INFO )
      IF (INFO.EQ.0)                                                    &
      CALL DGBCON( '1', N1, KL, KU, AB, LDAB, IPIV, ANORM, RCOND,       &
                  WORK, IWORK, INFO )

      DEALLOCATE(AB, work, iwork, ipiv)

      IF (info.EQ.0 .AND. rcond.NE.zero) THEN
         RCOND = one/RCOND
      ELSE
         RCOND = -one
      END IF

      END SUBROUTINE CheckConditionNumber_Serial

!#define SVD_ON
!#undef SVD_ON

      SUBROUTINE blk3d_factor(a, bm1, bp1, ipiv, mblk, nblocks)
      USE stel_constants, ONLY: one, zero
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER, INTENT(in) :: nblocks, mblk
      INTEGER, TARGET, INTENT(OUT) :: ipiv(mblk,nblocks)
      REAL(dp), TARGET, DIMENSION(mblk,mblk,nblocks),                   &
                     INTENT(INOUT) :: a, bm1, bp1
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
!      INTEGER :: ibuph, incnow, irecl, incbu, iunit=102, ndisk
      INTEGER :: k, k1, ier
      INTEGER, POINTER :: ipivot(:)
      REAL(dp), POINTER :: amat(:,:), bmat(:,:), cmat(:,:)
      REAL(dp), ALLOCATABLE, DIMENSION(:,:) :: temp
#if defined(SVD_ON)
      REAL(dp), ALLOCATABLE, DIMENSION(:)   :: w, work
      REAL(dp), ALLOCATABLE, DIMENSION(:,:) :: vt, u
      REAL(dp), PARAMETER :: small = 1.E-300_dp
      INTEGER :: nw
      INTEGER :: lwork
      CHARACTER(LEN=1), PARAMETER :: jobu='A', jobvt='A'
#endif
!-----------------------------------------------
!  modified (June, 2003, ORNL):         S. P. Hirshman
!-----------------------------------------------------------------------
!
!  this subroutine solves for the Q factors of a block-tridiagonal system of equations.
!
!-----------------------------------------------------------------------
!  INPUT
!  mblk                : block size
!  nblocks             : number of blocks
!  a                   : diagonal blocks
!  bp1, bm1            : lower, upper blocks (see equation below)
!
!  OUTPUT
!  ipiv                : pivot elements for kth block
!  a                   : a-1 LU factor blocks
!  bm1                 : q = a-1 * bm1 matrix
!
!  LOCAL VARIABLES
!  iunit               : unit number for block-tridiagonal solution disk file.
!
!  solutions are indexed in m-n fourier-space, legendre-space. the tri-diagonal
!  equation is:
!
!           bm1 * f(l-1) + a * f(l) + bp1 * f(l+1) = source(l)
!
!     GENERAL SOLUTION SCHEME APPLIED TO EACH BLOCK ROW (INDEX L)
!
!     1. Start from row N and solve for x(N) in terms of x(N-1):
!
!        x(N) = -q(N)*x(N-1) + r(N)
!
!        q(N) =  a(N)[-1] * bm1;    r(N) = a(N)[-1] * s(N)
!
!        where a(N)[-1] is the inverse of a(N)
!
!     2. Substitute for lth row to get recursion equation for q(l) and r(l):
!
!        x(l) = -q(l)*x(l-1) + r(l), in general, where:
!
!        q(l) = (a(l) - bp1(l)*q(l+1))[-1] * bm1(l)
!
!        qblk(l) == (a(l) - bp1(l) * q(l+1))[-1] on return
!
!        r(l) = (a(l) - bp1(l)*q(l+1))[-1] * (s(l) - bp1(l)*r(l+1))
!
!     3. At row l = 1, bm1(1) = 0 and get an equation for x(1) corresponding to q(1) = 0:
!
!        x(1) = r(1)
!
!     4. Finally, can back-solve for x(N) in terms of x(N-1) from eqn.(1) above
!
!
!     NUMERICAL IMPLEMENTATION (USING LAPACK ROUTINES)
!
!     IF USING SVD:
!
!     1. CALL dgesvd:   Perform SVD decomposition
!     2. CALL svbksb:   Solve Ax = bi, for i=1,mblk
!
!     OTHERWISE
!
!     1. CALL dgetrf:   Perform LU factorization of diagonal block (A) - faster than sgefa
!     2. CALL dgetrs:   With multiple (mblk) right-hand sides, to do block inversion
!                         operation, Ax = B  (stores result in B; here B is a matrix)
!
!      ndisk = mblk

!  main loop. load and process (backwards) block-rows nblocks to 1. 

#if defined(SVD_ON)
      lwork = 10*mblk
      ALLOCATE(vt(mblk,mblk), w(mblk),                                  &
               u(mblk,mblk), work(lwork), stat=ier)
      CALL ASSERT(ier.EQ.0,'Allocation error in blk3d_factor!')
#endif
      ipiv = 0
      ALLOCATE (temp(mblk,mblk), stat=ier)
      CALL ASSERT(ier.EQ.0, 'Allocation error in blk3d_factor!')

      BLOCKS: DO k = nblocks, 1, -1
!
!     Compute (and save) qblk(nblocks) = ablk(nblocks)[-1] * bml
!
         amat => a(:,:,k)

#if defined(SVD_ON)
!NOTE: v = vt coming out here
         CALL dgesvd (jobu, jobvt, mblk, mblk, amat, mblk, w, u, mblk,  &
                      vt, mblk, work, lwork, ier)
!Set SVD weights w to 0 for weights below the allowed threshold, so backsolver
!will compute pseudo-inverse
         WRITE (35, '(a,i4,a,1pe12.4)') 'Block: ',k, ' Condition #: ', SQRT(w(1)/w(mblk))
         PRINT '(a,i4,a,1pe12.4)',' BLOCK: ',k,' SVD COND #: ', SQRT(w(1)/w(mblk))
         DO nw = 2, mblk
            IF (w(nw) .LE. small*w(1)) w(nw) = 0
         END DO
	     CALL ASSERT(ier.EQ.0,'SVD ERROR')
!
!        STORE svd pseudo-inverse IN AMAT since u,vt,w will be deallocated at end
         CALL svdinv2 (amat, u, vt, w, mblk)
#else
         ipivot => ipiv(:,k)
         CALL dgetrf (mblk, mblk, amat, mblk, ipivot, ier)
	     CALL ASSERT(ier.EQ.0,'DGETRF ERROR')
#endif
         IF (k .eq. 1) EXIT

         bmat => bm1(:,:,k)

#if defined(SVD_ON)
!
!        Use (pseudo) inverse stored in AMAT = V*1/w*Ut 
!
         temp = bmat
         bmat = MATMUL(amat, temp)
#else
         CALL dgetrs('n', mblk, mblk, amat, mblk, ipivot,               &
                     bmat, mblk, ier)
         CALL ASSERT(ier.EQ.0,'dgetrs INFO != 0')
#endif
!         CALL wrdisk(iunit, ql, ndisk, incnow, ibuph, incbu, ier)
!         IF (ier .NE. 0) GOTO 302

!
!      Update effective diagonal "a" matrix. Use dgemm: faster AND doesn't overflow normal stack
!
         k1 = k-1 
         amat => bp1(:,:,k1)
         cmat => a(:,:,k1)
!         cmat = cmat - MATMUL(amat, bmat)
         CALL dgemm('N','N',mblk,mblk,mblk,-one,amat,mblk,              &
                    bmat, mblk, one, cmat, mblk)

      END DO BLOCKS

#if defined(SVD_ON)
      DEALLOCATE(vt, w, u, work)
#endif
!
!     COMPUTE TRANSPOSES HERE, SINCE REPEATEDLY CALLING MATMUL OPERATION
!     X*At IS FASTER THAN A*X DUE TO UNIT STRIDE
!

      DO k = 1, nblocks
         IF (k .NE. nblocks) THEN
            temp = TRANSPOSE(bp1(:,:,k))
            bp1(:,:,k) = temp
         END IF
         IF (k .NE. 1) THEN
            temp = TRANSPOSE(bm1(:,:,k))
            bm1(:,:,k) = temp
         END IF
      END DO

      DEALLOCATE (temp)

      END SUBROUTINE blk3d_factor


      SUBROUTINE block_precond(gc)
      USE timer_mod
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      REAL(dp), DIMENSION(mblk_size,ns), INTENT(INOUT) :: gc
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      REAL(dp), PARAMETER  ::  zero = 0
      INTEGER :: m, n, js, ntype, istat, irow
      REAL(dp) :: t1, error, error_sum, b_sum
      REAL(dp), ALLOCATABLE :: gc_s(:,:)
#if defined(MPI_OPT)
      REAL(dp) :: alpha0,beta0
      INTEGER  :: i,j,k,nn,kk,ia,ja,ix,jx,ic,jc,ib,jb
#endif
      REAL(dp) :: starttime, endtime, usedtime
      INTEGER :: PRECONDPASS=0
      INTEGER :: globrow

!-----------------------------------------------
      starttime = 0
      endtime = 0

!-----------------------------------------------
!
!     Applies 3D block-preconditioner to forces vector (gc)
!
      IF (PARSOLVER) THEN
         DO globrow = myStartBlock, myEndBlock
            CALL SetMatrixRHS(globrow, gc(1:mblk_size,globrow))
         END DO
        
         CALL BackwardSolve

         DO globrow = myStartBlock, myEndBlock
            CALL GetSolutionVector(globrow, gc(:,globrow))
         END DO

         CALL gather_array(gc)
      ELSE

         CALL second0(starttime)

         IF (l_backslv) THEN
            istat = 0
            IF (.NOT.ALLOCATED(gc_s)) THEN
               ALLOCATE (gc_s(mblk_size,ns), stat=istat)
            END IF
            CALL ASSERT(istat.EQ.0,'Allocation error0 in block_precond')
            gc_s = gc
         END IF

!     Solve Hx = gc, using Hessian (H) LU factors stored in block_... matrices
!     Store solution "x" in gc on output
         IF (LSCALAPACK) THEN
#if defined(MPI_OPT)
            ALLOCATE (tempp(ineedR), stat=istat)
            ia = 1
            ja = 1
            ib = 1
            jb = 1
            beta0 = 0.0d0
            alpha0 = 1.0d0
            mm0 = descR(3)
            nrhs0 = descR(4)
            call profstart('pdgeadd call:123')
            call pdgeadd('N', mm0, nrhs0, alpha0, gc, ia, ja, descR_all,       &
                         beta0, tempp, ib, jb, descR)
            call profend('pdgeadd call:123')
            ir = 1
            jr = 1
            call blacs_barrier(icontxt,'All')
            call profstart('solver call:123')
            call pdtrds(mblk_size,ns,dblkp,lblkp,ublkp,ipiv_blk,descA,         &
                        nrhs1,tempp,ir,jr,descR)
            call blacs_barrier(icontxt,'All')
            call profend('solver call:123')
            gc = 0
            ia = 1
            ja = 1
            ib = 1
            jb = 1
            beta0 = 0.0d0
            alpha0 = 1.0d0
            mm0 = descR(3)
            nrhs0 = descR(4)
!     descR_all(7)=-1
!     descR_all(8)=-1
            call profstart('pdgeadd2 call:123')
            call pdgeadd('N', mm0,nrhs0, alpha0, tempp,ia,ja,descR,            &
                         beta0,gc,ib,jb,descR_all)
            call dgsum2d( icontxt,'A', ' ', mm,1,gc,mm,-1,-1)
            call profend('pdgeadd2 call:123')

            DEALLOCATE (tempp, stat=istat)
#else
            CALL ASSERT(.FALSE.,'MPI_OPT must be true for LSCALAPACK!')
#endif
         ELSE       !NOT LSCALAPACK

            CALL second0(endtime)
            usedtime=endtime-starttime
            IF (SKSDBG) WRITE(TOFU,*)                                          &
               'HessianTag-2 : Time to BackwardSolve:',usedtime,               &
               'PRECONDPASS',PRECONDPASS,'in native run'
            IF (SKSDBG) CALL FLUSH(TOFU)

!        Serial solver
            CALL blk3d_slv(dblk, lblk, ublk, gc, ipiv_blk, mblk_size, ns)
         END IF     !END LSCALAPACK
      END IF !END OF IF(PARSOLVER) CONDITIONAL
      PRECONDPASS = PRECONDPASS + 1

      IF (l_backslv .AND. ALLOCATED(dblk_s) .AND. .NOT.l_linearize) THEN
         error_sum = 0;  b_sum = SQRT(SUM(gc_s*gc_s)/neqs)
         CALL ASSERT(ALLOCATED(gc_s),'gc_s is not allocated')

         WRITE (34, '(2(/,a))') ' BLK3D FACTORIZATION CHECK: Ax = b',   &
              '  IROW     M     N NTYPE         FORCE       DELTA-FORCE'
        
         DO js = startglobrow, endglobrow
            irow = 0
			WRITE (34,*) ' JS: ', js
            DO ntype = 1, ndims
               DO n = -ntor, ntor
                  DO m = 0, mpol
                     irow = irow+1
                     t1 = SUM(dblk_s(irow,:,js)*gc(:,js))
				     IF (js .LT. ns) THEN
                        t1 = t1 + SUM(ublk_s(irow,:,js)*gc(:,js+1))
                     END IF
			         IF (js .GT. 1) THEN
                        t1 = t1 + SUM(lblk_s(irow,:,js)*gc(:,js-1))
                     END IF

                     error = t1 - gc_s(irow,js)
                     error_sum = error_sum + error**2
                     WRITE(34,110) irow, m, n, ntype, gc_s(irow,js)/b_sum,     &
                                   error/b_sum
                  END DO
               END DO
            END DO
         END DO

         DEALLOCATE(dblk_s, lblk_s, ublk_s, gc_s, stat=istat)
         IF (lverbose) THEN
            PRINT '(3(a,1pe12.4))',' |b|: ', b_sum,                            &
                  ' |x| = ', SQRT(SUM(gc*gc)/neqs),                            &
                  ' Hessian Error: ', SQRT(error_sum/neqs)/b_sum
         END IF

      END IF
 110  FORMAT(4i6,1p,2e16.4)

 100  FORMAT(i6,1p,5e14.4)
 101  FORMAT(i6,1p,5e14.4,'  *')

      END SUBROUTINE block_precond

 
      SUBROUTINE apply_precond(gc)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      REAL(dp), DIMENSION(mblk_size,ns), INTENT(INOUT) :: gc
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      REAL(dp) :: ton, toff
!-----------------------------------------------
      CALL second0(ton)
!      IF (ltype) THEN
!         CALL diag_precond(gc)
!      ELSE
         CALL block_precond(gc)
!      END IF 
      CALL second0(toff)
      time_apply_precon = time_apply_precon+(toff-ton)

      END SUBROUTINE apply_precond


      SUBROUTINE blk3d_slv(ablk, qblk, bp1, source, ipiv, mblk, nblocks)
      USE stel_kinds
      USE stel_constants, ONLY: one
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER, INTENT(in) :: nblocks, mblk
      INTEGER, TARGET, INTENT(in) :: ipiv(mblk,nblocks)
      REAL(dp), TARGET, DIMENSION(mblk,mblk,nblocks), INTENT(IN) ::     &
                                     ablk, qblk, bp1
      REAL(dp), DIMENSION(mblk,nblocks), TARGET, INTENT(INOUT) :: source
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER, POINTER  :: ipivot(:)
!      INTEGER :: ibuph, incnow, irecl, incbu, iunit=102, ndisk
      INTEGER :: k, ier
      REAL(dp), POINTER :: amat(:,:), x1(:), source_ptr(:)
!-----------------------------------------------
!  modified (June, 2003, ORNL):         S. P. Hirshman
!-----------------------------------------------------------------------
!
!  this subroutine solves a block-tridiagonal system of equations, using 
!  the ABLK, QBLK factors from blk3d_factor,
!
!-----------------------------------------------------------------------
!  INPUT
!  mblk                : block size
!  nblocks             : number of blocks
!  bp1                 : upper blocks (see equation below)
!  ipiv                : pivot elements for kth block
!  ablk                : a-1 blocks
!  qblk                : q = a-1 * bm1
!  source              : input right side
!
!  OUTPUT
!  source              : Solution x of A x = source
! 
!  LOCAL VARIABLES
!  iunit               : unit number for block-tridiagonal solution disk file.
!
!  the tri-diagonal equation is:
!
!           bm1 * f(l-1) + a * f(l) + bp1 * f(l+1) = source(l)
!
!     GENERAL SOLUTION SCHEME APPLIED TO EACH BLOCK ROW (INDEX L)
!
!     1. Start from row N and solve for x(N) in terms of x(N-1):
!
!        x(N) = -q(N)*x(N-1) + r(N)
!
!        q(N) =  a(N)[-1] * bm1;    r(N) = a(N)[-1] * s(N)
!
!        where a(N)[-1] is the inverse of a(N)
!
!     2. Substitute for lth row to get recursion equation fo q(l) and r(l):
!
!        x(l) = -q(l)*x(l-1) + r(l), in general, where:
!
!        q(l) = (a(l) - bp1(l)*q(l+1))[-1] * bm1(l)
!
!        qblk(l) == (a(l) - bp1(l) * q(l+1))[-1] on return
!
!        r(l) = (a(l) - bp1(l)*q(l+1))[-1] * (s(l) - bp1(l)*r(l+1))
!
!     3. At row l = 1, bm1(1) = 0 and get an equation for x(1) corresponding to q(1) = 0:
!
!        x(1) = r(1)
!
!     4. Finally, can back-solve for x(N) in terms of x(N-1) from eqn.(1) above
!
!
!     NUMERICAL IMPLEMENTATION (USING LAPACK ROUTINES)
!
!     1. CALL dgetrs:   With single right hand side (source) to solve A x = b (b a vector)
!                       Faster than dgesl
!      ndisk = mblk

!  main loop. load and process (backwards) block-rows nblocks to 1. 
!  note: about equal time is spent in calling dgetrs and in performing
!  the two loop sums: on ibm-pc, 2 s (trs) vs 3 s (sums); on linux (logjam),
!  2.4 s (trs) vs 3 s (sums).

!
!     Back-solve for modified sources first
!
      BLOCKS: DO k = nblocks, 1, -1

         source_ptr => source(:,k)
         amat => ablk(:,:,k)
#if defined(SVD_ON)
         source_ptr = MATMUL(amat, source_ptr)
#else
         ipivot => ipiv(:,k);   
         CALL dgetrs('n', mblk, 1, amat, mblk, ipivot, source_ptr, mblk, ier)
		 CALL ASSERT(ier.EQ.0,'DGETRS INFO != 0')
#endif
         IF (k .eq. 1) EXIT
!
!        NOTE: IN BLK3D_FACTOR, BP1 AND BM1 WERE TRANSPOSED (AND STORED)
!        TO MAKE FIRST INDEX FASTEST VARYING IN THE FOLLOWING MATMUL OPS
!
         amat => bp1(:,:,k-1)
         x1 => source(:,k);  source_ptr => source(:,k-1)
         CALL dgemv('T',mblk,mblk,-one,amat,mblk,x1,1,one,source_ptr,1)

      END DO BLOCKS
!
!  forward (back-substitution) solution sweep for block-rows k = 2 to nblocks
!  now, source contains the solution vector
!
      DO k = 2, nblocks
!         CALL rddisk (iunit, ql, ndisk, incnow, ibuph, ier)
!         IF (ier .NE. 0) GOTO 303
!         ibuph = ibuph - incbu

         amat => qblk(:,:,k)
         x1 => source(:,k-1);  source_ptr => source(:,k)
!         source_ptr = source_ptr - MATMUL(x1,amat)  !USE THIS FORM IF TRANSPOSED qblk
!         source_ptr = source_ptr - MATMUL(amat,x1)  !UNTRANSPOSED FORM
         CALL dgemv('T',mblk,mblk,-one,amat,mblk,x1,1,one,source_ptr,1)

      END DO

      END SUBROUTINE blk3d_slv

      
      SUBROUTINE dealloc_hessian
      INTEGER :: istat

      IF (ALLOCATED(ublk)) DEALLOCATE(ublk, dblk, lblk, stat=istat)
      IF (ALLOCATED(ublkp))DEALLOCATE(dblkp,lblkp,ublkp,stat=istat)
      IF (ALLOCATED(mupar_norm)) DEALLOCATE(mupar_norm, stat=istat)
      IF (ALLOCATED(ipiv_blk)) DEALLOCATE (ipiv_blk)

      END SUBROUTINE dealloc_hessian

#if defined(MPI_OPT)
      SUBROUTINE blk3d_parallel(a,b,c,ap,bp,cp,mblk,mblkc,nblocks)
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      INTEGER :: ipart,jpart,k
      INTEGER :: ia,ja,inc1,ib,jb,inc2,nsm
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER, INTENT(IN)    :: nblocks, mblk, mblkc
        REAL(dp) :: aij2,aij
      REAL(dp), TARGET, DIMENSION(mblk,mblkc,nblocks),                  &
               INTENT(INOUT) ::  a, b, c
      REAL(dp), TARGET, DIMENSION(Locq,Locp,nblocks),                   &
               INTENT(INOUT) ::  ap
      REAL(dp), TARGET, DIMENSION(Locq,Locp,nblocks-1),                 &
               INTENT(INOUT) ::  bp, cp
      REAL(dp), POINTER :: amat(:,:), bmat(:,:), cmat(:,:)
      REAL(dp), POINTER :: amatp(:,:), bmatp(:,:), cmatp(:,:)

      do k=1,nblocks
        call profstart('change context:123')
        call blacs_barrier(descA(CTXT_),'All')
         amat => a(:,:,k)
         amatp => ap(:,:,k)
          call pdgemr2d(mblk,mblk,amat,1,1,descA_1xp,                   &
                       amatp,1,1,descA,  icontxt_global)
        if(k<nblocks)then
         bmat => b(:,:,k+1)
         bmatp => bp(:,:,k)
          call pdgemr2d(mblk,mblk,bmat,1,1,descA_1xp,                   &
                        bmatp,1,1,descA,  icontxt_global)
         cmat => c(:,:,k)
         cmatp => cp(:,:,k)
          call pdgemr2d(mblk,mblk,cmat,1,1,descA_1xp,                   &
                        cmatp,1,1,descA,  icontxt_global)
       endif
        call blacs_barrier(descA(CTXT_),'All')
        call profend('change context:123')
!        do ja=1,mblk
!          call pdelget( 'A',' ',aij, amatp,ja,ja,descA )
!          call pdelget( 'A',' ',aij2, amat,ja,ja,descA_1xp )
!              write(*,*) 'ja,ja,aij ',ja,ja,aij,aij2
!        enddo
       end do
      END SUBROUTINE blk3d_parallel
#endif

      END MODULE hessian
