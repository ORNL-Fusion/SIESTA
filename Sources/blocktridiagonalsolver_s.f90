!-------------------------------------------------------------------------------
!>  \brief Module contains subroutines for solver for block tri-diagonal matrices.
!!  \author  Kalyan S. Perumalla and Sudip K. Seal, ORNL 
!!  \date 2009-2017

!-------------------------------------------------------------------------------
MODULE blocktridiagonalsolver_s
USE v3_utilities, ONLY: assert
#if defined(MPI_OPT)
USE mpi_inc
#endif
USE descriptor_mod, ONLY: SIESTA_COMM
USE shared_data, ONLY: lequi1
IMPLICIT NONE
PRIVATE
!-------------------------------------------------------------------------------
! Precision settings
!-------------------------------------------------------------------------------
INTEGER, PARAMETER :: rprec = SELECTED_REAL_KIND(12,100)
INTEGER, PARAMETER :: iprec = SELECTED_INT_KIND(8)
INTEGER, PARAMETER :: cprec = KIND((1.0_rprec,1.0_rprec))
INTEGER, PARAMETER :: dp = rprec
INTEGER :: MPI_ERR

!-------------------------------------------------------------------------------
!>
!! Data associated with each row at each level
!<
!-------------------------------------------------------------------------------
TYPE LevelElement
  REAL(dp), ALLOCATABLE :: L(:,:), D(:,:), U(:,:), b(:,:)
  INTEGER, ALLOCATABLE :: pivot(:)
END TYPE LevelElement
!-------------------------------------------------------------------------------
!>
!! Solution of selected rows of interest to this rank
!<
!-------------------------------------------------------------------------------
TYPE SolutionElement
  REAL(dp), ALLOCATABLE :: x(:) !Vector of size M
END TYPE SolutionElement

!-------------------------------------------------------------------------------
!>
!! Data for each row at each level on this rank
!! The first dimension is the level number [1..L], L=#levels of forward solve
!! The 2nd dimension is the row number [1..K+g], K=#rows at level 1 on this rank
!! The +g is for incoming results from neighbors, 0<=g<=2
!<
!-------------------------------------------------------------------------------
TYPE (LevelElement), ALLOCATABLE :: lelement(:,:)

!-------------------------------------------------------------------------------
!>
!! Initial problem specification saved for verification at end of solution
!! The dimension is the row number [1..K], K=#rows at level 1 on this rank
!<
!-------------------------------------------------------------------------------
TYPE (LevelElement), ALLOCATABLE :: orig(:)

!-------------------------------------------------------------------------------
!>
!! The solution
!! The dimension is the global (level 1) row number [1..N]
!<
!-------------------------------------------------------------------------------
TYPE (SolutionElement), ALLOCATABLE :: selement(:)

!-------------------------------------------------------------------------------
INTEGER :: rank !<This MPI task's rank
INTEGER :: nranks !<Num of MPI tasks
INTEGER :: P !<Num of "master" tasks, adjusted such that P=Min(nranks,N)
INTEGER :: N !<Num of row blocks in input block tri-diag matrix
INTEGER :: M !<Size of each square sub-matrix block
INTEGER :: startglobrow !<Starting row number of this processor at first level
INTEGER :: endglobrow !<Ending row number (inclusive) at first level
INTEGER :: startrow1  !<1 row less than startglobrow
INTEGER :: endrow1    !<1 row greater than endglobrow
INTEGER :: nlevels !<Total number of levels in recursion
LOGICAL :: matdirtied !<Has a new matrix been set after ForwardSolve?
LOGICAL :: rhsdirtied !<Has a new RHS been set after ForwardSolve?

!-------------------------------------------------------------------------------
CHARACTER*100 :: kenvvar
CHARACTER*100 :: kenvval
LOGICAL :: LBCYLIC_INIT=.FALSE.
LOGICAL :: LPLB_INIT=.FALSE.
LOGICAL :: KPDBG !<Should debugging output be written? : SKS on Nov 9, 2011
INTEGER :: OFU !<Output file unit
INTEGER :: PFU !< Problem file unit
LOGICAL :: writeproblemfile !<Should save the random test case to a file?
LOGICAL :: writesolution !<Should dump solution in output?
LOGICAL :: usebarriers !<Should barriers be invoked between levels?
REAL :: membytes !<A running count of memory allocated so far
REAL :: dpsz !<Byte size of a double precision variable
REAL :: intsz !<Byte size of an integer variable
REAL :: ptrsz !<Byte size of a pointer variable
REAL(dp) :: ONE !<1.0
REAL(dp) :: ZERO !<0.0
LOGICAL(dp) :: L_COLSCALE=.FALSE.  !<set true if column scaling is applied

!-------------------------------------------------------------------------------
LOGICAL :: use_mpiwtime=.FALSE. !<Use MPI's timer function?
DOUBLE PRECISION :: loctimer1, loctimer2 !<Stopwatch snapshots to use locally
DOUBLE PRECISION :: mattimer1, mattimer2 !<Stopwatch snapshots for matrix ops
DOUBLE PRECISION :: globtimer1, globtimer2 !<Stopwatch snapshots to use globally
DOUBLE PRECISION :: timerfreq !<Timer frequency in Hz (count rate)
DOUBLE PRECISION :: tottime !<Total elapsed time
INTEGER :: totcount !<Total count
DOUBLE PRECISION :: totcommtime !<Total communication time
INTEGER :: totcommcount !<Total communication operation count
DOUBLE PRECISION :: totinvtime !<Total time to invert matrices
INTEGER :: totinvcount !<Total count of matrix inversions
DOUBLE PRECISION :: totmatmultime !<Total time for matrix-matrix multiplications
INTEGER :: totmatmulcount !<Total count of matrix-matrix multiplications
DOUBLE PRECISION :: totmatsoltime !<Total time for matrix-solve multiplications
INTEGER :: totmatsolcount !<Total count of matrix solutions

!-------------------------------------------------------------------------------
!>
!!  Unmodified diagonal element : SKS
!<
!-------------------------------------------------------------------------------
REAL(dp), ALLOCATABLE :: OrigDiagElement(:)
REAL(dp), ALLOCATABLE :: TopScaleFac(:), BotScaleFac(:)
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!>
!! BLACS/PBLAS options
!<
!-------------------------------------------------------------------------------
LOGICAL :: doblasonly
LOGICAL :: doblacscomm

!-------------------------------------------------------------------------------
!>
!! BLACS/PBLAS process grid information
!<
!-------------------------------------------------------------------------------
TYPE BlacsProcessGrid
  INTEGER :: myrow, mycol
  INTEGER :: nrows, ncols
  INTEGER :: blockszrows, blockszcols
  INTEGER, ALLOCATABLE :: map(:,:) !<Of size (1:nrows, 1:ncols)
END TYPE BlacsProcessGrid

!-------------------------------------------------------------------------------
!>
!! BLACS/PBLAS information
!<
!-------------------------------------------------------------------------------
TYPE BlacsParameters
  INTEGER :: iam
  INTEGER :: nprocs
  INTEGER :: maincontext
  INTEGER :: levelcontext
  TYPE(BlacsProcessGrid) :: pgrid
  INTEGER :: nbpp !<Min #PBLAS blocks per dim per proc (comp vs. comm tradeoff)
END TYPE BlacsParameters

TYPE(BlacsParameters) :: blacs

!-------------------------------------------------------------------------------
!>
!! Master-to-slave mapping
!<
!-------------------------------------------------------------------------------
TYPE MasterToSlaveMapping
  INTEGER :: masterrank
  INTEGER :: nslaves
  INTEGER, ALLOCATABLE :: slaveranks(:)
END TYPE MasterToSlaveMapping

!-------------------------------------------------------------------------------
!>
!! Level-specific PBLAS information
!<
!-------------------------------------------------------------------------------
TYPE PBLASLevelParameters
  LOGICAL :: ammaster
  TYPE(MasterToSlaveMapping) :: msmap

  INTEGER :: mpicomm, mpitag, mpierr
#if defined(MPI_OPT)
  INTEGER :: mpistatus(MPI_STATUS_SIZE)
#endif 

END TYPE PBLASLevelParameters

TYPE(PBLASLevelParameters) :: pblas

!-------------------------------------------------------------------------------
!>
!! Master-to-slave commands
!<
!-------------------------------------------------------------------------------
INTEGER, PARAMETER :: OP_NONE = 0
INTEGER, PARAMETER :: OP_DONE = 1
INTEGER, PARAMETER :: OP_DGEMM = 2
INTEGER, PARAMETER :: OP_DGEMV = 3
INTEGER, PARAMETER :: OP_DGETRF = 4
INTEGER, PARAMETER :: OP_DGETRS = 5

!-------------------------------------------------------------------------------
!>
!! Statistics (timing, etc.)
!<
!-------------------------------------------------------------------------------
TYPE TimeCount
  DOUBLE PRECISION :: tm
  INTEGER :: cnt
  DOUBLE PRECISION :: t1, t2 !<Stopwatch snapshots
END TYPE TimeCount
TYPE PBLASStats
  TYPE(TimeCount) :: wait, comm, comp
  TYPE(TimeCount) :: mm, trf, pmm, ptrf
  TYPE(TimeCount) :: mma, mmb, mmc, mmalpha, mmbeta, mmrc
  TYPE(TimeCount) :: extract, waitall
END TYPE PBLASStats

TYPE(PBLASStats) :: pstats

!-------------------------------------------------------------------------------
TYPE PBLASTempArray
  DOUBLE PRECISION, ALLOCATABLE :: temparray(:)
END TYPE PBLASTempArray

!-------------------------------------------------------------------------------
!PUBLIC METHODS
PUBLIC :: SetMatrixRowColL, SetMatrixRowColD, SetMatrixRowColU,         &
          GetMatrixRowColL, GetMatrixRowColD, GetMatrixRowColU,         &
          GetMatrixRHS, StoreDiagonal, Initialize, ForwardSolve,        &
          SetMatrixRHS, BackwardSolve, GetSolutionVector, PadSides,     &
          CheckSymmetry, ApplyParallelScaling, FindMinMax_TRI, GetRanks,&
          RefactorHessian, ParMatvec, Finalize, CheckConditionNumber
REAL(dp), PUBLIC :: MAXEIGEN_TRI  !<maximum eigenvalue from Gerschgorin theorem
REAL(dp), PUBLIC :: DMIN_TRI      !<minimum non-zero diagonal element, used to scale Levenberg parameter (SPH 08-31-16)

CONTAINS

!-------------------------------------------------------------------------------
!>
!! A convenience function to be able to change clock routine(s) easily
!<
!-------------------------------------------------------------------------------
SUBROUTINE BClockInit()
  !-----------------------------------------------------------------------------
  ! Local variables
  !-----------------------------------------------------------------------------
  INTEGER :: tempint

  IF ( use_mpiwtime ) THEN
      timerfreq = 1.0
  ELSE
      CALL SYSTEM_CLOCK(COUNT_RATE=tempint)
      timerfreq = tempint
  END IF
END SUBROUTINE BClockInit

!-------------------------------------------------------------------------------
!>
!! A convenience function to be able to change clock routine easily
!<
!-------------------------------------------------------------------------------
SUBROUTINE BSystemClock( ts )
  !-----------------------------------------------------------------------------
  ! Formal arguments
  !-----------------------------------------------------------------------------
  DOUBLE PRECISION, INTENT(INOUT) :: ts !<Time snapshot

  !-----------------------------------------------------------------------------
  ! Local variables
  !-----------------------------------------------------------------------------
  INTEGER :: tempint
  DOUBLE PRECISION :: tempdbl

  IF ( use_mpiwtime ) THEN
#if defined(MPI_OPT)
      tempdbl = MPI_WTIME()
      ts = tempdbl
#endif
  ELSE
      CALL SYSTEM_CLOCK( tempint )
      ts = tempint
  END IF
END SUBROUTINE BSystemClock

!-------------------------------------------------------------------------------
!>
!! A convenience function to be able to change or turn off easily
!<
!-------------------------------------------------------------------------------
SUBROUTINE FL( u )
  !-----------------------------------------------------------------------------
  ! Formal arguments
  !-----------------------------------------------------------------------------
  INTEGER, INTENT(IN) :: u !<Unit number to flush

  CALL FLUSH( u )
END SUBROUTINE FL

!-------------------------------------------------------------------------------
!>
!! Convenience routine to track allocated memory sizes
!<
!-------------------------------------------------------------------------------
SUBROUTINE ChargeMemory( bytes )
  !-----------------------------------------------------------------------------
  ! Formal arguments
  !-----------------------------------------------------------------------------
  REAL, INTENT(IN) :: bytes !<Number of bytes allocated

  !-----------------------------------------------------------------------------
  membytes = membytes + bytes
END SUBROUTINE ChargeMemory

!-------------------------------------------------------------------------------
!>
!! Convenience routine to accumulate timing values
!<
!-------------------------------------------------------------------------------
SUBROUTINE ChargeTime( tot, t2, t1, cnt )
  !-----------------------------------------------------------------------------
  ! Formal arguments
  !-----------------------------------------------------------------------------
  DOUBLE PRECISION, INTENT(INOUT) :: tot !<Total time
  DOUBLE PRECISION, INTENT(IN) :: t2 !<Timer value at end of segment being timed
  DOUBLE PRECISION, INTENT(IN) :: t1 !<Timer value at start of timed segment
  INTEGER, INTENT(INOUT) :: cnt !<Counter of number of samples

  !-----------------------------------------------------------------------------
  tot = tot + (REAL(t2-t1))/timerfreq
  cnt = cnt + 1
END SUBROUTINE ChargeTime

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!>
!! Initialize stats
!<
!-------------------------------------------------------------------------------
SUBROUTINE TimeCountInit( tc )
  !----------------------------------------------
  ! Formal arguments
  !----------------------------------------------
  TYPE(TimeCount), INTENT(INOUT) :: tc

  !----------------------------------------------
  tc%tm = 0
  tc%cnt = 0
END SUBROUTINE TimeCountInit

!-------------------------------------------------------------------------------
!>
!! Print stats
!<
!-------------------------------------------------------------------------------
SUBROUTINE TimeCountPrint( tc, msg )
  !----------------------------------------------
  ! Formal arguments
  !----------------------------------------------
  TYPE(TimeCount), INTENT(INOUT) :: tc
  CHARACTER(*), INTENT(IN) :: msg
  !----------------------------------------------
  ! Local Variables
  !----------------------------------------------
  DOUBLE PRECISION :: avg

  !----------------------------------------------
  avg = 0
  IF (tc%cnt .GT. 0) avg = tc%tm / tc%cnt
  IF(KPDBG) WRITE(OFU,'(A,I5.1,A,F8.4,A,F8.4,A)') msg,tc%cnt,' * ',avg,' sec = ',tc%tm,' sec'
END SUBROUTINE TimeCountPrint

!-------------------------------------------------------------------------------
!>
!! Init statistics
!<
!-------------------------------------------------------------------------------
SUBROUTINE PLBInitStats()
  CALL TimeCountInit( pstats%wait )
  CALL TimeCountInit( pstats%comm )
  CALL TimeCountInit( pstats%comp )
  CALL TimeCountInit( pstats%mm )
  CALL TimeCountInit( pstats%trf )
  CALL TimeCountInit( pstats%pmm )
  CALL TimeCountInit( pstats%ptrf )

  CALL TimeCountInit( pstats%mma )
  CALL TimeCountInit( pstats%mmb )
  CALL TimeCountInit( pstats%mmc )
  CALL TimeCountInit( pstats%mmalpha )
  CALL TimeCountInit( pstats%mmbeta )
  CALL TimeCountInit( pstats%mmrc )
  CALL TimeCountInit( pstats%extract )
  CALL TimeCountInit( pstats%waitall )
END SUBROUTINE PLBInitStats

!-------------------------------------------------------------------------------
!>
!! Print statistics so far
!<
!-------------------------------------------------------------------------------
SUBROUTINE PLBPrintStats()
  CALL TimeCountPrint( pstats%wait, 'PBLAS Wait ' )
  CALL TimeCountPrint( pstats%comm, 'PBLAS Comm ' )
  CALL TimeCountPrint( pstats%comp, 'PBLAS Comp ' )
  CALL TimeCountPrint( pstats%mm, 'PBLAS MM ' )
  CALL TimeCountPrint( pstats%trf, 'PBLAS TRF ' )
  CALL TimeCountPrint( pstats%pmm, 'PBLAS PMM ' )
  CALL TimeCountPrint( pstats%ptrf, 'PBLAS PTRF ' )

  CALL TimeCountPrint( pstats%mma, 'PBLAS MMA ' )
  CALL TimeCountPrint( pstats%mmb, 'PBLAS MMB ' )
  CALL TimeCountPrint( pstats%mmc, 'PBLAS MMC ' )
  CALL TimeCountPrint( pstats%mmalpha, 'PBLAS MMalpha ' )
  CALL TimeCountPrint( pstats%mmbeta, 'PBLAS MMbeta ' )
  CALL TimeCountPrint( pstats%mmrc, 'PBLAS MMRC ' )
  CALL TimeCountPrint( pstats%extract, 'PBLAS Extract ' )
  CALL TimeCountPrint( pstats%waitall, 'PBLAS Waitall ' )
END SUBROUTINE PLBPrintStats

!-------------------------------------------------------------------------------
!>
!! Initialize before foward or backward starts
!<
!-------------------------------------------------------------------------------
SUBROUTINE PLBInitialize
  !----------------------------------------------
  ! Local Variables
  !----------------------------------------------
  CHARACTER*100 :: envvar
  CHARACTER*100 :: envval

  LPLB_INIT=.TRUE.
  IF(KPDBG) WRITE(OFU,*) 'PLBInitialize Started'; CALL FL(OFU)

  !----------------------------------------------
  doblasonly = .TRUE.
  IF ( M .GE. 2048 ) THEN
    doblasonly = .TRUE.
    envvar = 'BLOCKTRI_BLASONLY'
    CALL GETENV( envvar, envval )
    IF ( envval .EQ. 'TRUE' ) THEN
      doblasonly = .TRUE.
      IF(KPDBG) WRITE(OFU,*) 'BLAS ONLY -- obeying env var ', envvar; CALL FL(OFU)
    END IF
  END IF
  IF(KPDBG) WRITE(OFU,*) 'doblasonly = ', doblasonly; CALL FL(OFU)

  !----------------------------------------------
  doblacscomm = .FALSE.
  envvar = 'BLOCKTRI_BLACSCOMM'
  CALL GETENV( envvar, envval )
  IF ( envval .EQ. 'TRUE' ) THEN
    doblacscomm = .TRUE.
    IF(KPDBG) WRITE(OFU,*) 'BLACS COMM -- obeying env var ', envvar; CALL FL(OFU)
  END IF
  IF(KPDBG) WRITE(OFU,*) 'doblacscomm = ', doblacscomm; CALL FL(OFU)

  !----------------------------------------------
  blacs%nbpp = 1
  envvar = 'BLOCKTRI_NBPP'
  CALL GETENV( envvar, envval )
  IF ( envval .NE. '' ) THEN
    READ( envval, *) blacs%nbpp
    IF(KPDBG) WRITE(OFU,*) 'NBPP -- obeying env var ', envvar; CALL FL(OFU)
  END IF
  IF(KPDBG) WRITE(OFU,*) 'NBPP = ', blacs%nbpp; CALL FL(OFU)

  !----------------------------------------------
  CALL PLBInitStats()

  !----------------------------------------------
  IF (doblasonly) THEN
    IF(KPDBG) WRITE(OFU,*) 'BLAS only (not using PBLAS)'; CALL FL(OFU)
  ELSE
#if defined(MPI_OPT)
    CALL BLACS_PINFO( blacs%iam, blacs%nprocs )
    IF(KPDBG) WRITE(OFU,*) 'BLACS_PINFO ', blacs%iam, ' ', blacs%nprocs; CALL FL(OFU)
    CALL BLACS_GET( 0, 0, blacs%maincontext )
    IF(KPDBG) WRITE(OFU,*) 'BLACS_GET ', blacs%maincontext; CALL FL(OFU)
    CALL BLACS_GRIDINIT( blacs%maincontext, 'R', 1, blacs%nprocs )
    IF(KPDBG) WRITE(OFU,*) 'BLACS_GRIDINIT'; CALL FL(OFU)
    CALL BLACS_BARRIER( blacs%maincontext, 'All' )
#endif
  END IF

  IF(KPDBG) WRITE(OFU,*) 'PLBInitialize Done'; CALL FL(OFU)
END SUBROUTINE PLBInitialize

!-------------------------------------------------------------------------------
!>
!! Finalize after foward or backward are done
!<
!-------------------------------------------------------------------------------
SUBROUTINE PLBFinalize

  IF(.NOT.LPLB_INIT) RETURN

  IF(KPDBG) WRITE(OFU,*) 'PLBFinalize Started'; CALL FL(OFU)
  IF (doblasonly) THEN
    IF(KPDBG) WRITE(OFU,*) 'BLAS only (not using PBLAS)'; CALL FL(OFU)
  ELSE
#if defined(MPI_OPT)
    CALL BLACS_BARRIER( blacs%maincontext, 'All' )
    IF(KPDBG) WRITE(OFU,*) 'BLACS_BARRIER'; CALL FL(OFU)
    CALL BLACS_EXIT(1)
    IF(KPDBG) WRITE(OFU,*) 'BLACS_EXIT nprocs= ', blacs%nprocs; CALL FL(OFU)
    CALL PLBPrintStats()
#endif
  END IF
  IF(KPDBG) WRITE(OFU,*) 'PLBFinalize Done'; CALL FL(OFU)
END SUBROUTINE PLBFinalize

!-------------------------------------------------------------------------------
!>
!! Initialize before any forward level starts
!<
!-------------------------------------------------------------------------------
SUBROUTINE PLBForwardInitialize
  !-----------------------------------------------------------------------------
  ! Nothing to do
END SUBROUTINE PLBForwardInitialize

!-------------------------------------------------------------------------------
!>
!! Finalize after all forward levels are done
!<
!-------------------------------------------------------------------------------
SUBROUTINE PLBForwardFinalize
  !-----------------------------------------------------------------------------
  ! Nothing to do
END SUBROUTINE PLBForwardFinalize

!-------------------------------------------------------------------------------
!>
!! Initialize before any backward level starts
!<
!-------------------------------------------------------------------------------
SUBROUTINE PLBBackwardInitialize
  !-----------------------------------------------------------------------------
  ! Nothing to do
END SUBROUTINE PLBBackwardInitialize

!-------------------------------------------------------------------------------
!>
!! Finalize after all backward levels are done
!<
!-------------------------------------------------------------------------------
SUBROUTINE PLBBackwardFinalize
  !-----------------------------------------------------------------------------
  ! Nothing to do
END SUBROUTINE PLBBackwardFinalize

!-------------------------------------------------------------------------------
!>
!! Initialize before the forward solve at a given level starts
!<
!-------------------------------------------------------------------------------
SUBROUTINE PLBForwardInitializeLevel( lvl, ammaster )
  !----------------------------------------------
  ! Formal arguments
  !----------------------------------------------
  INTEGER :: lvl
  LOGICAL :: ammaster
  !----------------------------------------------
  ! Local Variables
  !----------------------------------------------
  INTEGER :: I, J, K
  INTEGER :: maxslaves, actualslaves

  !----------------------------------------------
  !Short-circuit PBLAS to BLAS
  IF ( doblasonly ) THEN
    IF(KPDBG) WRITE(OFU,*) 'PLBForwardInitializeLevel BLAS only'; CALL FL(OFU)
    RETURN
  END IF

  IF(KPDBG) WRITE(OFU,*) 'PLBForwardInitializeLevel Started', ammaster; CALL FL(OFU)

#if defined(MPI_OPT)
  !----------------------------------------------
  pblas%ammaster = ammaster
  pblas%msmap%masterrank = -1
  pblas%msmap%nslaves = 0
  CALL DetermineMasterSlaveRanks()

  !----------------------------------------------
  !Set up BLACS process grid and a new context
  blacs%levelcontext = -1
  CALL BLACS_GET( blacs%maincontext, 10, blacs%levelcontext )
  IF(KPDBG) WRITE(OFU,*) 'Created new context'; CALL FL(OFU)
  IF(KPDBG) WRITE(OFU,*) 'Nslaves ',pblas%msmap%nslaves; CALL FL(OFU)

  blacs%pgrid%blockszrows = 64
  IF ( blacs%pgrid%blockszrows .GT. M ) blacs%pgrid%blockszrows = M
  blacs%pgrid%blockszcols = blacs%pgrid%blockszrows !<Portability => square
  IF(KPDBG) WRITE(OFU,*) 'Block NR=',blacs%pgrid%blockszrows; CALL FL(OFU)
  IF(KPDBG) WRITE(OFU,*) 'Block NC=',blacs%pgrid%blockszcols; CALL FL(OFU)

  maxslaves = (M*M) / (blacs%nbpp*blacs%nbpp* &
                       blacs%pgrid%blockszrows*blacs%pgrid%blockszcols)
  IF(KPDBG) WRITE(OFU,*) 'Max slaves ', maxslaves; CALL FL(OFU)
  IF ( maxslaves .LT. 1 ) maxslaves = 1
  IF(KPDBG) WRITE(OFU,*) 'Max slaves ', maxslaves; CALL FL(OFU)
  IF ( maxslaves .GT. pblas%msmap%nslaves ) maxslaves = pblas%msmap%nslaves
  IF(KPDBG) WRITE(OFU,*) 'Max slaves ', maxslaves; CALL FL(OFU)
  actualslaves = maxslaves
  IF(KPDBG) WRITE(OFU,*) ' Actual slaves ', actualslaves; CALL FL(OFU)

  blacs%pgrid%nrows = INT( SQRT( REAL( actualslaves ) ) )
  IF ( blacs%pgrid%nrows .LT. 1 ) blacs%pgrid%nrows = 1
  blacs%pgrid%ncols = INT( actualslaves / blacs%pgrid%nrows )
  IF(KPDBG) WRITE(OFU,*) 'NR=',blacs%pgrid%nrows,' NC=',blacs%pgrid%ncols; CALL FL(OFU)
  ALLOCATE( blacs%pgrid%map( 1 : blacs%pgrid%nrows, 1 : blacs%pgrid%ncols ) )
  K = 0
  DO I = 1, blacs%pgrid%nrows
    DO J = 1, blacs%pgrid%ncols
      blacs%pgrid%map(I,J) = pblas%msmap%slaveranks(K+1)
      K = K + 1
    END DO
  END DO
  IF(KPDBG) WRITE(OFU,*) 'NR*NC=',K; CALL FL(OFU)
  CALL BLACS_GRIDMAP( blacs%levelcontext, blacs%pgrid%map, &
                      blacs%pgrid%nrows, blacs%pgrid%nrows, blacs%pgrid%ncols )
  IF(KPDBG) WRITE(OFU,*) 'GridMap done'; CALL FL(OFU)
  CALL BLACS_GRIDINFO( blacs%levelcontext, &
                       blacs%pgrid%nrows, blacs%pgrid%ncols, &
                       blacs%pgrid%myrow, blacs%pgrid%mycol )
  IF(KPDBG) WRITE(OFU,*) 'GridInfo done'; CALL FL(OFU)
  IF(KPDBG) WRITE(OFU,*) 'Myrowcol ',blacs%pgrid%myrow,' ',blacs%pgrid%mycol; CALL FL(OFU)

  !----------------------------------------------
  pblas%mpicomm = SIESTA_COMM
  pblas%mpitag = 1234

  !----------------------------------------------
  CALL BLACS_BARRIER( blacs%maincontext, 'All' )
  IF(KPDBG) WRITE(OFU,*) 'BLACS_BARRIER'; CALL FL(OFU)

  !----------------------------------------------
  IF ( ammaster ) THEN
    IF(KPDBG) WRITE(OFU,*) 'PLBForwardInitializeLevel Master'; CALL FL(OFU)
  ELSE
    IF(KPDBG) WRITE(OFU,*) 'PLBForwardInitializeLevel Slave'; CALL FL(OFU)
    CALL SlaveService()
  END IF

  IF(KPDBG) WRITE(OFU,*) 'PLBForwardInitializeLevel Done', ammaster; CALL FL(OFU)
#endif

END SUBROUTINE PLBForwardInitializeLevel

!-------------------------------------------------------------------------------
!>
!! Finalize after the forward solve at a given level ends
!<
!-------------------------------------------------------------------------------
#if defined(MPI_OPT)
SUBROUTINE PLBForwardFinalizeLevel( lvl, ammaster )
  !----------------------------------------------
  ! Formal arguments
  !----------------------------------------------
  INTEGER :: lvl
  LOGICAL :: ammaster

  !----------------------------------------------
  !Short-circuit PBLAS to BLAS
  IF ( doblasonly ) THEN
    IF(KPDBG) WRITE(OFU,*) 'PLBForwardFinalizeLevel BLAS only'; CALL FL(OFU)
    RETURN
  END IF
  !----------------------------------------------
  !If I'm master, tell slaves that they're done
  IF ( ammaster .AND. (pblas%msmap%nslaves .GT. 1) ) THEN
    CALL MasterBcastNextOp( OP_DONE )
  END IF

  IF ( (blacs%pgrid%myrow .LT. 0) .OR. &
       (blacs%pgrid%myrow .GE. blacs%pgrid%nrows) ) THEN
    IF(KPDBG) WRITE(OFU,*) 'PLBForwardFinalizeLevel pariah !level-barrier';CALL FL(OFU)
  ELSE
    IF(KPDBG) WRITE(OFU,*) 'PLBForwardFinalizeLevel level-barrier'; CALL FL(OFU)
    CALL BLACS_BARRIER( blacs%levelcontext, 'All' )
    IF(KPDBG) WRITE(OFU,*) 'PLBForwardFinalizeLevel level grid exit'; CALL FL(OFU)
    CALL BLACS_GRIDEXIT( blacs%levelcontext )
  END IF

  IF(KPDBG) WRITE(OFU,*) 'PLBForwardFinalizeLevel main-barrier'; CALL FL(OFU)
  CALL BLACS_BARRIER( blacs%maincontext, 'All' )

  !----------------------------------------------
  pblas%msmap%masterrank = -1
  pblas%msmap%nslaves = 0
  DEALLOCATE( pblas%msmap%slaveranks )
  DEALLOCATE( blacs%pgrid%map )

  !----------------------------------------------
  CALL PLBPrintStats()

  IF(KPDBG) WRITE(OFU,*) 'PLBForwardFinalizeLevel ', ammaster; CALL FL(OFU)

END SUBROUTINE PLBForwardFinalizeLevel
#endif

!-------------------------------------------------------------------------------
!>
!! Initialize before the backward solve at a given level starts
!<
!-------------------------------------------------------------------------------
SUBROUTINE PLBBackwardInitializeLevel( lvl, ammaster )
  !----------------------------------------------
  ! Formal arguments
  !----------------------------------------------
  INTEGER :: lvl
  LOGICAL :: ammaster

  !----------------------------------------------
  ! Nothing to do
END SUBROUTINE PLBBackwardInitializeLevel

!-------------------------------------------------------------------------------
!>
!! Finalize after the backward solve at a given level ends
!<
!-------------------------------------------------------------------------------
SUBROUTINE PLBBackwardFinalizeLevel( lvl, ammaster )
  !----------------------------------------------
  ! Formal arguments
  !----------------------------------------------
  INTEGER :: lvl
  LOGICAL :: ammaster

  !----------------------------------------------
  !Nothing to do
END SUBROUTINE PLBBackwardFinalizeLevel

!-------------------------------------------------------------------------------
!>
!! Determine the ranks of master and slaves at this level
!<
!-------------------------------------------------------------------------------
#if defined(MPI_OPT)
SUBROUTINE DetermineMasterSlaveRanks()
  !----------------------------------------------
  ! Local Variables
  !----------------------------------------------
  LOGICAL, ALLOCATABLE :: recv_ammaster(:), assigned(:)
  TYPE(MasterToSlaveMapping), ALLOCATABLE :: allmsmaps(:)
  INTEGER :: totmasters, nslavespermaster, adjustednslavespermaster, tempmaster
  INTEGER :: mpi_err, I, J, masterindex

  IF(KPDBG) WRITE(OFU,*) 'DetermineMasterSlaveRanks started'; CALL FL(OFU)

  !----------------------------------------------
  !Initialize temporary data structures used for mapping
  ALLOCATE( recv_ammaster(1:nranks) )
  ALLOCATE( assigned(1:nranks) )
  DO I = 1, nranks
    assigned( i ) = .FALSE.
  END DO
  ALLOCATE( allmsmaps(1:nranks) )
  DO I = 1, nranks
    allmsmaps(I)%masterrank = -1
    allmsmaps(I)%nslaves = 0
  END DO

  !----------------------------------------------
  CALL MPI_Allgather(pblas%ammaster, 1, MPI_LOGICAL, &
                     recv_ammaster, 1, MPI_LOGICAL,  &
                     SIESTA_COMM, mpi_err)

  !----------------------------------------------
  !Figure out total number of masters
  totmasters = 0
  DO I = 1, nranks
    IF(KPDBG) WRITE(OFU,*) '    recv_ammaster ',I,' ',recv_ammaster(I); CALL FL(OFU)
    IF ( recv_ammaster(I) ) THEN
      totmasters = totmasters + 1
    END IF
  END DO

  IF ( totmasters .LE. 0 ) THEN
    IF(KPDBG) WRITE(OFU,*) 'Total masters must be greater than 0'; CALL FL(OFU)
    CALL ASSERT(.FALSE.,'blocktri #1')
  END IF

  nslavespermaster = (nranks / totmasters)
  IF(KPDBG) WRITE(OFU,*) 'Avg nslavespermaster',nslavespermaster; CALL FL(OFU)

  !----------------------------------------------
  !Assign slaves to each master
  tempmaster = 0
  masterloop: DO I = 1, nranks
    IF ( recv_ammaster(I) ) THEN

      tempmaster = tempmaster + 1

      !-------------------------------------------------------------
      !Give one more slave to earlier ranks, if not evenly divisible
      adjustednslavespermaster = nslavespermaster
      IF ( tempmaster .LE. MOD( nranks, totmasters ) ) THEN
        adjustednslavespermaster = adjustednslavespermaster + 1
      END IF

      IF(KPDBG) WRITE(OFU,*) 'Adjusted nslavespm',adjustednslavespermaster; CALL FL(OFU)
      allmsmaps(I)%masterrank = I-1
      ALLOCATE( allmsmaps(I)%slaveranks(1:adjustednslavespermaster) )

      !-------------------------------------------------------------
      !Add the master as first in its own slave rank set
      assigned(I) = .TRUE.
      allmsmaps(I)%nslaves = 1
      allmsmaps(I)%slaveranks(1) = I-1

      !-------------------------------------------------------------
      !Assign the next block of unassigned slaves
      slaveloop: DO J = 1, nranks
        IF ( allmsmaps(I)%nslaves .GE. adjustednslavespermaster ) THEN
          EXIT slaveloop
        END IF
        IF ( (.NOT. assigned(J)) .AND. (.NOT. recv_ammaster(J)) ) THEN
          assigned(J) = .TRUE.
          allmsmaps(J)%masterrank = I-1
          allmsmaps(I)%nslaves = allmsmaps(I)%nslaves + 1
          allmsmaps(I)%slaveranks(allmsmaps(I)%nslaves) = J-1
        END IF
      END DO slaveloop
    END IF
  END DO masterloop
  IF(KPDBG) WRITE(OFU,*) 'Computed all mappings'; CALL FL(OFU)

  !----------------------------------------------
  !Copy the mapping relevant to us (if I'm master, my own, else, my master's)
  masterindex = rank+1
  IF ( allmsmaps(masterindex)%masterrank .NE. rank ) THEN
    masterindex = allmsmaps(masterindex)%masterrank+1
  END IF
  pblas%msmap%masterrank = allmsmaps(masterindex)%masterrank
  pblas%msmap%nslaves = allmsmaps(masterindex)%nslaves
  ALLOCATE(pblas%msmap%slaveranks(1:allmsmaps(masterindex)%nslaves))
  pblas%msmap%slaveranks(:) = allmsmaps(masterindex)%slaveranks(:)

  IF(KPDBG) WRITE(OFU,*) 'Extracted my mappings'; CALL FL(OFU)
  IF(KPDBG) WRITE(OFU,*) 'My master rank',masterindex-1; CALL FL(OFU)
  IF(KPDBG) WRITE(OFU,*) 'Nslaves in my set',allmsmaps(masterindex)%nslaves; CALL FL(OFU)
  DO J = 1, allmsmaps(masterindex)%nslaves
    IF(KPDBG) WRITE(OFU,*) 'Slave',J,' ',allmsmaps(masterindex)%slaveranks(J);CALL FL(OFU)
  END DO

  !----------------------------------------------
  !Deallocations of temporary space
  DO I = 1, nranks
    IF ( recv_ammaster(I) ) THEN
      DEALLOCATE( allmsmaps(I)%slaveranks )
    END IF
  END DO
  DEALLOCATE( assigned )
  DEALLOCATE( allmsmaps )
  DEALLOCATE( recv_ammaster )

  IF(KPDBG) WRITE(OFU,*) 'DetermineMasterSlaveRanks done'; CALL FL(OFU)
END SUBROUTINE DetermineMasterSlaveRanks
#endif

!-------------------------------------------------------------------------------
!>
!! Inform the next operation to be done (calling rank is in master mode)
!<
!-------------------------------------------------------------------------------
SUBROUTINE MasterBcastNextOp( nextop )
  !----------------------------------------------
  ! Formal arguments
  !----------------------------------------------
  INTEGER, INTENT(IN) :: nextop
  !----------------------------------------------
  ! Local variables
  !----------------------------------------------
  INTEGER :: K, destrank
  REAL(dp) :: realnextop

  realnextop = REAL(nextop)
  IF(KPDBG) WRITE(OFU,*) 'MasterBcastNextOp started ', nextop; CALL FL(OFU)
  CALL MasterBcastValue( realnextop )
  IF(KPDBG) WRITE(OFU,*) 'MasterBcastNextOp done ', nextop; CALL FL(OFU)
END SUBROUTINE MasterBcastNextOp

!-------------------------------------------------------------------------------
!>
!! Determine the next operation to be done (calling rank is in slave mode)
!<
!-------------------------------------------------------------------------------
SUBROUTINE SlaveGetNextOp( nextop )
  !----------------------------------------------
  ! Formal arguments
  !----------------------------------------------
  REAL(dp) :: realnextop
  INTEGER, INTENT(INOUT) :: nextop

  IF(KPDBG) WRITE(OFU,*) 'SlaveGetNextOp started'; CALL FL(OFU)
  CALL SlaveReceiveValue( realnextop )
  nextop = INT( realnextop )
  IF(KPDBG) WRITE(OFU,*) 'SlaveGetNextOp done ', nextop; CALL FL(OFU)
END SUBROUTINE SlaveGetNextOp

!-------------------------------------------------------------------------------
!>
!! Slave's operation loop
!<
!-------------------------------------------------------------------------------
SUBROUTINE SlaveService()
  !----------------------------------------------
  ! Local variables
  !----------------------------------------------
  INTEGER :: nextop

  !----------------------------------------------
  IF ( (blacs%pgrid%myrow .LT. 0) .OR. &
       (blacs%pgrid%myrow .GE. blacs%pgrid%nrows) ) THEN
    IF(KPDBG) WRITE(OFU,*) 'SlaveService pariah falling out '; CALL FL(OFU)
  ELSE
    IF(KPDBG) WRITE(OFU,*) 'SlaveService started '; CALL FL(OFU)
    OpLoop: DO WHILE (.TRUE.)
      nextop = OP_NONE
      CALL SlaveGetNextOp( nextop )
      IF ( nextop .EQ. OP_DONE ) THEN
        EXIT OpLoop
      ELSE IF ( nextop .EQ. OP_DGEMM ) THEN
        CALL SlaveDGEMM()
      ELSE IF ( nextop .EQ. OP_DGETRF ) THEN
        CALL SlaveDGETRF()
      ELSE IF ( nextop .EQ. OP_DGETRS ) THEN
        CALL SlaveDGETRS()
      ELSE
        IF(KPDBG) WRITE(OFU,*) 'Bad Next Op', nextop; CALL FL(OFU)
        CALL ASSERT(.FALSE.,'blocktri #2')
      END IF
    END DO OpLoop
    IF(KPDBG) WRITE(OFU,*) 'SlaveService done '; CALL FL(OFU)
  END IF
END SUBROUTINE SlaveService

!-------------------------------------------------------------------------------
!>
!! Extracts from A the submatrix corresponding to process-grid element (pi,pj)
!<
!-------------------------------------------------------------------------------
SUBROUTINE ExtractSubMatrix( bszr, bszc, pnr, pnc, pi, pj, &
                             A, nrows, ncols, subA, subnrows, subncols )
  !----------------------------------------------
  ! Formal arguments
  !----------------------------------------------
  INTEGER, INTENT(IN) :: bszr, bszc
  INTEGER, INTENT(IN) :: pnr, pnc
  INTEGER, INTENT(IN) :: pi, pj
  REAL(dp), INTENT(IN) :: A(:,:)
  INTEGER, INTENT(IN) :: nrows, ncols
  REAL(dp), INTENT(OUT) :: subA(:)
  INTEGER, INTENT(IN) :: subnrows, subncols
  !----------------------------------------------
  ! Local variables
  !----------------------------------------------
  INTEGER :: I, J, K, Q, R
  INTEGER :: thisblkr, thisblkc

  IF(KPDBG) WRITE(OFU,*) 'ExtractSubMatrix NR=', subnrows, ' NC=', subncols; CALL FL(OFU)

  CALL BSystemClock(pstats%extract%t1)

  K = 0
  DO J = 1, ncols, bszc !<Col-major order
    thisblkc = (J-1) / bszc
    IF ( MOD(thisblkc, pnc) .EQ. pj-1 ) THEN
      DO R = 1, bszc
        IF ( J+R-1 .LE. ncols ) THEN
          DO I = 1, nrows, bszr
            thisblkr = (I-1) / bszr
            IF ( MOD(thisblkr, pnr) .EQ. pi-1 )  THEN
              DO Q = 1, bszr
                IF ( I+Q-1 .LE. nrows ) THEN
                  K = K + 1
                  subA(K) = A(I+Q-1,J+R-1)
                END IF
              END DO
            END IF
          END DO
        END IF
      END DO
    END IF
  END DO

  !----------------------------------------------
  IF ( K .NE. subnrows*subncols ) THEN
    IF(KPDBG) WRITE(OFU,*) 'Sanity check failed '; CALL FL(OFU)
    IF(KPDBG) WRITE(OFU,*) 'K=', K, ' subnr=', subnrows, ' subnc=', subncols; CALL FL(OFU)
    CALL ASSERT(.FALSE.,'blocktri #3')
  END IF

  CALL BSystemClock(pstats%extract%t2)
  CALL ChargeTime( pstats%extract%tm, pstats%extract%t2, pstats%extract%t1, pstats%extract%cnt )

  IF(KPDBG) WRITE(OFU,*) 'ExtractSubMatrix done K', K; CALL FL(OFU)

END SUBROUTINE ExtractSubMatrix

!-------------------------------------------------------------------------------
!>
!! Send the matrix by master (caller) to all slaves
!! The submatrix belonging to the caller (master is also one of the slaves) is
!! returned in the ssub* arguments (i.e., sent to self)
!<
!-------------------------------------------------------------------------------
SUBROUTINE MasterSendMatrix( A, nrows, ncols, ssubA, ssubnrows, ssubncols )
  !----------------------------------------------
  ! Formal arguments
  !----------------------------------------------
  REAL(dp), INTENT(IN) :: A(:,:)
  INTEGER, INTENT(IN) :: nrows, ncols
  REAL(dp), INTENT(OUT) :: ssubA(:)
  INTEGER, INTENT(IN) :: ssubnrows, ssubncols
  !----------------------------------------------
  ! Local variables
  !----------------------------------------------
  INTEGER :: I, J, K
  INTEGER :: aslaverank
  TYPE(PBLASTempArray), ALLOCATABLE :: subA(:,:)
  INTEGER :: subnrows, subncols
  INTEGER :: maxsubnrows, maxsubncols
  INTEGER :: bszr, bszc
  INTEGER :: pnr, pnc
#if defined(MPI_OPT)
  INTEGER, EXTERNAL :: NUMROC
  INTEGER :: mpistatus(MPI_STATUS_SIZE) !<MPI_Waitany status
#endif

  !----------------------------------------------
  INTEGER :: mpinisends !<Number of MPI_Isend requests initiated
  INTEGER, ALLOCATABLE :: mpireqs(:) !<MPI_Isend request handles
  INTEGER :: mpierr !<Generic use for MPI
  INTEGER, ALLOCATABLE :: mpistatuses(:,:) !<(nslaves,MPI_STATUS_SIZE) statuses
  INTEGER :: waitmatchindex !<MPI_Waitany match index

  !----------------------------------------------
  IF(KPDBG) WRITE(OFU,*) 'MasterSendMatrix started'; CALL FL(OFU)
  CALL BSystemClock(pstats%comm%t1)

#if defined(MPI_OPT)
  bszr = blacs%pgrid%blockszrows; bszc = blacs%pgrid%blockszcols
  pnr = blacs%pgrid%nrows; pnc = blacs%pgrid%ncols

  !----------------------------------------------
  mpinisends = 0
  ALLOCATE( mpireqs( pnr * pnc ) )
  ALLOCATE( mpistatuses( pnr * pnc, MPI_STATUS_SIZE ) )

  !----------------------------------------------
  maxsubnrows = NUMROC( nrows, bszr, 0, 0, pnr )
  maxsubncols = NUMROC( ncols, bszc, 0, 0, pnc )
  ALLOCATE( subA( 1 : pnr, 1 : pnc ) )

  DO I = 1, pnr
    DO J = 1, pnc
      aslaverank = blacs%pgrid%map(I,J)

      subnrows = NUMROC( nrows, bszr, I-1, 0, pnr )
      subncols = NUMROC( ncols, bszc, J-1, 0, pnc )

      IF(KPDBG) WRITE(OFU,*) 'MasterSendMatrix to ',I,' ',J,' ',aslaverank; CALL FL(OFU)

      IF ( I .EQ. 1 .AND. J .EQ. 1 ) THEN
        !Self (local)
        IF ( aslaverank .NE. rank ) THEN
          IF(KPDBG) WRITE(OFU,*) 'Inconsistency in slave rank of master'; CALL FL(OFU)
          CALL ASSERT(.FALSE.,'blocktri #4')
        END IF
        IF ( (ssubnrows .NE. subnrows) .OR. (ssubncols .NE. subncols) ) THEN
          IF(KPDBG) WRITE(OFU,*) 'Inconsistency in ssub dimensions'; CALL FL(OFU)
          IF(KPDBG) WRITE(OFU,*) 'SSNR ', ssubnrows, ' SSNC ', ssubncols; CALL FL(OFU)
          IF(KPDBG) WRITE(OFU,*) 'SNR  ', subnrows, ' SNC  ', subncols; CALL FL(OFU)
          CALL ASSERT(.FALSE.,'blocktri #5')
        END IF

        !Return a copy of this submatrix (will be like sent+received to/by self)
        IF(KPDBG) WRITE(OFU,*) 'MasterSendMatrix extracting self submatrix'; CALL FL(OFU)
        CALL ExtractSubMatrix(bszr, bszc, pnr, pnc, &
                              I, J, A, nrows, ncols, ssubA, subnrows, subncols)
        IF(KPDBG) WRITE(OFU,*) 'MasterSendMatrix kept self submatrix'; CALL FL(OFU)
      ELSE
        !A remote slave; send its corresponding matrix fragment
        IF(KPDBG) WRITE(OFU,*) 'MasterSendMatrix extracting submatrix',I,J; CALL FL(OFU)
        ALLOCATE( subA(I,J)%temparray( subnrows*subncols ) )
        CALL ExtractSubMatrix(bszr, bszc, pnr, pnc, &
                              I, J, A, nrows, ncols, subA(I,J)%temparray, &
                              subnrows,subncols)
        IF(KPDBG) WRITE(OFU,*) 'MasterSendMatrix extracted submatrix'; CALL FL(OFU)
        IF ( doblacscomm ) THEN
          CALL DGESD2D( blacs%levelcontext, subnrows, subncols, &
                        subA(I,J)%temparray, subnrows, I-1, J-1 )
        ELSE
          mpinisends = mpinisends + 1
          CALL MPI_Isend( subA(I,J)%temparray, subnrows*subncols, MPI_REAL8, &
            aslaverank, pblas%mpitag, pblas%mpicomm, mpireqs(mpinisends), mpierr )
        END IF
        IF(KPDBG) WRITE(OFU,*) 'MasterSendMatrix sent slave submatrix'; CALL FL(OFU)
        IF(KPDBG) WRITE(OFU,*) 'MasterSendMatrix deallocated temp submatrix'; CALL FL(OFU)
      END IF
    END DO
  END DO

  CALL BSystemClock(pstats%waitall%t1)
  CALL MPI_Waitall( mpinisends, mpireqs, mpistatuses, mpierr )
  CALL BSystemClock(pstats%waitall%t2)
  CALL ChargeTime( pstats%waitall%tm, pstats%waitall%t2, pstats%waitall%t1, pstats%waitall%cnt )

  DO I = 1, pnr
    DO J = 1, pnc
      IF ( .NOT. (I .EQ. 1 .AND. J .EQ. 1) ) THEN
        DEALLOCATE( subA(I,J)%temparray )
      END IF
    END DO
  END DO
  DEALLOCATE( subA )
  DEALLOCATE( mpireqs )
  DEALLOCATE( mpistatuses )

  CALL BSystemClock(pstats%comm%t2)
  CALL ChargeTime( pstats%comm%tm, pstats%comm%t2, pstats%comm%t1, pstats%comm%cnt )
#endif
  IF(KPDBG) WRITE(OFU,*) 'MasterSendMatrix done'; CALL FL(OFU)

END SUBROUTINE MasterSendMatrix

!-------------------------------------------------------------------------------
!>
!! Receive the fragment of matrix sent by master to slave (caller)
!<
!-------------------------------------------------------------------------------
SUBROUTINE SlaveReceiveMatrix( nrows, ncols, subA, subnrows, subncols )
  !----------------------------------------------
  ! Formal arguments
  !----------------------------------------------
#if defined(MPI_OPT)
  INTEGER :: mpistatus(MPI_STATUS_SIZE)
#endif
  INTEGER, INTENT(IN) :: nrows, ncols !<of global matrix
  REAL(dp), INTENT(OUT) :: subA(:)
  INTEGER, INTENT(IN) :: subnrows, subncols
  !----------------------------------------------
  ! Local variables
  !----------------------------------------------
  INTEGER :: masterrank
  INTEGER :: mpierr

  IF(KPDBG) WRITE(OFU,*) 'SlaveReceiveMatrix started ',subnrows,' ',subncols; CALL FL(OFU)

#if defined(MPI_OPT)
  masterrank = blacs%pgrid%map(1,1)

  CALL BSystemClock(pstats%comm%t1)

  IF ( doblacscomm ) THEN
    CALL DGERV2D( blacs%levelcontext, subnrows, subncols, subA, subnrows, 0, 0 )
  ELSE
    CALL MPI_Recv( subA, subnrows*subncols, MPI_REAL8, &
                   masterrank, pblas%mpitag, pblas%mpicomm, mpistatus, mpierr )
  END IF

  CALL BSystemClock(pstats%comm%t2)
  CALL ChargeTime( pstats%comm%tm, pstats%comm%t2, pstats%comm%t1, pstats%comm%cnt )
#endif

  IF(KPDBG) WRITE(OFU,*) 'SlaveReceiveMatrix done'; CALL FL(OFU)
END SUBROUTINE SlaveReceiveMatrix

!-------------------------------------------------------------------------------
!>
!! Send the given value by master (caller) to all slaves
!<
!-------------------------------------------------------------------------------
SUBROUTINE MasterBcastValue( val )
  !----------------------------------------------
  ! Formal arguments
  !----------------------------------------------
  REAL(dp), INTENT(IN) :: val

#if defined(MPI_OPT)
  CALL DGEBS2D( blacs%levelcontext, 'All', ' ', 1, 1, val, 1 );
#endif
  IF(KPDBG) WRITE(OFU,*) 'MasterBcastValue bcast to slaves'; CALL FL(OFU)

END SUBROUTINE MasterBcastValue

!-------------------------------------------------------------------------------
!>
!! Receive a value sent by master to slave (caller)
!<
!-------------------------------------------------------------------------------
SUBROUTINE SlaveReceiveValue( val )
  !----------------------------------------------
  ! Formal arguments
  !----------------------------------------------
  REAL(dp), INTENT(OUT) :: val

#if defined(MPI_OPT)
  CALL DGEBR2D( blacs%levelcontext, 'All', ' ', 1, 1, val, 1, 0, 0 )
#endif
  IF(KPDBG) WRITE(OFU,*) 'SlaveReceiveValue bcast from master'
  CALL FL(OFU)

END SUBROUTINE SlaveReceiveValue

!-------------------------------------------------------------------------------
!>
!! Incorporates into A the submatrix corresp. to process-grid element (pi,pj)
!<
!-------------------------------------------------------------------------------
SUBROUTINE InjectSubMatrix( bszr, bszc, pnr, pnc, pi, pj, &
                            A, nrows, ncols, subA, subnrows, subncols )
  !----------------------------------------------
  ! Formal arguments
  !----------------------------------------------
  INTEGER, INTENT(IN) :: bszr, bszc
  INTEGER, INTENT(IN) :: pnr, pnc
  INTEGER, INTENT(IN) :: pi, pj
  REAL(dp), INTENT(OUT) :: A(:,:)
  INTEGER, INTENT(IN) :: nrows, ncols
  REAL(dp), INTENT(IN) :: subA(:)
  INTEGER, INTENT(IN) :: subnrows, subncols
  !----------------------------------------------
  ! Local variables
  !----------------------------------------------
  INTEGER :: I, J, K, Q, R
  INTEGER :: thisblkr, thisblkc

  IF(KPDBG) WRITE(OFU,*) 'InjectSubMatrix NR=', subnrows, ' NC=', subncols; CALL FL(OFU)

  K = 0
  DO J = 1, ncols, bszc !<Col-major order
    thisblkc = (J-1) / bszc
    IF ( MOD(thisblkc, pnc) .EQ. pj-1 ) THEN
      DO R = 1, bszc
        IF ( J+R-1 .LE. ncols ) THEN
          DO I = 1, nrows, bszr
            thisblkr = (I-1) / bszr
            IF ( MOD(thisblkr, pnr) .EQ. pi-1 )  THEN
              DO Q = 1, bszr
                IF ( I+Q-1 .LE. nrows ) THEN
                  K = K + 1
                  A(I+Q-1,J+R-1) = subA(K)
                END IF
              END DO
            END IF
          END DO
        END IF
      END DO
    END IF
  END DO

  !----------------------------------------------
  IF ( K .NE. subnrows*subncols ) THEN
    IF(KPDBG) WRITE(OFU,*) 'Sanity check failed '; CALL FL(OFU)
    IF(KPDBG) WRITE(OFU,*) 'K=', K, ' subnr=', subnrows, ' subnc=', subncols; CALL FL(OFU)
    CALL ASSERT(.FALSE.,'blocktri #6')
  END IF

  IF(KPDBG) WRITE(OFU,*) 'InjectSubMatrix done K', K; CALL FL(OFU)

END SUBROUTINE InjectSubMatrix

!-------------------------------------------------------------------------------
!>
!! Incorporates into V the subvector corresp. to process-grid element (pi,0)
!<
!-------------------------------------------------------------------------------
SUBROUTINE InjectSubVector( bszr, pnr, pi, V, nrows, subV, subnrows )
  !----------------------------------------------
  ! Formal arguments
  !----------------------------------------------
  INTEGER, INTENT(IN) :: bszr
  INTEGER, INTENT(IN) :: pnr
  INTEGER, INTENT(IN) :: pi
  INTEGER, INTENT(OUT) :: V(:)
  INTEGER, INTENT(IN) :: nrows
  INTEGER, INTENT(IN) :: subV(:)
  INTEGER, INTENT(IN) :: subnrows
  !----------------------------------------------
  ! Local variables
  !----------------------------------------------
  INTEGER :: I, K, Q
  INTEGER :: thisblkr

  IF(KPDBG) WRITE(OFU,*) 'InjectSubVector NR=', subnrows; CALL FL(OFU)

  K = 0
  DO I = 1, nrows, bszr
    thisblkr = (I-1) / bszr
    IF ( MOD(thisblkr, pnr) .EQ. pi-1 ) THEN
      DO Q = 1, bszr
        IF ( I+Q-1 .LE. nrows ) THEN
          K = K + 1
          V(I+Q-1) = subV(K)
        END IF
      END DO
    END IF
  END DO

  !----------------------------------------------
  IF ( K .NE. subnrows ) THEN
    IF(KPDBG) WRITE(OFU,*) 'Sanity check failed '; CALL FL(OFU)
    IF(KPDBG) WRITE(OFU,*) 'K=', K, ' subnr=', subnrows; CALL FL(OFU)
    CALL ASSERT(.FALSE.,'blocktri #7')
  END IF

  IF(KPDBG) WRITE(OFU,*) 'InjectSubVector done K', K; CALL FL(OFU)

END SUBROUTINE InjectSubVector

!-------------------------------------------------------------------------------
!>
!! Receive the submatrices by master (caller) from all slaves
!! The submatrix belonging to the caller (master is also one of the slaves) is
!! given in the ssub* arguments (i.e., received from self)
!<
!-------------------------------------------------------------------------------
SUBROUTINE MasterRecvMatrix( A, nrows, ncols, ssubA, ssubnrows, ssubncols )
  !----------------------------------------------
  ! Formal arguments
  !----------------------------------------------
  REAL(dp), INTENT(OUT) :: A(:,:)
  INTEGER, INTENT(IN) :: nrows, ncols
  REAL(dp), INTENT(IN) :: ssubA(:)
  INTEGER, INTENT(IN) :: ssubnrows, ssubncols
  !----------------------------------------------
  ! Local variables
  !----------------------------------------------
  INTEGER :: I, J, K
  INTEGER :: aslaverank
  REAL(dp), ALLOCATABLE :: subA(:)
  INTEGER :: subnrows, subncols
  INTEGER :: bszr, bszc
  INTEGER :: pnr, pnc
#if defined(MPI_OPT)
  INTEGER, EXTERNAL :: NUMROC
#endif

  IF(KPDBG) WRITE(OFU,*) 'MasterRecvMatrix started'; CALL FL(OFU)
  CALL BSystemClock(pstats%comm%t1)

#if defined(MPI_OPT)
  bszr = blacs%pgrid%blockszrows; bszc = blacs%pgrid%blockszcols
  pnr = blacs%pgrid%nrows; pnc = blacs%pgrid%ncols

  DO I = 1, pnr
    DO J = 1, pnc
      aslaverank = blacs%pgrid%map(I,J)

      subnrows = NUMROC( nrows, bszr, I-1, 0, pnr )
      subncols = NUMROC( ncols, bszc, J-1, 0, pnc )

      IF(KPDBG) WRITE(OFU,*) 'MasterRecvMatrix from ',I,' ',J,' ',aslaverank; CALL FL(OFU)

      IF ( I .EQ. 1 .AND. J .EQ. 1 ) THEN
        !Self (local)
        IF ( aslaverank .NE. rank ) THEN
          IF(KPDBG) WRITE(OFU,*) 'Inconsistency in slave rank of master'; CALL FL(OFU)
          CALL ASSERT(.FALSE.,'blocktri #8')
        END IF
        IF ( (ssubnrows .NE. subnrows) .OR. (ssubncols .NE. subncols) ) THEN
          IF(KPDBG) WRITE(OFU,*) 'Inconsistency in ssub dimensions'; CALL FL(OFU)
          IF(KPDBG) WRITE(OFU,*) 'SSNR ', ssubnrows, ' SSNC ', ssubncols; CALL FL(OFU)
          IF(KPDBG) WRITE(OFU,*) 'SNR  ', subnrows, ' SNC  ', subncols; CALL FL(OFU)
          CALL ASSERT(.FALSE.,'blocktri #9')
        END IF

        !Use the given copy of this submatrix (like sent+received to/by self)
        IF(KPDBG) WRITE(OFU,*) 'MasterRecvMatrix injecting self submatrix'; CALL FL(OFU)
        CALL InjectSubMatrix(bszr, bszc, pnr, pnc, &
                              I, J, A, nrows, ncols, ssubA, subnrows, subncols)
        IF(KPDBG) WRITE(OFU,*) 'MasterRecvMatrix kept self submatrix'; CALL FL(OFU)
      ELSE
        !A remote slave; receive its corresponding matrix fragment
        ALLOCATE( subA( 1 : subnrows*subncols ) )
        IF(KPDBG) WRITE(OFU,*) 'MasterRecvMatrix receiving slave submatrix'; CALL FL(OFU)
        CALL BSystemClock(pstats%wait%t1)
        CALL DGERV2D( blacs%levelcontext, subnrows, subncols, &
                      subA, subnrows, I-1, J-1 )
        CALL BSystemClock(pstats%wait%t2)
        CALL ChargeTime( pstats%wait%tm, pstats%wait%t2, pstats%wait%t1, pstats%wait%cnt )
        IF(KPDBG) WRITE(OFU,*) 'MasterRecvMatrix injecting submatrix',I,J; CALL FL(OFU)
        CALL InjectSubMatrix( bszr, bszc, pnr, pnc, &
                              I, J, A, nrows, ncols, subA, subnrows, subncols )
        IF(KPDBG) WRITE(OFU,*) 'MasterRecvMatrix injected submatrix'; CALL FL(OFU)
        DEALLOCATE( subA ) !<We don't need it locally anymore
        IF(KPDBG) WRITE(OFU,*) 'MasterRecvMatrix deallocated temp submatrix'; CALL FL(OFU)
      END IF
    END DO
  END DO

  CALL BSystemClock(pstats%comm%t2)
  CALL ChargeTime( pstats%comm%tm, pstats%comm%t2, pstats%comm%t1, pstats%comm%cnt )
#endif

  IF(KPDBG) WRITE(OFU,*) 'MasterRecvMatrix done'; CALL FL(OFU)

END SUBROUTINE MasterRecvMatrix

!-------------------------------------------------------------------------------
!>
!! Send the fragment of matrix by slave (caller) to master
!<
!-------------------------------------------------------------------------------
SUBROUTINE SlaveSendMatrix( nrows, ncols, subA, subnrows, subncols )
  !----------------------------------------------
  ! Formal arguments
  !----------------------------------------------
  INTEGER, INTENT(IN) :: nrows, ncols !<of global matrix
  REAL(dp), INTENT(IN) :: subA(:)
  INTEGER, INTENT(IN) :: subnrows, subncols

  IF(KPDBG) WRITE(OFU,*) 'SlaveSendMatrix started ',subnrows,' ',subncols; CALL FL(OFU)

  CALL BSystemClock(pstats%comm%t1)
#if defined(MPI_OPT)
  CALL DGESD2D( blacs%levelcontext, subnrows, subncols, subA, subnrows, 0, 0 )
#endif
  CALL BSystemClock(pstats%comm%t2)
  CALL ChargeTime( pstats%comm%tm, pstats%comm%t2, pstats%comm%t1, pstats%comm%cnt )

  IF(KPDBG) WRITE(OFU,*) 'SlaveSendMatrix done'; CALL FL(OFU)
END SUBROUTINE SlaveSendMatrix

!-------------------------------------------------------------------------------
!>
!! Receive the subvector by master (caller) from all slaves in column 1 of pgrid
!! The subvector belonging to the caller (master is also one of the slaves) is
!! given in the ssub* arguments (i.e., received from self)
!<
!-------------------------------------------------------------------------------
SUBROUTINE MasterRecvVector( V, nrows, ssubV, ssubnrows )
  !----------------------------------------------
  ! Formal arguments
  !----------------------------------------------
  INTEGER, INTENT(OUT) :: V(:)
  INTEGER, INTENT(IN) :: nrows
  INTEGER, INTENT(IN) :: ssubV(:)
  INTEGER, INTENT(IN) :: ssubnrows
  !----------------------------------------------
  ! Local variables
  !----------------------------------------------
  INTEGER :: I, J, K
  INTEGER :: aslaverank
  INTEGER, ALLOCATABLE :: subV(:)
  INTEGER :: subnrows
  INTEGER :: bszr
  INTEGER :: pnr
#if defined(MPI_OPT)
  INTEGER, EXTERNAL :: NUMROC
#endif

  IF(KPDBG) WRITE(OFU,*) 'MasterRecvVector started'; CALL FL(OFU)
  CALL BSystemClock(pstats%comm%t1)

#if defined(MPI_OPT)
  bszr = blacs%pgrid%blockszrows
  pnr = blacs%pgrid%nrows

  DO I = 1, pnr
    DO J = 1, 1
      aslaverank = blacs%pgrid%map(I,J)

      subnrows = NUMROC( nrows, bszr, I-1, 0, pnr )

      IF(KPDBG) WRITE(OFU,*) 'MasterRecvVector from ',I,' ',J,' ',aslaverank; CALL FL(OFU)

      IF ( I .EQ. 1 .AND. J .EQ. 1 ) THEN
        !Self (local)
        IF ( aslaverank .NE. rank ) THEN
          IF(KPDBG) WRITE(OFU,*) 'Inconsistency in slave rank of master'; CALL FL(OFU)
          CALL ASSERT(.FALSE.,'blocktri #10')
        END IF
        IF ( ssubnrows .NE. subnrows ) THEN
          IF(KPDBG) WRITE(OFU,*) 'Inconsistency in ssub dimensions'; CALL FL(OFU)
          IF(KPDBG) WRITE(OFU,*) 'SSNR ', ssubnrows; CALL FL(OFU)
          IF(KPDBG) WRITE(OFU,*) 'SNR  ', subnrows; CALL FL(OFU)
          CALL ASSERT(.FALSE.,'blocktri #11')
        END IF

        !Use the given copy of this subvector (like sent+received to/by self)
        IF(KPDBG) WRITE(OFU,*) 'MasterRecvVector injecting self subvector'; CALL FL(OFU)
        CALL InjectSubVector(bszr, pnr, I, V, nrows, ssubV, subnrows)
        IF(KPDBG) WRITE(OFU,*) 'MasterRecvVector kept self subvector'; CALL FL(OFU)
      ELSE
        !A remote slave; receive its corresponding vector fragment
        ALLOCATE( subV( 1 : subnrows ) )
        IF(KPDBG) WRITE(OFU,*) 'MasterRecvVector receiving slave subvector'; CALL FL(OFU)
        CALL IGERV2D( blacs%levelcontext, subnrows, 1, subV, subnrows, I-1, 0 )
        IF(KPDBG) WRITE(OFU,*) 'MasterRecvVector injecting subvector',I,J; CALL FL(OFU)
        CALL InjectSubVector( bszr, pnr, I, V, nrows, subV, subnrows )
        IF(KPDBG) WRITE(OFU,*) 'MasterRecvVector injected subvector'; CALL FL(OFU)
        DEALLOCATE( subV ) !<We don't need it locally anymore
        IF(KPDBG) WRITE(OFU,*) 'MasterRecvVector deallocated temp subvector'; CALL FL(OFU)
      END IF
    END DO
  END DO

  CALL BSystemClock(pstats%comm%t2)
  CALL ChargeTime( pstats%comm%tm, pstats%comm%t2, pstats%comm%t1, pstats%comm%cnt )
#endif

  IF(KPDBG) WRITE(OFU,*) 'MasterRecvVector done'; CALL FL(OFU)

END SUBROUTINE MasterRecvVector

!-------------------------------------------------------------------------------
!>
!! Send the fragment of vector by slave (caller) to master
!! Only slaves who belong to zero'th column of pgrid send to master
!! Other slave simply do no-op for this routine
!<
!-------------------------------------------------------------------------------
SUBROUTINE SlaveSendVector( nrows, subV, subnrows )
  !----------------------------------------------
  ! Formal arguments
  !----------------------------------------------
  INTEGER, INTENT(IN) :: nrows !<of global matrix
  INTEGER, INTENT(IN) :: subV(:)
  INTEGER, INTENT(IN) :: subnrows
  INTEGER :: pj

  IF(KPDBG) WRITE(OFU,*) 'SlaveSendVector started ', subnrows; CALL FL(OFU)

  CALL BSystemClock(pstats%comm%t1)

  pj = blacs%pgrid%mycol
  IF ( pj .GT. 0 ) THEN
    IF(KPDBG) WRITE(OFU,*) 'SlaveSendVector skipping since mycol>0 ',pj; CALL FL(OFU)
  ELSE
#if defined(MPI_OPT)
    CALL IGESD2D( blacs%levelcontext, subnrows, 1, subV, subnrows, 0, 0 )
#endif
  END IF

  CALL BSystemClock(pstats%comm%t2)
  CALL ChargeTime( pstats%comm%tm, pstats%comm%t2, pstats%comm%t1, pstats%comm%cnt )

  IF(KPDBG) WRITE(OFU,*) 'SlaveSendVector done'; CALL FL(OFU)
END SUBROUTINE SlaveSendVector

!-------------------------------------------------------------------------------
!>
!! Encapsulates BLAS' DGEMM functionality
!<
!-------------------------------------------------------------------------------
SUBROUTINE PLBDGEMM( alpha, A, B, beta, C )
  !----------------------------------------------
  ! Formal arguments
  !----------------------------------------------
  REAL(dp), INTENT(IN)    :: alpha, beta
  REAL(dp), INTENT(IN)    :: A(:,:), B(:,:) !<Each matrix of size MxM
  REAL(dp), INTENT(INOUT) :: C(:,:) !<Matrix of size MxM
  !----------------------------------------------
  ! Local variables
  !----------------------------------------------
  REAL(dp), ALLOCATABLE :: subA(:), subB(:), subC(:)
  INTEGER :: descA(9), descB(9), descC(9)
  INTEGER :: lldA, lldB, lldC
  INTEGER :: nrows, ncols, subnrows, subncols
  INTEGER :: bszr, bszc, pnr, pnc, pi, pj
  LOGICAL :: useblas
  INTEGER :: ctxt, info
#if defined(MPI_OPT)
  INTEGER, EXTERNAL :: NUMROC
#endif

  IF(KPDBG) WRITE(OFU,*) 'MasterGEMM started'; CALL FL(OFU)

  !----------------------------------------------
  useblas = ( doblasonly ) .OR. &
            ( pblas%msmap%nslaves .EQ. 1) .OR. &
            ( M .LE. blacs%pgrid%blockszrows ) .OR. &
            ( M .LE. blacs%pgrid%blockszcols )

  IF ( useblas ) THEN
    IF(KPDBG) WRITE(OFU,*) 'BLAS DGEMM only (not using PBLAS)'; CALL FL(OFU)
    CALL DGEMM( 'N', 'N', M, M, M, alpha, A, M, B, M, beta, C, M )
  ELSE
#if defined(MPI_OPT)
    CALL BSystemClock(pstats%mm%t1)
    nrows = M; ncols = M
    bszr = blacs%pgrid%blockszrows; bszc = blacs%pgrid%blockszcols
    pnr = blacs%pgrid%nrows; pnc = blacs%pgrid%ncols
    pi = blacs%pgrid%myrow; pj = blacs%pgrid%mycol
    subnrows = NUMROC( nrows, bszr, pi, 0, pnr )
    subncols = NUMROC( ncols, bszc, pj, 0, pnc )

    IF(KPDBG) WRITE(OFU,*) 'MasterDGEMM allocating subABC'; CALL FL(OFU)
    ALLOCATE( subA(1:subnrows*subncols) )
    ALLOCATE( subB(1:subnrows*subncols) )
    ALLOCATE( subC(1:subnrows*subncols) )
    IF(KPDBG) WRITE(OFU,*) 'MasterDGEMM allocated subABC'; CALL FL(OFU)

    IF(KPDBG) WRITE(OFU,*) 'MasterDGEMM desciniting subABC'; CALL FL(OFU)
    ctxt = blacs%levelcontext
    lldA = subnrows; IF( lldA .LT. 1 ) lldA = 1
    lldB = lldA; lldC = lldA
    CALL DESCINIT(descA, nrows, ncols, bszr, bszc, 0, 0, ctxt, lldA, info )
    CALL DESCINIT(descB, nrows, ncols, bszr, bszc, 0, 0, ctxt, lldB, info )
    CALL DESCINIT(descC, nrows, ncols, bszr, bszc, 0, 0, ctxt, lldC, info )
    IF(KPDBG) WRITE(OFU,*) 'MasterDGEMM desciniting subABC'; CALL FL(OFU)

    IF(KPDBG) WRITE(OFU,*) 'MasterDGEMM sending OP_DGEMM'; CALL FL(OFU)
    CALL MasterBcastNextOp( OP_DGEMM )

    IF(KPDBG) WRITE(OFU,*) 'MasterDGEMM sending A'; CALL FL(OFU)
    CALL BSystemClock(pstats%mma%t1)
    CALL MasterSendMatrix( A, nrows, ncols, subA, subnrows, subncols )
    CALL BSystemClock(pstats%mma%t2)
    CALL ChargeTime( pstats%mma%tm, pstats%mma%t2, pstats%mma%t1, pstats%mma%cnt )

    IF(KPDBG) WRITE(OFU,*) 'MasterDGEMM sending B'; CALL FL(OFU)
    CALL BSystemClock(pstats%mmb%t1)
    CALL MasterSendMatrix( B, nrows, ncols, subB, subnrows, subncols )
    CALL BSystemClock(pstats%mmb%t2)
    CALL ChargeTime( pstats%mmb%tm, pstats%mmb%t2, pstats%mmb%t1, pstats%mmb%cnt )

    IF(KPDBG) WRITE(OFU,*) 'MasterDGEMM sending C'; CALL FL(OFU)
    CALL BSystemClock(pstats%mmc%t1)
    CALL MasterSendMatrix( C, nrows, ncols, subC, subnrows, subncols )
    CALL BSystemClock(pstats%mmc%t2)
    CALL ChargeTime( pstats%mmc%tm, pstats%mmc%t2, pstats%mmc%t1, pstats%mmc%cnt )

    IF(KPDBG) WRITE(OFU,*) 'MasterDGEMM sending alpha', alpha; CALL FL(OFU)
    CALL BSystemClock(pstats%mmalpha%t1)
    CALL MasterBcastValue( alpha )
    CALL BSystemClock(pstats%mmalpha%t2)
    CALL ChargeTime( pstats%mmalpha%tm, pstats%mmalpha%t2, pstats%mmalpha%t1, pstats%mmalpha%cnt )

    IF(KPDBG) WRITE(OFU,*) 'MasterDGEMM sending beta', beta; CALL FL(OFU)
    CALL BSystemClock(pstats%mmbeta%t1)
    CALL MasterBcastValue( beta )
    CALL BSystemClock(pstats%mmbeta%t2)
    CALL ChargeTime( pstats%mmbeta%tm, pstats%mmbeta%t2, pstats%mmbeta%t1, pstats%mmbeta%cnt )

    IF(KPDBG) WRITE(OFU,*) 'MasterDGEMM invoking PDGEMM'; CALL FL(OFU)
    CALL BSystemClock(pstats%comp%t1)
    CALL PDGEMM( 'N', 'N', M, M, M, alpha, &
                 subA, 1, 1, descA, &
                 subB, 1, 1, descB, &
                 beta, &
                 subC, 1, 1, descC )
    CALL BSystemClock(pstats%comp%t2)
    CALL ChargeTime( pstats%comp%tm, pstats%comp%t2, pstats%comp%t1, pstats%comp%cnt )
    CALL ChargeTime( pstats%pmm%tm, pstats%comp%t2, pstats%comp%t1, pstats%pmm%cnt )
    IF(KPDBG) WRITE(OFU,*) 'MasterDGEMM done PDGEMM'; CALL FL(OFU)

    IF(KPDBG) WRITE(OFU,*) 'MasterDGEMM receiving slave matrices'; CALL FL(OFU)
    CALL BSystemClock(pstats%mmrc%t1)
    CALL MasterRecvMatrix( C, nrows, ncols, subC, subnrows, subncols )
    CALL BSystemClock(pstats%mmrc%t2)
    CALL ChargeTime( pstats%mmrc%tm, pstats%mmrc%t2, pstats%mmrc%t1, pstats%mmrc%cnt )
    IF(KPDBG) WRITE(OFU,*) 'MasterDGEMM received slave matrices'; CALL FL(OFU)

    IF(KPDBG) WRITE(OFU,*) 'MasterDGEMM deallocating subABC'; CALL FL(OFU)
    DEALLOCATE( subA )
    DEALLOCATE( subB )
    DEALLOCATE( subC )
    IF(KPDBG) WRITE(OFU,*) 'MasterDGEMM deallocated subABC'; CALL FL(OFU)
    CALL BSystemClock(pstats%mm%t2)
    CALL ChargeTime( pstats%mm%tm, pstats%mm%t2, pstats%mm%t1, pstats%mm%cnt )
#endif
  END IF

  IF(KPDBG) WRITE(OFU,*) 'MasterGEMM done'; CALL FL(OFU)
END SUBROUTINE PLBDGEMM

!-------------------------------------------------------------------------------
!>
!! DGEMM support from slave
!<
!-------------------------------------------------------------------------------
SUBROUTINE SlaveDGEMM()
  !----------------------------------------------
  ! Local variables
  !----------------------------------------------
  REAL(dp), ALLOCATABLE :: subA(:), subB(:), subC(:)
  REAL(dp) :: alpha, beta
  INTEGER :: descA(9), descB(9), descC(9)
  INTEGER :: lldA, lldB, lldC
  INTEGER :: nrows, ncols, subnrows, subncols
  INTEGER :: bszr, bszc, pnr, pnc
  INTEGER :: pi, pj
  INTEGER :: ctxt, info
#if defined(MPI_OPT)
  INTEGER, EXTERNAL :: NUMROC
#endif

  CALL BSystemClock(pstats%mm%t1)

#if defined(MPI_OPT)
  nrows = M; ncols = M
  bszr = blacs%pgrid%blockszrows; bszc = blacs%pgrid%blockszcols
  pnr = blacs%pgrid%nrows; pnc = blacs%pgrid%ncols
  pi = blacs%pgrid%myrow; pj = blacs%pgrid%mycol

  subnrows = NUMROC( nrows, bszr, pi, 0, pnr )
  subncols = NUMROC( ncols, bszc, pj, 0, pnc )

  IF(KPDBG) WRITE(OFU,*) 'SlaveDGEMM allocating subABC'; CALL FL(OFU)
  ALLOCATE( subA(1:subnrows*subncols) )
  ALLOCATE( subB(1:subnrows*subncols) )
  ALLOCATE( subC(1:subnrows*subncols) )
  IF(KPDBG) WRITE(OFU,*) 'SlaveDGEMM allocated subABC'; CALL FL(OFU)

  IF(KPDBG) WRITE(OFU,*) 'SlaveDGEMM desciniting subABC'; CALL FL(OFU)
  ctxt = blacs%levelcontext
  lldA = subnrows; IF( lldA .LT. 1 ) lldA = 1
  lldB = lldA; lldC = lldA
  CALL DESCINIT(descA, nrows, ncols, bszr, bszc, 0, 0, ctxt, lldA, info )
  CALL DESCINIT(descB, nrows, ncols, bszr, bszc, 0, 0, ctxt, lldB, info )
  CALL DESCINIT(descC, nrows, ncols, bszr, bszc, 0, 0, ctxt, lldC, info )
  IF(KPDBG) WRITE(OFU,*) 'SlaveDGEMM desciniting subABC'; CALL FL(OFU)

  IF(KPDBG) WRITE(OFU,*) 'SlaveDGEMM receiving A'; CALL FL(OFU)
  CALL BSystemClock(pstats%mma%t1)
  CALL SlaveReceiveMatrix( nrows, ncols, subA, subnrows, subncols )
  CALL BSystemClock(pstats%mma%t2)
  CALL ChargeTime( pstats%mma%tm, pstats%mma%t2, pstats%mma%t1, pstats%mma%cnt )

  IF(KPDBG) WRITE(OFU,*) 'SlaveDGEMM receiving B'; CALL FL(OFU)
  CALL BSystemClock(pstats%mmb%t1)
  CALL SlaveReceiveMatrix( nrows, ncols, subB, subnrows, subncols )
  CALL BSystemClock(pstats%mmb%t2)
  CALL ChargeTime( pstats%mmb%tm, pstats%mmb%t2, pstats%mmb%t1, pstats%mmb%cnt )

  IF(KPDBG) WRITE(OFU,*) 'SlaveDGEMM receiving C'; CALL FL(OFU)
  CALL BSystemClock(pstats%mmc%t1)
  CALL SlaveReceiveMatrix( nrows, ncols, subC, subnrows, subncols )
  CALL BSystemClock(pstats%mmc%t2)
  CALL ChargeTime( pstats%mmc%tm, pstats%mmc%t2, pstats%mmc%t1, pstats%mmc%cnt )

  IF(KPDBG) WRITE(OFU,*) 'SlaveDGEMM receiving alpha'; CALL FL(OFU)
  CALL BSystemClock(pstats%mmalpha%t1)
  CALL SlaveReceiveValue( alpha )
  CALL BSystemClock(pstats%mmalpha%t2)
  CALL ChargeTime( pstats%mmalpha%tm, pstats%mmalpha%t2, pstats%mmalpha%t1, pstats%mmalpha%cnt )

  IF(KPDBG) WRITE(OFU,*) 'SlaveDGEMM receiving beta'; CALL FL(OFU)
  CALL BSystemClock(pstats%mmbeta%t1)
  CALL SlaveReceiveValue( beta )
  CALL BSystemClock(pstats%mmbeta%t2)
  CALL ChargeTime( pstats%mmbeta%tm, pstats%mmbeta%t2, pstats%mmbeta%t1, pstats%mmbeta%cnt )

  IF(KPDBG) WRITE(OFU,*) 'SlaveDGEMM invoking PDGEMM'; CALL FL(OFU)
  CALL BSystemClock(pstats%comp%t1)
  CALL PDGEMM( 'N', 'N', M, M, M, alpha, &
               subA, 1, 1, descA, &
               subB, 1, 1, descB, &
               beta, &
               subC, 1, 1, descC )
  CALL BSystemClock(pstats%comp%t2)
  CALL ChargeTime( pstats%comp%tm, pstats%comp%t2, pstats%comp%t1, pstats%comp%cnt )
  CALL ChargeTime( pstats%pmm%tm, pstats%comp%t2, pstats%comp%t1, pstats%pmm%cnt )
  IF(KPDBG) WRITE(OFU,*) 'SlaveDGEMM done PDGEMM'; CALL FL(OFU)

  IF(KPDBG) WRITE(OFU,*) 'SlaveDGEMM sending result matrix to master'; CALL FL(OFU)
  CALL BSystemClock(pstats%mmrc%t1)
  CALL SlaveSendMatrix( nrows, ncols, subC, subnrows, subncols )
  CALL BSystemClock(pstats%mmrc%t2)
  CALL ChargeTime( pstats%mmrc%tm, pstats%mmrc%t2, pstats%mmrc%t1, pstats%mmrc%cnt )
  IF(KPDBG) WRITE(OFU,*) 'SlaveDGEMM sent result matrix to master'; CALL FL(OFU)

  IF(KPDBG) WRITE(OFU,*) 'SlaveDGEMM deallocating subABC'; CALL FL(OFU)
  DEALLOCATE( subA )
  DEALLOCATE( subB )
  DEALLOCATE( subC )
  IF(KPDBG) WRITE(OFU,*) 'SlaveDGEMM deallocated subABC'; CALL FL(OFU)
#endif

  CALL BSystemClock(pstats%mm%t2)
  CALL ChargeTime( pstats%mm%tm, pstats%mm%t2, pstats%mm%t1, pstats%mm%cnt )

END SUBROUTINE SlaveDGEMM

!-------------------------------------------------------------------------------
!>
!! Encapsulates BLAS' DGEMV functionality
!<
!-------------------------------------------------------------------------------
SUBROUTINE PLBDGEMV( alpha, A, x, beta, y )
  !----------------------------------------------
  ! Formal arguments
  !----------------------------------------------
  REAL(dp), INTENT(IN)    :: alpha, beta
  REAL(dp), INTENT(IN)    :: A(:,:) !<Matrix of size MxM
  REAL(dp), INTENT(IN)    :: x(:) !<Vector of size Mx1
  REAL(dp), INTENT(INOUT) :: y(:) !<Vector of size Mx1

  !----------------------------------------------
  !Just do this locally (not worth asking slaves' help)
  CALL DGEMV( 'N', M, M, alpha, A, M, x, 1, beta, y, 1 )
END SUBROUTINE PLBDGEMV

!-------------------------------------------------------------------------------
!>
!! Encapsulates BLAS/LAPACK's DGETRF functionality
!<
!-------------------------------------------------------------------------------
SUBROUTINE PLBDGETRF( A, piv, info )
  !----------------------------------------------
  ! Formal arguments
  !----------------------------------------------
  REAL(dp), INTENT(INOUT) :: A(:,:) !<Matrix of size MxM
  INTEGER,  INTENT(OUT)   :: piv(:) !<Mx1
  INTEGER,  INTENT(INOUT) :: info
  !----------------------------------------------
  ! Local variables
  !----------------------------------------------
  REAL(dp), ALLOCATABLE :: subA(:)
  INTEGER, ALLOCATABLE :: subpiv(:)
  INTEGER :: descA(9)
  INTEGER :: lldA
  INTEGER :: nrows, ncols, subnrows, subncols
  INTEGER :: bszr, bszc, pnr, pnc, pi, pj
  INTEGER :: ctxt
  LOGICAL :: useblas
#if defined(MPI_OPT)
  INTEGER, EXTERNAL :: NUMROC
#endif

  IF(KPDBG) WRITE(OFU,*) 'MasterGETRF started'; CALL FL(OFU)

  !----------------------------------------------
  useblas = ( doblasonly ) .OR. &
            ( pblas%msmap%nslaves .EQ. 1) .OR. &
            ( M .LE. blacs%pgrid%blockszrows ) .OR. &
            ( M .LE. blacs%pgrid%blockszcols )

  !-----------------------------------------------------------------------------
  IF ( useblas ) THEN
    !IF(KPDBG) WRITE(OFU,*) 'BLAS DGETRF only (not using PBLAS)' ; CALL FL(OFU)
    IF(KPDBG) WRITE(OFU,*) 'BLAS DGETRF only (not using PBLAS) with M=',M ; CALL FL(OFU)
    CALL DGETRF( M, M, A, M, piv, info )
  ELSE
#if defined(MPI_OPT)
    CALL BSystemClock(pstats%trf%t1)
    nrows = M; ncols = M
    bszr = blacs%pgrid%blockszrows; bszc = blacs%pgrid%blockszcols
    pnr = blacs%pgrid%nrows; pnc = blacs%pgrid%ncols
    pi = blacs%pgrid%myrow; pj = blacs%pgrid%mycol
    subnrows = NUMROC( nrows, bszr, pi, 0, pnr )
    subncols = NUMROC( ncols, bszc, pj, 0, pnc )

    ctxt = blacs%levelcontext
    lldA = subnrows; IF( lldA .LT. 1 ) lldA = 1
    CALL DESCINIT(descA, nrows, ncols, bszr, bszc, 0, 0, ctxt, lldA, info )

    IF(KPDBG) WRITE(OFU,*) 'MasterDGETRF allocating subAPiv'; CALL FL(OFU)
    ALLOCATE( subA(1:subnrows*subncols) )
    ALLOCATE( subpiv(1:subnrows+bszr) )
    IF(KPDBG) WRITE(OFU,*) 'MasterDGETRF allocated subAPiv'; CALL FL(OFU)

    IF(KPDBG) WRITE(OFU,*) 'MasterDGETRF sending OP_GETRF'; CALL FL(OFU)
    CALL MasterBcastNextOp( OP_DGETRF )
    IF(KPDBG) WRITE(OFU,*) 'MasterDGETRF sending A'; CALL FL(OFU)
    CALL MasterSendMatrix( A, nrows, ncols, subA, subnrows, subncols )
    IF(KPDBG) WRITE(OFU,*) 'MasterDGETRF invoking PDGETRF'; CALL FL(OFU)
    CALL BSystemClock(pstats%comp%t1)
    CALL PDGETRF( nrows, ncols, subA, 1, 1, descA, subpiv, info )
    CALL BSystemClock(pstats%comp%t2)
    CALL ChargeTime( pstats%comp%tm, pstats%comp%t2, pstats%comp%t1, pstats%comp%cnt )
    CALL ChargeTime( pstats%ptrf%tm, pstats%comp%t2, pstats%comp%t1, pstats%ptrf%cnt )
    IF(KPDBG) WRITE(OFU,*) 'MasterDGETRF done PDGETRF'; CALL FL(OFU)

    IF(KPDBG) WRITE(OFU,*) 'MasterDGETRF receiving slave submatrices'; CALL FL(OFU)
    CALL MasterRecvMatrix( A, nrows, ncols, subA, subnrows, subncols )
    IF(KPDBG) WRITE(OFU,*) 'MasterDGETRF received slave submatrices'; CALL FL(OFU)
    IF(KPDBG) WRITE(OFU,*) 'MasterDGETRF receiving slave vectors'; CALL FL(OFU)
    CALL MasterRecvVector( piv, nrows, subpiv, subnrows )
    IF(KPDBG) WRITE(OFU,*) 'MasterDGETRF received slave vectors'; CALL FL(OFU)

    IF(KPDBG) WRITE(OFU,*) 'MasterDGETRF deallocating subAPiv'; CALL FL(OFU)
    DEALLOCATE( subpiv )
    DEALLOCATE( subA )
    IF(KPDBG) WRITE(OFU,*) 'MasterDGETRF deallocated subAPiv'; CALL FL(OFU)
    CALL BSystemClock(pstats%trf%t2)
    CALL ChargeTime( pstats%trf%tm, pstats%trf%t2, pstats%trf%t1, pstats%trf%cnt )
#endif
  END IF

  IF(KPDBG) WRITE(OFU,*) 'MasterDGETRF done'; CALL FL(OFU)
END SUBROUTINE PLBDGETRF

!-------------------------------------------------------------------------------
!>
!! DGETRF support from slave
!<
!-------------------------------------------------------------------------------
SUBROUTINE SlaveDGETRF()
  !----------------------------------------------
  ! Local variables
  !----------------------------------------------
  REAL(dp), ALLOCATABLE :: subA(:)
  INTEGER, ALLOCATABLE :: subpiv(:)
  INTEGER :: descA(9)
  INTEGER :: lldA
  INTEGER :: nrows, ncols, subnrows, subncols
  INTEGER :: bszr, bszc, pnr, pnc, pi, pj
  INTEGER :: ctxt, info
#if defined(MPI_OPT)
  INTEGER, EXTERNAL :: NUMROC
#endif

  CALL BSystemClock(pstats%trf%t1)

#if defined(MPI_OPT)
  nrows = M; ncols = M
  bszr = blacs%pgrid%blockszrows; bszc = blacs%pgrid%blockszcols
  pnr = blacs%pgrid%nrows; pnc = blacs%pgrid%ncols
  pi = blacs%pgrid%myrow; pj = blacs%pgrid%mycol
  subnrows = NUMROC( nrows, bszr, pi, 0, pnr )
  subncols = NUMROC( ncols, bszc, pj, 0, pnc )

  ctxt = blacs%levelcontext
  lldA = subnrows; IF( lldA .LT. 1 ) lldA = 1
  CALL DESCINIT(descA, nrows, ncols, bszr, bszc, 0, 0, ctxt, lldA, info )

  IF(KPDBG) WRITE(OFU,*) 'SlaveDGETRF allocating subAPiv'; CALL FL(OFU)
  ALLOCATE( subA(1:subnrows*subncols) )
  ALLOCATE( subpiv(1:subnrows+bszr) )
  IF(KPDBG) WRITE(OFU,*) 'SlaveDGETRF allocated subAPiv'; CALL FL(OFU)

  IF(KPDBG) WRITE(OFU,*) 'SlaveDGETRF receiving A submatrix'; CALL FL(OFU)
  CALL SlaveReceiveMatrix( nrows, ncols, subA, subnrows, subncols )
  IF(KPDBG) WRITE(OFU,*) 'SlaveDGETRF received A submatrix'; CALL FL(OFU)

  IF(KPDBG) WRITE(OFU,*) 'MasterDGETRF invoking PDGETRF'; CALL FL(OFU)
  CALL BSystemClock(pstats%comp%t1)
  CALL PDGETRF( nrows, ncols, subA, 1, 1, descA, subpiv, info )
  CALL BSystemClock(pstats%comp%t2)
  CALL ChargeTime( pstats%comp%tm, pstats%comp%t2, pstats%comp%t1, pstats%comp%cnt )
  CALL ChargeTime( pstats%ptrf%tm, pstats%comp%t2, pstats%comp%t1, pstats%ptrf%cnt )
  IF(KPDBG) WRITE(OFU,*) 'MasterDGETRF done PDGETRF'; CALL FL(OFU)

  IF(KPDBG) WRITE(OFU,*) 'SlaveDGETRF sending result matrix to master'; CALL FL(OFU)
  CALL SlaveSendMatrix( nrows, ncols, subA, subnrows, subncols )
  IF(KPDBG) WRITE(OFU,*) 'SlaveDGETRF sent result matrix to master'; CALL FL(OFU)
  CALL SlaveSendVector( nrows, subpiv, subnrows )
  IF(KPDBG) WRITE(OFU,*) 'SlaveDGETRF sent result vector to master'; CALL FL(OFU)

  IF(KPDBG) WRITE(OFU,*) 'SlaveDGETRF deallocating subAPiv'; CALL FL(OFU)
  DEALLOCATE( subpiv )
  DEALLOCATE( subA )
  IF(KPDBG) WRITE(OFU,*) 'SlaveDGETRF deallocated subAPiv'; CALL FL(OFU)

#endif

  CALL BSystemClock(pstats%trf%t2)
  CALL ChargeTime( pstats%trf%tm, pstats%trf%t2, pstats%trf%t1, pstats%trf%cnt )

END SUBROUTINE SlaveDGETRF

!-------------------------------------------------------------------------------
!>
!! Encapsulates BLAS/LAPACK's DGETRS functionality
!<
!-------------------------------------------------------------------------------
SUBROUTINE PLBDGETRS( nrhs, A, piv, B, info )
  !----------------------------------------------
  ! Formal arguments
  !----------------------------------------------
  INTEGER,  INTENT(IN) :: nrhs
  REAL(dp), INTENT(IN) :: A(:,:) !<Matrix of size MxM
  INTEGER,  INTENT(IN) :: piv(:) !<Mx1
  REAL(dp), INTENT(INOUT) :: B(:,:) !<Matrix of size Mxnrhs
  INTEGER :: info

  !----------------------------------------------
  !Just do this locally (not worth asking slaves' help)
  IF(KPDBG) WRITE(OFU,*) 'BLAS DGETRS only (not using PBLAS)'; CALL FL(OFU)
  CALL DGETRS( 'N', M, nrhs, A, M, piv, B, M, info )
END SUBROUTINE PLBDGETRS

!-------------------------------------------------------------------------------
!>
!! DGETRS support from slave
!<
!-------------------------------------------------------------------------------
SUBROUTINE SlaveDGETRS()

  IF(KPDBG) WRITE(OFU,*) 'SlaveDGETRS not implemented'; CALL FL(OFU)
  CALL ASSERT(.FALSE.,'blocktri #12') !To be completed

END SUBROUTINE SlaveDGETRS
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
!>
!! To be invoked, before solve, by user
!<
!-------------------------------------------------------------------------------
SUBROUTINE Initialize( inN, inM )
  !-----------------------------------------------------------------------------
  ! Formal arguments
  !-----------------------------------------------------------------------------
  INTEGER, INTENT(IN) :: inN !<Num of row blocks in input block tri-diag matrix
  INTEGER, INTENT(IN) :: inM !<Size of each square sub-matrix block

  !-----------------------------------------------------------------------------
  ! Local variables
  !-----------------------------------------------------------------------------
  INTEGER :: nrowsperrank !<Avg number of rows (+/-1) mapped to each processor
  INTEGER :: nspillrows !<Number of rows left over to smear after even mapping
  INTEGER :: mpierr !<MPI error code
  INTEGER :: globrow !<Loop variable
  INTEGER :: globrowoff !<Temporary variable
  LOGICAL :: flag   !< false if mpi NOT presently initialized
  REAL :: logNbase2 !<Temporary variable

  !-----------------------------------------------------------------------------
  LBCYLIC_INIT=.TRUE.
#if defined(MPI_OPT)
  CALL MPI_Initialized (flag, mpierr)
  IF (.not.flag) THEN
    CALL MPI_Init( mpierr )
  END IF
  use_mpiwtime = .TRUE.
  CALL MPI_Comm_rank( SIESTA_COMM, rank, mpierr )
  CALL MPI_Comm_size( SIESTA_COMM, nranks, mpierr)
#endif

  N = inN
  M = inM
  IF (N.EQ.0 .OR. M.EQ.0) return

  !-----------------------------------------------------------------------------
  ! Set output file unit, so we get output of each rank in its own file
  ! Extremely useful for debugging -- makes it so much easier to see progress!
  !-----------------------------------------------------------------------------
  OFU = rank+1000
  P = nranks
  IF ( P .GT. N ) THEN
    IF(KPDBG) WRITE(OFU,*) 'Detected extra ', P-N, ' ranks'; CALL FL(OFU)
    P = N
  END IF

  !-----------------------------------------------------------------------------
  !Set problem output file
  !-----------------------------------------------------------------------------
  PFU = OFU+nranks

  !-----------------------------------------------------------------------------
  ! Check to see if debugging output is to be written : SKS on Nov 9, 2011
  !-----------------------------------------------------------------------------
  KPDBG=.FALSE.; kenvvar='KPDBG'
  CALL GETENV(kenvvar,kenvval)
  IF (kenvval.NE.'') THEN
    IF (kenvval.EQ.'TRUE') THEN
      KPDBG=.TRUE.
    ELSE
      KPDBG=.FALSE.
    END IF
  END IF

  !-----------------------------------------------------------------------------
  IF(KPDBG) WRITE(OFU,*) 'Rank ', rank, ' NRanks ', nranks; CALL FL(OFU)
  IF(KPDBG) WRITE(OFU,*); CALL FL(OFU)
  IF(KPDBG) WRITE(OFU,*) '------ Initialization start ------'; CALL FL(OFU)
  IF(KPDBG) WRITE(OFU,*); CALL FL(OFU)

  !-----------------------------------------------------------------------------
  IF ( N .LT. 1 ) THEN
    IF(KPDBG) WRITE(OFU,*) 'Bad N', N; CALL FL(OFU)
    CALL ASSERT(.FALSE.,'blocktri #13')
  END IF

  !-----------------------------------------------------------------------------
  !Determine where we start and end in row list at level 1 on this processor
  !This processor's global row numbers span [startglobrow,endglobrow] inclusive
  !-----------------------------------------------------------------------------
  nrowsperrank = N/P !Number of rows evenly split across processors
  nspillrows = MOD(N, P) !Some left over after even split
  IF ( rank .LT. nspillrows ) THEN !The first few ranks get one extra row each
    startglobrow = rank*nrowsperrank + rank + 1
    endglobrow = startglobrow + nrowsperrank
  ELSE IF ( rank .LT. P ) THEN !The other active ranks
    startglobrow = rank*nrowsperrank + nspillrows + 1
    endglobrow = startglobrow + nrowsperrank - 1
  ELSE !The rest (unnecessary, inactive) ranks
    startglobrow = -1
    endglobrow = -2
  END IF
  
  !-----------------------------------------------------------------------------
  !Allocate array to store unmodified diagonal elements: SKS
  !-----------------------------------------------------------------------------
  IF (.NOT.ALLOCATED(OrigDiagElement)) ALLOCATE (OrigDiagElement((endglobrow-startglobrow+1)*M))
  
  !-----------------------------------------------------------------------------
  !No. of levels is 1+ceil(lg(N)); e.g., 6 for N=23
  !-----------------------------------------------------------------------------
  logNbase2 = LOG(REAL(N))/LOG(2.0)
  nlevels = 1 + INT(logNbase2)
  IF ( REAL(INT(logNbase2)) .NE. logNbase2 ) THEN
    nlevels = nlevels + 1
  END IF
  IF(KPDBG) WRITE(OFU,*) 'NLEVELS=', nlevels; CALL FL(OFU)

  !---------------------------------------------------------------------------
  matdirtied = .TRUE.
  rhsdirtied = .TRUE.

  !---------------------------------------------------------------------------
  !Not sure if barriers are needed between levels, but they can be optionally
  !turned on by setting usebarriers to TRUE
  !---------------------------------------------------------------------------
  usebarriers = .FALSE.

  !---------------------------------------------------------------------------
  writeproblemfile = .FALSE.
  writesolution = .FALSE.

  !---------------------------------------------------------------------------
  membytes = 0
  dpsz = 8 !8 bytes for double precision
  intsz = 4 !4 bytes for a single integer
  ptrsz = 8 !Assuming 64-bit addresses in the worst case
  ONE = 1
  ZERO = 0

  tottime = 0
  totcommtime = 0
  totinvtime = 0
  totinvcount = 0
  totmatmultime = 0
  totmatmulcount = 0
  totmatsoltime = 0
  totmatsolcount = 0
  CALL BClockInit()

  !---------------------------------------------------------------------------
  !The +2 in nlocrows+2 is for storing incoming values from neighbors, if any
  !The zero'th element stores the incoming element from top neighbor, the
  !last element stores the incoming element from the bottom neighbor.
  !---------------------------------------------------------------------------
  ! SKS 
  IF (.NOT.ALLOCATED(lelement)) ALLOCATE( lelement(nlevels, 0:endglobrow-startglobrow+2) )
  !ALLOCATE( lelement(nlevels, 0:endglobrow-startglobrow+2) )
  CALL ChargeMemory( nlevels*(endglobrow-startglobrow+3)*5*ptrsz )

  ! SKS
  IF (.NOT.ALLOCATED(selement)) ALLOCATE( selement(N) ) !Have a place for all x; only some will be allocated
  ! ALLOCATE( selement(N) ) !Have a place for all x; only some will be allocated
  CALL ChargeMemory( N*ptrsz )

  !---------------------------------------------------------------------------
  !Allocate (L,D,U,b,pivot) at level 1
  !---------------------------------------------------------------------------
  DO globrow = startglobrow, endglobrow, 1
    globrowoff = globrow-startglobrow+1
    !SKS start
    IF(.NOT.ALLOCATED(lelement(1,globrowoff)%L)) ALLOCATE( lelement(1, globrowoff)%L(M,M) )
    IF(.NOT.ALLOCATED(lelement(1,globrowoff)%D)) ALLOCATE( lelement(1, globrowoff)%D(M,M) )
    IF(.NOT.ALLOCATED(lelement(1,globrowoff)%U)) ALLOCATE( lelement(1, globrowoff)%U(M,M) )
    IF(.NOT.ALLOCATED(lelement(1,globrowoff)%b)) ALLOCATE( lelement(1, globrowoff)%b(M,1) )
    IF(.NOT.ALLOCATED(lelement(1,globrowoff)%pivot)) ALLOCATE( lelement(1, globrowoff)%pivot(M) )
    !SKS end
    !ALLOCATE( lelement(1, globrowoff)%L(M,M) )
    !ALLOCATE( lelement(1, globrowoff)%D(M,M) )
    !ALLOCATE( lelement(1, globrowoff)%U(M,M) )
    !ALLOCATE( lelement(1, globrowoff)%b(M,1) )
    !ALLOCATE( lelement(1, globrowoff)%pivot(M) )
    CALL ChargeMemory( (3*M*M+M)*dpsz + M*intsz )
    lelement(1, globrowoff)%L = 0
    lelement(1, globrowoff)%D = 0
    lelement(1, globrowoff)%U = 0
    lelement(1, globrowoff)%b = 0
    lelement(1, globrowoff)%pivot = 0
  END DO

  !---------------------------------------------------------------------------
  !Allocate (L,D,U,b) for saving the original problem
  !---------------------------------------------------------------------------
  !SKS
  IF (.NOT.ALLOCATED(orig)) ALLOCATE( orig(endglobrow-startglobrow+1) )
  ! ALLOCATE( orig(endglobrow-startglobrow+1) )
  DO globrow = startglobrow, endglobrow, 1
    globrowoff = globrow-startglobrow+1
    !SKS start
    IF(.NOT.ALLOCATED(orig(globrowoff)%L)) ALLOCATE( orig(globrowoff)%L(M,M) )
    IF(.NOT.ALLOCATED(orig(globrowoff)%D)) ALLOCATE( orig(globrowoff)%D(M,M) )
    IF(.NOT.ALLOCATED(orig(globrowoff)%U)) ALLOCATE( orig(globrowoff)%U(M,M) )
    IF(.NOT.ALLOCATED(orig(globrowoff)%b)) ALLOCATE( orig(globrowoff)%b(M,1) )
    !SKS end
    !ALLOCATE( orig(globrowoff)%L(M,M) )
    !ALLOCATE( orig(globrowoff)%D(M,M) )
    !ALLOCATE( orig(globrowoff)%U(M,M) )
    !ALLOCATE( orig(globrowoff)%b(M,1) )
    !-- DONT CHARGE MEMORY COST FOR THIS!!! --
    orig(globrowoff)%L = 0
    orig(globrowoff)%D = 0
    orig(globrowoff)%U = 0
    orig(globrowoff)%b = 0
  END DO

  !-----------------------------------------------------------------------------
  CALL PLBInitialize()

  !-----------------------------------------------------------------------------
  IF(KPDBG) WRITE(OFU,*) '   P=',P,'     N=',N,'    M=',M; CALL FL(OFU)
  IF(KPDBG) WRITE(OFU,*) 'SROW=',startglobrow,'  EROW=', endglobrow; CALL FL(OFU)
  IF(KPDBG) WRITE(OFU,*) '------ Initialization end ------'; CALL FL(OFU)
  IF(KPDBG) WRITE(OFU,*); CALL FL(OFU)
END SUBROUTINE Initialize

!-------------------------------------------------------------------------------
!SPH - Add this since now rank, nranks is PRIVATE to avoid conflict with nscalingtools
!      Must call Initialize first
SUBROUTINE GetRanks(irank, inranks)
INTEGER, INTENT(OUT) :: irank, inranks
irank = rank;  inranks = nranks
END SUBROUTINE GetRanks

!-------------------------------------------------------------------------------
!>
!! To be invoked, before solve, after initialize, by user, to set up the matrix
!<
!-------------------------------------------------------------------------------
SUBROUTINE SetMatrixRow(globrow, L, D, U)
  !-----------------------------------------------------------------------------
  ! Formal arguments
  !-----------------------------------------------------------------------------
  INTEGER, INTENT(IN) :: globrow !<Original/global block-row num in [1..N]
  REAL(dp), INTENT(IN) :: L(:,:), D(:,:), U(:,:)!<Vals in 1st L(or Nth U) == 0.0
  !-----------------------------------------------------------------------------
  ! Local variables
  !-----------------------------------------------------------------------------
  REAL(dp) :: val
  INTEGER :: i, j, globrowoff

  !-------------------------------------------
  ! Sanity checks on globrow
  !-------------------------------------------
  IF ((globrow .LT. 1) .OR. (globrow .GT. N)) THEN
     IF (KPDBG) THEN
        WRITE (OFU,*) 'SetMatrixRow: Bad input globrow ',globrow
        CALL FL(OFU)
     END IF
     CALL ASSERT(.FALSE.,'blocktri #14')
  END IF
  IF ((globrow .LT. startglobrow) .OR. (globrow .GT. endglobrow)) THEN
     IF (KPDBG) THEN
        WRITE (OFU,*) 'SetMatrixRow: Non-local globrow ',globrow
        CALL FL(OFU)
     END IF
     CALL ASSERT(.FALSE.,'blocktri #15')
  END IF

  globrowoff = globrow-startglobrow+1

  !-------------------------------------------
  ! Copy given blocks into allocated matrix
  !-------------------------------------------
  DO j = 1, M
     DO i = 1, M
        IF (globrow .EQ. 1) THEN
           val = 0
        ELSE
           val = L(i,j)
        END IF
        lelement(1, globrowoff)%L(i,j) = val

        val = D(i,j)
        lelement(1, globrowoff)%D(i,j) = val

        IF (globrow .EQ. N) THEN
           val = 0
        ELSE
           val = U(i,j)
        END IF
        lelement(1, globrowoff)%U(i,j) = val
     END DO
  END DO

  !-------------------------------------------
  !Save a copy of this original problem
  !-------------------------------------------
  orig(globrowoff)%L = lelement(1,globrowoff)%L
  orig(globrowoff)%D = lelement(1,globrowoff)%D
  orig(globrowoff)%U = lelement(1,globrowoff)%U

  !-------------------------------------------
  matdirtied = .TRUE.

END SUBROUTINE SetMatrixRow

!-------------------------------------------------------------------------------
!>
!! To be invoked, before solve, after initialize, by user, to set up the matrix
!<
!-------------------------------------------------------------------------------
SUBROUTINE SetMatrixRowL(globrow, L)
  !-----------------------------------------------------------------------------
  ! Formal arguments
  !-----------------------------------------------------------------------------
  INTEGER, INTENT(IN) :: globrow !<Original/global block-row num in [1..N]
  REAL(dp), INTENT(IN) :: L(:,:) !<Vals in 1st L == 0.0
  !-----------------------------------------------------------------------------
  ! Local variables
  !-----------------------------------------------------------------------------
  REAL(dp) :: val
  INTEGER :: i, j, globrowoff

  !-------------------------------------------
  ! Sanity checks on globrow
  !-------------------------------------------
  IF ((globrow .LT. 1) .OR. (globrow .GT. N)) THEN
     IF (KPDBG) THEN
        WRITE(OFU,*) 'SetMatrixRowL: Bad input globrow ',globrow
        CALL FL(OFU)
     END IF
     CALL ASSERT(.FALSE.,'blocktri #16')
  END IF
  IF ((globrow .LT. startglobrow) .OR. (globrow .GT. endglobrow)) THEN
     IF (KPDBG) THEN
        WRITE(OFU,*) 'SetMatrixRowL: Non-local globrow ',globrow
        CALL FL(OFU)
     END IF
     CALL ASSERT(.FALSE.,'blocktri #17')
  END IF

  globrowoff = globrow - startglobrow + 1

  !-------------------------------------------
  ! Copy given L block into allocated matrix
  !-------------------------------------------
  DO j = 1, M
     DO i = 1, M
        IF (globrow .EQ. 1) THEN
           val = 0
        ELSE
           val = L(i,j)
        END IF
        lelement(1, globrowoff)%L(i,j) = val
     END DO
  END DO

  !-------------------------------------------
  !Save a copy of this original problem
  !-------------------------------------------
  orig(globrowoff)%L = lelement(1,globrowoff)%L

  !-------------------------------------------
  matdirtied = .TRUE.

END SUBROUTINE SetMatrixRowL

!-------------------------------------------------------------------------------
!>
!! To be invoked, before solve, after initialize, by user, to set up the matrix
!<
!-------------------------------------------------------------------------------
SUBROUTINE SetMatrixRowD(globrow, D)
  !-----------------------------------------------------------------------------
  ! Formal arguments
  !-----------------------------------------------------------------------------
  INTEGER, INTENT(IN) :: globrow !<Original/global block-row num in [1..N]
  REAL(dp), INTENT(IN) :: D(:,:)
  !-----------------------------------------------------------------------------
  ! Local variables
  !-----------------------------------------------------------------------------
  REAL(dp) :: val
  INTEGER :: i, j, globrowoff

  !-------------------------------------------
  ! Sanity checks on globrow
  !-------------------------------------------
  IF ((globrow .LT. 1) .OR. (globrow .GT. N)) THEN
     IF (KPDBG) THEN
        WRITE(OFU,*) 'SetMatrixRowD: Bad input globrow ',globrow
        CALL FL(OFU)
     END IF
     CALL ASSERT(.FALSE.,'blocktri #18')
  END IF
  IF ((globrow .LT. startglobrow) .OR. (globrow .GT. endglobrow)) THEN
     IF (KPDBG) THEN
        WRITE(OFU,*) 'SetMatrixRowD: Non-local globrow ',globrow
        CALL FL(OFU)
     END IF
     CALL ASSERT(.FALSE.,'blocktri #19')
  END IF

  globrowoff = globrow - startglobrow + 1

  !-------------------------------------------
  ! Copy given D block into allocated matrix
  !-------------------------------------------
  DO j = 1, M
     DO i = 1, M
        val = D(i,j)
        lelement(1, globrowoff)%D(i,j) = val
     END DO
  END DO

  !-------------------------------------------
  !Save a copy of this original problem
  !-------------------------------------------
  orig(globrowoff)%D = lelement(1,globrowoff)%D

  !-------------------------------------------
  matdirtied = .TRUE.

END SUBROUTINE SetMatrixRowD

!-------------------------------------------------------------------------------
!>
!! To be invoked, before solve, after initialize, by user, to set up the matrix
!<
!-------------------------------------------------------------------------------
SUBROUTINE SetMatrixRowU(globrow, U)
  !-----------------------------------------------------------------------------
  ! Formal arguments
  !-----------------------------------------------------------------------------
  INTEGER, INTENT(IN) :: globrow !<Original/global block-row num in [1..N]
  REAL(dp), INTENT(IN) :: U(:,:)!<Vals in Nth U == 0.0
  !-----------------------------------------------------------------------------
  ! Local variables
  !-----------------------------------------------------------------------------
  REAL(dp) :: val
  INTEGER :: i, j, globrowoff

  !-------------------------------------------
  ! Sanity checks on globrow
  !-------------------------------------------
  IF ((globrow .LT. 1) .OR. (globrow .GT. N)) THEN
     IF (KPDBG) THEN
        WRITE(OFU,*) 'SetMatrixRowU: Bad input globrow ',globrow
        CALL FL(OFU)
     END IF
     CALL ASSERT(.FALSE.,'blocktri #20')
  END IF
  IF ((globrow .LT. startglobrow) .OR. (globrow .GT. endglobrow)) THEN
     IF (KPDBG) THEN
        WRITE(OFU,*) 'SetMatrixRowU: Non-local globrow ',globrow
        CALL FL(OFU)
     END IF
     CALL ASSERT(.FALSE.,'blocktri #21')
  END IF

  globrowoff = globrow - startglobrow + 1

  !-------------------------------------------
  ! Copy given U block into allocated matrix
  !-------------------------------------------
  DO j = 1, M
     DO i = 1, M
        IF (globrow .EQ. N) THEN
           val = 0
        ELSE
           val = U(i,j)
        END IF
        lelement(1, globrowoff)%U(i,j) = val
     END DO
  END DO

  !-------------------------------------------
  !Save a copy of this original problem
  !-------------------------------------------
  orig(globrowoff)%U = lelement(1,globrowoff)%U

  !-------------------------------------------
  matdirtied = .TRUE.

END SUBROUTINE SetMatrixRowU

!-------------------------------------------------------------------------------
!>
!! To be invoked, before solve, after initialize, by user, to set up the matrix
!<
!-------------------------------------------------------------------------------
SUBROUTINE SetMatrixRowColL(globrow, Lj, j)
  !-----------------------------------------------------------------------------
  ! Formal arguments
  !-----------------------------------------------------------------------------
  INTEGER, INTENT(IN) :: globrow !<Original/global block-row num in [1..N]
  REAL(dp), INTENT(IN) :: Lj(:) !<j'th colum of L at globrow; 1st L is always 0
  INTEGER, INTENT(IN) :: j !<column number of L at globrow that is being set
  !-----------------------------------------------------------------------------
  ! Local variables
  !-----------------------------------------------------------------------------
  REAL(dp) :: val
  INTEGER :: i, globrowoff

  !-------------------------------------------
  ! Sanity checks on globrow
  !-------------------------------------------
  IF ((globrow .LT. 1) .OR. (globrow .GT. N) ) THEN
     IF (KPDBG) THEN
        WRITE (OFU,*) 'SetMatrixRowColL: Bad input globrow ',globrow
        CALL FL(OFU)
     END IF
     CALL ASSERT(.FALSE.,'blocktri #22')
  END IF
  IF ((globrow .LT. startglobrow) .OR. (globrow .GT. endglobrow)) THEN
     IF (KPDBG) THEN
        WRITE (OFU,*) 'SetMatrixRowColL: Non-local globrow ',globrow
        CALL FL(OFU)
     END IF
     CALL ASSERT(.FALSE.,'blocktri #23')
  END IF
  IF (.NOT.((1 .LE. j) .AND. (j .LE. M))) THEN
     IF (KPDBG) THEN
        WRITE(OFU,*) 'SetMatrixRowColL: Bad j column ',j
        CALL FL(OFU)
     END IF
     CALL ASSERT(.FALSE.,'blocktri #24')
  END IF

  globrowoff = globrow - startglobrow + 1

  !-------------------------------------------
  ! Copy given L column into allocated matrix
  !-------------------------------------------
  DO i = 1, M
     IF (globrow .EQ. 1) THEN
        val = 0
     ELSE
        val = Lj(i)
     END IF
     lelement(1, globrowoff)%L(i,j) = val
  END DO

  !-------------------------------------------
  !Save a copy of this original problem
  !-------------------------------------------
  orig(globrowoff)%L(:,j) = lelement(1,globrowoff)%L(:,j)

  !-------------------------------------------
  matdirtied = .TRUE.

END SUBROUTINE SetMatrixRowColL

!-------------------------------------------------------------------------------
!>
!! To be invoked, before solve, after initialize, by user, to set up the matrix
!<
!-------------------------------------------------------------------------------
SUBROUTINE SetMatrixRowColD(globrow, Dj, j)
  !-----------------------------------------------------------------------------
  ! Formal arguments
  !-----------------------------------------------------------------------------
  INTEGER, INTENT(IN) :: globrow !<Original/global block-row num in [1..N]
  REAL(dp), INTENT(IN) :: Dj(:) !<j'th colum of D at globrow
  INTEGER, INTENT(IN) :: j !<column number of L at globrow that is being set
  !-----------------------------------------------------------------------------
  ! Local variables
  !-----------------------------------------------------------------------------
  REAL(dp) :: val
  INTEGER :: i, globrowoff

  !-------------------------------------------
  ! Sanity checks on globrow
  !-------------------------------------------
  IF ( (globrow .LT. 1) .OR. (globrow .GT. N) ) THEN
    IF(KPDBG) WRITE(OFU,*) 'SetMatrixRowColD: Bad input globrow ',globrow; CALL FL(OFU)
    CALL ASSERT(.FALSE.,'blocktri #25')
  END IF
  IF ( (globrow .LT. startglobrow) .OR. (globrow .GT. endglobrow) ) THEN
    IF(KPDBG) WRITE(OFU,*) 'SetMatrixRowColD: Non-local globrow ',globrow; CALL FL(OFU)
    CALL ASSERT(.FALSE.,'blocktri #26')
  END IF
  IF ( .NOT. ((1 .LE. j) .AND. (j .LE. M)) ) THEN
    IF(KPDBG) WRITE(OFU,*) 'SetMatrixRowColD: Bad j column ',j; CALL FL(OFU)
    CALL ASSERT(.FALSE.,'blocktri #27')
  END IF

  globrowoff = globrow-startglobrow+1

  !-------------------------------------------
  ! Copy given D column into allocated matrix
  !-------------------------------------------
  DO i = 1, M
    val = Dj(i)
    lelement(1, globrowoff)%D(i,j) = val
  END DO

  !-------------------------------------------
  !Save a copy of this original problem
  !-------------------------------------------
  orig(globrowoff)%D(:,j) = lelement(1,globrowoff)%D(:,j)

  !-------------------------------------------
  matdirtied = .TRUE.

END SUBROUTINE SetMatrixRowColD

!-------------------------------------------------------------------------------
!>
!! To be invoked, before solve, after initialize, by user, to set up the matrix
!<
!-------------------------------------------------------------------------------
SUBROUTINE SetMatrixRowColU(globrow, Uj, j)
USE descriptor_mod, ONLY: iam
  !-----------------------------------------------------------------------------
  ! Formal arguments
  !-----------------------------------------------------------------------------
  INTEGER, INTENT(IN)  :: globrow !<Original/global block-row num in [1..N]
  REAL(dp), INTENT(IN) :: Uj(:) !<j'th colum of U at globrow; Nth U is always 0
  INTEGER, INTENT(IN)  :: j !<column number of U at globrow that is being set
  !-----------------------------------------------------------------------------
  ! Local variables
  !-----------------------------------------------------------------------------
  REAL(dp) :: val
  INTEGER :: i, globrowoff

  !-------------------------------------------
  ! Sanity checks on globrow
  !-------------------------------------------
  IF ((globrow .LT. 1) .OR. (globrow .GT. N)) THEN
     IF (KPDBG) THEN
        WRITE (OFU,*) 'SetMatrixRowColU: Bad input globrow ', globrow
        CALL FL(OFU)
     END IF
     CALL ASSERT(.FALSE.,'blocktri #28')
  END IF
  IF ((globrow .LT. startglobrow) .OR. (globrow .GT. endglobrow)) THEN
     IF (KPDBG) THEN
        WRITE (OFU,*) 'SetMatrixRowColU: Non-local globrow ',globrow
        CALL FL(OFU)
     END IF
     CALL ASSERT(.FALSE.,'blocktri #29')
  END IF
  IF (.NOT.((1 .LE. j) .AND. (j .LE. M))) THEN
     IF (KPDBG) THEN
        WRITE (OFU,*) 'SetMatrixRowColL: Bad j column ',j
        CALL FL(OFU)
     END IF
     CALL ASSERT(.FALSE.,'blocktri #30')
  END IF

  globrowoff = globrow - startglobrow + 1

  !-------------------------------------------
  ! Copy given U column into allocated matrix
  !-------------------------------------------
  DO i = 1, M
     IF (globrow .EQ. N) THEN
        val = 0
     ELSE
        val = Uj(i)
     END IF
     lelement(1, globrowoff)%U(i,j) = val
  END DO

  !-------------------------------------------
  !Save a copy of this original problem
  !-------------------------------------------
  orig(globrowoff)%U(:,j) = lelement(1,globrowoff)%U(:,j)

  !-------------------------------------------
  matdirtied = .TRUE.

END SUBROUTINE SetMatrixRowColU

!-------------------------------------------------------------------------------
!>
!! To be invoked, before solve, after initialize, by user, to set up the matrix
!<
!-------------------------------------------------------------------------------
SUBROUTINE SetMatrixRHS(globrow, b)
  !-----------------------------------------------------------------------------
  ! Formal arguments
  !-----------------------------------------------------------------------------
  INTEGER, INTENT(IN) :: globrow !<Original/global block-row num in [1..N]
  REAL(dp), INTENT(IN) :: b(:) !<RHS column corr. to globrow
  !-----------------------------------------------------------------------------
  ! Local variables
  !-----------------------------------------------------------------------------
  INTEGER :: i, globrowoff

  !-------------------------------------------
  ! Sanity checks on globrow
  !-------------------------------------------
  IF ( (globrow .LT. 1) .OR. (globrow .GT. N) ) THEN
    IF(KPDBG) WRITE(OFU,*) 'SetMatrixRHS: Bad input globrow ',globrow; CALL FL(OFU)
    CALL ASSERT(.FALSE.,'blocktri #31')
  END IF
  IF ( (globrow .LT. startglobrow) .OR. (globrow .GT. endglobrow) ) THEN
    IF(KPDBG) WRITE(OFU,*) 'SetMatrixRHS: Non-local globrow ',globrow; CALL FL(OFU)
    CALL ASSERT(.FALSE.,'blocktri #32')
  END IF

  globrowoff = globrow - startglobrow + 1

  !-------------------------------------------
  ! Copy given values into allocated RHS
  !-------------------------------------------
  lelement(1, globrowoff)%b(1:M,1) = b(1:M)

  !-------------------------------------------
  !Save a copy of this original problem
  !-------------------------------------------
  orig(globrowoff)%b = lelement(1,globrowoff)%b

  !-------------------------------------------
  rhsdirtied = .TRUE.

END SUBROUTINE SetMatrixRHS

!-------------------------------------------------------------------------------
!>
!! Can be invoked after initialize, by user, to get a copy of the matrix
!<
!-------------------------------------------------------------------------------
SUBROUTINE GetMatrixRowColL( globrow, Lj, j )
  !-----------------------------------------------------------------------------
  ! Formal arguments
  !-----------------------------------------------------------------------------
  INTEGER, INTENT(IN) :: globrow !<Original/global block-row num in [1..N]
  REAL(dp), INTENT(OUT) :: Lj(:) !<j'th colum of L at globrow; 1st L is always 0
  INTEGER, INTENT(IN) :: j !<column number of L at globrow that is being set
  !-----------------------------------------------------------------------------
  ! Local variables
  !-----------------------------------------------------------------------------
  REAL(dp) :: val
  INTEGER :: i, globrowoff

  !-------------------------------------------
  ! Sanity checks on globrow
  !-------------------------------------------
  IF ( (globrow .LT. 1) .OR. (globrow .GT. N) ) THEN
    IF(KPDBG) WRITE(OFU,*) 'GetMatrixRowColL: Bad input globrow ',globrow; CALL FL(OFU)
    CALL ASSERT(.FALSE.,'blocktri #33')
  END IF
  IF ( (globrow .LT. startglobrow) .OR. (globrow .GT. endglobrow) ) THEN
    IF(KPDBG) WRITE(OFU,*) 'GetMatrixRowColL: Non-local globrow ',globrow; CALL FL(OFU)
    CALL ASSERT(.FALSE.,'blocktri #34')
  END IF
  IF ( .NOT. ((1 .LE. j) .AND. (j .LE. M)) ) THEN
    IF(KPDBG) WRITE(OFU,*) 'GetMatrixRowColL: Bad j column ',j; CALL FL(OFU)
    CALL ASSERT(.FALSE.,'blocktri #35')
  END IF

  globrowoff = globrow-startglobrow+1

  !-------------------------------------------
  ! Copy L from matrix
  !-------------------------------------------
  DO i = 1, M
    IF ( globrow .EQ. 1 ) THEN
      val = 0
    ELSE
      val = lelement(1, globrowoff)%L(i,j)
    END IF
    Lj(i) = val
  END DO

END SUBROUTINE GetMatrixRowColL

!-------------------------------------------------------------------------------
!>
!! Can be invoked after initialize, by user, to get a copy of the matrix
!<
!-------------------------------------------------------------------------------
SUBROUTINE GetMatrixRowColD( globrow, Dj, j )
  !-----------------------------------------------------------------------------
  ! Formal arguments
  !-----------------------------------------------------------------------------
  INTEGER, INTENT(IN) :: globrow !<Original/global block-row num in [1..N]
  REAL(dp), INTENT(OUT) :: Dj(:) !<j'th colum of D at globrow
  INTEGER, INTENT(IN) :: j !<column number of L at globrow that is being set
  !-----------------------------------------------------------------------------
  ! Local variables
  !-----------------------------------------------------------------------------
  REAL(dp) :: val
  INTEGER :: i, globrowoff

  !-------------------------------------------
  ! Sanity checks on globrow
  !-------------------------------------------
  IF ( (globrow .LT. 1) .OR. (globrow .GT. N) ) THEN
    IF(KPDBG) WRITE(OFU,*) 'GetMatrixRowColD: Bad input globrow ',globrow; CALL FL(OFU)
    CALL ASSERT(.FALSE.,'blocktri #36')
  END IF
  IF ( (globrow .LT. startglobrow) .OR. (globrow .GT. endglobrow) ) THEN
    IF(KPDBG) WRITE(OFU,*) 'GetMatrixRowColD: Non-local globrow ',globrow; CALL FL(OFU)
    CALL ASSERT(.FALSE.,'blocktri #37')
  END IF
  IF ( .NOT. ((1 .LE. j) .AND. (j .LE. M)) ) THEN
    IF(KPDBG) WRITE(OFU,*) 'GetMatrixRowColD: Bad j column ',j; CALL FL(OFU)
    CALL ASSERT(.FALSE.,'blocktri #38')
  END IF

  globrowoff = globrow-startglobrow+1

  !-------------------------------------------
  ! Copy D from matrix
  !-------------------------------------------
  DO i = 1, M
    val = lelement(1, globrowoff)%D(i,j)
    Dj(i) = val
  END DO

END SUBROUTINE GetMatrixRowColD

!-------------------------------------------------------------------------------
!>
!! Can be invoked after initialize, by user, to get a copy of the matrix
!<
!-------------------------------------------------------------------------------
SUBROUTINE GetMatrixRowColU( globrow, Uj, j )
  !-----------------------------------------------------------------------------
  ! Formal arguments
  !-----------------------------------------------------------------------------
  INTEGER, INTENT(IN) :: globrow !<Original/global block-row num in [1..N]
  REAL(dp), INTENT(OUT) :: Uj(:) !<j'th colum of L at globrow; Nth U is always 0
  INTEGER, INTENT(IN) :: j !<column number of U at globrow that is being set
  !-----------------------------------------------------------------------------
  ! Local variables
  !-----------------------------------------------------------------------------
  REAL(dp) :: val
  INTEGER :: i, globrowoff

  !-------------------------------------------
  ! Sanity checks on globrow
  !-------------------------------------------
  IF ( (globrow .LT. 1) .OR. (globrow .GT. N) ) THEN
    IF(KPDBG) WRITE(OFU,*) 'GetMatrixRowColU: Bad input globrow ',globrow; CALL FL(OFU)
    CALL ASSERT(.FALSE.,'blocktri #39')
  END IF
  IF ( (globrow .LT. startglobrow) .OR. (globrow .GT. endglobrow) ) THEN
    IF(KPDBG) WRITE(OFU,*) 'GetMatrixRowColU: Non-local globrow ',globrow; CALL FL(OFU)
    CALL ASSERT(.FALSE.,'blocktri #40')
  END IF
  IF ( .NOT. ((1 .LE. j) .AND. (j .LE. M)) ) THEN
    IF(KPDBG) WRITE(OFU,*) 'GetMatrixRowColL: Bad j column ',j; CALL FL(OFU)
    CALL ASSERT(.FALSE.,'blocktri #41')
  END IF

  globrowoff = globrow-startglobrow+1

  !-------------------------------------------
  ! Copy U from matrix
  !-------------------------------------------
  DO i = 1, M
    IF ( globrow .EQ. N ) THEN
      val = 0
    ELSE
      val = lelement(1, globrowoff)%U(i,j)
    END IF
    Uj(i) = val
  END DO

END SUBROUTINE GetMatrixRowColU

!-------------------------------------------------------------------------------
!>
!! Can be invoked after initialize, by user, to get a copy of the RHS
!<
!-------------------------------------------------------------------------------
SUBROUTINE GetMatrixRHS( globrow, b )
  !-----------------------------------------------------------------------------
  ! Formal arguments
  !-----------------------------------------------------------------------------
  INTEGER, INTENT(IN) :: globrow !<Original/global block-row num in [1..N]
  REAL(dp), INTENT(OUT) :: b(:) !<RHS column corr. to globrow
  !-----------------------------------------------------------------------------
  ! Local variables
  !-----------------------------------------------------------------------------
  REAL(dp) :: val
  INTEGER :: i, globrowoff

  !-------------------------------------------
  ! Sanity checks on globrow
  !-------------------------------------------
  IF ( (globrow .LT. 1) .OR. (globrow .GT. N) ) THEN
    IF(KPDBG) WRITE(OFU,*) 'GetMatrixRHS: Bad input globrow ',globrow; CALL FL(OFU)
    CALL ASSERT(.FALSE.,'blocktri #42')
  END IF
  IF ( (globrow .LT. startglobrow) .OR. (globrow .GT. endglobrow) ) THEN
    IF(KPDBG) WRITE(OFU,*) 'GetMatrixRHS: Non-local globrow ',globrow; CALL FL(OFU)
    CALL ASSERT(.FALSE.,'blocktri #43')
  END IF

  globrowoff = globrow-startglobrow+1

  !-------------------------------------------
  ! Copy RHS
  !-------------------------------------------
  DO i = 1, M
    val = lelement(1, globrowoff)%b(i,1)
    b(i) = val
  END DO

END SUBROUTINE GetMatrixRHS

!-------------------------------------------------------------------------------
!>
!! To be invoked, after solve, by user to get the solution vector of a local row
!<
!-------------------------------------------------------------------------------
SUBROUTINE GetSolutionVector( globrow, x )
  !-----------------------------------------------------------------------------
  ! Formal arguments
  !-----------------------------------------------------------------------------
  INTEGER, INTENT(IN) :: globrow !<Original/global block-row num in [1..N]
  REAL(dp), INTENT(OUT) :: x(:) !<Solution column corr. to globrow

  !-------------------------------------------
  ! Sanity checks on globrow
  !-------------------------------------------
  IF ( (globrow .LT. 1) .OR. (globrow .GT. N) ) THEN
    IF(KPDBG) WRITE(OFU,*) 'SetSolutionVector: Bad input globrow ',globrow; CALL FL(OFU)
    CALL ASSERT(.FALSE.,'blocktri #44')
  END IF
  IF ( (globrow .LT. startglobrow) .OR. (globrow .GT. endglobrow) ) THEN
    IF(KPDBG) WRITE(OFU,*) 'SetSolutionVector: Non-local globrow ',globrow; CALL FL(OFU)
    CALL ASSERT(.FALSE.,'blocktri #45')
  END IF

  !-------------------------------------------
  ! Copy solution into given vector
  !-------------------------------------------
  x = selement(globrow)%x

END SUBROUTINE GetSolutionVector

!-------------------------------------------------------------------------------
!>
!! To be invoked, before solve, after initialize, by user, to generate a test
!<
!-------------------------------------------------------------------------------
SUBROUTINE SetIdentityTestCase
  !-----------------------------------------------------------------------------
  ! Local variables
  !-----------------------------------------------------------------------------
  INTEGER :: i, j, globrow, globrowoff

  IF(KPDBG) WRITE(OFU,*) '------ Using Identity Test Case ------'; CALL FL(OFU)

  !-----------------------------------------
  !Initialize the arrays with identity
  !-----------------------------------------
  DO globrow = startglobrow, endglobrow, 1
    globrowoff = globrow-startglobrow+1
    DO i = 1, M
      DO j = 1, M
        lelement(1, globrowoff)%L(i,j) = ZERO
        IF ( i .EQ. j ) THEN
          lelement(1, globrowoff)%D(i,j) = ONE
        ELSE
          lelement(1, globrowoff)%D(i,j) = ZERO
        END IF
        lelement(1, globrowoff)%U(i,j) = ZERO
      END DO
      lelement(1, globrowoff)%b(i,1) = ONE
    END DO
    !-------------------------------------------
    !Save a copy of this original problem
    !-------------------------------------------
    orig(globrowoff)%L = lelement(1,globrowoff)%L
    orig(globrowoff)%D = lelement(1,globrowoff)%D
    orig(globrowoff)%U = lelement(1,globrowoff)%U
    orig(globrowoff)%b = lelement(1,globrowoff)%b
  END DO

  !-------------------------------------------
  matdirtied = .TRUE.
  rhsdirtied = .TRUE.

END SUBROUTINE SetIdentityTestCase

!-------------------------------------------------------------------------------
!>
!! To be invoked, before solve, after SetIdentityCase, by user,
!! to generate a test
!<
!-------------------------------------------------------------------------------
SUBROUTINE SetIdentityRHS
  !-----------------------------------------------------------------------------
  ! Local variables
  !-----------------------------------------------------------------------------
  INTEGER :: i, j, globrow, globrowoff

  IF(KPDBG) WRITE(OFU,*) '------ Using Identity Test Case RHS ------'; CALL FL(OFU)

  !-----------------------------------------
  !Initialize the arrays with identity
  !-----------------------------------------
  DO globrow = startglobrow, endglobrow, 1
    globrowoff = globrow-startglobrow+1
    DO i = 1, M
      lelement(1, globrowoff)%b(i,1) = ONE
    END DO
    !-------------------------------------------
    !Save a copy of this original RHS
    !-------------------------------------------
    orig(globrowoff)%b = lelement(1,globrowoff)%b
  END DO

  !-------------------------------------------
  rhsdirtied = .TRUE.

END SUBROUTINE SetIdentityRHS

!-------------------------------------------------------------------------------
!>
!! To be invoked, before solve, after initialize, by user, to generate a test
!<
!-------------------------------------------------------------------------------
SUBROUTINE SetRandomTestCase
  !-----------------------------------------------------------------------------
  ! Local variables
  !-----------------------------------------------------------------------------
  REAL(dp) :: randval, rmin, rmax
  INTEGER :: rngwidth
  INTEGER, ALLOCATABLE :: seeds(:)
  INTEGER :: i, j, globrow, globrowoff

  IF(KPDBG) WRITE(OFU,*) '------ Using Random Test Case ------'; CALL FL(OFU)

  !-----------------------------------------
  !Initialize random number seeds
  !-----------------------------------------
  CALL RANDOM_SEED( SIZE=rngwidth)
  ALLOCATE( seeds(rngwidth) )
  DO i = 1, rngwidth
    seeds(i) = i*(rank+100)*P
  END DO
  CALL RANDOM_SEED( PUT=seeds(1:rngwidth) )
  DEALLOCATE( seeds )

  !-----------------------------------------
  !Write the header information for problem
  !-----------------------------------------
  IF( writeproblemfile ) THEN !Write the dimensions N, NZ
    IF(KPDBG) WRITE(PFU,*) N*M; CALL FL(PFU) !<Print size of square matrix
    IF(KPDBG) WRITE(PFU,*) (3*N-2)*M*M; CALL FL(PFU) !<Print number of zeros
  END IF

  !-----------------------------------------
  !Initialize the arrays with random values
  !-----------------------------------------
  rmin = 0.01
  rmax = 1.00
  DO globrow = startglobrow, endglobrow, 1
    globrowoff = globrow-startglobrow+1
    DO i = 1, M
      DO j = 1, M
        IF ( globrow .EQ. 1 ) THEN
          randval = 0
        ELSE
          CALL RANDOM_NUMBER(randval); randval = rmin+(rmax-rmin)*randval
          IF( writeproblemfile ) THEN !<Print i, j, a[i,j]
            IF(KPDBG) WRITE(PFU,*) (globrow-1)*M+i, (globrow-2)*M+j, randval
          END IF
        END IF
        lelement(1, globrowoff)%L(i,j) = randval

        CALL RANDOM_NUMBER(randval); randval = rmin+(rmax-rmin)*randval
        lelement(1, globrowoff)%D(i,j) = randval
        IF( writeproblemfile ) THEN !<Print i, j, a[i,j]
          IF(KPDBG) WRITE(PFU,*) (globrow-1)*M+i, (globrow-1)*M+j, randval
        END IF

        IF ( globrow .EQ. N ) THEN
          randval = 0
        ELSE
          CALL RANDOM_NUMBER(randval); randval = rmin+(rmax-rmin)*randval
          IF( writeproblemfile ) THEN !<Print i, j, a[i,j]
            IF(KPDBG) WRITE(PFU,*) (globrow-1)*M+i, (globrow)*M+j, randval
          END IF
        END IF
        lelement(1, globrowoff)%U(i,j) = randval
      END DO
      CALL RANDOM_NUMBER(randval); randval = rmin+(rmax-rmin)*randval
      lelement(1, globrowoff)%b(i,1) = randval
    END DO

    !-------------------------------------------
    !Save a copy of this original problem
    !-------------------------------------------
    orig(globrowoff)%L = lelement(1,globrowoff)%L
    orig(globrowoff)%D = lelement(1,globrowoff)%D
    orig(globrowoff)%U = lelement(1,globrowoff)%U
    orig(globrowoff)%b = lelement(1,globrowoff)%b
  END DO

  IF (writeproblemfile) THEN
    DO globrow = startglobrow, endglobrow, 1
      globrowoff = globrow-startglobrow+1
      DO i = 1, M
        IF(KPDBG) WRITE(PFU,*) lelement(1,globrowoff)%b(i,1); CALL FL(PFU)
      END DO
    END DO
  END IF

  !-------------------------------------------
  matdirtied = .TRUE.
  rhsdirtied = .TRUE.

END SUBROUTINE SetRandomTestCase

!-------------------------------------------------------------------------------
!>
!! To be invoked, before solve, after SetRandomTestCase, by user, to
!! generate a an RHS
!<
!-------------------------------------------------------------------------------
SUBROUTINE SetRandomRHS( randseedoff )
  !-----------------------------------------------------------------------------
  ! Formal arguments
  !-----------------------------------------------------------------------------
  INTEGER, INTENT(IN) :: randseedoff
  !-----------------------------------------------------------------------------
  ! Local variables
  !-----------------------------------------------------------------------------
  REAL(dp) :: randval, rmin, rmax
  INTEGER :: rngwidth
  INTEGER, ALLOCATABLE :: seeds(:)
  INTEGER :: i, j, globrow, globrowoff

  IF(KPDBG) WRITE(OFU,*) '------ Using Random Test Case RHS ------'; CALL FL(OFU)

  !-----------------------------------------
  !Initialize random number seeds
  !-----------------------------------------
  CALL RANDOM_SEED( SIZE=rngwidth)
  ALLOCATE( seeds(rngwidth) )
  DO i = 1, rngwidth
    seeds(i) = i*(rank+100+34+randseedoff)*P
  END DO
  CALL RANDOM_SEED( PUT=seeds(1:rngwidth) )
  DEALLOCATE( seeds )

  !-----------------------------------------
  !Initialize the arrays with random values
  !-----------------------------------------
  rmin = 0.01
  rmax = 1.00
  DO globrow = startglobrow, endglobrow, 1
    globrowoff = globrow-startglobrow+1
    DO i = 1, M
      CALL RANDOM_NUMBER(randval); randval = rmin+(rmax-rmin)*randval
      lelement(1, globrowoff)%b(i,1) = randval
    END DO

    !-------------------------------------------
    !Save a copy of this original RHS
    !-------------------------------------------
    orig(globrowoff)%b = lelement(1,globrowoff)%b
  END DO

  !-------------------------------------------
  rhsdirtied = .TRUE.

END SUBROUTINE SetRandomRHS

!-------------------------------------------------------------------------------
!>
!! To be invoked, after solve, by user
!<
!-------------------------------------------------------------------------------
SUBROUTINE Finalize( do_mpifinalize )
  !-----------------------------------------------------------------------------
  ! Formal arguments
  !-----------------------------------------------------------------------------
  LOGICAL, INTENT(IN) :: do_mpifinalize !<Invoke MPI_Finalize?
  !-----------------------------------------------------------------------------
  ! Local variables
  !-----------------------------------------------------------------------------
  INTEGER :: level !<Loop variable
  INTEGER :: globrow !<Loop variable
  INTEGER :: mpierr !<MPI error code
  INTEGER :: tempinvcount, tempmulcount, tempsolcount !<Temp to avoid div by 0

  !-----------------------------------------------------------------------------
  IF(KPDBG) WRITE(OFU,*); CALL FL(OFU)
  IF(KPDBG) WRITE(OFU,*) '------ Finalizing start ------'; CALL FL(OFU)

  !-----------------------------------------------------------------------------
  CALL PLBFinalize()

  !-----------------------------------------------------------------------------
  ! Deallocate mats and hats, if any
  !-----------------------------------------------------------------------------
  IF ( ALLOCATED(lelement) ) THEN
    DO level = LBOUND(lelement,1), UBOUND(lelement,1)
      DO globrow = LBOUND(lelement,2), UBOUND(lelement,2)
        IF ( ALLOCATED(lelement(level,globrow)%L) ) THEN
          DEALLOCATE( lelement(level, globrow)%L )
        END IF
        IF ( ALLOCATED(lelement(level,globrow)%D) ) THEN
          DEALLOCATE( lelement(level, globrow)%D )
        END IF
        IF ( ALLOCATED(lelement(level,globrow)%U) ) THEN
          DEALLOCATE( lelement(level, globrow)%U )
        END IF
        IF ( ALLOCATED(lelement(level,globrow)%b) ) THEN
          DEALLOCATE( lelement(level, globrow)%b )
        END IF
        IF ( ALLOCATED(lelement(level,globrow)%pivot) ) THEN
          DEALLOCATE( lelement(level, globrow)%pivot )
        END IF
      END DO
    END DO
    DEALLOCATE(lelement)
  END IF
  IF ( ALLOCATED(orig) ) THEN
    DO globrow = LBOUND(orig,1), UBOUND(orig,1)
        IF ( ALLOCATED(orig(globrow)%L) ) THEN
          DEALLOCATE( orig(globrow)%L )
        END IF
        IF ( ALLOCATED(orig(globrow)%D) ) THEN
          DEALLOCATE( orig(globrow)%D )
        END IF
        IF ( ALLOCATED(orig(globrow)%U) ) THEN
          DEALLOCATE( orig(globrow)%U )
        END IF
        IF ( ALLOCATED(orig(globrow)%b) ) THEN
          DEALLOCATE( orig(globrow)%b )
        END IF
        IF ( ALLOCATED(orig(globrow)%pivot) ) THEN
          DEALLOCATE( orig(globrow)%pivot )
        END IF
    END DO
    DEALLOCATE(orig)
  END IF

  IF ( ALLOCATED(selement) ) THEN
    DO globrow = 1, N, 1
      IF ( ALLOCATED(selement(globrow)%x) ) THEN
        DEALLOCATE( selement(globrow)%x )
      END IF
    END DO
    DEALLOCATE(selement)
  END IF

  !-----------------------------------------------------------------------------
#if defined(MPI_OPT)
  IF ( usebarriers ) THEN
      IF(KPDBG) WRITE(OFU,*) 'Barrier in finalize'; CALL FL(OFU)
      CALL MPI_Barrier( SIESTA_COMM, mpierr )
      IF(KPDBG) WRITE(OFU,*) 'Done barrier in finalize'; CALL FL(OFU)
  END IF
  !-----------------------------------------------------------------------------
  IF ( do_mpifinalize ) THEN
      CALL MPI_Finalize( mpierr )
      use_mpiwtime = .FALSE.
  END IF
#endif

  !-----------------------------------------------------------------------------
  tempmulcount=totmatmulcount; IF ( tempmulcount .LE. 0 ) tempmulcount = 1
  tempinvcount=totinvcount; IF ( tempinvcount .LE. 0 ) tempinvcount = 1
  tempsolcount=totmatsolcount; IF ( tempsolcount .LE. 0 ) tempsolcount = 1
  IF(KPDBG) WRITE(OFU,*) 'N=', N, ' M=', M, ' P=', P, ' rank=', rank
  IF(KPDBG) WRITE(OFU,'(A,F6.1,A)') 'Memory        ', membytes/1e6, ' MB'
  IF(KPDBG) WRITE(OFU,'(A,F8.4,A)') 'Computation   ', tottime-totcommtime,' sec'
  IF(KPDBG) WRITE(OFU,'(A,F8.4,A)') 'Communication ', totcommtime,' sec'
  IF(KPDBG) WRITE(OFU,'(A,I5.1,A,F8.4,A,F8.4,A)') 'Matrix inv ',totinvcount, ' * ', &
                totinvtime/tempinvcount,' sec = ', totinvtime, ' sec'
  IF(KPDBG) WRITE(OFU,'(A,I5.1,A,F8.4,A,F8.4,A)') 'Matrix mul ',totmatmulcount, ' * ', &
                totmatmultime/tempmulcount, ' sec = ', totmatmultime, ' sec'
  IF(KPDBG) WRITE(OFU,'(A,I5.1,A,F8.4,A,F8.4,A)') 'Matrix sol ',totmatsolcount, ' * ', &
                totmatsoltime/tempsolcount, ' sec = ', totmatsoltime, ' sec'
  IF(KPDBG) WRITE(OFU,*) 'Finalized rank ', rank
  IF(KPDBG) WRITE(OFU,*) '------ Finalization end ------'; CALL FL(OFU)
END SUBROUTINE Finalize

!-------------------------------------------------------------------------------
!>
!! Is given integer even?
!<
!-------------------------------------------------------------------------------
LOGICAL FUNCTION ISEVEN( num )
  INTEGER, INTENT(IN) :: num !<A number
  ISEVEN = MOD(num,2) .EQ. 0
END FUNCTION ISEVEN

!-------------------------------------------------------------------------------
!>
!! Is given integer odd?
!<
!-------------------------------------------------------------------------------
LOGICAL FUNCTION ISODD( num )
  INTEGER, INTENT(IN) :: num !<A number
  ISODD = (.NOT. ISEVEN(num))
END FUNCTION ISODD

!-------------------------------------------------------------------------------
!>
!!Determines the global row number of a given local row number at a given level
!<
!-------------------------------------------------------------------------------
FUNCTION LR2GR( locrow, level ) result ( globrow )
  !-----------------------------------------------------------------------------
  ! Formal arguments
  !-----------------------------------------------------------------------------
  INTEGER, INTENT(IN) :: locrow !<local row number of global row, at given level
  INTEGER, INTENT(IN) :: level !<Level at which locrow's position is given
  INTEGER :: globrow !<Returned: Row number in original input (level 1)

  !-----------------------------------------------------------------------------
  ! Local variables
  !-----------------------------------------------------------------------------
  INTEGER :: i !<Loop variable

  !-----------------------------------------------------------------------------
  ! Invert the rule r=(r+1)/2 backward for all levels
  !-----------------------------------------------------------------------------
  globrow = locrow
  LevelLoop: DO i = level-1, 1, -1
    globrow = 2*globrow - 1
  END DO LevelLoop

  !-----------------------------------------------------------------------------
  ! DEBUGGING: Just double-check
  !-----------------------------------------------------------------------------
  IF ( ((2.0 ** (level-1)) * (locrow-1) + 1) .NE. globrow ) THEN
    CALL ASSERT(.FALSE.,'blocktri #46') !Consistency check failed
  END IF

  RETURN
END FUNCTION LR2GR

!-------------------------------------------------------------------------------
!>
!!Determines the local row number of "globrow" global row number when globrow
!!participates at a given level; returns zero if this globrow does not operate
!!at the given level
!<
!-------------------------------------------------------------------------------
FUNCTION GR2LR( globrow, level ) result ( locrow )
  !-----------------------------------------------------------------------------
  ! Formal arguments
  !-----------------------------------------------------------------------------
  INTEGER, INTENT(IN) :: globrow !<Row number in original input (level 1)
  INTEGER, INTENT(IN) :: level !<Level at which globrow's position is needed
  INTEGER :: locrow !<Returned: local row number of global row, at given level

  !-----------------------------------------------------------------------------
  ! Local variables
  !-----------------------------------------------------------------------------
  INTEGER :: i !<Loop variable

  !-----------------------------------------------------------------------------
  ! Apply the rule r=(r+1)/2 for all odd-numbered rows
  ! Any time r becomes even, give up
  !-----------------------------------------------------------------------------
  locrow = globrow
  LevelLoop: DO i = 1, level-1, 1
    IF ( ISEVEN(locrow) ) THEN !Even-numbered rows don't go any further in level
      locrow = 0
      EXIT LevelLoop !Stop looping
    END IF
    locrow = (locrow+1) / 2
  END DO LevelLoop

  RETURN
END FUNCTION GR2LR

!-------------------------------------------------------------------------------
!>
!!Determines the rank of the task holding the given global row (at level 1)
!<
!-------------------------------------------------------------------------------
FUNCTION GR2Rank( globrow ) result ( outrank )
  !-----------------------------------------------------------------------------
  ! Formal arguments
  !-----------------------------------------------------------------------------
  INTEGER, INTENT(IN) :: globrow !<Row number in original input (level 1)
  INTEGER :: outrank !<Returned: Rank holding the given global row

  !-----------------------------------------------------------------------------
  ! Local variables
  !-----------------------------------------------------------------------------
  INTEGER :: m !Average integral number of rows per processor
  INTEGER :: spill !Spill over beyond prefect loading of rows to processors

  IF ( (globrow .LT. 1) .OR. (globrow .GT. N) ) THEN
    outrank = -1
  ELSE
    m = N / P
    spill = MOD( N, P )
    IF ( globrow .LE. (m+1)*spill ) THEN
      outrank = (globrow-1)/(m+1)
    ELSE
      outrank = (globrow-1 - (m+1)*spill)/m + spill
    END IF
  END IF

  RETURN
END FUNCTION GR2Rank

!-------------------------------------------------------------------------------
!>
!!Determines the rank of the task holding the given local row (at a given level)
!<
!-------------------------------------------------------------------------------
FUNCTION LR2Rank( locrow, level ) result ( outrank )
  !-----------------------------------------------------------------------------
  ! Formal arguments
  !-----------------------------------------------------------------------------
  INTEGER, INTENT(IN) :: locrow !<Row number at a given level
  INTEGER, INTENT(IN) :: level !<Local row number is at this level
  INTEGER :: outrank !<Returned: Rank holding the given global row

  !-----------------------------------------------------------------------------
  ! Local variables
  !-----------------------------------------------------------------------------
  INTEGER :: globrow !Global row number of given locrow

  globrow = LR2GR( locrow, level )
  outrank = GR2Rank( globrow )

  RETURN
END FUNCTION LR2Rank

!-------------------------------------------------------------------------------
!>
!!Compute the matrix multiplications in the forward solve at some/any level
!!The row should be odd at the level,
!!The locrow should obey startlocrow <= locrow <= endlocrow
!<
!-------------------------------------------------------------------------------
#if defined(MPI_OPT)
SUBROUTINE ComputeForwardOddRowHats(locrow, level, startlocrow, endlocrow,bonly)
  !-----------------------------------------------------------------------------
  ! Formal arguments
  !-----------------------------------------------------------------------------
  INTEGER, INTENT(IN) :: locrow !<Row number at a given level
  INTEGER, INTENT(IN) :: level !<Local row number is at this level
  INTEGER, INTENT(IN) :: startlocrow !<Top-most local row number at given level
  INTEGER, INTENT(IN) :: endlocrow !<Bot-most local row number at given level
  LOGICAL, INTENT(IN) :: bonly !<Show compute/update only the b, and not L,D,U?
  !-----------------------------------------------------------------------------
  ! Local variables
  !-----------------------------------------------------------------------------
  INTEGER :: globrow !<Temporary variable
  INTEGER :: globrowoff !<Temporary variable
  INTEGER :: prevglobrowoff !<Temporary variable
  INTEGER :: nextglobrowoff !<Temporary variable

  !-------------------------------------------------------------
  ! Nothing to do if last level
  !-------------------------------------------------------------
  IF ( level .GE. nlevels ) THEN
    IF(KPDBG) WRITE(OFU,*) 'mat-mul last lvl: no op '; CALL FL(OFU)
    IF ( rank .NE. 0 ) THEN
      IF(KPDBG) WRITE(OFU,*) 'mat-mul last lvl: impossible ',rank; CALL FL(OFU)
      CALL ASSERT(.FALSE.,'blocktri #47')
    END IF
    RETURN
  END IF

  !-------------------------------------------------------------
  ! Sanity check
  !-------------------------------------------------------------
  IF ( ISEVEN( locrow ) ) THEN
    IF(KPDBG) WRITE(OFU,*) 'mat-mul: odd only ',globrow,' ',locrow; CALL FL(OFU); CALL ASSERT(.FALSE.,'blocktri #48')
  END IF
  IF ( .NOT. ( (startlocrow .LE. locrow) .AND. (locrow .LE. endlocrow) ) ) THEN
    IF(KPDBG) WRITE(OFU,*) 'mat-mul: locrow prob ',globrow,' ',locrow;CALL FL(OFU);CALL ASSERT(.FALSE.,'blocktri #49')
  END IF

  !-------------------------------------------------------------
  globrow = LR2GR( locrow, level )
  globrowoff = globrow - startglobrow + 1
  IF(KPDBG) WRITE(OFU,*) '  Compute mat-mul ',globrow,' ',locrow; CALL FL(OFU)

  !-------------------------------------------------------------
  ! Ensure there are mats at (thislevel, thisrow)
  !-------------------------------------------------------------
  IF ( .NOT. ALLOCATED(lelement(level,globrowoff)%D ) ) THEN
    IF(KPDBG) WRITE(OFU,*) '  Copying bad D'; CALL FL(OFU); CALL ASSERT(.FALSE.,'blocktri #50')
  END IF
  IF ( .NOT. ALLOCATED(lelement(level,globrowoff)%L ) ) THEN
    IF(KPDBG) WRITE(OFU,*) '  Copying bad L'; CALL FL(OFU); CALL ASSERT(.FALSE.,'blocktri #51')
  END IF
  IF ( .NOT. ALLOCATED(lelement(level,globrowoff)%U ) ) THEN
    IF(KPDBG) WRITE(OFU,*) '  Copying bad U'; CALL FL(OFU); CALL ASSERT(.FALSE.,'blocktri #52')
  END IF
  IF ( .NOT. ALLOCATED(lelement(level,globrowoff)%b ) ) THEN
    IF(KPDBG) WRITE(OFU,*) '  Copying bad b'; CALL FL(OFU); CALL ASSERT(.FALSE.,'blocktri #53')
  END IF

  !-------------------------------------------------------------
  ! Ensure there is space at (nextlevel, thisrow)
  !-------------------------------------------------------------
  IF ( .NOT. ALLOCATED(lelement(level+1,globrowoff)%D ) ) THEN
    IF(KPDBG) WRITE(OFU,*) '  Copying bad D +1'; CALL FL(OFU); CALL ASSERT(.FALSE.,'blocktri #54')
  END IF
  IF ( .NOT. ALLOCATED(lelement(level+1,globrowoff)%L ) ) THEN
    IF(KPDBG) WRITE(OFU,*) '  Copying bad L +1'; CALL FL(OFU); CALL ASSERT(.FALSE.,'blocktri #55')
  END IF
  IF ( .NOT. ALLOCATED(lelement(level+1,globrowoff)%U ) ) THEN
    IF(KPDBG) WRITE(OFU,*) '  Copying bad U +1'; CALL FL(OFU); CALL ASSERT(.FALSE.,'blocktri #56')
  END IF
  IF ( .NOT. ALLOCATED(lelement(level+1,globrowoff)%b ) ) THEN
    IF(KPDBG) WRITE(OFU,*) '  Copying bad b +1'; CALL FL(OFU); CALL ASSERT(.FALSE.,'blocktri #57')
  END IF

  !-------------------------------------------------------------
  ! Copy over the D and b as is first
  !-------------------------------------------------------------
  IF ( .NOT. bonly ) THEN
    lelement(level+1,globrowoff)%D = lelement(level,globrowoff)%D
  END IF
  lelement(level+1,globrowoff)%b = lelement(level,globrowoff)%b

  !-------------------------------------------------------------
  ! Check if this row has a top neighbor at this level
  !-------------------------------------------------------------
  IF ( LR2Rank(locrow-1,level) .LT. 0 ) THEN
    ! No rank before ours
    IF ( .NOT. bonly ) THEN
      lelement(level+1,globrowoff)%U = 0
      IF(KPDBG) WRITE(OFU,*) '  ZERO L'; CALL FL(OFU)
    END IF
  ELSE
    IF ( locrow .EQ. startlocrow ) THEN
      prevglobrowoff = 0
    ELSE IF ( locrow .GT. startlocrow ) THEN
      prevglobrowoff = LR2GR( locrow-1, level ) - startglobrow + 1
    ELSE
      IF(KPDBG) WRITE(OFU,*) 'mat-mul: impossible prev'; CALL FL(OFU); CALL ASSERT(.FALSE.,'blocktri #58')
    END IF

    IF ( .NOT. bonly ) THEN
      !------------------------------------------------
      ! D = D - L*Uhat(prev)
      !------------------------------------------------
      CALL BSystemClock(mattimer1)
      CALL PLBDGEMM( -ONE, lelement(level,globrowoff)%L, &
                        lelement(level,prevglobrowoff)%U, &
                        ONE, lelement(level+1,globrowoff)%D )
      CALL BSystemClock(mattimer2)
      CALL ChargeTime( totmatmultime, mattimer2, mattimer1, totmatmulcount )
      IF(KPDBG) WRITE(OFU,*) '  Dtop'; CALL FL(OFU)

      !------------------------------------------------
      ! L = -L*Lhat(prev)
      !------------------------------------------------
      CALL BSystemClock(mattimer1)
      CALL PLBDGEMM( -ONE, lelement(level,globrowoff)%L, &
                        lelement(level,prevglobrowoff)%L, &
                        ZERO, lelement(level+1,globrowoff)%L )
      CALL BSystemClock(mattimer2)
      CALL ChargeTime( totmatmultime, mattimer2, mattimer1, totmatmulcount )
      IF(KPDBG) WRITE(OFU,*) '  L'; CALL FL(OFU)
    END IF

    !------------------------------------------------
    ! b = b - L*bhat(prev)
    !------------------------------------------------
    CALL PLBDGEMV( -ONE, lelement(level,globrowoff)%L, &
                      lelement(level,prevglobrowoff)%b(:,1), &
                      ONE, lelement(level+1,globrowoff)%b(:,1) )
    IF(KPDBG) WRITE(OFU,*) '  btop'; CALL FL(OFU)
  END IF

  !-------------------------------------------------------------
  ! Check if this row has a bottom neighbor at this level
  !-------------------------------------------------------------
  IF ( LR2Rank(locrow+1,level) .LT. 0 ) THEN
    ! No rank after ours
    IF ( .NOT. bonly ) THEN
      lelement(level+1,globrowoff)%U = 0
      IF(KPDBG) WRITE(OFU,*) '  ZERO U'; CALL FL(OFU)
    END IF
  ELSE
    IF ( locrow .EQ. endlocrow ) THEN
      nextglobrowoff = endglobrow - startglobrow + 2
    ELSE IF ( locrow .LT. endlocrow ) THEN
      nextglobrowoff = LR2GR( locrow+1, level ) - startglobrow + 1
    ELSE
      IF(KPDBG) WRITE(OFU,*) 'mat-mul: impossible next'; CALL FL(OFU); CALL ASSERT(.FALSE.,'blocktri #59')
    END IF

    IF ( .NOT. bonly ) THEN
      !------------------------------------------------
      ! D = D - U*Lhat(next)
      !------------------------------------------------
      CALL BSystemClock(mattimer1)
      CALL PLBDGEMM( -ONE, lelement(level,globrowoff)%U, &
                        lelement(level,nextglobrowoff)%L, &
                        ONE, lelement(level+1,globrowoff)%D )
      CALL BSystemClock(mattimer2)
      CALL ChargeTime( totmatmultime, mattimer2, mattimer1, totmatmulcount )
      IF(KPDBG) WRITE(OFU,*) '  Dbot'; CALL FL(OFU)

      !------------------------------------------------
      ! U = -U*Uhat(next)
      !------------------------------------------------
      CALL BSystemClock(mattimer1)
      CALL PLBDGEMM( -ONE, lelement(level,globrowoff)%U, &
                        lelement(level,nextglobrowoff)%U, &
                        ZERO, lelement(level+1,globrowoff)%U )
      CALL BSystemClock(mattimer2)
      CALL ChargeTime( totmatmultime, mattimer2, mattimer1, totmatmulcount )
      IF(KPDBG) WRITE(OFU,*) '  U'; CALL FL(OFU)
    END IF

    !------------------------------------------------
    ! b = b - U*bhat(next)
    !------------------------------------------------
    CALL PLBDGEMV( -ONE, lelement(level,globrowoff)%U, &
                      lelement(level,nextglobrowoff)%b(:,1), &
                      ONE, lelement(level+1,globrowoff)%b(:,1) )
    IF(KPDBG) WRITE(OFU,*) '  bbot'; CALL FL(OFU)
  END IF
END SUBROUTINE ComputeForwardOddRowHats

!-------------------------------------------------------------------------------
!>
!!BCYCLIC forward phase; to be called after Initialize
!<
!-------------------------------------------------------------------------------
SUBROUTINE ForwardSolve
  !-----------------------------------------------------------------------------
  ! Local variables
  !-----------------------------------------------------------------------------
  INTEGER :: level !<Loop variable
  INTEGER :: locrow !<Loop variable for "local" row number
  INTEGER :: globrow !<Loop variable for "global" row number
  INTEGER :: globrowoff !<Temporary variable
  INTEGER :: startlocrow !<Starting row number of this processor at a level
  INTEGER :: endlocrow !<Ending row number (inclusive) at a given level
  INTEGER :: rowoffset !<Loop variable to generate row+/-1
  INTEGER :: nbrrank !<Temporary variable
  INTEGER :: nbrmpireqindx !<Index of next MPI_Irecv request
  INTEGER :: nbrmpireqcnt !<Number of MPI_Irecv/Isend requests initiated
  INTEGER :: nbrmpinirecvs !<Number of MPI_Irecv requests initiated
  INTEGER :: nbrmpinisends !<Number of MPI_Isend requests initiated
  INTEGER :: nbrmpireq(6) !<MPI_Irecv/Isend request handle for top or bot nbr
  INTEGER :: nbrmpierr(6) !<MPI_Irecv/Isend error for top or bot nbr
  INTEGER :: mpiwaiterr !<MPI_Waitall/any error
  INTEGER :: mpierr !<Generic use for MPI
  INTEGER :: nbrmpistatuses(MPI_STATUS_SIZE,6) !<MPI_Waitall statuses
  INTEGER :: nbrmpistatus(MPI_STATUS_SIZE) !<MPI_Waitany status
  INTEGER :: waitmatchindex !<MPI_Waitany match index
  INTEGER :: stag !<Tag number for Isends
  INTEGER :: reqi !<Looping variable for requests
  INTEGER :: blaserr !<Generic use for BLAS

  CALL BSystemClock(globtimer1)

  !-----------------------------------------------------------------------------
  IF(KPDBG) WRITE(OFU,*); CALL FL(OFU)
  IF(KPDBG) WRITE(OFU,*) '------ Forward solve start ------'; CALL FL(OFU)
  IF(KPDBG) WRITE(OFU,*); CALL FL(OFU)

  !-----------------------------------------------------------------------------
  CALL PLBForwardInitialize()

  !-----------------------------------------------------------------------------
  ! Forward solving loop in which even rows compute inverses and odd rows
  ! compute matrix-matrix products for new terms
  !-----------------------------------------------------------------------------
  ForwardSolveLoop: DO level = 1, nlevels, 1
    !---------------------------------------------------------------------------
    IF ( usebarriers ) THEN
      IF(KPDBG) WRITE(OFU,*) 'Barrier at forward level ', level; CALL FL(OFU)
      CALL MPI_Barrier( SIESTA_COMM, mpierr )
      IF(KPDBG) WRITE(OFU,*) 'Done barrier at forward level ', level; CALL FL(OFU)
    END IF

    !---------------------------------------------------------------------------
    !Determine start and end local rows at current level
    !---------------------------------------------------------------------------
    startlocrow = N+1 !Start with an invalid value (greater than endlocrow)
    endlocrow = -1 !Start with an invalid value (less than startlocrow)
    DO globrow = startglobrow, endglobrow, 1
      locrow = GR2LR( globrow, level ) !globrow's place at current level
      IF ( locrow .GT. 0 ) THEN !This global row participates at current level
        IF ( locrow .LT. startlocrow ) THEN !A valid, earlier local row
          startlocrow = locrow
        END IF
        IF ( locrow .GT. endlocrow ) THEN !A valid, later local row
          endlocrow = locrow
        END IF
      END IF
    END DO

    !---------------------------------------------------------------------------
    !We may have run out of work; see if we have a valid range of rows left
    !---------------------------------------------------------------------------
    IF ( startlocrow .GT. endlocrow ) THEN !No rows left at/above this level
      IF(KPDBG) WRITE(OFU,*) 'Nothing to do at forward level ', level; CALL FL(OFU)
      CALL PLBForwardInitializeLevel( level, .FALSE. )
      CALL PLBForwardFinalizeLevel( level, .FALSE. )
      CYCLE ForwardSolveLoop !Move on to next level; don't EXIT (for barriers)
    ELSE
      CALL PLBForwardInitializeLevel( level, .TRUE. )
    END IF

    IF(KPDBG) WRITE(OFU,*) '**** Forward level ', level, ' ****'; CALL FL(OFU)

    !---------------------------------------------------------------------------
    !Allocate memory for incoming nbr values, if any, at this level
    !---------------------------------------------------------------------------
    IF ( ISODD(startlocrow) .AND. &
     (LR2Rank(startlocrow-1,level) .GE. 0) ) THEN
      globrowoff = 0
      IF ( .NOT. ALLOCATED( lelement(level, globrowoff)%L ) ) THEN
        ALLOCATE( lelement(level, globrowoff)%L(M,M) )
        ALLOCATE( lelement(level, globrowoff)%U(M,M) )
        ALLOCATE( lelement(level, globrowoff)%b(M,1) )
        CALL ChargeMemory( (2*M*M+M)*dpsz )
      END IF
      lelement(level, globrowoff)%L = 0
      lelement(level, globrowoff)%U = 0
      lelement(level, globrowoff)%b = 0
    END IF
    IF ( ISODD(endlocrow) .AND. &
     (LR2Rank(endlocrow+1,level) .GE. 0) ) THEN
      globrowoff = endglobrow-startglobrow+2
      IF ( .NOT. ALLOCATED( lelement(level, globrowoff)%L ) ) THEN
        ALLOCATE( lelement(level, globrowoff)%L(M,M) )
        ALLOCATE( lelement(level, globrowoff)%U(M,M) )
        ALLOCATE( lelement(level, globrowoff)%b(M,1) )
        CALL ChargeMemory( (2*M*M+M)*dpsz )
      END IF
      lelement(level, globrowoff)%L = 0
      lelement(level, globrowoff)%U = 0
      lelement(level, globrowoff)%b = 0
    END IF

    !---------------------------------------------------------------------------
    !Allocate memory at next level for those rows that are odd at current level
    !---------------------------------------------------------------------------
    IF ( level .LT. nlevels ) THEN
      DO locrow = startlocrow, endlocrow, 1
        IF ( ISODD(locrow) ) THEN
          globrow = LR2GR( locrow, level )
          globrowoff = globrow-startglobrow+1
          IF ( .NOT. ALLOCATED( lelement(level+1, globrowoff)%L ) ) THEN
            ALLOCATE( lelement(level+1, globrowoff)%L(M,M) )
            ALLOCATE( lelement(level+1, globrowoff)%D(M,M) )
            ALLOCATE( lelement(level+1, globrowoff)%U(M,M) )
            ALLOCATE( lelement(level+1, globrowoff)%b(M,1) )
            CALL ChargeMemory( (2*M*M+M)*dpsz )
          END IF
          lelement(level+1, globrowoff)%L = 0
          lelement(level+1, globrowoff)%D = 0
          lelement(level+1, globrowoff)%U = 0
          lelement(level+1, globrowoff)%b = 0
        END IF
      END DO
    END IF

    !---------------------------------------------------------------------------
    !Reset requests
    !---------------------------------------------------------------------------
    DO nbrmpireqindx = 1, 6
      nbrmpireq(nbrmpireqindx) = MPI_REQUEST_NULL
    END DO

    !---------------------------------------------------------------------------
    !Pre-post the expected MPI receives via a non-blocking recv primitive
    !Pre-posting can help overlap communication with computing inverses
    !---------------------------------------------------------------------------
    nbrmpireqindx = 1
    nbrmpinirecvs = 0
    nbrrank = LR2Rank(startlocrow-1,level)
    IF ( ISODD(startlocrow) .AND. (nbrrank .GE. 0) ) THEN
      !-------------------------------------------------------------------------
      !Our top row at this level is odd and top row's previous is valid
      !We will get the previous (even-numbered) row's L-hat, U-hat, and b-hat
      !-------------------------------------------------------------------------
      globrowoff = 0
      IF(KPDBG) WRITE(OFU,*) '  Irecv ',startlocrow-1,' ',nbrrank,' '; CALL FL(OFU)

      IF ( .NOT. ALLOCATED(lelement(level,globrowoff)%L) ) THEN
        IF(KPDBG) WRITE(OFU,*)'  Irecv bad L'; CALL FL(OFU); CALL ASSERT(.FALSE.,'blocktri #60')
      END IF
      IF ( .NOT. ALLOCATED(lelement(level,globrowoff)%U) ) THEN
        IF(KPDBG) WRITE(OFU,*)'  Irecv bad U'; CALL FL(OFU); CALL ASSERT(.FALSE.,'blocktri #61')
      END IF
      IF ( .NOT. ALLOCATED(lelement(level,globrowoff)%b) ) THEN
        IF(KPDBG) WRITE(OFU,*)'  Irecv bad b'; CALL FL(OFU); CALL ASSERT(.FALSE.,'blocktri #62')
      END IF

      CALL BSystemClock(loctimer1)
      CALL MPI_Irecv( lelement(level, globrowoff)%L, M*M, MPI_REAL8, &
       nbrrank, 1, SIESTA_COMM, nbrmpireq(nbrmpireqindx), &
       nbrmpierr(nbrmpireqindx) )
      nbrmpireqindx = nbrmpireqindx + 1
      nbrmpinirecvs = nbrmpinirecvs + 1
      IF(KPDBG) WRITE(OFU,*) '  Irecv 1 top ',nbrrank,' ',nbrmpireqindx-1; CALL FL(OFU)

      CALL MPI_Irecv( lelement(level, globrowoff)%U, M*M, MPI_REAL8, &
       nbrrank, 2, SIESTA_COMM, nbrmpireq(nbrmpireqindx), &
       nbrmpierr(nbrmpireqindx))
      nbrmpireqindx = nbrmpireqindx + 1
      nbrmpinirecvs = nbrmpinirecvs + 1
      IF(KPDBG) WRITE(OFU,*) '  Irecv 2 top ',nbrrank,' ',nbrmpireqindx-1; CALL FL(OFU)

      CALL MPI_Irecv( lelement(level, globrowoff)%b, M, MPI_REAL8, &
       nbrrank, 3, SIESTA_COMM, nbrmpireq(nbrmpireqindx), &
       nbrmpierr(nbrmpireqindx))
      nbrmpireqindx = nbrmpireqindx + 1
      nbrmpinirecvs = nbrmpinirecvs + 1
      IF(KPDBG) WRITE(OFU,*) '  Irecv 3 top ',nbrrank,' ',nbrmpireqindx-1; CALL FL(OFU)
      CALL BSystemClock(loctimer2)
      CALL ChargeTime( totcommtime, loctimer2, loctimer1, totcommcount )
    END IF

    nbrrank = LR2Rank(endlocrow+1,level)
    IF ( ISODD(endlocrow) .AND. (nbrrank .GE. 0) ) THEN
      !-------------------------------------------------------------------------
      !Our bottom row at this level is odd and bottom row's next is valid
      !We will get the next (even-numbered) row's L-hat, U-hat, and b-hat
      !-------------------------------------------------------------------------
      globrowoff = endglobrow-startglobrow+2
      IF(KPDBG) WRITE(OFU,*) ' Irecv ',endlocrow+1,' ', nbrrank; CALL FL(OFU)

      IF ( .NOT. ALLOCATED(lelement(level,globrowoff)%L) ) THEN
        IF(KPDBG) WRITE(OFU,*)'Irecv bad L'; CALL FL(OFU); CALL ASSERT(.FALSE.,'blocktri #63')
      END IF
      IF ( .NOT. ALLOCATED(lelement(level,globrowoff)%U) ) THEN
        IF(KPDBG) WRITE(OFU,*)'Irecv bad U'; CALL FL(OFU); CALL ASSERT(.FALSE.,'blocktri #64')
      END IF
      IF ( .NOT. ALLOCATED(lelement(level,globrowoff)%b) ) THEN
        IF(KPDBG) WRITE(OFU,*)'Irecv bad b'; CALL FL(OFU); CALL ASSERT(.FALSE.,'blocktri #65')
      END IF

      CALL BSystemClock(loctimer1)
      CALL MPI_Irecv( lelement(level, globrowoff)%L, M*M, MPI_REAL8, &
       nbrrank, 4, SIESTA_COMM, nbrmpireq(nbrmpireqindx), &
       nbrmpierr(nbrmpireqindx))
      nbrmpireqindx = nbrmpireqindx + 1
      nbrmpinirecvs = nbrmpinirecvs + 1
      IF(KPDBG) WRITE(OFU,*) '  Irecv 4 bot ',nbrrank,' ',nbrmpireqindx-1; CALL FL(OFU)

      CALL MPI_Irecv( lelement(level, globrowoff)%U, M*M, MPI_REAL8, &
       nbrrank, 5, SIESTA_COMM, nbrmpireq(nbrmpireqindx), &
       nbrmpierr(nbrmpireqindx))
      nbrmpireqindx = nbrmpireqindx + 1
      nbrmpinirecvs = nbrmpinirecvs + 1
      IF(KPDBG) WRITE(OFU,*) '  Irecv 5 bot ',nbrrank,' ',nbrmpireqindx-1; CALL FL(OFU)

      CALL MPI_Irecv( lelement(level, globrowoff)%b, M, MPI_REAL8, &
       nbrrank, 6, SIESTA_COMM, nbrmpireq(nbrmpireqindx), &
       nbrmpierr(nbrmpireqindx))
      nbrmpireqindx = nbrmpireqindx + 1
      nbrmpinirecvs = nbrmpinirecvs + 1
      IF(KPDBG) WRITE(OFU,*) '  Irecv 6 bot ',nbrrank,' ',nbrmpireqindx-1; CALL FL(OFU)
      CALL BSystemClock(loctimer2)
      CALL ChargeTime( totcommtime, loctimer2, loctimer1, totcommcount )
    END IF

    !---------------------------------------------------------------------------
    !Compute inverses for even-numbered rows at this level
    !---------------------------------------------------------------------------
    nbrmpinisends = 0
    DO locrow = startlocrow, endlocrow, 1
      IF ( ISEVEN( locrow ) .OR. (level .EQ. nlevels) ) THEN
        globrow = LR2GR( locrow, level )
        globrowoff = globrow-startglobrow+1

        IF(KPDBG) WRITE(OFU,*) ' Compute even hats ',globrow,' ', locrow; CALL FL(OFU)

        !-----------------------------------------------------------------------
        !Sanity checks
        !-----------------------------------------------------------------------
        IF ( level .EQ. nlevels ) THEN
          IF ( (rank .NE. 0) .OR. &
               (startlocrow .NE. endlocrow) .OR. &
               (locrow .NE. 1) .OR. (globrow .NE. 1) ) THEN
            IF(KPDBG) WRITE(OFU,*) ' EVEN ERROR ',globrow,' ', locrow; CALL FL(OFU)
            CALL ASSERT(.FALSE.,'blocktri #66')
          END IF
        END IF
        IF ( .NOT. ALLOCATED(lelement(level,globrowoff)%D) ) THEN
          IF(KPDBG) WRITE(OFU,*)'Chats bad D'; CALL FL(OFU); CALL ASSERT(.FALSE.,'blocktri #67')
        END IF
        IF ( .NOT. ALLOCATED(lelement(1,globrowoff)%pivot) ) THEN
          IF(KPDBG) WRITE(OFU,*)'Chats bad pivot'; CALL FL(OFU); CALL ASSERT(.FALSE.,'blocktri #68')
        END IF
        IF ( .NOT. ALLOCATED(lelement(level,globrowoff)%L) ) THEN
          IF(KPDBG) WRITE(OFU,*)'Chats bad L'; CALL FL(OFU); CALL ASSERT(.FALSE.,'blocktri #69')
        END IF
        IF ( .NOT. ALLOCATED(lelement(level,globrowoff)%U) ) THEN
          IF(KPDBG) WRITE(OFU,*)'Chats bad U'; CALL FL(OFU); CALL ASSERT(.FALSE.,'blocktri #70')
        END IF
        IF ( .NOT. ALLOCATED(lelement(level,globrowoff)%b) ) THEN
          IF(KPDBG) WRITE(OFU,*)'Chats bad b'; CALL FL(OFU); CALL ASSERT(.FALSE.,'blocktri #71')
        END IF

        !-----------------------------------------------------------------------
        !Do LU factorization of D (in place into D itself)
        !-----------------------------------------------------------------------
        CALL BSystemClock(mattimer1)
        CALL PLBDGETRF( lelement(level,globrowoff)%D, &
                           lelement(1,globrowoff)%pivot, blaserr )
        CALL BSystemClock(mattimer2)
        CALL ChargeTime( totinvtime, mattimer2, mattimer1, totinvcount )
        IF (blaserr .NE. 0) THEN
          IF(KPDBG) WRITE(OFU,*) '  PLBDGETRF failed ', blaserr; CALL FL(OFU); CALL ASSERT(.FALSE.,'blocktri #72')
        END IF

        !-----------------------------------------------------------------------
        !Compute hats D^(-1)L and D^(-1)U if not last level
        !-----------------------------------------------------------------------
        IF ( level .LT. nlevels ) THEN
          CALL BSystemClock(mattimer1)
          CALL PLBDGETRS( M, lelement(level,globrowoff)%D, &
                             lelement(1,globrowoff)%pivot, &
                             lelement(level,globrowoff)%L, blaserr )
          CALL BSystemClock(mattimer2)
          CALL ChargeTime( totmatsoltime, mattimer2, mattimer1, totmatsolcount )
          IF (blaserr .NE. 0) THEN
            IF(KPDBG) WRITE(OFU,*) '   PLBDGETRS L failed'; CALL FL(OFU); CALL ASSERT(.FALSE.,'blocktri #73')
          END IF
          CALL BSystemClock(mattimer1)
          CALL PLBDGETRS( M, lelement(level,globrowoff)%D, &
                             lelement(1,globrowoff)%pivot, &
                             lelement(level,globrowoff)%U, blaserr )
          CALL BSystemClock(mattimer2)
          CALL ChargeTime( totmatsoltime, mattimer2, mattimer1, totmatsolcount )
          IF (blaserr .NE. 0) THEN
            IF(KPDBG) WRITE(OFU,*) '   PLBDGETRS U failed'; CALL FL(OFU); CALL ASSERT(.FALSE.,'blocktri #74')
          END IF
        END IF

        !-----------------------------------------------------------------------
        !Compute b-hats D^(-1)b
        !-----------------------------------------------------------------------
        CALL PLBDGETRS( 1, lelement(level,globrowoff)%D, &
                           lelement(1,globrowoff)%pivot, &
                           lelement(level,globrowoff)%b, blaserr )
        IF (blaserr .NE. 0) THEN
          IF(KPDBG) WRITE(OFU,*) 'PLBDGETRS b failed'; CALL FL(OFU); CALL ASSERT(.FALSE.,'blocktri #75')
        END IF

        !-----------------------------------------------------------------------
        !Send my Lhats, Uhats and bhats to neighbor row if that neighbor exists
        !and that neighbor row is hosted outside this rank
        !-----------------------------------------------------------------------
        DO rowoffset = -1, 1, 2
          nbrrank = LR2Rank(locrow+rowoffset, level)
          IF ( (nbrrank .GE. 0) .AND. (nbrrank .NE. rank) ) THEN
            CALL BSystemClock(loctimer1)
            IF ( .NOT. ALLOCATED(lelement(level,globrowoff)%L) ) THEN
              IF(KPDBG) WRITE(OFU,*)'ISend bad L'; CALL FL(OFU); CALL ASSERT(.FALSE.,'blocktri #76')
            END IF
            IF ( .NOT. ALLOCATED(lelement(level,globrowoff)%U) ) THEN
              IF(KPDBG) WRITE(OFU,*)'ISend bad U'; CALL FL(OFU); CALL ASSERT(.FALSE.,'blocktri #77')
            END IF
            IF ( .NOT. ALLOCATED(lelement(level,globrowoff)%b) ) THEN
              IF(KPDBG) WRITE(OFU,*)'ISend bad b'; CALL FL(OFU); CALL ASSERT(.FALSE.,'blocktri #78')
            END IF

            IF(KPDBG) WRITE(OFU,*)' ISend ',globrow,' ',locrow,' ',nbrrank;CALL FL(OFU)
            globrowoff = globrow-startglobrow+1
            stag = (((1-rowoffset))/2)*3
            CALL MPI_Isend( lelement(level,globrowoff)%L, M*M, MPI_REAL8, &
             nbrrank, stag+1, SIESTA_COMM, &
             nbrmpireq(nbrmpireqindx), nbrmpierr(nbrmpireqindx))
            nbrmpireqindx = nbrmpireqindx + 1
            nbrmpinisends = nbrmpinisends + 1
            IF(KPDBG) WRITE(OFU,*) '  Isend ',stag+1,' ',globrowoff; CALL FL(OFU)

            CALL MPI_Isend( lelement(level,globrowoff)%U, M*M, MPI_REAL8, &
             nbrrank, stag+2, SIESTA_COMM, &
             nbrmpireq(nbrmpireqindx), nbrmpierr(nbrmpireqindx))
            nbrmpireqindx = nbrmpireqindx + 1
            nbrmpinisends = nbrmpinisends + 1
            IF(KPDBG) WRITE(OFU,*) '  Isend ',stag+2,' ',globrowoff; CALL FL(OFU)

            CALL MPI_Isend( lelement(level,globrowoff)%b, M, MPI_REAL8, &
             nbrrank, stag+3, SIESTA_COMM, &
             nbrmpireq(nbrmpireqindx), nbrmpierr(nbrmpireqindx))
            nbrmpireqindx = nbrmpireqindx + 1
            nbrmpinisends = nbrmpinisends + 1
            IF(KPDBG) WRITE(OFU,*) '  Isend ',stag+3,' ',globrowoff; CALL FL(OFU)
            CALL BSystemClock(loctimer2)
            CALL ChargeTime( totcommtime, loctimer2, loctimer1, totcommcount )
          END IF
        END DO
      END IF
    END DO

    !---------------------------------------------------------------------------
    !Compute the matrix-matrix multiplications for non-boundary odd rows
    !---------------------------------------------------------------------------
    DO locrow = startlocrow, endlocrow, 1
      globrow = LR2GR( locrow, level )
      IF ( ISODD( locrow ) .AND. (level .NE. nlevels) ) THEN
        IF ( (locrow .NE. startlocrow) .AND. (locrow .NE. endlocrow) ) THEN
          IF(KPDBG) WRITE(OFU,*) ' Precomp mat-mul ',globrow,' ',locrow; CALL FL(OFU)
          CALL ComputeForwardOddRowHats( locrow, level, &
                                         startlocrow, endlocrow, .FALSE. )
        END IF
      END IF
    END DO

    !---------------------------------------------------------------------------
    !We have to wait for incoming even-numbered neighboring rows, before
    !computing inverses for the odd boundaries, if any.
    !Complete the send-receive of inverses, by doing a "wait all"
    !---------------------------------------------------------------------------
    nbrmpireqcnt = nbrmpireqindx - 1
    IF ( nbrmpireqcnt .GT. 0 ) THEN
      IF(KPDBG) WRITE(OFU,*) ' Wait ',nbrmpinirecvs,' ',nbrmpinisends; CALL FL(OFU)
      CALL BSystemClock(loctimer1)
      CALL MPI_Waitall( nbrmpireqcnt, nbrmpireq, nbrmpistatuses, mpiwaiterr )
      CALL BSystemClock(loctimer2)
      CALL ChargeTime( totcommtime, loctimer2, loctimer1, totcommcount )
    END IF

    !---------------------------------------------------------------------------
    !Now we can compute inverses for odd boundaries, if any
    !---------------------------------------------------------------------------
    IF ( ISODD(startlocrow) ) THEN
      globrow = LR2GR( startlocrow, level )
      IF(KPDBG) WRITE(OFU,*) ' Postcomp mat-mul ',globrow,' ',startlocrow;CALL FL(OFU)
      CALL ComputeForwardOddRowHats( startlocrow, level, &
                                     startlocrow, endlocrow, .FALSE. )
    END IF
    IF ( (startlocrow .NE. endlocrow) .AND. ISODD(endlocrow) ) THEN
      globrow = LR2GR( endlocrow, level )
      IF(KPDBG) WRITE(OFU,*) ' Postcomp mat-mul ',globrow,' ',endlocrow;CALL FL(OFU)
      CALL ComputeForwardOddRowHats( endlocrow, level, &
                                     startlocrow, endlocrow, .FALSE. )
    END IF
    CALL PLBForwardFinalizeLevel( level, .TRUE. )
  END DO ForwardSolveLoop

  !-----------------------------------------------------------------------------
  CALL PLBForwardFinalize()

  !-----------------------------------------------------------------------------
  matdirtied = .FALSE.
  rhsdirtied = .FALSE.

  !-----------------------------------------------------------------------------
  IF(KPDBG) WRITE(OFU,'(A,F6.1,A)') 'Memory so far: ', membytes/1e6, ' MB'; CALL FL(OFU)
  IF(KPDBG) WRITE(OFU,*) '------ Forward solve end ------'; CALL FL(OFU)

  CALL BSystemClock(globtimer2)
  CALL ChargeTime( tottime, globtimer2, globtimer1, totcount )

END SUBROUTINE ForwardSolve

!-------------------------------------------------------------------------------
!>
!!BCYCLIC forward phase to deal with a new b; to be called after a SetRHS that
!!may have been invoked after a ForwardSolve, before BackwardSolve.
!<
!-------------------------------------------------------------------------------
SUBROUTINE ForwardUpdateb
  !-----------------------------------------------------------------------------
  ! Local variables
  !-----------------------------------------------------------------------------
  INTEGER :: level !<Loop variable
  INTEGER :: locrow !<Loop variable for "local" row number
  INTEGER :: globrow !<Loop variable for "global" row number
  INTEGER :: globrowoff !<Temporary variable
  INTEGER :: startlocrow !<Starting row number of this processor at a level
  INTEGER :: endlocrow !<Ending row number (inclusive) at a given level
  INTEGER :: rowoffset !<Loop variable to generate row+/-1
  INTEGER :: nbrrank !<Temporary variable
  INTEGER :: nbrmpireqindx !<Index of next MPI_Irecv request
  INTEGER :: nbrmpireqcnt !<Number of MPI_Irecv/Isend requests initiated
  INTEGER :: nbrmpinirecvs !<Number of MPI_Irecv requests initiated
  INTEGER :: nbrmpinisends !<Number of MPI_Isend requests initiated
  INTEGER :: nbrmpireq(6) !<MPI_Irecv/Isend request handle for top or bot nbr
  INTEGER :: nbrmpierr(6) !<MPI_Irecv/Isend error for top or bot nbr
  INTEGER :: mpiwaiterr !<MPI_Waitall/any error
  INTEGER :: mpierr !<Generic use for MPI
  INTEGER :: nbrmpistatuses(MPI_STATUS_SIZE,6) !<MPI_Waitall statuses
  INTEGER :: nbrmpistatus(MPI_STATUS_SIZE) !<MPI_Waitany status
  INTEGER :: waitmatchindex !<MPI_Waitany match index
  INTEGER :: stag !<Tag number for Isends
  INTEGER :: reqi !<Looping variable for requests
  INTEGER :: blaserr !<Generic use for BLAS

  !-----------------------------------------------------------------------------
  IF(KPDBG) WRITE(OFU,*); CALL FL(OFU)
  IF(KPDBG) WRITE(OFU,*) '------ Forward updateb start ------'; CALL FL(OFU)
  IF(KPDBG) WRITE(OFU,*); CALL FL(OFU)

  !-----------------------------------------------------------------------------
  IF ( matdirtied ) THEN
    IF(KPDBG) WRITE(OFU,*); CALL FL(OFU)
    IF(KPDBG) WRITE(OFU,*) '  ------ Dirtied matrix ------'; CALL FL(OFU)
    IF(KPDBG) WRITE(OFU,*) '  ------ Forward solve, not updateb... ------'; CALL FL(OFU)
    CALL ForwardSolve()
    RETURN
  END IF

  CALL BSystemClock(globtimer1)

  !-----------------------------------------------------------------------------
  CALL PLBForwardInitialize()

  !-----------------------------------------------------------------------------
  ! Forward solving loop in which even rows have already computed inverses and
  ! odd rows recompute matrix-vector products with new rhs
  !-----------------------------------------------------------------------------
  ForwardSolveLoop: DO level = 1, nlevels, 1
    !---------------------------------------------------------------------------
    IF ( usebarriers ) THEN
      IF(KPDBG) WRITE(OFU,*) 'Barrier at forward updateb level ', level; CALL FL(OFU)
      CALL MPI_Barrier( SIESTA_COMM, mpierr )
      IF(KPDBG) WRITE(OFU,*) 'Done barrier at forward updateb level ', level; CALL FL(OFU)
    END IF

    !---------------------------------------------------------------------------
    !Determine start and end local rows at current level
    !---------------------------------------------------------------------------
    startlocrow = N+1 !Start with an invalid value (greater than endlocrow)
    endlocrow = -1 !Start with an invalid value (less than startlocrow)
    DO globrow = startglobrow, endglobrow, 1
      locrow = GR2LR( globrow, level ) !globrow's place at current level
      IF ( locrow .GT. 0 ) THEN !This global row participates at current level
        IF ( locrow .LT. startlocrow ) THEN !A valid, earlier local row
          startlocrow = locrow
        END IF
        IF ( locrow .GT. endlocrow ) THEN !A valid, later local row
          endlocrow = locrow
        END IF
      END IF
    END DO

    !---------------------------------------------------------------------------
    !We may have run out of work; see if we have a valid range of rows left
    !---------------------------------------------------------------------------
    IF ( startlocrow .GT. endlocrow ) THEN !No rows left at/above this level
      IF(KPDBG) WRITE(OFU,*) 'Nothing to do at forward updateb level ',level; CALL FL(OFU)
      CALL PLBForwardInitializeLevel( level, .FALSE. )
      CALL PLBForwardFinalizeLevel( level, .FALSE. )
      CYCLE ForwardSolveLoop !Move on to next level; don't EXIT (for barriers)
    ELSE
      CALL PLBForwardInitializeLevel( level, .TRUE. )
    END IF

    IF(KPDBG) WRITE(OFU,*) '**** Forward updateb level ', level, ' ****'; CALL FL(OFU)

    !---------------------------------------------------------------------------
    !Memory has already been allocated in ForwardSolve for incoming nbr values,
    !if any, at this level
    !---------------------------------------------------------------------------
    IF ( ISODD(startlocrow) .AND. &
     (LR2Rank(startlocrow-1,level) .GE. 0) ) THEN
      globrowoff = 0
      IF( .NOT. ( ALLOCATED( lelement(level, globrowoff)%L ) .AND. &
                  ALLOCATED( lelement(level, globrowoff)%U ) .AND. &
                  ALLOCATED( lelement(level, globrowoff)%b ) ) ) THEN
        IF(KPDBG) WRITE(OFU,*)'ForwardUpdateb +0 memory error'; CALL FL(OFU); CALL ASSERT(.FALSE.,'blocktri #79')
      END IF
      lelement(level, globrowoff)%b = 0
    END IF
    IF ( ISODD(endlocrow) .AND. &
     (LR2Rank(endlocrow+1,level) .GE. 0) ) THEN
      globrowoff = endglobrow-startglobrow+2
      IF ( .NOT. ( ALLOCATED( lelement(level, globrowoff)%L ) .AND. &
                   ALLOCATED( lelement(level, globrowoff)%U ) .AND. &
                   ALLOCATED( lelement(level, globrowoff)%b ) ) ) THEN
        IF(KPDBG) WRITE(OFU,*)'ForwardUpdateb +2 memory error'; CALL FL(OFU); CALL ASSERT(.FALSE.,'blocktri #80')
      END IF
      lelement(level, globrowoff)%b = 0
    END IF

    !---------------------------------------------------------------------------
    !Memory has already been allocated in ForwardSolve at next level for those
    !rows that are odd at current level
    !---------------------------------------------------------------------------
    IF ( level .LT. nlevels ) THEN
      DO locrow = startlocrow, endlocrow, 1
        IF ( ISODD(locrow) ) THEN
          globrow = LR2GR( locrow, level )
          globrowoff = globrow-startglobrow+1
          IF ( .NOT. ( ALLOCATED( lelement(level+1, globrowoff)%L ) .AND. &
                       ALLOCATED( lelement(level+1, globrowoff)%D ) .AND. &
                       ALLOCATED( lelement(level+1, globrowoff)%U ) .AND. &
                       ALLOCATED( lelement(level+1, globrowoff)%b ) ) ) THEN
            IF(KPDBG) WRITE(OFU,*)'ForwardUpdateb +1 memory error'; CALL FL(OFU); CALL ASSERT(.FALSE.,'blocktri #81')
          END IF
          lelement(level+1, globrowoff)%b = 0
        END IF
      END DO
    END IF

    !---------------------------------------------------------------------------
    !Reset requests
    !---------------------------------------------------------------------------
    DO nbrmpireqindx = 1, 6
      nbrmpireq(nbrmpireqindx) = MPI_REQUEST_NULL
    END DO

    !---------------------------------------------------------------------------
    !Pre-post the expected MPI receives via a non-blocking recv primitive
    !Pre-posting can help overlap communication with computing inverses
    !---------------------------------------------------------------------------
    nbrmpireqindx = 1
    nbrmpinirecvs = 0
    nbrrank = LR2Rank(startlocrow-1,level)
    IF ( ISODD(startlocrow) .AND. (nbrrank .GE. 0) ) THEN
      !-------------------------------------------------------------------------
      !Our top row at this level is odd and top row's previous is valid
      !We will get the previous (even-numbered) row's L-hat, U-hat, and b-hat
      !-------------------------------------------------------------------------
      globrowoff = 0
      IF(KPDBG) WRITE(OFU,*) '  Irecv ',startlocrow-1,' ',nbrrank,' '; CALL FL(OFU)

      IF ( .NOT. ALLOCATED(lelement(level,globrowoff)%L) ) THEN
        IF(KPDBG) WRITE(OFU,*)'  Irecv bad L'; CALL FL(OFU); CALL ASSERT(.FALSE.,'blocktri #82')
      END IF
      IF ( .NOT. ALLOCATED(lelement(level,globrowoff)%U) ) THEN
        IF(KPDBG) WRITE(OFU,*)'  Irecv bad U'; CALL FL(OFU); CALL ASSERT(.FALSE.,'blocktri #83')
      END IF
      IF ( .NOT. ALLOCATED(lelement(level,globrowoff)%b) ) THEN
        IF(KPDBG) WRITE(OFU,*)'  Irecv bad b'; CALL FL(OFU); CALL ASSERT(.FALSE.,'blocktri #84')
      END IF

      CALL BSystemClock(loctimer1)
      CALL MPI_Irecv( lelement(level, globrowoff)%b, M, MPI_REAL8, &
       nbrrank, 3, SIESTA_COMM, nbrmpireq(nbrmpireqindx), &
       nbrmpierr(nbrmpireqindx))
      nbrmpireqindx = nbrmpireqindx + 1
      nbrmpinirecvs = nbrmpinirecvs + 1
      IF(KPDBG) WRITE(OFU,*) '  Irecv 3 top ',nbrrank,' ',nbrmpireqindx-1; CALL FL(OFU)
      CALL BSystemClock(loctimer2)
      CALL ChargeTime( totcommtime, loctimer2, loctimer1, totcommcount )
    END IF

    nbrrank = LR2Rank(endlocrow+1,level)
    IF ( ISODD(endlocrow) .AND. (nbrrank .GE. 0) ) THEN
      !-------------------------------------------------------------------------
      !Our bottom row at this level is odd and bottom row's next is valid
      !We will get the next (even-numbered) row's L-hat, U-hat, and b-hat
      !-------------------------------------------------------------------------
      globrowoff = endglobrow-startglobrow+2
      IF(KPDBG) WRITE(OFU,*) ' Irecv ',endlocrow+1,' ', nbrrank; CALL FL(OFU)

      IF ( .NOT. ALLOCATED(lelement(level,globrowoff)%L) ) THEN
        IF(KPDBG) WRITE(OFU,*)'Irecv bad L'; CALL FL(OFU); CALL ASSERT(.FALSE.,'blocktri #85')
      END IF
      IF ( .NOT. ALLOCATED(lelement(level,globrowoff)%U) ) THEN
        IF(KPDBG) WRITE(OFU,*)'Irecv bad U'; CALL FL(OFU); CALL ASSERT(.FALSE.,'blocktri #86')
      END IF
      IF ( .NOT. ALLOCATED(lelement(level,globrowoff)%b) ) THEN
        IF(KPDBG) WRITE(OFU,*)'Irecv bad b'; CALL FL(OFU); CALL ASSERT(.FALSE.,'blocktri #87')
      END IF

      CALL BSystemClock(loctimer1)
      CALL MPI_Irecv( lelement(level, globrowoff)%b, M, MPI_REAL8, &
       nbrrank, 6, SIESTA_COMM, nbrmpireq(nbrmpireqindx), &
       nbrmpierr(nbrmpireqindx))
      nbrmpireqindx = nbrmpireqindx + 1
      nbrmpinirecvs = nbrmpinirecvs + 1
      IF(KPDBG) WRITE(OFU,*) '  Irecv 6 bot ',nbrrank,' ',nbrmpireqindx-1; CALL FL(OFU)
      CALL BSystemClock(loctimer2)
      CALL ChargeTime( totcommtime, loctimer2, loctimer1, totcommcount )
    END IF

    !---------------------------------------------------------------------------
    !Inverses have already been computed for even-numbered rows at this level
    !---------------------------------------------------------------------------
    nbrmpinisends = 0
    DO locrow = startlocrow, endlocrow, 1
      IF ( ISEVEN( locrow ) .OR. (level .EQ. nlevels) ) THEN
        globrow = LR2GR( locrow, level )
        globrowoff = globrow-startglobrow+1

        IF(KPDBG) WRITE(OFU,*) ' Computed even hats ',globrow,' ', locrow; CALL FL(OFU)

        !-----------------------------------------------------------------------
        !Sanity checks
        !-----------------------------------------------------------------------
        IF ( level .EQ. nlevels ) THEN
          IF ( (rank .NE. 0) .OR. &
               (startlocrow .NE. endlocrow) .OR. &
               (locrow .NE. 1) .OR. (globrow .NE. 1) ) THEN
            IF(KPDBG) WRITE(OFU,*) ' EVEN ERROR ',globrow,' ', locrow; CALL FL(OFU)
            CALL ASSERT(.FALSE.,'blocktri #88')
          END IF
        END IF
        IF ( .NOT. ALLOCATED(lelement(level,globrowoff)%D) ) THEN
          IF(KPDBG) WRITE(OFU,*)'Chats bad D'; CALL FL(OFU); CALL ASSERT(.FALSE.,'blocktri #89')
        END IF
        IF ( .NOT. ALLOCATED(lelement(1,globrowoff)%pivot) ) THEN
          IF(KPDBG) WRITE(OFU,*)'Chats bad pivot'; CALL FL(OFU); CALL ASSERT(.FALSE.,'blocktri #90')
        END IF
        IF ( .NOT. ALLOCATED(lelement(level,globrowoff)%L) ) THEN
          IF(KPDBG) WRITE(OFU,*)'Chats bad L'; CALL FL(OFU); CALL ASSERT(.FALSE.,'blocktri #91')
        END IF
        IF ( .NOT. ALLOCATED(lelement(level,globrowoff)%U) ) THEN
          IF(KPDBG) WRITE(OFU,*)'Chats bad U'; CALL FL(OFU); CALL ASSERT(.FALSE.,'blocktri #92')
        END IF
        IF ( .NOT. ALLOCATED(lelement(level,globrowoff)%b) ) THEN
          IF(KPDBG) WRITE(OFU,*)'Chats bad b'; CALL FL(OFU); CALL ASSERT(.FALSE.,'blocktri #93')
        END IF

        !-----------------------------------------------------------------------
        !LU factorization of D has already been performed in place into D itself
        !-----------------------------------------------------------------------

        !-----------------------------------------------------------------------
        !Compute b-hats D^(-1)b
        !-----------------------------------------------------------------------
        CALL PLBDGETRS( 1, lelement(level,globrowoff)%D, &
                           lelement(1,globrowoff)%pivot, &
                           lelement(level,globrowoff)%b, blaserr )
        IF (blaserr .NE. 0) THEN
          IF(KPDBG) WRITE(OFU,*) 'PLBDGETRS b failed'; CALL FL(OFU); CALL ASSERT(.FALSE.,'blocktri #94')
        END IF

        !-----------------------------------------------------------------------
        !Send my bhats to neighbor row if that neighbor exists
        !and that neighbor row is hosted outside this rank
        !-----------------------------------------------------------------------
        DO rowoffset = -1, 1, 2
          nbrrank = LR2Rank(locrow+rowoffset, level)
          IF ( (nbrrank .GE. 0) .AND. (nbrrank .NE. rank) ) THEN
            CALL BSystemClock(loctimer1)
            IF ( .NOT. ALLOCATED(lelement(level,globrowoff)%L) ) THEN
              IF(KPDBG) WRITE(OFU,*)'ISend bad L'; CALL FL(OFU); CALL ASSERT(.FALSE.,'blocktri #95')
            END IF
            IF ( .NOT. ALLOCATED(lelement(level,globrowoff)%U) ) THEN
              IF(KPDBG) WRITE(OFU,*)'ISend bad U'; CALL FL(OFU); CALL ASSERT(.FALSE.,'blocktri #96')
            END IF
            IF ( .NOT. ALLOCATED(lelement(level,globrowoff)%b) ) THEN
              IF(KPDBG) WRITE(OFU,*)'ISend bad b'; CALL FL(OFU); CALL ASSERT(.FALSE.,'blocktri #97')
            END IF

            IF(KPDBG) WRITE(OFU,*)' ISend ',globrow,' ',locrow,' ',nbrrank;CALL FL(OFU)
            globrowoff = globrow-startglobrow+1
            stag = (((1-rowoffset))/2)*3
            CALL MPI_Isend( lelement(level,globrowoff)%b, M, MPI_REAL8, &
             nbrrank, stag+3, SIESTA_COMM, &
             nbrmpireq(nbrmpireqindx), nbrmpierr(nbrmpireqindx))
            nbrmpireqindx = nbrmpireqindx + 1
            nbrmpinisends = nbrmpinisends + 1
            IF(KPDBG) WRITE(OFU,*) '  Isend ',stag+3,' ',globrowoff; CALL FL(OFU)
            CALL BSystemClock(loctimer2)
            CALL ChargeTime( totcommtime, loctimer2, loctimer1, totcommcount )
          END IF
        END DO
      END IF
    END DO

    !---------------------------------------------------------------------------
    !Compute the matrix-vector multiplications for non-boundary odd rows
    !---------------------------------------------------------------------------
    DO locrow = startlocrow, endlocrow, 1
      globrow = LR2GR( locrow, level )
      IF ( ISODD( locrow ) .AND. (level .NE. nlevels) ) THEN
        IF ( (locrow .NE. startlocrow) .AND. (locrow .NE. endlocrow) ) THEN
          IF(KPDBG) WRITE(OFU,*) ' Precomp mat-mul ',globrow,' ',locrow; CALL FL(OFU)
          CALL ComputeForwardOddRowHats( locrow, level, &
                                         startlocrow, endlocrow, .TRUE. )
        END IF
      END IF
    END DO

    !---------------------------------------------------------------------------
    !We have to wait for incoming even-numbered neighboring rows, before
    !computing inverses for the odd boundaries, if any.
    !Complete the send-receive of inverses, by doing a "wait all"
    !---------------------------------------------------------------------------
    nbrmpireqcnt = nbrmpireqindx - 1
    IF ( nbrmpireqcnt .GT. 0 ) THEN
      IF(KPDBG) WRITE(OFU,*) ' Wait ',nbrmpinirecvs,' ',nbrmpinisends; CALL FL(OFU)
      CALL BSystemClock(loctimer1)
      CALL MPI_Waitall( nbrmpireqcnt, nbrmpireq, nbrmpistatuses, mpiwaiterr )
      CALL BSystemClock(loctimer2)
      CALL ChargeTime( totcommtime, loctimer2, loctimer1, totcommcount )
    END IF

    !---------------------------------------------------------------------------
    !Now we can compute inverses for odd boundaries, if any
    !---------------------------------------------------------------------------
    IF ( ISODD(startlocrow) ) THEN
      globrow = LR2GR( startlocrow, level )
      IF(KPDBG) WRITE(OFU,*) ' Postcomp mat-mul ',globrow,' ',startlocrow;CALL FL(OFU)
      CALL ComputeForwardOddRowHats( startlocrow, level, &
                                     startlocrow, endlocrow, .TRUE. )
    END IF
    IF ( (startlocrow .NE. endlocrow) .AND. ISODD(endlocrow) ) THEN
      globrow = LR2GR( endlocrow, level )
      IF(KPDBG) WRITE(OFU,*) ' Postcomp mat-mul ',globrow,' ',endlocrow;CALL FL(OFU)
      CALL ComputeForwardOddRowHats( endlocrow, level, &
                                     startlocrow, endlocrow, .TRUE. )
    END IF
    CALL PLBForwardFinalizeLevel( level, .TRUE. )
  END DO ForwardSolveLoop

  !-----------------------------------------------------------------------------
  CALL PLBForwardFinalize()

  !-----------------------------------------------------------------------------
  rhsdirtied = .FALSE.

  !-----------------------------------------------------------------------------
  IF(KPDBG) WRITE(OFU,'(A,F6.1,A)') 'Memory so far: ', membytes/1e6, ' MB'; CALL FL(OFU)
  IF(KPDBG) WRITE(OFU,*) '------ updateb solve end ------'; CALL FL(OFU)

  CALL BSystemClock(globtimer2)
  CALL ChargeTime( tottime, globtimer2, globtimer1, totcount )

END SUBROUTINE ForwardUpdateb

!-------------------------------------------------------------------------------
!>
!!BCYCLIC backward phase; to be called after ForwardSolve
!<
!-------------------------------------------------------------------------------
SUBROUTINE BackwardSolve
  !-----------------------------------------------------------------------------
  ! Local variables
  !-----------------------------------------------------------------------------
  INTEGER :: level !<Loop variable
  INTEGER :: locrow !<Loop variable for "local" row number
  INTEGER :: globrow !<Loop variable for "global" row number
  INTEGER :: globrowoff !<Temporary variable
  INTEGER :: prevglobrow !<Temporary variable
  INTEGER :: nextglobrow !<Temporary variable
  INTEGER :: startlocrow !<Starting row number of this processor at a level
  INTEGER :: endlocrow !<Ending row number (inclusive) at a given level
  INTEGER :: rowoffset !<Loop variable to generate row+/-1
  INTEGER :: nbrrank !<Temporary variable
  INTEGER :: nbrmpireqindx !<Index of next MPI_Irecv request
  INTEGER :: nbrmpireqcnt !<Number of MPI_Irecv/Isend requests initiated
  INTEGER :: nbrmpinirecvs !<Number of MPI_Irecv requests initiated
  INTEGER :: nbrmpinisends !<Number of MPI_Isend requests initiated
  INTEGER :: nbrmpireq(2) !<MPI_Irecv/Isend request handle for top or bot nbr
  INTEGER :: nbrmpierr(2) !<MPI_Irecv/Isend error for top or bot nbr
  INTEGER :: mpiwaiterr !<MPI_Waitall/any error
  INTEGER :: mpierr !<Generic use for MPI
  INTEGER :: nbrmpistatuses(MPI_STATUS_SIZE,2) !<MPI_Waitall statuses
  INTEGER :: nbrmpistatus(MPI_STATUS_SIZE) !<MPI_Waitany status
  INTEGER :: waitmatchindex !<MPI_Waitany match index
  INTEGER :: stag !<Tag number for Isends
  INTEGER :: reqi !<Looping variable for requests
  INTEGER :: blaserr !<Generic use for BLAS

  CALL BSystemClock(globtimer1)

  !-----------------------------------------------------------------------------
  IF(KPDBG) WRITE(OFU,*); CALL FL(OFU)
  IF(KPDBG) WRITE(OFU,*) '------ Backward solve start ------'; CALL FL(OFU)

  !-----------------------------------------------------------------------------
  IF ( matdirtied ) THEN
    IF(KPDBG) WRITE(OFU,*); CALL FL(OFU)
    IF(KPDBG) WRITE(OFU,*) '  ------ Dirtied matrix; updating... ------'; CALL FL(OFU)
    CALL ForwardSolve()
  END IF
  IF ( rhsdirtied ) THEN
    IF(KPDBG) WRITE(OFU,*); CALL FL(OFU)
    IF(KPDBG) WRITE(OFU,*) '  ------ Dirtied RHS; updating... ------'; CALL FL(OFU)
    CALL ForwardUpdateb()
  END IF

  !-----------------------------------------------------------------------------
  CALL PLBBackwardInitialize()

  !-----------------------------------------------------------------------------
  ! Make sure we have the forward solve data
  !-----------------------------------------------------------------------------
  IF ( (.NOT. ALLOCATED(lelement)) .OR. (.NOT. ALLOCATED(selement)) ) THEN
    IF(KPDBG) WRITE(OFU,*) 'Forward not called before backward?'; CALL FL(OFU)
    RETURN
  END IF

  !-----------------------------------------------------------------------------
  ! Backward solving loop in which odd rows communicate their solution to even
  ! rows and even rows compute their solution
  !-----------------------------------------------------------------------------
  BackwardSolveLoop: DO level = nlevels, 1, -1
    !---------------------------------------------------------------------------
    IF ( usebarriers ) THEN
      IF(KPDBG) WRITE(OFU,*) 'Barrier at backward level ', level; CALL FL(OFU)
      CALL MPI_Barrier( SIESTA_COMM, mpierr )
      IF(KPDBG) WRITE(OFU,*) 'Done barrier at backward level ', level; CALL FL(OFU)
    END IF

    !---------------------------------------------------------------------------
    !Determine start and end local rows at current level
    !---------------------------------------------------------------------------
    startlocrow = N+1 !Start with an invalid value (greater than endlocrow)
    endlocrow = -1 !Start with an invalid value (less than startlocrow)
    DO globrow = startglobrow, endglobrow, 1
      locrow = GR2LR( globrow, level ) !globrow's place at current level
      IF ( locrow .GT. 0 ) THEN !This global row participates at current level
        IF ( locrow .LT. startlocrow ) THEN !A valid, earlier local row
          startlocrow = locrow
        END IF
        IF ( locrow .GT. endlocrow ) THEN !A valid, later local row
          endlocrow = locrow
        END IF
      END IF
    END DO

    !---------------------------------------------------------------------------
    !We may have nothing to do at this level; if so, just go back another level
    !---------------------------------------------------------------------------
    IF ( startlocrow .GT. endlocrow ) THEN !No more rows left at this level
      IF(KPDBG) WRITE(OFU,*) 'Nothing to do at backward level', level; CALL FL(OFU)
      CALL PLBBackwardInitializeLevel( level, .FALSE. )
      CALL PLBBackwardFinalizeLevel( level, .FALSE. )
      CYCLE BackwardSolveLoop !Move on to previous level
    ELSE
      CALL PLBBackwardInitializeLevel( level, .TRUE. )
    END IF

    IF(KPDBG) WRITE(OFU,*) '**** Backward level ', level, ' ****'; CALL FL(OFU)

    !---------------------------------------------------------------------------
    !Reset requests
    !---------------------------------------------------------------------------
    DO nbrmpireqindx = 1, 2
      nbrmpireq(nbrmpireqindx) = MPI_REQUEST_NULL
    END DO

    !---------------------------------------------------------------------------
    !Pre-post the expected MPI receives via a non-blocking recv primitive
    !Pre-posting can help overlap communication with computing inverses
    !---------------------------------------------------------------------------
    nbrmpireqindx = 1
    nbrmpinirecvs = 0
    nbrrank = LR2Rank(startlocrow-1,level)
    IF( ISEVEN(startlocrow) .AND. (nbrrank .GE. 0) ) THEN
      !-------------------------------------------------------------------------
      !Our top row at this level is even and top row's previous is valid
      !We will get the previous (odd-numbered) row's solution
      !-------------------------------------------------------------------------
      stag = 1
      globrow = LR2GR(startlocrow-1,level)
     IF(KPDBG) WRITE(OFU,*)' Irecv ',startlocrow-1,' ',globrow,' ',nbrrank;CALL FL(OFU)
      IF ( .NOT. ALLOCATED( selement(globrow)%x ) ) THEN
        ALLOCATE( selement(globrow)%x(M) )
        CALL ChargeMemory( M*dpsz )
      END IF
      CALL BSystemClock(loctimer1)
      CALL MPI_Irecv( selement(globrow)%x, M, MPI_REAL8, nbrrank, stag, &
       SIESTA_COMM, nbrmpireq(nbrmpireqindx), nbrmpierr(nbrmpireqindx) )
      nbrmpireqindx = nbrmpireqindx + 1
      nbrmpinirecvs = nbrmpinirecvs + 1
      CALL BSystemClock(loctimer2)
      CALL ChargeTime( totcommtime, loctimer2, loctimer1, totcommcount )
    END IF
    nbrrank = LR2Rank(endlocrow+1,level)
    IF ( ISEVEN(endlocrow) .AND. (nbrrank .GE. 0) ) THEN
      !-------------------------------------------------------------------------
      !Our bottom row at this level is even and bottom row's next is valid
      !We will get the next (odd-numbered) row's solution
      !-------------------------------------------------------------------------
      stag = 2
      globrow = LR2GR(endlocrow+1,level)
      IF(KPDBG) WRITE(OFU,*)' Irecv ',endlocrow+1,' ',globrow,' ',nbrrank;CALL FL(OFU)
      IF ( .NOT. ALLOCATED( selement(globrow)%x ) ) THEN
        ALLOCATE( selement(globrow)%x(M) )
        CALL ChargeMemory( M*dpsz )
      END IF
      CALL BSystemClock(loctimer1)
      CALL MPI_Irecv( selement(globrow)%x, M, MPI_REAL8, nbrrank, stag, &
       SIESTA_COMM, nbrmpireq(nbrmpireqindx), nbrmpierr(nbrmpireqindx) )
      nbrmpireqindx = nbrmpireqindx + 1
      nbrmpinirecvs = nbrmpinirecvs + 1
      CALL BSystemClock(loctimer2)
      CALL ChargeTime( totcommtime, loctimer2, loctimer1, totcommcount )
    END IF

    !---------------------------------------------------------------------------
    !Send our odd-numbered rows' solutions
    !---------------------------------------------------------------------------
    nbrmpinisends = 0
    DO locrow = startlocrow, endlocrow, 1
      IF ( ISODD( locrow ) ) THEN
        !-----------------------------------------------------------------------
        !Send my solution to neighbor row if that neighbor exists
        !and is hosted outside this rank; use non-blocking sends
        !-----------------------------------------------------------------------
        DO rowoffset = -1, 1, 2
          nbrrank = LR2Rank(locrow+rowoffset, level)
          IF ( (nbrrank .GE. 0) .AND. (nbrrank .NE. rank) ) THEN
            IF(KPDBG) WRITE(OFU,*) ' Isend ', locrow, ' ', nbrrank; CALL FL(OFU)
            globrow = LR2GR(locrow,level)
            IF ( .NOT. ALLOCATED(selement(globrow)%x) ) THEN
              IF(KPDBG) WRITE(OFU,*) 'ERROR BISEND AT LEVEL', level; CALL FL(OFU)
              IF(KPDBG) WRITE(OFU,*) locrow, ' ', globrow, ' ', nbrrank; CALL FL(OFU)
              CALL ASSERT(.FALSE.,'blocktri #98')
            END IF
            stag = (1-rowoffset)/2+1
            CALL BSystemClock(loctimer1)
            CALL MPI_Isend( selement(globrow)%x, M, MPI_REAL8, &
             nbrrank, stag, SIESTA_COMM, nbrmpireq(nbrmpireqindx), &
             nbrmpierr(nbrmpireqindx) )
            nbrmpireqindx = nbrmpireqindx + 1
            nbrmpinisends = nbrmpinisends + 1
            CALL BSystemClock(loctimer2)
            CALL ChargeTime( totcommtime, loctimer2, loctimer1, totcommcount )
          END IF
        END DO
      END IF
    END DO

    !---------------------------------------------------------------------------
    !Complete the send-receive of solutions, by doing an MPI "wait all"
    !---------------------------------------------------------------------------
    nbrmpireqcnt = nbrmpireqindx - 1
    IF ( nbrmpireqcnt .GT. 0 ) THEN
      IF(KPDBG) WRITE(OFU,*) ' Wait ',nbrmpinirecvs,' ',nbrmpinisends; CALL FL(OFU)
      CALL BSystemClock(loctimer1)
      CALL MPI_Waitall( nbrmpireqcnt, nbrmpireq, nbrmpistatuses, mpiwaiterr )
      CALL BSystemClock(loctimer2)
      CALL ChargeTime( totcommtime, loctimer2, loctimer1, totcommcount )
    END IF

    !---------------------------------------------------------------------------
    !Compute the solution for even rows at this level
    !---------------------------------------------------------------------------
    DO locrow = startlocrow, endlocrow, 1
      IF ( ISEVEN( locrow ) .OR. (level .EQ. nlevels) ) THEN
        globrow = LR2GR( locrow, level )
        globrowoff = globrow - startglobrow + 1
        IF ( ( level .EQ. nlevels ) .AND. (locrow .NE. 1) ) THEN
          IF(KPDBG) WRITE(OFU,*) 'ERROR ',level, ' ',globrow, ' ', locrow; CALL FL(OFU)
          CALL ASSERT(.FALSE.,'blocktri #99')
        END IF

        IF(KPDBG) WRITE(OFU,*) ' Compute solution ', globrow, ' ', locrow; CALL FL(OFU)

        IF ( .NOT. ALLOCATED( selement(globrow)%x ) ) THEN
          ALLOCATE( selement(globrow)%x(M) )
          CALL ChargeMemory( M*dpsz )
        END IF

        !----------------------------------------------------
        ! x = b-hat
        !----------------------------------------------------
        selement(globrow)%x = lelement(level,globrowoff)%b(:,1)

        !----------------------------------------------------
        ! x = x - L-hat*x(prev)
        !----------------------------------------------------
        nbrrank = LR2Rank(locrow-1,level)
        IF ( nbrrank .GE. 0 ) THEN
          prevglobrow = LR2GR(locrow-1,level)
          CALL PLBDGEMV( -ONE, lelement(level,globrowoff)%L, &
                            selement(prevglobrow)%x, &
                            ONE, selement(globrow)%x )
        END IF

        !----------------------------------------------------
        ! x = x - U-hat*x(next)
        !----------------------------------------------------
        nbrrank = LR2Rank(locrow+1,level)
        IF ( nbrrank .GE. 0 ) THEN
          nextglobrow = LR2GR(locrow+1,level)
          CALL PLBDGEMV( -ONE, lelement(level,globrowoff)%U, &
                            selement(nextglobrow)%x, &
                            ONE, selement(globrow)%x )
        END IF

        !----------------------------------------------------
        !Print the solution vector, if desired
        !----------------------------------------------------
        IF ( .FALSE. ) THEN
          IF(KPDBG) WRITE(OFU,*) ' x[',globrow,']=',selement(globrow)%x; CALL FL(OFU)
        END IF
      END IF
    END DO
    CALL PLBBackwardFinalizeLevel( level, .FALSE. )
  END DO BackwardSolveLoop

  !-----------------------------------------------------------------------------
  CALL PLBBackwardFinalize()

  !-----------------------------------------------------------------------------
  IF(KPDBG) WRITE(OFU,'(A,F6.1,A)') 'Memory so far: ', membytes/1e6, ' MB'; CALL FL(OFU)
  IF(KPDBG) WRITE(OFU,*) '------ Backward solve end ------'; CALL FL(OFU)

  CALL BSystemClock(globtimer2)
  CALL ChargeTime( tottime, globtimer2, globtimer1, totcount )

END SUBROUTINE BackwardSolve

!-------------------------------------------------------------------------------
!>
!!Verify the RMS error of solution after backward solve
!<
!-------------------------------------------------------------------------------
SUBROUTINE VerifySolution
  !-----------------------------------------------------------------------------
  ! Local variables
  !-----------------------------------------------------------------------------
  INTEGER :: i, k, globrow, globrowoff
  REAL(dp) :: term, totrmserr
  INTEGER :: nbrrank !<Rank of a neighboring row
  INTEGER :: nbrmpireqindx !<Index of next MPI_Irecv request
  INTEGER :: nbrmpireqcnt !<Number of MPI_Irecv/Isend requests initiated
  INTEGER :: nbrmpinirecvs !<Number of MPI_Irecv requests initiated
  INTEGER :: nbrmpinisends !<Number of MPI_Isend requests initiated
  INTEGER :: nbrmpireq(2) !<MPI_Irecv/Isend request handle for top or bot nbr
  INTEGER :: nbrmpierr(2) !<MPI_Irecv/Isend error for top or bot nbr
  INTEGER :: mpiwaiterr !<MPI_Waitall/any error
  INTEGER :: mpierr !<Generic use for MPI
  INTEGER :: nbrmpistatuses(MPI_STATUS_SIZE,2) !<MPI_Waitall statuses
  INTEGER :: nbrmpistatus(MPI_STATUS_SIZE) !<MPI_Waitany status
  INTEGER :: waitmatchindex !<MPI_Waitany match index

  !-----------------------------------------------------------------------------
  IF(KPDBG) WRITE(OFU,*); CALL FL(OFU)
  IF(KPDBG) WRITE(OFU,*) '------ Verifying solution ------'; CALL FL(OFU)

  !-----------------------------------------------------------------------------
  !Alocate memory, do Irecvs/Isends, for boundary odd/even rows resp.
  !-----------------------------------------------------------------------------
  nbrmpireqindx = 1
  nbrmpinirecvs = 0
  nbrmpinisends = 0
  nbrrank = GR2Rank(startglobrow-1)
  IF ( ISODD(startglobrow) .AND. (nbrrank .GE. 0) ) THEN
    IF ( .NOT. ALLOCATED( selement(startglobrow-1)%x ) ) THEN
        ALLOCATE( selement(startglobrow-1)%x(M) )
    END IF
    CALL MPI_Irecv( selement(startglobrow-1)%x, M, MPI_REAL8, nbrrank, 1, &
      SIESTA_COMM, nbrmpireq(nbrmpireqindx), nbrmpierr(nbrmpireqindx) )
    nbrmpireqindx = nbrmpireqindx + 1
    nbrmpinirecvs = nbrmpinirecvs + 1
  END IF
  nbrrank = GR2Rank(endglobrow+1)
  IF ( ISODD(endglobrow) .AND. (nbrrank .GE. 0) ) THEN
    IF ( .NOT. ALLOCATED( selement(endglobrow+1)%x ) ) THEN
        ALLOCATE( selement(endglobrow+1)%x(M) )
    END IF
    CALL MPI_Irecv( selement(endglobrow+1)%x, M, MPI_REAL8, nbrrank, 2, &
      SIESTA_COMM, nbrmpireq(nbrmpireqindx), nbrmpierr(nbrmpireqindx) )
    nbrmpireqindx = nbrmpireqindx + 1
    nbrmpinirecvs = nbrmpinirecvs + 1
  END IF
  nbrrank = GR2Rank(startglobrow-1)
  IF ( ISEVEN(startglobrow) .AND. (nbrrank .GE. 0) ) THEN
    CALL MPI_Isend( selement(startglobrow)%x, M, MPI_REAL8, nbrrank, 2, &
      SIESTA_COMM, nbrmpireq(nbrmpireqindx), nbrmpierr(nbrmpireqindx) )
    nbrmpireqindx = nbrmpireqindx + 1
    nbrmpinisends = nbrmpinisends + 1
  END IF
  nbrrank = GR2Rank(endglobrow+1)
  IF ( ISEVEN(endglobrow) .AND. (nbrrank .GE. 0) ) THEN
    CALL MPI_Isend( selement(endglobrow)%x, M, MPI_REAL8, nbrrank, 1, &
      SIESTA_COMM, nbrmpireq(nbrmpireqindx), nbrmpierr(nbrmpireqindx) )
    nbrmpireqindx = nbrmpireqindx + 1
    nbrmpinisends = nbrmpinisends + 1
  END IF
  nbrmpireqcnt = nbrmpireqindx - 1
  IF ( nbrmpireqcnt .GT. 0 ) THEN
    IF(KPDBG) WRITE(OFU,*) ' Wait ',nbrmpinirecvs,' ',nbrmpinisends; CALL FL(OFU)
    CALL MPI_Waitall( nbrmpireqcnt, nbrmpireq, nbrmpistatuses, mpiwaiterr )
  END IF

  !-----------------------------------------------------------------------------
  totrmserr = 0
  DO globrow = startglobrow, endglobrow, 1
    globrowoff = globrow - startglobrow + 1
    DO i = 1, M, 1
      term = 0
      IF ( globrow .GT. 1 ) THEN
        term = term + SUM( orig(globrowoff)%L(i,:) * selement(globrow-1)%x )
      END IF

      term = term + SUM( orig(globrowoff)%D(i,:) * selement(globrow)%x )

      IF (writesolution) THEN !< Write solution to output
        IF( i .EQ. 1 ) THEN
          DO k = 1, M
            IF(KPDBG) WRITE(OFU,*) 'X[',(globrow-1)*M+k,']=', selement(globrow)%x(k);CALL FL(OFU)
          END DO
        END IF
      END IF

      IF ( globrow .LT. N ) THEN
        term = term + SUM( orig(globrowoff)%U(i,:) * selement(globrow+1)%x )
      END IF

      totrmserr = totrmserr + (term - orig(globrowoff)%b(i,1))**2
    END DO
  END DO

  IF ( endglobrow .LT. startglobrow ) THEN
    totrmserr = 0
  ELSE
    totrmserr = SQRT( totrmserr / ((endglobrow-startglobrow+1) * M) )
  END IF
  IF(KPDBG) WRITE(OFU,'(A,E15.8E3)') 'TOTAL RMS ERROR = ', totrmserr; CALL FL(OFU)
  IF(KPDBG) WRITE(OFU,*) '------ Solution verified ------'; CALL FL(OFU)
END SUBROUTINE VerifySolution
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
SUBROUTINE CheckConditionNumber(nblkrow,bsize,anorm,rcond,info)
  
  INTEGER, INTENT(IN)   :: nblkrow, bsize
  INTEGER, INTENT(OUT)  :: info
  REAL(dp), INTENT(OUT) :: rcond, anorm
  
  
  INTEGER :: iblkrow, n1, i, j, offs, offsu, offk, kl, ku, ldab, icol,  &
             rowmin, rowmax
  REAL(dp), ALLOCATABLE :: Band(:,:), Acol(:), work(:)
  INTEGER, ALLOCATABLE, DIMENSION(:) :: iwork, ipiv


  n1=nblkrow*bsize
  IF (nranks.NE.1) THEN
    PRINT *, "Returning from CheckConditionNumber"
    info = -1
    RETURN
  END IF

  ku = 2*bsize-1
  kl = 2*bsize-1
  ldab = 2*kl+ku+1
  offk = kl+ku+1
  
  ALLOCATE(Acol(n1))
  ALLOCATE(Band(ldab, n1), ipiv(n1), iwork(n1), work(3*n1))
  
  offs = 0
  offsu = bsize
  anorm = 0
  Band = 0

!compute all cols (labelled by "j"), row-by-row Acol(1:n1)
  DO iblkrow = 1, nblkrow
    DO icol = 1, bsize
       j = offs+icol
       rowmin = offs-bsize+1; rowmax = offsu+bsize
       IF (iblkrow.GT.1)                                                &
       Acol(rowmin:offs)=orig(iblkrow-1)%U(:,icol)
       IF (iblkrow.LT.nblkrow)                                          &
       Acol(offsu+1:rowmax)=orig(iblkrow+1)%L(:,icol)
       Acol(offs+1:offsu)=orig(iblkrow)%D(:,icol)

       IF (iblkrow.EQ.1) THEN
          ROWMIN = offs+1
       ELSE IF (iblkrow.EQ.nblkrow) THEN
          ROWMAX = offsu
       END IF
       ANORM = MAX(ANORM, SUM(ABS(Acol(rowmin:rowmax))))
       
       DO i = MAX(1,j-ku,rowmin), MIN(n1,j+kl,rowmax)
          Band(offk+i-j,j) = Acol(i)
       END DO
       
    END DO
    offs = offsu
    offsu = offsu + bsize
  END DO


  DEALLOCATE (Acol)

  CALL DGBTRF( N1, N1, KL, KU, BAND, LDAB, IPIV, INFO )
  CALL DGBCON( '1', N1, KL, KU, BAND, LDAB, IPIV, ANORM, RCOND,          &
               WORK, IWORK, INFO )

  IF (info.eq.0 .and. rcond.ne.zero) THEN
     RCOND = one/RCOND
  ELSE
     RCOND = -one
  END IF

  DEALLOCATE(Band,work,iwork,ipiv)

END SUBROUTINE CheckConditionNumber


SUBROUTINE CheckSymmetry (asymIndx)

REAL(dp), INTENT (OUT) :: asymIndx                                       !On Processor 0 only
REAL(dp) :: partdiff, eij, eji, summ
REAL(dp), DIMENSION(2) :: temp
REAL(dp), ALLOCATABLE, DIMENSION(:,:) :: sendBuf, bottomLeftMatrix

INTEGER :: topNeighbor, bottomNeighbor

INTEGER :: mpistatus(MPI_STATUS_SIZE), mpi_err
INTEGER :: blkrow,blktype,i,j

ALLOCATE(bottomLeftMatrix(M,M), sendBuf(M,M), stat=i)
CALL ASSERT(i.EQ.0,'Allocation error in CheckSymmetry')
bottomLeftMatrix=0
sendBuf=0

IF(nranks.EQ.1) THEN
  topNeighbor=MPI_PROC_NULL
  bottomNeighbor=MPI_PROC_NULL
ELSE
  IF (rank.EQ.0) THEN
    topNeighbor=MPI_PROC_NULL
    bottomNeighbor=rank+1
  ELSE IF (rank.EQ.nranks-1) THEN
    topNeighbor=rank-1
    bottomNeighbor=MPI_PROC_NULL
    sendBuf=lelement(1,1)%L
  ELSE
    topNeighbor=rank-1
    bottomNeighbor=rank+1
    sendBuf=lelement(1,1)%L
  END IF
END IF

CALL MPI_Sendrecv(sendBuf,M*M,MPI_REAL8,topNeighbor,1,                  &
     bottomLeftMatrix,M*M,MPI_REAL8,bottomNeighbor,1,SIESTA_COMM,       &
     mpistatus,mpi_err)

DEALLOCATE(sendBuf)
IF (rank.EQ.nranks-1) bottomLeftMatrix=0

summ = 0
temp = 0

DO blkrow=1,endglobrow-startglobrow+1
  ! Compute contributions from the main diagonal blocks
  partdiff=0
  DO i=1,M
    eij=lelement(1,blkrow)%D(i,i)
    eji=lelement(1,blkrow)%D(i,i)
    temp(2)  = temp(2) + eij*eij
    DO j=i + 1,M
      eij=lelement(1,blkrow)%D(i,j)
      eji=lelement(1,blkrow)%D(j,i)
      partdiff = partdiff+(eij-eji)*(eij-eji)
      temp(2) = temp(2) + eij*eij + eji*eji
    END DO
  END DO
  temp(1) = temp(1) + 2*partdiff

  ! Compute contributions from off diagonal blocks
  IF (blkrow.NE.endglobrow-startglobrow+1) THEN
    DO i = 1, M
      DO j = 1, M
        eij=lelement(1,blkrow)%U(i,j)
        eji=lelement(1,blkrow+1)%L(j,i)
        temp(1) = temp(1) + (eij - eji)*(eij - eji)
        temp(2) = temp(2) + eij*eij + eji*eji
      END DO
    END DO
  ELSE IF (rank.NE.nranks-1) THEN
!Local bottom block row
    DO i = 1, M
      DO j = 1, M
        eij=lelement(1,blkrow)%U(i,j)
        eji=bottomLeftMatrix(j,i)
        temp(1) = temp(1) + (eij - eji)*(eij - eji)
        temp(2) = temp(2) + eij*eij + eji*eji
      END DO
    END DO
  END IF
END DO

CALL FLUSH(360+rank)

DEALLOCATE(bottomLeftMatrix)

IF (rank .eq. 0) THEN
   CALL MPI_Reduce(MPI_IN_PLACE,temp,2,MPI_REAL8,MPI_SUM,0,SIESTA_COMM,mpi_err)
ELSE
   CALL MPI_Reduce(temp,temp,2,MPI_REAL8,MPI_SUM,0,SIESTA_COMM,mpi_err)
END IF
IF (rank.EQ.0) THEN 
  asymIndx = SQRT(temp(1)/temp(2))
  !print *,'New asymIndx:',asymIndx
END IF

END SUBROUTINE CheckSymmetry
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
SUBROUTINE StoreDiagonal(blkrow, colnum, buf)

INTEGER, INTENT(IN) :: blkrow, colnum
REAL(dp), DIMENSION(1:M), INTENT(IN) :: buf
INTEGER :: indx

IF (blkrow .LT. startglobrow) THEN
   IF(KPDBG) WRITE(OFU,*) 'Error in indexing global block row index'; CALL FL(OFU)
END IF

indx = (blkrow - startglobrow)*M + colnum
OrigDiagElement(indx) = buf(colnum)

END SUBROUTINE StoreDiagonal
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
SUBROUTINE WriteBlocks(flag)

INTEGER, INTENT(IN) :: flag
!LOGICAL, INTENT(IN) :: flag
INTEGER :: row, col, localblkrow, globalrow

CHARACTER(100) :: fname
CHARACTER(50) :: ciam, cnprocs
INTEGER :: WOFU, istat

WRITE(ciam,*) rank; WRITE(cnprocs,*) nranks
ciam=ADJUSTL(ciam); cnprocs=ADJUSTL(cnprocs)
WOFU = 4*nranks+rank

fname='block-'//TRIM(ciam)//'-P-'//TRIM(cnprocs)//'.txt'
OPEN(UNIT=WOFU, FILE=fname, STATUS="REPLACE", ACTION="WRITE",&
&FORM="FORMATTED",POSITION="APPEND", IOSTAT=istat)
!OPEN(UNIT=WOFU, FILE=fname, STATUS="REPLACE", ACTION="WRITE",&
!&POSITION="REWIND", IOSTAT=istat)

IF(.TRUE.) THEN
DO row=1, (endglobrow-startglobrow+1)*M
  WRITE(WOFU,*) 'SavedDiag:',row,origDiagElement(row), flag
  CALL FL(WOFU)
END DO
END IF

IF(.TRUE.) THEN
DO globalrow=startglobrow,endglobrow
  localblkrow=globalrow-startglobrow+1
  
  WRITE(WOFU,*) 
  CALL FL(WOFU)
  DO row=1,M
    DO col=1,M
      WRITE(WOFU,*) globalrow, 'L', row, col, orig(localblkrow)%L(row,col), flag
      CALL FL(WOFU)
    END DO
  END DO

  WRITE(WOFU,*) 
  CALL FL(WOFU)
  DO row=1,M
    DO col=1,M
      WRITE(WOFU,*) globalrow, 'D', row, col, orig(localblkrow)%D(row,col), flag
      CALL FL(WOFU)
    END DO
  END DO

  WRITE(WOFU,*) 
  CALL FL(WOFU)
  DO row=1,M
    DO col=1,M
      WRITE(WOFU,*) globalrow, 'U', row, col, orig(localblkrow)%U(row,col), flag
      CALL FL(WOFU)
    END DO
  END DO
  WRITE(WOFU,*)'--------------------------------------------------------' 
  CALL FL(WOFU)

END DO
END IF

END SUBROUTINE WriteBlocks
!-------------------------------------------------------------------------------


!-------------------------------------------------------------------------------
SUBROUTINE DisplayBlocks(flag)

INTEGER, INTENT(IN) :: flag
INTEGER :: row, localblkrow, globalrow

IF(KPDBG) WRITE(OFU,*) '------- Here ------- ', 1; CALL FL(OFU)

IF(KPDBG) WRITE(OFU,*) ' '; CALL FL(OFU)
DO row=1, (endglobrow-startglobrow+1)*M
  IF(KPDBG) WRITE(OFU,*) 'OrigDiagElement:',OrigDiagElement(row),'-',flag
  CALL FL(OFU)
END DO
IF(KPDBG) WRITE(OFU,*) ' '; CALL FL(OFU)

IF(KPDBG) WRITE(OFU,*) '------- Here ------- ', 2; CALL FL(OFU)


DO globalrow=startglobrow,endglobrow
  localblkrow=globalrow-startglobrow+1
  
  IF(KPDBG) WRITE(OFU,*) ' '; CALL FL(OFU)

  DO row=1,M
    IF(KPDBG) WRITE(OFU,*) 'orig-L:',orig(localblkrow)%L(row,1:M),'-',flag
    CALL FL(OFU)
  END DO

  IF(KPDBG) WRITE(OFU,*) ' '; CALL FL(OFU)
  DO row=1,M
    IF(KPDBG) WRITE(OFU,*) 'orig-D:',orig(localblkrow)%D(row,1:M),'-',flag
    CALL FL(OFU)
  END DO

  IF(KPDBG) WRITE(OFU,*) ' '; CALL FL(OFU)
  DO row=1,M
    IF(KPDBG) WRITE(OFU,*) 'orig-U:',orig(localblkrow)%U(row,1:M),'-',flag
    CALL FL(OFU)
  END DO
  IF(KPDBG) WRITE(OFU,*) ' '; CALL FL(OFU)
END DO

IF(KPDBG) WRITE(OFU,*) '------- Here ------- ', 3; CALL FL(OFU)


END SUBROUTINE DisplayBlocks
!-------------------------------------------------------------------------------


!-------------------------------------------------------------------------------
SUBROUTINE ParMatVec (invec,outvec,outveclength)

INTEGER,INTENT(IN) :: outveclength 
REAL(dp), INTENT(IN) :: invec(1:N*M)
REAL(dp), INTENT(OUT) :: outvec(1:outveclength)

INTEGER :: col, row, offset
INTEGER :: globalrow, localblkrow, cnt, dcnt
REAL(dp) :: melement, velement, total

outvec=0; cnt=0; dcnt=0
DO globalrow=startglobrow,endglobrow

  localblkrow=globalrow-startglobrow+1
  
  IF (globalrow.EQ.1) THEN
    offset=0
  ELSE
    offset=(globalrow-2)*M
  END IF
  
  IF (globalrow.EQ.1) THEN
    DO row=1,M
      total=0
      DO col=1,2*M
        velement=invec(offset+col)
        IF (col.LE.M) THEN
          IF(row.EQ.col) THEN
            dcnt=dcnt+1
            melement=OrigDiagElement(dcnt)
          ELSE
            melement=orig(localblkrow)%D(row,col)
          END IF
        ELSE 
          melement=orig(localblkrow)%U(row,col-M)
        END IF
        total=total+melement*velement
      END DO
      cnt=cnt+1
      outvec(cnt)=total
    END DO
  ELSE IF (globalrow.EQ.N) THEN
    DO row=1,M
      total=0
      DO col=1,2*M
        velement=invec(offset+col)
        IF (col.LE.M) THEN 
          melement=orig(localblkrow)%L(row,col)
        ELSE
          IF (row.EQ.col-M) THEN
            dcnt=dcnt+1
            melement=OrigDiagElement(dcnt)
          ELSE
            melement=orig(localblkrow)%D(row,col-M)
          ENDIF
        END IF
        total=total+melement*velement
      END DO
      cnt=cnt+1
      outvec(cnt)=total
    END DO
  ELSE
    DO row=1,M
      total=0
      DO col=1,3*M
        velement=invec(offset+col)
        IF (col.LE.M) THEN
          melement=orig(localblkrow)%L(row,col)
        ELSE IF (col.GT.M.AND.col.LE.2*M) THEN 
          IF (row.EQ.col-M) THEN
            dcnt=dcnt+1
            melement=OrigDiagElement(dcnt)
          ELSE
            melement=orig(localblkrow)%D(row,col-M)
          END IF
        ELSE
          melement=orig(localblkrow)%U(row,col-2*M)
        END IF
        total=total+melement*velement
        IF (globalrow.EQ.2.AND.row.EQ.1) THEN
        END IF
      END DO
      cnt=cnt+1
      outvec(cnt)=total
    END DO
  END IF
END DO

END SUBROUTINE ParMatVec
!------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
SUBROUTINE GetColSum(cs)
!-------------------------------------------------------------------------------
REAL (dp), DIMENSION(M, N), INTENT(OUT) :: cs
!-------------------------------------------------------------------------------
INTEGER                                 :: js
REAL (dp)                               :: minScale
REAL (dp)                               :: eps
!-------------------------------------------------------------------------------
IF (lequi1) THEN
   minScale = 1.0E-5_dp
ELSE
   minScale = 1.0E-8_dp
END IF
CALL InitScaleFactors
cs = 0.0
DO js = startglobrow, endglobrow
   CALL GetScaleFactors(js,cs(:,js))
   eps = minScale*MAXVAL(cs(:,js))
   CALL ASSERT(eps .gt. zero, ' cs == 0 in GetColSum')
   cs(:,js) = one/SQRT(MAX(cs(:,js), eps))
END DO

CALL FinalizeScaleFactors

END SUBROUTINE GetColSum
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
SUBROUTINE InitScaleFactors

INTEGER :: top, bot, MPI_ERR, ic, j, k
INTEGER :: MPI_STAT(MPI_STATUS_SIZE)
REAL(dp), DIMENSION(:), ALLOCATABLE :: tmp1, tmp2
INTEGER :: numBlockRows

numBlockRows = endglobrow-startglobrow+1

top=rank-1; IF(rank.EQ.0) top=MPI_PROC_NULL
bot=rank+1; IF(rank.EQ.nranks-1) bot=MPI_PROC_NULL

ALLOCATE (TopScaleFac(M), BotScaleFac(M), tmp1(M), tmp2(M), stat=ic)
CALL ASSERT(ic.eq.0,'STAT != 0 IN InitScaleFactors')

DO ic=1, M
   IF (lequi1) THEN
      tmp1(ic) = SUM(ABS(orig(1)%L(:,ic)))
      tmp2(ic) = SUM(ABS(orig(numBlockRows)%U(:,ic)))
   ELSE
      tmp1(ic) = MAXVAL(ABS(orig(1)%L(:,ic)))
      tmp2(ic) = MAXVAL(ABS(orig(numBlockRows)%U(:,ic)))
   ENDIF
END DO

CALL MPI_Sendrecv(tmp1,M,MPI_REAL8,top,1,                               &
  BotScaleFac,M,MPI_REAL8,bot,1,SIESTA_COMM,MPI_STAT,MPI_ERR)

CALL MPI_Sendrecv(tmp2,M,MPI_REAL8,bot,1,                               &
  TopScaleFac,M,MPI_REAL8,top,1,SIESTA_COMM,MPI_STAT,MPI_ERR)

DEALLOCATE(tmp1, tmp2)

END SUBROUTINE InitScaleFactors
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
SUBROUTINE GetScaleFactors(js,scalevector)
!-------------------------------------------------------------------------------
INTEGER, INTENT(IN)   :: js
REAL(dp), INTENT(OUT) :: scalevector(M)
!-------------------------------------------------------------------------------
REAL(dp), ALLOCATABLE :: coltmp(:)
REAL(dp) :: eps
INTEGER  :: ic, ir
!-------------------------------------------------------------------------------

ALLOCATE (coltmp(M))

ir = js - startglobrow + 1
!COLUMN SUMS: 1-norm; MAXVAL: infinity norm
DO ic=1, M

   coltmp = ABS(orig(ir)%D(:,ic))
   IF (lequi1) THEN
      scalevector(ic) = SUM(coltmp)
   ELSE
      scalevector(ic) = MAXVAL(coltmp)
   END IF

   IF (js.GT.startglobrow) THEN
      coltmp = ABS(orig(ir-1)%U(:,ic))
      IF (lequi1) THEN
         scalevector(ic) = scalevector(ic) + SUM(coltmp)
         IF (js.EQ.endglobrow .AND. js.NE.N) THEN
            scalevector(ic) = scalevector(ic) + BotScaleFac(ic)
         END IF
      ELSE
         scalevector(ic) = MAX(scalevector(ic),MAXVAL(coltmp))
         IF (js.EQ.endglobrow .AND. js.NE.N) THEN
            scalevector(ic) = MAX(scalevector(ic),BotScaleFac(ic))
         END IF
      END IF
   END IF

   IF (js.LT.endglobrow) THEN
      coltmp = ABS(orig(ir+1)%L(:,ic))
      IF (lequi1) THEN
         IF (js.EQ.startglobrow .AND. js.NE.1) THEN
            scalevector(ic) = scalevector(ic) + TopScaleFac(ic)
         END IF
         scalevector(ic) = scalevector(ic) + SUM(coltmp)
      ELSE
         IF (js.EQ.startglobrow .AND. js.NE.1) THEN
            scalevector(ic) = MAX(scalevector(ic),TopScaleFac(ic))
         END IF
         scalevector(ic) = MAX(scalevector(ic),MAXVAL(coltmp))
      END IF
   END IF
END DO

DEALLOCATE (coltmp)
CALL ASSERT(ALL(scalevector.GT.zero),'SCALEVECTOR <= 0 IN GETSCALEFACTORS')

END SUBROUTINE GetScaleFactors
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
SUBROUTINE FinalizeScaleFactors
  DEALLOCATE (TopScaleFac,BotScaleFac)
END SUBROUTINE FinalizeScaleFactors
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
SUBROUTINE ApplyParallelScaling(lmp, colscale, lverbose)
!-------------------------------------------------------------------------------
REAL(dp), INTENT(IN)    :: lmp
REAL(dp), INTENT(INOUT) :: colscale(M,N)
LOGICAL, INTENT(IN)     :: lverbose
!-------------------------------------------------------------------------------
INTEGER               :: icount, istat
REAL(dp)              :: ton, toff, colcnd, colcnd1, temp(2)
REAL(dp), ALLOCATABLE :: colpad(:,:)
!-------------------------------------------------------------------------------
CALL second0(ton)

icount = 0
startrow1 = MAX(1, startglobrow-1); endrow1 = MIN(N, endglobrow+1)
CALL ASSERT(ALL(colscale(:,startrow1:endrow1).EQ.1),                    &
            'COLSCALE !=1 INITIALLY IN APPLYPARALLELSCALING')

ALLOCATE (colpad(M,N), stat=istat)
CALL ASSERT(istat.eq.0,' Allocation error in ApplyParallelScaling')

CALL GetColSum(colpad)
CALL PadSides(colpad, M, 1, 1)
CALL VectorCopyPar(colpad, colscale)
CALL ParallelScaling(colpad)

icount = icount + 1
!MAXIMUM COLUMN VARIATIONS (like an inverse "Condition number")
temp(1) = MAXVAL( colpad(:,startglobrow:endglobrow))
temp(2) = MAXVAL(-colpad(:,startglobrow:endglobrow))
CALL MPI_ALLREDUCE(MPI_IN_PLACE, temp, 2, MPI_REAL8, MPI_MAX,  SIESTA_COMM,    &
                   MPI_ERR)
colcnd1 = temp(1)
colcnd  = -temp(2)
                
IF (colcnd1 .NE. ZERO) THEN
   colcnd = ABS(colcnd/colcnd1)
ELSE 
   colcnd = 1
ENDIF
                
!IF (icount.LE.6 .AND. colcnd.LT.0.01_dp) GO TO 100                       !dgeequ logic

DEALLOCATE(colpad)

!Add levenberg to diagonals
CALL FindMinMax_TRI(lmp)

CALL second0(toff)
IF (rank.EQ.0 .and. lverbose) PRINT '(1x,3(a,1p,e10.3))',                      &
              'COLUMN-SCALING TIME: ', toff-ton, ' COLCND: ', colcnd,   &
              ' DMINL: ', DMIN_TRI*lmp

END SUBROUTINE ApplyParallelScaling
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
SUBROUTINE ParallelScaling(colsum)
!-------------------------------------------------------------------------------
REAL(dp), DIMENSION(M,N), INTENT(IN) :: colsum
!-------------------------------------------------------------------------------
INTEGER  :: js, icol, ir, dcnt
!-------------------------------------------------------------------------------

L_COLSCALE = .TRUE.

dcnt = 0
DO js = startglobrow, endglobrow
  ir=js-startglobrow+1
  DO icol = 1, M
    lelement(1,ir)%D(:,icol) = colsum(:,js)*orig(ir)%D(:,icol)*colsum(icol,js)
    orig(ir)%D(:,icol) = lelement(1,ir)%D(:,icol)
    dcnt=dcnt+1
    OrigDiagElement(dcnt) = orig(ir)%D(icol,icol)
    IF (js .GT. 1) THEN
       lelement(1,ir)%L(:,icol) = colsum(:,js)*orig(ir)%L(:,icol)*colsum(icol,js-1)
       orig(ir)%L(:,icol) = lelement(1,ir)%L(:,icol)
    END IF
    IF (js .LT. N) THEN
       lelement(1,ir)%U(:,icol) = colsum(:,js)*orig(ir)%U(:,icol)*colsum(icol,js+1)
       orig(ir)%U(:,icol) = lelement(1,ir)%U(:,icol)
    END IF
  END DO
END DO

END SUBROUTINE ParallelScaling
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
SUBROUTINE FindMinMax_TRI(lmp)
!SPH 091616: Find min diagonal element and max eigenvalue (Gerschgorin Theorem)
REAL(dp), PARAMETER  :: minThreshold=1.E-6_dp
REAL(dp), INTENT(IN) :: lmp
INTEGER  :: js, icol, ir, dcnt, imax_diag(1)
REAL(dp) :: eps, temp(2), minDiag, sign_diag = -1

dcnt=0
maxeigen_tri = 0
Dmin_TRI = HUGE(Dmin_TRI)
DO js = startglobrow, endglobrow
  ir=js-startglobrow+1

!Minimum allowable diagonal element after scaling
  minDiag = minThreshold*MAXVAL(ABS(OrigDiagElement(dcnt+1:dcnt+M)))
  IF (minDiag .EQ. ZERO) THEN
     minDiag = minThreshold
  ELSE
     imax_diag = MAXLOC(ABS(OrigDiagElement(dcnt+1:dcnt+M)))
     sign_diag = SIGN(one, OrigDiagElement(imax_diag(1)))
  END IF

  DO icol = 1, M
    dcnt = dcnt+1
    eps = OrigDiagElement(dcnt)

    IF (ABS(eps) .LE. minDiag) THEN                                      !Can happen for l_colscale=FALSE
       eps = sign_diag*minDiag
       OrigDiagElement(dcnt) = eps
       orig(ir)%D(icol,icol) = eps
       lelement(1,ir)%D(icol,icol) = eps
    ELSE
       IF (js .LT. N) Dmin_TRI = MIN(MAX(1.E-12_dp,ABS(eps)), Dmin_TRI)
    END IF

    maxeigen_tri = MAX(maxeigen_tri, SUM(ABS(orig(ir)%L(icol,:)) +      &
                   ABS(orig(ir)%D(icol,:)) + ABS(orig(ir)%U(icol,:))))
  END DO
END DO

temp(1)=-Dmin_TRI; temp(2)=maxeigen_tri
CALL MPI_ALLREDUCE(MPI_IN_PLACE,temp,2,MPI_REAL8, MPI_MAX, SIESTA_COMM,MPI_ERR)
Dmin_TRI=SIGN(-temp(1), eps)
maxeigen_tri=temp(2)


!ADD LEVENBERG PARAMETER LM*DMIN TO DIAGONAL [OR D*(1+LM) IF NO COL-SCALING]
CALL ResetDiagonal(lmp, .FALSE.)

END SUBROUTINE FindMinMax_TRI
!-------------------------------------------------------------------------------


!-------------------------------------------------------------------------------
SUBROUTINE ResetDiagonal(lmp, bReset)
REAL(dp), INTENT(IN) :: lmp
LOGICAL,  INTENT(IN) :: bReset
REAL(dp)             :: DminL, eps, eps0
INTEGER              :: dcnt, js, ir, icol

!ADD CONSTANT (SCALED BY DMIN) LEVENBERG PARAMETER TO DIAGONAL [NOT D*(1+LM)]
eps0 = SIGN(one, OrigDiagElement(1))
Dmin_TRI = SIGN(Dmin_TRI,eps0)
DminL = Dmin_TRI*ABS(lmp)

dcnt = 0

DO js = startglobrow, endglobrow
  ir=js-startglobrow+1
  DO icol = 1, M
    dcnt=dcnt+1
    eps = OrigDiagElement(dcnt)
    CALL ASSERT(SIGN(one,eps).EQ.eps0,'WRONG SIGN IN RESET DIAGONAL')
    IF (L_COLSCALE) THEN
       lelement(1,ir)%D(icol,icol) = eps + DminL
    ELSE
       lelement(1,ir)%D(icol,icol) = eps*(1 + ABS(lmp))
    END IF
    orig(ir)%D(icol,icol) = lelement(1,ir)%D(icol,icol)
  END DO

  IF (bReset) THEN
     lelement(1,ir)%L = orig(ir)%L  
     lelement(1,ir)%D = orig(ir)%D 
     lelement(1,ir)%U = orig(ir)%U
     lelement(1,ir)%pivot = 0
  END IF
END DO

END SUBROUTINE ResetDiagonal
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
SUBROUTINE RefactorHessian(lmp)
!SPH 08-31-16
!Routine to reset Hessian block prior to refactoring it with a new levenberg parameter (lmp)
!
REAL(dp), INTENT(IN) :: lmp

CALL ResetDiagonal(lmp, .TRUE.)
CALL ForwardSolve

END SUBROUTINE RefactorHessian
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
SUBROUTINE ComputeTranspose



END SUBROUTINE ComputeTranspose
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
SUBROUTINE VectorCopyPar (colsum, colscale)
!-----------------------------------------------
REAL(dp), INTENT(IN)    :: colsum(M,N)
REAL(dp), INTENT(INOUT) :: colscale(M,N)
!-----------------------------------------------
INTEGER               :: js
!-----------------------------------------------
DO js = startrow1, endrow1
   colscale(:,js) = colsum(:,js)*colscale(:,js)
END DO

END SUBROUTINE VectorCopyPar
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
SUBROUTINE PadSides(arrin, blksize, top, bot)
!-----------------------------------------------
INTEGER, INTENT(IN)     :: top, bot, blksize
REAL(dp), INTENT(INOUT) :: arrin(blksize, N)
!-----------------------------------------------
INTEGER                 :: MPI_STAT(MPI_STATUS_SIZE)
INTEGER                 :: M1, left, right, tag1=1, tag2=2,             &
                           tr, t1l, tl, t1r
!-----------------------------------------------
left=rank-1;  IF(rank.EQ.0) left=MPI_PROC_NULL
right=rank+1; IF(rank.EQ.nranks-1) right=MPI_PROC_NULL

tr = endglobrow;              tl = startglobrow
t1r= MIN(endglobrow+1, N);    t1l= MAX(1, startglobrow-1)

!Send top / Rec bottom elements
M1 = blksize*top
CALL MPI_Sendrecv(arrin(:,tr-(top-1):tr), M1, MPI_REAL8,                &
                  right, tag1, arrin(:,t1l-(top-1):t1l), M1, MPI_REAL8, &
                  left,  tag1, SIESTA_COMM, MPI_STAT, MPI_ERR)

!Send bottom / Rec top boundary elements
M1 = blksize*bot
CALL MPI_Sendrecv(arrin(:,tl:tl+(bot-1)), M1, MPI_REAL8,                &
                  left,  tag2, arrin(:,t1r:t1r+(bot-1)), M1, MPI_REAL8, &
                  right, tag2, SIESTA_COMM, MPI_STAT, MPI_ERR)
END SUBROUTINE PadSides

#endif
END MODULE blocktridiagonalsolver_s
!-------------------------------------------------------------------------------
