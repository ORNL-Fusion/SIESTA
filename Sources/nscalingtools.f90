!----------------------------------------------------------------------------
! Module to scale SIESTA with the number of surfaces (ns). [S. K. Seal, ORNL]
! Date: April 26, 2013
!----------------------------------------------------------------------------

!------------------------------------------------
MODULE nscalingtools
USE stel_kinds
USE v3_utilities, ONLY: assert
USE descriptor_mod, ONLY: SIESTA_COMM
#if defined(MPI_OPT)
USE blocktridiagonalsolver_s, ONLY: &
& SetMatrixRowColL, SetMatrixRowColD, SetMatrixRowColU, &
& GetMatrixRowColL, GetMatrixRowColD, GetMatrixRowColU, &
& GetMatrixRHS, StoreDiagonal
USE mpi_inc
#endif
IMPLICIT NONE
!------------------------------------------------
! User-controlled logicals
!------------------------------------------------
LOGICAL :: SKSDBG=.FALSE.
LOGICAL :: KPDBG=.FALSE.
LOGICAL :: PARSOLVER=.FALSE.
LOGICAL :: SYMMETRYCHECK=.TRUE.
LOGICAL :: PARGMRES=.FALSE.
LOGICAL :: BUFFERED=.FALSE.
LOGICAL :: PARFUNCTISL=.FALSE.
LOGICAL :: OUTPUT_TIMINGS=.FALSE.
!------------------------------------------------

LOGICAL :: INITREMAPDONE=.FALSE.
LOGICAL :: MAPPINGDONE=.FALSE.
LOGICAL :: FIRSTPARTITION_GE_2=.FALSE.
LOGICAL :: LASTPARTITION_GE_2=.FALSE.
CHARACTER(50) :: crank, cactiveranks
CHARACTER*100 :: envvar
CHARACTER*100 :: envval
CHARACTER, ALLOCATABLE :: remapBuffer(:)

INTEGER :: TOFU
INTEGER :: OpenStatus, nargs
INTEGER, PRIVATE :: N, M, K, KSQR

INTEGER, PARAMETER :: SAVEDIAG = 4
INTEGER, PARAMETER :: LOWER    = 3
INTEGER, PARAMETER :: DIAG     = 2
INTEGER, PARAMETER :: UPPER    = 1

!This is to stitch to blocktri variables of
!same name. Just USE nscalingtools.
INTEGER :: startglobrow, endglobrow, numBlocks
INTEGER, ALLOCATABLE, DIMENSION (:) :: bcyclicStartBlockProc, bcyclicEndBlockProc
INTEGER, ALLOCATABLE, DIMENSION (:), PRIVATE :: tagBlock, tagColInBlock
INTEGER, ALLOCATABLE, DIMENSION (:), TARGET :: rcounts, disp             !SPH added TARGET
INTEGER, ALLOCATABLE, DIMENSION (:) :: rcountsNS, dispNS                 !SPH 092116
INTEGER, PRIVATE :: a, b, c, flag, cnt 
INTEGER, PRIVATE :: icol, itype, block_row, irow, procID, ierror, localbrow
INTEGER :: activeranks, rank, nranks, leftproc, rightproc, totRecvs, nrecd
INTEGER :: tag, MPI_ERR, SYSMAXTAG, APPMAXTAG
INTEGER :: PACKSIZE

! HESSIAN MN-SPACE GATHER-related variables
INTEGER :: linblksize
INTEGER, ALLOCATABLE, DIMENSION(:) :: mnspcounts, mnspdisps
LOGICAL, PRIVATE :: SetUpTOMNSPAllGatherDONE=.FALSE.

#if defined(MPI_OPT)
INTEGER :: siesta_MPI_Status(MPI_STATUS_SIZE)
#endif 
!------------------------------------------------

CONTAINS

!--------------------------------------------------------------------------
SUBROUTINE  MyEnvVariables

!--------------------------------------------------------------------------
! Default values of environment values that users are allowed to change
!--------------------------------------------------------------------------
PARSOLVER=.TRUE.
PARFUNCTISL=.TRUE.
!--------------------------------------------------------------------------

!--------------------------------------------------------------------------
! Default values of environment values that users are NOT allowed to change
!--------------------------------------------------------------------------
OUTPUT_TIMINGS=.FALSE. !SKS: USERS MUST NOT CHANGE THIS VALUE: 04/26/2013
PARGMRES=.TRUE.        !SKS: USERS MUST NOT CHANGE THIS VALUE: 04/26/2013
BUFFERED=.FALSE.       !SKS: USERS MUST NOT CHANGE THIS VALUE: 04/26/2013
SYMMETRYCHECK=.TRUE.   !SKS: USERS MUST NOT CHANGE THIS VALUE: 04/26/2013
SKSDBG=.FALSE.         !SKS: USERS MUST NOT CHANGE THIS VALUE: 04/26/2013
KPDBG=.FALSE.          !SKS: USERS MUST NOT CHANGE THIS VALUE: 04/26/2013
!--------------------------------------------------------------------------

!--------------------------------------------------------------------------
! Read non-default values, if any
!--------------------------------------------------------------------------
envvar='PARSOLVER'
CALL GETENV(envvar,envval)
IF (envval.EQ.'ScaLAPACK' .OR. envval.EQ.'FALSE') THEN
  PARSOLVER=.FALSE.
  PARFUNCTISL=.FALSE.
  PARGMRES=.FALSE.
  OUTPUT_TIMINGS=.FALSE.
END IF
END SUBROUTINE MyEnvVariables
!--------------------------------------------------------------------------

!------------------------------------------------
SUBROUTINE SetOutputFile (iam, nprocs)

INTEGER :: iam, nprocs, istat
CHARACTER(100) :: fname, cfname
CHARACTER(50) :: ciam, cnprocs

WRITE(ciam,*) iam; WRITE(cnprocs,*) nprocs
ciam=ADJUSTL(ciam); cnprocs=ADJUSTL(cnprocs)
TOFU = 2*nprocs+iam

fname='sks-'//TRIM(ciam)//'-P-'//TRIM(cnprocs)//'.txt'
OPEN(UNIT=TOFU, FILE=fname, STATUS="REPLACE", ACTION="WRITE",&
&FORM="FORMATTED",POSITION="APPEND", IOSTAT=istat)

WRITE(TOFU,*)'SKSDBG:', SKSDBG; CALL FLUSH(TOFU)
WRITE(TOFU,*)'KPDBG:', KPDBG; CALL FLUSH(TOFU)

WRITE(TOFU,*)'SYMMETRYCHECK:', SYMMETRYCHECK; CALL FLUSH(TOFU)
WRITE(TOFU,*)'BUFFERED:', BUFFERED; CALL FLUSH(TOFU)
WRITE(TOFU,*)'PARGMRES:', PARGMRES; CALL FLUSH(TOFU)
WRITE(TOFU,*)'PARSOLVER:', PARSOLVER; CALL FLUSH(TOFU)
WRITE(TOFU,*)

END SUBROUTINE SetOutputFile
!------------------------------------------------

!------------------------------------------------
! This routine initializes all the necessary 
! variables used by all the routines in this 
! module
!------------------------------------------------
SUBROUTINE initRemap (mpol_in,ntor_in,blkrownum,numprocs,myrank)
USE shared_data, ONLY: ndims

INTEGER, INTENT(IN) :: mpol_in,ntor_in,blkrownum,numprocs,myrank 
INTEGER :: length
INTEGER :: i, FLAG, tmpint
INTEGER :: nBuffSize ! MRC Needed to compile with Gfortran 10.
!INTEGER*8 :: nBuffSize
!INTEGER, PARAMETER :: localprec = SELECTED_INT_KIND(24)
!INTEGER (localprec) :: nBuffSize

IF(INITREMAPDONE .OR. .NOT.PARSOLVER) RETURN

linblksize=(mpol_in+1)*(2*ntor_in+1)
M=ndims*linblksize
N=blkrownum

CALL ASSERT(M.NE.0, N.NE.0,'M, N = 0 IN initRemap')

nranks=numprocs
rank=myrank
nrecd=0;

IF(rank.GT.0) THEN 
  leftproc=rank-1
ELSE
  leftproc=MPI_PROC_NULL
END IF

IF(rank.LT.nranks-1) THEN
  rightproc=rank+1
ELSE
  rightproc=MPI_PROC_NULL
END IF

! To account for activeranks=P>N
activeranks=nranks
IF (activeranks.GT.N) activeranks=N

length=0; PACKSIZE=0
CALL MPI_Pack_size(2,MPI_INTEGER,SIESTA_COMM,length,MPI_ERR)
PACKSIZE=PACKSIZE+length
CALL MPI_Pack_size(M,MPI_REAL8,SIESTA_COMM,length,MPI_ERR)
PACKSIZE=PACKSIZE+length

!Each PACKSIZE packet has MPI_BSEND_OVERHEAD
IF (BUFFERED) THEN
  tmpint=PACKSIZE+MPI_BSEND_OVERHEAD
  nBuffsize=tmpint
  nBuffsize=nBuffsize*4
  nBuffsize=nBuffsize*M
  nBuffsize=nBuffsize*N
  nBuffsize=nBuffsize/MAX(nranks,1)
  IF (nBuffsize.LE.0) THEN
    CALL MPI_Barrier(SIESTA_COMM,MPI_ERR)
    IF (rank.EQ.0) THEN
      WRITE(7000+rank,*) '>>> N                  : ', N; CALL FLUSH(7000+rank)
      WRITE(7000+rank,*) '>>> M                  : ', M; CALL FLUSH(7000+rank)
      WRITE(7000+rank,*) '>>> PACKSIZE           : ', PACKSIZE; CALL FLUSH(7000+rank)
      WRITE(7000+rank,*) '>>> MPI_BSEND_OVERHEAD : ', MPI_BSEND_OVERHEAD; CALL FLUSH(7000+rank)
      WRITE(7000+rank,*) '>>> nBuffsize          : ', nBuffsize; CALL FLUSH(7000+rank)
    END IF
    CALL ASSERT(nBuffsize.GT.0,'initRemap nBuffsize')
  END IF

  ALLOCATE(remapBuffer(nBuffsize), stat=i)
  CALL ASSERT(i.eq.0,'MPI remapBuffer allocation failed in initRemap')
  CALL MPI_BUFFER_ATTACH(remapBuffer, nBuffsize, MPI_ERR)
END IF

IF (.NOT.ALLOCATED(bcyclicStartBlockProc)) ALLOCATE (bcyclicStartBlockProc(nranks))
IF (.NOT.ALLOCATED(bcyclicEndBlockProc)) ALLOCATE (bcyclicEndBlockProc(nranks))
INITREMAPDONE=.TRUE.

CALL bcyclicMapping
CALL SetUpTOMNSPAllGather
CALL computeAllGatherParameters

END SUBROUTINE initRemap
!------------------------------------------------

!------------------------------------------------
! This routine deallocates arrays, etc.
!------------------------------------------------
SUBROUTINE finalizeRemap

IF(ALLOCATED(bcyclicStartBlockProc)) DEALLOCATE(bcyclicStartBlockProc)
IF(ALLOCATED(bcyclicEndBlockProc)) DEALLOCATE(bcyclicEndBlockProc)
IF(ALLOCATED(remapBuffer)) DEALLOCATE(remapBuffer)

END SUBROUTINE finalizeRemap
!------------------------------------------------

!------------------------------------------------
! This computes the mapping of the block rows to 
! the processors. Computes the processors domain
! boundaries in terms of global block indices.
!------------------------------------------------
SUBROUTINE bcyclicMapping

!INTEGER, INTENT(OUT) :: startblock, endblock
INTEGER :: startblock, endblock
INTEGER :: lload, sload, myload
INTEGER :: numL, numS
INTEGER :: r, c
LOGICAL :: lnum, lnumg=.TRUE.

CALL ASSERT(INITREMAPDONE,                                              &
'Calling bcyclicMapping routine without calling initRemap')

IF (MAPPINGDONE) RETURN

lload=CEILING(REAL(N)/activeranks)
sload=FLOOR(REAL(N)/activeranks)

IF (lload.EQ.sload) THEN
  myload=lload
ELSE
  IF (rank.LT.MOD(N,activeranks)) THEN
    myload=lload
  ELSE
    myload=sload
  END IF
END IF

IF (sload.EQ.lload) THEN
  numS=0
  numL=rank
ELSE
  IF (myload.EQ.lload) THEN
    numL=rank
    numS=0
  ELSE
    numL=MOD(N,activeranks)
    numS=rank-numL
  END IF
END IF

IF (rank.LT.activeranks) THEN !active ranks
  startblock=numL*lload+numS*sload
  endblock=startblock+myload-1
ELSE                          !idle ranks
  startblock=-2
  endblock=-3
END IF

! Fortranized indices
startblock=startblock+1
endblock=endblock+1

CALL MPI_Allgather(startblock,1,MPI_INTEGER,bcyclicStartBlockProc,&
&1,MPI_INTEGER,SIESTA_COMM,MPI_ERR)

CALL MPI_Allgather(endblock,1,MPI_INTEGER,bcyclicEndBlockProc,&
&1,MPI_INTEGER,SIESTA_COMM,MPI_ERR)

startglobrow=bcyclicStartBlockProc(rank+1)
endglobrow=bcyclicEndBlockProc(rank+1)

FIRSTPARTITION_GE_2=.FALSE.
IF((bcyclicEndBlockProc(1)-bcyclicStartBlockProc(1)+1).GE.2) FIRSTPARTITION_GE_2=.TRUE.

LASTPARTITION_GE_2=.FALSE.
IF((bcyclicEndBlockProc(activeranks)-bcyclicStartBlockProc(activeranks)+1).GE.2) LASTPARTITION_GE_2=.TRUE.

numBlocks=bcyclicEndBlockProc(rank+1)-bcyclicStartBlockProc(rank+1)+1

!This may only occur on SOME processors, not all, so send to all
lnum = numBlocks.GE.2
CALL MPI_REDUCE(lnum, lnumg, 1, MPI_LOGICAL, MPI_LAND, 0, SIESTA_COMM,  &
                MPI_ERR)

CALL ASSERT(lnumg,'NS/nranks < 2...use smaller number of processors')

MAPPINGDONE=.TRUE.
END SUBROUTINE bcyclicMapping
!------------------------------------------------

!------------------------------------------------
! Compute AllGather vector variant parameters.
! This is to be called after bcyclicMapping.
!------------------------------------------------
SUBROUTINE computeAllGatherParameters

INTEGER :: i

CALL assert(MAPPINGDONE, 'computeAllGatherParameters')

IF(.NOT.ALLOCATED(rcounts)) ALLOCATE(rcounts(activeranks), rcountsNS(activeranks))
IF(.NOT.ALLOCATED(disp)) ALLOCATE(disp(activeranks),dispNS(activeranks))

DO i=1,activeranks
  rcountsNS(i)=(bcyclicEndBlockProc(i)-bcyclicStartBlockProc(i)+1)
  rcounts(i)  =(bcyclicEndBlockProc(i)-bcyclicStartBlockProc(i)+1)*M
END DO

disp(1)=0; dispNS(1)=0
DO i=2,activeranks
  dispNS(i)=dispNS(i-1)+rcountsNS(i-1)
  disp(i)=disp(i-1)+rcounts(i-1)
END DO

END SUBROUTINE computeAllGatherParameters
!------------------------------------------------

!------------------------------------------------
! A generic binary search routine
!------------------------------------------------
SUBROUTINE search(query, FOUND, location)

INTEGER, INTENT(IN)  :: query
LOGICAL, INTENT(OUT) :: FOUND
INTEGER, INTENT(OUT) :: location
INTEGER              :: p

CALL assert(MAPPINGDONE, 'search')

! To account for P>N
activeranks = activeranks
IF (activeranks.GT.N) THEN
   activeranks = N
END IF

! Dumb search algorithm
FOUND = .FALSE.
DO p = 1, activeranks
   IF ((bcyclicStartBlockProc(p) .LE. query) .AND. (query .LE. bcyclicEndBlockProc(p))) THEN
      FOUND = .TRUE.
      location = p
      EXIT
   END IF
END DO

END SUBROUTINE search 
!------------------------------------------------

!------------------------------------------------
SUBROUTINE send(columnData, blockRowNum, blockRowType, columnNum, procNum)

REAL(dp), DIMENSION(1:M), INTENT(IN) :: columnData
INTEGER, INTENT(IN)                  :: blockRowNum
INTEGER, INTENT(IN)                  :: blockRowType
INTEGER, INTENT(IN)                  :: columnNum
INTEGER, INTENT(OUT)                 :: procNum

CHARACTER, DIMENSION(PACKSIZE)       :: sendbuf
INTEGER                              :: positn
LOGICAL                              :: FOUND

positn=0
CALL MPI_Pack(blockRowNum, 1, MPI_INTEGER, sendBuf, PACKSIZE, positn, SIESTA_COMM, MPI_ERR)
CALL MPI_Pack(columnNum,   1, MPI_INTEGER, sendBuf, PACKSIZE, positn, SIESTA_COMM, MPI_ERR)
CALL MPI_Pack(columnData,  M, MPI_REAL8,   sendBuf, PACKSIZE, positn, SIESTA_COMM, MPI_ERR)

FOUND = .FALSE.
CALL search(blockRowNum, FOUND, procNum)
IF (FOUND) THEN
   IF (procNum-1 .EQ. rank) THEN
      IF (blockRowType .EQ. UPPER) THEN
         CALL SetMatrixRowColU(blockRowNum, columnData, columnNum)
      ELSE IF (blockRowType .EQ. SAVEDIAG) THEN
         CALL StoreDiagonal(blockRowNum, columnNum, columnData)
      ELSE IF (blockRowType .EQ. DIAG) THEN
         CALL SetMatrixRowColD(blockRowNum, columnData, columnNum)
      ELSE IF (blockRowType .EQ. LOWER) THEN
         CALL SetMatrixRowColL(blockRowNum, columnData, columnNum)
      ELSE
         CALL ASSERT(.FALSE., 'send error, blockRowType:')
      END IF
   ELSE 
      IF (BUFFERED) THEN
         CALL MPI_BSend(sendbuf, positn, MPI_PACKED, procNum-1, blockRowType, SIESTA_COMM, MPI_ERR)
      ELSE
         CALL MPI_Send(sendbuf, positn, MPI_PACKED, procNum-1, blockRowType, SIESTA_COMM, MPI_ERR)
      END IF
   END IF
ELSE
   CALL assert(MAPPINGDONE, 'send')
END IF

END SUBROUTINE send
!------------------------------------------------

!------------------------------------------------
SUBROUTINE receive(IPROBEFLAG)
USE blocktridiagonalsolver_s, ONLY:                                            &
& SetMatrixRowColL, SetMatrixRowColD, SetMatrixRowColU

LOGICAL, INTENT(IN)            :: IPROBEFLAG
LOGICAL                        :: FLAG
REAL(dp), DIMENSION(M)         :: coldata
CHARACTER, DIMENSION(PACKSIZE) :: recvbuf
INTEGER                        :: br
INTEGER                        :: ic
INTEGER                        :: it
INTEGER                        :: positn

FLAG = .TRUE.
IF (IPROBEFLAG) THEN
   CALL MPI_Iprobe(MPI_ANY_SOURCE, MPI_ANY_TAG, SIESTA_COMM, FLAG, siesta_MPI_Status, MPI_ERR)
END IF
IF (FLAG) THEN
   CALL MPI_Recv(recvbuf, PACKSIZE, MPI_PACKED, MPI_ANY_SOURCE, MPI_ANY_TAG, SIESTA_COMM, siesta_MPI_Status, MPI_ERR)
   nrecd = nrecd + 1

   it = siesta_MPI_Status(MPI_TAG)
   IF (it .NE. 1 .AND. it .NE. 2 .AND. it .NE. 3 .AND. it .NE. 4) THEN
      WRITE(TOFU,*) 'MPI_TAG:',it; FLUSH(TOFU)
      CALL ASSERT(.FALSE.,'receive 1 error')
   END IF
   positn = 0
   CALL MPI_Unpack(recvbuf, PACKSIZE, positn, br,      1, MPI_INTEGER, SIESTA_COMM, MPI_ERR)
   CALL MPI_Unpack(recvbuf, PACKSIZE, positn, ic,      1, MPI_INTEGER, SIESTA_COMM, MPI_ERR)
   CALL MPI_Unpack(recvbuf, PACKSIZE, positn, coldata, M, MPI_REAL8,   SIESTA_COMM, MPI_ERR)
   localbrow = br - bcyclicStartBlockProc(rank + 1) + 1

   IF (localbrow .LT. 1 .OR. localbrow .GT. numBlocks) THEN
      WRITE(TOFU,*) 'localbrow:',localbrow; FLUSH(TOFU)
      CALL ASSERT(.FALSE.,'receive 2 error')
   END IF

   IF (it .EQ. UPPER) THEN
      CALL SetMatrixRowColU(br, coldata, ic)
   ELSE IF (it .EQ. DIAG) THEN
      CALL SetMatrixRowColD(br, coldata, ic)
   ELSE IF (it .EQ. LOWER) THEN
      CALL SetMatrixRowColL(br, coldata, ic)
   ELSE IF (it .EQ. SAVEDIAG) THEN
      CALL StoreDiagonal(br, ic, coldata)
   ELSE
      WRITE(*,*)'Something wrong in ', rank,siesta_MPI_Status(MPI_TAG),br,ic,it
      CALL ASSERT(.FALSE.,'receive 3 error')
   END IF
END IF

END SUBROUTINE receive
!------------------------------------------------

!------------------------------------------------
! A utility function to gather all the partial 
! solution vectors and reconstruct the full 
! vector on each processor
!------------------------------------------------
SUBROUTINE GetFullSolution(invec,outvec) 

REAL(dp), DIMENSION(1:N*M), INTENT(IN) :: invec
REAL(dp), DIMENSION(N,M), INTENT(OUT)  :: outvec
INTEGER :: i, j, p
INTEGER :: numBlocks, cnt, indx, offset 

DO p=1, activeranks
  numBlocks=bcyclicEndBlockProc(p)-bcyclicStartBlockProc(p)+1
  DO i=1,numBlocks
    offset=(bcyclicStartBlockProc(p)-1)*M+i
    DO j=1,M
      outvec(bcyclicStartBlockProc(p)+i-1,j)=invec(offset+(j-1)*numBlocks)
    END DO
  END DO
END DO  

END SUBROUTINE GetFullSolution 
!------------------------------------------------

!-------------------------------------------------------------------------------
SUBROUTINE CheckPoint(ckpt, infname)

INTEGER, INTENT(IN) :: ckpt
CHARACTER(3) :: infname
      
CALL MPI_Barrier(SIESTA_COMM, MPI_ERR)
WRITE(90000+rank,*) 'CheckPoint: Location  ',infname, ' @ pt',ckpt,'in rank', rank
CALL FLUSH(90000+rank)

END SUBROUTINE CheckPoint
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
SUBROUTINE WriteTime(lfname, ofname, ltagname, tagname, usedtime, spass)

INTEGER, INTENT(IN) :: lfname, ltagname, spass
CHARACTER(len=lfname), INTENT(IN) :: ofname
CHARACTER(len=ltagname), INTENT(IN) :: tagname
DOUBLE PRECISION, INTENT(IN) :: usedtime

IF (spass.EQ.-1) THEN
  IF(SKSDBG) WRITE(TOFU,*) ofname," : ", tagname," : ", usedtime
  IF(SKSDBG) FLUSH(TOFU)
ELSE
  IF(SKSDBG) WRITE(TOFU,*) ofname," : ", tagname, " : ",usedtime, spass
  IF(SKSDBG) FLUSH(TOFU)
END IF

END SUBROUTINE WriteTime
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
SUBROUTINE SetUpTOMNSPAllGather()
INTEGER :: i

SetUpTOMNSPAllGatherDONE=.FALSE.

IF (.NOT.ALLOCATED(mnspcounts)) ALLOCATE (mnspcounts(nranks))
IF (.NOT.ALLOCATED(mnspdisps)) ALLOCATE (mnspdisps(nranks))

DO i=1,nranks
  mnspcounts(i)=(bcyclicEndBlockProc(i)-bcyclicStartBlockProc(i)+1)*linblksize
END DO

mnspdisps(1)=0
DO i=2,nranks
  mnspdisps(i)=mnspdisps(i-1)+mnspcounts(i-1)
END DO

SetUpTOMNSPAllGatherDONE=.TRUE.

END SUBROUTINE SetUpTOMNSPAllGather
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
SUBROUTINE GetTimes
USE timer_mod

! Get maximum
CALL MPI_Allreduce(total_time,total_time_max,1,MPI_REAL8,MPI_MAX,SIESTA_COMM,MPI_ERR)
CALL MPI_Allreduce(construct_hessian_time,construct_hessian_time_max,1,MPI_REAL8,MPI_MAX,SIESTA_COMM,MPI_ERR)
CALL MPI_Allreduce(asymmetry_check_time,asymmetry_check_time_max,1,MPI_REAL8,MPI_MAX,SIESTA_COMM,MPI_ERR)
CALL MPI_Allreduce(block_factorization_time,block_factorization_time_max,1,MPI_REAL8,MPI_MAX,SIESTA_COMM,MPI_ERR)
CALL MPI_Allreduce(hessian_funct_island_time,hessian_funct_island_time_max,1,MPI_REAL8,MPI_MAX,SIESTA_COMM,MPI_ERR)
CALL MPI_Allreduce(time_update_upperv,update_upperv_time_max,1,MPI_REAL8,MPI_MAX,SIESTA_COMM,MPI_ERR)
CALL MPI_Allreduce(time_init_state,init_state_time_max,1,MPI_REAL8,MPI_MAX,SIESTA_COMM,MPI_ERR)
CALL MPI_Allreduce(time_update_bfield,update_bfield_time_max,1,MPI_REAL8,MPI_MAX,SIESTA_COMM,MPI_ERR)
CALL MPI_Allreduce(time_update_pres,update_pres_time_max,1,MPI_REAL8,MPI_MAX,SIESTA_COMM,MPI_ERR)
CALL MPI_Allreduce(time_update_force,update_force_time_max,1,MPI_REAL8,MPI_MAX,SIESTA_COMM,MPI_ERR)
CALL MPI_Allreduce(cv_current_time,cv_currents_time_max,1,MPI_REAL8,MPI_MAX,SIESTA_COMM,MPI_ERR)
CALL MPI_Allreduce(get_force_harmonics_time,get_force_harmonics_time_max,1,MPI_REAL8,MPI_MAX,SIESTA_COMM,MPI_ERR)
CALL MPI_Allreduce(bhtobf_time,bhtobf_time_max,1,MPI_REAL8,MPI_MAX,SIESTA_COMM,MPI_ERR)
CALL MPI_Allreduce(tomnsp_time,tomnsp_time_max,1,MPI_REAL8,MPI_MAX,SIESTA_COMM,MPI_ERR)
CALL MPI_Allreduce(toijsp_time,toijsp_time_max,1,MPI_REAL8,MPI_MAX,SIESTA_COMM,MPI_ERR)
CALL MPI_Allreduce(to_full_mesh_time,to_full_mesh_time_max,1,MPI_REAL8,MPI_MAX,SIESTA_COMM,MPI_ERR)
CALL MPI_Allreduce(gmres_time,gmres_time_max,1,MPI_REAL8,MPI_MAX,SIESTA_COMM,MPI_ERR)
CALL MPI_Allreduce(gmres_wrap_time,gmres_wrap_time_max,1,MPI_REAL8,MPI_MAX,SIESTA_COMM,MPI_ERR)
CALL MPI_Allreduce(ParyAx_time,ParyAx_time_max,1,MPI_REAL8,MPI_MAX,SIESTA_COMM,MPI_ERR)
CALL MPI_Allreduce(sendrecv_time,sendrecv_time_max,1,MPI_REAL8,MPI_MAX,SIESTA_COMM,MPI_ERR)

! Get minimum
CALL MPI_Allreduce(time_total,total_time_min,1,MPI_REAL8,MPI_MIN,SIESTA_COMM,MPI_ERR)
CALL MPI_Allreduce(construct_hessian_time,construct_hessian_time_min,1,MPI_REAL8,MPI_MIN,SIESTA_COMM,MPI_ERR)
CALL MPI_Allreduce(asymmetry_check_time,asymmetry_check_time_min,1,MPI_REAL8,MPI_MAX,SIESTA_COMM,MPI_ERR)
CALL MPI_Allreduce(block_factorization_time,block_factorization_time_min,1,MPI_REAL8,MPI_MIN,SIESTA_COMM,MPI_ERR)
CALL MPI_Allreduce(hessian_funct_island_time,hessian_funct_island_time_min,1,MPI_REAL8,MPI_MIN,SIESTA_COMM,MPI_ERR)
CALL MPI_Allreduce(time_update_upperv,update_upperv_time_min,1,MPI_REAL8,MPI_MIN,SIESTA_COMM,MPI_ERR)
CALL MPI_Allreduce(time_init_state,init_state_time_min,1,MPI_REAL8,MPI_MIN,SIESTA_COMM,MPI_ERR)
CALL MPI_Allreduce(time_update_bfield,update_bfield_time_min,1,MPI_REAL8,MPI_MIN,SIESTA_COMM,MPI_ERR)
CALL MPI_Allreduce(time_update_pres,update_pres_time_min,1,MPI_REAL8,MPI_MIN,SIESTA_COMM,MPI_ERR)
CALL MPI_Allreduce(time_update_force,update_force_time_min,1,MPI_REAL8,MPI_MIN,SIESTA_COMM,MPI_ERR)
CALL MPI_Allreduce(cv_current_time,cv_currents_time_min,1,MPI_REAL8,MPI_MIN,SIESTA_COMM,MPI_ERR)
CALL MPI_Allreduce(get_force_harmonics_time,get_force_harmonics_time_min,1,MPI_REAL8,MPI_MIN,SIESTA_COMM,MPI_ERR)
CALL MPI_Allreduce(bhtobf_time,bhtobf_time_min,1,MPI_REAL8,MPI_MIN,SIESTA_COMM,MPI_ERR)
CALL MPI_Allreduce(tomnsp_time,tomnsp_time_min,1,MPI_REAL8,MPI_MIN,SIESTA_COMM,MPI_ERR)
CALL MPI_Allreduce(toijsp_time,toijsp_time_min,1,MPI_REAL8,MPI_MIN,SIESTA_COMM,MPI_ERR)
CALL MPI_Allreduce(to_full_mesh_time,to_full_mesh_time_min,1,MPI_REAL8,MPI_MIN,SIESTA_COMM,MPI_ERR)
CALL MPI_Allreduce(gmres_time,gmres_time_min,1,MPI_REAL8,MPI_MIN,SIESTA_COMM,MPI_ERR)
CALL MPI_Allreduce(gmres_wrap_time,gmres_wrap_time_min,1,MPI_REAL8,MPI_MIN,SIESTA_COMM,MPI_ERR)
CALL MPI_Allreduce(ParyAx_time,ParyAx_time_min,1,MPI_REAL8,MPI_MIN,SIESTA_COMM,MPI_ERR)
CALL MPI_Allreduce(sendrecv_time,sendrecv_time_min,1,MPI_REAL8,MPI_MIN,SIESTA_COMM,MPI_ERR)

WRITE (5000+rank, *)' Total runtime                                 :', time_total
WRITE (5000+rank, *)' *Time spent in init_data                      :', init_data_time
WRITE (5000+rank, *)' *Time spent in init_metric_elements           :', init_metric_elements_time
WRITE (5000+rank, *)' **Time spent in init_timers                   :', init_timers_time
WRITE (5000+rank, *)' **Time spent in read_wout_file                :', read_wout_file_time
WRITE (5000+rank, *)' **Time spent in Spline_Fourier_Modes          :', Spline_Fourier_Modes_time
WRITE (5000+rank, *)' **Time spent in Add_Ghost_Points              :', Add_Ghost_Points_time 
WRITE (5000+rank, *)' **Time spent in Spline_OneD_Array             :', Spline_OneD_Array_time
WRITE (5000+rank, *)' **Time spent in LoadRZL_VMEC                  :', LoadRZL_VMEC_time

WRITE (5000+rank, *)' *Time spent in test_fourier                   :', test_fourier_time
WRITE (5000+rank, *)' *Time spent in init_quantities                :', init_quantities_time
WRITE (5000+rank, *)' *Time spent in init_evolution                 :', init_evolution_time
WRITE (5000+rank, *)' *Time spent in converge_diagonal              :', converge_diagonal_time
WRITE (5000+rank, *)' *Time spent in comp_diag_elements             :', comp_diag_elements_time
WRITE (5000+rank, *)' *Time spent in converge_blocks                :', converge_blocks_time

WRITE (5000+rank, *)'   '
WRITE (5000+rank, *)' Time spent in converge_diagonal               :', converge_diagonal_time
WRITE (5000+rank, *)' *diag_evolve                                  :', diag_evolve_time
WRITE (5000+rank, *)' *diag_add_pert                                :', diag_add_pert_time
WRITE (5000+rank, *)' *Time spent in comp_diag_elements             :', comp_diag_elements_time

WRITE (5000+rank, *)'   '
WRITE (5000+rank, *)' Time spent in converge_blocks                 :', converge_blocks_time
WRITE (5000+rank, *)' *block_evolve                                 :', block_evolve_time
WRITE (5000+rank, *)' *block_add_pert                               :', block_add_pert_time

WRITE (5000+rank, *)'   '
WRITE (5000+rank, *)' Time spent in block_evolve                    :', block_evolve_time
WRITE (5000+rank, *)' *Time spent in Compute_Hessian_Blocks         :', compute_hessian_time
WRITE (5000+rank, *)' *Time spent in evolve_funct_island            :', evolve_funct_island_time
WRITE (5000+rank, *)' *Time spent in GMRES                          :', gmres_time
WRITE (5000+rank, *)' *Time spent in conj_grad                      :', conj_grad_time
WRITE (5000+rank, *)' *Time spent in evolve_restart_file_time       :', evolve_restart_file_time
WRITE (5000+rank, *)' *Time spent in evolve_add_resistive_E_time    :', evolve_add_resistive_E_time

WRITE (5000+rank, *)'   '
WRITE (5000+rank, *)' Time spent in Compute_Hessian_Blocks'
WRITE (5000+rank, *)' Time spent in Hessian_construction            :', construct_hessian_time
WRITE (5000+rank, *)' Time spent in Hessian_assymetry               :', asymmetry_check_time
WRITE (5000+rank, *)' Time spent in block_factorization             :', block_factorization_time

WRITE (5000+rank, *)'   '
WRITE (5000+rank, *)' Time spent in hess_funct_island               :', hessian_funct_island_time
WRITE (5000+rank, *)' *Time spent in update_upperv                  :', time_update_upperv
WRITE (5000+rank, *)' *Time spent in init_state                     :', time_init_state
WRITE (5000+rank, *)' *Time spent in update_bfield                  :', time_update_bfield
WRITE (5000+rank, *)' *Time spent in update_pres                    :', time_update_pres
WRITE (5000+rank, *)' *Time spent in update_force                   :', time_update_force

WRITE (5000+rank, *)'   '
WRITE (5000+rank, *)' Time spent in cv_currents                     :', cv_current_time
WRITE (5000+rank, *)' Time spent in get_force_harmonics             :', get_force_harmonics_time
WRITE (5000+rank, *)' Time spent in bhtobf                          :', bhtobf_time
WRITE (5000+rank, *)' Time spent in toijsp                          :', toijsp_time
WRITE (5000+rank, *)' Time spent in tomnsp                          :', tomnsp_time

WRITE (5000+rank, *)'   '
WRITE (5000+rank, *)' Time spent in GMRES                           :', gmres_time
WRITE (5000+rank, *)' Time spent in gmres_funct_island              :', gmres_funct_island_time
WRITE (5000+rank, *)' Time spent in gmres_init_dgmres               :', gmres_init_dgmres_time
WRITE (5000+rank, *)' Time spent in gmres_wrap                      :', gmres_wrap_time
WRITE (5000+rank, *)' Time spent in gmresr                          :', gmresr_time

WRITE (5000+rank, *)'   '
WRITE (5000+rank, *)' Time spent in drive_dgmres                    :', drive_dgmres_time
WRITE (5000+rank, *)' Time spent in ParyAx                          :', ParyAx_time
WRITE (5000+rank, *)' Time spent in dcopy                           :', dcopy_time
WRITE (5000+rank, *)' Time spent in apply_precond_time              :', time_apply_precon
WRITE (5000+rank, *)' Time spent in dgemv                           :', dgemv_time
WRITE (5000+rank, *)' Time spent in getnlforce                      :', getnlforce_time
WRITE (5000+rank, *)' Time spent in allgather                       :', gmres_wrap_allgather_time
WRITE (5000+rank, *)' Time spent in allreduce                       :', gmres_wrap_allreduce_time

END SUBROUTINE GetTimes
!-------------------------------------------------------------------------------

END MODULE nscalingtools
