      module prof_mod
        USE v3_utilities, ONLY: assert
        implicit none
        private

        integer, parameter :: maxroutines = 1024
        integer, parameter :: maxlevels = 1024
        integer, parameter :: maxstrlen = 127


        real, dimension(maxroutines) :: dictstart,dicttotal
        integer, dimension(maxroutines) :: dictcount
        integer :: nroutine,nlevels


        character(len=maxstrlen), dimension(maxroutines) :: dictname
        character(len=maxstrlen), dimension(maxlevels) :: lastroutine

        INTEGER, ALLOCATABLE, public :: scalcounts(:), scaldisps(:)
        INTEGER, ALLOCATABLE :: StartBlockProc(:), EndBlockProc(:)

#if defined(MPI_OPT)
        public :: profstart,profend,profstat,profinit,dclock
        public :: SetUpScalingAllGather

        contains

!=============================================
        real function dclock() 
        integer :: count,count_rate,count_max

        call system_clock(count,count_rate,count_max)
        if (count_rate.ne.0) then
           dclock = real(count)/real(count_rate)
        else
           dclock = 0.0
        endif

        end function dclock

!=============================================
        subroutine profinit()

!        integer :: i,ipos

        nroutine = 0

        dictname(:) = ' '
        dictstart(:) = 0.0
        dictcount(:) = 0.0
        dicttotal(:) = 0.0

        nlevels = 0
        lastroutine(:) = ' '

        end subroutine profinit

!=============================================
        subroutine profstart(rname)
        character(len=*),intent(in) :: rname


        character(len=maxstrlen) :: name
        logical :: found,isok
        integer :: i,j,ipos

        name = rname
        nlevels = nlevels + 1

        isok = (1 .le. nlevels).and.(nlevels .le. maxlevels)
        call assert( isok, '** profstart: invalid nlevels')
        
        lastroutine(nlevels) = name
        found = .false.
        do j=1,nroutine
           i = nroutine - j + 1
           if (dictname(i)(1:1).eq.name(1:1)) then
                found = (dictname(i) .eq. name)
                if (found) then
                   ipos = i
                   exit
                endif
           endif
        enddo

        if (.not.found) then
          nroutine = nroutine + 1
          isok = (nroutine .le. maxroutines)
          call assert(isok, '** profstart: nroutine > maxroutines')

          ipos = nroutine
          dictname(ipos) = name
          dictcount(ipos) = 0
          dicttotal(ipos) = 0.0
        endif

        dictstart(ipos) = dclock()
        dictcount(ipos) = dictcount(ipos) + 1

        return
        end subroutine profstart

!=============================================
        subroutine profend(rname)
        character(len=*),intent(in) :: rname

        character(len=maxstrlen) :: name
        integer :: i,j,ipos
        logical :: found,isok
        
        real :: tend

        name = rname
        tend = dclock()


        isok = (1.le.nlevels).and.(nlevels.le.maxlevels)
        call assert(isok, '** profend: invalid nlevels')

        isok = (name .eq. lastroutine(nlevels))
        if (.not.isok) then
          print*,'** profend name != lastroutine(',nlevels,') '
          print*,'name: ', name
          print*,'lastroutine(nlevels): ', lastroutine(nlevels)

          stop '** error ** '
        endif

        found = .false.
        do j=1,nroutine
           i = nroutine - j + 1

           if (dictname(i)(1:1) .eq. name(1:1)) then
                found = (dictname(i) .eq. name)
                if (found) then
                        ipos = i
                        exit
                endif
           endif
        enddo

        if (.not.found) then
                print*,'** profend: routine name not found '
                print*,'name: ',name
                stop '** error ** '
        endif

        dicttotal(ipos) = dicttotal(ipos) + (tend - dictstart(ipos));
        nlevels = nlevels - 1;
        
        return
        end subroutine profend

!=============================================
        subroutine profstat(outdev_in)
        implicit none
        integer, optional, intent(in):: outdev_in 
        character(len=maxstrlen) :: fname,fstr
        integer :: i, outdev

        if (nroutine .le. 0) return

        if (present(outdev_in)) then
           outdev = outdev_in
        else
           outdev = 16
        endif

        fname = 'profstat.dat'
        open(outdev, file=fname, form='formatted',                            &
     &       access='sequential',status='unknown')
        rewind(outdev)

        fstr = "(A20,' was called ',i10,' times, total ',f10.2,' secs')"
        do i=1,nroutine
          write(outdev,fstr) dictname(i), dictcount(i), dicttotal(i)
          write(*,fstr) dictname(i), dictcount(i), dicttotal(i)
        enddo

        close(outdev)

        IF (ALLOCATED(scalcounts)) DEALLOCATE (scalcounts, scaldisps)

        return
        end subroutine profstat
        
!Add SPH for col/row scaling
        SUBROUTINE SetUpScalingAllGather(mblk_size)
        USE descriptor_mod
        INTEGER, INTENT(IN) :: mblk_size
        INTEGER :: i, numroc
        EXTERNAL :: numroc, blacs_gridinfo

        IF (.NOT.ALLOCATED(scalcounts)) ALLOCATE (scalcounts(nprocs))
        IF (.NOT.ALLOCATED(scaldisps)) ALLOCATE (scaldisps(nprocs))

        CALL numrocMapping(iam, nprocs, mblk_size)

        DO i=1,nprocs
           scalcounts(i)=(EndBlockProc(i)-StartBlockProc(i)+1)
        END DO

        DEALLOCATE (StartBlockProc, EndBlockProc)

!Sanity consistency check        
!!!!!!!!!!!!!!!!!!!!!!!!!
        CALL blacs_gridinfo(icontxt_1xp,nprow,npcol,myrow,mycol)
        mb = mblk_size
        nb = 1

        rsrc = 0
        csrc = 0
        Locq = numroc( mblk_size, nb, mycol, csrc, npcol )
!        Locp = numroc( mblk_size, mb, myrow, rsrc, nprow )
        mblk_size2 = MAX(1,Locq)
        CALL ASSERT(scalcounts(iam+1).EQ.mblk_size2,                    &
        'scalcounts != mblk_size2 in SetupScalingAllGather')
!!!!!!!!!!!!!!!!!!!!!!!!!

        scaldisps(1)=0
        DO i=2,nprocs
           scaldisps(i)=scaldisps(i-1)+scalcounts(i-1)
        END DO

        END SUBROUTINE SetUpScalingAllGather

        SUBROUTINE numrocMapping(rank, activeranks, N)
        use descriptor_mod, ONLY: SIESTA_COMM
        use mpi_inc
        INTEGER, INTENT(IN) :: rank, activeranks, N
        INTEGER :: startblock, endblock
        INTEGER :: lload, sload, myload
        INTEGER :: numL, numS, mpi_err
        INTEGER :: r, c

        IF (.NOT.ALLOCATED(StartBlockProc)) ALLOCATE (StartBlockProc(activeranks))
        IF (.NOT.ALLOCATED(EndBlockProc)) ALLOCATE (EndBlockProc(activeranks))

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

        CALL MPI_Allgather(startblock,1,MPI_INTEGER,StartBlockProc,     &
                           1,MPI_INTEGER,SIESTA_COMM,MPI_ERR)

        CALL MPI_Allgather(endblock,1,MPI_INTEGER,EndBlockProc,         &
                           1,MPI_INTEGER,SIESTA_COMM,MPI_ERR)


        END SUBROUTINE numrocMapping
#endif
      end module prof_mod

