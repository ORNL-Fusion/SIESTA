      module ptrd_mod
        use stel_kinds
        USE v3_utilities, ONLY: assert
        implicit none

#if defined(MPI_OPT)
        contains

        subroutine pdtrdf(m,nblock,Amat,Bmat,Cmat,ipivmat, desc)
        implicit none
!
!       perform LU factorization of block tridiagonal system
!
!       [A1  C1       ]   [ L1          ] [ I U1        ]
!       [B1  A2 C2    ]   [ B1 L2       ] [   I  U2     ]
!       [    B2 A3 C3 ] = [    B2 L3    ] [      I   U3 ]
!       [       B4 A4 ]   [       B3 L4 ] [          I  ]
!
! 
! L1 = A1
! L1*U1 = C1 => U1 = L1\C1
! B1*U1 + L2 = A2 =>  L2 = A2 - B1*U1
! L2*U2 = C2 =>  U2 = L2\C2
!
! Lk = Ak - B_{k-1}*U_{k-1}
! Uk = Lk\Ck
! 
!
        REAL(rprec), parameter :: one = 1.0d0
        REAL(rprec), parameter :: zero = 0.0d0

        integer, intent(in) :: m, nblock
        REAL(rprec), target, dimension(:,:), intent(inout) :: Amat
        REAL(rprec), target, dimension(:,:), intent(inout) :: Bmat
        REAL(rprec), target, dimension(:,:), intent(inout) :: Cmat
        integer, dimension(:,:), target, intent(inout) :: ipivmat
        integer, dimension(:), intent(in) :: desc


        integer :: k, ia,ja, ib,jb, iu,ju, ic,jc
        integer :: info, nrhs, mm,nn
        integer, dimension(:), pointer :: ipiv
        REAL(rprec) :: alpha, beta
        REAL(rprec), dimension(:), pointer :: Ak, Bkm1, Ukm1
        REAL(rprec), dimension(:), pointer :: Lk, Ck


        nullify( Ak )
        nullify( Bkm1 )
        nullify( Ukm1 )
        nullify( Lk )
        nullify( Ck )


        do k=1,nblock
          if (k .ge. 2) then
!            -------------------------
!            Ak = Ak - B_{k-1}*U_{k-1}
!            -------------------------
             Ak => Amat(:,k)
             Bkm1 => Bmat(:,k-1)
             Ukm1 => Cmat(:,k-1)
             alpha = -one
             beta = one
             ib = 1
             jb = 1
             iu = 1
             ju = 1
             ia = 1
             ja = 1

             call pdgemm( 'N', 'N', m,m,m,                                    &
     &                   alpha, Bkm1,ib,jb,desc,                              &
     &                          Ukm1,iu,ju,desc,                              &
     &                   beta,  Ak, ia,ja, desc )
          endif

!         ------------
!         Uk = Lk \ Ck
!         ------------
          Lk => Amat(:,k)
          ipiv => ipivmat(:,k)
          ia = 1
          ja = 1
          info = 0
          mm = m
          nn = m
          call pdgetrf(mm,nn,Lk,ia,ja,desc,ipiv,info)
          call assert(info.eq.0,'pdtrdf: pdgetrf return info != 0')

          if (k .le. (nblock-1)) then
             nrhs = m
             ic = 1
             jc = 1
             Ck => Cmat(:,k)
             call pdgetrs('No',m,nrhs,Lk,ia,ja,desc,ipiv,               &
     &                    Ck,ic,jc,desc, info)
             call assert(info.eq.0,'pdtrdf: pdgetrs return info != 0')
          endif
          

        enddo
        return
        end subroutine pdtrdf

        subroutine pdtrds(m,nblock,Amat,Bmat,Cmat,ipivmat,desc,         &
     &                    nrhs, Rrhs_in,ir,jr,descRrhs_in)
        use descriptor_mod, ir1=>ir, jr1=>jr, mm1=>mm, info1=>info
        implicit none


!
!       use LU factorization of block tridiagonal system
!
!       [A1  C1       ]   [ L1          ] [ I U1        ]
!       [B1  A2 C2    ]   [ B1 L2       ] [   I  U2     ]
!       [    B2 A3 C3 ] = [    B2 L3    ] [      I   U3 ]
!       [       B4 A4 ]   [       B3 L4 ] [          I  ]
!
! 
! 
!
        integer, parameter :: idebug = 0
        REAL(rprec), parameter :: one = 1.0d0
        REAL(rprec), parameter :: zero = 0.0d0

        integer, intent(in) :: m, nblock
        REAL(rprec), dimension(:,:), target, intent(inout) :: Amat
        REAL(rprec), dimension(:,:), target, intent(inout) :: Bmat
        REAL(rprec), dimension(:,:), target, intent(inout) :: Cmat
        integer, dimension(:,:), target, intent(inout) :: ipivmat
        integer, dimension(:), intent(in) :: desc

        integer, intent(in) :: nrhs
        integer, intent(in) :: ir,jr
        REAL(rprec), dimension(:) :: Rrhs_in
        integer, dimension(:), intent(in) :: descRrhs_in


        integer :: k, info, mm,nn,kk
        integer :: ia,ja, irk, jrk, iy,jy, ix,jx, iu,ju, ib,jb
        REAL(rprec) :: alpha, beta

        REAL(rprec), pointer, dimension(:) :: Lk,Bkm1,Uk
        REAL(rprec), pointer, dimension(:) :: Rrhsk,Rrhskm1,Rrhskp1
        integer, pointer, dimension(:) :: ipiv


	integer :: inc1,inc2,iblock,irhs

	logical, parameter :: use_pdcopy = .false.
	REAL(rprec), dimension(:,:), target, allocatable :: Rrhs
	integer, dimension(DLEN_) :: descRrhs
    INTEGER           :: NUMROC
    EXTERNAL          :: NUMROC

!	-----------------------------------------------------
!	change storage to nblock copies of m by nrhs matrices
!	also to avoid alignment constraints in scalapack
!	-----------------------------------------------------
	icontxt = descRrhs_in(CTXT_)
	mb = descRrhs_in(MB_)
	nb = 1
	rsrc = 0
	csrc = 0

	call blacs_gridinfo( icontxt, nprow,npcol,myrow,mycol)
	Locp = numroc( m, mb, myrow,rsrc,nprow)
	Locq = numroc( nrhs,nb,mycol,csrc,npcol)
	lld = max(1,Locp)
	call descinit(descRrhs,m,nrhs,mb,nb,rsrc,csrc,icontxt,lld,info)
	call assert(info.eq.0,'pdtrds: descinit return info != 0')

	ineed = max(1,Locp*Locq)
	allocate( Rrhs(ineed,nblock),stat=ierr)
	call assert(ierr.eq.0,'pdtrds: alloc Rrhs,ierr != 0')
        Rrhs = 0.0d0

	nullify( Rrhsk )
	nullify( Rrhskm1 )
	nullify( Rrhskp1 )

!	---------
!	copy data
!	---------

	if (use_pdcopy) then

	do iblock=1,nblock
	do irhs=1,nrhs
	  ia = (ir-1) + 1 + (iblock-1)*m
	  ja = (jr-1) + irhs
	  inc1 = 1

	  ib = 1
	  jb = irhs
	  inc2 = 1

	  Rrhsk => Rrhs(:,iblock)
	  call pdcopy(m, Rrhs_in,ia,ja,descRrhs_in,inc1,                    &
     &                   Rrhsk,ib,jb,descRrhs,inc2 )
	enddo
	enddo

	else

	do iblock=1,nblock
	  irhs = 1
	  ia = (ir-1) + 1 + (iblock-1)*m
	  ja = (jr-1) + irhs
	  
	  ib = 1
	  jb = irhs
	  Rrhsk => Rrhs(:,iblock)
	  alpha = 1.0d0
	  beta = 0.0d0
	  call pdgeadd( 'N',m,nrhs,alpha,Rrhs_in,ia,ja,descRrhs_in,           &
     &            beta, Rrhsk,ib,jb,descRrhs )
	enddo

	endif

!   ----------------- 
!   L*U*x = r
!   (1) solve L*y = r
!   (2) solve U*x = y
!   ----------------- 
        nullify(Lk)
        nullify(Bkm1)
        nullify(Uk)
        nullify(ipiv)

        icontxt = desc(CTXT_)
	call blacs_gridinfo(icontxt, nprow,npcol,myrow,mycol)
        isroot = (myrow.eq.0).and.(mycol.eq.0)

!   --------------------------------
!  (1) solve L*y = r
!
!   [L1             ] [ y1 ]   [ r1 ]
!   [B1  L2         ] [ y2 ]   [ r2 ]
!   [    B2  L3     ] [ y3 ] = [ r3 ]
!   [        B3  L4 ] [ y4 ]   [ r4 ]
!
!
!   y1 = L1\r1
!   y2 = L2\( r2 - B1*y1 )
!   y3 = L3\( r3 - B2*y2 )
!   y4 = L4\( r4 - B3*y3 )
!
!   yk = Lk\( rk - B_{k-1}*y_{k-1} )
!   --------------------------------
          if (isroot .and. (idebug.ge.1)) then
             write(*,*) 'pdtrds 77: m,nblock,nrhs ',m,nblock,nrhs
             write(*,*) 'descRrhs_in(M_) ',descRrhs_in(M_)
             write(*,*) 'descRrhs_in(N_) ',descRrhs_in(N_)
             write(*,*) 'descRrhs_in(MB_) ',descRrhs_in(MB_)
             write(*,*) 'descRrhs_in(NB_) ',descRrhs_in(NB_)
          endif

        do k=1,nblock
          if (k .ge. 2) then
!           --------------------------
!           rk <- rk - B_{k-1}*y_{k-1}  
!           --------------------------
            Bkm1 => Bmat(:,k-1)

            alpha = -one
            beta = one
            mm = m
            nn = nrhs
            kk = m
            ib = 1
            jb = 1

            iy = 1
            jy = 1
            irk = 1
            jrk = 1
	    Rrhsk => Rrhs(:,k)
	    Rrhskm1 => Rrhs(:,k-1)

            call pdgemm( 'N', 'N', mm,nn,kk,                                 &
     &          alpha,  Bkm1, ib,jb, desc,                                   &
     &                  Rrhskm1, iy,jy, descRrhs,                            &
     &          beta,   Rrhsk, irk,jrk,descRrhs )
          endif

!         ------------
!         yk = Lk \ rk
!         ------------

          Lk => Amat(:,k)
          ia = 1
          ja = 1
          ipiv => ipivmat(:,k)

	  Rrhsk => Rrhs(:,k)
          jrk = 1
          irk = 1

          info = 0
          if (isroot .and. (idebug.ge.1)) then
             write(*,*) 'pdtrds 106: k,irk,jrk ',k,irk,jrk
          endif
          call pdgetrs( 'N', m,nrhs, Lk,ia,ja,desc, ipiv,               &
     &                  Rrhsk,irk,jrk,descRrhs,info)
          call assert(info.eq.0,'pdtrds: pdgetrs return info != 0')

        enddo


!  (2)  solve U*x = y
!
! [I  U1         ]   [ x1 ]   [ y1 ]
! [   I   U2     ]   [ x2 ]   [ y2 ]
! [       I   U3 ]   [ x3 ] = [ y3 ]
! [           I  ]   [ x4 ]   [ y4 ]
!
!
! x4 = y4
! x3 = y3 - U3*y4
! x2 = y2 - U2*y3
! x1 = y1 - U1*y2
!
! xk = yk - Uk*x_{k+1}
!

        do k=nblock-1,1,-1
!          --------------------
!          xk = yk - Uk*x_{k+1}
!          --------------------
           alpha = -one
           beta = one

           mm = m
           nn = nrhs
           kk = m
           Uk => Cmat(:,k)
           iu = 1
           ju = 1

           ix = 1
           jx = 1
           irk = 1
           jrk = 1

	   Rrhskp1 => Rrhs(:,k+1)
	   Rrhsk => Rrhs(:,k)

           call pdgemm( 'N','N', mm,nn,kk,                                     &
     &                 alpha,  Uk,iu,ju,desc,                                  &
     &                         Rrhskp1,ix,jx,descRrhs,                         &
     &                 beta,   Rrhsk,irk,jrk,descRrhs )
         enddo

!	----------------
!	copy results out
!	----------------

	if (use_pdcopy) then

	do iblock=1,nblock
	do irhs=1,nrhs
	  ia = (ir-1) + 1 + (iblock-1)*m
	  ja = (jr-1) + irhs
	  inc1 = 1

	  ib = 1
	  jb = irhs
	  inc2 = 1

	  Rrhsk => Rrhs(:,iblock)
          call pdcopy(m,Rrhsk,ib,jb,descRrhs,inc2,                          &
     &                  Rrhs_in,ia,ja,descRrhs_in,inc1)
	enddo
	enddo

	else

	do iblock=1,nblock
	  irhs = 1
	  ia = (ir-1) + 1 + (iblock-1)*m
	  ja = (jr-1) + irhs
	  
	  ib = 1
	  jb = irhs
	  Rrhsk => Rrhs(:,iblock)
	  alpha = 1.0d0
	  beta = 0.0d0

	  call pdgeadd('N',m,nrhs,alpha,Rrhsk,ib,jb,descRrhs,               &
     &                            beta,Rrhs_in,ia,ja,descRrhs_in)
	enddo

	endif
	deallocate( Rrhs, stat=ierr)
	call assert(ierr.eq.0,'pdtrds: dealloc Rhs ')
	  

        return
        end subroutine pdtrds

#if defined(NEED_TOOLS)      
      SUBROUTINE PDGETRS( TRANS, N, NRHS, A, IA, JA, DESCA, IPIV, B,    &
     &                    IB, JB, DESCB, INFO )
!
!  -- ScaLAPACK routine (version 1.7) --
!     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
!     and University of California, Berkeley.
!     May 1, 1997
!
!     .. Scalar Arguments ..
      CHARACTER          TRANS
      INTEGER            IA, IB, INFO, JA, JB, N, NRHS
!     ..
!     .. Array Arguments ..
      INTEGER            DESCA( * ), DESCB( * ), IPIV( * )
      REAL(dp)   A( * ), B( * )
!     ..
!
!  Purpose
!  =======
!
!  PDGETRS solves a system of distributed linear equations
!
!                   op( sub( A ) ) * X = sub( B )
!
!  with a general N-by-N distributed matrix sub( A ) using the LU
!  factorization computed by PDGETRF.
!  sub( A ) denotes A(IA:IA+N-1,JA:JA+N-1), op( A ) = A or A**T and
!  sub( B ) denotes B(IB:IB+N-1,JB:JB+NRHS-1).
!
!  Notes
!  =====
!
!  Each global data object is described by an associated description
!  vector.  This vector stores the information required to establish
!  the mapping between an object element and its corresponding process
!  and memory location.
!
!  Let A be a generic term for any 2D block cyclicly distributed array.
!  Such a global array has an associated description vector DESCA.
!  In the following comments, the character _ should be read as
!  "of the global array".
!
!  NOTATION        STORED IN      EXPLANATION
!  --------------- -------------- --------------------------------------
!  DTYPE_A(global) DESCA( DTYPE_ )The descriptor type.  In this case,
!                                 DTYPE_A = 1.
!  CTXT_A (global) DESCA( CTXT_ ) The BLACS context handle, indicating
!                                 the BLACS process grid A is distribu-
!                                 ted over. The context itself is glo-
!                                 bal, but the handle (the integer
!                                 value) may vary.
!  M_A    (global) DESCA( M_ )    The number of rows in the global
!                                 array A.
!  N_A    (global) DESCA( N_ )    The number of columns in the global
!                                 array A.
!  MB_A   (global) DESCA( MB_ )   The blocking factor used to distribute
!                                 the rows of the array.
!  NB_A   (global) DESCA( NB_ )   The blocking factor used to distribute
!                                 the columns of the array.
!  RSRC_A (global) DESCA( RSRC_ ) The process row over which the first
!                                 row of the array A is distributed.
!  CSRC_A (global) DESCA( CSRC_ ) The process column over which the
!                                 first column of the array A is
!                                 distributed.
!  LLD_A  (local)  DESCA( LLD_ )  The leading dimension of the local
!                                 array.  LLD_A >= MAX(1,LOCr(M_A)).
!
!  Let K be the number of rows or columns of a distributed matrix,
!  and assume that its process grid has dimension p x q.
!  LOCr( K ) denotes the number of elements of K that a process
!  would receive if K were distributed over the p processes of its
!  process column.
!  Similarly, LOCc( K ) denotes the number of elements of K that a
!  process would receive if K were distributed over the q processes of
!  its process row.
!  The values of LOCr() and LOCc() may be determined via a call to the
!  ScaLAPACK tool function, NUMROC:
!          LOCr( M ) = NUMROC( M, MB_A, MYROW, RSRC_A, NPROW ),
!          LOCc( N ) = NUMROC( N, NB_A, MYCOL, CSRC_A, NPCOL ).
!  An upper bound for these quantities may be computed by:
!          LOCr( M ) <= ceil( ceil(M/MB_A)/NPROW )*MB_A
!          LOCc( N ) <= ceil( ceil(N/NB_A)/NPCOL )*NB_A
!
!  This routine requires square block data decomposition ( MB_A=NB_A ).
!
!  Arguments
!  =========
!
!  TRANS   (global input) CHARACTER
!          Specifies the form of the system of equations:
!          = 'N':  sub( A )    * X = sub( B )  (No transpose)
!          = 'T':  sub( A )**T * X = sub( B )  (Transpose)
!          = 'C':  sub( A )**T * X = sub( B )  (Transpose)
!
!  N       (global input) INTEGER
!          The number of rows and columns to be operated on, i.e. the
!          order of the distributed submatrix sub( A ). N >= 0.
!
!  NRHS    (global input) INTEGER
!          The number of right hand sides, i.e., the number of columns
!          of the distributed submatrix sub( B ). NRHS >= 0.
!
!  A       (local input) REAL(dp) pointer into the local
!          memory to an array of dimension (LLD_A, LOCc(JA+N-1)).
!          On entry, this array contains the local pieces of the factors
!          L and U from the factorization sub( A ) = P!L!U; the unit
!          diagonal elements of L are not stored.
!
!  IA      (global input) INTEGER
!          The row index in the global array A indicating the first
!          row of sub( A ).
!
!  JA      (global input) INTEGER
!          The column index in the global array A indicating the
!          first column of sub( A ).
!
!  DESCA   (global and local input) INTEGER array of dimension DLEN_.
!          The array descriptor for the distributed matrix A.
!
!  IPIV    (local input) INTEGER array, dimension ( LOCr(M_A)+MB_A )
!          This array contains the pivoting information.
!          IPIV(i) -> The global row local row i was swapped with.
!          This array is tied to the distributed matrix A.
!
!  B       (local input/local output) REAL(dp) pointer into the
!          local memory to an array of dimension
!          (LLD_B,LOCc(JB+NRHS-1)).  On entry, the right hand sides
!          sub( B ). On exit, sub( B ) is overwritten by the solution
!          distributed matrix X.
!
!  IB      (global input) INTEGER
!          The row index in the global array B indicating the first
!          row of sub( B ).
!
!  JB      (global input) INTEGER
!          The column index in the global array B indicating the
!          first column of sub( B ).
!
!  DESCB   (global and local input) INTEGER array of dimension DLEN_.
!          The array descriptor for the distributed matrix B.
!
!  INFO    (global output) INTEGER
!          = 0:  successful exit
!          < 0:  If the i-th argument is an array and the j-entry had
!                an illegal value, then INFO = -(i!100+j), if the i-th
!                argument is a scalar and had an illegal value, then
!                INFO = -i.
!
!  =====================================================================
!
!     .. Parameters ..
      INTEGER            BLOCK_CYCLIC_2D, CSRC_, CTXT_, DLEN_, DTYPE_,   &
     &                   LLD_, MB_, M_, NB_, N_, RSRC_
      PARAMETER          ( BLOCK_CYCLIC_2D = 1, DLEN_ = 9, DTYPE_ = 1,   &
     &                     CTXT_ = 2, M_ = 3, N_ = 4, MB_ = 5, NB_ = 6,  &
     &                     RSRC_ = 7, CSRC_ = 8, LLD_ = 9 )
      REAL(dp)   ONE
      PARAMETER          ( ONE = 1.0D+0 )
!     ..
!     .. Local Scalars ..
      LOGICAL            NOTRAN
      INTEGER            IAROW, IBROW, ICOFFA, ICTXT, IROFFA, IROFFB,    &
     &                   MYCOL, MYROW, NPCOL, NPROW
!     ..
!     .. Local Arrays ..
      INTEGER            DESCIP( DLEN_ ), IDUM1( 1 ), IDUM2( 1 )
!     ..
!     .. External Subroutines ..
      EXTERNAL           BLACS_GRIDINFO, CHK1MAT, DESCSET, PCHK2MAT,     &
     &                   PDLAPIV, PDTRSM, PXERBLA
!     ..
!     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            INDXG2P, NUMROC
      EXTERNAL           INDXG2P, LSAME, NUMROC
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          ICHAR, MOD
!     ..
!     .. Executable Statements ..
!
!     Get grid parameters
!
      ICTXT = DESCA( CTXT_ )
      CALL BLACS_GRIDINFO( ICTXT, NPROW, NPCOL, MYROW, MYCOL )
!
!     Test the input parameters
!
      INFO = 0
      IF( NPROW.EQ.-1 ) THEN
         INFO = -(700+CTXT_)
      ELSE
         NOTRAN = LSAME( TRANS, 'N' )
         CALL CHK1MAT( N, 2, N, 2, IA, JA, DESCA, 7, INFO )
         CALL CHK1MAT( N, 2, NRHS, 3, IB, JB, DESCB, 12, INFO )
         IF( INFO.EQ.0 ) THEN
            IAROW = INDXG2P( IA, DESCA( MB_ ), MYROW, DESCA( RSRC_ ),   &
     &                       NPROW )
            IBROW = INDXG2P( IB, DESCB( MB_ ), MYROW, DESCB( RSRC_ ),   &
     &                       NPROW )
            IROFFA = MOD( IA-1, DESCA( MB_ ) )
            ICOFFA = MOD( JA-1, DESCA( NB_ ) )
            IROFFB = MOD( IB-1, DESCB( MB_ ) )
            IF( .NOT.NOTRAN .AND. .NOT.LSAME( TRANS, 'T' ) .AND. .NOT.  &
     &         LSAME( TRANS, 'C' ) ) THEN
               INFO = -1
            ELSE IF( IROFFA.NE.0 ) THEN
               INFO = -5
            ELSE IF( ICOFFA.NE.0 ) THEN
               INFO = -6
            ELSE IF( DESCA( MB_ ).NE.DESCA( NB_ ) ) THEN
               INFO = -(700+NB_)
!            ELSE IF( IROFFB.NE.0 .OR. IBROW.NE.IAROW ) THEN
!               INFO = -10
            ELSE IF( DESCB( MB_ ).NE.DESCA( NB_ ) ) THEN
               INFO = -(1200+NB_)
            ELSE IF( ICTXT.NE.DESCB( CTXT_ ) ) THEN
               INFO = -(1200+CTXT_)
            END IF
         END IF
         IF( NOTRAN ) THEN
            IDUM1( 1 ) = ICHAR( 'N' )
         ELSE IF( LSAME( TRANS, 'T' ) ) THEN
            IDUM1( 1 ) = ICHAR( 'T' )
         ELSE
            IDUM1( 1 ) = ICHAR( 'C' )
         END IF
         IDUM2( 1 ) = 1
         CALL PCHK2MAT( N, 2, N, 2, IA, JA, DESCA, 7, N, 2, NRHS, 3,    &
     &                  IB, JB, DESCB, 12, 1, IDUM1, IDUM2, INFO )
      END IF

      IF( INFO.NE.0 ) THEN
         CALL PXERBLA( ICTXT, 'PDGETRS', -INFO )
         RETURN
      END IF
!
!     Quick return if possible
!
      IF( N.EQ.0 .OR. NRHS.EQ.0 ) RETURN

      CALL DESCSET( DESCIP, DESCA( M_ ) + DESCA( MB_ )*NPROW, 1,        &
     &              DESCA( MB_ ), 1, DESCA( RSRC_ ), MYCOL, ICTXT,      &
     &              DESCA( MB_ ) + NUMROC( DESCA( M_ ), DESCA( MB_ ),   &
     &              MYROW, DESCA( RSRC_ ), NPROW ) )

      IF( NOTRAN ) THEN
!
!        Solve sub( A ) * X = sub( B ).
!
!        Apply row interchanges to the right hand sides.
!
         CALL PDLAPIV( 'Forward', 'Row', 'Col', N, NRHS, B, IB, JB,     &
     &                 DESCB, IPIV, IA, 1, DESCIP, IDUM1 )
!
!        Solve L*X = sub( B ), overwriting sub( B ) with X.
!
         CALL PDTRSM( 'Left', 'Lower', 'No transpose', 'Unit', N, NRHS, &
     &                ONE, A, IA, JA, DESCA, B, IB, JB, DESCB )
!
!        Solve U*X = sub( B ), overwriting sub( B ) with X.
!
         CALL PDTRSM( 'Left', 'Upper', 'No transpose', 'Non-unit', N,   &
     &                NRHS, ONE, A, IA, JA, DESCA, B, IB, JB, DESCB )
      ELSE
!
!        Solve sub( A )' * X = sub( B ).
!
!        Solve U'*X = sub( B ), overwriting sub( B ) with X.
!
         CALL PDTRSM( 'Left', 'Upper', 'Transpose', 'Non-unit', N, NRHS, &
     &                ONE, A, IA, JA, DESCA, B, IB, JB, DESCB )
!
!        Solve L'*X = sub( B ), overwriting sub( B ) with X.
!
         CALL PDTRSM( 'Left', 'Lower', 'Transpose', 'Unit', N, NRHS,    &
     &                ONE, A, IA, JA, DESCA, B, IB, JB, DESCB )
!
!        Apply row interchanges to the solution vectors.
!
         CALL PDLAPIV( 'Backward', 'Row', 'Col', N, NRHS, B, IB, JB,    &
     &                 DESCB, IPIV, IA, 1, DESCIP, IDUM1 )

      END IF

      RETURN

!     End of PDGETRS

      END SUBROUTINE PDGETRS
      
      SUBROUTINE PDELGET( SCOPE, TOP, ALPHA, A, IA, JA, DESCA )
      use descriptor_mod, desca_x=>desca
      implicit none
!
!  -- ScaLAPACK tools routine (version 1.7) --
!     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
!     and University of California, Berkeley.
!     May 1, 1997
!
!     .. Scalar Arguments ..
      CHARACTER*1        SCOPE, TOP
      INTEGER            IA, JA
      REAL(dp)   ALPHA
!     ..
!     .. Array arguments ..
      INTEGER            DESCA( * )
      REAL(dp)   A( * )
!     ..
!
!  Purpose
!  =======
!
!  PDELGET sets alpha to the distributed matrix entry A( IA, JA ).
!  The value of alpha is set according to the scope.
!
!  Notes
!  =====
!
!  Each global data object is described by an associated description
!  vector.  This vector stores the information required to establish
!  the mapping between an object element and its corresponding process
!  and memory location.
!
!  Let A be a generic term for any 2D block cyclicly distributed array.
!  Such a global array has an associated description vector DESCA.
!  In the following comments, the character _ should be read as
!  "of the global array".
!
!  NOTATION        STORED IN      EXPLANATION
!  --------------- -------------- --------------------------------------
!  DTYPE_A(global) DESCA( DTYPE_ )The descriptor type.  In this case,
!                                 DTYPE_A = 1.
!  CTXT_A (global) DESCA( CTXT_ ) The BLACS context handle, indicating
!                                 the BLACS process grid A is distribu-
!                                 ted over. The context itself is glo-
!                                 bal, but the handle (the integer
!                                 value) may vary.
!  M_A    (global) DESCA( M_ )    The number of rows in the global
!                                 array A.
!  N_A    (global) DESCA( N_ )    The number of columns in the global
!                                 array A.
!  MB_A   (global) DESCA( MB_ )   The blocking factor used to distribute
!                                 the rows of the array.
!  NB_A   (global) DESCA( NB_ )   The blocking factor used to distribute
!                                 the columns of the array.
!  RSRC_A (global) DESCA( RSRC_ ) The process row over which the first
!                                 row of the array A is distributed.
!  CSRC_A (global) DESCA( CSRC_ ) The process column over which the
!                                 first column of the array A is
!                                 distributed.
!  LLD_A  (local)  DESCA( LLD_ )  The leading dimension of the local
!                                 array.  LLD_A >= MAX(1,LOCr(M_A)).
!
!  Let K be the number of rows or columns of a distributed matrix,
!  and assume that its process grid has dimension p x q.
!  LOCr( K ) denotes the number of elements of K that a process
!  would receive if K were distributed over the p processes of its
!  process column.
!  Similarly, LOCc( K ) denotes the number of elements of K that a
!  process would receive if K were distributed over the q processes of
!  its process row.
!  The values of LOCr() and LOCc() may be determined via a call to the
!  ScaLAPACK tool function, NUMROC:
!          LOCr( M ) = NUMROC( M, MB_A, MYROW, RSRC_A, NPROW ),
!          LOCc( N ) = NUMROC( N, NB_A, MYCOL, CSRC_A, NPCOL ).
!  An upper bound for these quantities may be computed by:
!          LOCr( M ) <= ceil( ceil(M/MB_A)/NPROW )!MB_A
!          LOCc( N ) <= ceil( ceil(N/NB_A)/NPCOL )!NB_A
!
!  Arguments
!  =========
!
!  SCOPE   (global input) CHARACTER*1
!          The BLACS scope in which alpha is updated.
!          If SCOPE = 'R', alpha is updated only in the process row
!                          containing A( IA, JA ),
!          If SCOPE = 'C', alpha is updated only in the process column
!                          containing A( IA, JA ),
!          If SCOPE = 'A', alpha is updated in all the processes of the
!                          grid,
!          otherwise alpha is updated only in the process containing
!           A( IA, JA ).
!
!  TOP     (global input) CHARACTER!1
!          The topology to be used if broadcast is needed.
!
!  ALPHA   (global output) DOUBLE PRECISION, the scalar alpha.
!
!  A       (local input) REAL(dp) pointer into the local memory
!          to an array of dimension (LLD_A,!) containing the local
!          pieces of the distributed matrix A.
!
!  IA      (global input) INTEGER
!          The row index in the global array A indicating the first
!          row of sub( A ).
!
!  JA      (global input) INTEGER
!          The column index in the global array A indicating the
!          first column of sub( A ).
!
!  DESCA   (global and local input) INTEGER array of dimension DLEN_.
!          The array descriptor for the distributed matrix A.
!
!  =====================================================================
!
!     .. Parameters ..
!      INTEGER            BLOCK_CYCLIC_2D, CSRC_, CTXT_, DLEN_, DTYPE_,
!     $                   LLD_, MB_, M_, NB_, N_, RSRC_
!      PARAMETER          ( BLOCK_CYCLIC_2D = 1, DLEN_ = 9, DTYPE_ = 1,
!     $                     CTXT_ = 2, M_ = 3, N_ = 4, MB_ = 5, NB_ = 6,
!     $                     RSRC_ = 7, CSRC_ = 8, LLD_ = 9 )
      REAL(dp)   ZERO
      PARAMETER          ( ZERO = 0.0D+0 )
!     ..
!     .. Local Scalars ..
      INTEGER            IACOL, IAROW, ICTXT, IIA, IOFFA, JJA  !, MYCOL, MYROW, NPCOL, NPROW
!     ..
!     .. External Subroutines ..
      EXTERNAL           BLACS_GRIDINFO, DGEBR2D, DGEBS2D, INFOG2L
!     ..
!     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
!     ..
!     .. Executable Statements ..
!
!     Get grid parameters.
!
      ICTXT = DESCA( CTXT_ )
      CALL BLACS_GRIDINFO( ICTXT, NPROW, NPCOL, MYROW, MYCOL )
!
      CALL INFOG2L( IA, JA, DESCA, NPROW, NPCOL, MYROW, MYCOL, IIA, JJA,&
                    IAROW, IACOL )
!
      ALPHA = ZERO
!
      IF( LSAME( SCOPE, 'R' ) ) THEN
         IF( MYROW.EQ.IAROW ) THEN
            IF( MYCOL.EQ.IACOL ) THEN
               IOFFA = IIA+(JJA-1)*DESCA( LLD_ )
               CALL DGEBS2D( ICTXT, SCOPE, TOP, 1, 1, A( IOFFA ), 1 )
               ALPHA = A( IOFFA )
            ELSE
               CALL DGEBR2D( ICTXT, SCOPE, TOP, 1, 1, ALPHA, 1,         &
                             IAROW, IACOL )
            END IF
         END IF
      ELSE IF( LSAME( SCOPE, 'C' ) ) THEN
         IF( MYCOL.EQ.IACOL ) THEN
            IF( MYROW.EQ.IAROW ) THEN
               IOFFA = IIA+(JJA-1)*DESCA( LLD_ )
               CALL DGEBS2D( ICTXT, SCOPE, TOP, 1, 1, A( IOFFA ), 1 )
               ALPHA = A( IOFFA )
            ELSE
               CALL DGEBR2D( ICTXT, SCOPE, TOP, 1, 1, ALPHA, 1,         &
                             IAROW, IACOL )
            END IF
         END IF
      ELSE IF( LSAME( SCOPE, 'A' ) ) THEN
         IF( ( MYROW.EQ.IAROW ).AND.( MYCOL.EQ.IACOL ) ) THEN
            IOFFA = IIA+(JJA-1)*DESCA( LLD_ )
            CALL DGEBS2D( ICTXT, SCOPE, TOP, 1, 1, A( IOFFA ), 1 )
            ALPHA = A( IOFFA )
         ELSE
            CALL DGEBR2D( ICTXT, SCOPE, TOP, 1, 1, ALPHA, 1,            &
                          IAROW, IACOL )
         END IF
      ELSE
         IF( MYROW.EQ.IAROW .AND. MYCOL.EQ.IACOL )                      &
            ALPHA = A( IIA+(JJA-1)*DESCA( LLD_ ) )
      END IF
!
      RETURN
!
!     End of PDELGET
!
      END
      
      SUBROUTINE PDELSET( A, IA, JA, DESCA, ALPHA )
      use descriptor_mod, desca_x=>desca
      implicit none
!
!  -- ScaLAPACK tools routine (version 1.7) --
!     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
!     and University of California, Berkeley.
!     May 1, 1997
!
!     .. Scalar Arguments ..
      INTEGER            IA, JA
      REAL(dp)   ALPHA
!     ..
!     .. Array arguments ..
      INTEGER            DESCA( * )
      REAL(dp)   A( * )
!     ..
!
!  Purpose
!  =======
!
!  PDELSET sets the distributed matrix entry A( IA, JA ) to ALPHA.
!
!  Notes
!  =====
!
!  Each global data object is described by an associated description
!  vector.  This vector stores the information required to establish
!  the mapping between an object element and its corresponding process
!  and memory location.
!
!  Let A be a generic term for any 2D block cyclicly distributed array.
!  Such a global array has an associated description vector DESCA.
!  In the following comments, the character _ should be read as
!  "of the global array".
!
!  NOTATION        STORED IN      EXPLANATION
!  --------------- -------------- --------------------------------------
!  DTYPE_A(global) DESCA( DTYPE_ )The descriptor type.  In this case,
!                                 DTYPE_A = 1.
!  CTXT_A (global) DESCA( CTXT_ ) The BLACS context handle, indicating
!                                 the BLACS process grid A is distribu-
!                                 ted over. The context itself is glo-
!                                 bal, but the handle (the integer
!                                 value) may vary.
!  M_A    (global) DESCA( M_ )    The number of rows in the global
!                                 array A.
!  N_A    (global) DESCA( N_ )    The number of columns in the global
!                                 array A.
!  MB_A   (global) DESCA( MB_ )   The blocking factor used to distribute
!                                 the rows of the array.
!  NB_A   (global) DESCA( NB_ )   The blocking factor used to distribute
!                                 the columns of the array.
!  RSRC_A (global) DESCA( RSRC_ ) The process row over which the first
!                                 row of the array A is distributed.
!  CSRC_A (global) DESCA( CSRC_ ) The process column over which the
!                                 first column of the array A is
!                                 distributed.
!  LLD_A  (local)  DESCA( LLD_ )  The leading dimension of the local
!                                 array.  LLD_A >= MAX(1,LOCr(M_A)).
!
!  Let K be the number of rows or columns of a distributed matrix,
!  and assume that its process grid has dimension p x q.
!  LOCr( K ) denotes the number of elements of K that a process
!  would receive if K were distributed over the p processes of its
!  process column.
!  Similarly, LOCc( K ) denotes the number of elements of K that a
!  process would receive if K were distributed over the q processes of
!  its process row.
!  The values of LOCr() and LOCc() may be determined via a call to the
!  ScaLAPACK tool function, NUMROC:
!          LOCr( M ) = NUMROC( M, MB_A, MYROW, RSRC_A, NPROW ),
!          LOCc( N ) = NUMROC( N, NB_A, MYCOL, CSRC_A, NPCOL ).
!  An upper bound for these quantities may be computed by:
!          LOCr( M ) <= ceil( ceil(M/MB_A)/NPROW )*MB_A
!          LOCc( N ) <= ceil( ceil(N/NB_A)/NPCOL )*NB_A
!
!  Arguments
!  =========
!
!  A       (local output) REAL(dp) pointer into the local memory
!          to an array of dimension (LLD_A,*) containing the local
!          pieces of the distributed matrix A.
!
!  IA      (global input) INTEGER
!          The row index in the global array A indicating the first
!          row of sub( A ).
!
!  JA      (global input) INTEGER
!          The column index in the global array A indicating the
!          first column of sub( A ).
!
!  DESCA   (global and local input) INTEGER array of dimension DLEN_.
!          The array descriptor for the distributed matrix A.
!
!  ALPHA   (local input) DOUBLE PRECISION
!          The scalar alpha.
!
!  =====================================================================
!
!     .. Parameters ..
!      INTEGER            BLOCK_CYCLIC_2D, CSRC_, CTXT_, DLEN_, DTYPE_,
!     $                   LLD_, MB_, M_, NB_, N_, RSRC_
!      PARAMETER          ( BLOCK_CYCLIC_2D = 1, DLEN_ = 9, DTYPE_ = 1,
!     $                     CTXT_ = 2, M_ = 3, N_ = 4, MB_ = 5, NB_ = 6,
!     $                     RSRC_ = 7, CSRC_ = 8, LLD_ = 9 )
!     ..
!     .. Local Scalars ..
      INTEGER            IACOL, IAROW, IIA, JJA   !, MYCOL, MYROW, NPCOL,  NPROW
!     ..
!     .. External Subroutines ..
      EXTERNAL           BLACS_GRIDINFO, INFOG2L
!     ..
!     .. Executable Statements ..
!
!     Get grid parameters.
!
      CALL BLACS_GRIDINFO( DESCA( CTXT_ ), NPROW, NPCOL, MYROW, MYCOL )
!
      CALL INFOG2L( IA, JA, DESCA, NPROW, NPCOL, MYROW, MYCOL, IIA, JJA,&
                    IAROW, IACOL )
!
      IF( MYROW.EQ.IAROW .AND. MYCOL.EQ.IACOL )                         &
         A( IIA+(JJA-1)*DESCA( LLD_ ) ) = ALPHA
!
      RETURN
!
!     End of PDELSET
!
      END
#endif
     
#endif
        
	  end module ptrd_mod
