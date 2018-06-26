c######################################################################
c Libary module l_diag
c
c Contains subroutines for diagonalizing matrices
c "cdiag", "rdiag", "generaldiag"
c and related sorting subroutines
c "maxsort", "minsort", "maxabssort", "cmaxabssort", "cminsort"
c
c V1.1, 03.11.2009
c######################################################################


c----------------------------------------------------------------------
c Libary subroutine cDiag
c
c Diagonalizes a hermitian complex matrix.
c
c Input-variables:  matrix  - complex matrix (only lower triangular part used)
c                   maxdim1 - declared dimension of first index of matrix
c                   maxdim2 - declared dimension of second index of matrix
c                   n       - dimension of matrix
c		    dimEV1  - declared dimension of first index of ev
c		    dimEV2  - declared dimension of second index of ev
c                   copy    - if .false., matrix is overwriten
c Output-variables: ew      - eigenvalues (ew(j): j-th eigenvalue)
c                   ev      - orthonormal eigenvectors (ev(i,j): i-th component of j-th eigenvalue/-vector)
c----------------------------------------------------------------------

      subroutine cDiag 
     .(matrix, maxdim1, maxdim2, n, ew, ev, dimEV1, dimEV2, copy)

      implicit none
      integer      maxdim1,maxdim2,dimEV1,dimEV2, n
      complex*16   matrix(maxdim1,maxdim2)

      logical      copy
c      complex*16   a_copy(n,n)
      real*8	    abstol,dlamch
      integer      m
      real*8	    ew(n)
      complex*16   ev(dimEV1,dimEV2)
      integer	    ISUPPZ(2*n)
c      complex*16   work(41*n)
      complex*16   cxLWork,cxPhase
c      real*8	    rwork(24*n)
      real*8	    rLRWork
      integer      LWork,LRWork,LIWork,info
      integer	    i,j,k,kk
      complex*16, allocatable::a_copy(:,:),Work(:)
      real*8, allocatable::RWork(:)
      integer, allocatable::IWork(:)
      integer      istat
      
      if(maxdim1<n .or. maxdim2<n) then
        write (6,*) "Dimension-Error in cdiag (matrix)"
        stop      
      endif
      if(dimEV1<n .or. dimEV2<n) then
        write (6,*) "Dimension-Error in cdiag (ev)"
        stop      
      endif
      abstol=dlamch('S')
      if (copy) then
        allocate (a_copy(n,n),STAT=istat)
        if(istat>0) then
          write (6,*) "Allocation-Error in cdiag (a_copy)"
          stop
        endif
c 	copy lower trangular part (incl. diagonal)
        do i=1,n
          do j=1,i
            a_copy(i,j)=matrix(i,j)
          enddo
        enddo
      endif 
c Driver Routines for Solving Symmetric Eigenproblems: 
c Default Routine: zheevr 
c call ZHEEVR( JOBZ, RANGE, UPLO, N, A, LDA, VL, VU, IL, IU, ABSTOL, M, W, Z, LDZ, ISUPPZ, WORK, LWORK, RWORK, LRWORK, IWORK, LIWORK, INFO )
c          CHARACTER      JOBZ, RANGE, UPLO
c          INTEGER        IL, INFO, IU, LDA, LDZ, LIWORK, LRWORK, LWORK, M, N
c          DOUBLE         PRECISION ABSTOL, VL, VU
c          INTEGER        ISUPPZ( * ), IWORK( * )
c          DOUBLE         PRECISION RWORK( * ), W( * )
c          COMPLEX*16     A( LDA, * ), WORK( * ), Z( LDZ, * )
c JOBZ here: 'V': eigenvaluse+eigenvectors
c RANGE here: 'A':all eigenvalues
c UPLO 'U': use upper triangular part of A, 'L': use lower triangular part of A.
c by testing: 'L'--> first component of Eigenvector Z(1,i) is real
c by testing: 'U'--> last component of Eigenvector Z(N,i) is real
c N: Order of matrix A (N>=0) (A in M(NxN))
c A: Matrix A(lda,*>=max(1,n)), upper/lower (UPLO) part of triang. matrix
c A is hermite --> A_ij=comp. conj. (A_ji)--> upper/lower part of triag. is enough
c LDA: leading dimension of A: >=max(1,n)
c VL,VU: Not referenced (RANGE='A')
c IL,IU: Not referenced (RANGE='A')
c ABSTOL:high accuracy: dlamch('S')
c M: out: #eigenvalues, 0<=M<=N (RANGE='A': M=N)
c W: out: W(>=max(1,M)); eigenvalues in ascending order
c Z: out: Z(LDZ,*=max(1,M)): eigenvectors of A: Z(*,i)<-->W(i)
c LDZ: Leading dimension of z: >=max(1,n) if JOBZ='V'
c ISUPPZ: out: ISUPPZ(2*max(1,M))
c WORK:  complex*16 WORK(>=max(1,lwork))
c LWORK: dim of WORK: >=max(1,2n) (see app. note) (n<1000: 33*n; n>=1000: 41*n, tested on tc02,tc248,tc401)
c RWORK: real*8 RWORK(>=max(1,LRWORK))
c LRWORK: dim of RWORK: >=max(1,24n)
c IWORK: integer IWORK(max(1,LIWORK))
c LIWORK: >=max(1,10n)
c INFO: 0=OK, <0:INFO=-i, illegal value in i-th argument, >0: internal error
C
c     calculate optimal values for LWork,LRWork and LIWork
      call zHeEvr('V','A','L',n,matrix,maxdim1,0,0,0,0,abstol,
     .m,ew,ev,dimEV1,ISUPPZ,cxLWork,-1,rLRWork,-1,LIWork,-1,info)
      LWork=cxLWork
      LRWork=rLRWork
      allocate (RWork(LRWork),IWork(LIWork),STAT=istat)
      if(istat>0) then
        write (6,*) "Allocation-Error in cdiag (RWork,IWork)"
        stop
      endif
      allocate (Work(LWork),STAT=istat)
      if(istat>0) then
c	not enought memory??, try lowest allowed size
	LWork=2*n
	allocate (Work(LWork),STAT=istat)
      endif
      if(istat>0) then
        write (6,*) "Allocation-Error in cdiag (Work)"
        stop
      endif
      if (copy) then
        call zHeEvr('V','A','L',n,a_copy,n,0,0,0,0,abstol,
     .m,ew,ev,dimEV1,ISUPPZ,Work,LWork,RWork,LRWork,IWork,LIWork,info)
        deallocate (a_copy)
      else
        call zHeEvr('V','A','L',n,matrix,maxdim1,0,0,0,0,abstol,
     .m,ew,ev,dimEV1,ISUPPZ,Work,LWork,RWork,LRWork,IWork,LIWork,info)
      endif
      deallocate (Work,RWork,IWork)
      if (info.ne.0) then
        write (6,*) "Error ",info,
     ." in LaPack(zHeEvr) (diagonalising matrix)"
        stop
      endif
c phase convention for eigenvectors
      do k=1,m
        if (abs(dimag(ev(1,k))).le.abstol) then
          if (dble(ev(1,k)).lt.0) then
            do kk=1,n
              ev(kk,k)=-ev(kk,k)
            enddo
          endif
        else
c this part should never be used, except:
c internal calculation of LaPack routine changed, first component has an imag. part
c calculate phase factor cxPhase for (x,y)-->(r,0)
c (x,y)*cxPhase=(r,0) with r=abs((x,y)) (only rotate, don't change vector length)
c ->cxPhase=(r,0)/(x,y)
          cxPhase=abs(ev(1,k))/ev(1,k)   
          do kk=1,n
            ev(kk,k)=ev(kk,k)*cxPhase
          enddo
        endif
      enddo
      end 
c     end cDiag



c----------------------------------------------------------------------
c Libary subroutine rDiag
c
c Diagonalizes a real symmetric matrix.
c
c Input-variables:  matrix  - real symmetric matrix (only lower triangular part used)
c                   maxdim1 - declared dimension of first index of matrix
c                   maxdim2 - declared dimension of second index of matrix
c                   n       - dimension of matrix
c		    dimEV1  - declared dimension of first index of ev
c		    dimEV2  - declared dimension of second index of ev
c                   copy    - if .false., matrix is overwriten
c Output-variables: ew      - eigenvalues (ew(j): j-th eigenvalue)
c                   ev      - orthonormal eigenvectors (ev(i,j): i-th component of j-th eigenvalue/-vector)
c----------------------------------------------------------------------

      subroutine rDiag 
     .(matrix, maxdim1, maxdim2, n, ew, ev, dimEV1, dimEV2, copy,llo)

      implicit none
      integer      maxdim1,maxdim2,dimEV1,dimEV2, n
      real*8  	    matrix(maxdim1,maxdim2)
      integer llo ! defines convention

      logical      copy
      real*8	    abstol,dlamch
      integer      m
      real*8	    ew(n)
      real*8 	    ev(dimEV1,dimEV2)
      integer	    ISUPPZ(2*n)
      real*8	    rLWork
      integer      LWork,LIWork,info
      integer	    i,j,k,kk
      real*8, allocatable::a_copy(:,:),Work(:)
      integer, allocatable::IWork(:)
      integer      istat
      
      if(maxdim1<n .or. maxdim2<n) then
        write (6,*) "Dimension-Error in rdiag (matrix)"
        stop      
      endif
      if(dimEV1<n .or. dimEV2<n) then
        write (6,*) "Dimension-Error in rdiag (ev)"
        stop      
      endif
      abstol=dlamch('S')
c      abstol=0
      if (copy) then
        allocate (a_copy(n,n),STAT=istat)
        if(istat>0) then
          write (6,*) "Allocation-Error in rdiag (a_copy)"
          stop
        endif
c 	copy lower trangular part (incl. diagonal)
        do i=1,n
          do j=1,i
            a_copy(i,j)=matrix(i,j)
          enddo
        enddo
      endif 
c call DSYEVR( JOBZ, RANGE, UPLO, N, A, LDA, VL, VU, IL, IU, ABSTOL, M, W, Z, LDZ, ISUPPZ, WORK, LWORK, IWORK, LIWORK, INFO )
c           CHARACTER      JOBZ, RANGE, UPLO
c           INTEGER        IL, INFO, IU, LDA, LDZ, LIWORK, LWORK, M, N
c           DOUBLE         PRECISION ABSTOL, VL, VU
c           INTEGER        ISUPPZ( * ), IWORK( * )
c           DOUBLE         PRECISION A( LDA, * ), W( * ), WORK( * ), Z( LDZ, * )
c JOBZ,RANGE,UPLO,N,A,LDA,VL,VU,IL,IU,ABSTOL,M,W,Z,LDZ,ISUPPZ,IWORK,LIWORK,INFO see above (zHeEvr)
c WORK: real*8 WORK(>=max(1,LWORK))
c LWORK: dim of WORK: >=max(1,26n) (n<1100: 26*n; n>=1100: 33*n, tested on tc02,tc401,tc248)
C LWORK=-1 or LIWORK=-1: calc. optimal values for LWORK, LIWORK, stored in WORK,IWORK
c     calculate optimal values for LWork and LIWork
      call dSyEvr('V','A','L',n,matrix,maxdim1,0,0,0,0,abstol,
     .m,ew,ev,dimEV1,ISUPPZ,rLWork,-1,LIWork,-1,info)
      LWork=rLWork
      allocate (IWork(LIWork),STAT=istat)
      if(istat>0) then
        write (6,*) "Allocation-Error in rdiag (IWork)"
        stop
      endif
      allocate (Work(LWork),STAT=istat)
      if(istat>0) then
c	not enought memory??, try lowest allowed size
	LWork=26*n
	allocate (Work(LWork),STAT=istat)
      endif
      if(istat>0) then
        write (6,*) "Allocation-Error in rdiag (Work)"
        stop
      endif
      if (copy) then
        call dSyEvr('V','A','L',n,a_copy,n,0,0,0,0,abstol,
     .m,ew,ev,dimEV1,ISUPPZ,Work,LWork,IWork,LIWork,info)
        deallocate (a_copy)
      else
        call dSyEvr('V','A','L',n,matrix,maxdim1,0,0,0,0,abstol,
     .m,ew,ev,dimEV1,ISUPPZ,Work,LWork,IWork,LIWork,info)
      endif
      deallocate (Work,IWork)
      if (info.ne.0) then
         write (6,*) "Error ",info,
     ." in LaPack(dSyEvr) (diagonalising matrix)"
         stop
      endif
c phase convention for eigenvectors
      if (llo.eq.1) then
      do k=1,m
         if (ev(1,k).lt.0) then
            do kk=1,n
               ev(kk,k)=-ev(kk,k)
            enddo
         endif
      enddo
      endif
      end 
c     end rDiag


c----------------------------------------------------------------------
c Libary subroutine generalDiag
c
c Diagonalizes a complex matrix.
c
c Input-variables:  matrix  - complex matrix
c                   maxdim1 - declared dimension of first index of matrix
c                   maxdim2 - declared dimension of second index of matrix
c                   n       - dimension of matrix
c		    dimEV1  - declared dimension of first index of ev
c		    dimEV2  - declared dimension of second index of ev
c                   copy    - if .false., matrix is overwriten
c Output-variables: ew      - eigenvalues (ew(j): j-th eigenvalue)
c                   ev      - orthonormal eigenvectors (ev(i,j): i-th component of j-th eigenvalue/-vector)
c----------------------------------------------------------------------

      subroutine generalDiag
     .(matrix, maxdim1, maxdim2, n, ew, ev, dimEV1, dimEV2, copy)

      implicit none
      integer      maxdim1,maxdim2,dimEV1,dimEV2, n
      complex*16   matrix(maxdim1,maxdim2)

      logical      copy
c      complex*16   a_copy(n,n)
      real*8	    abstol,dlamch
      integer       m
      complex*16    ew(n)
      complex*16    ev(dimEV1,dimEV2)
      integer	    ISUPPZ(2*n)
      complex*16    cxLWork,cxPhase
      real*8	    RWork(2*n)
      real*8	    rLRWork
      integer       LWork,info
      integer	    i,j,k,kk
      complex*16, allocatable::a_copy(:,:),Work(:)
      integer      istat
      
      if(maxdim1<n .or. maxdim2<n) then
        write (6,*) "Dimension-Error in generalDiag (matrix)"
        stop      
      endif
      if(dimEV1<n .or. dimEV2<n) then
        write (6,*) "Dimension-Error in generalDiag (ev)"
        stop      
      endif
      abstol=dlamch('S')
      if (copy) then
        allocate (a_copy(n,n),STAT=istat)
        if(istat>0) then
          write (6,*) "Allocation-Error in generalDiag (a_copy)"
          stop
        endif
c 	copy full matrix
        do i=1,n
          do j=1,n
            a_copy(i,j)=matrix(i,j)
          enddo
        enddo
      endif 

c     calculate optimal values for LWork,LRWork and LIWork
      call zGeEv('N','V',n,matrix,maxdim1,ew,ev,dimEV1,ev,dimEV1,
     .     cxLWork,-1,RWork,info)
      LWork=cxLWork
      allocate (Work(LWork),STAT=istat)
      if(istat>0) then
c	not enought memory??, try lowest allowed size
         LWork=2*n
         allocate (Work(LWork),STAT=istat)
      endif
      if(istat>0) then
        write (6,*) "Allocation-Error in generalDiag (Work)"
        stop
      endif
      if (copy) then
         call zGeEv('N','V',n,a_copy,n,ew,ev,dimEV1,ev,dimEV1,
     .        Work,LWork,RWork,info)
         deallocate (a_copy)
      else
         call zGeEv('N','V',n,matrix,maxdim1,ew,ev,dimEV1,ev,dimEV1,
     .        Work,LWork,RWork,info)
      endif
      deallocate(Work)
      if (info.ne.0) then
        write (6,*) "Error ",info,
     ." in LaPack(zGeEv) (diagonalising matrix)"
        stop
      endif
c phase convention for eigenvectors
      do k=1,m
        if (abs(dimag(ev(1,k))).le.abstol) then
          if (dble(ev(1,k)).lt.0) then
            do kk=1,n
              ev(kk,k)=-ev(kk,k)
            enddo
          endif
        else
c this part should never be used, except:
c internal calculation of LaPack routine changed, first component has an imag. part
c calculate phase factor cxPhase for (x,y)-->(r,0)
c (x,y)*cxPhase=(r,0) with r=abs((x,y)) (only rotate, do not change vector length)
c ->cxPhase=(r,0)/(x,y)
          cxPhase=abs(ev(1,k))/ev(1,k)   
          do kk=1,n
            ev(kk,k)=ev(kk,k)*cxPhase
          enddo
        endif
      enddo
      end 
c     end generalDiag

c--------------------------------------------------------------------
c Libary subroutine maxsort
c
c Resorts vectors "vec" to decreasing order of corresponding value "w"
c
c Input-variables:  w     - values
c                   vec   - vectors
c                   dim   - dimension of vector space
c                   zahl  - number of vectors
c Output-variables: w     - sorted values
c                   vec   - sorted eigenvectors
c----------------------------------------------------------------------

c     subroutine maxsort (w, vec, dim, zahl)
c     
c     implicit none 
c     integer     dim, zahl, i, j, n
c     real*8      w(*), x
c     complex*16  vec(dim,zahl), z

c     do i=1,zahl
c        do j=i+1,zahl
c           if (w(j).gt.w(i)) then
c              x=w(i)
c              w(i)=w(j)
c              w(j)=x
c              do n=1,dim
c                 z=vec(n,i)
c                 vec(n,i)=vec(n,j)
c                 vec(n,j)=z
c              enddo
c           endif
c        enddo
c     enddo
c     end         

c--------------------------------------------------------------------
c Libary subroutine minsort
c
c Resorts vectors "vec" to increasing order of corresponding value "w"
c
c Input-variables:  w     - values
c                   vec   - vectors
c                   dim   - dimension of vector space
c                   zahl  - number of vectors
c Output-variables: w     - sorted values
c                   vec   - sorted eigenvectors
c----------------------------------------------------------------------

c     subroutine minsort (w, vec, dim, zahl)
c     
c     implicit none 
c     integer     dim, zahl, i, j, n
c     real*8      w(*), x
c     complex*16  vec(dim,zahl), z

c     do i=1,zahl
c        do j=i+1,zahl
c           if (w(j).lt.w(i)) then
c              x=w(i)
c              w(i)=w(j)
c              w(j)=x
c              do n=1,dim
c                 z=vec(n,i)
c                 vec(n,i)=vec(n,j)
c                 vec(n,j)=z
c              enddo
c           endif
c        enddo
c     enddo
c     end         

c--------------------------------------------------------------------
c Libary subroutine maxabssort
c
c Resorts vectors "vec" to decreasing order of corresponding absolute
c value of "w"
c
c Input-variables:  w     - values
c                   vec   - vectors
c                   dim   - dimension of vector space
c                   zahl  - number of vectors
c Output-variables: w     - sorted values
c                   vec   - sorted eigenvectors
c----------------------------------------------------------------------

c     subroutine maxabssort (w, vec, dim, zahl)
c     
c     implicit none 
c     integer     dim, zahl, i, j, n
c     real*8      w(*), x
c     complex*16  vec(dim,zahl), z

c     do i=1,zahl
c        do j=i+1,zahl
c           if (abs(w(j)).gt.abs(w(i))) then
c              x=w(i)
c              w(i)=w(j)
c              w(j)=x
c              do n=1,dim
c                 z=vec(n,i)
c                 vec(n,i)=vec(n,j)
c                 vec(n,j)=z
c              enddo
c           endif
c        enddo
c     enddo
c     end         


c--------------------------------------------------------------------
c Libary subroutine cmaxabssort
c
c Resorts vectors "vec" to decreasing order of corresponding absolute
c value of "w"
c
c Input-variables:  w     - values
c                   vec   - vectors
c                   dim   - dimension of vector space
c                   zahl  - number of vectors
c Output-variables: w     - sorted values
c                   vec   - sorted eigenvectors
c----------------------------------------------------------------------

      subroutine cmaxabssort (w, vec, dim, zahl)
      
      implicit none 
      integer     dim, zahl, i, j, n
      complex*16  vec(dim,zahl), w(*), z

      do i=1,zahl
         do j=i+1,zahl
            if (abs(w(j)).gt.abs(w(i))) then
               z=w(i)
               w(i)=w(j)
               w(j)=z
               do n=1,dim
                  z=vec(n,i)
                  vec(n,i)=vec(n,j)
                  vec(n,j)=z
               enddo
            endif
         enddo
      enddo
      end         

c--------------------------------------------------------------------
c Libary subroutine cminsort
c
c Resorts vectors "vec" to increasing order of corresponding real
c value of "w"
c
c Input-variables:  w     - values
c                   vec   - vectors
c                   dim   - dimension of vector space
c                   zahl  - number of vectors
c Output-variables: w     - sorted values
c                   vec   - sorted eigenvectors
c----------------------------------------------------------------------

      subroutine cminsort (w, vec, dim, zahl)
      
      implicit none 
      integer     dim, zahl, i, j, n
      complex*16  vec(dim,zahl), w(*), z

      do i=1,zahl
         do j=i+1,zahl
            if (dreal(w(j)).lt.dreal(w(i))) then
               z=w(i)
               w(i)=w(j)
               w(j)=z
               do n=1,dim
                  z=vec(n,i)
                  vec(n,i)=vec(n,j)
                  vec(n,j)=z
               enddo
            endif
         enddo
      enddo
      end         



