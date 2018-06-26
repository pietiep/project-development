      subroutine mydiag(N,W,A)
*     .. Parameters ..
      implicit none
      INTEGER     N,i,j
      INTEGER     LWMAX
      PARAMETER   (LWMAX= 1000)
*     .. Local Scalars ..
      INTEGER           INFO, LWORK
*     .. Local Arrays ..

      complex *16  A(N,N), WORK(N*N),Temp(N,N)
      real*8 W(N),RWORK(3*N)

*
*     ..Executable Statements..
!      WRITE(*,*)'DSYEV Results'
      LWORK = N*N
*        Solve the symmetric eigenvalue problem
*


!        CALL zheev( 'Vectors', 'Upper', N, A, N, W, WORK, LWORK
!     $             ,RWORK,
!     $              INFO)

      call cDiag(A,N,N,N,W,Temp,N,N,.false.)

      A=Temp

*
         IF (INFO.GT.0) THEN
            WRITE(*,*)'Algorith failed'
            STOP
         END IF

c phasen convention for eigenvectors

!      do i=1,N
!          do j=1,N
!           if (real(A(j,i)).lt.0) then
!            A(j,i)=-A(j,i) 
!            endif
!           enddo
!           enddo      
*
*           Print solution
*
!      CALL r_PRINT_MATRIX('Eigenvalues', 1, N, W)
!      CALL PRINT_MATRIX('Eigenvectors', N, N, A)

      END


      SUBROUTINE PRINT_MATRIX( DESC, M, N, A )
      CHARACTER*(*)    DESC
      INTEGER          M, N
      complex*16 A(M, N )
*
      INTEGER          I, J
*
      WRITE(*,*)
      WRITE(*,*) DESC
      DO I = 1, M
         WRITE(*,9998) (real(A( I, J )), J = 1, N )
      END DO
*
 9998 FORMAT( 11(:,1X,F6.2) )
      RETURN
      end subroutine

      SUBROUTINE r_PRINT_MATRIX( DESC, M, N, A )
      CHARACTER*(*)    DESC
      INTEGER          M, N
      real*8 A(M, N )
*
      INTEGER          I, J
*
      WRITE(*,*)
      WRITE(*,*) DESC
      DO I = 1, M
         WRITE(*,9998) (A( I, J ), J = 1, N )
      END DO
*
 9998 FORMAT( 11(:,1X,F6.2) )
      RETURN
      end subroutine
      
      SUBROUTINE c_PRINT_MATRIX( DESC, M, N, A )
      CHARACTER*(*)    DESC
      INTEGER          M, N
      complex*16 A(M, N )
*
      INTEGER          I, J
*
      WRITE(*,*)
      WRITE(*,*) DESC
      DO I = 1, M
         WRITE(*,9998) (A( I, J ), J = 1, N )
      END DO
*
 9998 FORMAT( 11(:,1X,F6.2) )
      RETURN
      end subroutine

      SUBROUTINE PRINT_Vector_MATRIX( DESC, M, N, A,B )
      CHARACTER*(*)    DESC
      INTEGER          M, N, I ,J
      complex*16 A(M, N )
      real*8  B(M)
      complex*16 C(M,N+1)

      do I=1,M
       C(I,1)=B(I) 
      enddo

      do J=1,N
      do I=1,M
      C(I,J+1)=A(I,J)
      enddo
      enddo
      
*
*
      WRITE(*,*)
      WRITE(*,*) DESC
      DO I = 1, M 
         WRITE(*,9998) (real(C( I, J )), J = 1, N+1 )
      
      END DO
*
 9998 FORMAT( 11(:,1X,F6.2) )
      RETURN

      end subroutine


      SUBROUTINE compl_PRINT_MATRIX( DESC, M, N, A )
      CHARACTER*(*)    DESC
      INTEGER          M, N
      COMPLEX*16       A(M, N )
*
      INTEGER          I, J
*
      WRITE(*,*)
      WRITE(*,*) DESC
      DO I = 1, M
         WRITE(*,9998) ( A( I, J ), J = 1, N )
      END DO
*
 9998 FORMAT( 11(:,1X,'(',F6.2,',',F6.2,')') )
      RETURN
      END subroutine






