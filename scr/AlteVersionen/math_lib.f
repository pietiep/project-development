      subroutine vec_norm(psi,anzahlba,norm)
      implicit none
      
      integer i,j,anzahlba
      complex*16 psi(anzahlba),norm,dotpro
      
      norm=0
      
      do i=1,anzahlba
      dotpro = dotpro + (abs(psi(i))**2)
      enddo
      norm=sqrt(dotpro)
      end subroutine
     

 
      subroutine skalpro(a,b,anzahlba,dotpro)
      implicit none

      integer i,j
      integer anzahlba
      complex*16 a(anzahlba), b(anzahlba),dotpro

      dotpro=0
      do i=1,anzahlba
      dotpro=dotpro+conjg(a(i))*b(i)
      enddo

      end subroutine      
      

      subroutine matrixvec(a,b,anzahlba,prodmavec)
      implicit none
      integer i,j,anzahlba
      complex*16 a(anzahlba,anzahlba),b(anzahlba),prodmavec(anzahlba)
      

      do i=1,anzahlba
      prodmavec(i)=0
      do j=1,anzahlba
      prodmavec(i)=prodmavec(i)+a(i,j)*b(j)
      enddo
      enddo
      end subroutine

      subroutine matrixvec2(a,b,anzahlba,m,prodmavec)
      implicit none
      integer i,j,anzahlba,m
      complex*16 a(anzahlba,m),b(m),prodmavec(anzahlba)
      

      do i=1,anzahlba
      prodmavec(i)=0
      do j=1,m
      prodmavec(i)=prodmavec(i)+a(i,j)*b(j)
      enddo
      enddo
      end subroutine


      subroutine vecvec(a,b,anzahlba,vecneu)
      implicit none

      integer i,anzahlba
      complex*16 b(anzahlba),vecneu(anzahlba)
      complex*16 a(anzahlba)

      
      do i=1,anzahlba
      vecneu(i)=a(i)*b(i)
      enddo
      
      end subroutine
    
      subroutine vecadd(a,b,anzahlba,vecneu) 
      implicit none
      
      integer i,anzahlba
      complex*16 b(anzahlba),vecneu(anzahlba)
      complex*16 a(anzahlba)
      
      do i=1,anzahlba
      vecneu(i)=a(i)+b(i)
      enddo
      
      
      end subroutine      

      subroutine matrix_adj_vec(a,b,anzahlba,prodmavec)
      implicit none
      integer i,j,anzahlba
      complex*16 a(anzahlba,anzahlba),b(anzahlba),prodmavec(anzahlba)

      do i=1,anzahlba
      prodmavec(i)=0
      do j=1,anzahlba
      prodmavec(i)=prodmavec(i)+conjg(a(j,i))*b(j)
      enddo
      enddo
      end subroutine
      

      subroutine gramschmidt(wtilde,alpha,v1,beta,v0,m,w)
      implicit none

      integer j,m
      complex*16 v1(m),v0(m),beta
      complex*16 wtilde(m),w(m),alpha
      
      do j=1,m

      w(j)=wtilde(j)-alpha*v1(j)-beta*v0(j)

      enddo

      end subroutine

      subroutine copyvec(a,b,N,skalar)
      implicit none

      integer j,N
      complex*16 a(N),b(N),skalar

      do j=1,N
      a(j)=b(j)/skalar
      enddo

      end subroutine

      subroutine skal_vec(a,b,skal1,skal2,N)
      implicit none
      integer j,N

      complex*16 a(N),b(N),skal1,skal2
      
      do j=1,N
      a(j)=b(j)*skal1*skal2
      enddo
      
      end subroutine

      subroutine matr_mul(a,b,N,l,m,pro_ma)
      implicit none

      integer i,j,k,m,N,l
      complex*16 a(N,l),b(l,m),pro_ma(N,m)
      

      do i=1,m
      do k=1,N
      pro_ma(k,i)=0
      do j=1,l
!      pro_ma(k,i)=pro_ma(k,i)+conjg(a(k,j))*b(j,i) ! NEU (DEBUGGED)
      pro_ma(k,i)=pro_ma(k,i)+a(k,j)*b(j,i)
      enddo
      enddo
      enddo

      end subroutine 

      subroutine adjun(a,N,m,b)

      implicit none
      
      complex*16 a(N,m),b(m,N)
      integer i, j, N, m

      do i=1,N
      do j=1,m
      b(j,i)=conjg(a(i,j))
      enddo
      enddo

      end subroutine 

