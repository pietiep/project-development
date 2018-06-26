      subroutine lanczos(v1,krylow_h,krylow_ew,
     $ T_kin,V,trafo,m,N,xeigen)
      implicit none

      integer j,i,m,N
      complex*16 V(N),T_kin(N,N),trafo(N,N)
      complex*16 v1(N,m),v1_DVR(N,m),v0(N),beta
      complex*16 wtilde(N),w(N),alpha, wtilde_DVR(N)
      complex*16 d(m),e(m),krylow_h(m,m),norm
      real*8 krylow_ew(m),xeigen(N)

      v0=0
      beta=0

******************************************
*     Lanczos
******************************************
!      call PRINT_MATRIX('ekin',N,N,T_kin)

      do j=1,m-1

      call matrix_adj_vec(trafo,v1(1,j),N,v1_DVR(1,j))
!      call PRINT_MATRIX('start DVR',1,N,v1_DVR(1,j))
!      call PRINT_Vector_MATRIX('plot start',N,1,v1_DVR,xeigen)
      
      call hamilton(T_kin,V,trafo,v1_DVR(1,j),N,wtilde_DVR)
!      call PRINT_MATRIX('H|Psi> in DVR',1,N,wtilde_DVR)

      call matrixvec(trafo,wtilde_DVR,N,wtilde)      
!      call PRINT_MATRIX('H|Psi> in FBR',1,N,wtilde)

      call skalpro(wtilde,v1(1,j),N,alpha)

!      call PRINT_MATRIX('<Psi|H|Psi> in FBR',1,1,alpha)

      !Koeffizient auf der Diagonalen

      d(j)=alpha

      call gramschmidt(wtilde,alpha,v1(1,j),beta,v0,N,w)

      !Nichtnormierter Vektor aus Gram-Schmidt
!      call PRINT_MATRIX('w in FBR',1,N,w)

      call vec_norm(w,N,beta)

      !Norm entspricht Koeffzienten der Off-Diagonalen (Herl. S.320
      !Tannor 11.14)
      

!      call PRINT_MATRIX('beta in FBR',1,1,beta)
      e(j+1)=beta

!      v0(1,j)=v1(1,j)
      call copyvec(v0,v1(1,j),N,1d0)
!      call PRINT_MATRIX('v0',1,N,v0)
!      v1(1,j)=w(1,j)/beta
      call copyvec(v1(1,j+1),w,N,beta)
!      call PRINT_MATRIX('v1',1,N,v1)

      enddo
!      call PRINT_MATRIX('v1 FBR',1,N,v1(1,m))
      call matrix_adj_vec(trafo,v1(1,m),N,v1_DVR(1,m))
!      call PRINT_MATRIX('v1 DVR',1,N,v1_DVR(1,m))
     
!      call PRINT_MATRIX('V(x-r0)',1,N,V)
!      call r_PRINT_MATRIX('xeigen',1,N,xeigen)
!      call PRINT_MATRIX('Tkin',1,N,T_kin)
!      call PRINT_MATRIX('trafo',1,N,trafo)
!      call PRINT_MATRIX('trans_trafo',1,N,Transpose(trafo))


 
      call hamilton(T_kin,V,trafo,v1_DVR(1,m),N,wtilde_DVR)
!      call PRINT_MATRIX('H|Psi> in DVR',1,N,wtilde_DVR)

!      call skalpro(wtilde_DVR,v1_DVR,N,alpha)

      call matrixvec(trafo,wtilde_DVR,N,wtilde)      
!      call PRINT_MATRIX('H|Psi> in FBR',1,N,wtilde)
      
      call skalpro(wtilde,v1(1,m),N,alpha)
!      call PRINT_MATRIX('<Psi|H|Psi> in FBR',1,1,alpha)
      
      d(j)=alpha

!     call PRINT_MATRIX('Krylowvektor in FBR',N,m,v1)

!      call gramschmidt(wtilde(1,m),alpha,v1(1,m),beta,v0(1,m),N,w(1,m))
!
!      call vec_norm(w(1,m),N,beta)
!      call PRINT_MATRIX('beta in FBR',1,1,beta)
!
!!      v0(1,m)=v1(1,m)
!      call copyvec(v0(1,m),v1(1,m),N,1d0)
!      call PRINT_MATRIX('v0',1,N,v0(1,m))
!      call copyvec(v1(1,j),w(1,j),N,beta)
!      call PRINT_MATRIX('v1',1,N,v1(1,m))
!      v1(1,m)=w(1,m)/beta

      do i=1,m
      do j=1,m
      if(j.eq.i) then
      krylow_h(j,i)=d(j)
      else
      krylow_h(j,i)=0
      endif
      enddo
      enddo


      do i=1,m
      do j=1,m
      if (i.eq.(j+1)) then
      krylow_h(j,i)=e(j+1) !da, e(1)=0
      else if (j.eq.(i+1)) then
      krylow_h(j,i)=e(j)
      endif
      enddo
      enddo
    

      
!      call PRINT_MATRIX("Hn",m,m,krylow_h)
      
      call mydiag(m,krylow_ew,krylow_h,)


!      call r_PRINT_MATRIX('krylow_ew im lanczos',1,m,krylow_ew)      
!      call PRINT_MATRIX('Eigenzustaende in Krylowdarstellung'
!     $ ,m,m,krylow_h)
     

       end subroutine


       
     
