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

      do j=1,m-1

      call matrix_adj_vec(trafo,v1(1,j),N,v1_DVR(1,j))
      
!     H|Psi> 
      call hamilton(T_kin,V,trafo,v1_DVR(1,j),N,wtilde_DVR)

!     H|Psi> in FBR
      call matrixvec(trafo,wtilde_DVR,N,wtilde)      

!    alpha = <Psi|H|Psi> in FBR
      call skalpro(wtilde,v1(1,j),N,alpha)


      !Koeffizient auf der Diagonalen

      d(j)=alpha

      !Nichtnormierter Vektor aus Gram-Schmidt
      call gramschmidt(wtilde,alpha,v1(1,j),beta,v0,N,w)


      call vec_norm(w,N,beta)

      !Norm entspricht Koeffzienten der Off-Diagonalen (Herl. S.320
      !Tannor 11.14)
      

!     Koeffizient der Offdiagonalen
      e(j+1)=beta

      call copyvec(v0,v1(1,j),N,1d0)
      call copyvec(v1(1,j+1),w,N,beta)

      enddo
      call matrix_adj_vec(trafo,v1(1,m),N,v1_DVR(1,m))
     
 
      call hamilton(T_kin,V,trafo,v1_DVR(1,m),N,wtilde_DVR)
!     H|Psi> in DVR


      call matrixvec(trafo,wtilde_DVR,N,wtilde)      
!     wtilde=H|Psi> in FBR
      
      call skalpro(wtilde,v1(1,m),N,alpha)
!     <Psi|H|Psi> = alpha
      
      d(j)=alpha


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
    
      
      call mydiag(m,krylow_ew,krylow_h,)

       end subroutine


       
     
