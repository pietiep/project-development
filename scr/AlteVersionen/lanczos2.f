      subroutine lanczos(T_kin,V,trafo,m,N)
      implicit none

      integer j,i,m,N
      complex*16 V(N),T_kin(N,N),trafo(N,N)
      complex*16 v1(N),v1_DVR(N),v0(N),beta,vv1(m,N)
      complex*16 wtilde(N),w(N),alpha, wtilde_DVR(N)
      complex*16 d(m),e(m),krylow_h(m,m)
      real*8 krylow_ew(m)

      v0=0
      beta=0
*****************************************
*      Startwellenfunktion v1
*****************************************
    
      v1=0
      v1(1)=1/(sqrt(2d0))
      v1(3)=1/(sqrt(2d0))
      call PRINT_MATRIX('start FBR',1,N,v1)

!      call matrix_adj_vec(trafo,v1,N,v0)
!      do i=1,N
!      v1(i)=v0(i)
!      v0(i)=0
!      enddo
!      
!      call PRINT_MATRIX('start DVR',1,N,v1)
!
!      call matrixvec(trafo,v1,N,v0)
!
!      call PRINT_MATRIX('start FBR',1,N,v0)
!
!      call hamilton(T_kin,V,trafo,v1,N,wtilde)
!
!      call matrixvec(trafo,wtilde,N,v0)
!
!      call PRINT_MATRIX('H|Psi> in FBR',1,N,v0)
!      Stop
******************************************
*     Lanczos
******************************************
!      call PRINT_MATRIX('ekin',N,N,T_kin)

      do j=1,m-1

      call matrix_adj_vec(trafo,v1,N,v1_DVR)
      call PRINT_MATRIX('start DVR',1,N,v1_DVR)
      
      call hamilton(T_kin,V,trafo,v1_DVR,N,wtilde_DVR)
      call PRINT_MATRIX('H|Psi> in DVR',1,N,wtilde_DVR)

      call matrixvec(trafo,wtilde_DVR,N,wtilde)      
      call PRINT_MATRIX('H|Psi> in FBR',1,N,wtilde)

      call skalpro(wtilde,v1,N,alpha)

      call PRINT_MATRIX('<Psi|H|Psi> in FBR',1,1,alpha)

      !Koeffizient auf der Diagonalen

      d(j)=alpha

      call gramschmidt(wtilde,alpha,v1,beta,v0,N,w)

      !Nichtnormierter Vektor aus Gram-Schmidt
      call PRINT_MATRIX('w in FBR',1,N,w)

      call vec_norm(w,N,beta)

      !Norm entspricht Koeffzienten der Off-Diagonalen (Herl. S.320
      !Tannor 11.14)
      

      call PRINT_MATRIX('beta in FBR',1,1,beta)

      e(j+1)=beta

      v0=v1

      v1=w/beta

      enddo

      call PRINT_MATRIX('v1 FBR',1,N,v1)
      call matrix_adj_vec(trafo,v1,N,v1_DVR)
      call PRINT_MATRIX('v1 DVR',1,N,v1_DVR)
            
      call hamilton(T_kin,V,trafo,v1_DVR,N,wtilde_DVR)
      call PRINT_MATRIX('H|Psi> in DVR',1,N,wtilde_DVR)

      call matrixvec(trafo,wtilde_DVR,N,wtilde)      
      call PRINT_MATRIX('H|Psi> in FBR',1,N,wtilde)
     
      call skalpro(wtilde,v1,N,alpha)
      call PRINT_MATRIX('<Psi|H|Psi> in FBR',1,1,alpha)

      d(j)=alpha

      call gramschmidt(wtilde,alpha,v1,beta,v0,N,w)

      call vec_norm(w,N,beta)
      call PRINT_MATRIX('beta in FBR',1,1,beta)

      v0=v1
      
      v1=w/beta

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
    

      
      call PRINT_MATRIX("krylow_h",m,m,krylow_h)
      
      call mydiag(m,krylow_ew,krylow_h,)
      do i=1,m
      print*, krylow_ew(i)
      enddo
      

      end subroutine


       
     
