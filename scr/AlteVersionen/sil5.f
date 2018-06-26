      subroutine sil(v1,krylow_ev,krylow_ew,
     $ T_kin,V,trafo,m,N,dt,tphase,psi,xeigen)

      implicit none

      integer j,i,m,N
      complex*16 V(N),T_kin(N,N),trafo(N,N)
      complex*16 v1(N,m),norm,tphase,tphase_dt
      complex*16 krylow_ev(m,m),kry_dvr_ma(N,m),krylow_ev_dagger(m,m)
      real*8 krylow_ew(m),xeigen(N),dt
      complex*16 kry_FBR(N,m),kry_FBR_dagger(m,N)
      complex*16 phasen_vec(m),phasen_vecII(m), phasen_ew_vec(m)
      complex*16 psi_t(N),psi(N),psi_t_DVR(N),z
      integer steps
      real*8 total_t,pi
      parameter (pi=3.14159265359)

   
*****************************************
*     Time - Steps
*****************************************
      total_t=1
      steps=10
      dt=total_t/steps       
      
*****************************************
*      Startwellenfunktion v1
*****************************************
      do j=1,N
      do i=1,m 
      v1(j,i)=0
      enddo
      enddo
      
      do j=1,m
      v1(j,1)=psi(j)
      enddo
!      open(3,file='psi')
!      do i=1,N
!      read(3,*) v1(i,1)
!      print*, 'bla',psi_t(i)
!      enddo
      
      
!      call compl_ PRINT_MATRIX('psi_FBR in startwelle'
!     $ ,N,1,v1)
!      do j=1,m
!      v1(j,1)=1.d0
!      enddo
!      call vec_norm(v1,N,norm)
!      v1=v1/norm
!!      call PRINT_MATRIX('start',N,1,v1)

********************************************
*     First Iteration with Startwavefunction
********************************************



      call lanczos(v1,krylow_ev,krylow_ew,  T_kin,V,trafo,m,N
     $,xeigen )
!      call r_PRINT_MATRIX('ew im sil',1,m,krylow_ew)
!      call compl_PRINT_MATRIX('1 ev im sil',m,m,krylow_ev)
      call krytra(v1,krylow_ev,N,m,trafo,kry_dvr_ma,kry_FBR)
!      call compl_PRINT_MATRIX('2 ev im sil',m,m,krylow_ev)
     
      
!      print*,'dt',dt 
      
!      dt=-0.5
      do j=1,m
      phasen_ew_vec(j)=
     .      conjg(krylow_ev(1,j))*exp(dt*tphase*krylow_ew(j))
!      print*, dt,'dt'
!      do i=1,m
!      write(*,*) krylow_ev(i,j)
      enddo
!      enddo
      
!      call compl_PRINT_MATRIX('2 ev im sil',m,m,krylow_ev)
!      call r_PRINT_MATRIX('dt',1,1,dt)
!      call compl_PRINT_MATRIX('tphase',1,1,tphase)
!      call r_PRINT_MATRIX('ew',m,1,krylow_ew)
!      call compl_PRINT_MATRIX('exp(..)',m,1,exp(dt*tphase*krylow_ew))
!      
!      call compl_ PRINT_MATRIX('phasen_ew_vec',m,1,phasen_ew_vec)
!
!      call compl_PRINT_MATRIX('ev'  
!     $,m,m,krylow_ev)
      print*,''


      call matrixvec2(kry_FBR,phasen_ew_vec,N,m,psi_t) 
      
!      call compl_PRINT_MATRIX('Psi(t)',N,1,(psi_t))
      open(3,file='psi')
      do i=1,N
      write(3,*) psi_t(i)
      enddo
      close (3)
      
      
      call matrix_adj_vec(trafo,psi_t,N,psi_t_DVR)
      open(2,file='psi_DVR')
      write(2,*) ''
      write(2,*) '# 1'
      do i=1,N
      write(2,"(4e20.9)")xeigen(i), abs(psi_t_DVR(i)),
     $ real(psi_t_DVR(i)),aimag(psi_t_DVR(i))
      enddo

      
!      call compl_PRINT_MATRIX('Psi(t)_DVR',N,1,(psi_t_DVR))
!      call r_PRINT_MATRIX('abs(Psi(t))',N,1,abs(psi_t))
      print*,''
      print*,''
      print*,''
!      do j=1,N
!      print*,xeigen(j),psi_t_DVR(j)
!      enddo

*****************************************
*      Further time step Iterations
*****************************************
      do i=2,steps 
      
      v1=0

! Ueberschreibe v1 mit vorherigen Psi_t
      open(3,file='psi')
      do j=1,N
      read(3,*) v1(j,1)
      enddo
      close (3)

      write(*,*) i
!      call compl_ PRINT_MATRIX('psi_FBR in startwelle'
!     $ ,N,1,v1)
      
      krylow_ev=0
      krylow_ew=0 
     
      call lanczos(v1,krylow_ev,krylow_ew,  T_kin,V,trafo,m,N
     $,xeigen )
!      call r_PRINT_MATRIX('ew im sil',1,m,krylow_ew)
      call krytra(v1,krylow_ev,N,m,trafo,kry_dvr_ma,kry_FBR)
     
      
      do j=1,m
      phasen_ew_vec(j)=krylow_ev(1,j)*exp(dt*tphase*krylow_ew(j))
      enddo
     
!      call compl_ PRINT_MATRIX('phasen_ew_vec',m,1,phasen_ew_vec)

!      call compl_PRINT_MATRIX('ev'  
!     $,m,m,krylow_ev)
      print*,''


      call matrixvec2(kry_FBR,phasen_ew_vec,N,m,psi_t) 
      
!      call compl_PRINT_MATRIX('Psi(t)',N,1,(psi_t))
!      call r_PRINT_MATRIX('abs(Psi(t))',N,1,abs(psi_t))
     
      open(3,file='psi')
      do j=1,N
      write(3,*) psi_t(j)
      enddo
      close (3)

      call matrix_adj_vec(trafo,psi_t,N,psi_t_DVR)
      open(2,file='psi_DVR')
      write(2,*) ''
      write(2,*) '#',i
      do j=1,N
      write(2,"(4e20.9)")xeigen(j),abs(psi_t_DVR(j)), real(psi_t_DVR(i))
     $,aimag(psi_t_DVR(i)) 
      enddo
      
      print*,''
      print*,''
      print*,''
 
      enddo
      end subroutine

      
      
