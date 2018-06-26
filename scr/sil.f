      subroutine sil(name,v1,krylow_ev,krylow_ew,
     $ T_kin,V,trafo,m,N,dt,
     $ tphase,total_t,psi,xeigen)

      implicit none

      character*20 name
      integer j,i,m,N
      complex*16 V(N),T_kin(N,N),trafo(N,N)
      complex*16 v1(N,m),norm,tphase,tphase_dt
      complex*16 krylow_ev(m,m),kry_dvr_ma(N,m),krylow_ev_dagger(m,m)
      real*8 krylow_ew(m),xeigen(N),dt
      complex*16 kry_FBR(N,m),kry_FBR_dagger(m,N)
      complex*16 phasen_vec(m),phasen_vecII(m), phasen_ew_vec(m)
      complex*16 psi_t(N),psi(N),psi_DVR(N),z
      integer steps
      real*8 total_t,pi
      parameter (pi=3.14159265359)

*****************************************
*     Time - Steps
*****************************************
      steps=int(total_t/dt)
      
*****************************************
*      Startwellenfunktion v1
*****************************************
      do j=1,N
      do i=1,m 
      v1(j,i)=0
      enddo
      enddo
      
*****************************************
*      Further time step Iterations
*****************************************
      do i=1,steps 
      
      v1=0
! Ueberschreibe v1 mit vorherigen Psi_t
      do j=1,N
      v1(j,1)=psi(j)
      enddo
      
      krylow_ev=0
      krylow_ew=0 
     
      call lanczos(v1,krylow_ev,krylow_ew,  T_kin,V,trafo,m,N
     $,xeigen )
      call krytra(v1,krylow_ev,N,m,trafo,kry_dvr_ma,kry_FBR)
      
      do j=1,m
      phasen_ew_vec(j)=krylow_ev(1,j)*exp(dt*tphase*krylow_ew(j))
      enddo

      call matrixvec2(kry_FBR,phasen_ew_vec,N,m,psi) 
      call matrix_adj_vec(trafo,psi,N,psi_DVR)
      call writepsi(psi_DVR,N)
 
      enddo

      end subroutine
