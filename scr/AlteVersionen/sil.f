      subroutine sil(v1,krylow_ev,krylow_ew,
     $ T_kin,V,trafo,m,N,dt,tphase,psi)

      implicit none

      integer j,i,m,N
      complex*16 V(N),T_kin(N,N),trafo(N,N)
      complex*16 v1(N,m),norm,tphase,tphase_dt
      complex*16 krylow_ev(m,m),kry_dvr_ma(N,m),krylow_ev_dagger(m,m)
      real*8 krylow_ew(m),xeigen(N),dt
      complex*16 kry_FBR(N,m),kry_FBR_dagger(m,N)
      complex*16 phasen_vec(m),phasen_vecII(m), phasen_ew_vec(m)
      complex*16 psi_t(N),psi(N)
      integer steps
      real*8 total_t

*****************************************
*      Startwellenfunktion v1
*****************************************
      v1=0
      
!      do j=1,m
!      v1(j,1)=psi(j)
!      enddo
      
      do j=1,m
      v1(j,1)=1.d0
      enddo
      call vec_norm(v1,N,norm)
      v1=v1/norm
!      call PRINT_MATRIX('start',N,1,v1)

*****************************************

      call lanczos(v1,krylow_ev,krylow_ew,  T_kin,V,trafo,m,N
     $,xeigen )
      call r_PRINT_MATRIX('ew im sil',1,m,krylow_ew)
      call krytra(v1,krylow_ev,N,m,trafo,kry_dvr_ma,kry_FBR)
      call adjun(krylow_ev,m,m,krylow_ev_dagger)
      
      dt=1.d0      

      do j=1,m
      phasen_ew_vec(j)=krylow_ev(1,j)*exp(dt*tphase*krylow_ew(j))
      enddo
     
       do j=1,m
      print*, (phasen_ew_vec(j)), 'phasen_ew_vec'
      enddo
      
      call PRINT_MATRIX('phasen_ew_vec',m,1,phasen_ew_vec)

      call PRINT_MATRIX('ev'  
     $,m,m,krylow_ev)
!      tphase=(0.d0,-1.d0)
      dt=1.d0
      tphase_dt=tphase*dt
      print*,''
!      print*, dt,tphase_dt, 'tphase_dt'
      print*,''
      do j=1,m
      phasen_vec(j)=exp(tphase_dt*krylow_ew(j))
      print*, (phasen_vec(j)), 'phasen_vec'
      enddo
      print*,''

      call vecvec(phasen_vec,krylow_ev_dagger(1,1),m,phasen_vecII)
      do j=1,m
      print*, (phasen_vecII(j)), 'phasen_vecII'
      enddo
      print*,''
!      call PRINT_MATRIX('Phasenvek',m,1,phasen_vec)
!      call PRINT_MATRIX('EIGENV transponiert',m,m,krylow_ev_dagger)
!      call PRINT_MATRIX('Eigenv mit Phasenvekt',m,1,phasen_vecII)

      call matrixvec2(kry_FBR,phasen_vecII,N,m,psi_t) 
      
      call r_PRINT_MATRIX('Psi(t)',N,1,abs(psi_t))
      print*,''

      do j=1,N
!      print*, abs(psi_t(j)), 'psi_t'
      enddo

      STOP
*****************************************
* Zeitloops
*****************************************
      
      total_t=2.d0
      steps=2
      dt=total_t/steps       
      
      do i=2, steps 
     
      v1=0
      do j=1,N
      v1(j,1)=psi_t(j)
      enddo
      
      
      call lanczos(v1,krylow_ev,krylow_ew,  T_kin,V,trafo,m,N
     $,xeigen )
      call r_PRINT_MATRIX('ew im sil',1,m,krylow_ew)
      call krytra(v1,krylow_ev,N,m,trafo,kry_dvr_ma,kry_FBR)
      call adjun(krylow_ev,m,m,krylow_ev_dagger)
      call PRINT_MATRIX('ev'  
     $,m,m,krylow_ev)
      tphase_dt=tphase*dt
      print*,''
      print*,''
      do j=1,m
      phasen_vec(j)=exp(tphase_dt*krylow_ew(j))
!      print*, (phasen_vec(j)), 'phasen_vec'
      enddo
      print*,''

      call vecvec(phasen_vec,krylow_ev_dagger(1,1),m,phasen_vecII)
      do j=1,m
!      print*, (phasen_vecII(j)), 'phasen_vecII'
      enddo
      print*,''

      call matrixvec2(kry_FBR,phasen_vecII,N,m,psi_t) 
      call r_PRINT_MATRIX('Psi(t)bla',N,1,abs(psi_t))
      print*,''

      enddo
 

      end subroutine

      
      
