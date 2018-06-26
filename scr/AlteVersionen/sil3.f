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
*     Time - Steps
*****************************************
      total_t=2.d0
      steps=2
      dt=total_t/steps       

*****************************************
*      Startwellenfunktion v1
*****************************************
      v1=0
      
      do j=1,N
      v1(j,1)=psi(j)
      enddo
      
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
      call r_PRINT_MATRIX('ew im sil',1,m,krylow_ew)
      call krytra(v1,krylow_ev,N,m,trafo,kry_dvr_ma,kry_FBR)
     
      
      print*,'' 
      do j=1,m
      phasen_ew_vec(j)=krylow_ev(1,j)*exp(dt*tphase*krylow_ew(j))
      enddo
     
      
      call compl_ PRINT_MATRIX('phasen_ew_vec',m,1,phasen_ew_vec)

      call compl_PRINT_MATRIX('ev'  
     $,m,m,krylow_ev)
      print*,''
      print*,''
      print*,''


      call matrixvec2(kry_FBR,phasen_ew_vec,N,m,psi_t) 
      
      call compl_PRINT_MATRIX('Psi(t)',N,1,(psi_t))
      call r_PRINT_MATRIX('abs(Psi(t))',N,1,abs(psi_t))
      print*,''


      STOP
*****************************************
*      Further time step Iterations
*****************************************
      
      v1=0
      
      do j=1,m
      v1(j,1)=psi_t(j)
      enddo
      
 

      end subroutine

      
      
