      
      subroutine initmatrix(btype,anzahlba, bdef,V,trafo,xeigen,T)
      implicit none
      
      integer btype, anzahlba
      real*8  bdef(3)
      complex*16  T(anzahlba,anzahlba)
      complex*16  f_x,y
      real*8      g_x
      complex*16  trafo(anzahlba,anzahlba)
      real*8 xeigen(anzahlba)
      complex*16 V(anzahlba)
      complex*16 psi(anzahlba)
      integer i,j
     
********************************************
*          Transformations - Matrix
********************************************



!Basistype fuer DVR
 
      if(btype.eq.1) then   
         do i = 1, anzahlba
         do j = 1, anzahlba
         trafo(j,i)= f_x(j,i,btype,bdef) 
      !f(x)= zum Auff"ullen der Matrixelemente
         enddo
         enddo
      endif

      
      print *
      print *, "trafo"
      print *
     
       
      call mydiag(anzahlba,xeigen,trafo)     
      do i=1,anzahlba
      xeigen(i)=xeigen(i)!ich will das Potential nicht mit verschieben  -bdef(1)                   
      enddo
      
      

******************************************
*          Impuls**2 in harmonischer Basis
******************************************
         
         do i = 1, anzahlba
         do j = 1, anzahlba
         T(j,i)= g_x(j,i,bdef)
         enddo
         enddo
      
      end subroutine

*****************************************
*     H|PSI>
*****************************************
      subroutine hamilton(t,v,trafo,psi,N,h_psi)
      implicit none
      
      integer N,i
      complex*16 t(N,N), v(N), trafo(N,N)
      complex*16 psi(N),psi_v(N),psi_trafo(N),psi_kin(N)
      complex*16 psi_kin_back(N)
      complex*16 h_psi(N) 
      complex*16 skalpsi
      real*8      test

!      |psi>
      call vecvec(v,psi,N,psi_v)
      
      
      call matrixvec(trafo,psi,N,psi_trafo)
      

!     T|psi>
      call matrixvec(t,psi_trafo,N,psi_kin)

!     Trafo
      call matrix_adj_vec(trafo,psi_kin,N,psi_kin_back)
      
!      
      
      call vecadd(psi_v,psi_kin_back,N,h_psi)
      
      

      call skalpro(psi,h_psi,N,skalpsi)
    
      end subroutine
      
      


      complex*16 function f_x(j,i,btype,bdef)
      
      implicit none
      integer j,i
      integer readpsi, routine
      integer btype, anzahlba
      real*8 bdef(3)

!Harmonische Oszillator-Funktion Ort
      j=j-1
      i=i-1
      
      if(btype.eq.1) then
      if(j.eq.(i+1)) then
      f_x=1.0d0*j
      else if (i.eq.(j+1)) then
      f_x=1.0d0*i
      else 
      f_x=0
      endif
      
      f_x=sqrt(f_x/(2.0d0*bdef(2)*bdef(3))) !2/2 * 1/m * 1/w
      endif
      j=j+1
      i=i+1
    
      endfunction 

      real*8 function g_x(j,i,bdef)
      
      implicit none
      integer j,i
      integer readpsi, routine
      integer anzahlba
      real*8 bdef(3)

!Harmonische Oszillator-Funktion Impuls
      j=j-1
      i=i-1
      
 
      if(j.eq.i) then
      g_x=2.0d0*i+1
      else if (i.eq.(j+2)) then
      g_x=-sqrt(1.0d0*i*(i-1))
      else if (j.eq.(i+2)) then
      g_x=-sqrt(1.0d0*j*(j-1))
      else 
      g_x=0
      endif


      
      g_x=g_x*bdef(2)/(4.0d0)
      j=j+1
      i=i+1
      endfunction        
