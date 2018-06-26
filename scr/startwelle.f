****************************************    
*  Wellenfunktioninitialisierung  Psi  *
****************************************      


      subroutine wavefkt(anzahlba,xeigen,bdef,trafo,psi_FBR)
      implicit none
      
      integer anzahlba
      real*8  bdef(3)
      real*8    xeigen(anzahlba)
      complex*16 norm
      integer i
      real*8 pi
      complex*16 trafo(anzahlba,anzahlba)
      parameter (pi=3.14159265358)
      complex*16 psi_DVR(anzahlba), psi_FBR(anzahlba)
      
      


      do i=1,anzahlba
      psi_DVR(i)= exp(-0.5d0*bdef(2)*bdef(3)*(xeigen(i)-bdef(1))**2)
      enddo

      call vec_norm(psi_DVR,anzahlba,norm)
      
      do i=1,anzahlba
      psi_DVR(i)=psi_DVR(i)/norm
      enddo
      
      call matrixvec(trafo, psi_DVR, anzahlba,psi_FBR)
      end subroutine

      
      subroutine readpsi(psi, anzahlba)
      implicit none

      integer     anzahlba, i
      complex*16  psi(anzahlba)


      do i=1,anzahlba
      read(3) psi(i)
      enddo

      end subroutine


      subroutine writepsi(psi, anzahlba)
      implicit none

      integer     anzahlba, i
      complex*16  psi(anzahlba)


      do i=1,anzahlba
      write(3) psi(i)
      enddo

      end subroutine

