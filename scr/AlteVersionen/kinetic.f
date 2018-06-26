      subroutine kin_energy(T,T_kin,bdef,N)

      implicit none
      integer N
      real*8 bdef(3)
      complex*16  T(N,N)
      complex*16  T_kin(N,N)

      T_kin=T/(2.0d0*bdef(2))  

      end subroutine

