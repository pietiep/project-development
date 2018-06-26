      subroutine potential(x,V,bdef)

      implicit none
      real*8 bdef(3)
      real*8 x
      complex*16 V

      V= 0.5*bdef(3)*(bdef(2)*x)**2-0.02*x**4

      end subroutine

