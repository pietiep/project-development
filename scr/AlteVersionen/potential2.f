      subroutine potential(x,V,bdef)

      implicit none
      real*8 bdef(3)
      real*8 x
      complex*16 V

      if (x.le.-5) then
      V=400
      else if (x.ge.5) then
      V=400
      else
      V= 0.5*bdef(3)*(bdef(2)*x)**2
      endif

      end subroutine

