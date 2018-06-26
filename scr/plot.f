      subroutine plotwfk(dim,name2,x)
      implicit none
      
      integer     i,j,dim
      integer     max
      parameter   (max=10000)
      real*8      x(max)
      complex*16  psi(dim)
      character*100 name,number
      character*20  name2

      open(10,file=trim(name2)//'/psi',form='unformatted')

      j=10000
      do while(.true.)
     
      do i=1,dim
      read(10,end=999) psi(i)
      enddo

      write(number,*) j
      number=adjustL(number)
      name=trim(name2)//'/psi.'//trim(number)

      open(12,file=trim(name))
      do i=1,dim
      write(12,'(20f20.10)') x(i), abs(psi(i))**2, real(psi(i)),
     .      aimag(psi(i))
      enddo
      close(12)

      j=j+1
      enddo


999   CONTINUE
      end
