
      implicit none
      
      integer     i,j,dim
      integer     max
      parameter   (max=10000)
      real*8      psi(max),x(max), real_psi(max), imag_psi(max)
      real*8      psi1(max), psi2(max)
      character*100 name,number

      open(3,file='psi_DVR')
      write(6,*) 'Number of grid points?'
      read(5,*)   dim

      j=10000
      do while(.true.)
     
      read(3,*,end=999)
      read(3,*,end=999)
     
      do i=1,dim
      read(3,*) x(i), psi(i), psi1(i), psi2(i)
      enddo

      write(number,*) j
      number=adjustL(number)
      name='psi.'//trim(number)

      open(12,file=trim(name))
      do i=1,dim
      write(12,'(20f20.10)') x(i), psi(i), psi1(i), psi2(i)
      enddo

      close(12)


      j=j+1
      enddo


999   CONTINUE
      end
