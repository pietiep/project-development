
      implicit none
      
      integer     i,j,dim
      integer     max
      parameter   (max=10000)
      real*8      psi(max),x(max), real_psi(max), imag_psi(max)
      character*100 name,number

      open(3,file='psi_DVR_reim')
      write(6,*) 'Number of grid points?'
      read(5,*)   dim

      j=1
      do while(.true.)
     
      read(3,*,end=999)
      read(3,*,end=999)
     
      do i=1,dim
      read(3,*) x(i), real_psi(i), imag_psi(i)
      enddo

      write(number,*) j
      number=adjustL(number)
      name='psireim.'//trim(number)

      open(12,file=trim(name))
      do i=1,dim
      write(12,*) x(i), real_psi(i), imag_psi(i)
      enddo

      close(12)


      j=j+1
      enddo


999   CONTINUE
      end
