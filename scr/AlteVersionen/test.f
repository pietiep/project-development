      subroutine krytra(kry_vec,kry_eigen,N,m,trafo,kry_dvr_ma)
      implicit none
      
      integer N,m
      
      complex*16 kry_vec(N,m),kry_eigen(m,m),kry_FBR(N,m),
     $ kry_dvr(N),kry_dvr_ma(N,m),trafo(N,N)

      call matr_mul(kry_FBR,kry_vec,kry_eigen,N,m,m)
     
       call PRINT_MATRIX('Eigenzustaende in FBR'
     $,N,m,kry_FBR)

      call matrix_adj_vec(trafo,kry_FBR(1,1),N,kry_dvr)

!      call PRINT_MATRIX('Eigenzustaende in DVR'
!     $,N,1,kry_dvr)

      call matr_mul(kry_dvr_ma,Transpose(trafo),kry_FBR,N,N,m)

!      call PRINT_MATRIX('Eigenzustaende in DVR'
!     $,N,m,kry_dvr_ma)


      end subroutine
