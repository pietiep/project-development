      subroutine krytra(kry_vec,kry_eigen,N,m,trafo,kry_dvr_ma,
     $ kry_FBR)
      implicit none
      
      integer N,m, i
      
      complex*16 kry_vec(N,m),kry_eigen(m,m),kry_FBR(N,m),
     $ kry_dvr(N),kry_dvr_ma(N,m),trafo(N,N),trafo_dagger(N,N)


      call matr_mul(kry_vec,kry_eigen,N,m,m,kry_FBR)
     
!      call compl_PRINT_MATRIX('Eigenzustaende in FBR'
!     $,N,m,kry_FBR)
!      call r_PRINT_MATRIX('Eigenzustaende in FBR'
!     $,N,m,abs(kry_FBR))

      call matrix_adj_vec(trafo,kry_FBR(1,1),N,kry_dvr)

!      call PRINT_MATRIX('Eigenzustaende in DVR'
!     $,N,1,kry_dvr)

      call adjun(trafo,N,N,trafo_dagger)
!      call compl_PRINT_MATRIX('trafo',N,N,trafo)
!      call compl_PRINT_MATRIX('trafo_dagger',N,N,trafo_dagger)
      
      call matr_mul(trafo_dagger,kry_FBR,N,N,m,kry_dvr_ma)

!      call PRINT_MATRIX('Eigenzustaende in DVR'
!     $,N,m,kry_dvr_ma)


      end subroutine
