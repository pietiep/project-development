      program ewpd
      implicit none

      character*10 name

      integer readpsi
      integer routine
      integer btype,anzahlba
      real*8 bdef(3)
      integer i, j, m, typ
      complex*16, allocatable :: V(:)
      complex*16, allocatable :: T(:,:)
      complex*16, allocatable :: T_kin(:,:)
      real*8, allocatable :: xeigen(:)
      complex*16, allocatable :: trafo(:,:) 
      complex*16, allocatable :: psi(:)
      complex*16, allocatable :: h_psi(:)
      complex*16, allocatable :: v1(:,:)
      complex*16, allocatable :: kry_trafo(:,:)
      real*8, allocatable :: krylow_ew(:)
      complex*16, allocatable :: kry_dvr_ma(:,:)

      ! Auslesen der Konfigurationsdatei
      call read_system(readpsi,routine,btype,anzahlba,bdef,m,typ)
      print *, typ 
       ! Speicher Allokieren
      allocate(V(anzahlba))
      allocate(xeigen(anzahlba))
      allocate(T(anzahlba,anzahlba))
      allocate(T_kin(anzahlba,anzahlba))
      allocate(trafo(anzahlba,anzahlba))
      allocate(psi(anzahlba))
      allocate(h_psi(anzahlba))
      allocate(v1(anzahlba,m))
      allocate(kry_trafo(m,m))
      allocate(krylow_ew(m))
      allocate(kry_dvr_ma(anzahlba,m))
!      print *, "readpsi"
!      print *, readpsi
!      print *, "routine"
!      print *, routine
!      print *, "btype"
!      print *, btype
!      print *, "anzahlba"
!      print *, anzahlba
!      write(6,*)bdef(1:3)
       
!      print *, "Gitterpunkte" 

      ! Aufbau des DVRs, des Potentialvektors,...
      call initmatrix(btype,anzahlba, bdef,V,trafo,xeigen,T_kin)

!       call kin_energy(T,T_kin,bdef,anzahlba)
       !Wellenfunktion initiallisieren
      call wavefkt(anzahlba,xeigen,bdef,psi)

       
       !Potential-Vektor 
      do i=1,anzahlba
      call potential(xeigen(i),V(i),bdef)
      enddo

!      call hamilton(T_kin,V,trafo,psi,anzahlba,h_psi)
! Krylowvektoren v1(N,m) in DVR. Eigenvektoren kry_Trafo in
! Krylowdarstellung
!      call lanczos(v1,kry_trafo,krylow_ew,  T_kin,V,trafo,m,anzahlba
!     $,xeigen )
!      call krytra(v1,kry_trafo,anzahlba,m,trafo,kry_dvr_ma)

***********************************************************
* Short Iterativ Lanczos
***********************************************************

      call sil(v1,kry_trafo,krylow_ew, T_kin,V,trafo,m,anzahlba)

***********************************************************





!      call PRINT_MATRIX('trafo',anzahlba,anzahlba,trafo)
!      call PRINT_Vector_MATRIX('plot',anzahlba,m,kry_dvr_ma,xeigen)


!      call r_PRINT_MATRIX('Gitter',anzahlba,1,xeigen)      




      end
