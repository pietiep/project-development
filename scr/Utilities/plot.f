      program plot
      implicit none

      character*20 name !Name of the result files

      integer readpsi
      integer routine
      integer btype,anzahlba
      real*8 bdef(3),dt,total_t
      complex*16 tphase
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
      call read_system(name,readpsi,routine,btype,anzahlba,
     $ bdef,m,typ,dt,total_t)
      
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
      open(3,file=trim(name)//'/psi',form='unformatted')

     
      ! Aufbau des DVRs, des Potentialvektors,...
      call initmatrix(btype,anzahlba, bdef,V,trafo,xeigen,T_kin)

      ! plot wavefunction
      call plotwfk(anzahlba,name,xeigen)

      end
