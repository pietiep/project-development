
        subroutine initmatrix(Trafo,btype,bdef,hermite)
        implicit none
        
        integer btype(2),bdef(3)
        real*8  Trafo(btype(2),btype(2))
        real*8  f_x(btype(2),btype(2))
        real*8  hermite(btype(2),btype(2))
        integer i,j

        if(btype(1).eq.1) then       !Basistype fuer DVR
                do i = 0, btype(2)
                  do j = 0, btype(2)
                hermite(j,i)= f_x(j,i) !f(x)= zum Auff"ullen der Matrixelemente
                enddo
                enddo
                endif
        end subroutine




        real*8 function f_x(j,i)
        
        implicit none
        integer j,i
        integer readpsi, routine
        integer btype(2)
        real*8 bdef(3)
        call read_system(readpsi,routine,btype,bdef) 

!Harmonische Oszillator-Funktion
        
        if(btype(1).eq.1) then
        f_x=0.0d0
        elseif(j.eq.(i+1)) then
        f_x=1.0d0*(i+1)
        elseif (j.eq.(i-1)) then
        f_x=1.0d0*i
        endif
 
        f_x=sqrt(f_x/(2.0d0*bdef(2)*bdef(3))) !1/2 * 1/m * 1/w

        endfunction        
