        subroutine read_system(name2,readpsi,routine,
     $  btype,anzahlba,bdef,m,typ,dt,total_t)
        implicit none
        
        real*8 bdef(3),dt, total_t

        integer i
        character*20 name, name2

        integer readpsi
        integer routine
        integer btype, anzahlba, m, typ   
         

        call getarg(1,name)
        open(10,file=name)
        read(10,*,err=109) name2    ! directory
        read(10,*,err=109) readpsi ! read wavefunction
        read(10,*,err=109) routine ! type of routine
        read(10,*,err=109) btype
        read(10,*,err=109) anzahlba
        read(10,*,err=109) bdef(1), bdef(2), bdef(3) ! r0, freq, mass 
        read(10,*,err=109) m
        read(10,*,err=109) typ
        read(10,*,err=109) dt
        read(10,*,err=109) total_t
        close(10)

        name=adjustl(name)
        return
               
109     write(6,*) "Cannot read input file. First Argument represents
     . the input file name."
        end subroutine read_system


