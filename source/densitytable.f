      subroutine densitytable(Zix,Nix)
c
c +---------------------------------------------------------------------
c | Author: Arjan Koning and Marieke Duijvestijn
c | Date  : July 5, 2004
c | Task  : Tabulated level densities
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      logical          lexist
      character*4      denchar
      character*90     denfile
      integer          Zix,Nix,Z,A,ibar,iloop,ia,nex,J
      double precision ld2j1(0:numJ)
c
c *********** Tabulated level densities from Goriely *******************
c
c Zix       : charge number index for residual nucleus
c Nix       : neutron number index for residual nucleus
c ZZ,Z      : charge number of residual nucleus
c AA,A,ia   : mass number of residual nucleus 
c ldexist   : flag for existence of level density table
c denfile   : level density parameter file
c ldtable   : level density from table
c ldtottable: total level density from table
c
      Z=ZZ(Zix,Nix,0)
      A=AA(Zix,Nix,0)
      denchar='z   '
      write(denchar(2:4),'(i3.3)') Z
      ibar=0
      do 10 iloop=0,2
        ldexist(Zix,Nix,iloop)=.false.
        if(iloop.eq.0) 
     +    denfile=path(1:lenpath)//'density/ground/goriely/'//denchar
        if(iloop.eq.1) denfile=path(1:lenpath)//
     +    'density/fission/goriely/inner/'//denchar
        if(iloop.eq.2) denfile=path(1:lenpath)//
     +    'density/fission/goriely/outer/'//denchar   
        inquire (file=denfile,exist=lexist)
        if (lexist) then
          open (unit=2,status='old',file=denfile)
   20     read(2,'(/31x,i3//)',end=10) ia
          if (A.ne.ia) then
            do 30 nex=1,56
              read(2,'()')
   30       continue
            goto 20
          else    
            ldexist(Zix,Nix,ibar)=.true.
            do 40 nex=1,55
              read(2,'(24x,e9.2,9x,30e9.2)',err=100) 
     +          ldtottable(Zix,Nix,nex,ibar),(ld2j1(J),J=0,29)
              do 50 J=0,29
                ldtable(Zix,Nix,nex,J,ibar)=ld2j1(J)
   50         continue
   40       continue
            do 60 J=0,29
              ldtable(Zix,Nix,0,J,ibar)=0.
   60       continue
            ibar=ibar+1
          endif
        endif
   10 continue
      return
  100 write(*,'("TALYS-error: Wrong level density table for",$)') 
      write(*,'(" Z=",i3," A=",i3)') Z,A
      stop    
      end
Copyright (C) 2004  A.J. Koning, S. Hilaire and M.C. Duijvestijn
