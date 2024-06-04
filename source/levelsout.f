      subroutine levelsout(Zix,Nix)
c
c +---------------------------------------------------------------------
c | Author: Arjan Koning 
c | Date  : October 13, 2004
c | Task  : Output of discrete levels 
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      integer Zix,Nix,Z,N,A,type,i,k
c
c *************************** Discrete levels **************************
c
c Zix        : charge number index for residual nucleus
c Nix        : neutron number index for residual nucleus
c ZZ,Z       : charge number of residual nucleus
c NN,N       : neutron number of residual nucleus
c AA,A       : mass number of residual nucleus
c nuc        : symbol of nucleus
c parskip    : logical to skip outgoing particle
c parname    : name of particle
c S          : separation energy per particle
c nlev       : number of levels for nucleus
c edis       : energy of level
c jdis       : spin of level
c cparity    : parity of level (character)
c parlev     : parity of level 
c tau        : lifetime of state in seconds
c jassign    : flag for assignment of spin
c passign    : flag for assignment of parity
c ENSDF      : string from original ENSDF discrete level file  
c branchratio: gamma-ray branching ratio to level
c bassign    : flag for assignment of branching ratio
c
      Z=ZZ(Zix,Nix,0)
      N=NN(Zix,Nix,0)
      A=AA(Zix,Nix,0)       
      write(*,'(/"NUCLEAR STRUCTURE INFORMATION FOR Z=",i3," N=",i3,$)')
     +   Z,N
      write(*,'(" (",i3,a2,") ")') A,nuc(Z)
      write(*,'(/"Separation energies:")') 
      write(*,'(/"Particle        S         "/)')
      do 10 type=1,6
        if (parskip(type)) goto 10
        write(*,'(a8,3x,f9.5)') parname(type),S(Zix,Nix,type)
   10 continue
      write(*,'(/"Discrete levels of Z=",i3," N=",i3,$)') Z,N
      write(*,'(" (",i3,a2,") ")') A,nuc(Z)
      write(*,'(/"Number Energy  Spin Parity  Branching  Ratio (%)",$)')
      write(*,'(" Lifetime(sec) Assignment        ENSDF"/)')
      do 20 i=0,nlev(Zix,Nix)
        write(*,'(i3,4x,f7.4,1x,f4.1,3x,a1,$)') i,edis(Zix,Nix,i),
     +    jdis(Zix,Nix,i),cparity(parlev(Zix,Nix,i))
        write(*,'(24(" "),$)')
        if (tau(Zix,Nix,i).ne.0.) then
          write(*,'(2x,1p,e10.3,$)') tau(Zix,Nix,i)
        else
          write(*,'(12(" "),$)')
        endif
        write(*,'(7x,2a1,a18)') jassign(Zix,Nix,i),passign(Zix,Nix,i),
     +    ENSDF(Zix,Nix,i)
        do 40 k=0,i
          if (branchratio(Zix,Nix,i,k).ne.0.) then
            write(*,'(30x,"--->",i3,2x,f8.4,18x,a1)') k,
     +        branchratio(Zix,Nix,i,k)*100.,bassign(Zix,Nix,i,k)
          endif
   40   continue
   20 continue
      return 
      end
Copyright (C) 2004  A.J. Koning, S. Hilaire and M.C. Duijvestijn
