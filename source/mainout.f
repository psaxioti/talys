      subroutine mainout
c
c +---------------------------------------------------------------------
c | Author: Arjan Koning 
c | Date  : December 15, 2006    
c | Task  : Main output
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      integer Zcomp,Ncomp,Zix,Nix,type,i
c
c *************************** Code and version *************************
c
      write(*,'(/"    TALYS-0.72      (Version: December 15, 2006)"/)')
      write(*,'(" Copyright (C) 2004  A.J. Koning, S. Hilaire ",
     +  "and M.C. Duijvestijn")')
      write(*,'(24x," NRG          CEA              NRG"/)')
      write(*,'(" Dimensions - Cross sections: mb, Energies: MeV, ",
     +  "Angles: degrees")')
c
c ***************** Write input file and default parameters ************
c
c inputout: subroutine to write the input parameters
c
      call inputout
c
c ********************** Main nuclear parameters ***********************
c
c Zcomp  : charge number index for compound nucleus
c Ncomp  : neutron number index for compound nucleus
c Zindex : charge number index for residual nucleus
c Nindex : neutron number index for residual nucleus
c k0     : index of incident particle
c parname: name of particle
c parmass: mass of particle in a.m.u.
c Atarget: mass number of target nucleus
c Starget: symbol of target nucleus
c tarmass: mass of target nucleus
c Ltarget: excited level of target
c edis   : energy of level
c jdis   : spin of level
c cparity: parity of level (character)
c parlev : parity of level 
c tau    : lifetime of state in seconds
c parskip: logical to skip outgoing particle
c numinc : number of incident energies
c eninc  : incident energy in MeV
c parsym : symbol of particle
c Q      : Q-value for target nucleus 
c
      Zcomp=0
      Ncomp=0
      Zix=Zindex(Zcomp,Ncomp,k0)
      Nix=Nindex(Zcomp,Ncomp,k0)
      write(*,'(/" ########## BASIC REACTION PARAMETERS ##########"/)')
      write(*,'(" Projectile           : ",a8,4x, 
     +  "Mass in a.m.u.      : ",f10.6)') parname(k0),parmass(k0)
      write(*,'(" Target               : ",i3,a2,7x, 
     +  "Mass in a.m.u.      : ",f10.6)') Atarget,Starget,tarmass
      if (Ltarget.ne.0) then
        write(*,'(/" Excited target level : Number  Energy  ",
     +    "Spin Parity Lifetime(sec)")')
        write(*,'(24x,i3,4x,f7.4,2x,f4.1,3x,a1,4x,1p,e10.3)') 
     +    Ltarget,edis(Zix,Nix,Ltarget),jdis(Zix,Nix,Ltarget),
     +    cparity(parlev(Zix,Nix,Ltarget)),tau(Zix,Nix,Ltarget)
      endif
      write(*,'(/" Included channels:")')
      do 10 type=-1,6
        if (parskip(type)) goto 10
        write(*,'(21x,a8)') parname(type)
   10 continue
      if (numinc.eq.1) then
        write(*,'(/," 1 incident energy (LAB):"/)') 
      else
        write(*,'(/,1x,i3," incident energies (LAB):"/)') numinc
      endif
      do 20 i=1,numinc
        if (eninc(i).lt.0.001) then
          write(*,'(1x,1p,e10.3)') eninc(i)
        else
          write(*,'(1x,f10.3)') eninc(i)
        endif
   20 continue
      write(*,'(/" Q-values for binary reactions:"/)') 
      do 30 type=0,6
        if (parskip(type)) goto 30
        write(*,'(" Q(",a1,",",a1,"):",f9.5)') parsym(k0),parsym(type),
     +    Q(type)
   30 continue
c
c * Write nuclear structure parameters for target and compound nucleus *
c
c flaglevels   : flag for output of discrete level information
c levelsout    : subroutine for output of discrete levels
c flagdensity  : flag for output of level densities
c densityout   : subroutine for output of level density parameters
c flagfisout   : flag for output of fission information
c fissionparout: subroutine for output for fission parameters
c strucwrite   : flag for output of nuclear structure info 
c
      if (flaglevels) call levelsout(Zix,Nix)
      if (flagdensity) call densityout(Zix,Nix)
      if (flagfisout) call fissionparout(Zix,Nix)
      strucwrite(Zix,Nix)=.true.
      if (parskip(0)) return
      if (k0.ne.0) then
        if (flaglevels) call levelsout(Zcomp,Ncomp)
        if (flagdensity) call densityout(Zcomp,Ncomp)
        if (flagfisout) call fissionparout(Zcomp,Ncomp)
        strucwrite(Zcomp,Ncomp)=.true.
      endif
      return 
      end
Copyright (C) 2004  A.J. Koning, S. Hilaire and M.C. Duijvestijn
