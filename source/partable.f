      subroutine partable(Zix,Nix)
c
c +---------------------------------------------------------------------
c | Author: Arjan Koning
c | Date  : August 26, 2004
c | Task  : Write model parameters per nucleus to separate file
c +---------------------------------------------------------------------
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      integer Zix,Nix,Z,A,ibar,l
c
c ****************************** Z and A of nucleus ********************
c
c Zix : charge number index for residual nucleus
c Nix : neutron number index for residual nucleus
c ZZ,Z: charge number of residual nucleus
c AA,A: mass number of residual nucleus
c nuc : symbol of nucleus
c
      Z=ZZ(Zix,Nix,0)
      A=AA(Zix,Nix,0)
      write(11,'("#")') 
      write(11,'("# Parameters for ",i3,a2)')  A,nuc(Z)
c
c ********************** Level density parameters **********************
c
c alev        : level density parameter
c gammald     : gamma-constant for asymptotic level density parameter
c pair        : total pairing correction
c D0          : experimental s-wave resonance spacing in eV
c ibar        : fission barrier
c nfisbar     : number of fission barrier parameters
c deltaW      : shell correction in nuclear mass 
c T           : nuclear temperature
c E0          : constant of temperature formula
c Exmatch     : matching point for Ex
c Ntop        : highest discrete level for temperature matching
c Nlow        : lowest discrete level for temperature matching
c Krotconstant: normalization constant for rotational enhancement
c g           : single-particle level density parameter
c gp          : single-particle proton level density parameter
c gn          : single-particle neutron level density parameter
c
      write(11,'("#")') 
      write(11,'("# Level density")')
      write(11,'("#")') 
      write(11,'("a       ",2i4,f10.5)') Z,A,alev(Zix,Nix)
      write(11,'("gammald ",2i4,f10.5)') Z,A,gammald(Zix,Nix)
      write(11,'("pair    ",2i4,f10.5)') Z,A,pair(Zix,Nix)
      if (D0(Zix,Nix).ne.0.) write(11,'("D0      ",2i4,f10.2)') Z,A,
     +  0.001*D0(Zix,Nix)
      do 10 ibar=0,nfisbar(Zix,Nix)
        write(11,'("deltaW  ",2i4,f10.5,i4)') Z,A,deltaW(Zix,Nix,ibar),
     +    ibar
        write(11,'("T       ",2i4,f10.5,i4)') Z,A,T(Zix,Nix,ibar),ibar
        write(11,'("E0      ",2i4,f10.5,i4)') Z,A,E0(Zix,Nix,ibar),ibar
        write(11,'("Exmatch ",2i4,f10.5,i4)') Z,A,Exmatch(Zix,Nix,ibar),
     +    ibar
        write(11,'("Ntop    ",2i4,2i4)') Z,A,Ntop(Zix,Nix,ibar),ibar
        write(11,'("Nlow    ",2i4,2i4)') Z,A,Nlow(Zix,Nix,ibar),ibar
        write(11,'("Krotconstant ",2i4,f10.5,i4)') Z,A,
     +    Krotconstant(Zix,Nix,ibar),ibar
   10 continue
      write(11,'("g       ",2i4,f10.5)') Z,A,g(Zix,Nix)
      write(11,'("gp      ",2i4,f10.5)') Z,A,gp(Zix,Nix)
      write(11,'("gn      ",2i4,f10.5)') Z,A,gn(Zix,Nix)
c
c ************************ Gamma-ray parameters ************************
c
c gamgam: experimental total radiative width in eV                 
c gammax: number of l-values for gamma multipolarity 
c sgr   : strength of GR     
c egr   : energy of GR
c ggr   : width of GR
c
      write(11,'("#")') 
      write(11,'("# Gamma-ray")')
      write(11,'("#")') 
      write(11,'("gamgam  ",2i4,f10.5)') Z,A,gamgam(Zix,Nix)
      do 110 l=1,gammax                      
        write(11,'("sgr     ",2i4,f8.3," M",i1)') Z,A,
     +    sgr(Zix,Nix,0,l,1),l
        write(11,'("egr     ",2i4,f8.3," M",i1)') Z,A,
     +    egr(Zix,Nix,0,l,1),l
        write(11,'("ggr     ",2i4,f8.3," M",i1)') Z,A,
     +    ggr(Zix,Nix,0,l,1),l
        write(11,'("sgr     ",2i4,f8.3," E",i1)') Z,A,
     +    sgr(Zix,Nix,1,l,1),l
        write(11,'("egr     ",2i4,f8.3," E",i1)') Z,A,
     +    egr(Zix,Nix,1,l,1),l
        write(11,'("ggr     ",2i4,f8.3," E",i1)') Z,A,
     +    ggr(Zix,Nix,1,l,1),l
  110 continue
c
c ************************** Fission parameters ************************
c
c flagfission: flag for fission 
c nfisbar    : number of fission barrier parameters
c fbarrier   : height of fission barrier
c fwidth     : width of fission barrier             
c widthc2    : width of class2 states
c Rtransmom  : normalization constant for moment of inertia for 
c              transition states
c Rclass2mom : normalization constant for moment of inertia for 
c              class 2 states
c
      if (flagfission) then
        write(11,'("#")') 
        write(11,'("# Fission parameters")')
        write(11,'("#")') 
        do 210 ibar=1,nfisbar(Zix,Nix) 
          write(11,'("fisbar  ",2i4,f10.5,i3)') Z,A,
     +      fbarrier(Zix,Nix,ibar),ibar
          write(11,'("fishw   ",2i4,f10.5,i3)') Z,A,
     +      fwidth(Zix,Nix,ibar),ibar
          if (ibar.lt.nfisbar(Zix,Nix)) 
     +      write(11,'("class2width ",2i4,f10.5,i3)') Z,A,
     +      widthc2(Zix,Nix,ibar),ibar
          write(11,'("Rtransmom ",2i4,f10.5,i3)') Z,A,
     +      Rtransmom(Zix,Nix,ibar),ibar
          write(11,'("Rclass2mom ",2i4,f10.5,i3)') Z,A,
     +      Rclass2mom(Zix,Nix,ibar),ibar
  210   continue
      endif
      write(11,'("#---------------------------------------------")') 
      return 
      end
Copyright (C) 2004  A.J. Koning, S. Hilaire and M.C. Duijvestijn
