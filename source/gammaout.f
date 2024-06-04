      subroutine gammaout(Zcomp,Ncomp)
c
c +---------------------------------------------------------------------
c | Author: Arjan Koning 
c | Date  : July 7, 2004
c | Task  : Output of gamma-ray strength functions, transmission 
c |         coefficients and cross sections
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      integer Zcomp,Ncomp,Z,N,A,l,nen
      real    e,fstrength
c
c ***** Gamma-ray strength functions and transmission coefficients *****
c
c Zcomp    : charge number index for compound nucleus
c Ncomp    : neutron number index for compound nucleus
c ZZ,Z     : charge number of residual nucleus
c NN,N     : neutron number of residual nucleus
c AA,A     : mass number of residual nucleus
c nuc      : symbol of nucleus
c gamgam   : total radiative width
c dgamgam  : uncertainty in gamgam
c D0       : experimental s-wave resonance spacing in eV
c dD0      : uncertainty in D0
c Dth      : theoretical s-wave resonance spacing
c S0       : s-wave strength function
c dS0      : uncertainty in S0
c swaveth  : theoretical strength function for s-wave
c gnorm    : gamma normalization factor
c gammax   : number of l-values for gamma multipolarity
c sgr      : strength of GR
c ngr      : number of GR
c egr      : energy of GR
c ggr      : width of GR
c kgr      : constant for gamma-ray strength function
c ebegin   : first energy point of energy grid 
c eend     : last energy point of energy grid 
c egrid,e  : outgoing energy grid 
c fstrength: gamma-ray strength function
c Tjl      : transmission coefficients as a function of particle type, 
c            energy, spin and l-value
c
      Z=ZZ(Zcomp,Ncomp,0)
      N=NN(Zcomp,Ncomp,0)
      A=AA(Zcomp,Ncomp,0)    
      write(*,'(/"########## GAMMA STRENGTH FUNCTIONS, TRANSMIS",$)')
      write(*,'("SION COEFFICIENTS AND CROSS SECTIONS ##########")')
      write(*,'(/"Gamma-ray information for Z=",i3," N=",i3,$)') Z,N
      write(*,'(" (",i3,a2,") "/)') A,nuc(Z)   
      write(*,'("S-wave strength function parameters:"/)')
      write(*,'("Exp. total radiative width=",f10.5," eV +/-",f8.5)') 
     +  gamgam(Zcomp,Ncomp),dgamgam(Zcomp,Ncomp)
      write(*,'("Exp. D0                   =",f10.2," eV +/-",f8.2,$)') 
     +  D0(Zcomp,Ncomp),dD0(Zcomp,Ncomp)
      write(*,'(" Theor. D0                   =",f15.2," eV")') Dth
      write(*,'("Exp. S-wave strength func.=",f10.5,"E-4 +/-",$)') 
     +  S0(Zcomp,Ncomp)
      write(*,'(f8.5,$)') dS0(Zcomp,Ncomp)
      write(*,'(" Theor. S-wave strength func.=",f15.5,"E-4")') 
     +  1.e4*swaveth
      write(*,'("Normalization factor      =",f10.5)') gnorm
      do 10 l=1,gammax
        write(*,'(/"Normalized gamma-ray strength functions and ",$)')
        write(*,'("transmission coefficients for l=",i2,/)') l
        write(*,'(" Giant resonance parameters :"/)')
        write(*,'(" sigma0(M",i1,") =",f8.3,$)') l,
     +    sgr(Zcomp,Ncomp,0,l,1)
        write(*,'("       sigma0(E",i1,") =",f8.3,$)') 
     +    l,sgr(Zcomp,Ncomp,1,l,1)
        if (ngr(Zcomp,Ncomp,1,l).eq.2) then
          write(*,'(" and ",f8.3)') sgr(Zcomp,Ncomp,1,l,2)
        else
          write(*,'()')
        endif
        write(*,'("      E(M",i1,") =",f8.3,$)') l,
     +    egr(Zcomp,Ncomp,0,l,1)
        write(*,'("            E(E",i1,") =",f8.3,$)') 
     +    l,egr(Zcomp,Ncomp,1,l,1)
        if (ngr(Zcomp,Ncomp,1,l).eq.2) then
          write(*,'(" and ",f8.3)') egr(Zcomp,Ncomp,1,l,2)
        else
          write(*,'()')
        endif
        write(*,'("  gamma(M",i1,") =",f8.3,$)') l,
     +    ggr(Zcomp,Ncomp,0,l,1)
        write(*,'("        gamma(E",i1,") =",f8.3,$)') 
     +    l,ggr(Zcomp,Ncomp,1,l,1)
        if (ngr(Zcomp,Ncomp,1,l).eq.2) then
          write(*,'(" and ",f8.3)') ggr(Zcomp,Ncomp,1,l,2)
        else
          write(*,'()')
        endif
        write(*,'("      k(M",i1,") =",1p,e14.5,$)') 
     +    l,kgr(Zcomp,Ncomp,0,l)
        write(*,'("      k(E",i1,") =",1p,e14.5/)') 
     +    l,kgr(Zcomp,Ncomp,1,l)
        write(*,'("    E       f(M",i1,")        f(E",i1,")",$)') l,l
        write(*,'("        T(M",i1,")        T(E",i1,")"/)') l,l
        do 20 nen=ebegin(0),eend(0)
          e=egrid(nen)
          write(*,'(f7.3,1p,4e13.5)') e,fstrength(Zcomp,Ncomp,e,0,l)
     +      *gnorm,fstrength(Zcomp,Ncomp,e,1,l)*gnorm,Tjl(0,nen,0,l),
     +      Tjl(0,nen,1,l)
   20   continue
   10 continue
c
c **************** Cross sections for inverse channels *****************
c
c xsreac: reaction cross section
c
      write(*,'(/"Photoabsorption cross sections"/)')
      write(*,'("   E      reaction"/)')
      do 110 nen=ebegin(0),eend(0)
        write(*,'(f7.3,1p,e12.4)') egrid(nen),xsreac(0,nen)
  110 continue
      return
      end
Copyright (C) 2004  A.J. Koning, S. Hilaire and M.C. Duijvestijn
