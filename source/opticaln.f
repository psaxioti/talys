      subroutine opticaln(Zix,Nix,eopt)
c
c +---------------------------------------------------------------------
c | Author: Arjan Koning 
c | Date  : July 7, 2004
c | Task  : Optical potential for neutrons
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      character*72 optmodfile
      integer      Zix,Nix,nen,i,mw,md
      real         eopt,elow,eup,eint,vloc(19),f
c
c ************************ Calculate parameters ************************
c
c Zix          : charge number index for residual nucleus
c Nix          : neutron number index for residual nucleus
c eopt         : incident energy
c optmodfile   : file with optical model parameters
c optmod       : file with optical model parameters
c eomp         : energies on optical model file
c omplines     : number of lines on optical model file
c elow,eup,eint: help variables
c vloc         : interpolated optical model parameters 
c vomp         : optical model parameters from file
c v,rv,av      : real volume potential, radius, diffuseness
c vd,rvd,avd   : real surface potential, radius, diffuseness
c w,rw,aw      : imaginary volume potential, radius, diffuseness
c wd,rwd,awd   : imaginary surface potential, radius, diffuseness
c vso,rvso,avso: real spin-orbit potential, radius, diffuseness
c wso,rwso,awso: imaginary spin-orbit potential, radius, diffuseness
c rc           : Coulomb radius
c
c 1. In case of an optical model file, we interpolate between the
c    tabulated values.
c
      optmodfile=optmod(Zix,Nix,1)
      if (optmodfile(1:1).ne.' ') then
        if (eopt.lt.eomp(1,1).or.eopt.gt.eomp(1,omplines(1))) goto 100
        do 20 nen=1,omplines(1)-1
          elow=eomp(1,nen)
          eup=eomp(1,nen+1)
          if (elow.le.eopt.and.eopt.le.eup) then
            eint=(eopt-elow)/(eup-elow)
            do 30 i=1,19
              vloc(i)=vomp(1,nen,i)+eint*(vomp(1,nen+1,i)-vomp(1,nen,i))
   30       continue
            v=vloc(1)
            rv=vloc(2)
            av=vloc(3)
            w=vloc(4)
            rw=vloc(5)
            aw=vloc(6)
            vd=vloc(7)
            rvd=vloc(8)
            avd=vloc(9)
            wd=vloc(10)
            rwd=vloc(11)
            awd=vloc(12)
            vso=vloc(13)
            rvso=vloc(14)
            avso=vloc(15)
            wso=vloc(16)
            rwso=vloc(17)
            awso=vloc(18)
            rc=vloc(19)
            return
          endif
   20   continue
      endif
c
c 2. The general energy-dependent form of the optical potential
c    using parameters per nucleus or the global optical model,
c    both from subroutine omppar.
c
c f        : E-Ef
c ef       : Fermi energy
c v1,v2,v3 : components for V
c w1,w2    : components for W
c d1,d2,d3 : components for Wd
c mw,md    : powers for W and Wd
c vso1,vso2: components for Vso
c wso1,wso2: components for Wso
c
  100 f=eopt-ef(Zix,Nix,1)
      v=v1(Zix,Nix,1)*(1.-v2(Zix,Nix,1)*f+v3(Zix,Nix,1)*f**2
     +  -7.e-9*f**3)
      rv=rv0(Zix,Nix,1)
      av=av0(Zix,Nix,1)
      mw=2
      w=w1(Zix,Nix,1)*f**mw/(f**mw+w2(Zix,Nix,1)**mw)
      rw=rv
      aw=av
      vd=0.
      rvd=rvd0(Zix,Nix,1)
      avd=avd0(Zix,Nix,1)
      md=2
      wd=d1(Zix,Nix,1)*f**md*exp(-d2(Zix,Nix,1)*f)/
     +  (f**md+d3(Zix,Nix,1)**md)
      rwd=rvd
      awd=avd
      vso=vso1(Zix,Nix,1)*exp(-vso2(Zix,Nix,1)*f)
      rvso=rvso0(Zix,Nix,1)
      avso=avso0(Zix,Nix,1)
      wso=wso1(Zix,Nix,1)*f**2/(f**2+wso2(Zix,Nix,1)**2)
      rwso=rvso
      awso=avso
      rc=rc0(Zix,Nix,1)
      return
      end
Copyright (C) 2004  A.J. Koning, S. Hilaire and M.C. Duijvestijn
