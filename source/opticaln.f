      subroutine opticaln(Zix,Nix,eopt)
c
c +---------------------------------------------------------------------
c | Author: Arjan Koning 
c | Date  : February 7, 2006
c | Task  : Optical potential for neutrons
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      character*72 optmodfile
      integer      Zix,Nix,nen,Z,A,i,mw,md
      real         eopt,elow,eup,eint,vloc(19),f,v1loc,v2loc,v3loc,
     +             v4loc,w1loc,w2loc,d1loc,d2loc,d3loc,vso1loc,vso2loc,
     +             wso1loc,wso2loc
c
c ************************ Calculate parameters ************************
c
c Zix          : charge number index for residual nucleus
c Nix          : neutron number index for residual nucleus
c eopt         : incident energy
c numNph       : maximal number of neutrons away from the initial 
c                compound nucleus for multiple pre-equilibrium emission
c numZph       : maximal number of protons away from the initial 
c                compound nucleus for multiple pre-equilibrium emission
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
      optmodfile='                                                     '
      if (Zix.le.numZph.and.Nix.le.numNph) optmodfile=optmod(Zix,Nix,1)
      if (optmodfile(1:1).ne.' ') then
        if (eopt.lt.eomp(Zix,Nix,1,1).or.
     +    eopt.gt.eomp(Zix,Nix,1,omplines(1))) goto 100
        do 20 nen=1,omplines(1)-1
          elow=eomp(Zix,Nix,1,nen)
          eup=eomp(Zix,Nix,1,nen+1)
          if (elow.le.eopt.and.eopt.le.eup) then
            eint=(eopt-elow)/(eup-elow)
            do 30 i=1,19
              vloc(i)=vomp(Zix,Nix,1,nen,i)+
     +          eint*(vomp(Zix,Nix,1,nen+1,i)-vomp(Zix,Nix,1,nen,i))
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
c ZZ,Z         : charge number of residual nucleus
c AA,A         : mass number of residual nucleus
c ompglobal    : flag for use of global optical model
c soukhovitskii: subroutine for global optical model parameters for 
c                actinides by Soukhovitskii et al.
c f            : E-Ef
c v1loc.....   : help variables
c v1adjust..   : adjustable factors for OMP (default 1.)
c ef           : Fermi energy
c v1,v2,v3     : components for V
c w1,w2        : components for W
c d1,d2,d3     : components for Wd
c mw,md        : powers for W and Wd
c vso1,vso2    : components for Vso
c wso1,wso2    : components for Wso
c
  100 Z=ZZ(Zix,Nix,0)
      A=AA(Zix,Nix,0)
      if (Z.ge.90.and.ompglobal(Zix,Nix,1)) then
        call soukhovitskii(1,Z,A,eopt)
      else
        f=eopt-ef(Zix,Nix,1)
        v1loc=v1adjust(1)*v1(Zix,Nix,1)
        v2loc=v2adjust(1)*v2(Zix,Nix,1)
        v3loc=v3adjust(1)*v3(Zix,Nix,1)
        v4loc=v4adjust(1)*7.e-9
        v=v1loc*(1.-v2loc*f+v3loc*f**2-v4loc*f**3)
        rv=rvadjust(1)*rv0(Zix,Nix,1)
        av=avadjust(1)*av0(Zix,Nix,1)
        mw=2
        w1loc=w1adjust(1)*w1(Zix,Nix,1)
        w2loc=w2adjust(1)*w2(Zix,Nix,1)
        w=w1loc*f**mw/(f**mw+w2loc**mw)
        rw=rv
        aw=av
        vd=0.
        rvd=rvdadjust(1)*rvd0(Zix,Nix,1)
        avd=avdadjust(1)*avd0(Zix,Nix,1)
        md=2
        d1loc=d1adjust(1)*d1(Zix,Nix,1)
        d2loc=d2adjust(1)*d2(Zix,Nix,1)
        d3loc=d3adjust(1)*d3(Zix,Nix,1)
        wd=d1loc*f**md*exp(-d2loc*f)/(f**md+d3loc**md)
        rwd=rvd
        awd=avd
        vso1loc=vso1adjust(1)*vso1(Zix,Nix,1)
        vso2loc=vso2adjust(1)*vso2(Zix,Nix,1)
        vso=vso1loc*exp(-vso2loc*f)
        rvso=rvsoadjust(1)*rvso0(Zix,Nix,1)
        avso=avsoadjust(1)*avso0(Zix,Nix,1)
        wso1loc=wso1adjust(1)*wso1(Zix,Nix,1)
        wso2loc=wso2adjust(1)*wso2(Zix,Nix,1)
        wso=wso1loc*f**2/(f**2+wso2loc**2)
        rwso=rvso
        awso=avso
        rc=rcadjust(1)*rc0(Zix,Nix,1)
      endif
      return
      end
Copyright (C) 2004  A.J. Koning, S. Hilaire and M.C. Duijvestijn
