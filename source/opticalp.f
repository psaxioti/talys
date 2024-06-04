      subroutine opticalp(Zix,Nix,eopt)
c
c +---------------------------------------------------------------------
c | Author: Arjan Koning 
c | Date  : July 7, 2004
c | Task  : Optical potential for protons
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      character*72 optmodfile
      integer      Zix,Nix,nen,i,Z,A,mw,md
      real         eopt,elow,eup,eint,vloc(19),f,Vc,vcoul
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
c
c 1. In case of an optical model file, we interpolate between the
c    tabulated values.
c
      optmodfile=optmod(Zix,Nix,2)
      if (optmodfile(1:1).ne.' ') then
        if (eopt.lt.eomp(2,1).or.eopt.gt.eomp(2,omplines(2))) goto 100
        do 20 nen=1,omplines(2)-1
          elow=eomp(2,nen)
          eup=eomp(2,nen+1)
          if (elow.le.eopt.and.eopt.le.eup) then
            eint=(eopt-elow)/(eup-elow)
            do 30 i=1,19
              vloc(i)=vomp(2,nen,i)+eint*(vomp(2,nen+1,i)-vomp(2,nen,i))
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
c    using parameters per nucleus or a global optical model,
c    both from subroutine omppar.
c
c f        : E-Ef
c ef       : Fermi energy
c ompglobal: flag for use of global optical model
c ZZ,Z     : charge number of residual nucleus  
c AA,A     : mass number of residual nucleus  
c Vc       : Coulomb constant
c vcoul    : Coulomb term of V 
c v1,v2,v3 : components for V
c rc,rc0   : Coulomb radius
c onethird : 1/3
c w1,w2    : components for W
c d1,d2,d3 : components for Wd
c mw,md    : powers for W and Wd
c vso1,vso2: components for Vso
c wso1,wso2: components for Wso
c
  100 f=eopt-ef(Zix,Nix,2)
      if (ompglobal(Zix,Nix,2)) then
        Z=ZZ(Zix,Nix,0)
        A=AA(Zix,Nix,0)
        Vc=1.73/rc0(Zix,Nix,2)*Z/(A**onethird)
        vcoul=Vc*v1(Zix,Nix,2)*(v2(Zix,Nix,2)-2.*v3(Zix,Nix,2)*f+
     +    3.*7.e-9*f**2)
      else
        vcoul=0.
      endif
      v=v1(Zix,Nix,2)*(1.-v2(Zix,Nix,2)*f+v3(Zix,Nix,2)*f**2
     +  -7.e-9*f**3)+vcoul
      rv=rv0(Zix,Nix,2)
      av=av0(Zix,Nix,2)
      mw=2
      w=w1(Zix,Nix,2)*f**mw/(f**mw+w2(Zix,Nix,2)**mw)
      rw=rv
      aw=av
      vd=0.
      rvd=rvd0(Zix,Nix,2)
      avd=avd0(Zix,Nix,2)
      md=2
      wd=d1(Zix,Nix,2)*f**md*exp(-d2(Zix,Nix,2)*f)/
     +  (f**md+d3(Zix,Nix,2)**md)
      rwd=rvd
      awd=avd
      vso=vso1(Zix,Nix,2)*exp(-vso2(Zix,Nix,2)*f)
      rvso=rvso0(Zix,Nix,2)
      avso=avso0(Zix,Nix,2)
      wso=wso1(Zix,Nix,2)*f**2/(f**2+wso2(Zix,Nix,2)**2)
      rwso=rvso
      awso=avso
      rc=rc0(Zix,Nix,2)
      return
      end                     
Copyright (C) 2004  A.J. Koning, S. Hilaire and M.C. Duijvestijn
