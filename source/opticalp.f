      subroutine opticalp(Zix,Nix,eopt)
c
c +---------------------------------------------------------------------
c | Author: Arjan Koning
c | Date  : June 29, 2011
c | Task  : Optical potential for protons
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      character*72 optmodfile
      integer      Zix,Nix,nen,i,Z,A,mw,md
      real         eopt,elow,eup,eint,vloc(19),f,Vc,vcoul,v1loc,v2loc,
     +             v3loc,v4loc,v5loc,w1loc,w2loc,w3,w4,d1loc,d2loc,
     +             d3loc,vso1loc,vso2loc,wso1loc,wso2loc,factor(12)
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
c
c 1. In case of an optical model file, we interpolate between the
c    tabulated values.
c
      optmodfile='                                                     '
      if (Zix.le.numZph.and.Nix.le.numNph) optmodfile=optmod(Zix,Nix,2)
      if (optmodfile(1:1).ne.' ') then
        if (eopt.lt.eomp(Zix,Nix,2,1).or.
     +    eopt.gt.eomp(Zix,Nix,2,omplines(Zix,Nix,2))) goto 100
        do 20 nen=1,omplines(Zix,Nix,2)-1
          elow=eomp(Zix,Nix,2,nen)
          eup=eomp(Zix,Nix,2,nen+1)
          if (elow.le.eopt.and.eopt.le.eup) then
            eint=(eopt-elow)/(eup-elow)
            do 30 i=1,19
              vloc(i)=vomp(Zix,Nix,2,nen,i)+
     +          eint*(vomp(Zix,Nix,2,nen+1,i)-vomp(Zix,Nix,2,nen,i))
   30       continue
            v=v1adjust(2)*vloc(1)
            rv=rvadjust(2)*vloc(2)
            av=avadjust(2)*vloc(3)
            w=w1adjust(2)*vloc(4)
            rw=rwadjust(2)*vloc(5)
            aw=awadjust(2)*vloc(6)
            vd=d1adjust(2)*vloc(7)
            rvd=rvdadjust(2)*vloc(8)
            avd=avdadjust(2)*vloc(9)
            wd=d1adjust(2)*vloc(10)
            rwd=rwdadjust(2)*vloc(11)
            awd=awdadjust(2)*vloc(12)
            vso=vso1adjust(2)*vloc(13)
            rvso=rvsoadjust(2)*vloc(14)
            avso=avsoadjust(2)*vloc(15)
            wso=wso1adjust(2)*vloc(16)
            rwso=rwsoadjust(2)*vloc(17)
            awso=awsoadjust(2)*vloc(18)
            rc=rcadjust(2)*vloc(19)
            goto 200
          endif
   20   continue
      endif
c
c 2. The general energy-dependent form of the optical potential
c    using parameters per nucleus or a global optical model,
c    both from subroutine omppar.
c
c ZZ,Z         : charge number of residual nucleus
c AA,A         : mass number of residual nucleus
c ompglobal    : flag for use of global optical model
c soukhovitskii: subroutine for global optical model parameters for
c                actinides by Soukhovitskii et al.
c f            : E-Ef
c ef           : Fermi energy
c v1loc.....   : help variables
c v1adjust..   : adjustable factors for OMP (default 1.)
c Vc           : Coulomb constant
c vcoul        : Coulomb term of V
c v1,v2,v3     : components for V
c rc,rc0       : Coulomb radius
c onethird     : 1/3
c w1,w2        : components for W
c d1,d2,d3     : components for Wd
c mw,md        : powers for W and Wd
c vso1,vso2    : components for Vso
c wso1,wso2    : components for Wso
c
  100 Z=ZZ(Zix,Nix,0)
      A=AA(Zix,Nix,0)
      if (Z.ge.90.and.ompglobal(Zix,Nix,2)) then
        call soukhovitskii(2,Z,A,eopt)
      else
        f=eopt-ef(Zix,Nix,2)
        rc=rcadjust(2)*rc0(Zix,Nix,2)
        v1loc=v1adjust(2)*v1(Zix,Nix,2)
        v2loc=v2adjust(2)*v2(Zix,Nix,2)
        v3loc=v3adjust(2)*v3(Zix,Nix,2)
        v4loc=v4adjust(2)*7.e-9
c
c Extension to 600 MeV only included if the maximum incident energy
c is above 200 MeV. Otherwise, stick to standard OMP. In this way,
c small discontinuities are avoided.
c
        if (enincmax.gt.200.) then
            v5loc=1.74e-11-5.1e-14*A
          else
            v5loc=0.
        endif
c
c Coulomb term
c
        if (ompglobal(Zix,Nix,2)) then
          Vc=1.73/rc*Z/(A**onethird)
          vcoul=Vc*v1loc*(v2loc-2.*v3loc*f+3.*v4loc*f**2+4.*v5loc*f**3)
        else
          vcoul=0.
        endif
        v=v1loc*(1.-v2loc*f+v3loc*f**2-v4loc*f**3-v5loc*f**4)+vcoul
        rv=rvadjust(2)*rv0(Zix,Nix,2)
        av=avadjust(2)*av0(Zix,Nix,2)
        mw=2
        w1loc=w1adjust(2)*w1(Zix,Nix,2)
        w2loc=w2adjust(2)*w2(Zix,Nix,2)
        if (enincmax.gt.200.) then
          w3=20.
        else
          w3=0.
        endif
        w4=600.
        w=w1loc*f**mw/(f**mw+w2loc**mw)+w3*f**4/(f**4+w4**4)
        rw=rwadjust(2)*rv0(Zix,Nix,2)
        aw=awadjust(2)*av0(Zix,Nix,2)
        vd=0.
        rvd=rvdadjust(2)*rvd0(Zix,Nix,2)
        avd=avdadjust(2)*avd0(Zix,Nix,2)
        md=2
        d1loc=d1adjust(2)*d1(Zix,Nix,2)
        d2loc=d2adjust(2)*d2(Zix,Nix,2)
        d3loc=d3adjust(2)*d3(Zix,Nix,2)
        wd=d1loc*f**md*exp(-d2loc*f)/(f**md+d3loc**md)
        rwd=rwdadjust(2)*rvd0(Zix,Nix,2)
        awd=awdadjust(2)*avd0(Zix,Nix,2)
        vso1loc=vso1adjust(2)*vso1(Zix,Nix,2)
        vso2loc=vso2adjust(2)*vso2(Zix,Nix,2)
        vso=vso1loc*exp(-vso2loc*f)
        rvso=rvsoadjust(2)*rvso0(Zix,Nix,2)
        avso=avsoadjust(2)*avso0(Zix,Nix,2)
        wso1loc=wso1adjust(2)*wso1(Zix,Nix,2)
        wso2loc=wso2adjust(2)*wso2(Zix,Nix,2)
        wso=wso1loc*f**2/(f**2+wso2loc**2)
        rwso=rwsoadjust(2)*rvso0(Zix,Nix,2)
        awso=awsoadjust(2)*avso0(Zix,Nix,2)
      endif
c
c Possible additional energy-dependent adjustment of the geometry
c
c ompadjustF: logical for local OMP adjustment
c adjustF   : subroutine for local optical model geometry adjustment
c factor    : Woods-Saxon multiplication factor
c
  200 if (ompadjustF(2)) then
        call adjustF(2,eopt,factor)
        rv=factor(1)*rv
        av=factor(2)*av
        rw=factor(3)*rw
        aw=factor(4)*aw
        rvd=factor(5)*rvd
        avd=factor(6)*avd
        rwd=factor(7)*rwd
        awd=factor(8)*awd
        rvso=factor(9)*rvso
        avso=factor(10)*avso
        rwso=factor(11)*rwso
        awso=factor(12)*awso
      endif
      return
      end
Copyright (C) 2004  A.J. Koning, S. Hilaire and M.C. Duijvestijn
