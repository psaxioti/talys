      subroutine opticalt(Zix,Nix,eopt)
c
c +---------------------------------------------------------------------
c | Author: Arjan Koning 
c | Date  : December 10, 2009
c | Task  : Optical potential for tritons
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      character*72 optmodfile
      integer      Zix,Nix,nen,i
      real         eopt,elow,eup,eint,vloc(19),vn,rvn,avn,wn,rwn,awn,
     +             vdn,rvdn,avdn,wdn,vp,rvp,avp,wp,rwp,awp,vdp,rvdp,
     +             avdp,wdp,vson,rvson,avson,wson,vsop,rvsop,avsop,wsop,
     +             factor(6)
c
c ******************* Optical model input file *************************
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
c 1 In case of an optical model file, we interpolate between the
c   tabulated values.
c
      optmodfile='                                                     '
      if (Zix.le.numZph.and.Nix.le.numNph) optmodfile=optmod(Zix,Nix,4)
      if (optmodfile(1:1).ne.' ') then
        if (eopt.lt.eomp(Zix,Nix,4,1).or.
     +    eopt.gt.eomp(Zix,Nix,4,omplines(Zix,Nix,4))) goto 100
        do 20 nen=1,omplines(Zix,Nix,4)-1
          elow=eomp(Zix,Nix,4,nen)
          eup=eomp(Zix,Nix,4,nen+1)
          if (elow.le.eopt.and.eopt.le.eup) then
            eint=(eopt-elow)/(eup-elow)
            do 30 i=1,19
              vloc(i)=vomp(Zix,Nix,4,nen,i)+
     +          eint*(vomp(Zix,Nix,4,nen+1,i)-vomp(Zix,Nix,4,nen,i))
   30       continue
            v=v1adjust(4)*vloc(1)
            rv=rvadjust(4)*vloc(2)
            av=avadjust(4)*vloc(3)
            w=w1adjust(4)*vloc(4)
            rw=rvadjust(4)*vloc(5)
            aw=avadjust(4)*vloc(6)
            vd=d1adjust(4)*vloc(7)
            rvd=rvdadjust(4)*vloc(8)
            avd=avdadjust(4)*vloc(9)
            wd=d1adjust(4)*vloc(10)
            rwd=rvdadjust(4)*vloc(11)
            awd=avdadjust(4)*vloc(12)
            vso=vso1adjust(4)*vloc(13)
            rvso=rvsoadjust(4)*vloc(14)
            avso=avsoadjust(4)*vloc(15)
            wso=wso1adjust(4)*vloc(16)
            rwso=rvsoadjust(4)*vloc(17)
            awso=avsoadjust(4)*vloc(18)
            rc=rcadjust(4)*vloc(19)
            goto 200
          endif
   20   continue
      endif                                    
c
c 2. The general energy-dependent form of the optical potential
c    using parameters per nucleus or a global optical model,
c    both from subroutine omppar.
c
c We use the Watanabe method (S. Watanabe, Nucl. Phys. 8, 484 (1958),
c see also D.G. Madland, Proceedings of a Specialists' Meeting on 
c preequilibrium nuclear reactions, Semmering, Austria, February 10-12 
c 1988, p. 103) to make a triton potential out of the proton and 
c neutron potential. The simplified formula is:
c                V(t,E)=2*V(n,E/3)+V(p,E/3) 
c and similarly for Wd and W. Furthermore,
c              VSO(t,E)=VSO(n,E)/3
c              WSO(t,E)=WSO(n,E)/3
c We take a similar weighting for the geometry parameters.  
c
c 1. Neutron potential for E/3.
c
c onethird : 1/3
c eopt     : incident energy (help variable)
c opticaln : optical potential for neutrons 
c vn,wn....: optical model parameters for neutrons   
c
  100 eopt=onethird*eopt
      call opticaln(Zix,Nix,eopt)
      vn=v
      rvn=rv
      avn=av
      wn=w
      rwn=rw
      awn=aw
      vdn=vd
      rvdn=rvd
      avdn=avd
      wdn=wd
c
c 2. Proton potential for E/3.
c
c opticalp : optical potential for protons 
c vp,wp....: optical model parameters for protons   
c
      call opticalp(Zix,Nix,eopt)
      vp=v
      rvp=rv
      avp=av
      wp=w
      rwp=rw
      awp=aw
      vdp=vd
      rvdp=rvd
      avdp=avd
      wdp=wd
c
c 3. Another 2 calls to opticaln for eopt=E to determine Vso and Wso.
c
c vson,wson: Vso and Wso for neutrons
c vsop,wsop: Vso and Wso for protons
c
      eopt=eopt*3.
      call opticaln(Zix,Nix,eopt)
      vson=vso
      rvson=rvso
      avson=avso
      wson=wso
      call opticalp(Zix,Nix,eopt)
      vsop=vso
      rvsop=rvso
      avsop=avso
      wsop=wso
c
c 4. Final potential depths: construct V, W, Wd, Vso and Wso.
c
c v1adjust..: adjustable factors for OMP (default 1.)
c
      v=v1adjust(4)*(2.*vn+vp)
      rv=rvadjust(4)*(2.*rvn+rvp)*onethird
      av=avadjust(4)*(2.*avn+avp)*onethird
      w=w1adjust(4)*(2.*wn+wp)
      rw=rvadjust(4)*(2.*rwn+rwp)*onethird
      aw=avadjust(4)*(2.*awn+awp)*onethird
      vd=d1adjust(4)*(2.*vdn+vdp)
      rvd=rvdadjust(4)*(2.*rvdn+rvdp)*onethird
      avd=avdadjust(4)*(2.*avdn+avdp)*onethird
      wd=d1adjust(4)*(2.*wdn+wdp)
      rwd=rvd
      awd=avd
      vso=vso1adjust(4)*0.5*(vson+vsop)*onethird
      rvso=rvsoadjust(4)*(2.*rvson+rvsop)*onethird
      avso=avsoadjust(4)*(2.*avson+avsop)*onethird
      wso=wso1adjust(4)*0.5*(wson+wsop)*onethird
      rwso=rvso
      awso=avso
c
c Possible additional energy-dependent adjustment of the geometry
c
c ompadjustF: logical for local OMP adjustment
c adjustF   : subroutnie for local optical model geometry adjustment
c factor    : Woods-Saxon multiplication factor
c
  200 if (ompadjustF(4)) then
        call adjustF(4,eopt,factor)
        rv=factor(1)*rv
        av=factor(2)*av
        rwd=factor(3)*rwd
        awd=factor(4)*awd
        rvso=factor(5)*rvso
        avso=factor(6)*avso
      endif
      return
      end
Copyright (C) 2004  A.J. Koning, S. Hilaire and M.C. Duijvestijn
