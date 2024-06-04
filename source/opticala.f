      subroutine opticala(Zix,Nix,eopt)
c
c +---------------------------------------------------------------------
c | Author: Arjan Koning 
c | Date  : December 8, 2005
c | Task  : Optical potential for alpha
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      character*72 optmodfile
      integer      Zix,Nix,i,nen
      real         eopt,elow,eup,eint,vloc(19),vn,rvn,avn,wn,rwn,awn,
     +             vdn,rvdn,avdn,wdn,rwdn,awdn,vp,rvp,avp,wp,rwp,awp,
     +             vdp,rvdp,avdp,wdp,rwdp,awdp
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
      if (Zix.le.numZph.and.Nix.le.numNph) optmodfile=optmod(Zix,Nix,6)
      if (optmodfile(1:1).ne.' ') then
        if (eopt.lt.eomp(Zix,Nix,6,1).or.
     +    eopt.gt.eomp(Zix,Nix,6,omplines(6))) goto 100
        do 20 nen=1,omplines(6)-1
          elow=eomp(Zix,Nix,6,nen)
          eup=eomp(Zix,Nix,6,nen+1)
          if (elow.le.eopt.and.eopt.le.eup) then
            eint=(eopt-elow)/(eup-elow)
            do 30 i=1,19
              vloc(i)=vomp(Zix,Nix,6,nen,i)+
     +          eint*(vomp(Zix,Nix,6,nen+1,i)-vomp(Zix,Nix,6,nen,i))
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
c We use the Watanabe method (S. Watanabe, Nucl. Phys. 8, 484 (1958),
c see also D.G. Madland, Proceedings of a Specialists' Meeting on 
c preequilibrium nuclear reactions, Semmering, Austria, February 10-12 
c 1988, p. 103) to make an alpha potential out of the proton and 
c neutron potential. The simplified formula is:
c                V(a,E)=2*V(n,E/4)+2*V(p,E/4) 
c and similarly for Wd and W. Furthermore,
c              VSO(a,E)=0.
c              WSO(a,E)=0.
c We take a similar weighting for the geometry parameters. 
c
c 1. Neutron potential for E/4.
c
c eopt     : incident energy (help variable)
c opticaln : optical potential for neutrons
c vn,wn,wdn: V, W and Wd for neutrons
c
  100 eopt=0.25*eopt
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
      rwdn=rwd
      awdn=awd         
c
c 2. Proton potential for E/4.
c
c opticalp : optical potential for protons
c vp,wp,wdp: V, W and Wd for protons
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
      rwdp=rwd
      awdp=awd        
c
c 3. Final potential depths: construct V, W, Wd, Vso and Wso.
c
c v1adjust..: adjustable factors for OMP (default 1.)
c
      v=v1adjust(6)*(2.*vn+2.*vp)
      rv=rvadjust(6)*0.5*(rvn+rvp)
      av=avadjust(6)*0.5*(avn+avp)
      w=w1adjust(6)*(2.*wn+2.*wp)
      rw=rvadjust(6)*0.5*(rwn+rwp)
      aw=avadjust(6)*0.5*(awn+awp)
      vd=d1adjust(6)*(2.*vdn+2.*vdp)
      rvd=rvdadjust(6)*0.5*(rvdn+rvdp)
      avd=avdadjust(6)*0.5*(avdn+avdp)   
      wd=d1adjust(6)*(2.*wdn+2.*wdp)
      rwd=rvdadjust(6)*0.5*(rwdn+rwdp)
      awd=avdadjust(6)*0.5*(awdn+awdp)  
      vso=0.
      wso=0.
      eopt=eopt*4.
      return
      end
Copyright (C) 2004  A.J. Koning, S. Hilaire and M.C. Duijvestijn
