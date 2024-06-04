      subroutine opticala(Zix,Nix,eopt)
c
c +---------------------------------------------------------------------
c | Author: Arjan Koning and Stephane Goriely
c | Date  : May 11, 2011
c | Task  : Optical potential for alpha
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      character*72 optmodfile
      integer      Zix,Nix,i,nen
      real         eopt,elow,eup,eint,vloc(19),e,vn,rvn,avn,wn,rwn,awn,
     +             vdn,rvdn,avdn,wdn,rwdn,awdn,vp,rvp,avp,wp,rwp,awp,
     +             vdp,rvdp,avdp,wdp,rwdp,awdp,factor(12)
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
     +    eopt.gt.eomp(Zix,Nix,6,omplines(Zix,Nix,6))) goto 100
        do 20 nen=1,omplines(Zix,Nix,6)-1
          elow=eomp(Zix,Nix,6,nen)
          eup=eomp(Zix,Nix,6,nen+1)
          if (elow.le.eopt.and.eopt.le.eup) then
            eint=(eopt-elow)/(eup-elow)
            do 30 i=1,19
              vloc(i)=vomp(Zix,Nix,6,nen,i)+
     +          eint*(vomp(Zix,Nix,6,nen+1,i)-vomp(Zix,Nix,6,nen,i))
   30       continue
            v=v1adjust(6)*vloc(1)
            rv=rvadjust(6)*vloc(2)
            av=avadjust(6)*vloc(3)
            w=w1adjust(6)*vloc(4)
            rw=rwadjust(6)*vloc(5)
            aw=awadjust(6)*vloc(6)
            vd=d1adjust(6)*vloc(7)
            rvd=rvdadjust(6)*vloc(8)
            avd=avdadjust(6)*vloc(9)
            wd=d1adjust(6)*vloc(10)
            rwd=rwdadjust(6)*vloc(11)
            awd=awdadjust(6)*vloc(12)
            vso=vso1adjust(6)*vloc(13)
            rvso=rvsoadjust(6)*vloc(14)
            avso=avsoadjust(6)*vloc(15)
            wso=wso1adjust(6)*vloc(16)
            rwso=rwsoadjust(6)*vloc(17)
            awso=awsoadjust(6)*vloc(18)
            rc=rcadjust(6)*vloc(19)
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
c e        : incident energy (help variable)
c opticaln : optical potential for neutrons
c vn,wn,wdn: V, W and Wd for neutrons
c
  100 e=0.25*eopt
      call opticaln(Zix,Nix,e)
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
      call opticalp(Zix,Nix,e)
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
      rw=rwadjust(6)*0.5*(rwn+rwp)
      aw=awadjust(6)*0.5*(awn+awp)
      vd=d1adjust(6)*(2.*vdn+2.*vdp)
      rvd=rvdadjust(6)*0.5*(rvdn+rvdp)
      avd=avdadjust(6)*0.5*(avdn+avdp)
      wd=d1adjust(6)*(2.*wdn+2.*wdp)
      rwd=rwdadjust(6)*0.5*(rwdn+rwdp)
      awd=awdadjust(6)*0.5*(awdn+awdp)
      vso=0.
      wso=0.
c
c S. Goriely: inclusion of the alpha OMP of Mc Fadden & Satchler
c for alphaomp=2, folding model for alphaomp=3,4,5
c alphaomp=1 --> Watanabe potential
c alphaomp=2 --> Global OMP of L. McFadden, G.R. Satchler,
c Nucl. Phys. 84 (1966) 177
c alphaomp=3-5 --> Global folding alpha OMP of P. Demetriou, C. Grama,
c S. Goriely, Nucl Phys A707, 253 (2002)
c
c Overwrite some of the previous values.
c
      if (alphaomp.eq.2) then
        v=185.0
        rv=1.40
        av=0.52
        w=25.0
        rw=rv
        aw=av
        vd=0.
        wd=0.
        vso=0.
        wso=0.
      endif
c
c Possible additional energy-dependent adjustment of the geometry
c
c ompadjustF: logical for local OMP adjustment
c adjustF   : subroutine for local optical model geometry adjustment
c factor    : Woods-Saxon multiplication factor
c
  200 if (ompadjustF(6)) then
        call adjustF(6,eopt,factor)
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
