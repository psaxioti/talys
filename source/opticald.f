      subroutine opticald(Zix,Nix,eopt)
c
c +---------------------------------------------------------------------
c | Author: Arjan Koning 
c | Date  : July 7, 2004
c | Task  : Optical potential for deuterons
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      character*72 optmodfile
      integer      Zix,Nix,nen,i
      real         eopt,elow,eup,eint,vloc(19),vn,rvn,avn,wn,rwn,awn,
     +             vdn,rvdn,avdn,wdn,rwdn,awdn,vp,rvp,avp,wp,rwp,awp,
     +             vdp,rvdp,avdp,wdp,rwdp,awdp,vson,rvson,avson,wson,
     +             rwson,awson,vsop,rvsop,avsop,wsop,rwsop,awsop
c
c ******************* Optical model input file *************************
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
c    tabulated values
c
      optmodfile=optmod(Zix,Nix,3)
      if (optmodfile(1:1).ne.' ') then
        if (eopt.lt.eomp(3,1).or.eopt.gt.eomp(3,omplines(3))) goto 100
        do 20 nen=1,omplines(3)-1
          elow=eomp(3,nen)
          eup=eomp(3,nen+1)
          if (elow.le.eopt.and.eopt.le.eup) then
            eint=(eopt-elow)/(eup-elow)
            do 30 i=1,19
              vloc(i)=vomp(3,nen,i)+eint*(vomp(3,nen+1,i)-vomp(3,nen,i))
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
c 1988, p. 103) to make a deuteron potential out of the proton and 
c neutron potential. The simplified formula is:
c                V(d,E)=V(n,E/2)+V(p,E/2) 
c and similarly for Wd and W. VSO and Wso are as in the nucleon 
c potentials.
c We take a similar weighting for the geometry parameters.
c
c 1. Neutron potential for E/2.
c
c eopt     : incident energy (help variable)
c opticaln : optical potential for neutrons 
c vn,wn....: optical model parameters for neutrons
c
  100 eopt=0.5*eopt
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
c 2. Proton potential for E/2.
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
      rwdp=rwd
      awdp=awd
c
c 3. Another 2 calls to opticaln for eopt=E to determine Vso and Wso.
c
c vson,wson: Vso and Wso for neutrons
c vsop,wsop: Vso and Wso for protons
c
      eopt=eopt*2.
      call opticaln(Zix,Nix,eopt)
      vson=vso
      rvson=rvso
      avson=avso
      wson=wso
      rwson=rwso
      awson=awso
      call opticalp(Zix,Nix,eopt)
      vsop=vso
      rvsop=rvso
      avsop=avso
      wsop=wso
      rwsop=rwso
      awsop=awso
c
c 4. Final potential depths: construct V, W, Wd, Vso and Wso.
c
      v=vn+vp
      rv=0.5*(rvn+rvp)
      av=0.5*(avn+avp)
      w=wn+wp
      rw=0.5*(rwn+rwp)
      aw=0.5*(awn+awp)
      vd=vdn+vdp
      rvd=0.5*(rvdn+rvdp)
      avd=0.5*(avdn+avdp)
      wd=wdn+wdp
      rwd=0.5*(rwdn+rwdp)
      awd=0.5*(awdn+awdp)
      vso=0.5*(vson+vsop)
      rvso=0.5*(rvson+rvsop)
      avso=0.5*(avson+avsop)
      wso=0.5*(wson+wsop)
      rwso=0.5*(rwson+rwsop)
      awso=0.5*(awson+awsop)
      return
      end
Copyright (C) 2004  A.J. Koning, S. Hilaire and M.C. Duijvestijn
