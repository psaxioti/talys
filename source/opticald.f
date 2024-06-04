      subroutine opticald(Zix,Nix,eopt)
c
c +---------------------------------------------------------------------
c | Author: Arjan Koning
c | Date  : May 11, 2011
c | Task  : Optical potential for deuterons
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      character*72 optmodfile
      integer      Zix,Nix,nen,i,Z,A,N
      real         eopt,elow,eup,eint,vloc(19),e,vn,rvn,avn,wn,rwn,awn,
     +             vdn,rvdn,avdn,wdn,vp,rvp,avp,wp,rwp,awp,vdp,rvdp,
     +             avdp,wdp,vson,rvson,avson,wson,vsop,rvsop,avsop,wsop,
     +             A13,mu,summu,beta,asym,factor(12)
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
c 1. In case of an optical model file, we interpolate between the
c    tabulated values
c
      optmodfile='                                                     '
      if (Zix.le.numZph.and.Nix.le.numNph) optmodfile=optmod(Zix,Nix,3)
      if (optmodfile(1:1).ne.' ') then
        if (eopt.lt.eomp(Zix,Nix,3,1).or.
     +    eopt.gt.eomp(Zix,Nix,3,omplines(Zix,Nix,3))) goto 100
        do 20 nen=1,omplines(Zix,Nix,3)-1
          elow=eomp(Zix,Nix,3,nen)
          eup=eomp(Zix,Nix,3,nen+1)
          if (elow.le.eopt.and.eopt.le.eup) then
            eint=(eopt-elow)/(eup-elow)
            do 30 i=1,19
              vloc(i)=vomp(Zix,Nix,3,nen,i)+
     +          eint*(vomp(Zix,Nix,3,nen+1,i)-vomp(Zix,Nix,3,nen,i))
   30       continue
            v=v1adjust(3)*vloc(1)
            rv=rvadjust(3)*vloc(2)
            av=avadjust(3)*vloc(3)
            w=w1adjust(3)*vloc(4)
            rw=rwadjust(3)*vloc(5)
            aw=awadjust(3)*vloc(6)
            vd=d1adjust(3)*vloc(7)
            rvd=rvdadjust(3)*vloc(8)
            avd=avdadjust(3)*vloc(9)
            wd=d1adjust(3)*vloc(10)
            rwd=rwdadjust(3)*vloc(11)
            awd=awdadjust(3)*vloc(12)
            vso=vso1adjust(3)*vloc(13)
            rvso=rvsoadjust(3)*vloc(14)
            avso=avsoadjust(3)*vloc(15)
            wso=wso1adjust(3)*vloc(16)
            rwso=rwsoadjust(3)*vloc(17)
            awso=awsoadjust(3)*vloc(18)
            rc=rcadjust(3)*vloc(19)
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
c 1988, p. 103) to make a deuteron potential out of the proton and
c neutron potential. The simplified formula is:
c                V(d,E)=V(n,E/2)+V(p,E/2)
c and similarly for Wd and W. VSO and Wso are as in the nucleon
c potentials.
c We take a similar weighting for the geometry parameters.
c
c 1. Neutron potential for E/2.
c
c e        : incident energy (help variable)
c opticaln : optical potential for neutrons
c vn,wn....: optical model parameters for neutrons
c
  100 e=0.5*eopt
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
c
c 2. Proton potential for E/2.
c
c opticalp : optical potential for protons
c vp,wp....: optical model parameters for protons
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
c
c 3. Another 2 calls to opticaln for eopt=E to determine Vso and Wso.
c
c vson,wson: Vso and Wso for neutrons
c vsop,wsop: Vso and Wso for protons
c
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
      v=v1adjust(3)*(vn+vp)
      rv=rvadjust(3)*0.5*(rvn+rvp)
      av=avadjust(3)*0.5*(avn+avp)
      w=w1adjust(3)*(wn+wp)
      rw=rwadjust(3)*0.5*(rwn+rwp)
      aw=awadjust(3)*0.5*(awn+awp)
      vd=d1adjust(3)*(vdn+vdp)
      rvd=rvdadjust(3)*0.5*(rvdn+rvdp)
      avd=avdadjust(3)*0.5*(avdn+avdp)
      wd=d1adjust(3)*(wdn+wdp)
      rwd=rwdadjust(3)*0.5*(rvdn+rvdp)
      awd=awdadjust(3)*0.5*(avdn+avdp)
      vso=vso1adjust(3)*0.5*(vson+vsop)
      rvso=rvsoadjust(3)*0.5*(rvson+rvsop)
      avso=avsoadjust(3)*0.5*(avson+avsop)
      wso=wso1adjust(3)*0.5*(wson+wsop)
      rwso=rwsoadjust(3)*0.5*(rvson+rvsop)
      awso=awsoadjust(3)*0.5*(avson+avsop)
c
c Alternative options for deuteron OMP
c
c deuteronomp=1 --> Watanabe potential
c deuteronomp=2 --> Daehnick potential
c deuteronomp=3 --> Bojowald potential
c deuteronomp=4 --> Han potential
c deuteronomp=5 --> Haixia An potential
c
c Overwrite some of the previous values.
c
c ZZ,Z    : charge number of residual nucleus
c AA,A    : mass number of residual nucleus
c onethird: 1/3
c magic   : magic numbers
c beta,mu : help variables
c
      if (deuteronomp.ge.2) then
        Z=ZZ(Zix,Nix,0)
        A=AA(Zix,Nix,0)
        N=A-Z
        A13=real(A)**onethird
        asym=(N-Z)/real(A)
c
c deuteronomp=2 --> Global OMP of W.W. Daehnick, J.D. Childs, Z. Vrcelj,
c Phys. Rev C21, 2253 (1980): relativistic case
c
        if (deuteronomp.eq.2) then
          v=88.0+0.88*Z/A13-0.283*eopt
          rv=1.17
          av=0.717+0.0012*eopt
          beta=-(0.01*eopt)**2
          w=(12+0.031*eopt)*(1.-exp(beta))
          rw=1.376-0.01*sqrt(eopt)
          summu=0.
          do i=1,8
            mu=(0.5*(magic(i)-N))**2
            summu=summu+exp(-mu)
          enddo
          aw=0.52+0.07*A13-0.04*summu
          vd=0.
          rvd=rw
          avd=aw
          wd=(12.+0.031*eopt)*exp(beta)
          rwd=rvd
          awd=avd
          vso=5.0
          rvso=1.04
          avso=0.60
          wso=0.37*A13-0.03*eopt
          rwso=0.80
          awso=0.25
          rc=1.30
        endif
c
c deuteronomp=3 --> Global OMP of J. Bojowald et al,
c Phys. Rev. C38, 1153 (1988).
c
        if (deuteronomp.eq.3) then
          v=81.32-0.24*eopt+1.43*Z/A13
          rv=1.18
          av=0.636+0.035*A13
          w=max(0.,0.132*(eopt-45.))
          rw=1.27
          aw=0.768+0.021*A13
          vd=0.
          rvd=rw
          avd=aw
          wd=max(0.,7.80+1.04*A13-0.712*w)
          rwd=rvd
          awd=avd
          vso=6.
          rvso=0.78+0.038*A13
          avso=rvso
          wso=0.
          rwso=rvso
          awso=avso
          rc=1.30
        endif
c
c deuteronomp=4 --> Global OMP of Y. Han, Y. Shi and Q. Shen,
c Phys. Rev. C74, 044615 (2006).
c
        if (deuteronomp.eq.4) then
          v=82.18-0.148*eopt-0.000886*eopt*eopt-34.811*asym+1.058*Z/A13
          rv=1.174
          av=0.809
          w=max(0.,-4.916+0.0555*eopt+4.42e-5*eopt*eopt+35.*asym)
          rw=1.563
          aw=0.700+0.045*A13
          vd=0.
          rvd=1.328
          avd=0.465+0.045*A13
          wd=20.968-0.0794*eopt-43.398*asym
          rwd=rvd
          awd=avd
          vso=3.703
          rvso=1.234
          avso=0.813
          wso=-0.206
          rwso=rvso
          awso=avso
          rc=1.698
        endif
c
c deuteronomp=5 --> Global OMP of Haixia An and Chonghai Cai,
c Phys. Rev. C73, 054605 (2006).
c
        if (deuteronomp.eq.5) then
          v=91.85-0.249*eopt-0.000116*eopt*eopt+0.642*Z/A13
          rv=1.152-0.00776/A13
          av=0.719+0.0126*A13
          w=1.104+0.0622*eopt
          rw=1.305+0.0997/A13
          aw=0.855-0.100*A13
          vd=0.
          rvd=1.334+0.152/A13
          avd=0.531+0.062*A13
          wd=10.83-0.0306*eopt
          rwd=rvd
          awd=avd
          vso=3.557
          rvso=0.972
          avso=1.011
          wso=0.
          rwso=rvso
          awso=avso
          rc=1.303
        endif
c
c Possible adjustment of parameters
c
        v=v1adjust(3)*v
        rv=rvadjust(3)*rv
        av=avadjust(3)*av
        w=w1adjust(3)*w
        rw=rwadjust(3)*rw
        aw=awadjust(3)*aw
        vd=d1adjust(3)*vd
        rvd=rvdadjust(3)*rvd
        avd=avdadjust(3)*avd
        wd=d1adjust(3)*wd
        rwd=rwdadjust(3)*rwd
        awd=awdadjust(3)*awd
        vso=vso1adjust(3)*vso
        rvso=rvsoadjust(3)*rvso
        avso=avsoadjust(3)*avso
        wso=wso1adjust(3)*wso
        rwso=rwsoadjust(3)*rwso
        awso=awsoadjust(3)*awso
      endif
c
c Possible additional energy-dependent adjustment of the geometry
c
c ompadjustF: logical for local OMP adjustment
c adjustF   : subroutine for local optical model geometry adjustment
c factor    : Woods-Saxon multiplication factor
c
  200 if (ompadjustF(3)) then
        call adjustF(3,eopt,factor)
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
