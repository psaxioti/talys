      subroutine multipreeq(Zcomp,Ncomp,nex)
c
c +---------------------------------------------------------------------
c | Author: Arjan Koning 
c | Date  : August 1, 2008
c | Task  : Multiple preequilibrium model
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      logical          surfwell
      integer          Zcomp,Ncomp,nex,ip,ih,ZNcomp,itype,p,h,type,Zix,
     +                 Nix,nen,nexout,parity,J
      real             sumfeed,gs,ignatyuk,Rfactor,omegaph,phdens,dEx,
     +                 Eex,Eo,term,omegap1h,EoplusS,omega1p,proba,
     +                 Tswave,Pescape,Exm,Exmin,factor
      double precision summpe,feedph,sumph,sumtype
c
c mpreeqmode=1: Multiple emission exciton model
c mpreeqmode=2: Transmission coefficients for particle emission
c
c Multiple preequilibrium emission model 2 is adopted from Chadwick
c and Young, Phys. Rev. C50 (1994) p. 996.
c
c ************ Check presence of multiple pre-equilibrium **************
c
c Zcomp  : charge number index for compound nucleus
c Ncomp  : neutron number index for compound nucleus  
c sumfeed: help variable
c ip     : particle number of mother nucleus
c maxpar : maximal particle number
c ih     : hole number of mother nucleus
c xspopph: population cross section per particle-hole configuration
c
c There must be excited particles and holes present for multiple
c pre-equilibrium emission to occur.
c
      if (Zcomp.eq.0.and.Ncomp.eq.0) return
      sumfeed=0.
      do 10 ip=1,maxpar
        do 10 ih=1,maxpar
          sumfeed=sumfeed+xspopph(Zcomp,Ncomp,nex,ip,ih)
   10 continue
      if (sumfeed.le.1.e-10) return
c
c ***************** Multiple preequilibrium emission *******************
c
c Ecomp         : total energy of composite system
c Exinc         : excitation energy of entrance bin
c gs,g          : single-particle level density parameter
c flaggshell    : flag for energy dependence of single particle level
c                 density parameter g
c ignatyuk      : function for energy dependent level density 
c                 parameter a
c alev          : level density parameter
c surfwell      : flag for surface effects in finite well 
c Rfactor,Rblann: Blann's factor
c ZNcomp        : help variable
c
      Ecomp=Exinc
      gs=g(Zcomp,Ncomp)
      if (flaggshell) gs=gs*ignatyuk(Zcomp,Ncomp,Ecomp,0)/
     +  alev(Zcomp,Ncomp)
      surfwell=.false.
      Rfactor=0.5
      ZNcomp=Zcomp+Ncomp
c
c For secondary pre-equilibrium emission, determine the type of the
c first emitted particle.
c
c itype : help variable
c summpe: multiple preequilibrium emission for excitation energy bin
c
      if (ZNcomp.eq.1) then
        if (Ncomp.eq.1) then
          itype=1
        else
          itype=2
        endif
      endif
      summpe=0.
c
c 110 and 120: Loops over all possible particle-hole excitations of
c the mother bin.
c
c feedph      : feeding term from previous particle-hole calculation
c sumph       : multiple preequilibrium flux for particle-hole pair
c mpreeqmode  : designator for multiple pre-equilibrium model
c omegaph     : particle-hole state density for compound system
c phdens      : particle-hole state density
c Efermi      : depth of Fermi well
c p0          : initial particle number
c h0          : initial hole number
c p           : particle number
c h           : hole number
c n           : exciton number
c emissionrate: subroutine for emission rate
c lifetime    : subroutine for calculation of lifetime of exciton state
c
      do 110 ip=1,maxpar
        do 120 ih=1,maxpar
          feedph=xspopph(Zcomp,Ncomp,nex,ip,ih)
          if (feedph.le.1.e-10) goto 120
          sumph=0.
          if (mpreeqmode.eq.2) then
            omegaph=phdens(Zcomp,Ncomp,ip,ih,gs,Exinc,Efermi,surfwell)
            omegaph=max(omegaph,1.)
          else
            p0=ip
            h0=ih
            do 130 p=p0,maxpar
              h=h0+p-p0
              if (h.gt.maxpar) goto 130
              call emissionrate(Zcomp,Ncomp,p,h)
              call lifetime(Zcomp,Ncomp,p,h)
  130       continue
          endif
c
c 140: Loop over neutron and proton type of next emitted particle
c
c parskip    : logical to skip outgoing particle
c Zindex,Zix : charge number index for residual nucleus
c Nindex,Nix : neutron number index for residual nucleus  
c numZph     : maximal number of protons away from the initial
c              compound nucleus for multiple pre-equilibrium emission
c numNph     : maximal number of neutrons away from the initial
c              compound nucleus for multiple pre-equilibrium emission   
c sumtype    : help variable
c deltaEx,dEx: excitation energy bin for population arrays
c
          do 140 type=1,2
            if (parskip(type)) goto 140
            Zix=Zindex(Zcomp,Ncomp,type)
            Nix=Nindex(Zcomp,Ncomp,type)           
            if (Zix.gt.numZph.or.Nix.gt.numNph) goto 140
            if (ZNcomp.eq.1) Rfactor=Rblann(itype,type,ip)
            sumtype=0.
            dEx=deltaEx(Zix,Nix)
c
c 150: Loop over residual excitation energy bins
c
c Nlast     : last discrete level
c nexmax    : maximum excitation energy bin for residual nucleus
c Ex,Eex    : excitation energy
c Eo        : outgoing energy
c S         : separation energy per particle
c locate    : subroutine to find value in ordered table
c egrid     : energies of basic energy grid in MeV
c ebegin    : first energy point of energy grid
c eend      : last energy point for energy grid of each particle
c term      : help variable
c omegap1h  : particle-hole state density for residual system
c EoplusS   : outgoing energy+S 
c omega1p   : particle-hole state density for continuum particle
c proba     : probability to find continuum particle
c Tswave    : transmission coefficient for s-wave
c Tjl       : transmission coefficients as a function of particle 
c             type, energy, spin and l-value     
c Pescape   : escape probability
c Exm       : maximal attainable energy
c dExinc    : excitation energy bin for mother nucleus
c Exmin     : lower bound of energy bin
c factor    : help variable
c tauexc    : mean lifetime
c wemission : emission rate
c mcontrib  : contribution to emission spectrum
c mpecontrib: contribution to multiple pre-equilibrium emission spectrum
c xspopex   : population cross section summed over spin and parity
c
            do 150 nexout=Nlast(Zix,Nix,0)+1,nexmax(type)
              Eex=Ex(Zix,Nix,nexout)
              Eo=Exinc-Eex-S(Zcomp,Ncomp,type)
              call locate(egrid,ebegin(type),eend(type),Eo,nen)
              term=0.
              if (mpreeqmode.eq.2) then
                gs=g(Zix,Nix)
                if (flaggshell) gs=g(Zix,Nix)*
     +            ignatyuk(Zix,Nix,Eex,0)/alev(Zix,Nix)
                omegap1h=phdens(Zix,Nix,ip-1,ih,gs,Eex,Efermi,surfwell)
                EoplusS=Exinc-Eex
                omega1p=phdens(Zix,Nix,1,0,gs,EoplusS,Efermi,surfwell)
                proba=omega1p*omegap1h/omegaph/ip*Rfactor
                Tswave=Tjl(type,nen,1,0)
                Pescape=proba*Tswave
                if (nexout.eq.nexmax(type)) then
                  Exm=Exinc+0.5*dExinc-S(Zcomp,Ncomp,type)
                  Exmin=Ex(Zix,Nix,nexout)-0.5*dEx
                  dEx=Exm-Exmin
                endif
                term=feedph*Pescape*dEx
                xspopph(Zix,Nix,nexout,ip-1,ih)=
     +            xspopph(Zix,Nix,nexout,ip-1,ih)+term
              else
                do 160 p=p0,maxpar
                  h=h0+p-p0
                  if (h.gt.maxpar) goto 160
                  factor=feedph*tauexc(p,h)*wemission(type,p,h,nen)
     +              *Rfactor*dEx
                  xspopph(Zix,Nix,nexout,p-1,h)=
     +              xspopph(Zix,Nix,nexout,p-1,h)+factor
                  term=term+factor
  160           continue
              endif
              sumtype=sumtype+term
              mcontrib(type,nex,nexout)=mcontrib(type,nex,nexout)+term
              mpecontrib(type,nex,nexout)=mpecontrib(type,nex,nexout)+
     +          term
              xspopex(Zix,Nix,nexout)=xspopex(Zix,Nix,nexout)+term
c
c 170: Adopt spin distribution from primary pre-equilibrium.
c
c parity: parity
c maxJ  : maximal J-value
c xspop : population cross section
c RnJ   : spin distribution for particle-hole states
c RnJsum: (2J+1)*sum over spin distributions
c
              do 170 parity=-1,1,2
                do 170 J=0,maxJ(Zix,Nix,nexout)
                  xspop(Zix,Nix,nexout,J,parity)=
     +              xspop(Zix,Nix,nexout,J,parity)+
     +              0.5*(2*J+1)*RnJ(2,J)/RnJsum(2)*term
  170         continue
  150       continue
c
c Add total multiple pre-equilibrium contributions
c
c xspopnuc : population cross section per nucleus 
c xspartial: emitted cross section flux per energy bin 
c xsfeed   : cross section from compound to residual nucleus
c xsmpe    : multiple-preequilibrium cross section per energy bin
c xsmpetot : total multiple-preequilibrium cross section 
c mulpreZN : logical for multiple pre-equilibrium per nucleus
c
            sumph=sumph+sumtype
            summpe=summpe+sumtype
            xspopnuc(Zix,Nix)=xspopnuc(Zix,Nix)+sumtype
            xspartial(type,nex)=xspartial(type,nex)+sumtype
            xsfeed(Zcomp,Ncomp,type)=xsfeed(Zcomp,Ncomp,type)+sumtype
            xsmpe(type,nex)=xsmpe(type,nex)+sumtype
            xsmpetot(type)=xsmpetot(type)+sumtype
            if (sumtype.ne.0.) mulpreZN(Zix,Nix)=.true.
  140     continue
c
c Flux that is not emitted during a particular stage, is always 
c transferred to the next stage.
c
          if (mpreeqmode.eq.2) then
            if (ip.le.maxpar-1.and.ih.le.maxpar-1)
     +        xspopph(Zcomp,Ncomp,nex,ip+1,ih+1)=
     +        xspopph(Zcomp,Ncomp,nex,ip+1,ih+1)+feedph-sumph
          endif
  120   continue
  110 continue
c
c ************************ Normalization *******************************
c
c Dmulti: depletion factor for multiple preequilibrium
c
      Dmulti(nex)=summpe/xspopex(Zcomp,Ncomp,nex)
      xspopex(Zcomp,Ncomp,nex)=xspopex(Zcomp,Ncomp,nex)-summpe
      return
      end
Copyright (C) 2004  A.J. Koning, S. Hilaire and M.C. Duijvestijn
