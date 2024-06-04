      subroutine gammanorm(Zcomp,Ncomp)
c
c +---------------------------------------------------------------------
c | Author: Stephane Hilaire and Arjan Koning 
c | Date  : September 13, 2005
c | Task  : Normalization of gamma ray strength functions 
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      integer          Zcomp,Ncomp,NL,nex,nexgam,J2b,J2e,J2,nexout,
     +                 Pprimebeg,Pprimeend,Irspin2beg,Irspin2end,Pprime,
     +                 Irspin2,l2beg,l2end,l2,l,modl,irad,nen
      real             Sn,Exgam(0:10*numex),dEx,gamsum,Exmid,Egamma,
     +                 Exmin,Explus,dE1,dE2,Rspin,Tgamma,fstrength,
     +                 gammaxs
      double precision density,rho1,rho2,rho3,r1log,r2log,r3log,rho
c
c We normalize the gamma-ray strength functions by imposing the 
c condition that the transmission coefficients integrated from zero up 
c to neutron separation energy are equal to the ratio of the 
c experimental mean gamma width and mean level spacing for s-wave 
c neutrons.
c
c *************** Initialization of excitation energies ****************
c
c Zcomp       : charge number index for compound nucleus
c Ncomp       : neutron number index for compound nucleus
c swaveth     : theoretical strength function for s-wave
c gnorm       : gamma normalization factor
c Sn          : neutron separation energy
c S           : separation energy per particle
c Nlast,NL    : last discrete level
c edis        : energy of level
c nexgam,numex: maximum excitation energy bin for gamma normalization
c Exgam       : excitation energy for gamma normalization
c dEx         : excitation energy bin 
c
c If the gamma normalization factor has been given in the input, we
c do not need to calculate it.
c
      swaveth=0.
      if (gnorm.ne.-1..or.gnorm.eq.1.) goto 210
      Sn=S(Zcomp,Ncomp,1)
      NL=Nlast(Zcomp,Ncomp,0)
      do 10 nex=0,NL
        if (edis(Zcomp,Ncomp,nex).gt.Sn) then
          nexgam=nex-1
          goto 30
        endif
        Exgam(nex)=edis(Zcomp,Ncomp,nex)
   10 continue
      nexgam=10*numex
      dEx=(Sn-Exgam(NL))/(nexgam-NL)
      do 20 nex=NL+1,nexgam
        Exgam(nex)=Exgam(NL)+(nex-NL)*dEx
   20 continue
c
c ********************** Loop over quantum numbers *********************
c
c J2b       : 2 * start of J summation
c targetspin: spin of target
c parspin   : spin of incident particle
c J2e       : 2 * end of J summation
c gamsum    : sum over gamma-ray transmission coefficients
c
c Specify the summation boundaries for the neutron channel.
c
   30 J2b=int(abs(2.*(targetspin-parspin(1))))
      J2e=int(2.*(targetspin+parspin(1)))
      gamsum=0.
c
c 110: Sum over total angular momentum J (J2) of compound nucleus
c
c J2: 2 * J
c
      do 110 J2=J2b,J2e,2
c
c 120: Sum over outgoing excitation energies
c
c Exmid : help variable
c Egamma: gamma energy
c
        do 120 nexout=0,nexgam
          Exmid=Exgam(nexout)
          Egamma=Sn-Exmid
          if (Egamma.le.0.) goto 120
c
c Set begin and end energies for level density integration
c
c Exmin,Explus: help variables
c dE1,dE2     : help variables
c
          if (nexout.gt.NL) then
            Exmin=max(Exmid-0.5*dEx,0.)
            Explus=Exmid+0.5*dEx
            dE1=Exmid-Exmin
            dE2=Explus-Exmid
          endif
c
c Initialization of summations
c
c Pprimebeg : start of residual parity summation
c parlev    : parity of level
c Pprimeend : end of residual parity summation
c Irspin2beg: 2 * start of residual spin summation
c jdis      : spin of level
c Irspin2end: 2 * end of residual spin summation
c rho       : integrated level density
c gammax    : number of l-values for gamma multipolarity
c numJ      : maximal J-value
c
c For discrete states, the begin and end points of the residual
c spin/parity summation are both set equal to the residual discrete
c level spin/parity.
c    
          if (nexout.le.NL) then
            Pprimebeg=parlev(Zcomp,Ncomp,nexout)
            Pprimeend=Pprimebeg  
            Irspin2beg=int(2.*jdis(Zcomp,Ncomp,nexout))
            Irspin2end=Irspin2beg
            rho=1.
          else
c
c For the continuum, the begin and end points of the residual
c spin/parity summation are set to the maximally accessible values.
c      
            Pprimebeg=-1
            Pprimeend=1   
            Irspin2beg=mod(J2,2)
            Irspin2end=J2+2*gammax
            Irspin2end=min(Irspin2end,2*numJ)
          endif
c
c 130: Sum over residual parity
c
c pardif : difference between residual and compound nucleus parity
c targetP: parity of target
c
          do 130 Pprime=Pprimebeg,Pprimeend,2 
            pardif=abs(targetP-Pprime)/2
c
c 140: Sum over residual spin
c
c Irspin2       : 2 * residual spin
c Rspin         : residual spin
c rho1,rho2,rho3: help variables
c density       : level density
c ldmodel       : level density model
c r1log,..      : help variables   
c
c rho is the level density integrated over the bin (Explus - Exmin).
c Note that we use logarithmic integration for more precision.
c
            do 140 Irspin2=Irspin2beg,Irspin2end,2
              Rspin=0.5*Irspin2
              if (nexout.gt.NL) then
                rho1=real(density(Zcomp,Ncomp,Exmin,Rspin,Pprime,0,
     +            ldmodel))+1.e-30
                rho2=real(density(Zcomp,Ncomp,Exmid,Rspin,Pprime,0,
     +            ldmodel))
                rho3=real(density(Zcomp,Ncomp,Explus,Rspin,Pprime,0,
     +            ldmodel))+1.e-30
                r1log=log(rho1)
                r2log=log(rho2)
                r3log=log(rho3)
                if (r2log.ne.r1log.and.r2log.ne.r3log) then
                  rho=(rho1-rho2)/(r1log-r2log)*dE1
     +              +(rho2-rho3)/(r2log-r3log)*dE2
                else
                  rho=rho2*(dE1+dE2)
                endif
              endif
c
c 150: Sum over l of outgoing channel
c
c l2beg: 2 * start of l-summation
c l2end: 2 * end l-summation
c
              l2beg=max(abs(J2-Irspin2),2)
              l2end=min(J2+Irspin2,2*gammax)
              do 150 l2=l2beg,l2end,2
c
c modl: help variable
c irad: variable to indicate M(=0) or E(=1) radiation
c
c Multipole radiation selection rules
c (irad=0: M-transition, irad=1: E-transition)
c
                l=l2/2
                modl=mod(l,2)
                if (pardif.eq.modl) then
                  irad=1
                else
                  irad=0
                endif 
c
c Create sum for normalization
c
c Tgamma   : gamma transmission coefficient
c twopi    : 2.*pi
c fstrength: gamma ray strength function
c
                Tgamma=twopi*(Egamma**(l2+1))*
     +            fstrength(Zcomp,Ncomp,Egamma,irad,l)
                gamsum=gamsum+rho*Tgamma
  150         continue
  140       continue
  130     continue
  120   continue
  110 continue
c
c ************** Normalization of transmission coefficients ***********
c
c gamgam: total radiative width
c Dtheo : theoretical s-wave resonance spacing
c ebegin: first energy point of energy grid
c eend  : last energy point of energy grid
c egrid : outgoing energy grid
c lmax  : maximal l-value for transmission coefficients
c Tjl   : transmission coefficients as a function of particle type, 
c         energy, spin and l-value
c
      if (gamsum.ne.0.) then
        swaveth=gamgam(Zcomp,Ncomp)/Dtheo(Zcomp,Ncomp)
        gnorm=twopi*swaveth/gamsum
      else
        gnorm=1.
      endif
  210 do 220 nen=ebegin(0),eend(0)
        Egamma=egrid(nen)
        lmax(0,nen)=gammax
        do 230 l=1,gammax
          do 230 irad=0,1
            Tgamma=twopi*(Egamma**(2*l+1))*
     +        fstrength(Zcomp,Ncomp,Egamma,irad,l)
            Tjl(0,nen,irad,l)=Tgamma*gnorm
  230   continue
c
c Photo-absorption cross section
c
c xsreac : reaction cross section   
c gammaxs: function for gamma ray cross sections 
c
        xsreac(0,nen)=gammaxs(Zcomp,Ncomp,Egamma)
  220 continue
      return
      end
Copyright (C) 2004  A.J. Koning, S. Hilaire and M.C. Duijvestijn
