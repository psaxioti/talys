      subroutine preeqtotal
c
c +---------------------------------------------------------------------
c | Author: Arjan Koning
c | Date  : August 14, 2004
c | Task  : Total pre-equilibrium cross sections
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      integer type,NL,p,nen,i,ppi,pnu,parity,J
      real    Elast,frac,norm
c
c ************************ Total pre-equilibrium ***********************
c
c parskip     : logical to skip outgoing particle
c Nlast,NL    : last discrete level
c parZ        : charge number of particle
c parN        : neutron number of particle
c Elast,frac  : help variables
c eoutdis     : outgoing energy of discrete state reaction
c Etop        : top of outgoing energy bin
c nendisc     : last discrete bin
c p           : particle number
c p0          : initial particle number
c maxpar      : maximal particle number   
c ebegin      : first energy point of energy grid 
c xssteptot   : preequilibrium cross section per particle type and stage
c xsstep      : preequilibrium cross section per particle type, stage
c               and outgoing energy
c deltaE      : energy bin around outgoing energies
c xspreeqtot  : preequilibrium cross section per particle type
c eend        : last energy point of energy grid 
c xspreeqtotps: preequilibrium cross section per particle type for
c               pickup and stripping
c xspreeqps   : preequilibrium cross section per particle type and
c               outgoing energy for pickup and stripping
c xspreeqtotki: preequilibrium cross section per particle type for
c               knockout and inelastic       
c xspreeqki   : preequilibrium cross section per particle type and
c               outgoing energy for knockout and inelastic              
c xspreeqsum  : total preequilibrium cross section summed over particles
c
c The pre-equilibrium spectra and spectra per exciton number are summed
c to total pre-equilibrium cross sections. Special care is taken for the
c continuum bin with the highest outgoing energy, i.e. the one that
c overlaps with the energy corresponding to the last discrete state.
c
      do 10 type=0,6
        if (parskip(type)) goto 10
        NL=Nlast(parZ(type),parN(type),0)
        Elast=eoutdis(type,NL)
        if (Elast.gt.0.) then
          frac=Etop(nendisc(type))-Elast
        else
          frac=0.
        endif
        do 20 p=p0,maxpar
          do 30 nen=ebegin(type),nendisc(type)
            xssteptot(type,p)=xssteptot(type,p)+xsstep(type,p,nen)*
     +        deltaE(nen)
   30     continue
          xssteptot(type,p)=xssteptot(type,p)-
     +      xsstep(type,p,nendisc(type))*frac     
          xspreeqtot(type)=xspreeqtot(type)+xssteptot(type,p)
   20   continue
        do 40 nen=ebegin(type),eend(type)
          xspreeqtotps(type)=xspreeqtotps(type)+xspreeqps(type,nen)*
     +      deltaE(nen)
          xspreeqtotki(type)=xspreeqtotki(type)+xspreeqki(type,nen)*
     +      deltaE(nen)
   40   continue
        xspreeqtotps(type)=xspreeqtotps(type)-
     +    xspreeqps(type,nendisc(type))*frac
        xspreeqtotki(type)=xspreeqtotki(type)-
     +    xspreeqki(type,nendisc(type))*frac
        xspreeqtot(type)=xspreeqtot(type)+xspreeqtotps(type)+
     +    xspreeqtotki(type)
        xspreeqsum=xspreeqsum+xspreeqtot(type)
   10 continue
c
c ************************* Unitarity condition ************************
c
c In line with unitarity, the summed direct + pre-equilibrium cross 
c section may not exceed the reaction cross section. In these cases,
c we normalize the results.
c
c xsflux        : cross section flux
c xsreacinc     : reaction cross section for incident channel  
c xsdirdiscsum  : total direct cross section
c xsgrsum       : sum over giant resonance cross sections
c norm,preeqnorm: preequilibrium normalization factor
c xspreeqdiscsum: total preequilibrium cross section for discrete states
c xspreeqdisc   : preequilibrium cross section for discrete state
c xspreeqdisctot: preequilibrium cross section summed over discrete
c                 states
c xspreeq       : preequilibrium cross section per particle type and 
c                 outgoing energy
c flag2comp     : flag for two-component pre-equilibrium model
c ppi           : proton particle number
c ppi0          : initial proton number
c pnu           : neutron particle number
c pnu0          : initial neutron number         
c xsstep2       : two-component preequilibrium cross section
c flagpespin    : flag for pre-equilibrium spin distribution or compound
c                 spin distribution for pre-equilibrium cross section  
c parity        : parity
c maxJph        : maximal spin for particle-hole states
c xspreeqJP     : preequilibrium cross section per particle type,
c                 outgoing energy, spin and parity
c
      xsflux=xsreacinc-xsdirdiscsum-xsgrsum
      preeqnorm=0.
      if (xspreeqsum+xspreeqdiscsum.gt.xsflux) then
        norm=xsflux/(xspreeqsum+xspreeqdiscsum)
        preeqnorm=norm
        xspreeqdiscsum=xspreeqdiscsum*norm
        xspreeqsum=xsflux-xspreeqdiscsum
        do 110 type=0,6
          if (parskip(type)) goto 110
          xspreeqtot(type)=xspreeqtot(type)*norm
          xspreeqtotps(type)=xspreeqtotps(type)*norm
          xspreeqtotki(type)=xspreeqtotki(type)*norm
          xspreeqdisctot(type)=xspreeqdisctot(type)*norm
          do 120 p=p0,maxpar
            xssteptot(type,p)=xssteptot(type,p)*norm
  120     continue
          do 130 i=0,Nlast(parZ(type),parN(type),0)
            xspreeqdisc(type,i)=xspreeqdisc(type,i)*norm
  130     continue
          do 140 nen=ebegin(type),eend(type)
            xspreeq(type,nen)=xspreeq(type,nen)*norm
            xspreeqps(type,nen)=xspreeqps(type,nen)*norm
            xspreeqki(type,nen)=xspreeqki(type,nen)*norm
            do 150 p=p0,maxpar
              xsstep(type,p,nen)=xsstep(type,p,nen)*norm
  150       continue
            if (flag2comp) then
              do 160 ppi=ppi0,maxpar
                do 160 pnu=pnu0,maxpar
                  xsstep2(type,ppi,pnu,nen)=xsstep2(type,ppi,pnu,nen)*
     +              norm
  160         continue
            endif
            if (flagpespin) then
              do 170 parity=-1,1,2
                do 170 J=0,maxJph
                  xspreeqJP(type,nen,J,parity)=
     +              xspreeqJP(type,nen,J,parity)*norm
  170         continue
            endif
  140     continue
  110   continue
      endif
c
c **** Add discrete pre-equilibrium contribution to discrete states ****
c
c xsdirdisctot: direct cross section summed over discrete states       
c xsdirdisc   : direct cross section for discrete state 
c
      xsdirdiscsum=xsdirdiscsum+xspreeqdiscsum
      do 210 type=0,6
        if (parskip(type)) goto 210
        xsdirdisctot(type)=xsdirdisctot(type)+xspreeqdisctot(type)
        do 220 i=0,Nlast(parZ(type),parN(type),0)  
          xsdirdisc(type,i)=xsdirdisc(type,i)+xspreeqdisc(type,i)
  220   continue
  210 continue
      return
      end
Copyright (C) 2004  A.J. Koning, S. Hilaire and M.C. Duijvestijn
