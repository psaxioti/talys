      subroutine normalization
c
c +---------------------------------------------------------------------
c | Author: Arjan Koning 
c | Date  : September 23, 2009
c | Task  : Normalize cross sections to experimental or evaluated data
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      integer mt,nen,type,idc,ident,npart,ia,ih,it,id,ip,in,i1,i2,nex,
     +        type2
      real    Efac,xsadd,xsdif(nummt),add,ratio
c
c ************************ Adjustment factors **************************
c
c Set incident energy dependent adjustment factors (purely for
c fitting purposes).
c
c Nrescue: number of energies for adjustment factors
c Crescue: adjustment factor for this incident energy
c Einc   : incident energy in MeV
c Erescue: energy grid for adjustment factors
c frescue: adjustment factor
c
      do 10 mt=1,nummt
        if (Nrescue(mt).eq.0) goto 10
        Crescue(mt)=1.
        if (Einc.le.Erescue(mt,1)) then
          Crescue(mt)=frescue(mt,1)
          goto 10
        endif
        if (Einc.ge.Erescue(mt,Nrescue(mt))) then
          Crescue(mt)=frescue(mt,Nrescue(mt))
          goto 10
        endif
        do  20 nen=1,Nrescue(mt)-1
          if (Einc.gt.Erescue(mt,nen).and.Einc.le.Erescue(mt,nen+1)) 
     +      then
            Efac=(Einc-Erescue(mt,nen))/
     +        (Erescue(mt,nen+1)-Erescue(mt,nen))
            Crescue(mt)=frescue(mt,nen)+
     +        Efac*(frescue(mt,nen+1)-frescue(mt,nen))
            goto 10
          endif
   20   continue
   10 continue
c
c *************************** Normalization ****************************
c
c xsadd      : total difference in cross section after normalization
c xsnonel    : non-elastic cross section
c idc,ident  : help variables
c idnum      : counter for exclusive channel
c idchannel  : identifier for exclusive channel
c xsdif      : difference in cross section after normalization
c xschannel  : channel cross section
c flagfission: flag for fission
c xsfistot   : total fission cross section
c
c Application of normalization from rescue files. For each partial 
c cross section the applied difference is stored, so that later it can 
c be redistributed over correlated channels.
c
      xsadd=0.
      do 110 type=0,6
        mt=101+type
        if (type.eq.0) mt=102
        if (type.eq.1) mt=4
        if (type.eq.0) ident=0
        if (type.eq.1) ident=100000
        if (type.eq.2) ident=10000
        if (type.eq.3) ident=1000
        if (type.eq.4) ident=100
        if (type.eq.5) ident=10
        if (type.eq.6) ident=1
        do 120 idc=0,idnum
          if (idchannel(idc).eq.ident.and.Crescue(mt).ne.1.and.
     +      Crescue(mt).ne.0.) then
            xsdif(mt)=xschannel(idc)*(1./Crescue(mt)-1.)
            xsadd=xsadd+xsdif(mt)
            goto 110
          endif
  120   continue
  110 continue
      do 130 idc=0,idnum
        if (idchannel(idc).eq.200000.and.Crescue(16).ne.1.and.
     +    Crescue(16).ne.0.) then
          xsdif(16)=xschannel(idc)*(1./Crescue(16)-1.)
          xsadd=xsadd+xsdif(16)
          goto 140
        endif
  130 continue
  140 if (flagfission.and.Crescue(18).ne.1..and.Crescue(18).ne.0.) then
        xsdif(18)=xsfistot*(1./Crescue(18)-1.)
        xsadd=xsadd+xsdif(18)
      endif
      if (Crescue(1).ne.1..and.Crescue(1).ne.0.)
     +  xsdif(1)=xstotinc*(1./Crescue(1)-1.)
      if (Crescue(2).ne.1..and.Crescue(2).ne.0.)
     +  xsdif(2)=xselastot*(1./Crescue(2)-1.)
      if (Crescue(3).ne.1..and.Crescue(3).ne.0.)
     +  xsdif(3)=xsnonel*(1./Crescue(3)-1.)
c
c Normalization of channel cross sections and/or spectra, by the 
c rescue factors.
c
c npart        : number of particles in outgoing channel
c maxchannel   : maximal number of outgoing particles in individual 
c                channel description (e.g. this is 3 for (n,2np))
c numia,....   : maximal number of ejectile in channel description  
c factor       : factor to determine redistribution of cross section
c add          : help variable
c ratio        : adjustment ratio
c xsgamchannel : gamma channel cross section
c xsgamdischan : discrete gamma channel cross section
c xschannelsp  : channel cross section spectra
c xsdisc       : total cross section for discrete state
c xscompdisc   : compound cross section for discrete state     
c xsdirdisc    : direct cross section for discrete state 
c xsexclcont   : exclusive single channel cross section for continuum
c xsngn        : total (projectile,gamma-ejectile) cross section
c xsexclusive  : exclusive single channel cross section
c xsdisctot    : total cross section summed over discrete states
c xsdirdisctot : direct cross section summed over discrete states
c xscompdisctot: compound cross section summed over discrete states  
c xscompcont   : compound cross section for continuum     
c xsdircont    : direct cross section for continuum
c xsconttot    : total cross section for continuum
c xsdirect     : total direct cross section
c xsbinary     : cross section from initial compound to residual nucleus
c xsparticle   : total particle production cross section
c multiplicity : particle multiplicity
c xsfischannel : fission channel cross section
c
      do 210 npart=0,maxchannel
      do 210 ia=0,numia 
      do 210 ih=0,numih
      do 210 it=0,numit
      do 210 id=0,numid
      do 210 ip=0,numip
      do 210 in=0,numin
        if (in+ip+id+it+ih+ia.ne.npart) goto 210
        ident=100000*in+10000*ip+1000*id+100*it+10*ih+ia
        do 220 idc=0,idnum
          if (xschannel(idc).le.0.) goto 220
          if (idchannel(idc).eq.ident) then
            add=0.
            type=-1
            if (ident.eq.0) type=0
            if (ident.eq.100000) type=1
            if (ident.eq.10000) type=2
            if (ident.eq.1000) type=3
            if (ident.eq.100) type=4
            if (ident.eq.10) type=5
            if (ident.eq.1) type=6
            if (ident.eq.200000) type=7
            if (type.eq.0.and.Crescue(102).ne.1.) add=xsdif(102)
            if (type.eq.1.and.Crescue(4).ne.1.) add=xsdif(4)
            if (type.eq.2.and.Crescue(103).ne.1.) add=xsdif(103)
            if (type.eq.3.and.Crescue(104).ne.1.) add=xsdif(104)
            if (type.eq.4.and.Crescue(105).ne.1.) add=xsdif(105)
            if (type.eq.5.and.Crescue(106).ne.1.) add=xsdif(106)
            if (type.eq.6.and.Crescue(107).ne.1.) add=xsdif(107)
            if (type.eq.7.and.Crescue(16).ne.1.) add=xsdif(16)
            if (add.eq.0.) goto 220
            ratio=1.+add/xschannel(idc)
            xschannel(idc)=xschannel(idc)*ratio
            xsgamchannel(idc)=xsgamchannel(idc)*ratio
            do 230 i1=0,numlev
              do 230 i2=0,numlev
                xsgamdischan(idc,i1,i2)=xsgamdischan(idc,i1,i2)*ratio
  230       continue
            if (flagspec) then
              do 240 type2=0,6
                do 250 nen=0,numen
                  xschannelsp(idc,type2,nen)=xschannelsp(idc,type2,nen)*
     +              ratio
  250           continue
  240         continue
            endif
            if (type.ge.0.and.type.le.6) then
              do 260 nex=0,numlev
                xsdisc(type,nex)=xsdisc(type,nex)*ratio
                xscompdisc(type,nex)=xscompdisc(type,nex)*ratio
                xsdirdisc(type,nex)=xsdirdisc(type,nex)*ratio
  260         continue
              xsexclcont(type)=xsexclcont(type)*ratio
              xsngn(type)=xsngn(type)*ratio
              xsexclusive(type)=xsexclusive(type)*ratio
              xsdisctot(type)=xsdisctot(type)*ratio
              xsdirdisctot(type)=xsdirdisctot(type)*ratio
              xscompdisctot(type)=xscompdisctot(type)*ratio
              xscompcont(type)=xscompcont(type)*ratio
              xsdircont(type)=xsdircont(type)*ratio
              xsconttot(type)=xsconttot(type)*ratio
              xsdirect(type)=xsdirect(type)*ratio
              xsbinary(type)=xsbinary(type)*ratio
              xsparticle(type)=xsparticle(type)+add
              multiplicity(type)=xsparticle(type)/xsreacinc
            endif
            if (type.eq.7) then
              xsparticle(1)=xsparticle(1)+2.*add
              multiplicity(1)=xsparticle(1)/xsreacinc
            endif
          endif
  220   continue
  210 continue
      if (flagfission) then
        add=xsdif(18)
        if (xsfistot.gt.0.) then
          ratio=1.+add/xsfistot
          xsfistot=xsfistot*ratio
          do 270 idc=0,idnum
            xsfischannel(idc)=xsfischannel(idc)*ratio
  270     continue
        endif
      endif
c
c Put difference in the elastic (or total) cross section
c
c xselastot  : total elastic cross section (shape + compound)
c xstotadjust: total cross section adjustment
c xseladjust : elastic cross section adjustment
c xsnonel    : non-elastic cross section
c xsnonadjust: nonelastic cross section adjustment
c channelsum : sum over exclusive channel cross sections
c
      if (Crescue(2).ne.1..and.Crescue(2).ne.0.) then
        xseladjust(nin)=xsdif(2)
        xselastot=xselastot+xseladjust(nin)
        xstotadjust(nin)=xsdif(2)+xsadd
        xstotinc=xstotinc+xstotadjust(nin)
      else
        xstotadjust(nin)=xsdif(1)
        xstotinc=xstotinc+xstotadjust(nin)
        xseladjust(nin)=xsdif(1)-xsadd
        xselastot=xselastot+xseladjust(nin)
      endif
      xsnonadjust(nin)=xsadd
      xsnonel=xsnonel+xsnonadjust(nin)
      channelsum=channelsum+xsadd
      return
      end
Copyright (C) 2004  A.J. Koning, S. Hilaire and M.C. Duijvestijn
