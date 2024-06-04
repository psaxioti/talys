      subroutine cascade(Zcomp,Ncomp,nex)
c
c +---------------------------------------------------------------------
c | Author: Arjan Koning 
c | Date  : May 22, 2007
c | Task  : Gamma-ray cascade
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      integer          Zcomp,Ncomp,nex,J,parity,k,Jres,Pres
      double precision xsJP,intens,xsgamma
c
c ******************* Gamma-ray cascade ********************************
c
c Zcomp      : charge number index for compound nucleus
c Ncomp      : neutron number index for compound nucleus
c nex        : excitation energy bin of compound nucleus
c jdis,J,Jres: spin of level
c parity,Pres: parity
c parlev     : parity of level 
c xsJP,xspop : population cross section for spin and parity
c intens     : total gamma intensity
c branchratio: gamma-ray branching ratio to level
c xspopex    : population cross section summed over spin and parity
c xspartial  : emitted cross section flux per energy bin         
c numZchan   : maximal number of outgoing proton units in individual
c              channel description
c numNchan   : maximal number of outgoing neutron units in individual
c              channel description
c mcontrib   : contribution to emission spectrum
c
      J=int(jdis(Zcomp,Ncomp,nex))
      parity=parlev(Zcomp,Ncomp,nex)
      xsJP=xspop(Zcomp,Ncomp,nex,J,parity)
      do 10 k=0,nex
        Jres=int(jdis(Zcomp,Ncomp,k))
        Pres=parlev(Zcomp,Ncomp,k)
        intens=xsJP*branchratio(Zcomp,Ncomp,nex,k)
        xspop(Zcomp,Ncomp,k,Jres,Pres)=
     +    xspop(Zcomp,Ncomp,k,Jres,Pres)+intens
        xspopex(Zcomp,Ncomp,k)=xspopex(Zcomp,Ncomp,k)+intens
        xspopex(Zcomp,Ncomp,nex)=xspopex(Zcomp,Ncomp,nex)-intens
        xspartial(0,nex)=xspartial(0,nex)+intens
        if (Zcomp.le.numZchan.and.Ncomp.le.numNchan)
     +    mcontrib(0,nex,k)=intens
c
c ************ Storage of discrete gamma line intensities **************
c
c flaggamdis  : flag for output of discrete gamma-ray intensities
c flagelectron: flag for application of electron conversion coefficient
c xsgamma     : help variable
c conv        : conversion coefficient
c xsgamdis    : discrete gamma-ray cross section
c
        if (flaggamdis) then
          if (flagelectron) then
            xsgamma=intens/(1.+conv(Zcomp,Ncomp,nex,k))
          else
            xsgamma=intens
          endif
          xsgamdis(Zcomp,Ncomp,nex,k)=xsgamma
        endif
   10 continue
      return
      end
Copyright (C) 2004  A.J. Koning, S. Hilaire and M.C. Duijvestijn
