      subroutine gammanorm(Zcomp,Ncomp)
c
c +---------------------------------------------------------------------
c | Author: Stephane Hilaire and Arjan Koning 
c | Date  : September 11, 2007
c | Task  : Normalization of gamma ray strength functions 
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      integer Zcomp,Ncomp,nen,l,irad
      real    Egamma,Tgamma,fstrength,gammaxs
c
c We normalize the gamma-ray strength functions by imposing the 
c condition that the transmission coefficients integrated from zero up 
c to neutron separation energy are equal to the ratio of the 
c experimental mean gamma width and mean level spacing for s-wave 
c neutrons.
c
c ************** Normalization of transmission coefficients ***********
c
c Zcomp    : charge number index for compound nucleus
c Ncomp    : neutron number index for compound nucleus
c gnorm    : gamma normalization factor
c gamgam   : total radiative width
c gamgamth : theoretical total radiative width
c ebegin   : first energy point of energy grid
c eend     : last energy point of energy grid
c Egamma   : gamma energy
c egrid    : outgoing energy grid
c lmax     : maximal l-value for transmission coefficients
c gammax   : number of l-values for gamma multipolarity
c Tgamma   : gamma transmission coefficient
c twopi    : 2.*pi
c fstrength: gamma ray strength function
c Tjl      : transmission coefficients as a function of particle type, 
c            energy, spin and l-value
c
      if (gnorm.eq.-1.) then
        if (gamgamth(Zcomp,Ncomp).ne.0.) then
          gnorm=gamgam(Zcomp,Ncomp)/gamgamth(Zcomp,Ncomp)
        else
          gnorm=1.
        endif
      endif
      do 10 nen=ebegin(0),eend(0)
        Egamma=egrid(nen)
        lmax(0,nen)=gammax
        do 20 l=1,gammax
          do 20 irad=0,1
            Tgamma=twopi*(Egamma**(2*l+1))*
     +        fstrength(Zcomp,Ncomp,Egamma,irad,l)
            Tjl(0,nen,irad,l)=Tgamma*gnorm
   20   continue
c
c Photo-absorption cross section
c
c xsreac : reaction cross section   
c gammaxs: function for gamma ray cross sections 
c
        xsreac(0,nen)=gammaxs(Zcomp,Ncomp,Egamma)
   10 continue
      return
      end
Copyright (C) 2004  A.J. Koning, S. Hilaire and M.C. Duijvestijn
