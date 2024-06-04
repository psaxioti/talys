      subroutine tgamma(Zcomp,Ncomp)
c
c +---------------------------------------------------------------------
c | Author: Stephane Hilaire and Arjan Koning
c | Date  : April 4, 2022
c | Task  : Photon transmission coefficients
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      integer      Zcomp,Ncomp,nen,l,irad
      real         Egamma,fstrength,gammaxs
c
c ************** Normalization of transmission coefficients ***********
c
c Zcomp       : charge number index for compound nucleus
c Ncomp       : neutron number index for compound nucleus
c ebegin      : first energy point of energy grid
c eend        : last energy point of energy grid
c Einc        : incident energy in MeV
c Egamma      : gamma energy
c egrid       : outgoing energy grid
c lmax        : maximal l-value for transmission coefficients
c gammax      : number of l-values for gamma multipolarity
c twopi       : 2.*pi
c fstrength   : gamma ray strength function
c Tjl         : transmission coefficients as a function of
c               particle type, energy, spin and l-value
c
      do 10 nen=ebegin(0),eend(0)
        Egamma=egrid(nen)
        lmax(0,nen)=gammax
        do 20 l=1,gammax
          do 20 irad=0,1
            Tjl(0,nen,irad,l)=twopi*(Egamma**(2*l+1))*
     +        fstrength(Zcomp,Ncomp,Einc,Egamma,irad,l)
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
Copyright (C)  2022 A.J. Koning, S. Hilaire and S. Goriely
