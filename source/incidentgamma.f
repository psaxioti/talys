      subroutine incidentgamma
c
c +---------------------------------------------------------------------
c | Author: Arjan Koning 
c | Date  : October 25, 2004
c | Task  : Incident photons
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      integer l
      real    gammaxs,fstrength,Tgamma
c
c **** Photo-absorption cross section and transmission coefficients ****
c
c lmaxinc  : maximal l-value for transmission coefficients for
c            incident channel
c gammax   : number of l-values for gamma multipolarity
c xsreacinc: reaction cross section for incident channel 
c gammaxs  : function for gamma ray cross sections 
c Einc     : incident energy in MeV
c Tgamma   : gamma transmission coefficient
c twopi    : 2.*pi
c fstrength: gamma ray strength function
c Tjlinc   : transmission coefficients as a function of spin and l for
c            the incident channel 
c
      lmaxinc=gammax
      xsreacinc=gammaxs(0,0,Einc)
      do 10 l=1,gammax
        Tgamma=twopi*(Einc**(2*l+1))*fstrength(0,0,Einc,1,l)
        Tjlinc(0,l)=Tgamma
   10 continue       
      return
      end
Copyright (C) 2004  A.J. Koning, S. Hilaire and M.C. Duijvestijn
