      subroutine gamma(Zcomp,Ncomp)
c
c +---------------------------------------------------------------------
c | Author: Arjan Koning
c | Date  : April 4, 2022
c | Task  : Gamma cross section and transmission coefficients
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      integer Zcomp,Ncomp
c
c ********** Gamma cross section and transmission coefficients *********
c
c Zcomp    : charge number index for compound nucleus
c Ncomp    : neutron number index for compound nucleus
c tgamma   : subroutine for photon transmission coefficients
c gammaout : subroutine for output of gamma-ray strength functions,
c            transmission coefficients and cross sections
c
      call tgamma(Zcomp,Ncomp)
      if (flaggamma) call gammaout(Zcomp,Ncomp)
      return
      end
Copyright (C)  2013 A.J. Koning, S. Hilaire and S. Goriely
