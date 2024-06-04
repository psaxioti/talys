      subroutine basicxs(Zcomp,Ncomp)
c
c +---------------------------------------------------------------------
c | Author: Arjan Koning
c | Date  : July 7, 2004
c | Task  : Basic cross sections and transmission coefficients
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      integer Zcomp,Ncomp
c
c ******************* ECIS calculations and output *********************
c
c Zcomp       : charge number index for compound nucleus
c Ncomp       : neutron number index for compound nucleus
c basicinitial: subroutine for initialization of arrays for basic cross 
c               sections
c inverse     : subroutine for ECIS calculation of total, reaction and
c               elastic cross sections and transmission coefficients 
c               for outgoing energy grid.
c parinclude  : logical to include outgoing particle
c gamma       : subroutine for gamma cross section and transmission 
c               coefficients
c
c The transmission coefficients and inverse reaction cross sections for
c the outgoing energy grid need to be calculated only once, for the
c maximal incident energy.
c
      call basicinitial
      call inverse(Zcomp,Ncomp)
      if (parinclude(0)) call gamma(Zcomp,Ncomp)
      return
      end
Copyright (C) 2004  A.J. Koning, S. Hilaire and M.C. Duijvestijn
