      subroutine basicinitial
c
c +---------------------------------------------------------------------
c | Author: Arjan Koning
c | Date  : January 24, 2023
c | Task  : Initialization of arrays for basic cross sections
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      integer nen
c
c *************************** Initialization ***************************
c
c numen     : maximum number of outgoing energies
c lmax      : maximal l-value for transmission coefficients
c gammax    : number of l-values for gamma multipolarity
c xstot     : total cross section (neutrons only) for
c xsreac    : reaction cross section
c xsopt     : optical model reaction cross section
c xselas    : total elastic cross section (neutrons only)
c numl      : maximum l-value (set in talys.cmb)
c Tl        : transmission coefficients as a function of particle type,
c             energy and l-value (averaged over spin)
c Tjl       : transmission coefficients as a function of particle type,
c             energy, spin and l-value
c
      lmax=0
      do nen=0,numen
        lmax(0,nen)=gammax
      enddo
      xstot=0.
      xsreac=0.
      xsopt=0.
      xselas=0.
      Tl=0.
      Tjl=0.
      return
      end
Copyright (C)  2023 A.J. Koning, S. Hilaire and S. Goriely
