      subroutine astroinit
c
c +---------------------------------------------------------------------
c | Author: Stephane Hilaire and Stephane Goriely
c | Date  : January 24, 2023
c | Task  : Initialization of astrophysics quantities
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
c
c **** Initialization of arrays of astrophysics interest ***************
c
c maxNastro   : maximal number of neutrons away from the initial
c               compound nucleus for astrophysical calculations
c maxZastro   : maximal number of protons away from the initial
c               compound nucleus for astrophysical calculations
c numT        : number of temperatures
c numinc      : number of incident energies
c xsastro     : cross section for astrophysical calculation
c xsastroex   : cross section for astrophysical calculation to a given excited state
c rateastro   : thermonuclear reaction rate factor
c rateastroex : thermonuclear reaction rate factor to a given excited state
c macsastro   : Maxwellian-averaged thermonuclear reaction cross section
c partf       : integrated partition function
c rateastrofis: thermonuclear reaction rate factor for fission
c macsastro   : thermonuclear reaction cross section
c macsastroex : thermonuclear reaction cross section to a given excited state
c macsastrofis: thermonuclear reaction cross section for fission
c xsastrofis  : astrophysical fission cross section
c
      maxZastro=numZastro
      maxNastro=numNastro
      xsastro=0.
      xsastroex=0.
      rateastro=0.
      macsastro=0.
      rateastroex=0.
      macsastroex=0.
      partf=0.
      rateastrofis=0.
      rateastroracap=0.
      macsastrofis=0.
      macsastroracap=0.
      xsastrofis=0.
      return
      end
Copyright (C)  2013 A.J. Koning, S. Hilaire and S. Goriely
