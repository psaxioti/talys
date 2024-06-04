      subroutine optical(Zix,Nix,kopt,eopt)
c
c +---------------------------------------------------------------------
c | Author: Arjan Koning
c | Date  : July 7, 2004
c | Task  : Determination of optical potential
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      implicit none
      integer  Zix,Nix,kopt
      real     eopt
c
c ********************** Call optical model module *********************
c
c Zix          : charge number index for residual nucleus
c Nix          : neutron number index for residual nucleus
c kopt         : index for fast particle
c eopt         : incident energy
c opticaln,....: optical model for neutrons, protons, deuterons, etc.
c
      if (kopt.eq.1) call opticaln(Zix,Nix,eopt)
      if (kopt.eq.2) call opticalp(Zix,Nix,eopt)
      if (kopt.eq.3) call opticald(Zix,Nix,eopt)
      if (kopt.eq.4) call opticalt(Zix,Nix,eopt)
      if (kopt.eq.5) call opticalh(Zix,Nix,eopt)
      if (kopt.eq.6) call opticala(Zix,Nix,eopt)
      return
      end
Copyright (C) 2004  A.J. Koning, S. Hilaire and M.C. Duijvestijn
