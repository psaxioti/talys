      function ignatyuk(Zix,Nix,Eex,ibar)
c
c +---------------------------------------------------------------------
c | Author: Arjan Koning
c | Date  : July 5, 2004
c | Task  : Energy dependent level density parameter a
c +---------------------------------------------------------------------
c
c ******************* Declarations and common blocks *******************
c
      include "talys.cmb"
      integer Zix,Nix,ibar,A
      real    ignatyuk,Eex,U,damp,fU,aldlim,aldlow,expo,qfermi
c
c *********************** Level density formula ************************
c
c Zix     : charge number index for residual nucleus
c Nix     : neutron number index for residual nucleus
c ignatyuk: function for energy dependent level density parameter a
c Eex     : excitation energy
c ibar    : fission barrier
c U       : excitation energy minus pairing energy
c pair    : total pairing correction
c fU,damp : help variables
c gammald : gamma-constant for asymptotic level density parameter
c deltaW  : shell correction in nuclear mass 
c aldlim  : asymptotic level density parameter
c alimit  : asymptotic level density parameter
c
c Formalism from Ignatyuk et al. Sov. Jour. Nuc. Phys. 21 (1975), 255. 
c
      U=Eex-pair(Zix,Nix)
c
c 1. For very low Eex, i.e. U < 0, we use the first order Taylor 
c    expansion
c
      if (U.le.0.) then
        damp=(1.+deltaW(Zix,Nix,ibar)*gammald(Zix,Nix))
      else
c
c 2. Higher energies 
c
        fU=1.-exp(-gammald(Zix,Nix)*U)
        damp=1.+fU*deltaW(Zix,Nix,ibar)/U
      endif
      aldlim=alimit(Zix,Nix)
      if (damp.le.0.) damp=1.
c
c Fermi distribution for asymptotic level density parameter. This takes
c the damping of collective effects into account in a phenomenological 
c way. The level density parameters are then equal to A/13 in the 
c high-energy limit rather than A/8.
c
c flagcolldamp: flag for damping of collective effects
c AA,A        : mass number of residual nucleus
c aldlow      : lower limit of a
c expo        : help variable
c Ufermi      : parameter for Fermi distribution
c cfermi      : width of fermi distribution
c qfermi      : Fermi distribution
c
      if (flagcolldamp) then
        A=AA(Zix,Nix,0)
        aldlow=A/13.
        expo=(U-Ufermi)/cfermi
        qfermi=0.
        if (expo.gt.-80.) qfermi=1./(1.+exp(-expo))
        aldlim=aldlow*qfermi+aldlim*(1.-qfermi)
      endif
c
c Final level density parameter (the value is at least 1.)
c
      ignatyuk=max(aldlim*damp,1.)
      return
      end
Copyright (C) 2004  A.J. Koning, S. Hilaire and M.C. Duijvestijn
