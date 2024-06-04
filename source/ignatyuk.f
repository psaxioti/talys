      function ignatyuk(Zix,Nix,Eex,ibar)
c
c +---------------------------------------------------------------------
c | Author: Arjan Koning
c | Date  : May 10, 2006
c | Task  : Energy dependent level density parameter a
c +---------------------------------------------------------------------
c
c ******************* Declarations and common blocks *******************
c
      include "talys.cmb"
      integer Zix,Nix,ibar
      real    ignatyuk,Eex,U,damp,fU
c
c *********************** Level density formula ************************
c
c Zix     : charge number index for residual nucleus
c Nix     : neutron number index for residual nucleus
c ignatyuk: function for energy dependent level density parameter a
c Eex     : excitation energy
c ibar    : fission barrier
c U       : excitation energy minus pairing energy
c delta   : energy shift
c fU,damp : help variables
c gammald : gamma-constant for asymptotic level density parameter
c deltaW  : shell correction in nuclear mass 
c alimit  : asymptotic level density parameter
c
c Formalism from Ignatyuk et al. Sov. Jour. Nuc. Phys. 21 (1975), 255. 
c
      U=Eex-delta(Zix,Nix,ibar)
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
      ignatyuk=max(alimit(Zix,Nix)*damp,1.)
      return
      end
Copyright (C) 2004  A.J. Koning, S. Hilaire and M.C. Duijvestijn
