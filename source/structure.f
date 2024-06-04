      subroutine structure(Zix,Nix)
c
c +---------------------------------------------------------------------
c | Author: Arjan Koning
c | Date  : October 14, 2004
c | Task  : Nuclear structure parameters
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      integer Zix,Nix
c
c ********** Calculate and read various nuclear parameters *************
c
c Zix         : charge number index for residual nucleus
c Nix         : neutron number index for residual nucleus   
c levels      : subroutine for discrete levels
c flagendf    : flag for information for ENDF-6 file
c primary     : flag to designate primary (binary) reaction
c gammadecay  : subroutine for scheme for discrete gamma decay  
c deformpar   : subroutine for deformation parameters
c parinclude  : logical to include outgoing particle
c resonancepar: subroutine for s-wave resonance parameters
c gammapar    : subroutine for gamma ray parameters
c flagompall  : flag for new optical model calculation for all residual
c               nuclides
c omppar      : subroutine for optical model parameters
c flagfission : flag for fission
c fissionpar  : subroutine for fission parameters
c densitypar  : subroutine for level density parameters
c densitymatch: subroutine for level density matching solution
c ldmodel     : level density model
c densitytable: subroutine for tabulated level densities
c D0theory    : subroutine for theoretical calculation of D0
c flagpartable: flag for output of model parameters on separate file
c partable    : subroutine to write model parameters per nucleus to 
c               separate file
c
c All the nuclear structure info is read and/or calculated for the 
c nucleus under consideration.
c
      call levels(Zix,Nix)
      if (flagendf.and.primary) call gammadecay(Zix,Nix)
      call deformpar(Zix,Nix)
      if (parinclude(0)) then
        call resonancepar(Zix,Nix)
        call gammapar(Zix,Nix)
      endif
      if ((Zix.le.2.and.Nix.le.2).or.flagompall) call omppar(Zix,Nix)
      if (flagfission) call fissionpar(Zix,Nix)
      call densitypar(Zix,Nix)
      call densitymatch(Zix,Nix)
      if (ldmodel.eq.3) call densitytable(Zix,Nix)
      call D0theory(Zix,Nix)
      if (flagpartable) call partable(Zix,Nix)
      return
      end
Copyright (C) 2004  A.J. Koning, S. Hilaire and M.C. Duijvestijn
