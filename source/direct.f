      subroutine direct
c
c +---------------------------------------------------------------------
c | Author: Arjan Koning 
c | Date  : July 12, 2004
c | Task  : Calculation of direct inelastic cross sections 
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
c
c ****************** Direct cross section calculation ******************
c
c k0        : index for incident particle  
c Ltarget   : excited level of target
c directecis: subroutine for ECIS calculation of direct cross section
c directread: subroutine to read ECIS results for direct cross section 
c flaggiant : flag for collective contribution from giant resonances
c giant     : subroutine for giant resonance contribution     
c flagdirect: flag for output of direct reaction results
c directout : subroutine for output of direct reaction cross sections
c
      if (k0.eq.0) return
      if (Ltarget.ne.0) return
      call directecis
      call directread
      if (flaggiant) call giant
      if (flagdirect) call directout
      return
      end
Copyright (C) 2004  A.J. Koning, S. Hilaire and M.C. Duijvestijn
