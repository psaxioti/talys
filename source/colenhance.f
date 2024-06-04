      subroutine colenhance(Zix,Nix,Eex,ald,ibar,Krot,Kvib)
c
c +---------------------------------------------------------------------
c | Author: Marieke Duijvestijn
c | Date  : October 23, 2004
c | Task  : Collective enhancement
c +---------------------------------------------------------------------
c
c ******************* Declarations and common blocks *******************
c
      include "talys.cmb"
      integer Zix,Nix,ibar,A
      real    Eex,ald,Krot,Kvib,U,ald0,ignatyuk,Tempgs,Tempbf,damprotgs,
     +        expo,damprotbf,dampdefgs,b2,b2gs,spincutgs,spincutbf,
     +        Krotgs,Krotbf
c
c **************** Collective enhancement for fission  *****************
c
c Zix        : charge number index for residual nucleus
c Nix        : neutron number index for residual nucleus
c Eex        : excitation energy
c ald,ald0   : level density parameter
c ibar       : fission barrier number, zero for states on ground state
c Krot       : rotational enhancement factor
c Kvib       : vibrational enhancement
c U          : excitation energy minus pairing energy
c pair       : total pairing correction              
c ldmodel    : level density model
c flagfission: flag for fission
c AA,A       : mass number of residual nucleus
c ignatyuk   : function for energy dependent level density parameter a
c Tempgs     : nuclear temperature for ground state
c Tempbf     : nuclear temperature for fission barrier
c
      Krot=1.
      Kvib=1.
      if (ldmodel.eq.1.and..not.flagfission) return
      U=Eex-pair(Zix,Nix)
      if (U.le.0.) return
      A=AA(Zix,Nix,0)
      if (ibar.eq.0) then
        ald0=ald
      else
        ald0=ignatyuk(Zix,Nix,Eex,0) 
      endif
      Tempgs=sqrt(U/ald0)
      Tempbf=Tempgs
c
c Calculation of Krot
c
c damprotgs   : rotational damping function for ground state
c expo        : exponent
c Ufermi      : energy of Fermi distribution for damping of ground-state
c             : rotational effects
c cfermi      : width of Fermi distribution for damping of ground-state
c             : rotational effects
c damprotbf   : rotational damping function for fission barrier
c Ufermibf    : energy of Fermi distribution for damping of barrier
c             : rotational effects
c cfermibf    : width of Fermi distribution for damping of barrier
c             : rotational effects      
c dampdefgs   : deformation damping function for ground state
c b2gs,b2     : help variables
c beta2       : deformation parameter
c spincutgs   : spin-cutoff parameter squared (perpendicular projection)
c               for ground state
c Krotconstant: normalization constant for rotational enhancement
c Krotgs      : rotational enhancement factor for ground state
c spincutbf   : spin-cutoff parameter squared (perpendicular projection)
c               for fission barrier
c Irigid0     : undeformed rigid body value of moment of inertia
c Irigid      : rigid body value of moment of inertia
c axtype      : type of axiality of barrier (1: axial, 2: tri-axial)
c pi2         : pi**2
c twothird    : 2/3
c Krotbf      : rotational enhancement factor for fission barrier
c
c
      damprotgs=0.
      expo=(U-Ufermi)/cfermi
      if (expo.gt.-80.) damprotgs=1./(1.+exp(-expo))
      damprotbf=0.
      expo=(U-Ufermibf)/cfermibf
      if (expo.gt.-80.) damprotbf=1./(1.+exp(-expo))
      dampdefgs=1./(1.+exp((Tempgs-2.0)/0.5))
      b2=beta2(Zix,Nix,ibar)
      b2gs=beta2(Zix,Nix,0)*dampdefgs
      if (b2gs.eq.0.) b2gs=0.012*dampdefgs
      spincutgs=Krotconstant(Zix,Nix,0)*Irigid0(Zix,Nix)*
     +  (1.+b2gs/3.)*Tempgs
      spincutgs=max(spincutgs,1.)
      Krotgs=spincutgs*(1.-damprotgs)+damprotgs
      if (ibar.gt.0) then
        spincutbf=Krotconstant(Zix,Nix,ibar)*Irigid(Zix,Nix,ibar)*
     +    Tempbf
        if (axtype(Zix,Nix,ibar).eq.2) then
          ald=ignatyuk(Zix,Nix,Eex,ibar) 
          spincutbf=spincutbf*2*sqrt(twopi)*
     +      sqrt((6./pi2)*ald*0.24*(A**twothird)*(1.-2.*b2/3.)*Tempbf)
        endif
        spincutbf=max(spincutbf,1.)
        Krotbf=spincutbf*(1.-damprotbf)+damprotbf
        if (ibar.eq.2) Krotbf=2.*Krotbf
      endif
c
c ldmodel 1: collective effects for the ground state are included 
c   implicitly in the intrinsic level density, and collective effects 
c   on the barrier are determined relative to the ground state.
c
c ldmodel 2: both ground state and barrier collective effects are 
c   included explicitly
c
      if (ibar.eq.0) then
        Krot=Krotgs
        if (ldmodel.eq.1) Krot=Krot/spincutgs
      else
        Krot=Krotbf
        if (ldmodel.eq.1) Krot=Krot/Krotgs
        Krot=max(Krot,1.)
      endif
c
c Calculation of Kvib     
c
      if (ldmodel.eq.2) Kvib=exp(0.0555*(A**twothird)*(Tempgs**(4./3.)))
      return
      end
Copyright (C) 2004  A.J. Koning, S. Hilaire and M.C. Duijvestijn
