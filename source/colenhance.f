      subroutine colenhance(Zix,Nix,Eex,ald,ibar,Krot,Kvib,Kcoll)
c
c +---------------------------------------------------------------------
c | Author: Marieke Duijvestijn and Arjan Koning
c | Date  : August 15, 2006
c | Task  : Collective enhancement
c +---------------------------------------------------------------------
c
c ******************* Declarations and common blocks *******************
c
      include "talys.cmb"
      integer Zix,Nix,ibar,A
      real    Eex,ald,Krot,Kvib,Kcoll,U,aldgs,ignatyuk,Temp,damper,expo,
     +        spincutgs,spincutbf,spincut
c
c **************** Collective enhancement for fission  *****************
c
c Zix      : charge number index for residual nucleus
c Nix      : neutron number index for residual nucleus
c Eex      : excitation energy
c ald,aldgs: level density parameter
c ibar     : fission barrier number, zero for states on ground state
c Krot     : rotational enhancement factor
c Kvib     : vibrational enhancement
c Kcoll    : total collective enhancement
c flagcol  : flag for collective enhancement of level density
c ldmodel  : level density model
c Exmatch  : matching energy for Ex
c U        : excitation energy minus pairing energy
c delta    : energy shift
c pair     : total pairing correction              
c Pshift   : adjustable pairing shift
c ignatyuk : function for energy dependent level density parameter a
c Temp     : nuclear temperature 
c
      Krot=1.
      Kvib=1.
      Kcoll=1.
      if (.not.flagcol) return
      if (ldmodel.eq.1.and.Eex.lt.Exmatch(Zix,Nix,ibar)) return
      U=Eex-delta(Zix,Nix,ibar)
      if (U.le.0.) return
      if (ibar.eq.0) then
        aldgs=ald
      else
        aldgs=ignatyuk(Zix,Nix,Eex,0) 
      endif
      Temp=sqrt(U/aldgs)
c
c Calculation of Kvib     
c
c AA,A    : mass number of residual nucleus
c twothird: 2/3
c
      A=AA(Zix,Nix,0)
      Kvib=exp(0.0555*(A**twothird)*(Temp**(4./3.)))
c
c Calculation of damping function and Krot
c
c damper      : energy damping function 
c expo        : exponent
c Ufermi      : energy of Fermi distribution for damping of ground-state
c             : rotational effects
c cfermi      : width of Fermi distribution for damping of ground-state
c             : rotational effects
c spincutgs   : spin-cutoff parameter squared (perpendicular projection)
c               for ground state
c Krotconstant: normalization constant for rotational enhancement
c Irigid      : rigid body value of moment of inertia
c
c Ground state
c
      if (ibar.eq.0) then
        damper=0.
        expo=(U-Ufermi)/cfermi
        if (expo.le.80.) damper=1./(1.+exp(expo))
        spincutgs=Krotconstant(Zix,Nix,0)*Irigid(Zix,Nix,0)*Temp
        Krot=max(spincutgs,1.)
c
c Fission barrier
c
c Ufermibf : energy of Fermi distribution for damping of barrier
c          : rotational effects
c cfermibf : width of Fermi distribution for damping of barrier
c          : rotational effects      
c spincutbf: spin-cutoff parameter squared (perpendicular projection)
c            for fission barrier
c axtype   : type of axiality of barrier (1: axial, 2: tri-axial)
c twopi    : 2.*pi
c spincut  : spin cutoff factor
c
      else
        damper=0.
        expo=(U-Ufermibf)/cfermibf
        if (expo.le.80.) damper=1./(1.+exp(expo))
        spincutbf=Krotconstant(Zix,Nix,ibar)*Irigid(Zix,Nix,ibar)*
     +    Temp
        if (axtype(Zix,Nix,ibar).eq.2) then
          ald=ignatyuk(Zix,Nix,Eex,ibar) 
          spincutbf=spincutbf*2*sqrt(twopi)*
     +      sqrt(spincut(Zix,Nix,ald,Eex,ibar))
        endif
        Krot=max(spincutbf,1.)
        if (ibar.eq.2) Krot=2.*Krot
      endif
      Kcoll=1.+(Krot*Kvib-1.)*damper
      Kcoll=max(Kcoll,1.)
      return
      end
Copyright (C) 2004  A.J. Koning, S. Hilaire and M.C. Duijvestijn
