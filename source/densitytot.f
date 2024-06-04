      function densitytot(Zix,Nix,Eex,ibar,ldmod)
c
c +---------------------------------------------------------------------
c | Author: Arjan Koning
c | Date  : December 18, 2007
c | Task  : Total level density 
c +---------------------------------------------------------------------
c
c ******************* Declarations and common blocks *******************
c
      include "talys.cmb"
      integer          Zix,Nix,ibar,ldmod,nex2
      real             Eex,P,ald,ignatyuk,Krot,Kvib,Kcoll,Eshift
      double precision densitytot,dens,gilcam,bsfgmodel,superfluid,
     +                 eb,ee,ldb,lde,ldtab,cctable
c
c *********************** Total level density **************************
c
c densitytot : total level density
c Zix        : charge number index for residual nucleus
c Nix        : neutron number index for residual nucleus
c Eex        : excitation energy
c ibar       : fission barrier number, zero for states on ground state
c ldmod      : level density model
c ldexist    : flag for existence of level density table
c
      densitytot=0.
      if (Eex.lt.0.) goto 100
      if (ldmod.le.3.or..not.ldexist(Zix,Nix,ibar)) then
c
c ald       : level density parameter 
c ignatyuk  : function for energy dependent level density parameter a 
c colenhance: subroutine for collective enhancement
c Krot      : rotational enhancement factor
c Kvib      : vibrational enhancement factor
c Kcoll     : total collective enhancement
c delta,P   : energy shift
c dens      : total level density
c gilcam    : Gilbert-Cameron level density formula
c bsfgmodel : Back-shifted Fermi gas level density formula
c superfluid: Superfluid model level density formula
c
c Kcoll will be determined.
c
        ald=ignatyuk(Zix,Nix,Eex,ibar)
        call colenhance(Zix,Nix,Eex,ald,ibar,Krot,Kvib,Kcoll)
        P=delta(Zix,Nix,ibar)
c
c 1. Gilbert and Cameron
c
        if (ldmod.eq.1.or.(ldmodel.ge.4.and..not.ldexist(Zix,Nix,ibar)))
     +    dens=gilcam(Zix,Nix,ald,Eex,P,ibar)
c
c 2. Back-shifted Fermi gas
c
        if (ldmod.eq.2) dens=bsfgmodel(Zix,Nix,ald,Eex,P,ibar)
c
c 3. Superfluid model
c
        if (ldmod.eq.3) dens=superfluid(Zix,Nix,ald,Eex,P,ibar)
        densitytot=Kcoll*dens
      else
c
c 4. Tabulated level densities     
c
c Eshift             : shifted excitation energy
c ctable,ptable      : constant to adjust tabulated level densities
c Edensmax           : maximum energy on level density table
c locate             : subroutine to find value in ordered table 
c edens              : energy grid for tabulated level densities
c nendens            : number of energies for level density grid
c eb,ee,ldb,lde,ldtab: help variables
c ldtottable         : total level density from table
c
        Eshift=Eex-ptable(Zix,Nix,ibar)
        if (Eshift.le.0.) goto 100
        if (Eshift.le.Edensmax) then
          call locate(edens,0,nendens,Eshift,nex2) 
          eb=edens(nex2)
          ee=edens(nex2+1)
          ldb=ldtottable(Zix,Nix,nex2,ibar)
          lde=ldtottable(Zix,Nix,nex2+1,ibar)
          ldtab=ldb+(Eshift-eb)/(ee-eb)*(lde-ldb)
        else
          eb=edens(nendens-1)
          ee=edens(nendens)
          ldb=log(ldtottable(Zix,Nix,nendens-1,ibar))
          lde=log(ldtottable(Zix,Nix,nendens,ibar))
          ldtab=exp(ldb+(Eshift-eb)/(ee-eb)*(lde-ldb))
        endif
        cctable=exp(ctable(Zix,Nix,ibar)*sqrt(Eshift))
        densitytot=cctable*ldtab
      endif
  100 densitytot=max(densitytot,1.d-30)
      return
      end
Copyright (C) 2004  A.J. Koning, S. Hilaire and M.C. Duijvestijn
