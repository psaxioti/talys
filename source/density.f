      function density(Zix,Nix,Eex,Rspin,parity,ibar,ldmod)
c
c +---------------------------------------------------------------------
c | Author: Arjan Koning
c | Date  : December 18, 2007
c | Task  : Level density 
c +---------------------------------------------------------------------
c
c ******************* Declarations and common blocks *******************
c
c Note that the parity is not used explicitly for analytical level 
c densities. For those, an equidistant parity distribution is assumed.
c
      include "talys.cmb"
      integer          Zix,Nix,parity,ibar,ldmod,jj,nex2
      real             Eex,Rspin,ald,ignatyuk,spindis,Eshift
      double precision density,densitytot,eb,ee,ldb,lde,ldtab,cctable
c
c ************************** Level density *****************************
c
c density   : level density
c Zix       : charge number index for residual nucleus
c Nix       : neutron number index for residual nucleus
c Eex       : excitation energy
c Rspin     : spin
c parity    : parity
c ibar      : fission barrier number, zero for states on ground state
c ldmod     : level density model
c ldexist   : flag for existence of level density table
c ald       : level density parameter 
c ignatyuk  : function for energy dependent level density parameter a 
c densitytot: total level density
c pardis    : parity distribution 
c spindis   : Wigner spin distribution
c
c 1. Gilbert and Cameron
c 2. Back-shifted Fermi gas
c 3. Superfluid model
c
      density=0.
      if (Eex.lt.0.) goto 100
      if (ldmod.le.3.or..not.ldexist(Zix,Nix,ibar)) then
        ald=ignatyuk(Zix,Nix,Eex,ibar)
        density=densitytot(Zix,Nix,Eex,ibar,ldmod)*pardis*
     +    spindis(Zix,Nix,Eex,ald,Rspin,ibar)
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
c ldtable            : level density from table
c expo               : exponent
c cctable            : constant to adjust tabulated level densities
c
        jj=min(29,int(Rspin))
        Eshift=Eex-ptable(Zix,Nix,ibar)
        if (Eshift.le.0.) goto 100
        if (Eshift.le.Edensmax) then
          call locate(edens,0,nendens,Eshift,nex2) 
          eb=edens(nex2)
          ee=edens(nex2+1)
          ldb=ldtable(Zix,Nix,nex2,jj,parity,ibar)
          lde=ldtable(Zix,Nix,nex2+1,jj,parity,ibar)
          if (ldb.le.1..or.lde.le.1.) then
            ldtab=ldb+(Eshift-eb)/(ee-eb)*(lde-ldb)
          else
            ldtab=exp(log(ldb)+(Eshift-eb)/(ee-eb)*
     +       (log(lde)-log(ldb)))
          endif
        else
          eb=edens(nendens-1)
          ee=edens(nendens)
          ldb=log(ldtable(Zix,Nix,nendens-1,jj,parity,ibar))
          lde=log(ldtable(Zix,Nix,nendens,jj,parity,ibar))
          ldtab=exp(ldb+(Eshift-eb)/(ee-eb)*(lde-ldb))
        endif
        cctable=exp(dble(ctable(Zix,Nix,ibar)*sqrt(Eshift)))
        density=cctable*ldtab
      endif
  100 density=max(density,1.d-30)
      return
      end
Copyright (C) 2004  A.J. Koning, S. Hilaire and M.C. Duijvestijn
