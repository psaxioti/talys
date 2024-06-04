      function density(Zix,Nix,Eex,Rspin,ibar,ldmod)
c
c +---------------------------------------------------------------------
c | Author: Arjan Koning
c | Date  : October 15, 2004
c | Task  : Total level density 
c +---------------------------------------------------------------------
c
c ******************* Declarations and common blocks *******************
c
c Although the parity is not an argument of this function, it is 
c included in the level density value, i.e. the final level density is
c per parity. An equidistant parity distribution is assumed.
c
      include "talys.cmb"
      integer          Zix,Nix,ibar,ldmod,jj,nex2
      real             Eex,Rspin,ald,ignatyuk,Krot,Kvib,spindis
      double precision density,gilcam,eb,ee,ldb,lde,ldtab,expo,ctable
c
c *********************** Total level density **************************
c
c density : level density
c Zix     : charge number index for residual nucleus
c Nix     : neutron number index for residual nucleus
c Eex     : excitation energy
c Rspin   : spin
c ibar    : fission barrier number, zero for states on ground state
c ldmod   : level density model
c ldexist : flag for existence of level density table
c ald     : level density parameter 
c ignatyuk: function for energy dependent level density parameter a 
c
c 1. Gilbert and Cameron: effective
c 2. Gilbert and Cameron: explicit collective enhancement
c
      density=0.
      if (Eex.lt.0.) goto 100
      if (ldmod.le.2.or..not.ldexist(Zix,Nix,ibar)) then
        ald=ignatyuk(Zix,Nix,Eex,ibar)
c
c Enhancement only in the Fermi gas region
c Krot (and Kvib) will be determined.
c
c colenhance: subroutine for collective enhancement
c Kvib      : vibrational enhancement factor
c Krot      : rotational enhancement factor
c gilcam    : Gilbert-Cameron level density formula
c spindis   : Wigner spin distribution
c pardis    : parity distribution
c
        call colenhance(Zix,Nix,Eex,ald,ibar,Krot,Kvib)
        density=Krot*Kvib*gilcam(Zix,Nix,Eex,ald,ibar)*
     +    spindis(Zix,Nix,Eex,ald,Rspin,ibar)*pardis
      endif
c
c 3. Tabulated level densities     
c
c locate             : subroutine to find value in ordered table 
c edens              : energy grid for tabulated level densities
c eb,ee,ldb,lde,ldtab: help variables
c ldtable            : level density from table
c expo               : exponent
c c1table,c2table    : constant to adjust tabulated level densities
c
      if (ldmod.eq.3.and.ldexist(Zix,Nix,ibar)) then
        jj=min(29,int(Rspin))
        if (Eex.le.150.) then
          call locate(edens,0,55,Eex,nex2) 
          eb=edens(nex2)
          ee=edens(nex2+1)
          ldb=ldtable(Zix,Nix,nex2,jj,ibar)
          lde=ldtable(Zix,Nix,nex2+1,jj,ibar)
          ldtab=ldb+(Eex-eb)/(ee-eb)*(lde-ldb)
        else
          eb=edens(54)
          ee=edens(55)
          ldb=log(ldtable(Zix,Nix,54,jj,ibar))
          lde=log(ldtable(Zix,Nix,55,jj,ibar))
          ldtab=exp(ldb+(Eex-eb)/(ee-eb)*(lde-ldb))
        endif
        expo=c2table(Zix,Nix)*sqrt(Eex)
        ctable=c1table(Zix,Nix)*exp(expo)
        density=ctable*ldtab*pardis
      endif
  100 density=max(density,1.d-30)
      return
      end
Copyright (C) 2004  A.J. Koning, S. Hilaire and M.C. Duijvestijn
