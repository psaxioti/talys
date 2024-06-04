      function preeqpair(Zix,Nix,n,E)
c
c +---------------------------------------------------------------------
c | Author: Marieke Duijvestijn
c | Date  : August 9, 2004
c | Task  : Pairing effects in exciton model
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      integer Zix,Nix,n
      real    preeqpair,E,gsp,gsn,gs,pair0,Tcrit,ncrit,pairex
c
c **************************** Fu formula ******************************
c
c preeqpair: pre-equilibrium pairing energy 
c Zix      : charge number index for residual nucleus
c Nix      : neutron number index for residual nucleus
c n        : exciton number
c E        : excitation energy
c pair     : pairing energy
c flag2comp: flag for two-component pre-equilibrium model 
c gp,gsp   : single-particle proton level density parameter  
c gn,gsn   : single-particle neutron level density parameter
c gs       : single-particle level density parameter
c pair0    : ground-state pairing gap
c Tcrit    : critical temperature
c ncrit    : average quasi-particle number at the Tcrit
c pairex   : excited state pairing gap
c
      preeqpair=0.
      if (pair(Zix,Nix).eq.0.) return
      if (flag2comp) then
        gsp=gp(Zix,Nix)
        gsn=gn(Zix,Nix)
        gs=gsp+gsn
      else
        gs=g(Zix,Nix)         
      endif
      pair0=sqrt(pair(Zix,Nix)/(0.25*gs))
      Tcrit=2.*pair0/3.5
      ncrit=2.*gs*Tcrit*log(2.)
      if (E/pair(Zix,Nix).ge.(0.716+2.44*(n/ncrit)**2.17)) then
        pairex=pair0*(0.996-1.76*(n/ncrit)**1.6/
     +    (E/pair(Zix,Nix)**0.68))
      else
        pairex=0.
      endif
      preeqpair=pair(Zix,Nix)-0.25*gs*(pairex**2.)
      return
      end
Copyright (C) 2004  A.J. Koning, S. Hilaire and M.C. Duijvestijn
