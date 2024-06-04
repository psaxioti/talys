      function funcfismode(Zix,Nix,epar)
c
c +---------------------------------------------------------------------
c | Author: Marieke Duijvestijn
c | Date  : September 9, 2004
c | Task  : Transmission coefficient per fission mode
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      integer          Zix,Nix
      real             epar,temps(9),ald,tmp,bft,hwt,ignatyuk
      double precision density,funcfismode,expo
      data temps /0.0,0.3,0.6,0.9,1.2,1.6,2.0,2.5,3.0/
c
c **********************************************************************
c
      ald=ignatyuk(Zix,Nix,epar,0)
      Tmp=sqrt(Epar/ald)
      call splint(temps,bf,bfsplin,numtemp,tmp,bft)
      call splint(temps,hw,hwsplin,numtemp,tmp,hwt)
      expo=2.*pi*(bft+epar-excfis)/hwt
      if (expo.le.80.) then
        funcfismode=1./(1.+exp(expo))
      else
        funcfismode=exp(-80.)/(exp(-80.)+exp(expo-80.))
      endif
      if (epar.gt.0.) 
     +  funcfismode=funcfismode*density(Zix,Nix,epar,0.,0,1)
      return
      end
Copyright (C) 2004  A.J. Koning, S. Hilaire and M.C. Duijvestijn
