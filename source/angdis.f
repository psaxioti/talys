      subroutine angdis
c
c +---------------------------------------------------------------------
c | Author: Arjan Koning 
c | Date  : October 14, 2004
c | Task  : Calculation of angular distributions for discrete states
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      integer iang,LL,type,nex
      real    ang,x,leg(0:2*numJ,0:numang),plegendre
c
c ************* Create arrays with Legendre polynomials ****************
c
c nangle   : number of angles
c ang      : angle
c angle    : angle in degrees
c deg2rad  : conversion factor for degrees to radians     
c x        : help variable
c J2end    : 2 * end of J summation
c leg      : Legendre polynomial
c plegendre: function for calculation of Legendre polynomial
c
      do 10 iang=0,nangle
        ang=angle(iang)*deg2rad
        x=cos(ang)
        do 20 LL=0,J2end,2
          leg(LL,iang)=plegendre(LL,x)
   20   continue
   10 continue
c
c ******** Calculation of compound and total angular distribution ******
c
c parskip : logical to skip outgoing particle
c Nlast   : last discrete level
c parZ    : charge number of particle
c parN    : neutron number of particle
c cleg    : compound nucleus Legendre coefficient
c compad  : compound angular distribution
c discad  : discrete state angular distribution
c directad: direct angular distribution
c
      do 110 type=0,6
        if (parskip(type)) goto 110
        do 120 nex=0,Nlast(parZ(type),parN(type),0)
c
c For some unknown reason, compound inelastic inelastic scattering 
c to 0+ states can have a L=2 Legendre coefficient that is larger
c than the L=0 Legendre coefficient. This only happens for high 
c incident energies, typically above 10 MeV, i.e. for cases that 
c compound inelastic scattering to indvidual states is negligible. 
c Nevertheless, to avoid complaints by ENDF6 checking codes we put 
c the L=2 Legendre coefficient equal to the L=0 coefficient for these
c rare cases.
c
          if (type.ne.k0.or.nex.gt.0) 
     +      cleg(type,nex,2)=min(cleg(type,nex,2),cleg(type,nex,0))
          do 130 LL=0,J2end,2
            do 140 iang=0,nangle
              compad(type,nex,iang)=compad(type,nex,iang)+
     +          (2*LL+1)*cleg(type,nex,LL)*leg(LL,iang)
  140       continue
  130     continue
          do 150 iang=0,nangle
            discad(type,nex,iang)=directad(type,nex,iang)+
     +        compad(type,nex,iang)
  150     continue
  120   continue
  110 continue
c
c ************ Total Legendre coefficients and normalization ***********
c
c tleg     : total Legendre coefficient
c dleg     : direct reaction Legendre coefficient
c xspopex0 : binary population cross section
c tlegnor  : total Legendre coefficient normalized to 1
c k0       : index of incident particle
c xselastot: total elastic cross section (shape + compound)
c
      do 210 type=0,6
        if (parskip(type)) goto 210
        do 220 nex=0,Nlast(parZ(type),parN(type),0)
          do 230 LL=0,J2end
            tleg(type,nex,LL)=dleg(type,nex,LL)+cleg(type,nex,LL)
            if (xspopex0(type,nex).ge.1.e-20) tlegnor(type,nex,LL)=
     +        tleg(type,nex,LL)/xspopex0(type,nex)
  230     continue
  220   continue
  210 continue
      if (k0.eq.1) then
        do 240 LL=0,J2end
          tlegnor(k0,0,LL)=tleg(k0,0,LL)/xselastot
  240   continue
      endif
c
c ***************************** Recoils ********************************
c
c flagrecoil  : flag for calculation of recoils
c angdisrecoil: subroutine for recoil angular distributions for 
c               discrete states
c
      if (flagrecoil) call angdisrecoil
      return
      end
Copyright (C) 2004  A.J. Koning, S. Hilaire and M.C. Duijvestijn
