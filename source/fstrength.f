      function fstrength(Zcomp,Ncomp,Egamma,irad,l)
c
c +---------------------------------------------------------------------
c | Author: Arjan Koning
c | Date  : July 7, 2004
c | Task  : Gamma ray strength functions 
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      integer Zcomp,Ncomp,irad,l,i
      real    fstrength,Egamma,sgr1,egr1,ggr1,kgr1,egr2,ggr2,Egam2,Tnuc,
     +        ggredep0,ggredep,enum,denom,factor1,factor2
c
c ************************* Strength functions *************************
c
c fstrength: gamma-ray strength function
c Egamma   : gamma energy
c irad     : variable to indicate M(=0) or E(=1) radiation
c ngr      : number of GR
c Zcomp    : charge number index for compound nucleus
c Ncomp    : neutron number index for compound nucleus
c sgr,sgr1 : strength of GR
c egr,egr1 : energy of GR
c ggr,ggr1 : width of GR
c kgr,kgr1 : constant for gamma-ray strength function
c egr2,ggr2: help variables
c Egam2    : help variable
c strength : strength function of Kopecky-Uhl (1) or Brink-Axel (2) 
c Tnuc     : nuclear temperature
c S        : separation energy per particle
c alev     : level density parameter
c ggredep0 : energy dependent damping width at zero gamma energy
c ggredep  : energy dependent damping width
c twopi    : 2.*pi
c enum     : enumerator of Lorentzian
c denom    : denominator of Lorentzian
c factor1-2: help variables
c egrid    : outgoing energy grid
c ebegin   : first energy point of energy grid  
c
      fstrength=0.
      do 10 i=1,ngr(Zcomp,Ncomp,irad,l)
        sgr1=sgr(Zcomp,Ncomp,irad,l,i)
        egr1=egr(Zcomp,Ncomp,irad,l,i)
        ggr1=ggr(Zcomp,Ncomp,irad,l,i)
        kgr1=kgr(Zcomp,Ncomp,irad,l)
        egr2=egr1**2
        ggr2=ggr1**2
        Egam2=Egamma**2
c
c 1. Kopecky-Uhl generalized Lorentzian.
c
        if (strength.eq.1.and.l.eq.1.and.irad.eq.1) then
          Tnuc=sqrt(max(S(Zcomp,Ncomp,1)-Egamma,0.)/alev(Zcomp,Ncomp))
          ggredep0=ggr1*twopi**2*Tnuc**2/egr2
          ggredep=ggredep0+ggr1*Egam2/egr2
          enum=ggredep*Egamma
          denom=(Egam2-egr2)**2+Egam2*ggredep**2
          factor1=enum/denom
          factor2=0.7*ggredep0/(egr1**3)
          fstrength=fstrength+kgr1*sgr1*ggr1*(factor1+factor2)
        else
c
c 2. Brink-Axel standard Lorentzian.
c
          enum=0.
          if (Egamma.gt.egrid(ebegin(0))) then
            enum=ggr2*Egamma**(3-2*l)
            denom=(Egam2-egr2)**2+Egam2*ggr2
            fstrength=fstrength+kgr1*sgr1*enum/denom
          endif
        endif
   10 continue
      return
      end
Copyright (C) 2004  A.J. Koning, S. Hilaire and M.C. Duijvestijn
