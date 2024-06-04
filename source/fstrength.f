      function fstrength(Zcomp,Ncomp,Egamma,irad,l)
c
c +---------------------------------------------------------------------
c | Author: Arjan Koning
c | Date  : October 12, 2007
c | Task  : Gamma ray strength functions 
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      integer Zcomp,Ncomp,irad,l,i,nen
      real    fstrength,Egamma,sgr1,egr1,ggr1,kgr1,egr2,ggr2,Egam2,e,
     +        Tnuc,ggredep0,ggredep,enum,denom,factor1,factor2,
     +        Eq(0:numgamqrpa),eb,ee,gamb,game,f1
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
c
c Models for E1 gamma-ray strength function:
c
c 1. Kopecky-Uhl
c 2. Brink-Axel 
c 3. Goriely HFBCS
c 4. Goriely HFB
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
c qrpaexist: flag for existence of tabulated QRPA strength functions
c k0       : index of incident particle
c Einc     : incident energy
c Tnuc     : nuclear temperature
c delta    : energy shift
c S        : separation energy per particle
c alev     : level density parameter
c ggredep0 : energy dependent damping width at zero gamma energy
c ggredep  : energy dependent damping width
c twopi    : 2.*pi
c enum     : enumerator of Lorentzian
c denom    : denominator of Lorentzian
c factor1-2: help variables
c
        if ((strength.eq.1.or.(strength.ge.3.and.
     +    .not.qrpaexist(Zcomp,Ncomp))).and.l.eq.1.and.irad.eq.1) then
          if (k0.gt.0.or.Egamma.ne.Einc) then
            e=max(Einc+S(Zcomp,Ncomp,1)-delta(Zcomp,Ncomp,0)-Egamma,0.)
            Tnuc=sqrt(e/alev(Zcomp,Ncomp))
          else
            Tnuc=0.
          endif
          ggredep0=ggr1*twopi**2*Tnuc**2/egr2
          ggredep=ggredep0+ggr1*Egam2/egr2
          enum=ggredep*Egamma
          denom=(Egam2-egr2)**2+Egam2*ggredep**2
          factor1=enum/denom
          factor2=0.7*ggredep0/(egr1**3)
          fstrength=fstrength+kgr1*sgr1*ggr1*(factor1+factor2)
        endif
c
c 2. Brink-Axel standard Lorentzian.
c
        if (strength.eq.2.or.l.ne.1.or.irad.ne.1) then
          enum=0.
          if (Egamma.gt.0.001) then
            enum=ggr2*Egamma**(3-2*l)
            denom=(Egam2-egr2)**2+Egam2*ggr2
            fstrength=fstrength+kgr1*sgr1*enum/denom
          endif
        endif
c
c 3+4. Tabulated QRPA strength functions
c
c locate    : subroutine to find value in ordered table
c numgamqrpa: number of energies for QRPA strength function
c eqrpa,Eq  : energy grid for QRPA strength function
c fe1qrpa   : tabulated QRPA strength function
c eb,ee,....: help variables
c
        if (strength.ge.3.and.qrpaexist(Zcomp,Ncomp).and.l.eq.1.and.
     +    irad.eq.1) then
          do 110 nen=0,numgamqrpa
            Eq(nen)=eqrpa(Zcomp,Ncomp,nen)
  110     continue
          if (Egamma.le.Eq(numgamqrpa)) then
            call locate(Eq,0,numgamqrpa,Egamma,nen)
            eb=Eq(nen)
            ee=Eq(nen+1)
            gamb=fe1qrpa(Zcomp,Ncomp,nen)
            game=fe1qrpa(Zcomp,Ncomp,nen+1)
            if (gamb.gt.0..and.game.gt.0.) then
              f1=log10(gamb)+(Egamma-eb)/(ee-eb)*
     +          (log10(game)-log10(gamb))
              fstrength=10.**f1
            else
              fstrength=gamb+(Egamma-eb)/(ee-eb)*(game-gamb)
            endif
          else
            eb=Eq(numgamqrpa-1)
            ee=Eq(numgamqrpa)
            gamb=fe1qrpa(Zcomp,Ncomp,numgamqrpa-1)
            game=fe1qrpa(Zcomp,Ncomp,numgamqrpa)
            if (gamb.gt.0..and.game.gt.0.) then
              f1=log10(gamb)+(Egamma-eb)/(ee-eb)*
     +          (log10(game)-log10(gamb))
              fstrength=10.**f1
            else
              fstrength=gamb+(Egamma-eb)/(ee-eb)*(game-gamb)
            endif
          endif
        endif
   10 continue
      return
      end
Copyright (C) 2004  A.J. Koning, S. Hilaire and M.C. Duijvestijn
