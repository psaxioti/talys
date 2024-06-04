      subroutine egridastro
c
c +---------------------------------------------------------------------
c | Author  : Stephane Goriely
c | Date    : May 17, 2009
c | Task    : Calculate default incident energy grid for astrophysical 
c |           rate
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      integer          nen,Zix,Nix,neg,neg1,neg2
      real             Temp,dTgrid,Teps
      double precision emin,emax,eg1,eg2,deg1,deg2,de1,de2,deg,e,de,
     +                 Q0,Qn,Qp,Qa,qmin,qmax,am,acm,b1,b2,b3,kt
c
c ******** Set temperature grid for astrophysical calculations *********
c 
c T9       : Temperature grid in 10**9 K
c Temp,Teps: temperatures of basic grid in 10**9 K
c dTgrid   : temperature increment
c
      T9(1)=0.0001
      Temp=0.0001
      dTgrid=0.0004
      nen=1
   10 Temp=Temp+dTgrid
      Teps=Temp+1.e-4
      nen=nen+1
      if (nen.gt.numT) goto 100
      T9(nen)=Temp
      if (Teps.gt.0.0005) dTgrid=0.0005
      if (Teps.gt.0.001) dTgrid=0.004
      if (Teps.gt.0.005) dTgrid=0.005
      if (Teps.gt.0.01) dTgrid=0.04
      if (Teps.gt.0.05) dTgrid=0.05
      if (Teps.gt.0.05) dTgrid=0.05
      if (Teps.gt.0.30) dTgrid=0.1
      if (Teps.gt.1.) dTgrid=0.50
      if (Teps.gt.4.) dTgrid=1.
      if (Teps.gt.10.) goto 100
      goto 10
c
c ************************ incident energy grid ************************
c
c Zix        : charge number index for target nucleus
c Nix        : neutron number index for target nucleus
c redumass,am: reduced mass
c emax,emin  : max and minimum energies
c eg1,eg2    : Gamow energies at T9(min) and T9(max)
c kt         : energy kT expressed in MeV corresponding to a 
c              temperature T9=1
c
  100 kt=0.0862d0
      emin=1.d-12
      if (k0.gt.1) emin=1.d-3
      emax=50.
      neg=50
      neg1=5 
      neg2=10
      Zix=parZ(k0)
      Nix=parN(k0)
      acm=1./specmass(Zix,Nix,k0)
      am=redumass(Zix,Nix,k0)
      Q0=S(0,0,k0)
      Qn=S(0,0,1)
      Qp=S(0,0,2)
      Qa=S(0,0,6)
      qmin=min(Qn,Qp,Qa)
      qmax=max(Qn,Qp,Qa)
      if (k0.eq.0) then
        eg1=max(0.d0,qmin-4.*kt*T9(numT))
        eg2=qmax+4.*kt*T9(numT)
        eg2=max(eg2,20.d0)
        b1=Qn
        b2=Qp
        b3=Qa
      elseif (k0.eq.1) then
        eg1=max(kt*T9(1)/4.,0.001d0)
        eg2=kt*T9(numT)*4.
        b1=-999.
        b2=Qp-Q0
        b3=Qa-Q0
      elseif (k0.gt.1) then
        eg1=0.122*am**(1./3.)*(Ztarget*parZ(k0)*T9(1))**(2./3.)
        deg1=0.237*am**(1./6.)*(Ztarget*parZ(k0))**(1./3.)*
     &    T9(1)**(5./6.)
        eg1=eg1-2.*deg1
        eg2=0.122*am**(1./3.)*(Ztarget*parZ(k0)*T9(numT))**(2./3.)
        deg2=0.237*am**(1./6.)*(Ztarget*parZ(k0))**(1./3.)*
     &    T9(numT)**(5./6.)
        eg2=eg2+2.*deg2
        b1=Qn-Q0
        b2=Qp-Q0
        if (k0.eq.2.and.k0.ne.6) b2=Qa-Q0
        b3=-999.
        eg2=max(eg2,b1,b2)
      endif
      if (eg1.lt.0.) eg1=0.
      nen=0
      e=emin 
      de1=(eg1-emin)/float(neg1)
      deg=(eg2-eg1)/float(neg)
      de2=(emax-eg2)/float(neg2)
  110 continue 
      if (e.gt.emax) goto 200
      if (e.lt.eg1) then
        de=de1
      elseif (e.lt.eg2) then
        de=deg
      else
        de=de2
      endif
      if (k0.eq.1) then
         de=0.000002d00
         if (e.ge.0.00001d00) de=0.000010d00
         if (e.ge.0.0001d00) de=0.00010d00
         if (e.ge.0.001d00) de=0.0010d00
         if (e.ge.0.010d00) de=0.0040d00
         if (e.ge.0.100d00) de=0.0500d00
         if (e.ge.1.000d00) de=0.2500d00
         if (e.ge.5.000d00) de=0.5d00
         if (e.ge.10.00d00) de=2.d00
         if (e.ge.20.00d00) de=5.d00
      endif
      if (abs(e-b1).lt.de.or.abs(e-b2).lt.de.or.
     &  abs(e-b3).lt.de) de=de/10.
      e=e+de
      nen=nen+1
      if (nen.gt.numenin+2) then
        write(*,*) 'Astro-warning: too many energy points'
        nen=numenin+2
        goto 200
      endif
      eninc(nen)=e*acm
      goto 110
  200 numinc=nen 
c
c The minimum and maximum value of all the incident energies is
c determined.
c
c enincmin: minimum incident energy
c enincmax: maximum incident energy
c
      enincmin=eninc(1)
      enincmax=eninc(numinc)
      return
      end
