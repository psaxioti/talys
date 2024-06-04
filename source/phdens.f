      function phdens(p,h,gs,Eex,Ewell,surfwell)
c
c +---------------------------------------------------------------------
c | Author: Arjan Koning
c | Date  : August 9, 2004
c | Task  : Particle-hole state density
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      logical surfwell
      integer p,h,n1
      real    phdens,gs,Eex,Ewell,Ap,fac1,factor,finitewell
c
c ***** State density of Betak and Dobes, Z. Phys. A279 (1976) 319. ****
c
c phdens    : particle-hole state density
c p         : particle number
c h         : hole number
c gs        : single-particle level density parameter
c Eex       : excitation energy
c Ewell     : depth of potential well
c surfwell  : flag for surface effects in finite well
c Apauli,Ap : Pauli blocking correction energy
c n1        : exciton number-1
c fac1      : help variable
c factor    : help variable
c nfac      : n!
c finitewell: function for correction for finite well depth 
c
c In general, the finite depth of the hole is included. If the 
c uncorrected state density is required, it should be specified by 
c Ewell=0. (which thus actually means Ewell=infinity).
c
      phdens=0.
      if (p.lt.0.or.h.lt.0) return
      if (p+h.eq.0) return
      Ap=Apauli(p,h)
      if (Ap.ge.Eex) return
      n1=p+h-1
      fac1=nfac(p)*nfac(h)*nfac(n1)
      factor=gs**(p+h)/fac1
      phdens=factor*(Eex-Ap)**n1
      phdens=phdens*finitewell(p,h,Eex,Ewell,surfwell) 
      if (phdens.lt.1.e-10) phdens=0.
      return
      end
Copyright (C) 2004  A.J. Koning, S. Hilaire and M.C. Duijvestijn
