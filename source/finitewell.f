      function finitewell(p,h,Eex,Ewell,surfwell)
c
c +---------------------------------------------------------------------
c | Author: Marieke Duijvestijn and Arjan Koning
c | Date  : September 2, 2004
c | Task  : Correction for finite well depth
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      logical surfwell
      integer p,h,n,jbin,jwell,k
      real    finitewell,Eex,Ewell,widthdis,wtsum,fwtsum,Ewellj,x,wtinv,
     +        wt,fwell,factor
c
c ********************** Finite well correction ************************
c
c finitewell: correction function for finite well depth
c p         : particle number
c h         : hole number
c Eex       : excitation energy
c Ewell     : average depth of potential well
c surfwell  : flag for surface effects in finite well  
c n         : exciton number
c widthdis  : width of the weighting distribution for averaging the 
c             finite well depth correction function 
c Efermi    : maximal depth of Fermi well 
c jbin      : number of bins used in averaging procedure equals 2*jbin+1
c wtsum     : sum over all contributions to the weight
c fwtsum    : sum over all weighted contributions to the finite well 
c             depth correction function
c jwell     : help variable in loop over potential well depth bins
c Ewellj    : potential well depth used in loop
c x         : help variable
c wtinv,wt  : (inverse) weight
c fwell     : help variable
c factor    : help variable
c sgn       : +1 for even argument, -1 for odd argument
c ncomb     : n!/(n!(n-k)!)
c
c The general finite well correction comes from Betak and Dobes, 
c Z. Phys. A279 (1976) 319. For the first p-h interaction, surface
c effects can be included according to Kalbach, Phys Rev C32, 
c p. 1157 (1985).
c
      finitewell=1.
      n=p+h
c
c Kalbach surface effects
c
      if (surfwell.and.Ewell.lt.(Efermi-0.5)) then
        widthdis=Ewell*(Efermi-Ewell)/(2.*Efermi)
        if (p.eq.0.and.h.eq.1) then
          if (Eex.gt.1.16*Efermi) then
            finitewell=0.
            return
          endif
          if (Eex.lt.Ewell) then
            finitewell=1.
            return
          endif
          finitewell=1./(1.+exp((Eex-Ewell)/widthdis))
        else
          jbin=4
          wtsum=0.
          fwtsum=0.
          do 20 jwell=-jbin,jbin
            Ewellj=Ewell+jwell*widthdis
            x=(Ewellj-Ewell)/widthdis
            if (x.le.80.) then
              wtinv=(1.+exp(x))*(1.+exp(-x))
              wt=1./wtinv
            else
              wt=0.
            endif
            if (Ewellj.gt.Efermi.or.Ewellj.lt.0.) goto 20
            fwell=0.
            do 10 k=0,h
              if (Eex-k*Ewellj.gt.0.) fwell=fwell+
     +          sgn(k)*ncomb(h,k)*((Eex-k*Ewellj)/Eex)**(n-1)
 10         continue
            wtsum=wtsum+wt
            fwtsum=fwtsum+wt*fwell
 20       continue
          if (wtsum.gt.0.) finitewell=fwtsum/wtsum  
        endif
      else
c
c Betak and Dobes surface effects
c
        if (Eex.le.Ewell) return
        if (h.eq.1.and.n.eq.1) return
        do 30 k=1,h
          if (Eex-k*Ewell.gt.0.) then
            factor=((Eex-k*Ewell)/Eex)**(n-1)
            finitewell=finitewell+sgn(k)*ncomb(h,k)*factor
          endif
 30     continue
      endif
      return
      end
Copyright (C) 2004  A.J. Koning, S. Hilaire and M.C. Duijvestijn
