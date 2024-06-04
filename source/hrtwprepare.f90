      subroutine hrtwprepare(tjl,numtjl,st,tav,sv,v,w,numtr,numinc,
     +  WFCfactor)
c
c +---------------------------------------------------------------------
c | Author: Stephane Hilaire
c | Date  : January 24, 2007
c | Task  : Preparation of HRTW width fluctuation correction
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      implicit none
      integer          numtjl,numtr,numinc,niter,i,ni,WFCfactor
      real             sv,v(numtr),w(numtr),alpha,gamma,beta
      double precision tjl(0:5,numtr),st,tav(numtr),st2,factor,t,f,va,
     +                 ft,gt,denom,Ta
c
c ************************** HRTW calculation **************************
c
c We use the model of HRTW (1975), revised in 1980.
c
c Initialization and average transmission coefficients
c
c niter : number of iterations
c numtjl: number of transmission coefficients
c tav   : average transmission coefficients
c tjl   : transmission coefficients
c
      niter=60
      do 10 i=1,numtr
        tav(i)=0.
   10 continue
      do 20 i=1,numtjl+1
        if (tjl(0,i).gt.0.) tav(i)=tjl(1,i)/tjl(0,i)
   20 continue
c
c Calculation of sum over t**2
c
c st2: help variable
c
      st2=0.
      do 30 i=numinc+1,numtjl+1
        st2=st2+tjl(2,i)
   30 continue
c
c Initialisation for HRTW calculation
c
c factor,t,f,sv: help variables
c v,w          : variables for final HRTW calculation
c
c WFCfactor: 1 Original HRTW
c            2 Ernebjerg and Herman, AIP Conference Proceedings
c              Volume 769 Issue 1 Pages 1233-1236,
c              American Institute of Physics, 2005
c
      if (WFCfactor.eq.1) then
        factor=real(4.*st2/(st*st+3.*st2))
        do 40 i=1,numtjl+1
          t=tav(i)
          if (t.lt.st) then
            f=factor*(1.+t/st)
            w(i)=real(1.+2./(1.+t**f)+87.*(((t-st2/st)/st)**2)*
     +        (t/st)**5)
          else
            w(i)=3.
          endif
   40   continue
      else
        alpha=0.139
        beta=15.247
        gamma=4.081
        do 45 i=1,numtjl+1
          Ta=tav(i)
          fT=alpha/(1.-Ta**beta)
          gT=1.+gamma*Ta*(1.-Ta)
          denom=1.+fT*(st**gT)
          va=2.-1./denom
          w(i)=1.+2./va
   45   continue
      endif
      do 50 i=1,numtjl+1
        v(i)=real(tav(i)/(1.+tav(i)/st*(w(i)-1.)))
   50 continue
c
c Loop over iterations
c
c ni: counter
c
      do 60 ni=1,niter
        sv=0.
c
c Sum over v
c
        do 70 i=numinc+1,numtjl+1
          sv=sv+v(i)*tjl(0,i)
   70   continue
c
c Determination of new v's
c
        do 80 i=1,numtjl+1
          v(i)=real(tav(i)/(1.+v(i)/sv*(w(i)-1.)))
   80   continue
   60 continue
      return
      end
Copyright (C)  2013 A.J. Koning, S. Hilaire and S. Goriely
