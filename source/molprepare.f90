      subroutine molprepare(tjl,numtjl,st,npmold,xmo,wmo,tav,vnu,
     +  product,numtr,numinc,WFCfactor)
c
c +---------------------------------------------------------------------
c | Author: Stephane Hilaire
c | Date  : September 20, 2022
c | Task  : Preparation of Moldauer width fluctuation correction
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      implicit none
      integer          numtjl,npmold,numtr,numinc,i,im,WFCfactor
      real             xmo(npmold),wmo(npmold),vnu(numtr),
     +                 product(npmold),x,expo,capt,factor,eps,yy,x1,
     +                 prod,fxmsqrt,fxmold,alpha,beta,gamma
      double precision tjl(0:5,numtr),st,tav(numtr),tgamma,Ta,fT,gT,
     +                 denom,alpha1,beta1,gamma1,delta1,f
c
c **************** Preparation of Moldauer integral ********************
c
c tgamma: gamma transmission coefficient
c numtjl: number of transmission coefficients
c tav   : average transmission coefficients
c vnu   : number of degrees of freedom
c tjl   : transmission coefficients
c
c Initialisation and average transmission coefficients
c
      tgamma=0.
      do 10 i=1,numtr
        tav(i)=0.
        vnu(i)=1.
   10 continue
      do 20 i=1,numtjl
        if (tjl(0,i).gt.0.) tav(i)=tjl(1,i)/tjl(0,i)
   20 continue
      tgamma=tjl(1,numtjl+1)
c
c Calculation of number of degrees of freedom
c
c st : denominator of compound nucleus formula
c
c WFCfactor: 1 Original Moldauer, Nucl. Phys., A344, 185(1980)
c            2 Ernebjerg and Herman, AIP Conference Proceedings 
c              Volume 769 Issue 1 Pages 1233-1236, 
c              American Institute of Physics, 2005
c            3 Kawano and Talou, Nuclear Data Sheets 118 (2014) 183–186
c
      if (WFCfactor.eq.1) then
        do 30 i=1,numtjl
          Ta=tav(i)
          vnu(i)=min(1.78+real(((Ta**1.212)-0.78)*exp(-0.228*st)),2.)
   30   continue
      endif
      if (WFCfactor.eq.2) then
        alpha = 0.177
        beta = 20.337
        gamma = 3.148
        do i = 1,numtjl
          Ta = tav(i)
          fT = alpha/(1.-Ta**beta)
          gT = 1.+gamma*Ta*(1.-Ta)
          denom = 1.+fT*(st**gT)
          vnu(i) = min(2.-1./denom,2.)
          vnu(i) = max(vnu(i),1.)
        enddo
      endif
      if (WFCfactor.eq.3) then
        do i = 1,numtjl
          Ta = tav(i)
          alpha1 = 0.0287892*Ta+0.245856
          beta1 = 1.+2.5*Ta*(1.-Ta)*exp(-2.*st)
          gamma1 = Ta*Ta-(st-2.*Ta)**2 
          if (gamma1 > 0. .and. st < 2*Ta) then
            delta1 = sqrt(gamma1)/Ta
          else
            delta1 = 1.
          endif
          f = alpha1*(st+Ta)/(1.-Ta)*beta1*delta1
          denom = 1.+f
          vnu(i) = 2.-1./denom
        enddo
      endif
      vnu(numtjl+1)=1.
c
c Loop over integration points
c
c x          : help variables
c capt       : help variables
c
      do 40 im=1,npmold
        x=xmo(im)
c
c Special case for capture
c
        expo=real(tgamma*x/st)
        if (expo.gt.80.) then
          capt=0.
        else
          capt=exp(-expo)
        endif
c
c Calculation of product over tjl
c
c factor,eps,yy,x1   : help variables
c prod               : help variable
c wmo   : help variable
c fxmold: help variable
c fxmsqrt: help variable
c product            : product used in final Moldauer calculation
c
        factor=0.
        do 50 i=numinc+1,numtjl
          eps=real(2.*tav(i)*x/(st*vnu(i)))
          yy=real(-vnu(i)*0.5*tjl(0,i))
          if (eps.le.1.e-30) goto 50
          if (eps.lt.1.e-5) then
            x1=eps-0.5*eps*eps+(eps**3)/3.-0.25*(eps**4)+0.2*(eps**5)
          else
            x1=log(1.+eps)
          endif
          factor=factor+yy*x1
   50   continue
        prod=exp(factor)
        fxmsqrt=wmo(im)*exp(0.5*x)
        fxmold=fxmsqrt*prod*fxmsqrt
        product(im)=fxmold*capt
   40 continue
      return
      end
Copyright (C)  2013 A.J. Koning, S. Hilaire and S. Goriely
