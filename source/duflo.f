      subroutine duflo(nn,nz,ebin)
c
c +---------------------------------------------------------------------
c | Author: Duflo-Zuker (adapted by Arjan Koning)
c | Date  : July 1, 2004
c | Task  : Analytical mass formula of Duflo-Zuker
c +---------------------------------------------------------------------
c
c
c ****************** Declarations and common blocks ********************
c
      real xmag(6),zmag(5),a(21)
c
c ************************ Duflo-Zuker formula *************************
c
c nn  : neutron number
c nz  : proton number 
c ebin: total binding energy
c
      data zmag/14.,28.,50.,82.,114./
      data xmag/14.,28.,50.,82.,126.,184./
      data a/16.178,18.422,120.146,202.305,12.454,0.73598,5.204,1.0645,
     +  1.4206,0.0548,0.1786,.6181,.0988,.0265,-.1537,.3113,-.6650,
     +  -.0553,-.0401,.1774,.4523/
      x=nn
      z=nz
      t=abs(x-z)*.5
      v=x+z
      s=v**(2./3.)
      u=v**(1./3.)
c
c E0: macroscopic part of the binding energy
c
      E0=a(1)*v-a(2)*s-a(3)*t*t/v+a(4)*t*t/u/v-a(5)*t/u-a(6)*z*z/u
      esh=0.
      esh1=0.
      do 10 k=2,5
        f1=zmag(k-1)
        f2=zmag(k)
        dfz=f2-f1
        if (z.ge.f1.and.z.lt.f2) then
          roz=(z-f1)/dfz
          pz=roz*(roz-1)*dfz
          do 20 l=2,6
            f3=xmag(l-1)
            f4=xmag(l)
            dfn=f4-f3
            if (x.ge.f3.and.x.lt.f4) then
              ron=(x-f3)/dfn
              pn=ron*(ron-1)*dfn
              esh=(pn+pz)*a(8)+a(10)*pn*pz
              xx=2.*ron-1.
              zz=2.*roz-1.
              txxx=pn*xx
              tzzz=pz*zz
              txxz=pn*zz
              tzzx=pz*xx
              kl=l-k
              if (kl.eq.0) esh1=a(k+10)*(txxx+tzzz)+a(k+15)*(txxz+tzzx)
              if (kl.eq.1)
     +          esh1=a(k+11)*txxx-a(k+16)*txxz+a(k+10)*tzzz-a(k+15)*tzzx
              if (kl.eq.2)
     +          esh1=a(k+12)*txxx+a(k+17)*txxz+a(k+10)*tzzz+a(k+15)*tzzx
              edef=a(9)*(pn+pz)+a(11)*pn*pz
              if (esh.lt.edef) esh=edef
            endif
   20     continue
        endif
   10 continue
      ebin=e0+esh+esh1
      nn2=nn/2
      nz2=nz/2
      nn2=2*nn2
      nz2=2*nz2
      if (nn2.ne.nn) ebin=ebin-a(7)/u
      if (nz2.ne.nz) ebin=ebin-a(7)/u
      return
      end
