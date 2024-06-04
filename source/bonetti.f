      subroutine bonetti
c
c +---------------------------------------------------------------------
c | Author: Marieke Duijvestijn 
c | Date  : July 7, 2004
c | Task  : Determination of effective absorption optical potential
c +---------------------------------------------------------------------
c
c R. Bonetti, S. Galbiati, and L. Milazzo-Colli, Lett. Nuov. Cim. 18, 
c 557 (1977)
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb" 
      real    dr,emax,e,sum1,sum2,radv,radd,rr,fwsvol,fwssurf,term1,
     +        term2
      integer Zix,Nix,nrbins,nenend,k,nen,i
c
c ********** Calculation of volume absorptive optical potential ********
c
c Zix,Zindex : charge number index for residual nucleus
c Nix,Nindex : neutron number index for residual nucleus
c k0         : index of incident particle
c nrbins     : number of continuum excitation energy bins
c dr         : integration bin width 
c emax       : maximum excitation energy
c enincmax   : maximum incident energy
c S          : separation energy
c nenend     : help variable
c e          : energy 
c optical    : subroutine for determination of optical potential 
c radv       : imaginary volume nuclear radius 
c radd       : imaginary surface nuclear radius 
c Atarget    : mass number of target nucleus
c onethird   : 1/3
c w,rw,aw    : imaginary volume potential, radius, diffuseness
c wd,rwd,awd : imaginary surface potential, radius, diffuseness
c rr         : running variable in integration over the radius 
c fwsvol     : imaginary volume optical potential form factor 
c fwssurf    : derivative of the imaginary surface optical potential
c              form factor 
c term1,term2: help variables 
c sum1,sum2  : help variables 
c wvol       : absorption part of the optical potential averaged over
c              the volume
c
      Zix=Zindex(0,0,k0)
      Nix=Nindex(0,0,k0)
      nrbins=50
      dr=20./nrbins
      emax=enincmax+S(0,0,k0)
      nenend=int(emax+1.)
      do 10 k=1,2
        do 20 nen=-200,10*nenend
          e=0.1*real(nen)
          call optical(Zix,Nix,k,e)
          radv=rw*Atarget**onethird
          radd=rwd*Atarget**onethird
          sum1=0.
          sum2=0.
          do 30 i=1,nrbins
            rr=(i-0.5)*dr
            fwsvol=1./(1.+exp((rr-radv)/aw))
            fwssurf=-exp((rr-radd)/awd)/(awd*(1.+exp((rr-radd)/awd))**2)
            term2=fwsvol*(rr**2)*dr
            term1=term2*(w*fwsvol-4.*awd*wd*fwssurf)
            sum1=sum1+term1
            sum2=sum2+term2
 30       continue
          wvol(k,nen)=sum1/sum2
 20     continue
 10   continue
      return
      end
Copyright (C) 2004  A.J. Koning, S. Hilaire and M.C. Duijvestijn
