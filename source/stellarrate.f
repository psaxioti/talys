      subroutine stellarrate
c
c +---------------------------------------------------------------------
c | Author  : Stephane Goriely and Stephane Hilaire
c | Date    : May 19, 2009
c | Task    : Calculate reaction rate for a Maxwell-Boltzmann 
c |           distribution
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      integer          Zix,Nix,Ncomp,Zcomp,nloop,iloop,nen,i
      real             fact0,am,e,fex,ee0,ee1
      double precision fact,sum,xs,term,dE,rate
c
c ************************ Partition Function **************************
c
c fact0      : reaction rate factor for photons
c Zix        : charge number index for target nucleus
c parZ       : charge number of particle
c k0         : index of incident particle
c Nix        : neutron number index for target nucleus
c parN       : neutron number of particle
c redumass,am: reduced mass
c partfunc   : subroutine to calculate partition function
c
c Initialization
c
      fact0=3.951776e+17
      Zix=parZ(k0)
      Nix=parN(k0)
      am=redumass(Zix,Nix,k0)
c
c Calculation of the T-dependent partition function
c
      call partfunc
c
c Calculation of thermonuclear reaction rate factor
c
c Ncomp       : neutron number index for compound nucleus
c maxN        : maximal number of neutron away from initial compound 
c               nucleus
c Zcomp       : proton number index for compound nucleus
c maxZ        : maximal number of protons away from initial compound 
c               nucleus
c flagfission : flag for fission
c numT        : number of temperatures 
c T9          : Temperature grid in 10**9 K
c fact        : reaction rate factor for particles
c sum,fex     : help variables
c numinc      : number of incident energies
c eninc,e     : incident energy in MeV
c xs          : cross section expressed in barn 
c xsastro     : cross section for astrophysical calculation
c xsastrofis  : astrophysical fission cross section
c term        : help variable
c dE          : incident energy bin
c flagastrogs : flag for calculation of astrophysics reaction rate with 
c               target in ground state only
c rateastro   : thermonuclear reaction rate factor
c rateastrofis: thermonuclear reaction rate factor for fission
c partf       : integrated partition function
c
      do 10 Ncomp=0,maxN
        do 20 Zcomp=0,maxZ
          nloop=1
          if (flagfission.and.Ncomp.eq.0.and.Zcomp.eq.0) nloop=2
          do 30 iloop=1,nloop
            do 40 i=1,numT
              rateastrofis(i)=0.
              fact=3.7335e+10/(sqrt(am*(T9(i)**3)))
              sum=0.
              do 50 nen=1,numinc
                e=real(eninc(nen)*specmass(Zix,Nix,k0))
                if (nen.lt.numinc) then
                  ee1=real(eninc(nen+1)*specmass(Zix,Nix,k0))
                else
                  ee1=e
                endif
                if (nen.gt.1) then
                  ee0=real(eninc(nen-1)*specmass(Zix,Nix,k0))
                else
                  ee0=e
                endif
                if (iloop.eq.1) then
                  xs=xsastro(Zcomp,Ncomp,nen)/1000.
                else
                  xs=xsastrofis(nen)/1000.
                endif
                fex=11.605*e/T9(i)
                if (fex.gt.80.) goto 50
                if (k0.ne.0) then
                  term=fact*e*xs*exp(-fex)
                else
                  term=fact0*e**2*xs/(exp(fex)-1.)
                endif
                dE=(ee1-ee0)/2.
                sum=sum+term*dE
  50          continue
              if (flagastrogs) then
                rate=sum
              else
                rate=sum/partf(i)
              endif
              if (iloop.eq.1) then
                rateastro(Zcomp,Ncomp,i)=rate
              else
                rateastrofis(i)=rate
              endif
  40        continue
  30      continue
  20    continue
  10  continue
      return
      end
Copyright (C) 2005  A.J. Koning, S. Hilaire and M.C. Duijvestijn
