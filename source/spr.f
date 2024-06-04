      subroutine spr
c
c +---------------------------------------------------------------------
c | Author: Arjan Koning
c | Date  : January 13, 2009
c | Task  : S, P and R' resonance parameters
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      integer l,J,odd
      real    Efac,Rpot,r2k2,gamjl,gamgamjl
c
c **************** S, P and R' resonance parameters ********************
c
c Rprime     : potential scattering radius
c xselasinc  : total elastic cross section (neutrons only) for incident
c              channel
c fourpi     : 4.*pi
c Efac       : factor with incident energy
c Einc       : incident energy in MeV
c twopi      : 2.*pi
c Rpot       : standard value for R
c Atarget    : mass number of target nucleus
c onethird   : 1/3
c Sstrength  : s,p,d,etc-wave strength function
c Tlinc      : transmission coefficients as a function of l for the 
c              incident channel, averaged over spin
c r2k2       : help variable
c wavenum    : wave number
c flagurr    : flag for output of unresolved resonance parameters
c numinclow  : number of incident energies below Elow
c Ztarget    : charge number of target nucleus    
c nuc        : symbol of nucleus
c jdis       : spin of level
c xscaptherm : thermal capture cross section
c S          : neutron separation energy
c odd        : odd (1) or even (0) nucleus
c dtheory    : subroutine for theoretical calculation of average neutron
c              spacings
c numJ       : number of J-values
c Djltheo    : theoretical resonance spacing per J,l value
c gamjl      : average neutron width per J,l value
c Dltheo     : theoretical resonance spacing per l value
c gamgamjl   : experimental total radiative width in eV per J,l value
c strength   : model for E1 gamma-ray strength function
c gnorm      : gamma normalization factor
c numinc     : number of incident energies
c flagendf   : flag for information for ENDF-6 file
c flagendfdet: flag for detailed ENDF-6 information per channel
c nin        : counter for incident energy
c
      Rprime=10.*sqrt(0.001*xselasinc/fourpi)
      Efac=1./(sqrt(1.e6*Einc)*twopi)
      Rpot=1.35*Atarget**onethird
      r2k2=Rpot*Rpot*wavenum*wavenum
      Sstrength(0)=Tlinc(0)*Efac
      Sstrength(1)=Tlinc(1)*Efac*(1.+r2k2)/r2k2
      Sstrength(2)=Tlinc(2)*Efac*(9.+3.*r2k2+r2k2*r2k2)/(r2k2*r2k2)
      if (flagurr) then
        if (nin.eq.numinclow+1) then
          open (unit=21,status='unknown',file='urr.dat')
          write(21,'("#")')
          write(21,'("# Resonance parameters for Z=",i4," A=",i4," (",
     +      i3,a2,") Target spin=",f4.1)') Ztarget,Atarget,Atarget,
     +      nuc(Ztarget),jdis(0,1,0)
          write(21,'("# Thermal capture cross section=",1pe12.5," mb",
     +      "   Sn=",e12.5," MeV")') xscaptherm,S(0,0,1)
          write(21,'("#")')
        endif
        write(21,'("#  Einc[MeV]=",1pe12.5)') Einc
        write(21,'("# Rprime[fm]=",1pe12.5)') Rprime
        write(21,'("#   J  l    D(J,l)[eV]    D(l)[eV]     S(l)",
     +    "    Gn(J,l)[eV] Gg(J,l)[eV]")')
        odd=mod(Atarget+1,2)
        do 10 l=0,2
          call dtheory(0,0,l,Einc)
          do 20 J=0,numJ
            if (Djltheo(J,l).eq.0.) goto 20
            gamjl=Djltheo(J,l)*Sstrength(l)
            if (strength.eq.1) then
              gamgamjl=Djltheo(J,l)*Sgam(J,l)*gnorm
            else
              gamgamjl=Djltheo(J,l)*Sgam(J,l)
            endif
            write(21,'(f5.1,i3,2x,1p,5e12.5)') J+0.5*odd,l,Djltheo(J,l),
     +        Dltheo(l),Sstrength(l),gamjl,gamgamjl
   20     continue
   10   continue
        if (nin.eq.numinc) close (unit=21)
      endif
      if (flagendf.and.flagendfdet.and.
     +  (Einc.le.0.1.or.nin.eq.numinclow+1)) then
        open (unit=1,status='unknown',file='spr.opt')
        write(1,'(2i4,3f8.4)') Atarget,Ztarget,Sstrength(0)*1.e4,
     +    Sstrength(1)*1.e4,Rprime
        close (unit=1)
      endif
      return
      end
Copyright (C) 2004  A.J. Koning, S. Hilaire and M.C. Duijvestijn
