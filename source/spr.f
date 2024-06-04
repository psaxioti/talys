      subroutine spr
c
c +---------------------------------------------------------------------
c | Author: Arjan Koning
c | Date  : October 19, 2006
c | Task  : S, P and R' resonance parameters
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      real Rpot,r2k2
c
c **************** S, P and R' resonance parameters ********************
c
c swave      : S-wave strength function
c Tlinc      : transmission coefficients as a function of l for the 
c              incident channel, averaged over spin
c Einc       : incident energy in MeV
c twopi      : 2.*pi
c Rpot       : standard value for R
c Atarget    : mass number of target nucleus
c onethird   : 1/3
c wavenum    : wave number
c r2k2       : help variable
c pwave      : P-wave strength function
c Rprime     : potential scattering radius
c xselasinc  : total elastic cross section (neutrons only) for incident
c              channel
c fourpi     : 4.*pi
c flagendf   : flag for information for ENDF-6 file
c flagendfdet: flag for detailed ENDF-6 information per channel
c Ztarget    : charge number of target nucleus    
c
      swave=Tlinc(0)/(sqrt(1.e6*Einc)*twopi)
      Rpot=1.35*Atarget**onethird
      r2k2=Rpot*Rpot*wavenum*wavenum
      pwave=Tlinc(1)/(sqrt(1.e6*Einc)*twopi)*(1.+r2k2)/r2k2
      Rprime=10.*sqrt(0.001*xselasinc/fourpi)
      if (flagendf.and.flagendfdet) then
        open (unit=1,status='unknown',file='spr.opt')
        write(1,'(2i4,3f8.4)') Atarget,Ztarget,swave*1.e4,pwave*1.e4,
     +    Rprime
        close (unit=1)
      endif
      return
      end
Copyright (C) 2004  A.J. Koning, S. Hilaire and M.C. Duijvestijn
