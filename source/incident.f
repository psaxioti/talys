      subroutine incident
c
c +---------------------------------------------------------------------
c | Author: Arjan Koning 
c | Date  : June 12, 2007
c | Task  : Main settings and basic cross sections for incident energy
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      character*1 yesno
c
c ************ Setting of some parameters for incident channel *********
c
c flagpreeq : flag for pre-equilibrium calculation
c flagmulpre: flag for multiple pre-equilibrium calculation 
c flagmain  : flag for main output
c Einc      : incident energy in MeV
c yesno     : function to assign y or n to logical value
c flagwidth : flag for width fluctuation calculation
c
c Multiple pre-equilibrium emission is always set off if there is no
c primary pre-equilibrium emission.
c
      if (.not.flagpreeq) flagmulpre=.false.
c
c Write the energy dependent flags to the output file.
c
      if (flagmain) then
        if (Einc.ge.0.001) then
          write(*,'(/," ########## RESULTS FOR E=",f9.5,
     +      " ##########"/)') Einc       
        else
          write(*,'(/," ########## RESULTS FOR E=",1p,e12.5,
     +      " ##########"/)') Einc        
        endif
        write(*,'(" Energy dependent input flags"/)')
        write(*,'(" Width fluctuations (flagwidth)      : ",a1)') 
     +    yesno(flagwidth)
        write(*,'(" Preequilibrium (flagpreeq)          : ",a1)') 
     +    yesno(flagpreeq)
        write(*,'(" Multiple preequilibrium (flagmulpre): ",a1)') 
     +    yesno(flagmulpre)
      endif
c
c *** Calculation of total, reaction, elastic cross section, 
c     transmission coefficients and elastic angular distribution for 
c     incident energy ***
c
c k0           : index for incident particle  
c coullimit    : energy limit for charged particle OMP calculation
c incidentecis : subroutine for ECIS calculation for incident energy
c incidentread : subroutine to read ECIS results for incident energy
c incidentnorm : subroutine for normalization of reaction cross sections
c                and transmission coefficients for incident channel
c flaginitpop  : flag for initial population distribution
c incidentgamma: subroutine for incident photons
c nin          : counter for incident energy
c numinclow    : number of incident energies below Elow
c spr          : subroutine for S, P and R' resonance parameters
c flaginverse  : flag for output of transmission coefficients and 
c                inverse reaction cross sections
c incidentout  : subroutine for reaction output for incident channel
c
      if (k0.gt.1.and.Einc.lt.coullimit(k0)) return
      if (k0.gt.0) then
        call incidentecis
        call incidentread
        call incidentnorm
      else
        if (.not.flaginitpop) call incidentgamma
      endif
      if (k0.eq.1.and.(Einc.le.0.1.or.nin.eq.numinclow+1)) call spr
      if (flaginverse) call incidentout
      return
      end
Copyright (C) 2004  A.J. Koning, S. Hilaire and M.C. Duijvestijn
