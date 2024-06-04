      subroutine talysreaction
c
c +---------------------------------------------------------------------
c | Author: Arjan Koning 
c | Date  : October 14, 2004
c | Task  : Reaction models
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
c
c ************************* Reaction calculation ***********************
c
c Basic cross sections (optical model) and initialisation
c
c flagompall  : flag for new optical model calculation for all residual
c               nuclei 
c basicxs     : subroutine for basic cross sections and transmission
c               coefficients        
c enincmax    : maximum incident energy
c epreeq      : on-set incident energy for preequilibrium calculation
c preeqinit   : subroutine for initialization of general 
c               pre-equilibrium parameters
c excitoninit : subroutine for initialization of exciton model 
c               parameters
c flagcomp    : flag for compound nucleus calculation
c compoundinit: subroutine for initialization of compound model 
c               parameters
c
      if (.not.flagompall) call basicxs(0,0)
      if (enincmax.ge.epreeq) then
        call preeqinit
        call excitoninit
      endif
      if (flagcomp) call compoundinit
c
c Loop over incident energies
c
c Initialisation
c
c nin        : counter for incident energies 
c numinc     : number of incident energies 
c eninc,Einc : incident energy in MeV
c energies   : subroutine for energies  
c reacinitial: subroutine for initialization of arrays for various 
c              cross sections
c eninclow   : minimal incident energy for nuclear model calculations
c
      do 10 nin=1,numinc
        if (flagompall) call basicxs(0,0)
        Einc=eninc(nin)
        call energies
        call reacinitial
        if (Einc.lt.eninclow) goto 10
c
c Optical model
c
c incident  : subroutine for main settings and basic cross sections for 
c             incident energy
c exgrid    : subroutine to set excitation energy grid
c flagrecoil: flag for calculation of recoils
c recoilinit: subroutine for calculation of initial recoil velocity
c             and direction
c lmaxinc   : maximal l-value for transmission coefficients
c xsreacinc : reaction cross section for incident channel
c xseps     : limit for cross sections
c
        call incident
        call exgrid(0,0)
        if (flagrecoil) call recoilinit
        if (lmaxinc.eq.-1.or.xsreacinc.lt.xseps) goto 20
c
c Direct reactions
c
c direct: subroutine for calculation of direct inelastic cross sections
c
        call direct
c
c Pre-equilibrium reactions
c
c flagpreeq : flag for pre-equilibrium calculation
c preeq     : subroutine for preequilibrium reactions
c population: subroutine for processing of pre-equilibrium spectra 
c             into population bins
c
        if (flagpreeq) then
          call preeq
          call population
        endif
c
c Binary compound reactions
c
c compnorm  : subroutine for normalization of compound nucleus 
c             cross section
c comptarget: subroutine for compound reaction for initial compound
c             nucleus
c
        if (flagcomp) then
          call compnorm
          call comptarget
        endif
c
c Collecting binary reaction results
c
c binary : subroutine for binary reaction results
c flagang: flag for output of angular distributions
c flagddx: flag for output of double-differential cross sections
c angdis : subroutine for calculation of angular distributions for 
c          discrete states
c
        call binary
        if (flagang.or.flagddx.or.flagrecoil) call angdis
c
c Multiple emission               
c
c multiple: subroutine for multiple emission
c
        call multiple
c
c Exclusive channels
c
c flagchannels: flag for exclusive channels calculation
c channels    : subroutine for exclusive reaction channels
c
        if (flagchannels) call channels
c
c Collecting total cross sections, spectra, angular distributions, etc.
c
c totalxs    : subroutine for total cross sections
c flagspec   : flag for output of spectra
c spectra    : subroutine for creation of particle spectra
c flagfission: flag for fission     
c flagmassdis: flag for calculation of fission fragment mass yields
c massdis    : subroutine for fission fragment yields
c residual   : subroutine for residual production cross sections
c totalrecoil: subroutine for recoil results
c numinclow  : number of incident energies below Elow   
c thermal    : subroutine for estimate of thermal cross sections
c output     : subroutine for output
c
        call totalxs
        if (flagspec.or.flagddx) call spectra
        if (flagfission.and.flagmassdis) call massdis
        call residual
        if (flagrecoil) call totalrecoil
        if (nin.eq.numinclow+1.and.numinclow.gt.0) call thermal
   20   call output
   10 continue
c
c Final output
c
c finalout: subroutine for output of final results
c flagendf: flag for information for ENDF-6 file
c endf    : subroutine for cross sections and information for 
c           ENDF-6 file
c flagmain: flag for main output
c timer   : subroutine for output of execution time
c
      call finalout
      if (flagendf) call endf
      if (flagmain) call timer
      return
      end
Copyright (C) 2004  A.J. Koning, S. Hilaire and M.C. Duijvestijn
