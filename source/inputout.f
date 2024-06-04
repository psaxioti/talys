      subroutine inputout
c
c +---------------------------------------------------------------------
c | Author: Arjan Koning 
c | Date  : December 2, 2004
c | Task  : Write input parameters
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      character*1  yesno
      character*12 sysstring,rotstring
      integer      i,type
c
c ************************** User input file ***************************
c
c nlines: number of input lines
c inline: input line
c
      write(*,'(/"########## USER INPUT ##########")')
      write(*,'(/"USER INPUT FILE"/)')
      do 10 i=1,nlines
        write(*,'(a80)') inline(i)
   10 continue
c
c ********* All possible input parameters including defaults ***********
c
      write(*,'(/"USER INPUT FILE + DEFAULTS"/)')
      write(*,'("Keyword           Value   Variable     Explanation"/)')
c
c 1. Four main keywords
c
c ptype0    : type of incident particle
c Starget   : symbol of target nucleus
c Atarget   : mass number of target nucleus
c numinc    : number of incident energies
c eninc     : incident energy in MeV
c energyfile: file with incident energies
c
      write(*,'("#"/"# Four main keywords"/"#")')
      write(*,'("projectile          ",a1,"     ptype0",$)') ptype0
      write(*,'("       type of incident particle")') 
      write(*,'("element            ",a2,"     Starget",$)') Starget
      write(*,'("      symbol of target nucleus")')
      write(*,'("mass              ",i3,"     mass        ",$)') Atarget
      write(*,'(" mass number of target nucleus")')
      if (numinc.eq.1) then
        write(*,'("energy            ",f7.3," eninc",$)') eninc(1)
        write(*,'("        incident energy in MeV")') 
      else
        write(*,'("energy ",a14,"     energyfile",$)') energyfile
        write(*,'("   file with incident energies")') 
      endif
c
c 2. Basic physical and numerical parameters
c
c outtype     : type of outgoing particles
c maxZ        : maximal number of protons away from the initial compound
c               nucleus
c maxN        : maximal number of neutrons away from the initial 
c               compound nucleus
c nbins       : number of continuum excitation energy bins
c segment     : number of segments to divide emission energy grid
c nlevmax     : maximum number of included discrete levels for target
c nlevmaxres  : maximum number of included discrete levels for residual
c               nucleus  
c parsym      : symbol of particle    
c nlevbin     : number of excited levels for binary nucleus 
c parname     : name of particle     
c Ltarget     : excited level of target
c isomer      : definition of isomer in seconds
c transpower  : power for transmission coefficient limit
c transeps    : absolute limit for transmission coefficient
c xseps       : limit for cross sections
c popeps      : limit for population cross section per nucleus
c Rfiseps     : ratio for limit for fission cross section per nucleus 
c eninclow    : minimal incident energy for nuclear model calculations
c nangle      : number of angles
c nanglecont  : number of angles for continuum
c nanglerec   : number of recoil angles
c maxenrec    : number of recoil energies 
c yesno       : function to assign y or n to logical value
c flagchannels: flag for exclusive channels calculation
c maxchannel  : maximal number of outgoing particles in individual
c               channel description (e.g. this is 3 for (n,2np))
c flagrel     : flag for relativistic kinematics
c flagrecoil  : flag for calculation of recoils     
c flaglabddx  : flag for calculation of DDX in LAB system
c flagrecoilav: flag for average velocity in recoil calculation
c flagEchannel: flag for channel energy for emission spectrum
c
      write(*,'("#"/"# Basic physical and numerical parameters"/"#")')
      write(*,'("ejectiles ",7(1x,a1),"  outtype      ",$)') 
     +  (outtype(type),type=0,6)
      write(*,'("outgoing particles")')
      write(*,'("maxz              ",i3,"     maxZ         ",$)') maxZ
      write(*,'("maximal number of protons away from the initial",$)')
      write(*,'(" compound nucleus")')
      write(*,'("maxn              ",i3,"     maxN         ",$)') maxN
      write(*,'("maximal number of neutrons away from the initial",$)')
      write(*,'(" compound nucleus")')
      write(*,'("bins              ",i3,"     nbins        ",$)') nbins
      write(*,'("number of continuum excitation energy bins")')
      write(*,'("segment           ",i3,"     segment     ",$)') segment
      write(*,'(" number of segments to divide emission energy grid")')
      write(*,'("maxlevelstar      ",i3,"     nlevmax",$)') nlevmax
      write(*,'("      maximum number of included discrete levels",$)')
      write(*,'(" for target")')
      write(*,'("maxlevelsres      ",i3,"     nlevmaxres",$)') 
     +  nlevmaxres
      write(*,'("   maximum number of included discrete levels",$)')
      write(*,'(" for residual nucleus")')
      do 20 type=0,6
        write(*,'("maxlevelsbin ",a1,"    ",i3,"     nlevbin   ",$)') 
     +    parsym(type),nlevbin(type)
        write(*,'("   maximum number of included discrete levels",$)')
        write(*,'(" for ",a8," channel")') parname(type)
   20 continue
      write(*,'("ltarget           ",i3,"     ltarget",$)') Ltarget
      write(*,'("      excited level of target")')
      write(*,'("isomer          ",1p,e9.2," isomer ",$)') isomer 
      write(*,'("      definition of isomer in seconds")')
      write(*,'("transpower        ",i3,"     transpower",$)') 
     +  transpower
      write(*,'("   power for transmission coefficient limit")')
      write(*,'("transeps        ",1p,e9.2," transeps",$)') transeps
      write(*,'("     limit for transmission coefficient")')
      write(*,'("xseps           ",1p,e9.2," xseps   ",$)') xseps
      write(*,'("     limit for cross sections")')
      write(*,'("popeps          ",1p,e9.2," popeps  ",$)') popeps
      write(*,'("     limit for population cross section per nucleus")')
      write(*,'("Rfiseps         ",1p,e9.2," Rfiseps ",$)') Rfiseps
      write(*,'("     ratio for limit for fission cross section",$)')
      write(*,'(" per nucleus")')
      write(*,'("elow            ",1p,e9.2," elow    ",$)') eninclow
      write(*,'("     minimal incident energy for nuclear model",$)')
      write(*,'(" calculations")')
      write(*,'("angles            ",i3,"     nangle",$)') nangle
      write(*,'("       number of angles")')
      write(*,'("anglescont        ",i3,"     nanglecont",$)') 
     +  nanglecont
      write(*,'("   number of angles for continuum")')
      write(*,'("anglesrec         ",i3,"     nanglerec ",$)') 
     +  nanglerec
      write(*,'("   number of recoil angles")')
      write(*,'("maxenrec          ",i3,"     maxenrec  ",$)') 
     +  maxenrec
      write(*,'("   number of recoil energies")')
      write(*,'("channels            ",a1,"     flagchannels",$)') 
     +  yesno(flagchannels)
      write(*,'(" flag for exclusive channels calculation")')
      write(*,'("maxchannel         ",i2,"     maxchannel",$)') 
     +  maxchannel
      write(*,'("   maximal number of outgoing particles in",$)')
      write(*,'(" individual channel description")')
      write(*,'("relativistic        ",a1,"     flagrel",$)') 
     +  yesno(flagrel)
      write(*,'("      flag for relativistic kinematics")')
      write(*,'("recoil              ",a1,"     flagrecoil",$)') 
     +  yesno(flagrecoil)
      write(*,'("   flag for calculation of recoils")')
      write(*,'("labddx              ",a1,"     flaglabddx",$)') 
     +  yesno(flaglabddx)
      write(*,'("   flag for calculation of DDX in LAB system")')
      write(*,'("recoilaverage       ",a1,"     flagrecoilav",$)') 
     +  yesno(flagrecoilav)
      write(*,'(" flag for average velocity in recoil calculation")')
      write(*,'("channelenergy       ",a1,"     flagEchannel",$)') 
     +  yesno(flagEchannel)
      write(*,'(" flag for channel energy for emission spectrum")')
c
c 3. Optical model
c
c flaglocalomp: flag for local (y) or global (n) optical model
c flagompall  : flag for new optical model calculation for all residual
c               nuclei  
c flagautorot : flag for automatic rotational coupled channels
c               calculations for A > 150
c flagspher   : flag to force spherical optical model   
c flagstate   : flag for optical model potential for each excited state
c maxband     : highest vibrational level added to rotational model 
c maxrot      : number of included excited rotational levels
c sysstring   : help variable
c flagsys     : flag for reaction cross section from systematics
c rotstring   : help variable
c flagrot     : flag for use of rotational optical model per
c               outgoing particle, if available
c core        : even-even core for weakcoupling (-1 or 1)  
c flagecissave: flag for saving ECIS input and output files    
c flageciscalc: flag for new ECIS calculation for outgoing particles
c               and energy grid
c flaginccalc : flag for new ECIS calculation for incident channel
c
      write(*,'("#"/"# Optical model"/"#")')
      write(*,'("localomp            ",a1,"     flaglocalomp",$)') 
     +  yesno(flaglocalomp)
      write(*,'(" flag for local (y) or global (n) optical model")')
      write(*,'("optmodall           ",a1,"     flagompall  ",$)') 
     +  yesno(flagompall)
      write(*,'(" flag for new optical model calculation for ",$)')
      write(*,'("all residual nuclei")')
      write(*,'("autorot             ",a1,"     flagautorot ",$)') 
     +  yesno(flagautorot)
      write(*,'(" flag for automatic rotational coupled channels ",$)')
      write(*,'("calculations for A > 150")')
      write(*,'("spherical           ",a1,"     flagspher   ",$)') 
     +  yesno(flagspher)
      write(*,'(" flag to force spherical optical model")')
      write(*,'("statepot            ",a1,"     flagstate   ",$)') 
     +  yesno(flagstate)
      write(*,'(" flag for optical model potential for each",$)')
      write(*,'(" excited state")')
      write(*,'("maxband           ",i3,"     maxband     ",$)') 
     +  maxband
      write(*,'(" highest vibrational band added to rotational model")')
      write(*,'("maxrot            ",i3,"     maxrot      ",$)') 
     +  maxrot
      write(*,'(" number of included excited rotational levels")')
      sysstring='            '
      i=-1
      do 110 type=1,6
        if (flagsys(type)) then
          i=i+2
          write(sysstring(i:i),'(a1)') parsym(type)
        endif
  110 continue
      write(*,'("sysreaction  ",a12," sysreaction  ",$)') sysstring
      write(*,'("particles with reaction cross section from ",$)')
      write(*,'("systematics")')
      rotstring='            '
      i=-1
      do 120 type=1,6
        if (flagrot(type)) then
          i=i+2
          write(rotstring(i:i),'(a1)') parsym(type)
        endif
  120 continue
      write(*,'("rotational   ",a12," rotational   ",$)') rotstring
      write(*,'("particles with possible rotational optical model")')
      write(*,'("core              ",i3,"     core   ",$)') core
      write(*,'("      even-even core for weakcoupling (-1 or 1)")')
      write(*,'("ecissave            ",a1,"     flagecissave",$)') 
     +  yesno(flagecissave)
      write(*,'(" flag for saving ECIS input and output files")')
      write(*,'("eciscalc            ",a1,"     flageciscalc",$)') 
     +  yesno(flageciscalc)
      write(*,'(" flag for new ECIS calculation for outgoing ",$)')
      write(*,'("particles and energy grid")')
      write(*,'("inccalc             ",a1,"     flaginccalc ",$)') 
     +  yesno(flaginccalc)
      write(*,'(" flag for new ECIS calculation for incident channel")')
c
c 4. Compound nucleus
c
c enincmin    : minimum incident energy
c ewfc        : off-set incident energy for width fluctuation 
c               calculation
c enincmax    : maximum incident energy
c flagwidth   : flag for width fluctuation calculation
c wmode       : designator for width fluctuation model
c flagcomp    : flag for compound nucleus calculation
c flagfullhf  : flag for full spin dependence of transmission
c               coefficients         
c flageciscomp: flag for compound nucleus calculation by ECIS
c
      write(*,'("#"/"# Compound nucleus"/"#")')
      if (numinc.gt.1.and.enincmin.lt.ewfc.and.enincmax.ge.ewfc) then
        write(*,'("widthfluc         ",f7.3," ewfc         ",$)') ewfc
        write(*,'("off-set incident energy for width fluctuation",$)')
        write(*,'(" calculation")')
      else
        write(*,'("widthfluc           ",a1,"     flagwidth  ",$)') 
     +    yesno(flagwidth)
        write(*,'("  flag for width fluctuation calculation")')
      endif
      write(*,'("widthmode          ",i2,"     wmode      ",$)') wmode
      write(*,'("  designator for width fluctuation model")')
      write(*,'("compound            ",a1,"     flagcomp     ",$)') 
     +  yesno(flagcomp)
      write(*,'("flag for compound nucleus model")')
      write(*,'("fullhf              ",a1,"     flagfullhf   ",$)') 
     +  yesno(flagfullhf)
      write(*,'("flag for full spin dependence of transmission ",$)')
      write(*,'("coefficients")')
      write(*,'("eciscompound        ",a1,"     flageciscomp",$)') 
     +  yesno(flageciscomp)
      write(*,'(" flag for compound nucleus calculation by ECIS")')
c
c 5. Gamma emission    
c
c gammax      : number of l-values for gamma multipolarity
c strength    : strength function of Kopecky-Uhl (1) or Brink-Axel (2)
c flagelectron: flag for application of electron conversion coefficient
c
      write(*,'("#"/"# Gamma emission"/"#")')
      write(*,'("gammax             ",i2,"     gammax",$)') gammax
      write(*,'("       number of l-values for gamma multipolarity")')
      write(*,'("strength           ",i2,"     strength",$)') strength
      write(*,'("     strength function of Kopecky-Uhl (1) or ",$)')
      write(*,'("Brink-Axel (2)")')
      write(*,'("electronconv        ",a1,"     flagelectron",$)') 
     +  yesno(flagelectron)
      write(*,'(" flag for application of electron conversion",$)') 
      write(*,'(" coefficient")')
c
c 6. Pre-equilibrium   
c
c epreeq      : on-set incident energy for preequilibrium calculation
c flagpreeq   : flag for pre-equilibrium calculation
c preeqmode   : designator for pre-equilibrium model
c flagmulpre  : flag for multiple pre-equilibrium calculation 
c mpreeqmode  : designator for multiple pre-equilibrium model 
c emulpre     : on-set incident energy for multiple preequilibrium
c flagpespin  : flag for pre-equilibrium spin distribution or compound
c               spin distribution for pre-equilibrium cross section   
c flaggiant0  : flag for collective contribution from giant resonances
c flagsurface : flag for surface effects in exciton model
c flagpecomp  : flag for Kalbach complex particle emission model
c flag2comp   : flag for two-component pre-equilibrium model
c flagecisdwba: flag for new ECIS calculation for DWBA for MSD
c flagonestep : flag for continuum one-step direct only
c
      write(*,'("#"/"# Pre-equilibrium"/"#")')
      if (numinc.gt.1.and.enincmin.lt.epreeq.and.enincmax.ge.epreeq) 
     +  then
        write(*,'("preequilibrium    ",f7.3," epreeq   ",$)') epreeq
        write(*,'("    on-set incident energy for preequilibrium",$)')
        write(*,'(" calculation")')
      else
        write(*,'("preequilibrium      ",a1,"     flagpreeq  ",$)') 
     +    yesno(flagpreeq)
        write(*,'("  flag for pre-equilibrium calculation")')
      endif
      write(*,'("preeqmode          ",i2,"     preeqmode",$)') preeqmode
      write(*,'("    designator for pre-equilibrium model")')
      if (numinc.gt.1.and.enincmin.lt.emulpre.and.enincmax.ge.emulpre) 
     +  then
        write(*,'("multipreeq        ",f7.3," emulpre",$)') emulpre
        write(*,'("      on-set incident energy for multiple",$)')
        write(*,'(" preequilibrium")')
      else
        write(*,'("multipreeq          ",a1,"     flagmulpre   ",$)') 
     +    yesno(flagmulpre)
        write(*,'("flag for multiple pre-equilibrium calculation")')
      endif
      write(*,'("mpreeqmode         ",i2,"     mpreeqmode",$)') 
     +  mpreeqmode
      write(*,'("   designator for multiple pre-equilibrium model")')
      write(*,'("preeqspin           ",a1,"     flagpespin",$)') 
     +  yesno(flagpespin)
      write(*,'("   flag for pre-equilibrium spin distribution")')
      write(*,'("giantresonance      ",a1,"     flaggiant    ",$)') 
     +  yesno(flaggiant0)
      write(*,'("flag for collective contribution from giant",$)')
      write(*,'(" resonances")')
      write(*,'("preeqsurface        ",a1,"     flagsurface",$)') 
     +  yesno(flagsurface)
      write(*,'("  flag for surface effects in exciton model")')
      write(*,'("preeqcomplex        ",a1,"     flagpecomp",$)') 
     +  yesno(flagpecomp)
      write(*,'("   flag for Kalbach complex particle emission model")')
      write(*,'("twocomponent        ",a1,"     flag2comp",$)') 
     +  yesno(flag2comp)
      write(*,'("    flag for two-component pre-equilibrium model")')
      write(*,'("ecisdwba            ",a1,"     flagecisdwba",$)') 
     +  yesno(flagecisdwba)
      write(*,'(" flag for new ECIS calculation for DWBA for MSD")')
      write(*,'("onestep             ",a1,"     flagonestep ",$)') 
     +  yesno(flagonestep)
      write(*,'(" flag for continuum one-step direct only")')
c
c 7. Level densities   
c
c ldmodel     : level density model
c flagasys    : flag for all level density parameters a from 
c               systematics
c flagcolldamp: flag for damping of collective effects
c flaggshell  : flag for energy dependence of single particle level
c               density parameter g

c
      write(*,'("#"/"# Level densities"/"#")')
      write(*,'("ldmodel            ",i2,"     ldmodel    ",$)') ldmodel
      write(*,'("  level density model")')
      write(*,'("asys                ",a1,"     flagasys    ",$)') 
     +  yesno(flagasys)
      write(*,'(" flag for all level density parameters a from ",$)')
      write(*,'("systematics")')
      write(*,'("colldamp            ",a1,"     flagcolldamp",$)') 
     +  yesno(flagcolldamp)
      write(*,'(" flag for damping of collective effects")')
      write(*,'(" flag for coupling of single particle level ",$)')
      write(*,'("density parameter g to level density parameter a")')
      write(*,'("gshell              ",a1,"     flaggshell  ",$)') 
     +  yesno(flaggshell)
      write(*,'(" flag for energy dependence of single particle",$)')
      write(*,'(" level density parameter g")')
c
c 8. Fission           
c
c flagfission: flag for fission
c fismodel   : fission model
c fismodelalt: alternative fission model for default barriers
c flagclass2 : flag for class2 states in fission
c flagmassdis: flag for calculation of fission fragment mass yields
c flagffevap : flag for calculation of particle evaporation from
c              fission fragment mass yields 
c
      write(*,'("#"/"# Fission"/"#")')
      write(*,'("fission             ",a1,"     flagfission",$)') 
     +  yesno(flagfission)
      write(*,'("  flag for fission")')
      write(*,'("fismodel           ",i2,"     fismodel   ",$)') 
     +  fismodel  
      write(*,'("  fission model")')
      write(*,'("fismodelalt        ",i2,"     fismodelalt ",$)') 
     +  fismodelalt
      write(*,'(" alternative fission model for default barriers")')
      write(*,'("class2              ",a1,"     flagclass2 ",$)') 
     +  yesno(flagclass2)
      write(*,'("  flag for class2 states in fission")')
      write(*,'("massdis             ",a1,"     flagmassdis",$)') 
     +  yesno(flagmassdis)
      write(*,'("  flag for calculation of fission fragment",$)')
      write(*,'(" mass yields")')
      write(*,'("ffevaporation       ",a1,"     flagffevap  ",$)') 
     +  yesno(flagffevap)
      write(*,'(" flag for calculation of particle evaporation",$)')
      write(*,'(" from fission fragment mass yields")')
c
c 9. Output
c
c flagmain    : flag for main output
c flagbasic   : flag for output of basic information and results
c flagpop     : flag for output of population
c flagcheck   : flag for output of numerical checks
c flaglevels  : flag for output of discrete level information
c flagdensity : flag for output of level densities
c flagoutomp  : flag for output of optical model parameters
c flagdirect  : flag for output of direct reaction results
c flaginverse : flag for output of transmission coefficients and inverse
c               reaction cross sections
c flagtransen : flag for output of transmission coefficients per energy
c flagoutecis : flag for output of ECIS results
c flaggamma   : flag for output of gamma-ray information
c flagpeout   : flag for output of pre-equilibrium results
c flagfisout  : flag for output of fission information
c flagdisc    : flag for output of discrete state cross sections
c flagspec    : flag for output of spectra
c flagadd     : flag for addition of discrete states to spectra
c flagaddel   : flag for addition of elastic peak to spectra
c flagang     : flag for output of angular distributions
c flaglegendre: flag for output of Legendre coefficients
c ddxmode     : mode for double-differential cross sections: 0: None,
c               1: Angular distributions, 2: Spectra per angle, 3: Both
c flagoutdwba : flag for output of DWBA cross sections for MSD
c flaggamdis  : flag for output of discrete gamma-ray intensities
c flagexc     : flag for output of excitation functions
c flagendf    : flag for information for ENDF-6 file
c flagpartable: flag for output of model parameters on separate file
c
      write(*,'("#"/"# Output"/"#")')
      write(*,'("outmain             ",a1,"     flagmain     ",$)') 
     +  yesno(flagmain)
      write(*,'("flag for main output")')
      write(*,'("outbasic            ",a1,"     flagbasic    ",$)') 
     +  yesno(flagbasic)
      write(*,'("flag for output of basic information and results")')
      write(*,'("outpopulation       ",a1,"     flagpop",$)') 
     +  yesno(flagpop)
      write(*,'("      flag for output of population")')
      write(*,'("outcheck            ",a1,"     flagcheck",$)') 
     +  yesno(flagcheck)
      write(*,'("    flag for output of numerical checks")')
      write(*,'("outlevels           ",a1,"     flaglevels",$)') 
     +  yesno(flaglevels)
      write(*,'("   flag for output of discrete level information")')
      write(*,'("outdensity          ",a1,"     flagdensity",$)') 
     +  yesno(flagdensity)
      write(*,'("  flag for output of level densities")')
      write(*,'("outomp              ",a1,"     flagoutomp ",$)') 
     +  yesno(flagoutomp)
      write(*,'("  flag for output of optical model parameters")')
      write(*,'("outdirect           ",a1,"     flagdirect ",$)') 
     +  yesno(flagdirect)
      write(*,'("  flag for output of direct reaction results")')
      write(*,'("outinverse          ",a1,"     flaginverse",$)') 
     +  yesno(flaginverse)
      write(*,'("  flag for output of transmission coefficients",$)')
      write(*,'(" and inverse reaction cross sections")')
      write(*,'("outtransenergy      ",a1,"     flagtransen",$)') 
     +  yesno(flagtransen)
      write(*,'("  flag for output of transmission coefficients",$)')
      write(*,'(" per energy")')
      write(*,'("outecis             ",a1,"     flagoutecis ",$)') 
     +  yesno(flagoutecis)
      write(*,'(" flag for output of ECIS results")')
      write(*,'("outgamma            ",a1,"     flaggamma",$)') 
     +  yesno(flaggamma)
      write(*,'("    flag for output of gamma-ray information")')
      write(*,'("outpreequilibrium   ",a1,"     flagpeout",$)') 
     +  yesno(flagpeout)
      write(*,'("    flag for output of pre-equilibrium results ")')
      write(*,'("outfission          ",a1,"     flagfisout",$)') 
     +  yesno(flagfisout)
      write(*,'("   flag for output of fission information")')
      write(*,'("outdiscrete         ",a1,"     flagdisc",$)') 
     +  yesno(flagdisc)
      write(*,'("     flag for output of discrete state cross",$)')
      write(*,'(" sections")')
      write(*,'("outspectra          ",a1,"     flagspec",$)') 
     +  yesno(flagspec)
      write(*,'("     flag for output of double-differential cross",$)')
      write(*,'(" sections")')
      write(*,'("adddiscrete         ",a1,"     flagadd ",$)') 
     +  yesno(flagadd)
      write(*,'("     flag for addition of discrete states to",$)')
      write(*,'(" spectra")')
      write(*,'("addelastic          ",a1,"     flagaddel",$)') 
     +  yesno(flagaddel)
      write(*,'("    flag for addition of elastic peak to spectra")')
      write(*,'("outangle            ",a1,"     flagang",$)') 
     +  yesno(flagang)
      write(*,'("      flag for output of angular distributions")')
      write(*,'("outlegendre         ",a1,"     flaglegendre",$)') 
     +  yesno(flaglegendre)
      write(*,'(" flag for output of Legendre coefficients")')
      write(*,'("ddxmode             ",i1,"     ddxmode      ",$)') 
     +  ddxmode
      write(*,'("mode for double-differential cross sections")')
      write(*,'("outdwba             ",a1,"     flagoutdwba  ",$)') 
     +  yesno(flagoutdwba)
      write(*,'("flag for output of DWBA cross sections for MSD")')
      write(*,'("outgamdis           ",a1,"     flaggamdis  ",$)') 
     +  yesno(flaggamdis)
      write(*,'(" flag for output of discrete gamma-ray ",$)')
      write(*,'("intensities")')
      write(*,'("outexcitation       ",a1,"     flagexc",$)') 
     +  yesno(flagexc)
      write(*,'("      flag for output of excitation functions")')
      write(*,'("endf                ",a1,"     flagendf",$)') 
     +  yesno(flagendf)
      write(*,'("     flag for information for ENDF-6 file")')
      write(*,'("partable            ",a1,"     flagpartable",$)') 
     +  yesno(flagpartable)
      write(*,'(" flag for output of model parameters on",$)')
      write(*,'(" separate file")')
      return
      end
Copyright (C) 2004  A.J. Koning, S. Hilaire and M.C. Duijvestijn
