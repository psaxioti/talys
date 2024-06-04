      subroutine input5
c
c +---------------------------------------------------------------------
c | Author: Arjan Koning
c | Date  : January 24, 2023
c | Task  : Read input for fifth set of variables
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      logical      flagassign,wtadjust
      character*1  ch
      character*132 word(40),key,value,cval,line
      integer      Zix,Nix,ibar,irad,Z,N,oddZ,oddN,lval,type,mt,is,igr,
     +             i,class,iz,ia,ival,type2,omptype,nr,k,iword,ilev0,
     +             nbr,ilev1,istat,A,Aact,l
      real         alphaldall,betaldall,gammashell1all,br,sum,val,
     +             Pshiftconstantall,Vf0,rfiscor
c
c ************** Defaults for fifth set of input variables *************
c
c k0            : index for incident particle
c flagendf      : flag for information for ENDF-6 file
c eninclow      : minimal incident energy for nuclear model calculations
c Atarget       : mass number of target nucleus
c strength      : model for E1 gamma-ray strength function
c flagffruns    : flag to denote that run is for fission fragment
c Rspincut      : adjustable constant (global) for spin cutoff factor
c spincutmodel  : model for spin cutoff factor for ground state
c shellmodel    : model for shell correction energies
c kvibmodel     : model for vibrational enhancement
c Cnubar1       : adjustable parameter for nubar constant value
c Cnubar2       : adjustable parameter for nubar energy slope
c Tmadjust      : adjustable parameter for PFNS temperature
c Fsadjust      : adjustable parameter for PFNS scission fraction
c gammashell2   : gamma-constant for asymptotic level density parameter
c flagcol       : flag for collective enhancement of level density
c flagcolldamp  : flag for damping of collective effects in effective
c                 level density (without explicit collective
c                 enhancement). Only used for Bruyeres-le-Chatel
c                 (Pascal Romain) fission model
c pairconstant  : constant for pairing energy systematics
c Ufermi        : energy of Fermi distribution for damping of
c                 rotational effects
c cfermi        : width of Fermi distribution for damping of
c                 rotational effects
c Kph           : constant for single-particle level density parameter
c                 (g=A/Kph)
c M2constant    : overall constant for matrix element in exciton model
c M2limit       : constant for asymptotic value for matrix element
c M2shift       : constant for energy shift for matrix element
c Rspincutpreeq : adjustable constant (global) for preequilibrium spin 
c                 cutoff factor
c Rpinu,Rnupi.  : ratio for two-component matrix element
c Esurf0        : well depth for surface interaction
c Rgamma        : adjustable parameter for pre-equilibrium gamma decay
c elwidth       : width of elastic peak in MeV
c xscaptherm    : thermal capture cross section
c xsptherm      : thermal (n,p) cross section
c xsalphatherm  : thermal (n,a) cross section
c Zix           : charge number index for residual nucleus
c numZ          : maximal number of protons away from the initial
c                 compound nucleus
c Nix           : neutron number index for residual nucleus
c numN          : maximal number of neutrons away from the initial
c                 compound nucleus
c ldmodel       : level density model
c alphad        : alpha-constant for asymptotic level density parameter
c betald        : beta-constant for asymptotic level density parameter
c gammashell1   : gamma-constant for asymptotic level density parameter
c Pshiftconstant: global constant for pairing shift
c alev          : level density parameter
c alimit        : asymptotic level density parameter
c gammald       : gamma-constant for asymptotic level density parameter
c pair          : total pairing correction
c ibar          : fission barrier
c numbar        : number of fission barriers
c Pshift        : adjustable pairing shift
c Pshiftadjust  : adjustable correction to pairing shift
c deltaW        : shell correction in nuclear mass
c Exmatch       : matching point for Ex
c T             : nuclear temperature
c E0            : constant of temperature formula
c Nlow          : lowest discrete level for temperature matching
c Ntop          : highest discrete level for temperature matching
c Krotconstant  : normalization constant for rotational enhancement
c s2adjust      : adjustable constant (Z,A,barrier-dependent) for spin
c                 cutoff parameter
c ctable,ptable : constant to adjust tabulated level densities
c ctableadjust  : adjustable correction to ctable
c ptableadjust  : adjustable correction to ptable
c cglobal       : global constant to adjust tabulated level densities
c pglobal       : global constant to adjust tabulated level densities
c g             : single-particle level density parameter
c gp            : single-particle proton level density parameter
c gn            : single-particle neutron level density parameter
c gamgam        : total radiative width in eV
c D0            : experimental s-wave resonance spacing in eV
c Risomer       : adjustable correction to level branching ratios
c etable,ftable : constant to adjust tabulated strength functions
c etableadjust..: adjustable correction to tabulated strength functions
c ldadjust      : logical for energy-dependent level density adjustment
c gamadjust     : logical for energy-dependent gamma adjustment
c fisadjust     : logical for energy-dependent fission adjustment
c irad          : variable to indicate M(=0) or E(=1) radiation
c lval          : multipolarity
c numgam        : maximum number of l-values for gamma multipolarity
c igr           : giant resonance
c egr           : energy of GR
c ggr           : width of GR
c sgr           : strength of GR
c epr           : energy of PR
c gpr           : width of PR
c tpr           : strength of PR
c egradjust.....: adjustable factors for giant resonance parameters
c                 (default 1.)
c upbend        : properties of the low-energy upbend of given multipolarity
c fiso          : correction factor for isospin forbidden transitions
c fisom         : correction factor for isospin forbidden transitions
c                 in multiple emission
c fisominit     : initial value for fisom
c aadjust,...   : adjustable factors for level density parameters
c                 (default 1.)
c axtype        : type of axiality of barrier
c                 1: axial symmetry
c                 2: left-right asymmetry
c                 3: triaxial and left-right symmetry
c                 4: triaxial no left-right symmetry
c                 5: no symmetry
c fbarrier      : height of fission barrier
c fwidth        : width of fission barrier
c bdamp         : fission partial damping parameter
c fbaradjust,...: adjustable factors for fission parameters
c                 (default 1.)
c betafiscor    : adjustable factor for fission path width
c Zinit         : charge number of initial compound nucleus
c Ninit         : neutron number of initial compound nucleus
c Aact          : mass number for actinide
c Vf0           : vfiscor at mass 240
c vfiscoradjust : adjustable correction to vfiscor
c betafiscoradjust: adjustable correction to betafiscor
c rfiscor       : slope of mass dependence
c vfiscor       : adjustable factor for fission path height
c                 transition states
c Rclass2mom    : normalization constant for moment of inertia for
c                 class 2 states
c widthc2       : width of class2 states
c Ninit         : neutron number of initial compound nucleus
c hbtransfile   : file with head band transition states
c class2file    : file with class 2 transition states
c levelfile     : discrete level file
c deformfile    : deformation parameter file
c Exlfile       : tabulated strength function file
c densfile      : tabulated level density file
c optmodfileN   : optical model parameter file for neutrons
c optmodfileP   : optical model parameter file for protons
c radialfile    : radial matter density file
c ompenergyfile : file with energies for OMP calculation (ENDF files
c                 only)
c yieldfile     : file with fission fragment yields
c radialmodel   : model for radial matter densities (JLM OMP only)
c breakupmodel  : model for break-up reaction: 1. Kalbach 2. Avrigeanu
c massnucleus   : mass of nucleus in amu as read from user input file
c massexcess    : mass excess in MeV as read from user input file
c massdir       : directory with mass tables
c beta2         : deformation parameter
c msdbins       : number of energy points for DWBA calculation for MSD
c Emsdmin       : minimal outgoing energy for MSD calculation
c Cstrip        : adjustable parameter for stripping/pick-up reactions
c Cknock        : adjustable parameter for knockout reactions
c Cbreak        : adjustable parameter for breakup reactions
c nbranch       : number of branching levels
c branchlevel   : level to which branching takes place
c branchratio   : gamma-ray branching ratio to level
c v1adjust....  : adjustable factors for OMP (default 1.)
c preeqadjust   : logical for energy-dependent pre-eq. adjustment
c Nadjust       : number of adjustable parameters
c adjustkey     : keyword for local adjustment
c adjustfile    : file for local adjustment
c adjustpar     : local adjustment parameters
c adjustix      : local adjustment index
c nenadjust     : number of tabulated energies of local adjustment
c Eadjust       : tabulated energy of local adjustment
c Dadjust       : tabulated depth of local adjustment
c ompadjustF    : logical for local OMP adjustment
c ompadjustN    : number of energy ranges for local OMP adjustment
c ompadjustE1   : start energy of local OMP adjustment
c ompadjustE2   : end energy of local OMP adjustment
c ompadjustD    : depth of local OMP adjustment
c ompadjusts    : variance of local OMP adjustment
c ompadjustp    : logical for energy-dependent OMP adjustment
c Ejoin         : joining energy for high energy OMP
c Vinfadjust    : adjustable factor for high energy limit of
c                 real central potential
c adjustTJ      : logical for energy-dependent Tlj adjustment
c jlmmode       : option for JLM imaginary potential normalization
c flagrescue    : flag for final rescue: normalization to data
c rescuefile    : file with incident energy dependent adjustment factors
c grescue       : global multiplication factor for incident energy
c                 dependent adjustment factors
c flaglabddx    : flag for calculation of DDX in LAB system
c aradialcor    : adjustable parameter for shape of DF alpha potential
c adepthcor     : adjustable parameter for depth of DF alpha potential
c nanglerec     : number of recoil angles
c RprimeU       : potential scattering radius
c nTmax         : effective number of temperatures for Maxwellian
c                 average
c astroT9       : temperature, in 10^9 K, for Maxwellian average
c astroE        : energy, in MeV, for Maxwellian average
c Ebeam         : incident energy in MeV for isotope production
c Eback         : lower end of energy range in MeV for isotope
c                 production
c radiounit     : unit for radioactivity: Bq, kBq, MBq, Gbq,
c                 mCi, Ci or kCi
c yieldunit     : unit for isotope yield: num (number),
c                 mug (micro-gram), mg, g, or kg
c Ibeam         : beam current in mA for isotope production
c Area          : target area in cm^2 for isotope production
c Tirrad        : irradiation time per unit
c unitTirrad    : irradiation time unit (y,d,h,m,s)
c Tcool         : cooling time per unit
c unitTcool     : cooling time unit (y,d,h,m,s)
c rhotarget     : target material density
c flagriplomp   : flag for RIPL OMP
c riplomp       : RIPL OMP number
c
      if (k0.eq.1.and.flagendf) then
        eninclow=0.
      else
        eninclow=1.e-6
      endif
c
c Advice of S. Goriely: no gamma normalization for A < 40.
c
      Rspincut=1.
      spincutmodel=1
      shellmodel=1
      kvibmodel=2
      gammashell2=0.
      pairconstant=12.
      Kph=15.
      M2constant=1.
      M2limit=1.
      M2shift=1.
      Rpipi=1.
      Rnunu=1.5
      Rpinu=1.
      Rnupi=1.
      Rspincutpreeq=1.
      Esurf0=-1.
      Rgamma=2.
      elwidth=0.5
      xscaptherm=0.
      xsptherm=0.
      xsalphatherm=0.
      alphaldall=-99.
      betaldall=-99.
      gammashell1all=-99.
      Pshiftconstantall=-99.
      Cnubar1=1.
      Cnubar2=1.
      Tmadjust=1.
      Fsadjust=1.
      if (Atarget <= fislim) then
        Cbarrier = 0.85
      else
        Cbarrier = 1.20
      endif
      do Zix=0,numZ
        do Nix=0,numN
          if (ldmodel(Zix,Nix).eq.1.or.ldmodel(Zix,Nix).ge.4) then
            if (flagcol(Zix,Nix)) then
              alphald(Zix,Nix)=0.0207305
              betald(Zix,Nix)=0.229537
              gammashell1(Zix,Nix)=0.473625
              Pshiftconstant(Zix,Nix)=0.
            else
              alphald(Zix,Nix)=0.0692559
              betald(Zix,Nix)=0.282769
              gammashell1(Zix,Nix)=0.433090
              Pshiftconstant(Zix,Nix)=0.
            endif
            if (flagcolldamp) then
              alphald(Zix,Nix)=0.0666
              betald(Zix,Nix)=0.258
              gammashell1(Zix,Nix)=0.459
              Pshiftconstant(Zix,Nix)=0.
            endif
          endif
          if (ldmodel(Zix,Nix).eq.2) then
            if (flagcol(Zix,Nix)) then
              alphald(Zix,Nix)=0.0381563
              betald(Zix,Nix)=0.105378
              gammashell1(Zix,Nix)=0.546474
              Pshiftconstant(Zix,Nix)=0.743229
            else
              alphald(Zix,Nix)=0.0722396
              betald(Zix,Nix)=0.195267
              gammashell1(Zix,Nix)=0.410289
              Pshiftconstant(Zix,Nix)=0.173015
            endif
          endif
          if (ldmodel(Zix,Nix).eq.3) then
            if (flagcol(Zix,Nix)) then
              alphald(Zix,Nix)=0.0357750
              betald(Zix,Nix)=0.135307
              gammashell1(Zix,Nix)=0.699663
              Pshiftconstant(Zix,Nix)=-0.149106
            else
              alphald(Zix,Nix)=0.110575
              betald(Zix,Nix)=0.0313662
              gammashell1(Zix,Nix)=0.648723
              Pshiftconstant(Zix,Nix)=1.13208
            endif
          endif
        enddo
      enddo
      alev=0.
      alimit=0.
      gammald=-1.
      pair=1.e-20
      Ufermi=30.
      cfermi=5.
      Pshift=1.e-20
      Pshiftadjust=0.
      deltaW=0.
      Exmatch=0.
      T=0.
      E0=1.e-20
      Tadjust=1.
      E0adjust=1.
      Exmatchadjust=1.
      Nlow=-1
      Ntop=-1
      s2adjust=1.
      Krotconstant=1.
      ctable=cglobal
      ptable=pglobal
      ctableadjust=0.
      ptableadjust=0.
      g=0.
      gp=0.
      gn=0.
      gamgam=0.
      D0=0.
      Risomer=1.
      ldadjust=.false.
      gamadjust=.false.
      fisadjust=.false.
      adjustTJ=.false.
      TJadjust=1.
      etable=0.
      ftable=1.
      wtable=1.
      wtadjust=.false.
      etableadjust=0.
      ftableadjust=1.
      wtableadjust=1.
      upbend=0.
      egr=0.
      ggr=0.
      sgr=0.
      epr=0.
      gpr=0.
      tpr=0.
      egradjust=1.
      ggradjust=1.
      sgradjust=1.
      epradjust=1.
      gpradjust=1.
      tpradjust=1.
      do Zix=0,numZ
        do Nix=0,numN
          if (k0.le.1) then
            if (strength.eq.8) then
              if (ldmodel(Zix,Nix).eq.1) then
                if (flagcol(Zix,Nix)) then
                  wtable(Zix,Nix,1,1)=1.052
                else
                  wtable(Zix,Nix,1,1)=1.076
                endif
              endif
              if (ldmodel(Zix,Nix).eq.2) then
                if (flagcol(Zix,Nix)) then
                  wtable(Zix,Nix,1,1)=0.942
                else
                  wtable(Zix,Nix,1,1)=0.943
                endif
              endif
              if (ldmodel(Zix,Nix).eq.3) then
                if (flagcol(Zix,Nix)) then
                  wtable(Zix,Nix,1,1)=0.905
                else
                  wtable(Zix,Nix,1,1)=0.911
                endif
              endif
              if (ldmodel(Zix,Nix).eq.4) wtable(Zix,Nix,1,1)=0.918
              if (ldmodel(Zix,Nix).eq.5) wtable(Zix,Nix,1,1)=1.017
              if (ldmodel(Zix,Nix).eq.6) wtable(Zix,Nix,1,1)=0.942
            endif
            if (strength.eq.9) then
              if (ldmodel(Zix,Nix).eq.1) then
                if (flagcol(Zix,Nix)) then
                  wtable(Zix,Nix,1,1)=1.048
                else
                  wtable(Zix,Nix,1,1)=1.081
                endif
              endif
              if (ldmodel(Zix,Nix).eq.2) then
                if (flagcol(Zix,Nix)) then
                  wtable(Zix,Nix,1,1)=0.911
                else
                  wtable(Zix,Nix,1,1)=0.934
                endif
              endif
              if (ldmodel(Zix,Nix).eq.3) then
                if (flagcol(Zix,Nix)) then
                  wtable(Zix,Nix,1,1)=0.884
                else
                  wtable(Zix,Nix,1,1)=0.919
                endif
              endif
              if (ldmodel(Zix,Nix).eq.4) wtable(Zix,Nix,1,1)=0.921
              if (ldmodel(Zix,Nix).eq.5) wtable(Zix,Nix,1,1)=1.021
              if (ldmodel(Zix,Nix).eq.6) wtable(Zix,Nix,1,1)=0.936
            endif
          endif
          if (strengthM1.eq.8.or.strengthM1.eq.10) then
            if (Ainit.ge.105) then
              upbend(Zix,Nix,0,1,1)=1.e-8
              upbend(Zix,Nix,0,1,3)=0.
            else
              upbend(Zix,Nix,0,1,1)=3.e-8
              upbend(Zix,Nix,0,1,3)=4.
            endif
          endif
          if (strengthM1.eq.3) then
            upbend(Zix,Nix,0,1,1)=3.5e-8
            upbend(Zix,Nix,0,1,3)=6.
          endif
          upbend(Zix,Nix,0,1,2)=0.8
          if (strength.eq.8) then
            upbend(Zix,Nix,1,1,1)=1.e-10
            upbend(Zix,Nix,1,1,2)=3.
          else
            upbend(Zix,Nix,1,1,1)=0.
            upbend(Zix,Nix,1,1,2)=0.
          endif
        enddo
      enddo
      fiso=-1.
      fisom=-1.
      fisominit=1.
      aadjust=1.
      gnadjust=1.
      gpadjust=1.
      gadjust=1.
      gamgamadjust=1.
      axtype=1
      fbarrier=0.
      fwidth=0.
      bdamp=0.01
      fbaradjust=1.
      fwidthadjust=1.
      bdampadjust=1.
      Rtransmom=1.
      Rclass2mom=1.
      widthc2=0.2
      Ufermi=45.
      cfermi=5.
      betafiscor=1.
      do Zix=0,numZ
        do Nix=0,numN
          Z=Zinit-Zix
          N=Ninit-Nix
          oddZ=mod(Z,2)
          oddN=mod(N,2)
          Vf0=0.86
          if (oddZ.eq.0.and.oddN.eq.0) Vf0=0.83
          if (oddZ.eq.1.and.oddN.eq.0) Vf0=0.91
          if (oddZ.eq.0.and.oddN.eq.1) Vf0=0.86
          if (oddZ.eq.1.and.oddN.eq.1) Vf0=0.85
          A=Z+N
          Aact=max(min(A,255),225)
          rfiscor=0.005
          vfiscor(Zix,Nix)=Vf0-rfiscor*(Aact-240)
          if (Ninit-Nix.gt.144.or.fismodel.eq.5) axtype(Zix,Nix,1)=3
          if (fismodel.lt.5) axtype(Zix,Nix,2)=2
          Rtransmom(Zix,Nix,1)=0.6
          Ufermi(Zix,Nix,0)=30.
          cfermi(Zix,Nix,0)=5.
          do ibar=1,numbar
            Ufermi(Zix,Nix,ibar)=45.
            cfermi(Zix,Nix,ibar)=5.
          enddo
        enddo
      enddo
      betafiscoradjust=1.
      vfiscoradjust=1.
      fismodelx=fismodel
      hbtransfile='                                       '
      clas2file='                                         '
      densfile='                                          '
      Exlfile='                                '
      levelfile='                                               '
      deformfile='                                              '
      optmodfileN='                                            '
      optmodfileP='                                            '
      radialfile='                                              '
      ompenergyfile='                                            '
      yieldfile='                                            '
      radialmodel=2
      breakupmodel=1
      massdir='                                                        '
      massnucleus=0.
      massexcess=0.
      beta2=0.
      do Zix=0,numZ+4
        do Nix=0,numN+4
          beta2(Zix,Nix,1)=0.6
          beta2(Zix,Nix,2)=0.8
          beta2(Zix,Nix,3)=1.
        enddo
      enddo
      msdbins=6
      Emsdmin=0.
      Cstrip=1.
      Cknock=1.
      Cbreak=1.
c
c Extra setting for problems with Kalbach model for (alpha,p) reactions
c
      if (k0.eq.6) Cbreak=0.
      branchlevel=0
      branchratio=0.
      nbranch=0
      v1adjust=1.
      v2adjust=1.
      v3adjust=1.
      v4adjust=1.
      rvadjust=1.
      avadjust=1.
      w1adjust=1.
      w2adjust=1.
      w3adjust=1.
      w4adjust=1.
      rwadjust=-1.
      awadjust=-1.
      rvdadjust=-1.
      avdadjust=-1.
      d1adjust=1.
      d2adjust=1.
      d3adjust=1.
      rwdadjust=1.
      awdadjust=1.
      vso1adjust=1.
      vso2adjust=1.
      rvsoadjust=1.
      avsoadjust=1.
      wso1adjust=1.
      wso2adjust=1.
      rwsoadjust=-1.
      awsoadjust=-1.
      rcadjust=1.
      Vinfadjust=1.
      ompadjustN=0
      ompadjustE1=0.
      ompadjustE2=0.
      ompadjustD=1.
      ompadjusts=1.
      ompadjustF=.false.
      ompadjustp=.false.
      Ejoin=200.
      riplomp=0
      preeqadjust=.false.
      Nadjust=0
      adjustkey='                                                 '
      adjustfile='                                                '
      nenadjust=0
      adjustix=0
      adjustpar=0.
      Eadjust=0.
      Dadjust=0.
      flagrescue=.false.
      rescuefile='                                          '
      grescue=1.
      jlmmode=0
      lvadjust=1.
      lwadjust=1.
      lv1adjust=1.
      lw1adjust=1.
      lvsoadjust=1.
      lwsoadjust=1.
      aradialcor=1.
      adepthcor=1.
      if (flaglabddx) then
        nanglerec=numangrec
      else
        nanglerec=1
      endif
      RprimeU=0.
      nTmax=numT
      astroT9=0.
      astroE=0.
      Ebeam=-1.
      Eback=-1.
      Ibeam=1.
      radiounit='gbq'
      yieldunit='num'
      Area=1.
      Tirrad=0
      unitTirrad=' '
      Tcool=0
      unitTcool=' '
      Tirrad(1)=1
      unitTirrad(1)='d'
      Tcool(1)=1
      unitTcool(1)='d'
      rhotarget=-1.
      flagriplomp=.false.
      if (.not.flagsoukhoinp.and.Atarget.gt.fislim) then
        if ((Ztarget.ge.90.and.Ztarget.le.97.and.Atarget.ge.228.and.
     +    Atarget.le.249).or.flagriplrisk) then
          flagriplrisk=.true.
          flagriplomp=.true.
          flagsoukho=.false.
          riplomp(1)=2408
        endif
      endif
c
c **************** Read fifth set of input variables *******************
c
c nlines     : number of input lines
c getkeywords: subroutine to retrieve keywords and values from input
c              line
c inline     : input line
c word       : words on input line
c key        : keyword
c value      : value or string
c ch         : character
c
c The keyword is identified and the corresponding values are read.
c Erroneous input is immediately checked. The keywords and number of
c values on each line are retrieved from the input.
c
      do i=1,nlines
        line = inline(i)
        call getkeywords(line,word)
        key=word(1)
        value=word(2)
        ch=word(2)(1:1)
        Zix=0
        Nix=0
        type=0
        lval=0
        ibar=0
        igr=1
c
c Test for keywords
c
c Here, the various model parameters can be set to overrule the default
c values. Most default values will be computed later on, since they
c require more computation (e.g. level density parameters).
c
c Each keyword is characterized by a certain order of parameter
c and value input. They are distinguished by different classes.
c
c Classes:
c
c 1: keyword Z A real-value [optional: local adjustment]
c 2: keyword Z A integer-value
c 3: keyword Z A real-value barrier [optional: local adjustment]
c 4: keyword Z A integer-value barrier [optional: local adjustment]
c 5: keyword Z A real-value rad-type l-val [optional: local adjustment]
c 6: keyword particle-type real-value [optional: local adjustment]
c 7: keyword particle-type integer-value
c 8: keyword value [optional: local adjustment]
c 9: keyword Z filename
c
c class     : input class
c getvalues : subroutine to assign values to keywords
c iz        : charge number
c ia        : mass number
c val       : real value
c ival      : integer value
c cval      : character value
c flagassign: flag to assign value or not
c ilev0     : counter for level
c ilev1     : counter for level
c iword     : word counter
c nbr       : number of branches
c
        if (key.eq.'elow') then
          read(value,*,end=1000,err=1000) eninclow
          cycle
        endif
        if (key.eq.'massnucleus') then
          class=1
          call getvalues(class,word,Zix,Nix,type,
     +      ibar,irad,lval,igr,val,ival,cval,flagassign)
          if (flagassign) massnucleus(Zix,Nix)=val
          cycle
        endif
        if (key.eq.'massexcess') then
          class=1
          call getvalues(class,word,Zix,Nix,type,
     +      ibar,irad,lval,igr,val,ival,cval,flagassign)
          if (flagassign) massexcess(Zix,Nix)=val
          cycle
        endif
        if (key.eq.'massdir') then
          read(value,*,end=1000,err=1000) massdir
          cycle
        endif
        if (key.eq.'a') then
          class=1
          call getvalues(class,word,Zix,Nix,type,
     +      ibar,irad,lval,igr,val,ival,cval,flagassign)
          if (flagassign) alev(Zix,Nix)=val
          cycle
        endif
        if (key.eq.'alimit') then
          class=1
          call getvalues(class,word,Zix,Nix,type,
     +      ibar,irad,lval,igr,val,ival,cval,flagassign)
          if (flagassign) alimit(Zix,Nix)=val
          cycle
        endif
        if (key.eq.'gammald') then
          class=1
          call getvalues(class,word,Zix,Nix,type,
     +      ibar,irad,lval,igr,val,ival,cval,flagassign)
          if (flagassign) gammald(Zix,Nix)=val
          cycle
        endif
        if (key.eq.'pair') then
          class=1
          call getvalues(class,word,Zix,Nix,type,
     +      ibar,irad,lval,igr,val,ival,cval,flagassign)
          if (flagassign) pair(Zix,Nix)=val
          cycle
        endif
        if (key.eq.'pshift') then
          class=3
          call getvalues(class,word,Zix,Nix,type,
     +      ibar,irad,lval,igr,val,ival,cval,flagassign)
          if (flagassign) pshift(Zix,Nix,ibar)=val
          cycle
        endif
        if (key.eq.'pshiftadjust') then
          class=3
          call getvalues(class,word,Zix,Nix,type,
     +      ibar,irad,lval,igr,val,ival,cval,flagassign)
          if (flagassign) pshiftadjust(Zix,Nix,ibar)=val
          cycle
        endif
        if (key.eq.'deltaw') then
          class=3
          call getvalues(class,word,Zix,Nix,type,
     +      ibar,irad,lval,igr,val,ival,cval,flagassign)
          if (flagassign) deltaW(Zix,Nix,ibar)=val
          cycle
        endif
        if (key.eq.'exmatch') then
          class=3
          call getvalues(class,word,Zix,Nix,type,
     +      ibar,irad,lval,igr,val,ival,cval,flagassign)
          if (flagassign) Exmatch(Zix,Nix,ibar)=val
          cycle
        endif
        if (key.eq.'t') then
          class=3
          call getvalues(class,word,Zix,Nix,type,
     +      ibar,irad,lval,igr,val,ival,cval,flagassign)
          if (flagassign) T(Zix,Nix,ibar)=val
          cycle
        endif
        if (key.eq.'e0') then
          class=3
          call getvalues(class,word,Zix,Nix,type,
     +      ibar,irad,lval,igr,val,ival,cval,flagassign)
          if (flagassign) E0(Zix,Nix,ibar)=val
          cycle
        endif
        if (key.eq.'nlow') then
          class=4
          call getvalues(class,word,Zix,Nix,type,
     +      ibar,irad,lval,igr,val,ival,cval,flagassign)
          if (flagassign) Nlow(Zix,Nix,ibar)=ival
          cycle
        endif
        if (key.eq.'ntop') then
          class=4
          call getvalues(class,word,Zix,Nix,type,
     +      ibar,irad,lval,igr,val,ival,cval,flagassign)
          if (flagassign) Ntop(Zix,Nix,ibar)=ival
          cycle
        endif
        if (key.eq.'s2adjust') then
          class=3
          call getvalues(class,word,Zix,Nix,type,
     +      ibar,irad,lval,igr,val,ival,cval,flagassign)
          if (flagassign) s2adjust(Zix,Nix,ibar)=val
          cycle
        endif
        if (key.eq.'krotconstant') then
          class=3
          call getvalues(class,word,Zix,Nix,type,
     +      ibar,irad,lval,igr,val,ival,cval,flagassign)
          if (flagassign) then
            Krotconstant(Zix,Nix,ibar)=val
            ldadjust(Zix,Nix)=.true.
          endif
          cycle
        endif
        if (key.eq.'ufermi') then
          class=3
          call getvalues(class,word,Zix,Nix,type,
     +      ibar,irad,lval,igr,val,ival,cval,flagassign)
          if (flagassign) then
            Ufermi(Zix,Nix,ibar)=val
            ldadjust(Zix,Nix)=.true.
          endif
          cycle
        endif
        if (key.eq.'cfermi') then
          class=3
          call getvalues(class,word,Zix,Nix,type,
     +      ibar,irad,lval,igr,val,ival,cval,flagassign)
          if (flagassign) then
            cfermi(Zix,Nix,ibar)=val
            ldadjust(Zix,Nix)=.true.
          endif
          cycle
        endif
        if (key.eq.'ctable') then
          class=3
          call getvalues(class,word,Zix,Nix,type,
     +      ibar,irad,lval,igr,val,ival,cval,flagassign)
          if (flagassign) then
            ctable(Zix,Nix,ibar)=val
            ldadjust(Zix,Nix)=.true.
          endif
          cycle
        endif
        if (key.eq.'ctableadjust') then
          class=3
          call getvalues(class,word,Zix,Nix,type,
     +      ibar,irad,lval,igr,val,ival,cval,flagassign)
          if (flagassign) then
            ctableadjust(Zix,Nix,ibar)=val
            ldadjust(Zix,Nix)=.true.
          endif
          cycle
        endif
        if (key.eq.'ptable') then
          class=3
          call getvalues(class,word,Zix,Nix,type,
     +      ibar,irad,lval,igr,val,ival,cval,flagassign)
          if (flagassign) ptable(Zix,Nix,ibar)=val
          cycle
        endif
        if (key.eq.'ptableadjust') then
          class=3
          call getvalues(class,word,Zix,Nix,type,
     +      ibar,irad,lval,igr,val,ival,cval,flagassign)
          if (flagassign) ptableadjust(Zix,Nix,ibar)=val
          cycle
        endif
        if (key.eq.'g') then
          class=1
          call getvalues(class,word,Zix,Nix,type,
     +      ibar,irad,lval,igr,val,ival,cval,flagassign)
          if (flagassign) g(Zix,Nix)=val
          cycle
        endif
        if (key.eq.'gp') then
          class=1
          call getvalues(class,word,Zix,Nix,type,
     +      ibar,irad,lval,igr,val,ival,cval,flagassign)
          if (flagassign) gp(Zix,Nix)=val
          cycle
        endif
        if (key.eq.'gn') then
          class=1
          call getvalues(class,word,Zix,Nix,type,
     +      ibar,irad,lval,igr,val,ival,cval,flagassign)
          if (flagassign) gn(Zix,Nix)=val
          cycle
        endif
        if (key.eq.'risomer') then
          class=1
          call getvalues(class,word,Zix,Nix,type,
     +      ibar,irad,lval,igr,val,ival,cval,flagassign)
          if (flagassign) Risomer(Zix,Nix)=val
          cycle
        endif
        if (key.eq.'egr') then
          class=5
          call getvalues(class,word,Zix,Nix,type,
     +      ibar,irad,lval,igr,val,ival,cval,flagassign)
          if (flagassign) then
            egr(Zix,Nix,irad,lval,igr)=val
            gamadjust(Zix,Nix)=.true.
          endif
          cycle
        endif
        if (key.eq.'ggr') then
          class=5
          call getvalues(class,word,Zix,Nix,type,
     +      ibar,irad,lval,igr,val,ival,cval,flagassign)
          if (flagassign) then
            ggr(Zix,Nix,irad,lval,igr)=val
            gamadjust(Zix,Nix)=.true.
          endif
          cycle
        endif
        if (key.eq.'sgr') then
          class=5
          call getvalues(class,word,Zix,Nix,type,
     +      ibar,irad,lval,igr,val,ival,cval,flagassign)
          if (flagassign) then
            sgr(Zix,Nix,irad,lval,igr)=val
            gamadjust(Zix,Nix)=.true.
          endif
          cycle
        endif
        if (key.eq.'epr') then
          class=5
          call getvalues(class,word,Zix,Nix,type,
     +      ibar,irad,lval,igr,val,ival,cval,flagassign)
          if (flagassign) then
            epr(Zix,Nix,irad,lval,igr)=val
            gamadjust(Zix,Nix)=.true.
          endif
          cycle
        endif
        if (key.eq.'gpr') then
          class=5
          call getvalues(class,word,Zix,Nix,type,
     +      ibar,irad,lval,igr,val,ival,cval,flagassign)
          if (flagassign) then
            gpr(Zix,Nix,irad,lval,igr)=val
            gamadjust(Zix,Nix)=.true.
          endif
          cycle
        endif
        if (key.eq.'spr') then
          class=5
          call getvalues(class,word,Zix,Nix,type,
     +      ibar,irad,lval,igr,val,ival,cval,flagassign)
          if (flagassign) then
            tpr(Zix,Nix,irad,lval,igr)=val
            gamadjust(Zix,Nix)=.true.
          endif
          cycle
        endif
        if (key.eq.'egradjust') then
          class=5
          call getvalues(class,word,Zix,Nix,type,
     +      ibar,irad,lval,igr,val,ival,cval,flagassign)
          if (flagassign) then
            egradjust(Zix,Nix,irad,lval,igr)=val
            gamadjust(Zix,Nix)=.true.
          endif
          cycle
        endif
        if (key.eq.'ggradjust') then
          class=5
          call getvalues(class,word,Zix,Nix,type,
     +      ibar,irad,lval,igr,val,ival,cval,flagassign)
          if (flagassign) then
            ggradjust(Zix,Nix,irad,lval,igr)=val
            gamadjust(Zix,Nix)=.true.
          endif
          cycle
        endif
        if (key.eq.'sgradjust') then
          class=5
          call getvalues(class,word,Zix,Nix,type,
     +      ibar,irad,lval,igr,val,ival,cval,flagassign)
          if (flagassign) then
            sgradjust(Zix,Nix,irad,lval,igr)=val
            gamadjust(Zix,Nix)=.true.
          endif
          cycle
        endif
        if (key.eq.'epradjust') then
          class=5
          call getvalues(class,word,Zix,Nix,type,
     +      ibar,irad,lval,igr,val,ival,cval,flagassign)
          if (flagassign) then
            epradjust(Zix,Nix,irad,lval,igr)=val
            gamadjust(Zix,Nix)=.true.
          endif
          cycle
        endif
        if (key.eq.'gpradjust') then
          class=5
          call getvalues(class,word,Zix,Nix,type,
     +      ibar,irad,lval,igr,val,ival,cval,flagassign)
          if (flagassign) then
            gpradjust(Zix,Nix,irad,lval,igr)=val
            gamadjust(Zix,Nix)=.true.
          endif
          cycle
        endif
        if (key.eq.'spradjust') then
          class=5
          call getvalues(class,word,Zix,Nix,type,
     +      ibar,irad,lval,igr,val,ival,cval,flagassign)
          if (flagassign) then
            tpradjust(Zix,Nix,irad,lval,igr)=val
            gamadjust(Zix,Nix)=.true.
          endif
          cycle
        endif
        if (key.eq.'upbendc') then
          class=5
          call getvalues(class,word,Zix,Nix,type,
     +      ibar,irad,lval,igr,val,ival,cval,flagassign)
          if (flagassign) upbend(Zix,Nix,irad,lval,1)=val
          cycle
        endif
        if (key.eq.'upbende') then
          class=5
          call getvalues(class,word,Zix,Nix,type,
     +      ibar,irad,lval,igr,val,ival,cval,flagassign)
          if (flagassign) upbend(Zix,Nix,irad,lval,2)=val
          cycle
        endif
        if (key.eq.'upbendf') then
          class=5
          call getvalues(class,word,Zix,Nix,type,
     +      ibar,irad,lval,igr,val,ival,cval,flagassign)
          if (flagassign) upbend(Zix,Nix,irad,lval,3)=val
          cycle
        endif
        if (key.eq.'fiso') then
          class=6
          call getvalues(class,word,Zix,Nix,type,
     +      ibar,irad,lval,igr,val,ival,cval,flagassign)
          if (flagassign) fiso(type)=val
          cycle
        endif
        if (key.eq.'fisom') then
          class=6
          call getvalues(class,word,Zix,Nix,type,
     +      ibar,irad,lval,igr,val,ival,cval,flagassign)
          if (flagassign) then
            fisom(type)=val
            fisominit(type)=val
          endif
          cycle
        endif
        if (key.eq.'gamgam') then
          class=1
          call getvalues(class,word,Zix,Nix,type,
     +      ibar,irad,lval,igr,val,ival,cval,flagassign)
          if (flagassign) gamgam(Zix,Nix)=val
          cycle
        endif
        if (key.eq.'d0') then
          class=1
          call getvalues(class,word,Zix,Nix,type,
     +      ibar,irad,lval,igr,val,ival,cval,flagassign)
          if (flagassign) D0(Zix,Nix)=val*1000.
          cycle
        endif
        if (key.eq.'aadjust') then
          class=1
          call getvalues(class,word,Zix,Nix,type,
     +      ibar,irad,lval,igr,val,ival,cval,flagassign)
          if (flagassign) aadjust(Zix,Nix)=val
          cycle
        endif
        if (key.eq.'tadjust') then
          class=3
          call getvalues(class,word,Zix,Nix,type,
     +      ibar,irad,lval,igr,val,ival,cval,flagassign)
          if (flagassign) Tadjust(Zix,Nix,ibar)=val
          cycle
        endif
        if (key.eq.'e0adjust') then
          class=3
          call getvalues(class,word,Zix,Nix,type,
     +      ibar,irad,lval,igr,val,ival,cval,flagassign)
          if (flagassign) E0adjust(Zix,Nix,ibar)=val
          cycle
        endif
        if (key.eq.'exmatchadjust') then
          class=3
          call getvalues(class,word,Zix,Nix,type,
     +      ibar,irad,lval,igr,val,ival,cval,flagassign)
          if (flagassign) Exmatchadjust(Zix,Nix,ibar)=val
          cycle
        endif
        if (key.eq.'gnadjust') then
          class=1
          call getvalues(class,word,Zix,Nix,type,
     +      ibar,irad,lval,igr,val,ival,cval,flagassign)
          if (flagassign) gnadjust(Zix,Nix)=val
          cycle
        endif
        if (key.eq.'gpadjust') then
          class=1
          call getvalues(class,word,Zix,Nix,type,
     +      ibar,irad,lval,igr,val,ival,cval,flagassign)
          if (flagassign) gpadjust(Zix,Nix)=val
          cycle
        endif
        if (key.eq.'gadjust') then
          class=1
          call getvalues(class,word,Zix,Nix,type,
     +      ibar,irad,lval,igr,val,ival,cval,flagassign)
          if (flagassign) gadjust(Zix,Nix)=val
          cycle
        endif
        if (key.eq.'gamgamadjust') then
          class=1
          call getvalues(class,word,Zix,Nix,type,
     +      ibar,irad,lval,igr,val,ival,cval,flagassign)
          if (flagassign) gamgamadjust(Zix,Nix)=val
          cycle
        endif
        if (key.eq.'etable') then
          class=5
          call getvalues(class,word,Zix,Nix,type,
     +      ibar,irad,lval,igr,val,ival,cval,flagassign)
          if (flagassign) then
            etable(Zix,Nix,irad,lval)=val
            gamadjust(Zix,Nix)=.true.
          endif
          cycle
        endif
        if (key.eq.'etableadjust') then
          class=5
          call getvalues(class,word,Zix,Nix,type,
     +      ibar,irad,lval,igr,val,ival,cval,flagassign)
          if (flagassign) then
            etableadjust(Zix,Nix,irad,lval)=val
            gamadjust(Zix,Nix)=.true.
          endif
          cycle
        endif
        if (key.eq.'ftable') then
          class=5
          call getvalues(class,word,Zix,Nix,type,
     +      ibar,irad,lval,igr,val,ival,cval,flagassign)
          if (flagassign) then
            ftable(Zix,Nix,irad,lval)=val
            gamadjust(Zix,Nix)=.true.
          endif
          cycle
        endif
        if (key.eq.'ftableadjust') then
          class=5
          call getvalues(class,word,Zix,Nix,type,
     +      ibar,irad,lval,igr,val,ival,cval,flagassign)
          if (flagassign) then
            ftableadjust(Zix,Nix,irad,lval)=val
            gamadjust(Zix,Nix)=.true.
          endif
          cycle
        endif
        if (key.eq.'wtable') then
          class=5
          call getvalues(class,word,Zix,Nix,type,
     +      ibar,irad,lval,igr,val,ival,cval,flagassign)
          if (flagassign) then
            wtable(Zix,Nix,irad,lval)=val
            gamadjust(Zix,Nix)=.true.
            wtadjust=.true.
          endif
          cycle
        endif
        if (key.eq.'wtableadjust') then
          class=5
          call getvalues(class,word,Zix,Nix,type,
     +      ibar,irad,lval,igr,val,ival,cval,flagassign)
          if (flagassign) then
            wtableadjust(Zix,Nix,irad,lval)=val
            gamadjust(Zix,Nix)=.true.
            wtadjust=.true.
          endif
          cycle
        endif
        if (key.eq.'fisbar') then
          ibar=1
          class=3
          call getvalues(class,word,Zix,Nix,type,
     +      ibar,irad,lval,igr,val,ival,cval,flagassign)
          if (flagassign) then
            fbarrier(Zix,Nix,ibar)=val
            fisadjust(Zix,Nix)=.true.
          endif
          cycle
        endif
        if (key.eq.'fishw') then
          ibar=1
          class=3
          call getvalues(class,word,Zix,Nix,type,
     +      ibar,irad,lval,igr,val,ival,cval,flagassign)
          if (flagassign) then
            fwidth(Zix,Nix,ibar)=val
            fisadjust(Zix,Nix)=.true.
          endif
          cycle
        endif
        if (key.eq.'bdamp') then
          ibar=1
          class=3
          call getvalues(class,word,Zix,Nix,type,
     +      ibar,irad,lval,igr,val,ival,cval,flagassign)
          if (flagassign) then
            bdamp(Zix,Nix,ibar)=val
            fisadjust(Zix,Nix)=.true.
          endif
          cycle
        endif
        if (key.eq.'fisbaradjust') then
          ibar=1
          class=3
          call getvalues(class,word,Zix,Nix,type,
     +      ibar,irad,lval,igr,val,ival,cval,flagassign)
          if (flagassign) then
            fbaradjust(Zix,Nix,ibar)=val
            fisadjust(Zix,Nix)=.true.
          endif
          cycle
        endif
        if (key.eq.'fishwadjust') then
          ibar=1
          class=3
          call getvalues(class,word,Zix,Nix,type,
     +      ibar,irad,lval,igr,val,ival,cval,flagassign)
          if (flagassign) then
            fwidthadjust(Zix,Nix,ibar)=val
            fisadjust(Zix,Nix)=.true.
          endif
          cycle
        endif
        if (key.eq.'bdampadjust') then
          ibar=1
          class=3
          call getvalues(class,word,Zix,Nix,type,
     +      ibar,irad,lval,igr,val,ival,cval,flagassign)
          if (flagassign) then
            bdampadjust(Zix,Nix,ibar)=val
            fisadjust(Zix,Nix)=.true.
          endif
          cycle
        endif
        if (key.eq.'rtransmom') then
          ibar=1
          class=3
          call getvalues(class,word,Zix,Nix,type,
     +      ibar,irad,lval,igr,val,ival,cval,flagassign)
          if (flagassign) Rtransmom(Zix,Nix,ibar)=val
          cycle
        endif
        if (key.eq.'rclass2mom') then
          ibar=1
          class=3
          call getvalues(class,word,Zix,Nix,type,
     +      ibar,irad,lval,igr,val,ival,cval,flagassign)
          if (flagassign) Rclass2mom(Zix,Nix,ibar)=val
          cycle
        endif
        if (key.eq.'class2width') then
          ibar=1
          class=3
          call getvalues(class,word,Zix,Nix,type,
     +      ibar,irad,lval,igr,val,ival,cval,flagassign)
          if (flagassign) widthc2(Zix,Nix,ibar)=val
          cycle
        endif
        if (key.eq.'levelfile') then
          class=10
          call getvalues(class,word,Zix,Nix,type,
     +      ibar,irad,lval,igr,val,ival,cval,flagassign)
          if (flagassign) levelfile(Zix)=cval
          cycle
        endif
        if (key.eq.'deformfile') then
          class=10
          call getvalues(class,word,Zix,Nix,type,
     +      ibar,irad,lval,igr,val,ival,cval,flagassign)
          if (flagassign) deformfile(Zix)=cval
          cycle
        endif
        if (key.eq.'e1file') then
          class=11
          call getvalues(class,word,Zix,Nix,type,
     +      ibar,irad,lval,igr,val,ival,cval,flagassign)
          if (flagassign) Exlfile(Zix,Nix,1,1)=cval
          cycle
        endif
        if (key.eq.'m1file') then
          class=11
          call getvalues(class,word,Zix,Nix,type,
     +      ibar,irad,lval,igr,val,ival,cval,flagassign)
          if (flagassign) Exlfile(Zix,Nix,0,1)=cval
          cycle
        endif
        if (key.eq.'densfile') then
          class=11
          call getvalues(class,word,Zix,Nix,type,
     +      ibar,irad,lval,igr,val,ival,cval,flagassign)
          if (flagassign) densfile(Zix,Nix)=cval
          cycle
        endif
        if (key.eq.'hbtransfile') then
          class=11
          call getvalues(class,word,Zix,Nix,type,
     +      ibar,irad,lval,igr,val,ival,cval,flagassign)
          if (flagassign) hbtransfile(Zix,Nix)=cval
          cycle
        endif
        if (key.eq.'class2file') then
          class=11
          call getvalues(class,word,Zix,Nix,type,
     +      ibar,irad,lval,igr,val,ival,cval,flagassign)
          if (flagassign) clas2file(Zix,Nix)=cval
          cycle
        endif
        if (key.eq.'optmodfilen') then
          class=10
          call getvalues(class,word,Zix,Nix,type,
     +      ibar,irad,lval,igr,val,ival,cval,flagassign)
          if (flagassign) optmodfilen(Zix)=cval
          cycle
        endif
        if (key.eq.'optmodfilep') then
          class=10
          call getvalues(class,word,Zix,Nix,type,
     +      ibar,irad,lval,igr,val,ival,cval,flagassign)
          if (flagassign) optmodfilep(Zix)=cval
          cycle
        endif
        if (key.eq.'radialfile') then
          class=10
          call getvalues(class,word,Zix,Nix,type,
     +      ibar,irad,lval,igr,val,ival,cval,flagassign)
          if (flagassign) radialfile(Zix)=cval
          cycle
        endif
        if (key.eq.'ompenergyfile') then
          ompenergyfile=value
          cycle
        endif
        if (key.eq.'yieldfile') then
          yieldfile=value
          cycle
        endif
        if (key.eq.'beta2') then
          class=3
          call getvalues(class,word,Zix,Nix,type,
     +      ibar,irad,lval,igr,val,ival,cval,flagassign)
          if (flagassign) beta2(Zix,Nix,ibar)=val
          cycle
        endif
        if (key.eq.'axtype') then
          ibar=1
          class=4
          call getvalues(class,word,Zix,Nix,type,
     +      ibar,irad,lval,igr,val,ival,cval,flagassign)
          if (flagassign) axtype(Zix,Nix,ibar)=ival
          cycle
        endif
        if (key.eq.'vfiscor') then
          class=1
          call getvalues(class,word,Zix,Nix,type,
     +      ibar,irad,lval,igr,val,ival,cval,flagassign)
          if (flagassign) vfiscor(Zix,Nix)=val
          cycle
        endif
        if (key.eq.'vfiscoradjust') then
          class=1
          call getvalues(class,word,Zix,Nix,type,
     +      ibar,irad,lval,igr,val,ival,cval,flagassign)
          if (flagassign) vfiscoradjust(Zix,Nix)=val
          cycle
        endif
        if (key.eq.'cbarrier') then
          class=9
          call getvalues(class,word,Zix,Nix,type,
     +      ibar,irad,lval,igr,val,ival,cval,flagassign)
          if (flagassign) Cbarrier=val
          cycle
        endif
        if (key.eq.'betafiscor') then
          class=1
          call getvalues(class,word,Zix,Nix,type,
     +      ibar,irad,lval,igr,val,ival,cval,flagassign)
          if (flagassign) betafiscor(Zix,Nix)=val
          cycle
        endif
        if (key.eq.'betafiscoradjust') then
          class=1
          call getvalues(class,word,Zix,Nix,type,
     +      ibar,irad,lval,igr,val,ival,cval,flagassign)
          if (flagassign) betafiscoradjust(Zix,Nix)=val
          cycle
        endif
        if (key.eq.'branch') then
          read(word(2),*,end=1000,err=1000) iz
          read(word(3),*,end=1000,err=1000) ia
          read(word(4),*,end=1000,err=1000) ilev0
          read(word(5),*,end=1000,err=1000) nbr
          Zix=Zinit-iz
          Nix=Ninit-ia+iz
          if (Zix.lt.0.or.Zix.gt.numZ.or.Nix.lt.0.or.Nix.gt.numN)
     +      goto 1040
          if (ilev0.lt.0..or.ilev0.gt.numlev) goto 1060
          if (nbr.lt.0..or.nbr.gt.numlev) goto 1060
          iword=5
          sum=0.
          do k=1,nbr
            iword=iword+1
            read(word(iword),*,end=1000,err=1000) ilev1
            if (ilev1.lt.0..or.ilev1.gt.numlev) goto 1060
            iword=iword+1
            read(word(iword),*,end=1000,err=1000) br
            if (br.lt.0.) goto 1070
            branchlevel(Zix,Nix,ilev0,k)=ilev1
            branchratio(Zix,Nix,ilev0,k)=br
            sum=sum+br
          enddo
          if (sum.gt.0.) then
            do k=1,nbr
              branchratio(Zix,Nix,ilev0,k)=branchratio(Zix,Nix,ilev0,k)
     +          /sum
            enddo
          endif
          nbranch(Zix,Nix,ilev0)=nbr
          cycle
        endif
        if (key.eq.'cstrip') then
          class=6
          call getvalues(class,word,Zix,Nix,type,
     +      ibar,irad,lval,igr,val,ival,cval,flagassign)
          if (flagassign) Cstrip(type)=val
          cycle
        endif
        if (key.eq.'cknock') then
          class=6
          call getvalues(class,word,Zix,Nix,type,
     +      ibar,irad,lval,igr,val,ival,cval,flagassign)
          if (flagassign) Cknock(type)=val
          cycle
        endif
        if (key.eq.'cbreak') then
          class=6
          call getvalues(class,word,Zix,Nix,type,
     +      ibar,irad,lval,igr,val,ival,cval,flagassign)
          if (flagassign) Cbreak(type)=val
          cycle
        endif
        if (key.eq.'v1adjust') then
          class=6
          call getvalues(class,word,Zix,Nix,type,
     +      ibar,irad,lval,igr,val,ival,cval,flagassign)
          if (flagassign) then
            v1adjust(type)=val
            ompadjustp(type)=.true.
          endif
          cycle
        endif
        if (key.eq.'v2adjust') then
          class=6
          call getvalues(class,word,Zix,Nix,type,
     +      ibar,irad,lval,igr,val,ival,cval,flagassign)
          if (flagassign) then
            v2adjust(type)=val
            ompadjustp(type)=.true.
          endif
          cycle
        endif
        if (key.eq.'v3adjust') then
          class=6
          call getvalues(class,word,Zix,Nix,type,
     +      ibar,irad,lval,igr,val,ival,cval,flagassign)
          if (flagassign) then
            v3adjust(type)=val
            ompadjustp(type)=.true.
          endif
          cycle
        endif
        if (key.eq.'v4adjust') then
          class=6
          call getvalues(class,word,Zix,Nix,type,
     +      ibar,irad,lval,igr,val,ival,cval,flagassign)
          if (flagassign) then
            v4adjust(type)=val
            ompadjustp(type)=.true.
          endif
          cycle
        endif
        if (key.eq.'rvadjust') then
          class=6
          call getvalues(class,word,Zix,Nix,type,
     +      ibar,irad,lval,igr,val,ival,cval,flagassign)
          if (flagassign) then
            rvadjust(type)=val
            ompadjustp(type)=.true.
          endif
          cycle
        endif
        if (key.eq.'avadjust') then
          class=6
          call getvalues(class,word,Zix,Nix,type,
     +      ibar,irad,lval,igr,val,ival,cval,flagassign)
          if (flagassign) then
            avadjust(type)=val
            ompadjustp(type)=.true.
          endif
          cycle
        endif
        if (key.eq.'w1adjust') then
          class=6
          call getvalues(class,word,Zix,Nix,type,
     +      ibar,irad,lval,igr,val,ival,cval,flagassign)
          if (flagassign) then
            w1adjust(type)=val
            ompadjustp(type)=.true.
          endif
          cycle
        endif
        if (key.eq.'w2adjust') then
          class=6
          call getvalues(class,word,Zix,Nix,type,
     +      ibar,irad,lval,igr,val,ival,cval,flagassign)
          if (flagassign) then
            w2adjust(type)=val
            ompadjustp(type)=.true.
          endif
          cycle
        endif
        if (key.eq.'w3adjust') then
          class=6
          call getvalues(class,word,Zix,Nix,type,
     +      ibar,irad,lval,igr,val,ival,cval,flagassign)
          if (flagassign) then
            w3adjust(type)=val
            ompadjustp(type)=.true.
          endif
          cycle
        endif
        if (key.eq.'w4adjust') then
          class=6
          call getvalues(class,word,Zix,Nix,type,
     +      ibar,irad,lval,igr,val,ival,cval,flagassign)
          if (flagassign) then
            w4adjust(type)=val
            ompadjustp(type)=.true.
          endif
          cycle
        endif
        if (key.eq.'rwadjust') then
          class=6
          call getvalues(class,word,Zix,Nix,type,
     +      ibar,irad,lval,igr,val,ival,cval,flagassign)
          if (flagassign) then
            rwadjust(type)=val
            ompadjustp(type)=.true.
          endif
          cycle
        endif
        if (key.eq.'awadjust') then
          class=6
          call getvalues(class,word,Zix,Nix,type,
     +      ibar,irad,lval,igr,val,ival,cval,flagassign)
          if (flagassign) then
            awadjust(type)=val
            ompadjustp(type)=.true.
          endif
          cycle
        endif
        if (key.eq.'rvdadjust') then
          class=6
          call getvalues(class,word,Zix,Nix,type,
     +      ibar,irad,lval,igr,val,ival,cval,flagassign)
          if (flagassign) then
            rvdadjust(type)=val
            ompadjustp(type)=.true.
          endif
          cycle
        endif
        if (key.eq.'avdadjust') then
          class=6
          call getvalues(class,word,Zix,Nix,type,
     +      ibar,irad,lval,igr,val,ival,cval,flagassign)
          if (flagassign) avdadjust(type)=val
          ompadjustp(type)=.true.
          cycle
        endif
        if (key.eq.'rwdadjust') then
          class=6
          call getvalues(class,word,Zix,Nix,type,
     +      ibar,irad,lval,igr,val,ival,cval,flagassign)
          if (flagassign) then
            rwdadjust(type)=val
            ompadjustp(type)=.true.
          endif
          cycle
        endif
        if (key.eq.'awdadjust') then
          class=6
          call getvalues(class,word,Zix,Nix,type,
     +      ibar,irad,lval,igr,val,ival,cval,flagassign)
          if (flagassign) then
            awdadjust(type)=val
            ompadjustp(type)=.true.
          endif
          cycle
        endif
        if (key.eq.'d1adjust') then
          class=6
          call getvalues(class,word,Zix,Nix,type,
     +      ibar,irad,lval,igr,val,ival,cval,flagassign)
          if (flagassign) then
            d1adjust(type)=val
            ompadjustp(type)=.true.
          endif
          cycle
        endif
        if (key.eq.'d2adjust') then
          class=6
          call getvalues(class,word,Zix,Nix,type,
     +      ibar,irad,lval,igr,val,ival,cval,flagassign)
          if (flagassign) then
            d2adjust(type)=val
            ompadjustp(type)=.true.
          endif
          cycle
        endif
        if (key.eq.'d3adjust') then
          class=6
          call getvalues(class,word,Zix,Nix,type,
     +      ibar,irad,lval,igr,val,ival,cval,flagassign)
          if (flagassign) then
            d3adjust(type)=val
            ompadjustp(type)=.true.
          endif
          cycle
        endif
        if (key.eq.'rvsoadjust') then
          class=6
          call getvalues(class,word,Zix,Nix,type,
     +      ibar,irad,lval,igr,val,ival,cval,flagassign)
          if (flagassign) then
            rvsoadjust(type)=val
            ompadjustp(type)=.true.
          endif
          cycle
        endif
        if (key.eq.'avsoadjust') then
          class=6
          call getvalues(class,word,Zix,Nix,type,
     +      ibar,irad,lval,igr,val,ival,cval,flagassign)
          if (flagassign) then
            avsoadjust(type)=val
            ompadjustp(type)=.true.
          endif
          cycle
        endif
        if (key.eq.'vso1adjust') then
          class=6
          call getvalues(class,word,Zix,Nix,type,
     +      ibar,irad,lval,igr,val,ival,cval,flagassign)
          if (flagassign) then
            vso1adjust(type)=val
            ompadjustp(type)=.true.
          endif
          cycle
        endif
        if (key.eq.'vso2adjust') then
          class=6
          call getvalues(class,word,Zix,Nix,type,
     +      ibar,irad,lval,igr,val,ival,cval,flagassign)
          if (flagassign) then
            vso2adjust(type)=val
            ompadjustp(type)=.true.
          endif
          cycle
        endif
        if (key.eq.'rwsoadjust') then
          class=6
          call getvalues(class,word,Zix,Nix,type,
     +      ibar,irad,lval,igr,val,ival,cval,flagassign)
          if (flagassign) then
            rwsoadjust(type)=val
            ompadjustp(type)=.true.
          endif
          cycle
        endif
        if (key.eq.'awsoadjust') then
          class=6
          call getvalues(class,word,Zix,Nix,type,
     +      ibar,irad,lval,igr,val,ival,cval,flagassign)
          if (flagassign) then
            awsoadjust(type)=val
            ompadjustp(type)=.true.
          endif
          cycle
        endif
        if (key.eq.'wso1adjust') then
          class=6
          call getvalues(class,word,Zix,Nix,type,
     +      ibar,irad,lval,igr,val,ival,cval,flagassign)
          if (flagassign) then
            wso1adjust(type)=val
            ompadjustp(type)=.true.
          endif
          cycle
        endif
        if (key.eq.'wso2adjust') then
          class=6
          call getvalues(class,word,Zix,Nix,type,
     +      ibar,irad,lval,igr,val,ival,cval,flagassign)
          if (flagassign) then
            wso2adjust(type)=val
            ompadjustp(type)=.true.
          endif
          cycle
        endif
        if (key.eq.'rcadjust') then
          class=6
          call getvalues(class,word,Zix,Nix,type,
     +      ibar,irad,lval,igr,val,ival,cval,flagassign)
          if (flagassign) then
            rcadjust(type)=val
            ompadjustp(type)=.true.
          endif
          cycle
        endif
        if (key.eq.'ejoin') then
          class=6
          call getvalues(class,word,Zix,Nix,type,
     +      ibar,irad,lval,igr,val,ival,cval,flagassign)
          if (flagassign) Ejoin(type)=val
          cycle
        endif
        if (key.eq.'vinfadjust') then
          class=6
          call getvalues(class,word,Zix,Nix,type,
     +      ibar,irad,lval,igr,val,ival,cval,flagassign)
          if (flagassign) Vinfadjust(type)=val
          cycle
        endif
        if (key.eq.'tjadjust') then
          class=12
          call getvalues(class,word,Zix,Nix,type,
     +      ibar,irad,lval,igr,val,ival,cval,flagassign)
          if (flagassign) then
            TJadjust(Zix,Nix,type)=val
            adjustTJ(Zix,Nix,type)=.true.
          endif
          cycle
        endif
        if (key.eq.'rvadjustf'.or.key.eq.'avadjustf'.or.
     +    key.eq.'rwadjustf'.or.key.eq.'awadjustf'.or.
     +    key.eq.'rvdadjustf'.or.key.eq.'avdadjustf'.or.
     +    key.eq.'rwdadjustf'.or.key.eq.'awdadjustf'.or.
     +    key.eq.'rvsoadjustf'.or.key.eq.'avsoadjustf'.or.
     +    key.eq.'rwsoadjustf'.or.key.eq.'awsoadjustf') then
          if (key.eq.'rvadjustf') omptype=1
          if (key.eq.'avadjustf') omptype=2
          if (key.eq.'rwadjustf') omptype=3
          if (key.eq.'awadjustf') omptype=4
          if (key.eq.'rvdadjustf') omptype=5
          if (key.eq.'avdadjustf') omptype=6
          if (key.eq.'rwdadjustf') omptype=7
          if (key.eq.'awdadjustf') omptype=8
          if (key.eq.'rvsoadjustf') omptype=9
          if (key.eq.'avsoadjustf') omptype=10
          if (key.eq.'rwsoadjustf') omptype=11
          if (key.eq.'awsoadjustf') omptype=12
          if (key.eq.'rcadjustf') omptype=13
          do type=1,6
            if (ch.eq.parsym(type)) then
              type2=type
              goto 170
            endif
          enddo
          goto 1000
  170     ompadjustF(type2)=.true.
          ompadjustN(type2,omptype)=ompadjustN(type2,omptype)+1
          nr=ompadjustN(type2,omptype)
          read(word(3),*,end=1000,err=1000)
     +      ompadjustE1(type2,omptype,nr)
          read(word(4),*,end=1000,err=1000)
     +      ompadjustE2(type2,omptype,nr)
          read(word(5),*,end=1000,err=1000) ompadjustD(type2,omptype,nr)
          read(word(6),*,iostat=istat) ompadjusts(type2,omptype,nr)
          if (istat.lt.0) goto 1000
          if (istat.gt.0) cycle
          cycle
        endif
        if (key.eq.'rescuefile') then
          read(word(2),*,end=1000,err=1000) mt
          if (mt.lt.1.or.mt.gt.nummt) goto 1050
          is=-1
          do k=1,71
            if (word(3)(k:k+1).eq.'_g') then
              is=0
              exit
            endif
            if (word(3)(k:k+1).eq.'_m') then
              is=1
              exit
            endif
            if (word(3)(k:k+1).eq.'_n') then
              is=2
              exit
            endif
            if (word(3)(k:k+1).eq.'_o') then
              is=3
              exit
            endif
            if (word(3)(k:k+1).eq.'_p') then
              is=4
              exit
            endif
            if (word(3)(k:k+1).eq.'_q') then
              is=5
              exit
            endif
            if (word(3)(k:k+1).eq.'_r') then
              is=6
              exit
            endif
            if (word(3)(k:k+1).eq.'_s') then
              is=7
              exit
            endif
            if (word(3)(k:k+1).eq.'_t') then
              is=8
              exit
            endif
            if (word(3)(k:k+1).eq.'_u') then
              is=9
              exit
            endif
            if (word(3)(k:k+1).eq.'_v') then
              is=10
              exit
            endif
          enddo
          rescuefile(mt,is)=word(3)
          val=1.
          read(word(4),*,end=230,err=1000) val
  230     grescue(mt,is)=val
          flagrescue=.true.
          cycle
        endif
        if (key.eq.'jlmmode') then
          read(value,*,end=1000,err=1000) jlmmode
          cycle
        endif
        if (key.eq.'lvadjust') then
          class=9
          call getvalues(class,word,Zix,Nix,type,
     +      ibar,irad,lval,igr,val,ival,cval,flagassign)
          if (flagassign) then
            lvadjust=val
            ompadjustp(1)=.true.
          endif
          cycle
        endif
        if (key.eq.'lwadjust') then
          class=9
          call getvalues(class,word,Zix,Nix,type,
     +      ibar,irad,lval,igr,val,ival,cval,flagassign)
          if (flagassign) then
            lwadjust=val
            ompadjustp(1)=.true.
          endif
          cycle
        endif
        if (key.eq.'lv1adjust') then
          class=9
          call getvalues(class,word,Zix,Nix,type,
     +      ibar,irad,lval,igr,val,ival,cval,flagassign)
          if (flagassign) then
            lv1adjust=val
            ompadjustp(1)=.true.
          endif
          cycle
        endif
        if (key.eq.'lw1adjust') then
          class=9
          call getvalues(class,word,Zix,Nix,type,
     +      ibar,irad,lval,igr,val,ival,cval,flagassign)
          if (flagassign) then
            lw1adjust=val
            ompadjustp(1)=.true.
          endif
          cycle
        endif
        if (key.eq.'lvsoadjust') then
          class=9
          call getvalues(class,word,Zix,Nix,type,
     +      ibar,irad,lval,igr,val,ival,cval,flagassign)
          if (flagassign) then
            lvsoadjust=val
            ompadjustp(1)=.true.
          endif
          cycle
        endif
        if (key.eq.'lwsoadjust') then
          class=9
          call getvalues(class,word,Zix,Nix,type,
     +      ibar,irad,lval,igr,val,ival,cval,flagassign)
          if (flagassign) then
            lwsoadjust=val
            ompadjustp(1)=.true.
          endif
          cycle
        endif
        if (key.eq.'aradialcor') then
          class=9
          call getvalues(class,word,Zix,Nix,type,
     +      ibar,irad,lval,igr,val,ival,cval,flagassign)
          if (flagassign) then
            aradialcor=val
            ompadjustp(6)=.true.
          endif
          cycle
        endif
        if (key.eq.'adepthcor') then
          class=9
          call getvalues(class,word,Zix,Nix,type,
     +      ibar,irad,lval,igr,val,ival,cval,flagassign)
          if (flagassign) then
            adepthcor=val
            ompadjustp(6)=.true.
          endif
          cycle
        endif
        if (key.eq.'spincutmodel') then
          read(value,*,end=1000,err=1000) spincutmodel
          cycle
        endif
        if (key.eq.'shellmodel') then
          read(value,*,end=1000,err=1000) shellmodel
          cycle
        endif
        if (key.eq.'kvibmodel') then
          read(value,*,end=1000,err=1000) kvibmodel
          cycle
        endif
        if (key.eq.'radialmodel') then
          read(value,*,end=1000,err=1000) radialmodel
          cycle
        endif
        if (key.eq.'breakupmodel') then
          read(value,*,end=1000,err=1000) breakupmodel
          cycle
        endif
        if (key.eq.'riplomp') then
          if (ch.eq.'n'.and.trim(word(3)).eq.'') then
            flagriplomp=.false.
          else
            class=7
            call getvalues(class,word,Zix,Nix,type,
     +        ibar,irad,lval,igr,val,ival,cval,flagassign)
            if (flagassign) riplomp(type)=ival
            flagriplomp=.true.
          endif
          cycle
        endif
        if (key.eq.'rspincut') then
          class=9
          call getvalues(class,word,Zix,Nix,type,
     +      ibar,irad,lval,igr,val,ival,cval,flagassign)
          if (flagassign) Rspincut=val
          cycle
        endif
        if (key.eq.'rspincutpreeq') then
          class=9
          call getvalues(class,word,Zix,Nix,type,
     +      ibar,irad,lval,igr,val,ival,cval,flagassign)
          if (flagassign) Rspincutpreeq=val
          cycle
        endif
        if (key.eq.'alphald') then
          read(value,*,end=1000,err=1000) alphaldall
          cycle
        endif
        if (key.eq.'betald') then
          read(value,*,end=1000,err=1000) betaldall
          cycle
        endif
        if (key.eq.'gammashell1') then
          read(value,*,end=1000,err=1000) gammashell1all
          cycle
        endif
        if (key.eq.'gammashell2') then
          read(value,*,end=1000,err=1000) gammashell2
          cycle
        endif
        if (key.eq.'pairconstant') then
          read(value,*,end=1000,err=1000) pairconstant
          cycle
        endif
        if (key.eq.'pshiftconstant') then
          read(value,*,end=1000,err=1000) Pshiftconstantall
          cycle
        endif
        if (key.eq.'kph') then
          read(value,*,end=1000,err=1000) Kph
          cycle
        endif
        if (key.eq.'m2constant') then
          class=9
          call getvalues(class,word,Zix,Nix,type,
     +      ibar,irad,lval,igr,val,ival,cval,flagassign)
          if (flagassign) then
            M2constant=val
            preeqadjust=.true.
          endif
          cycle
        endif
        if (key.eq.'m2limit') then
          read(value,*,end=1000,err=1000) M2limit
          cycle
        endif
        if (key.eq.'m2shift') then
          read(value,*,end=1000,err=1000) M2shift
          cycle
        endif
        if (key.eq.'rpipi') then
          read(value,*,end=1000,err=1000) Rpipi
          cycle
        endif
        if (key.eq.'rnunu') then
          read(value,*,end=1000,err=1000) Rnunu
          cycle
        endif
        if (key.eq.'rpinu') then
          read(value,*,end=1000,err=1000) Rpinu
          cycle
        endif
        if (key.eq.'rnupi') then
          read(value,*,end=1000,err=1000) Rnupi
          cycle
        endif
        if (key.eq.'esurf') then
          read(value,*,end=1000,err=1000) Esurf0
          cycle
        endif
        if (key.eq.'rgamma') then
          read(value,*,end=1000,err=1000) Rgamma
          cycle
        endif
        if (key.eq.'tmadjust') then
          class=9
          call getvalues(class,word,Zix,Nix,type,
     +      ibar,irad,lval,igr,val,ival,cval,flagassign)
          if (flagassign) Tmadjust=val
          cycle
        endif
        if (key.eq.'fsadjust') then
          class=9
          call getvalues(class,word,Zix,Nix,type,
     +      ibar,irad,lval,igr,val,ival,cval,flagassign)
          if (flagassign) Fsadjust=val
          cycle
        endif
        if (key.eq.'cnubar1') then
          read(value,*,end=1000,err=1000) Cnubar1
          cycle
        endif
        if (key.eq.'cnubar2') then
          read(value,*,end=1000,err=1000) Cnubar2
          cycle
        endif
        if (key.eq.'msdbins')  then
          read(value,*,end=1000,err=1000) msdbins
          cycle
        endif
        if (key.eq.'emsdmin')  then
          read(value,*,end=1000,err=1000) Emsdmin
          cycle
        endif
        if (key.eq.'elwidth')  then
          read(value,*,end=1000,err=1000) elwidth
          cycle
        endif
        if (key.eq.'xscaptherm')  then
          read(value,*,end=1000,err=1000) xscaptherm(-1)
          cycle
        endif
        if (key.eq.'xsptherm')  then
          read(value,*,end=1000,err=1000) xsptherm(-1)
          cycle
        endif
        if (key.eq.'xsalphatherm')  then
          read(value,*,end=1000,err=1000) xsalphatherm(-1)
          cycle
        endif
        if (key.eq.'anglesrec')  then
          read(value,*,end=1000,err=1000) nanglerec
          cycle
        endif
        if (key.eq.'rprime')  then
          read(value,*,end=1000,err=1000) RprimeU
          cycle
        endif
        if (key.eq.'astrot')  then
          read(value,*,end=1000,err=1000) astroT9
          nTmax=1
          flagastro=.true.
          cycle
        endif
        if (key.eq.'astroe')  then
          read(value,*,end=1000,err=1000) astroE
          nTmax=1
          flagastro=.true.
          flagastrogs=.true.
          cycle
        endif
        if (key.eq.'ebeam') then
          read(value,*,end=1000,err=1000) Ebeam
          cycle
        endif
        if (key.eq.'eback') then
          read(value,*,end=1000,err=1000) Eback
          cycle
        endif
        if (key.eq.'radiounit') then
          read(value,*,end=1000,err=1000) radiounit
          cycle
        endif
        if (key.eq.'yieldunit') then
          read(value,*,end=1000,err=1000) yieldunit
          cycle
        endif
        if (key.eq.'ibeam') then
          read(value,*,end=1000,err=1000) Ibeam
          cycle
        endif
        if (key.eq.'area') then
          read(value,*,end=1000,err=1000) Area
          cycle
        endif
        if (key.eq.'tirrad') then
          Tirrad(1)=0
          unitTirrad(1)=' '
          do k=1,5
            read(word(2*k),'(i9)',iostat=istat) Tirrad(k)
            if (istat.lt.0) goto 1000
            if (istat.gt.0) cycle
            read(word(2*k+1),'(a1)',end=1000,err=1000) unitTirrad(k)
          enddo
          cycle
        endif
        if (key.eq.'rho') then
          read(value,*,end=1000,err=1000) rhotarget
          cycle
        endif
        if (key.eq.'tcool') then
          Tcool(1)=0
          unitTcool(1)=' '
          do k=1,5
            read(word(2*k),'(i9)',iostat=istat) Tcool(k)
            if (istat.lt.0) goto 1000
            if (istat.gt.0) cycle
            read(word(2*k+1),'(a1)',end=1000,err=1000) unitTcool(k)
          enddo
          cycle
        endif
        cycle
 1040   write(*,'(" TALYS-warning: Z,N index out of range,",
     +    " keyword ignored: ",a)') trim(line)
      enddo
c
c Set level density parameters per nucleus
c
c alphaldall : variable for level density
c betaldall  : variable for level density
c gammashell1all: variable for level density
c Pshiftconstantall: variable for level density
c
      if (alphaldall.ne.-99.) alphald=alphaldall
      if (betaldall.ne.-99.) betald=betaldall
      if (gammashell1all.ne.-99.) gammashell1=gammashell1all
      if (Pshiftconstantall.ne.-99.) Pshiftconstant=Pshiftconstantall
c
c Apply consistent OMP adjustment factors for Koning-Delaroche and other
c potentials.
c
      do type=1,6
        if ((type.eq.3.and.deuteronomp.ge.2).or.(type.eq.6.and.
     +    alphaomp.ge.2)) then
          if (rwadjust(type).eq.-1.) rwadjust(type)=1.
          if (awadjust(type).eq.-1.) awadjust(type)=1.
          if (rvdadjust(type).eq.-1.) rvdadjust(type)=1.
          if (avdadjust(type).eq.-1.) avdadjust(type)=1.
          if (rwsoadjust(type).eq.-1.) rwsoadjust(type)=1.
          if (awsoadjust(type).eq.-1.) awsoadjust(type)=1.
        else
          if (rwadjust(type).eq.-1.) rwadjust(type)=rvadjust(type)
          if (awadjust(type).eq.-1.) awadjust(type)=avadjust(type)
          if (rvdadjust(type).eq.-1.) rvdadjust(type)=rwdadjust(type)
          if (avdadjust(type).eq.-1.) avdadjust(type)=awdadjust(type)
          if (rwsoadjust(type).eq.-1.) rwsoadjust(type)=rvsoadjust(type)
          if (awsoadjust(type).eq.-1.) awsoadjust(type)=avsoadjust(type)
        endif
      enddo
c
c In case of external PSF files, no default normalization of width
c
      if (.not.wtadjust) then
        do Zix=0,numZ
          do Nix=0,numN
            do irad=0,1
              do l=1,numgam
                if (Exlfile(Zix,Nix,irad,l)(1:1).ne.' ') 
     +            wtable(Zix,Nix,irad,l)=1.
              enddo
            enddo
          enddo
        enddo
      endif
      return
c
c Error and warning messages
c
 1000 write(*,'(" TALYS-error: Wrong input: ",a)') trim(line)
      stop
 1050 write(*,'(" TALYS-error: 1 <= MT number <= ",i4,
     +  ", MT index out of range: ",a)') nummt,trim(line)
      stop
 1060 write(*,'(" TALYS-error: 0 <= level number <= ",i4,
     +  ", ilev0, ilev1 or nbr index out of range: ",a)')
     +  numlev,trim(line)
      stop
 1070 write(*,'(" TALYS-error: 0 <= branching ratio",
     +  ", br index out of range: ",a)') trim(line)
      stop
      end
Copyright (C)  2023 A.J. Koning, S. Hilaire and S. Goriely
