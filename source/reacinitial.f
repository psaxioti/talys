      subroutine reacinitial
c
c +---------------------------------------------------------------------
c | Author: Arjan Koning
c | Date  : January 19, 2023
c | Task  : Initialization of arrays for various cross sections
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
c
c *************** Initialization of primary reaction *******************
c
c primary: flag to designate primary (binary) reaction
c
      primary=.true.
c
c *************** Initialize incident reaction arrays ******************
c
c lmaxinc     : maximal l-value for transmission coefficients for
c               incident channel
c xstotinc    : total cross section (neutrons only) for incident channel
c xsreacinc   : reaction cross section for incident channel
c xsoptinc    : optical model reaction cross section for incident
c               channel
c xsracape    : direct radiative capture cross section
c xselasinc   : total elastic cross section (neutrons only) for
c               incident channel
c xsinitpop   : initial population cross section
c numl        : maximal number of l-values
c Sstrength   : s,p,d,etc-wave strength function
c strengthl   : l neutron strength function
c strengthlj  : (l,j) neutron strength function
c sigurrs     : (l,j) scattering cross section for URR
c sigurrc     : (l,j) capture cross section for URR
c sigurrf     : (l,j) fission cross section for URR
c Tjlinc      : transmission coefficients as a function of spin
c               and l for the incident channel
c Turr        : (l,j) transmission coefficient for URR calculation
c Purrlj      : (l,j) parity for URR calculation
c xsbinarylj  : (l,j) cross section for URR calculation
c nulj        : (l,j) number of degrees of freedom for URR calculation
c urrwidth    : channel width in URR
c Tlinc       : transmission coefficients as a function of l for
c               the incident channel, averaged over spin
c numlev2     : maximum number of included discrete levels
c dleg        : direct reaction Legendre coefficient
c directad    : direct angular distribution
c ruth        : elastic/Rutherford ratio
c dorigin     : origin of direct cross section (Direct or Preeq)
c xsdirdisc   : direct cross section for discrete state
c xsdirdisctot: direct cross section summed over discrete states
c xsdirdiscsum: total direct cross section
c
      lmaxinc=0
      xstotinc=0.
      xsreacinc=0.
      xsracape=0.
      xsracapedisc=0.
      xsracapecont=0.
      xsoptinc=0.
      xselasinc=0.
      xsinitpop=0.
      Sstrength=0.
      strengthl=0.
      Tjlinc=0.
      strengthlj=0.
      sigurrs=0.
      sigurrc=0.
      sigurrf=0.
      Turrljinc=0.
      Purrlj=1
      Turrlj=0.
      xsbinarylj=0.
      nulj=0
      urrwidth=0.
      Tlinc=0.
      spot=0.
      dleg=0.
      directad=0.
      ruth=0.
      elasni=0.
      dorigin='      '
      xsdirdisc=0.
      xsdirdisctot=0.
      xsdirdiscsum=0.
c
c ********************** Initialization of energies ********************
c
c Nix    : neutron number index for residual nucleus
c numN   : maximal number of neutrons away from the initial
c          compound nucleus
c Zix    : charge number index for residual nucleus
c numZ   : maximal number of protons away from the initial
c          compound nucleus
c Exmax0 : maximum excitation energy for residual nucleus (including
c          negative energies)
c Exmax  : maximum excitation energy for residual nucleus
c maxex  : maximum excitation energy bin for compound nucleus
c numex  : maximal number of excitation energies
c Ex     : excitation energy
c Etotal : total energy of compound system (target + projectile)
c Exinc  : excitation energy of entrance bin
c deltaEx: excitation energy bin for population arrays
c
      Exmax0=0.
      Exmax=0.
      maxex=0
      Ex=0.
      Exmax0(0,0)=Etotal
      Exmax(0,0)=Etotal
      Exinc=Etotal
      deltaEx=0.
c
c ****************** Initialize giant resonance arrays *****************
c
c xscollconttot: total collective cross section in the continuum
c xsgrcoll     : giant resonance cross section
c eoutgr       : outgoing energy
c flagddx      : flag for output of double-differential cross sections
c numangcont   : maximum number of angles for continuum
c grcollad     : giant resonance angular distribution
c numen        : maximum number of outgoing energies
c xsgrstate    : smoothed giant resonance cross section per state
c xsgr         : smoothed giant resonance cross section
c xsgrad       : smoothed giant resonance angular distribution
c xsgrtot      : total smoothed giant resonance cross section
c collcontad   : collective angular distribution in the continuum
c xsgrsum      : sum over giant resonance cross sections
c xscollcont   : collective cross section in the continuum
c
      xsgrcoll=0.
      eoutgr=0.
      if (flagddx) grcollad=0.
      xsgrstate=0.
      xsgr=0.
      xsgrad=0.
      xscollconttot=0.
      xsgrtot=0.
      collcontad=0.
      xsgrsum=0.
      xscollcont=0.
      xscollcontJP=0.
c
c *************** Initialize pre-equilibrium arrays ********************
c
c numparx       : maximum number of particles
c xsstep        : preequilibrium cross section per particle type, stage
c                 and outgoing energy
c xsstep2       : two-component preequilibrium cross section
c wemission     : emission rate
c wemission2    : two-component emission rate
c xspreeq       : preequilibrium cross section per particle type and
c                 outgoing energy
c xsmpreeq      : multiple pre-equilibrium emission spectrum
c xspreeqps     : preequilibrium cross section per particle type and
c                 outgoing energy for pickup and stripping
c xspreeqki     : preequilibrium cross section per particle type and
c                 outgoing energy for knockout and inelastic
c xspreeqbu     : preequilibrium cross section per particle type and
c                 outgoing energy for breakup
c parity        : parity
c J             : total angular momentum
c numJph        : maximal spin for particle-hole states
c xspreeqJP     :   preequilibrium cross section per particle type,
c                 outgoing energy, spin and parity
c xspreeqtot    : preequilibrium cross section per particle type
c xspreeqtotps  : preequilibrium cross section per particle type for
c                 pickup and stripping
c xspreeqtotki  : preequilibrium cross section per particle type for
c                 knockout and inelastic
c xspreeqtotbu  : preequilibrium cross section per particle type for
c                 breakup
c xsBUnuc       : nucleon breakup cross section
c xsEB          : elastic breakup cross section
c xsBF          : nucleon inelastic breakup cross section
c xssteptot     : preequilibrium cross section per particle type and
c                 stage
c numJ          : maximal J-value
c xscompad      : compound emission angular distribution
c xsbinemisad   : angular distribution for emission from first compound
c                 nucleus
c flagang       : flag for output of angular distributions
c xspreeqad     : preequilibrium angular distribution per particle type
c xsmpreeqad    : multiple preequilibrium angular distribution
c xspreeqdisc   : preequilibrium cross section for discrete state
c xspreeqdisctot: preequilibrium cross section summed over discrete
c                 states
c xspreeqdiscsum: total preequilibrium cross section for discrete states
c flagmulpre    : flag for multiple pre-equilibrium calculation
c numNph        : maximal number of neutrons away from the initial
c                 compound nucleus for multiple pre-equilibrium emission
c numZph        : maximal number of protons away from the initial
c                 compound nucleus for multiple pre-equilibrium emission
c xspopph       : population cross section per particle-hole
c                 configuration
c PP2           : total strength
c Spre          : time-integrated strength of two-component exciton
c                 state
c xspopph2      : population cross section per two-component
c                 particle-hole configuration
c xspreeqsum    : total preequilibrium cross section summed over
c                 particles
c Esurf         : well depth for surface interaction
c n3            : counter
c n4            : counter
c
      xsstep=0.
      xsstep2=0.
      wemission=0.
      wemission2=0.
      xspreeq=0.
      xsmpreeq=0.
      xspreeqps=0.
      xspreeqki=0.
      xspreeqbu=0.
      xspreeqJP=0.
      xspreeqtot=0.
      xspreeqtotps=0.
      xspreeqtotki=0.
      xspreeqtotbu=0.
      xsBUnuc=0.
      xsEB=0.
      xsBF=0.
      xssteptot=0.
      if (flagang.or.flagddx) then
        xscompad=0.
        xsbinemisad=0.
        xspreeqad=0.
        xsmpreeqad=0.
      endif
      xspreeqdisc=0.
      xspreeqdisctot=0.
      xspreeqdiscsum=0.
      if (flagmulpre) then
        xspopph=0.
        PP2=0.
        Spre=0.
        xspopph2=0.
      endif
      xspreeqsum=0.
      Esurf=0.
c
c Multi-step direct arrays
c
c preeqmode : designator for pre-equilibrium model
c numenmsd  : maximum number of energy points for DWBA calculation for
c             MSD
c numJmsd   : maximal spin for MSD
c xsdwin    : DWBA cross section as a function of incident energy,
c             outgoing energy and angular momentum
c xsdw      : DWBA angular distribution as a function of incident
c             energy, outgoing energy, angular momentum and angle
c Emsd      : MSD energy grid
c xscont    : continuum one-step direct cross section for MSD
c xscont1   : continuum one-step direct cross section for MSD
c             (unnormalized)
c xscontad  : continuum one-step direct angular distribution for MSD
c xscontad1 : continuum one-step direct angular distribution for MSD
c             (unnormalized)
c msdstep1  : continuum one-step direct cross section (unnormalized)
c msdstepad1: continuum one-step direct angular distribution
c             (unnormalized)
c msdstep0  : n-step cross section for MSD
c msdstep   : continuum n-step direct cross section
c msdstepad0: n-step angular distribution for MSD
c msdstepad : continuum n-step direct angular distribution
c
      if (preeqmode.eq.4) then
        xsdwin=0.
        xsdw=0.
        Emsd=0.
        xscont=0.
        xscont1=0.
        xscontad=0.
        xscontad1=0.
        msdstep1=0.
        msdstepad1=0.
        msdstep=0.
        msdstepad=0.
        msdstep0=0.
        msdstepad0=0.
      endif
c
c ********* Initialization of primary compound and binary arrays *******
c
c J2beg      : 2 * start of J summation
c J2end      : 2 * end of J summation
c xspopnuc   : population cross section per nucleus
c xspopdir   : direct population cross section per nucleus
c xspoppreeq : preequilibrium population cross section per nucleus
c xspopcomp  : compound population cross section per nucleus
c Fdir       : direct population fraction per nucleus
c Fpreeq     : preequilibrium population fraction per nucleus
c Fcomp      : compound population fraction per nucleus
c xsBFnuc    : inelastic breakup enhancement brought by breakup neutrons
c              and protons interacting with the same Atarget and leading
c              to the same residual nucleus, (Z,N), populated by the
c              deuteron interaction process
c xspopnucT  : total population cross section per nucleus including the
c              inelastic breakup enhancement
c xspopex    : population cross section summed over spin and parity
c xsisoBU    : population cross section summed over spin and parity
c              including inelastic breakup enhancemeny            
c ENHratio   : breakup nucleons enhancing reaction cross
c              sections, PRC 89,044613 Eq. (7),
c              n + Atarget  sig(n,Z,A,Eout)/sig_Total(n,Enout);
c              p + Atarget  sig(p,Z,A,Eout)/sig_Reaction(p,Epout)
c xspopex    : population cross section summed over spin and parity
c maxJ       : maximal J-value
c xspop      : population cross section
c xscompel   : compound elastic cross section
c xselastot  : total elastic cross section (shape + compound)
c xsnonel    : non-elastic cross section
c xscompnonel: total compound non-elastic cross section
c xsbinary   : cross section from initial compound to residual nucleus
c xsdisctot  : total cross section summed over discrete states
c xsdircont  : direct cross section for continuum
c xsdirect   : total direct cross section
c xsconttot  : total cross section for continuum
c xscompound : total compound cross section
c xscompcont : compound cross section for continuum
c lmaxhf     : maximal l-value for transmission coefficients
c contrib    : contribution to emission spectrum
c feedbinary : feeding from first compound nucleus
c flagspec   : flag for output of spectra
c binemis    : emission spectra from initial compound nucleus
c xscomp     : compound emission spectrum
c xsbinemis  : cross section for emission from first compound nucleus
c xsemis     : emission spectrum from compound nucleus
c xsngnspec  : total (projectile,gamma-ejectile) spectrum
c numlev     : maximum number of included discrete levels
c xsracappop : population cross section for radiative capture
c xspopex0   : binary population cross section for discrete states
c preeqpopex : pre-equilibrium population cross section summed over
c              spin and parity
c preeqpop   : pre-equilibrium population cross section
c xsdisc     : total cross section for discrete state
c xscompdisc : compound cross section for discrete state
c numang     : maximum number of angles
c compad     : compound angular distribution
c discad     : discrete state angular distribution
c cleg       : compound nucleus Legendre coefficient
c tleg       : total Legendre coefficient
c tlegnor    : total Legendre coefficient normalized to 1
c cleg0      : Legendre coefficient normalized to the first one
c transjl    : array for width fluctuation calculation
c
      J2beg=0
      J2end=0
      xspopnuc=0.
      xspopdir=0.
      xspoppreeq=0.
      xspopcomp=0.
      Fdir=0.
      Fpreeq=0.
      Fcomp=0.
      xsBFnuc=0.
      xspopnucT=0.
      xspopnucP=0.
      CNterm=0.
      xsracappopex=0.
      xspopex=0.
      xsisoBU=0.
      preeqpopex=0.
      maxJ=numJ
      xspopexP=0.
      xsracappop=0.
      xspop=0.
      preeqpop=0.
      xscompel=0.
      xselastot=0.
      xsnonel=0.
      xscompnonel=0.
      xsbinary=0.
      xsdisctot=0.
      xsdircont=0.
      xsdirect=0.
      xsconttot=0.
      xscompound=0.
      xscompcont=0.
      Eaveragebin=0.
      Eaverage=0.
      lmaxhf=0
      contrib=0.
      feedbinary=0.
      if (flagspec) then
        binemis=0.
        xscomp=0.
        xsbinemis=0.
        xsemis=0.
        xsngnspec=0.
      endif
      xspopex0=0.
      xsdisc=0.
      xscompdisc=0.
      if (flagang.or.flagddx) then
        compad=0.
        discad=0.
      endif
      cleg=0.
      tleg=0.
      tlegnor=0.
      cleg0=0.
      transjl=0.
c
c ************* Initialization of multiple mission arrays **************
c
c xsngnsum    : sum over total (projectile,gamma-ejectile) cross
c               sections
c idnum       : counter for exclusive channel
c flagchannels: flag for exclusive channels calculation
c channelsum  : sum over exclusive channel cross sections
c xsabs       : absorption cross section
c numNchan    : maximal number of outgoing neutron units in individual
c               channel description
c numZchan    : maximal number of outgoing proton units in individual
c               channel description
c feedexcl    : feeding terms from compound excitation energy bin to
c               residual excitation energy bin
c popexcl     : population cross section of bin just before decay
c flagfission : flag for fission
c fisfeedex   : fission contribution from excitation energy bin
c tfis        : fission transmission coefficients
c xsfeed      : cross section from compound to residual nucleus
c xsngn       : total (projectile,gamma-ejectile) cross section
c xsgamdis    : discrete gamma-ray cross section
c xsgamdistot : total discrete gamma-ray cross section
c
      xsngnsum=0.
      idnum=-1
      if (flagchannels) then
        channelsum=0.
        xsabs=0.
        feedexcl=0.
        popexcl=0.
      endif
      if (flagfission) then
        denfis=0.
        gamfis=0.
        taufis=0.
        tfisA=0.
        rhofisA=0.
        tfisdown=0.
        tfisup=0.
        fisfeedex=0.
        fisfeedJP=0.
      endif
      tfis=0.
      xsfeed=0.
      Eaveragemul=0.
      xsngn=0.
      xsgamdis=0.
      xsgamdistot=0.
c
c ************* Initialization of channels and total arrays ************
c
c nin          : counter for incident energy
c idnumfull    : flag to designate maximum number of exclusive channels
c opennum      : total number of open channels
c numin,....   : maximal number of ejectile in channel description
c chanopen     : flag to open channel with first non-zero cross section
c numen2       : maximum number of outgoing energies
c xssumoutad   : angular distribution summed over mechanisms
c xsdiscoutad  : smoothed angular distribution for discrete state
c xscompoutad  : compound emission angular distribution
c xspreeqoutad : preequilibrium angular distribution per particle type
c xsmpreeqoutad: multiple preequilibrium angular distribution
c xsbranch     : branching ratio for isomeric cross section
c xsparcheck   : total particle production cross section
c xsspeccheck  : total particle production spectra
c
      if (nin.eq.1) then
        idnumfull=.false.
        opennum=-1
        chanopen=.false.
      endif
      if (flagang.or.flagddx) then
        xssumoutad=0.
        xsdiscoutad=0.
        xscompoutad=0.
        xspreeqoutad=0.
        xsmpreeqoutad=0.
      endif
      xsbranch=0.
      xsparcheck=0.
      xsspeccheck=0.
      espec=0.
      xsmpreeqout=0.
      xssumout=0.
      xscompout=0.
      xsdiscout=0.
      xspreeqout=0.
      xspreeqpsout=0.
      xspreeqkiout=0.
      xspreeqbuout=0.
      preeqratio=0.
      buratio=0.
c
c ************* Initialization of total cross section arrays ***********
c
c nexmax       : maximum excitation energy bin for residual nucleus
c xsexclusive  : exclusive single channel cross section
c xsexclcont   : exclusive single channel cross section for continuum
c multiplicity : particle multiplicity
c xsparticle   : total particle production cross section
c xsfistot     : total fission cross section
c maxA         : maximal number of nucleons away from the initial
c                compound nucleus
c maxZ         : maximal number of protons away from the initial
c                compound nucleus
c maxN         : maximal number of protons away from the initial
c                compound nucleus
c xsresprod    : total residual production (= reaction) cross section
c xsmassprod   : residual production cross section per mass unit
c xsmassprodT  : residual production cross section per mass unit
c                including inelastic breakup enhancement
c xstot6       : total cross section (neutrons only) for ENDF-6 file
c xsreac6      : reaction cross section for ENDF-6 file
c xsnon6       : non-elastic cross section for ENDF-6 file
c xsopt6       : optical model reaction cross section for ENDF-6 file
c xselassh6    : shape elastic cross section (neutrons only) for ENDF-6
c                file
c xselas6      : total elastic cross section (neutrons only) for ENDF-6
c                file
c xscompel6    : compound elastic cross section
c e6           : energies of ENDF-6 energy grid in MeV
c fxsbinary    : cross section from initial compound to residual nucleus
c fxsexclusive : exclusive single channel cross section
c fxsdisctot   : total cross section summed over discrete states
c fxsexclcont  : exclusive single channel cross section for continuum
c fxsngn       : total (projectile,gamma-ejectile) cross section
c fxsdisc      : total cross section for discrete state
c fxsdirdisc   : direct cross section for discrete state
c fxscompdisc  : compound cross section for discrete state
c fxspopnuc    : population cross section per nucleus
c fxspopex     : population cross section summed over spin and parity
c fxsisoBU     : population cross section summed over spin and parity
c                including inelastic breakup enhancement
c fxsbranch    : branching ratio for isomeric cross section
c idchannel    : identifier for exclusive channel
c fxschannel   : channel cross section
c fxsgamchannel: gamma channel cross section
c fxsratio     : ratio of exclusive cross section over residual
c                production cross section (for exclusive gamma ray
c                intensities)
c fxsgamdischan: discrete gamma channel cross section
c fxschaniso   : channel cross section per isomer
c fexclbranch  : exclusive channel yield per isomer
c fxsgamdischan: discrete gamma channel cross section
c fxsnonel     : non-elastic cross section
c fxselastot   : total elastic cross section (shape + compound)
c fxstotinc    : total cross section (neutrons only) for incident
c                channel
c fxscompel    : compound elastic cross section
c fxselasinc   : total elastic cross section (neutrons only) for
c                incident channel
c fxsreacinc   : reaction cross section for incident channel
c fxscompnonel : total compound non-elastic cross section
c fxsdirdiscsum: total direct cross section
c fxspreeqsum  : total preequilibrium cross section summed over
c                particles
c fxsracape    : direct capture cross section
c fisstring    : string for exclusive fission reaction channel
c
      nexmax=-1
      xsexclusive=0.
      xsexclcont=0.
      multiplicity=0.
      xsparticle=0.
      xsfistot=0.
      if (.not.flagffruns) xsfistot0=0.
      if (.not.flagrpruns) xspopnuc0=0.
      maxA=maxZ+maxN
      xsresprod=0.
      xsmassprod=0.
      xsmassprodT=0.
      if (nin.eq.1) then
        xstot6=0.
        xsreac6=0.
        xsnon6=0.
        xsnonel6=0.
        xsopt6=0.
        xselassh6=0.
        xselas6=0.
        xscompel6=0.
        e6=0.
        fxsbinary=0.
        fxsexclusive=0.
        fxsdisctot=0.
        fxsexclcont=0.
        fxsngn=0.
        fnubar=0.
        fxsdisc=0.
        fxsdirdisc=0.
        fxscompdisc=0.
        fxsngn=0.
        fxsreacinc=0.
        fxselastot=0.
        fxsnonel=0.
        fxscompel=0.
        fxselasinc=0.
        fxstotinc=0.
        fxscompnonel=0.
        fxsdirdiscsum=0.
        fxspreeqsum=0.
        fxsracape=0.
        fxspopnuc=0.
        fxspopex=0.
        fxsisoBU=0.
        fxsbranch=0.
        idchannel=-1
        fxschannel=0.
        fxsgamchannel=0.
        fxsratio=0.
      endif
      fxsgamdischan=0.
      fxschaniso=0.
      fexclbranch=0.
      fisstring='                  '
      return
      end
Copyright (C)  2023 A.J. Koning, S. Hilaire and S. Goriely
