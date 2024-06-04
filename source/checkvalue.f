      subroutine checkvalue
c
c +---------------------------------------------------------------------
c | Author: Arjan Koning
c | Date  : December 18, 2007
c | Task  : Check for errors in values
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      logical lexist
      integer type,i,Zix,Nix,A,irad,lval,igr,ibar,fax0
      real    egr0,ggr0,sgr0,value,fbar0,fhw0,fR0
c
c All parameters need to fall within certain ranges. These ranges are 
c specified in this subroutine and in the manual.
c
c ******************* Check for wrong input variables ******************
c
c 1. Check of values for four main keywords.
c
c ptype0  : type of incident particle
c parsym  : symbol of particle
c numelem : number of elements   
c Starget : symbol of target nucleus
c nuc     : symbol of nucleus
c Atarget : mass number of target nucleus
c nummass : number of masses
c Zinit   : charge number of initial compound nucleus
c Ninit   : neutron number of initial compound nucleus
c enincmin: minimum incident energy
c enincmax: maximum incident energy
c
      do 10 type=0,6
        if (ptype0.eq.parsym(type)) goto 20
   10 continue
      if (ptype0.eq.'0') goto 20
      write(*,'(" TALYS-error: Wrong symbol for projectile: ",a1)') 
     +  ptype0
      stop
   20 do 30 i=3,numelem
        if (Starget.eq.nuc(i)) goto 40
   30 continue
      write(*,'(" TALYS-error: Wrong symbol for element: ",a2)') Starget
      stop
   40 if (Atarget.le.5.or.Atarget.gt.nummass) then
        write(*,'(" TALYS-error: 5 < Target mass <= ",i3)') nummass
        stop
      endif
      if (Zinit.le.2) then
        write(*,'(" TALYS-error: Target Z > 2")')
        stop
      endif
      if (Ninit.le.2) then
        write(*,'(" TALYS-error: Target N > 2")')
        stop
      endif
      if (Zinit.gt.numelem) then
        write(*,'(" TALYS-error: combination of projectile ",
     +    "+ element > ",i3)') numelem
        stop
      endif
      if (enincmin.lt.1.e-11.or.enincmax.ge.250.) then
        write(*,'(" TALYS-error: 1.e-5 eV <= Incident energy",
     +    " < 250 MeV")')
        stop
      endif
c
c 2. Check of values for basic physical and numerical parameters
c
c maxZ,numZ    : maximal number of protons away from the initial
c                compound nucleus
c maxN,numN    : maximal number of neutrons away from the initial
c                compound nucleus
c nbins,numbins: number of continuum excitation energy bins
c numbins      : maximal number of continuum excitation energy bins
c ptype0       : type of incident particle
c segment      : number of segments to divide emission energy grid
c flagastro    : flag for calculation of astrophysics reaction rate
c nlevmax      : maximum number of included discrete levels for target
c nlevmaxres   : maximum number of included discrete levels for residual
c                nucleus  
c numlev       : maximum number of included discrete levels 
c nlevbin      : number of excited levels for binary nucleus
c nlev         : number of excited levels for nucleus  
c Zix          : charge number index for residual nucleus
c numZ         : maximal number of protons away from initial compound 
c                nucleus
c Nix          : neutron number index for residual nucleus           
c numN         : maximal number of neutrons away from initial compound 
c                nucleus
c A            : mass number
c Ainit        : mass number of initial compound nucleus
c massnucleus  : mass of nucleus in amu as read from user input file
c massexcess   : mass excess in MeV as read from user input file
c Ltarget      : excited level of target
c isomer       : definition of isomer in seconds
c core         : even-even core for weakcoupling (-1 or 1)  
c transpower   : power for transmission coefficient limit
c transeps     : absolute limit for transmission coefficient
c xseps        : limit for cross sections
c popeps       : limit for population cross section per nucleus
c Rfiseps      : ratio for limit for fission cross section per nucleus 
c eninclow     : minimal incident energy for nuclear model calculations
c nangle       : number of angles
c numang       : maximum number of angles
c nanglecont   : number of angles for continuum
c numangcont   : maximum number of angles for continuum
c nanglerec    : number of recoil angles   
c numangrec    : maximum number of recoil angles   
c maxenrec     : number of recoil energies 
c numenrec     : maximum number of recoil energies 
c maxchannel   : maximal number of outgoing particles in individual
c                channel description (e.g. this is 3 for (n,2np))
c massmodel    : model for theoretical nuclear mass
c
      if (maxZ.lt.0.or.maxZ.gt.numZ-2) then
        write(*,'(" TALYS-error: 0 <= maxZ <=",i3)') numZ-2
        stop
      endif
      if (maxN.lt.0.or.maxN.gt.numN-2) then
        write(*,'(" TALYS-error: 0 <= maxN <=",i3)') numN-2
        stop
      endif
      if (nbins.lt.2.or.nbins.gt.numbins) then
        write(*,'(" TALYS-error: 2 <= bins <=",i3)') numbins
        stop
      endif
      if (segment.le.0.or.segment.gt.4) then
        write(*,'(" TALYS-error: 1 <= segment <= 4")') 
        stop
      else
        if (segment.gt.1.and.enincmax.gt.100.) then
          write(*,'(" TALYS-error: segment = 1",
     +      " for incident energy of ",f7.3," MeV")') enincmax
          stop
        endif
        if (segment.gt.2.and.enincmax.gt.40.) then
          write(*,'(" TALYS-error: 1 <= segment = 2",
     +      " for incident energy of ",f7.3," MeV")') enincmax
          stop
        endif
        if (segment.gt.3.and.enincmax.gt.20.) then
          write(*,'(" TALYS-error: 1 <= segment = 3",
     +      " for incident energy of ",f7.3," MeV")') enincmax
          stop
        endif
        if (segment.gt.1.and.flagastro) then
          write(*,'(" TALYS-error: segment = 1",
     +      " for astrophysical calculations")')
          stop
        endif
      endif
      if (nlevmax.lt.0.or.nlevmax.gt.numlev) then
        write(*,'(" TALYS-error: 0 <= maxlevelstar <=",i3)') numlev
        stop
      endif
      if (nlevmaxres.lt.0.or.nlevmaxres.gt.numlev) then
        write(*,'(" TALYS-error: 0 <= maxlevelsres <=",i3)') numlev
        stop
      endif
      do 50 type=0,6
        if (nlevbin(type).lt.0.or.nlevbin(type).gt.numlev) then
          write(*,'(" TALYS-error: 0 <= maxlevelsbin <=",i3)') numlev
          stop
        endif
   50 continue
      do 60 Zix=0,numZ
        do 70 Nix=0,numN
          if (nlev(Zix,Nix).lt.0.or.nlev(Zix,Nix).gt.numlev) then
            write(*,'(" TALYS-error: 0 <= nlevels <=",i3)') numlev
            stop
          endif
          A=Ainit-Zix-Nix
          if (massnucleus(Zix,Nix).ne.0.and.
     +      (massnucleus(Zix,Nix).lt.real(A)-1..or.
     +      massnucleus(Zix,Nix).gt.real(A)+1.)) then
            write(*,'(" TALYS-error: A-1. <= massnucleus <= A+1.")')
            stop
          endif
          if (massexcess(Zix,Nix).ne.0.and.
     +      (massexcess(Zix,Nix).lt.-1000..or.
     +      massexcess(Zix,Nix).gt.1000.)) then
            write(*,'(" TALYS-error: -1000. <= massexcess <= 1000.")')
            stop
          endif
   70   continue
   60 continue
      if (Ltarget.lt.0.or.Ltarget.gt.numlev) then
        write(*,'(" TALYS-error: 0 <= Ltarget <=",i3)') numlev
        stop
      endif
      if (isomer.lt.0..or.isomer.gt.1.e38) then
        write(*,'(" TALYS-error: 0. <= isomer <= 1.e38")') 
        stop
      endif
      if (core.ne.-1.and.core.ne.1) then
        write(*,'(" TALYS-error: core = -1 or 1")')
        stop
      endif
      if (transpower.lt.2.or.transpower.gt.20) then
        write(*,'(" TALYS-error: 2 <= transpower <= 20")') 
        stop
      endif
      if (transeps.lt.0..or.transeps.gt.1.) then
        write(*,'(" TALYS-error: 0. <= transeps <= 1.")') 
        stop
      endif
      if (xseps.lt.0..or.xseps.gt.1000.) then
        write(*,'(" TALYS-error: 0. <= xseps <= 1000.")') 
        stop
      endif
      if (popeps.lt.0..or.popeps.gt.1000.) then
        write(*,'(" TALYS-error: 0. <= popeps <= 1000.")') 
        stop
      endif
      if (Rfiseps.lt.0..or.Rfiseps.gt.1.) then
        write(*,'(" TALYS-error: 0. <= Rfiseps <= 1.")') 
        stop
      endif
      if (eninclow.ne.0..and.(eninclow.lt.1.e-6.or.eninclow.gt.1.)) 
     +  then
        write(*,'(" TALYS-error: 1.e-6 <= Elow <= 1.")') 
        stop
      endif
      if (nangle.lt.1..or.nangle.gt.numang) then
        write(*,'(" TALYS-error: 1 <= angles <=",i3)') numang 
        stop
      endif
      if (nanglecont.lt.1..or.nanglecont.gt.numangcont) then
        write(*,'(" TALYS-error: 1 <= anglescont <=",i3)') numangcont
        stop
      endif
      if (nanglerec.lt.1..or.nanglerec.gt.numangrec) then
        write(*,'(" TALYS-error: 1 <= anglesrec <=",i3)') numangrec
        stop
      endif
      if (maxenrec.lt.1..or.maxenrec.gt.numenrec) then
        write(*,'(" TALYS-error: 1 <= maxenrec <=",i3)') numenrec
        stop
      endif
      if (maxchannel.lt.1.or.maxchannel.gt.8) then
        write(*,'(" TALYS-error: 1 <= maxchannel <= 8")')
        stop
      endif
      if (massmodel.lt.1.or.massmodel.gt.3) then
        write(*,'(" TALYS-error: 1 <= massmodel <= 3")')
        stop
      endif
c
c 3. Check of values of optical model
c
c numNph     : maximal number of neutrons away from the initial 
c              compound nucleus for multiple pre-equilibrium emission
c numZph     : maximal number of protons away from the initial 
c              compound nucleus for multiple pre-equilibrium emission
c optmod     : file with optical model parameters
c optmodfileN: optical model parameter file for neutrons
c optmodfileP: optical model parameter file for protons
c radialfile : radial matter density file
c
      do 110 Zix=0,numZph
        do 120 Nix=0,numNph
          do 130 type=1,6
            if (optmod(Zix,Nix,type)(1:1).eq.' ') goto 130
            inquire (file=optmod(Zix,Nix,type),exist=lexist)
            if (.not.lexist) then
              write(*,'(" TALYS-error: Non-existent optical model",
     +          " file: ",a72)') optmod(Zix,Nix,type)
              stop
            endif
  130     continue
  120   continue
        if (optmodfileN(Zix)(1:1).ne.' ') then
          inquire (file=optmodfileN(Zix),exist=lexist)
          if (.not.lexist) then
            write(*,'(" TALYS-error: Non-existent optical model",
     +        " file: ",a72)') optmodfileN(Zix)
            stop
          endif
        endif
        if (optmodfileP(Zix)(1:1).ne.' ') then
          inquire (file=optmodfileP(Zix),exist=lexist)
          if (.not.lexist) then
            write(*,'(" TALYS-error: Non-existent optical model",
     +        " file: ",a72)') optmodfileP(Zix)
            stop
          endif
        endif
        if (radialfile(Zix)(1:1).ne.' ') then
          inquire (file=radialfile(Zix),exist=lexist)
          if (.not.lexist) then
            write(*,'(" TALYS-error: Non-existent radial file: ",a72)')
     +        radialfile(Zix)
            stop
          endif
        endif
c
c Check other parameter input files
c
c levelfile  : discrete level file
c deformfile : deformation parameter file
c hbtransfile: file with head band transition states
c class2file : file with class 2 transition states
c alphaomp   : alpha optical model (1=normal, 2= McFadden-Satchler)
c radialmodel: model for radial matter densities (JLM OMP only)
c
        if (levelfile(Zix)(1:1).ne.' ') then
          inquire (file=levelfile(Zix),exist=lexist)
          if (.not.lexist) then
            write(*,'(" TALYS-error: Non-existent level file: ",a72)')
     +        levelfile(Zix)
            stop
          endif
        endif
        if (deformfile(Zix)(1:1).ne.' ') then
          inquire (file=deformfile(Zix),exist=lexist)
          if (.not.lexist) then
            write(*,'(" TALYS-error: Non-existent deformation ",
     +        "parameter file: ",a72)') deformfile(Zix)
            stop
          endif
        endif
        do 140 Nix=0,numN
          if (hbtransfile(Zix,Nix)(1:1).ne.' ') then
            inquire (file=hbtransfile(Zix,Nix),exist=lexist)
            if (.not.lexist) then
              write(*,'(" TALYS-error: Non-existent head band ",
     +          "transition state file: ",a72)') 
     +          hbtransfile(Zix,Nix)
              stop
            endif
          endif
          if (class2file(Zix,Nix)(1:1).ne.' ') then
            inquire (file=class2file(Zix,Nix),exist=lexist)
            if (.not.lexist) then
              write(*,'(" TALYS-error: Non-existent class 2 ",
     +          "transition state file: ",a72)') class2file(Zix,Nix)
              stop
            endif
          endif
  140   continue
  110 continue
      if (alphaomp.lt.1.or.alphaomp.gt.2) then
        write(*,'(" TALYS-error: 1 <= alphaomp <= 2")')
        stop
      endif
      if (radialmodel.lt.1.or.radialmodel.gt.2) then
        write(*,'(" TALYS-error: 1 <= radialmodel <= 2")')
        stop
      endif
c
c Check adjustable OMP parameters
c
c v1adjust..: adjustable factors for OMP (default 1.)
c
      do 150 type=1,6
        if (v1adjust(type).lt.0.2.or.v1adjust(type).gt.5.) then
          write(*,'(" TALYS-error: 0.2 <= v1adjust <= 5.")')
          stop
        endif
        if (v2adjust(type).lt.0.2.or.v2adjust(type).gt.5.) then
          write(*,'(" TALYS-error: 0.2 <= v2adjust <= 5.")')
          stop
        endif
        if (v3adjust(type).lt.0.2.or.v3adjust(type).gt.5.) then
          write(*,'(" TALYS-error: 0.2 <= v3adjust <= 5.")')
          stop
        endif
        if (v4adjust(type).lt.0.2.or.v4adjust(type).gt.5.) then
          write(*,'(" TALYS-error: 0.2 <= v4adjust <= 5.")')
          stop
        endif
        if (rvadjust(type).lt.0.5.or.rvadjust(type).gt.2.) then
          write(*,'(" TALYS-error: 0.5 <= rvadjust <= 2.")')
          stop
        endif
        if (avadjust(type).lt.0.5.or.avadjust(type).gt.2.) then
          write(*,'(" TALYS-error: 0.5 <= avadjust <= 2.")')
          stop
        endif
        if (w1adjust(type).lt.0.2.or.w1adjust(type).gt.5.) then
          write(*,'(" TALYS-error: 0.2 <= w1adjust <= 5.")')
          stop
        endif
        if (w2adjust(type).lt.0.2.or.w2adjust(type).gt.5.) then
          write(*,'(" TALYS-error: 0.2 <= w2adjust <= 5.")')
          stop
        endif
        if (d1adjust(type).lt.0.2.or.d1adjust(type).gt.5.) then
          write(*,'(" TALYS-error: 0.2 <= d1adjust <= 5.")')
          stop
        endif
        if (d2adjust(type).lt.0.2.or.d2adjust(type).gt.5.) then
          write(*,'(" TALYS-error: 0.2 <= d2adjust <= 5.")')
          stop
        endif
        if (d3adjust(type).lt.0.2.or.d3adjust(type).gt.5.) then
          write(*,'(" TALYS-error: 0.2 <= d3adjust <= 5.")')
          stop
        endif
        if (rvdadjust(type).lt.0.5.or.rvdadjust(type).gt.2.) then
          write(*,'(" TALYS-error: 0.5 <= rvdadjust <= 2.")')
          stop
        endif
        if (avdadjust(type).lt.0.5.or.avdadjust(type).gt.2.) then
          write(*,'(" TALYS-error: 0.5 <= avdadjust <= 2.")')
          stop
        endif
        if (vso1adjust(type).lt.0.2.or.vso1adjust(type).gt.5.) then
          write(*,'(" TALYS-error: 0.2 <= vso1adjust <= 5.")')
          stop
        endif
        if (vso2adjust(type).lt.0.2.or.vso2adjust(type).gt.5.) then
          write(*,'(" TALYS-error: 0.2 <= vso2adjust <= 5.")')
          stop
        endif
        if (wso1adjust(type).lt.0.2.or.wso1adjust(type).gt.5.) then
          write(*,'(" TALYS-error: 0.2 <= wso1adjust <= 5.")')
          stop
        endif
        if (wso2adjust(type).lt.0.2.or.wso2adjust(type).gt.5.) then
          write(*,'(" TALYS-error: 0.2 <= wso2adjust <= 5.")')
          stop
        endif
        if (rvsoadjust(type).lt.0.5.or.rvsoadjust(type).gt.2.) then
          write(*,'(" TALYS-error: 0.5 <= rvsoadjust <= 2.")')
          stop
        endif
        if (avsoadjust(type).lt.0.5.or.avsoadjust(type).gt.2.) then
          write(*,'(" TALYS-error: 0.5 <= avsoadjust <= 2.")')
          stop
        endif
        if (rcadjust(type).lt.0.5.or.rcadjust(type).gt.2.) then
          write(*,'(" TALYS-error: 0.5 <= rcadjust <= 2.")')
          stop
        endif
  150 continue
      if (lvadjust.lt.0.5.or.lvadjust.gt.1.5) then
        write(*,'(" TALYS-error: 0.2 <= lvadjust <= 5.")')
        stop
      endif
      if (lwadjust.lt.0.5.or.lwadjust.gt.1.5) then
        write(*,'(" TALYS-error: 0.5 <= lwadjust <= 1.5")')
        stop
      endif
      if (lv1adjust.lt.0.5.or.lv1adjust.gt.1.5) then
        write(*,'(" TALYS-error: 0.5 <= lv1adjust <= 1.5")')
        stop
      endif
      if (lw1adjust.lt.0.5.or.lw1adjust.gt.1.5) then
        write(*,'(" TALYS-error: 0.5 <= lw1adjust <= 1.5")')
        stop
      endif
      if (lvsoadjust.lt.0.5.or.lvsoadjust.gt.1.5) then
        write(*,'(" TALYS-error: 0.5 <= lvsoadjust <= 1.5")')
        stop
      endif
      if (lwsoadjust.lt.0.5.or.lwsoadjust.gt.1.5) then
        write(*,'(" TALYS-error: 0.5 <= lwsoadjust <= 1.5")')
        stop
      endif
c
c Check direct reaction parameters
c
c maxband   : highest vibrational band added to rotational model
c maxrot    : number of included excited rotational levels
c k0        : index for incident particle
c flaggiant0: flag for collective contribution from giant resonances
c
      if (maxband.lt.0.or.maxband.gt.100) then
        write(*,'(" TALYS-error: 0 <= maxband <= 100")')
        stop
      endif
      if (maxrot.lt.0.or.maxrot.gt.10) then
        write(*,'(" TALYS-error: 0 <= maxrot <= 10")')
        stop
      endif
      if (k0.eq.0.and.flaggiant0) then
        write(*,'(" TALYS-error: No giant resonance sumrules for",
     +    " photonuclear reactions")')
        stop
      endif
c
c 4. Check of values for compound nucleus
c
c ewfc        : off-set incident energy for width fluctuation 
c               calculation
c wmode       : designator for width fluctuation model
c flageciscomp: flag for compound nucleus calculation by ECIS 
c
      if ((ewfc.lt.0..or.ewfc.gt.20.).and.ewfc.ne.-1.) then
        write(*,'(" TALYS-error: 0. <= ewfc <= 20.")') 
        stop
      endif
      if (wmode.lt.0.or.wmode.gt.3) then
        write(*,'(" TALYS-error: 0 <= wmode <= 3")')
        stop
      endif
      if (k0.eq.0.and.ewfc.gt.0.) then
        write(*,'(" TALYS-error: No width fluctuations for",
     +    " photonuclear reactions")')
        stop
      endif
      if (k0.eq.0..and.flageciscomp) then
        write(*,'(" TALYS-error: No compound calculation by",
     +    " ECIS for incident photons")')
        stop
      endif
      if (enincmax.gt.20..and.flageciscomp) then
        write(*,'(" TALYS-error: No compound calculation by",
     +    " ECIS for E > 20 MeV")')
        stop
      endif
c
c 5. Check of values for gamma emission  
c
c gammax  : number of l-values for gamma multipolarity
c strength: model for E1 gamma-ray strength function
c egr,egr0: energy of GR
c irad    : variable to indicate M(=0) or E(=1) radiation
c lval    : multipolarity              
c igr     : giant resonance
c ggr,ggr0: width of GR
c sgr,sgr0: strength of GR    
c gamgam  : experimental total radiative width in eV
c D0      : experimental s-wave resonance spacing in eV
c S0      : s-wave strength function          
c gnorm   : gamma normalization factor
c etable  : constant to adjust tabulated strength functions
c ftable  : constant to adjust tabulated strength functions
c
      if (gammax.lt.1.or.gammax.gt.6) then
        write(*,'(" TALYS-error: 1 <= gammax <= 6")')
        stop
      endif
      if (strength.lt.1.or.strength.gt.4) then
        write(*,'(" TALYS-error: 1 <= strength <= 4")')
        stop
      endif
      if (strength.ge.3) then
        inquire (file=path(1:lenpath)//'gamma/hfb/z050',exist=lexist)
        if (.not.lexist) then
          write(*,'(" TALYS-error: Microscopic HFB tables are not ",
     +    " installed: download the full TALYS package from ",
     +    "www.talys.eu")')
          stop
        endif
      endif
      do 210 Zix=0,numZ
        do 220 Nix=0,numN
          do 230 irad=0,1        
            do 240 lval=1,gammax
              do 250 igr=1,2
                egr0=egr(Zix,Nix,irad,lval,igr)
                if (egr0.ne.0..and.(egr0.lt.1..or.egr0.gt.100.)) then
                  write(*,'(" TALYS-error: 1.<= energy of GR <=100.")')
                  stop
                endif
                ggr0=ggr(Zix,Nix,irad,lval,igr)
                if (ggr0.ne.0..and.(ggr0.lt.1..or.ggr0.gt.100.)) then
                  write(*,'(" TALYS-error: 1.<= width of GR <=100.")')
                  stop
                endif
                sgr0=sgr(Zix,Nix,irad,lval,igr)
                if (sgr0.lt.0..or.sgr0.gt.10000.) then
                  write(*,'(" TALYS-error: 0.<= strength of GR",
     +              " <=10000.")')
                  stop
                endif
  250         continue
  240       continue
  230     continue
          if (gamgam(Zix,Nix).ne.0.) then
            if (gamgam(Zix,Nix).lt.0..or.gamgam(Zix,Nix).gt.10.) then
              write(*,'(" TALYS-error: 0. < gamgam <= 10.")')
              stop
            endif
          endif
          if (D0(Zix,Nix).ne.0.) then
            if (D0(Zix,Nix).lt.1.e-6.or.D0(Zix,Nix).gt.1.e7) then
              write(*,'(" TALYS-error: 1.e-6 <= D0 <= 1.e7")')
              stop
            endif
          endif
          if (S0(Zix,Nix).lt.0..or.S0(Zix,Nix).gt.10.) then
            write(*,'(" TALYS-error: 0. <= S0 <= 10.")')
            stop
          endif
          if (etable(Zix,Nix).lt.-10..or.etable(Zix,Nix).gt.10.) then
            write(*,'(" TALYS-error: -10. <= etable <= 10.")')
            stop
          endif
          if (ftable(Zix,Nix).lt.0.1.or.ftable(Zix,Nix).gt.10.) then
            write(*,'(" TALYS-error: 0.1 <= ftable <= 10.")')
            stop
          endif
  220   continue
  210 continue
      if ((gnorm.le.0..or.gnorm.gt.1000.).and.gnorm.ne.-1.) then
        write(*,'(" TALYS-error: 0. < gnorm <= 100.")') 
        stop
      endif
c
c 6. Check of values for pre-equilibrium 
c
c epreeq    : on-set incident energy for preequilibrium calculation
c preeqmode : designator for pre-equilibrium model
c mpreeqmode: designator for multiple pre-equilibrium model 
c emulpre   : on-set incident energy for multiple preequilibrium
c pairmodel : model for preequilibrium pairing energy
c M2constant: constant for matrix element in exciton model
c M2limit   : constant for asymptotical value for matrix element
c M2shift   : constant for energy shift for matrix element      
c Rpinu,....: ratio for two-component matrix element
c Esurf0    : well depth for surface interaction
c Rgamma    : adjustable parameter for pre-equilibrium gamma decay  
c msdbins   : number of energy points for DWBA calculation for MSD
c numenmsd  : maximum number of energy points for DWBA calculation for 
c             MSD
c Emsdmin   : minimal outgoing energy for MSD calculation
c elwidth   : width of elastic peak in MeV
c flagpecomp: flag for Kalbach complex particle emission model 
c Cstrip    : adjustable parameter for stripping/pick-up reactions
c Cknock    : adjustable parameter for knockout reactions
c
      if ((epreeq.lt.0..or.epreeq.gt.250.).and.epreeq.ne.-1.) then
        write(*,'(" TALYS-error: 0. <= epreeq < 250.")') 
        stop
      endif
      if (preeqmode.lt.1.or.preeqmode.gt.4) then
        write(*,'(" TALYS-error: 1 <= preeqmode <= 4")')
        stop
      endif
      if (preeqmode.gt.3.and.ptype0.eq.'g') then
        write(*,'(" TALYS-error: preeqmode <= 3 for incident photons")')
        stop
      endif
      if (mpreeqmode.lt.1.or.mpreeqmode.gt.2) then
        write(*,'(" TALYS-error: 1 <= mpreeqmode <= 2")')
        stop
      endif
      if (emulpre.lt.0..or.emulpre.gt.250.) then
        write(*,'(" TALYS-error: 0. <= emulpre < 250.")') 
        stop
      endif
      if (pairmodel.lt.1.or.pairmodel.gt.2) then
        write(*,'(" TALYS-error: 1 <= pairmodel <= 2")')
        stop
      endif
      if (M2constant.lt.0..or.M2constant.gt.100.) then
        write(*,'(" TALYS-error: 0. <= M2constant <= 100.")') 
        stop
      endif
      if (M2limit.lt.0..or.M2limit.gt.100.) then
        write(*,'(" TALYS-error: 0. <= M2limit <= 100.")') 
        stop
      endif
      if (M2shift.lt.0..or.M2shift.gt.100.) then
        write(*,'(" TALYS-error: 0. <= M2shift <= 100.")') 
        stop
      endif
      if (Rpipi.lt.0..or.Rpipi.gt.100.) then
        write(*,'(" TALYS-error: 0. <= Rpipi <= 100.")') 
        stop
      endif
      if (Rnunu.lt.0..or.Rnunu.gt.100.) then
        write(*,'(" TALYS-error: 0. <= Rnunu <= 100.")') 
        stop
      endif
      if (Rpinu.lt.0..or.Rpinu.gt.100.) then
        write(*,'(" TALYS-error: 0. <= Rpinu <= 100.")') 
        stop
      endif
      if (Rnupi.lt.0..or.Rnupi.gt.100.) then
        write(*,'(" TALYS-error: 0. <= Rnupi <= 100.")') 
        stop
      endif
      if (Esurf0.ne.-1..and.(Esurf0.lt.0..or.Esurf0.gt.38.)) then
        write(*,'(" TALYS-error: 0. <= Esurf <= 38.")') 
        stop
      endif
      if (Rgamma.lt.0..or.Rgamma.gt.100.) then
        write(*,'(" TALYS-error: 0. <= Rgamma <= 100.")') 
        stop
      endif
      if ((msdbins.lt.2.or.msdbins.gt.numenmsd/2-1).and.msdbins.ne.0) 
     +  then
        write(*,'(" TALYS-error: 2 <= msdbins <=",i3)') numenmsd/2-1
        stop
      endif
      if (Emsdmin.lt.0.) then
        write(*,'(" TALYS-error: 0. <= E-in")') 
        stop
      endif
      if (elwidth.lt.1.e-6.or.elwidth.gt.100.) then
        write(*,'(" TALYS-error: 1.e-6 <= elwidth <= 100.")') 
        stop
      endif
      if (k0.eq.0.and.flagpecomp) then
        write(*,'(" TALYS-error: No pick-up and knock-out ",
     +    "mechanism for photonuclear reactions")')
        stop
      endif
      do 260 type=0,6
        if (Cstrip(type).lt.0..or.Cstrip(type).gt.100.) then
          write(*,'(" TALYS-error: 0. <= Cstrip <= 100.")')
          stop
        endif
        if (Cknock(type).lt.0..or.Cknock(type).gt.100.) then
          write(*,'(" TALYS-error: 0. <= Cknock <= 100.")')
          stop
        endif
  260 continue
c
c 7. Check of values for level densities 
c
c ldmodel     : level density model
c spincutmodel: model for spin cutoff factor for ground state
c shellmodel  : model for shell correction energies
c kvibmodel   : model for vibrational enhancement
c alev        : level density parameter
c alimit      : asymptotic level density parameter
c gammald     : gamma-constant for asymptotic level density parameter
c ibar        : fission barrier
c numbar      : number of fission barriers
c deltaW      : shell correction in nuclear mass
c Nlow        : lowest discrete level for temperature matching
c Ntop        : highest discrete level for temperature matching
c E0          : constant of temperature formula
c beta2       : deformation parameter
c Krotconstant: normalization constant for rotational enhancement
c T           : nuclear temperature
c Exmatch     : matching point for Ex 
c Pshift      : adjustable pairing shift
c pair        : pairing energy
c ctable      : constant to adjust tabulated level densities
c ptable      : constant to adjust tabulated level densities
c cglobal     : global constant to adjust tabulated level densities
c pglobal     : global constant to adjust tabulated level densities
c g           : single-particle level density parameter
c gp          : single-particle proton level density parameter
c gn          : single-particle neutron level density parameter   
c
      if (ldmodel.lt.1.or.ldmodel.gt.5) then
        write(*,'(" TALYS-error: 1 <= ldmodel <= 5")')
        stop
      endif
      if (ldmodel.ge.4) then
        inquire (file=path(1:lenpath)//
     +    'density/ground/hilaire/z050.tab',exist=lexist)
        if (.not.lexist) then
          write(*,'(" TALYS-error: Microscopic HFB tables are not ",
     +    " installed: download the full TALYS package from ",
     +    "www.talys.eu")')
          stop
        endif
      endif
      if (spincutmodel.lt.1.or.spincutmodel.gt.2) then
        write(*,'(" TALYS-error: 1 <= spincutmodel <= 2")')
        stop
      endif
      if (shellmodel.lt.1.or.shellmodel.gt.2) then
        write(*,'(" TALYS-error: 1 <= shellmodel <= 2")')
        stop
      endif
      if (kvibmodel.lt.1.or.kvibmodel.gt.2) then
        write(*,'(" TALYS-error: 1 <= kvibmodel <= 2")')
        stop
      endif
      do 310 Zix=0,numZ
        do 320 Nix=0,numN
          if (alev(Zix,Nix).ne.0.) then
            if (alev(Zix,Nix).lt.1.or.alev(Zix,Nix).gt.100.) then
              write(*,'(" TALYS-error: 1. <= a <= 100.")')
              stop
            endif
          endif
          if (alimit(Zix,Nix).ne.0.) then
            if (alimit(Zix,Nix).lt.1.or.alimit(Zix,Nix).gt.100.) then
              write(*,'(" TALYS-error: 1. <= alimit <= 100.")')
              stop
            endif
          endif
          if (gammald(Zix,Nix).ne.-1.) then
            if (gammald(Zix,Nix).lt.0..or.gammald(Zix,Nix).gt.1.) then
              write(*,'(" TALYS-error: 0. <= gammald <= 1.")')
              stop
            endif
          endif
          do 330 ibar=0,numbar
            if (deltaW(Zix,Nix,ibar).ne.0.) then
              if (deltaW(Zix,Nix,ibar).lt.-20..or.
     +          deltaW(Zix,Nix,ibar).gt.20.) then
                write(*,'(" TALYS-error: -20. <= deltaW <= 20.")')
                stop
              endif
            endif
            if (Nlow(Zix,Nix,ibar).ne.-1) then
              if (Nlow(Zix,Nix,ibar).lt.0.or.
     +          Nlow(Zix,Nix,ibar).gt.200) then
                write(*,'(" TALYS-error: 0 <= Nlow <= 200")')
                stop
              endif
            endif
            if (Ntop(Zix,Nix,ibar).ne.-1) then
              if (Ntop(Zix,Nix,ibar).lt.0.or.
     +          Ntop(Zix,Nix,ibar).gt.200) then
                write(*,'(" TALYS-error: 0 <= Ntop <= 200")')
                stop
              endif
            endif
            if (Nlow(Zix,Nix,ibar).ne.-1.and.Ntop(Zix,Nix,ibar).ne.-1.
     +        and.Nlow(Zix,Nix,ibar).ge.Ntop(Zix,Nix,ibar)) then
              write(*,'(" TALYS-error: Ntop <= Nlow")')
              stop
            endif
            if (E0(Zix,Nix,ibar).ne.0.) then
              if (E0(Zix,Nix,ibar).lt.-10..or.
     +          E0(Zix,Nix,ibar).gt.10.) then
                write(*,'(" TALYS-error: -10. <= E0 <= 10.")')
                stop
              endif
            endif
            if (beta2(Zix,Nix,ibar).lt.-0.5.or.
     +        beta2(Zix,Nix,ibar).ge.1.5) then
              write(*,'(" TALYS-error: -0.5 <= beta2 < 1.5")')
              stop
            endif
            if (Krotconstant(Zix,Nix,ibar).lt.0.01.or.
     +        Krotconstant(Zix,Nix,ibar).gt.100.) then
              write(*,'(" TALYS-error: 0.01 <= Krotconstant <= 100.")')
              stop
            endif
            if (T(Zix,Nix,ibar).ne.0.) then
              if (T(Zix,Nix,ibar).lt.1.e-3.or.
     +          T(Zix,Nix,ibar).gt.10.) then
                write(*,'(" TALYS-error: 1.e-3 <= T <= 10.")')
                stop
              endif
            endif
            if (Exmatch(Zix,Nix,ibar).ne.0.) then
              if (Exmatch(Zix,Nix,ibar).lt.0.1.or.
     +          Exmatch(Zix,Nix,ibar).gt.20.) then
                write(*,'(" TALYS-error: 0.1 <= Exmatch <= 20.")')
                stop
              endif
            endif
            if (Pshift(Zix,Nix,ibar).lt.-10..or.
     +        Pshift(Zix,Nix,ibar).gt.10.) then
              write(*,'(" TALYS-error: -10. <= Pshift <= 10.")')
              stop
            endif
            if (ctable(Zix,Nix,ibar).ne.0.) then
              if (ctable(Zix,Nix,ibar).lt.-10..or.ctable(Zix,Nix,ibar).
     +          gt.10.) then
                write(*,'(" TALYS-error: -10. <= ctable <= 10.")')
                stop
              endif
            endif
            if (ptable(Zix,Nix,ibar).ne.0.) then
              if (ptable(Zix,Nix,ibar).lt.-10..or.ptable(Zix,Nix,ibar).
     +          gt.10.) then
                write(*,'(" TALYS-error: -10. <= ptable <= 10.")')
                stop
              endif
            endif
  330     continue
          if (pair(Zix,Nix).lt.-10..or.pair(Zix,Nix).gt.10.) then
            write(*,'(" TALYS-error: -10. <= pair <= 10.")')
            stop
          endif
          if (cglobal.ne.0.) then
            if (cglobal.lt.-10..or.cglobal.gt.10.) then
              write(*,'(" TALYS-error: -10. <= cglobal <= 10.")')
              stop
            endif
          endif
          if (pglobal.ne.0.) then
            if (pglobal.lt.-10..or.pglobal.gt.10.) then
              write(*,'(" TALYS-error: -10. <= pglobal <= 10.")')
              stop
            endif
          endif
          if (g(Zix,Nix).ne.0.) then
            if (g(Zix,Nix).lt.0.1.or.g(Zix,Nix).gt.100.) then
              write(*,'(" TALYS-error: 0.1 <= g <= 100.")')
              stop
            endif
          endif
          if (gp(Zix,Nix).ne.0.) then
            if (gp(Zix,Nix).lt.0.1.or.gp(Zix,Nix).gt.100.) then
              write(*,'(" TALYS-error: 0.1 <= gp <= 100.")')
              stop
            endif
          endif
          if (gn(Zix,Nix).ne.0.) then
            if (gn(Zix,Nix).lt.0.1.or.gn(Zix,Nix).gt.100.) then
              write(*,'(" TALYS-error: 0.1 <= gn <= 100.")')
              stop
            endif
          endif
  320   continue
  310 continue
c
c There are many input possibilities for the energy dependent level
c density parameter of the Ignatyuk formula. The required parameters
c are alev, alimit, gammald and deltaW. The Ignatyuk formula implies
c that they can not all be given at the same time in the input file.
c
c alphad        : alpha-constant for asymptotic level density parameter
c betald        : beta-constant for asymptotic level density parameter
c gammashell1   : gamma-constant for asymptotic level density parameter
c gammashell2   : gamma-constant for asymptotic level density parameter
c pairconstant  : constant for pairing energy systematics
c Pshiftconstant: global constant for pairing shift
c Ufermi        : energy of Fermi distribution for damping of 
c                 ground-state rotational effects
c cfermi        : width of Fermi distribution for damping of 
c               : ground-state rotational effects
c Ufermibf      : energy of Fermi distribution for damping of barrier
c               : rotational effects
c cfermibf      : width of Fermi distribution for damping of barrier
c               : rotational effects           
c Kph           : constant for single-particle level density parameter
c                 (g=A/Kph)          
c Rspincut      : adjustable constant for spin cutoff factor
c
      do 340 Zix=0,numZ
        do 350 Nix=0,numN
          do 360 ibar=0,numbar
            if (alev(Zix,Nix).ne.0.and.deltaW(Zix,Nix,ibar).ne.0..and.
     +        alimit(Zix,Nix).ne.0.and.gammald(Zix,Nix).ne.-1.) then
              write(*,'(" TALYS-error: Level density conflict - a,",
     +          " deltaW, alimit and gammald are ALL given",
     +          " in the input for Z=",i3," A=",i3,
     +          " fission barrier=",i3)') Zinit-Zix,Ainit-Zix-Nix,ibar
              stop
            endif
  360     continue
  350   continue
  340 continue
      if (alphald.lt.0.01.or.alphald.gt.0.2) then
        write(*,'(" TALYS-error: 0.01 <= alphald <= 0.2")') 
        stop
      endif
      if (betald.lt.-0.5.or.betald.gt.0.5) then
        write(*,'(" TALYS-error: -0.5 <= betald <= 0.5")') 
        stop
      endif
      if (betald.lt.0..and.abs(betald).gt.alphald) then
        write(*,'(" TALYS-error: if betald < 0, |betald| < alphald")') 
        stop
      endif
      if (gammashell1.lt.0..or.gammashell1.gt.1.) then
        write(*,'(" TALYS-error: 0. <= gammashell1 <= 1.")') 
        stop
      endif
      if (gammashell2.lt.0..or.gammashell2.gt.0.2) then
        write(*,'(" TALYS-error: 0. <= gammashell2 <= 0.2")') 
        stop
      endif
      if (pairconstant.lt.0..or.pairconstant.gt.30.) then
        write(*,'(" TALYS-error: 0. <= pairconstant <= 30.")') 
        stop
      endif
      if (Pshiftconstant.lt.-5..or.Pshiftconstant.gt.5.) then
        write(*,'(" TALYS-error: -5. <= Pshiftconstant <= 5.")') 
        stop
      endif
      if (Ufermi.lt.0..or.Ufermi.gt.1000.) then
        write(*,'(" TALYS-error: 0. <= Ufermi <= 1000.")') 
        stop
      endif
      if (cfermi.le.0..or.cfermi.gt.1000.) then
        write(*,'(" TALYS-error: 0. <= cfermi <= 1000.")') 
        stop
      endif
      if (Ufermibf.lt.0..or.Ufermibf.gt.1000.) then
        write(*,'(" TALYS-error: 0. <= Ufermibf <= 1000.")') 
        stop
      endif
      if (cfermibf.le.0..or.cfermibf.gt.1000.) then
        write(*,'(" TALYS-error: 0. <= cfermibf <= 1000.")') 
        stop
      endif
      if (Kph.lt.1..or.Kph.gt.100.) then
        write(*,'(" TALYS-error: 1. <= Kph <= 100.")') 
        stop
      endif
      if (Rspincut.le.0..or.Rspincut.gt.10.) then
        write(*,'(" TALYS-error: 0. < Rspincut <= 10.")') 
        stop
      endif
c
c 8. Check of values for fission
c
c flagfission: flag for fission
c flagfisout : flag for output of fission information 
c fismodel   : fission model  
c fismodelalt: alternative fission model for default barriers
c axtype     : type of axiality of barrier 
c                 1: axial symmetry
c                 2: left-right asymmetry
c                 3: triaxial and left-right symmetry
c                 4: triaxial no left-right symmetry
c                 5: no symmetry
c fbarrier   : height of fission barrier
c fwidth     : width of fission barrier   
c Rtransmom  : normalization constant for moment of inertia for 
c              transition states
c Rclass2mom : normalization constant for moment of inertia for 
c              class 2 states
c widthc2    : width of class2 states
c betafiscor : adjustable factor for fission path width
c vfiscor    : adjustable factor for fission path height
c
      if ((flagfission.or.flagfisout).and.Atarget.le.56) then
        write(*,'(" TALYS-error: Fission not allowed for A <= 56")')
        stop
      endif
      if (fismodel.lt.1.or.fismodel.gt.5) then
        write(*,'(" TALYS-error: 1 <= fismodel <= 5")')
        stop
      endif
      if (fismodelalt.lt.3.or.fismodelalt.gt.4) then
        write(*,'(" TALYS-error: 3 <= fismodelalt <= 4")')
        stop
      endif
      do 410 Zix=0,numZ
        do 420 Nix=0,numN
          do 430 ibar=1,numbar
            fax0=axtype(Zix,Nix,ibar)
            if (fax0.ne.0.and.(fax0.lt.1.or.fax0.gt.5)) then
              write(*,'(" TALYS-error: 1 <= type of axiality <= 5")')
              stop
            endif
            fbar0=fbarrier(Zix,Nix,ibar)
            if (fbar0.ne.0..and.(fbar0.lt.0..or.fbar0.gt.100.)) then
              write(*,'(" TALYS-error: 0.<= fission barrrier <=100.")')
              stop
            endif
            fhw0=fwidth(Zix,Nix,ibar)
            if (fhw0.ne.0..and.(fhw0.lt.0.01.or.fhw0.gt.10.)) then
              write(*,'(" TALYS-error: 0.01 <= fission width <= 10.")')
              stop
            endif
            fR0=Rtransmom(Zix,Nix,ibar)
            if (fR0.ne.0..and.(fR0.lt.0.1.or.fR0.gt.10.)) then
              write(*,'(" TALYS-error: 0.1 <= Rtransmom <= 10.")')
              stop
            endif
            fR0=Rclass2mom(Zix,Nix,ibar)
            if (fR0.ne.0..and.(fR0.lt.0.1.or.fR0.gt.10.)) then
              write(*,'(" TALYS-error: 0.1 <= Rclass2mom <= 10.")')
              stop
            endif
  430     continue
          fhw0=betafiscor(Zix,Nix)
          if (fhw0.lt.0.1.or.fhw0.gt.10.) then
            write(*,'(" TALYS-error: 0.1 <= betafiscor <=10.")')
            stop
          endif
          fbar0=vfiscor(Zix,Nix)
          if (fbar0.lt.0.1.or.fbar0.gt.10.) then
            write(*,'(" TALYS-error: 0.1 <= vfiscor <=10.")')
            stop
          endif
  420   continue
  410 continue
c
c 9. Check of values for output
c
c ddxmode  : mode for double-differential cross sections: 0: None,
c            1: Angular distributions, 2: Spectra per angle, 3: Both
c ddxecount: counter for double-differential cross section files
c fileddxe : designator for double-differential cross sections on
c            separate file: angular distribution
c enincmax : maximum incident energy
c ddxacount: counter for double-differential cross section files
c fileddxa : designator for double-differential cross sections on
c            separate file: spectrum per angle    
c
      if (ddxmode.lt.0.or.ddxmode.gt.3) then
        write(*,'(" TALYS-error: 0 <= ddxmode <= 3")')
        stop
      endif
      do 510 type=0,6
        do 520 i=1,ddxecount(type)
          value=fileddxe(type,i)
          if (value.lt.0.or.value.gt.enincmax) then
            write(*,'(" TALYS-error: 0. <= fileddxe <=",f7.3)') enincmax
            stop
          endif
  520   continue
        do 530 i=1,ddxacount(type)
          value=fileddxa(type,i)
          if (value.lt.0.or.value.gt.180.) then
            write(*,'(" TALYS-error: 0. <= fileddxa <= 180.")')
            stop
          endif
  530   continue
  510 continue
      return
      end
Copyright (C) 2004  A.J. Koning, S. Hilaire and M.C. Duijvestijn
