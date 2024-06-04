      subroutine input5
c
c +---------------------------------------------------------------------
c | Author: Arjan Koning 
c | Date  : December 7, 2009
c | Task  : Read input for fifth set of variables
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      character*1  ch
      character*80 word(40),key,value
      integer      Zix,Nix,ibar,irad,Z,N,oddZ,oddN,lval,type,mt,igr,i,
     +             iz,ia,ivalue,type2,omptype,nr
      real         val
c
c ************** Defaults for fifth set of input variables *************
c
c k0            : index for incident particle
c flagendf      : flag for information for ENDF-6 file
c eninclow      : minimal incident energy for nuclear model calculations
c Atarget       : mass number of target nucleus
c strength      : model for E1 gamma-ray strength function
c gnorm         : gamma normalization factor
c Rspincut      : adjustable constant (global) for spin cutoff factor
c spincutmodel  : model for spin cutoff factor for ground state
c shellmodel    : model for shell correction energies
c kvibmodel     : model for vibrational enhancement
c gammashell2   : gamma-constant for asymptotic level density parameter 
c flagcol       : flag for collective enhancement of level density
c flagcolldamp  : flag for damping of collective effects in effective
c                 level density (without explicit collective 
c                 enhancement). Only used for Bruyeres-le-Chatel 
c                 (Pascal Romain) fission model
c alphad        : alpha-constant for asymptotic level density parameter
c betald        : beta-constant for asymptotic level density parameter
c gammashell1   : gamma-constant for asymptotic level density parameter 
c pairconstant  : constant for pairing energy systematics
c ldmodel       : level density model
c Pshiftconstant: global constant for pairing shift
c Ufermi        : energy of Fermi distribution for damping of 
c                 ground-state rotational effects
c cfermi        : width of Fermi distribution for damping of 
c                 ground-state rotational effects
c Ufermibf      : energy of Fermi distribution for damping of barrier
c                 rotational effects
c cfermibf      : width of Fermi distribution for damping of barrier
c                 rotational effects              
c Kph           : constant for single-particle level density parameter 
c                 (g=A/Kph)
c M2constant    : overall constant for matrix element in exciton model
c M2limit       : constant for asymptotical value for matrix element
c M2shift       : constant for energy shift for matrix element
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
c alev          : level density parameter
c alimit        : asymptotic level density parameter
c gammald       : gamma-constant for asymptotic level density parameter
c pair          : total pairing correction
c ibar          : fission barrier
c numbar        : number of fission barriers
c Pshift        : adjustable pairing shift
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
c cglobal       : global constant to adjust tabulated level densities
c pglobal       : global constant to adjust tabulated level densities
c g             : single-particle level density parameter
c gp            : single-particle proton level density parameter
c gn            : single-particle neutron level density parameter
c gamgam        : experimental total radiative width in eV
c D0            : experimental s-wave resonance spacing in eV
c S0            : s-wave strength function
c etable,ftable : constant to adjust tabulated strength functions
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
c fiso          : correction factor for isospin forbidden transitions
c aadjust....   : adjustable factors for level density parameters 
c                 (default 1.)
c axtype        : type of axiality of barrier 
c                 1: axial symmetry
c                 2: left-right asymmetry
c                 3: triaxial and left-right symmetry
c                 4: triaxial no left-right symmetry
c                 5: no symmetry
c fbarrier      : height of fission barrier
c fwidth        : width of fission barrier 
c betafiscor    : adjustable factor for fission path width
c Zinit         : charge number of initial compound nucleus
c Ninit         : neutron number of initial compound nucleus
c vfiscor       : adjustable factor for fission path height
c fismodel      : fission model
c Rtransmom     : normalization constant for moment of inertia for 
c                 transition states
c Rclass2mom    : normalization constant for moment of inertia for 
c                 class 2 states
c widthc2       : width of class2 states
c Ninit         : neutron number of initial compound nucleus
c hbtransfile   : file with head band transition states
c class2file    : file with class 2 transition states
c levelfile     : discrete level file   
c deformfile    : deformation parameter file   
c optmodfileN   : optical model parameter file for neutrons
c optmodfileP   : optical model parameter file for protons
c radialfile    : radial matter density file
c ompenergyfile : file with energies for OMP calculation (ENDF files 
c                 only)
c radialmodel   : model for radial matter densities (JLM OMP only)
c massnucleus   : mass of nucleus in amu as read from user input file
c massexcess    : mass excess in MeV as read from user input file
c beta2         : deformation parameter
c msdbins       : number of energy points for DWBA calculation for MSD
c Emsdmin       : minimal outgoing energy for MSD calculation
c Cstrip        : adjustable parameter for stripping/pick-up reactions
c Cknock        : adjustable parameter for knockout reactions
c v1adjust....  : adjustable factors for OMP (default 1.)
c ompadjustF    : logical for local OMP adjustment
c ompadjustN    : number of energy ranges for local OMP adjustment
c ompadjustE1   : start energy of local OMP adjustment
c ompadjustE2   : end energy of local OMP adjustment
c ompadjustD    : depth of local OMP adjustment
c ompadjusts    : variance of local OMP adjustment
c jlmmode       : option for JLM imaginary potential normalization
c flagrescue    : flag for final rescue: normalization to data
c rescuefile    : file with incident energy dependent adjustment factors
c grescue       : global multiplication factor for incident energy 
c                 dependent adjustment factors
c flaglabddx    : flag for calculation of DDX in LAB system
c nanglerec     : number of recoil angles
c
      if (k0.eq.1.and.flagendf) then
        eninclow=0.
      else
        eninclow=1.e-6
      endif
c
c Advice of S. Goriely: no gamma normalization for A < 40.
c
      if (k0.ne.1.or.Atarget.lt.40.or.strength.eq.3.or.strength.eq.4) 
     +  then
        gnorm=1.
      else
        gnorm=-1.
      endif
      Rspincut=1.
      spincutmodel=1
      shellmodel=1
      kvibmodel=2
      gammashell2=0.
c
c Level density systematics 
c
c ldmodel 1: Gilbert and Cameron
c ldmodel 2: Back-shifted Fermi gas
c ldmodel 3: Superfluid model
c
      if (ldmodel.eq.1.or.ldmodel.ge.4) then
        if (flagcol) then
          alphald=0.0207305
          betald=0.229537
          gammashell1=0.473625
          Pshiftconstant=0.
        else
          alphald=0.0692559
          betald=0.282769
          gammashell1=0.433090
          Pshiftconstant=0.
        endif                              
        if (flagcolldamp) then
          alphald=0.0666
          betald=0.258
          gammashell1=0.459
          Pshiftconstant=0.
        endif
      endif                              
      if (ldmodel.eq.2) then
        if (flagcol) then
          alphald=0.0381563
          betald=0.105378
          gammashell1=0.546474
          Pshiftconstant=0.743229
        else
          alphald=0.0722396
          betald=0.195267
          gammashell1=0.410289
          Pshiftconstant=0.173015
        endif                              
      endif                              
      if (ldmodel.eq.3) then
        if (flagcol) then
          alphald=0.0357750
          betald=0.135307
          gammashell1=0.699663
          Pshiftconstant=-0.149106
        else
          alphald=0.110575
          betald=0.0313662
          gammashell1=0.648723
          Pshiftconstant=1.13208
        endif                              
      endif                              
      pairconstant=12.
      Ufermi=30.
      cfermi=5.
      Ufermibf=45.
      cfermibf=5.
      Kph=15.
      M2constant=1.
      M2limit=1.
      M2shift=1.
      Rpipi=1.
      Rnunu=1.5
      Rpinu=1.
      Rnupi=1.
      Esurf0=-1.
      Rgamma=2.
      elwidth=0.5
      xscaptherm=0.
      xsptherm=0.
      xsalphatherm=0.
      do 10 Zix=0,numZ
        do 20 Nix=0,numN
          alev(Zix,Nix)=0.
          alimit(Zix,Nix)=0.
          gammald(Zix,Nix)=-1.
          pair(Zix,Nix)=1.e-20
          do 30 ibar=0,numbar
            Pshift(Zix,Nix,ibar)=1.e-20
            deltaW(Zix,Nix,ibar)=0.
            Exmatch(Zix,Nix,ibar)=0.
            T(Zix,Nix,ibar)=0.
            E0(Zix,Nix,ibar)=1.e-20
            Nlow(Zix,Nix,ibar)=-1
            Ntop(Zix,Nix,ibar)=-1
            s2adjust(Zix,Nix,ibar)=1.
            Krotconstant(Zix,Nix,ibar)=1.
            ctable(Zix,Nix,ibar)=cglobal
            ptable(Zix,Nix,ibar)=pglobal
   30     continue     
          g(Zix,Nix)=0.
          gp(Zix,Nix)=0.
          gn(Zix,Nix)=0.
          gamgam(Zix,Nix)=0.
          D0(Zix,Nix)=0.
          S0(Zix,Nix)=0.
          etable(Zix,Nix)=0.
          ftable(Zix,Nix)=1.
          do 40 irad=0,1
            do 40 lval=1,numgam
              do 40 igr=1,2
                egr(Zix,Nix,irad,lval,igr)=0.
                ggr(Zix,Nix,irad,lval,igr)=0.
                sgr(Zix,Nix,irad,lval,igr)=0.
                epr(Zix,Nix,irad,lval)=0.
                gpr(Zix,Nix,irad,lval)=0.
                tpr(Zix,Nix,irad,lval)=0.
   40     continue     
          do type=0,6
            fiso(Zix,Nix,type)=1.
          enddo
          if (Zix.eq.0.and.Nix.eq.0) then
            if (Zinit.eq.Ninit) then
              if (k0.eq.1) fiso(Zix,Nix,k0)=2.
              if (k0.eq.2) fiso(Zix,Nix,k0)=2.
              if (k0.eq.6) fiso(Zix,Nix,k0)=5.
            endif
            if (Zinit.eq.Ninit-1.or.Zinit.eq.Ninit+1) then
              if (k0.eq.1) fiso(Zix,Nix,k0)=1.5
              if (k0.eq.2) fiso(Zix,Nix,k0)=1.5
              if (k0.eq.6) fiso(Zix,Nix,k0)=1.5
            endif
          endif
          aadjust(Zix,Nix)=1.
          gnadjust(Zix,Nix)=1.
          gpadjust(Zix,Nix)=1.
          gamgamadjust(Zix,Nix)=1.
          do 50 ibar=1,numbar
            axtype(Zix,Nix,ibar)=1
            fbarrier(Zix,Nix,ibar)=0.
            fwidth(Zix,Nix,ibar)=0.
            Rtransmom(Zix,Nix,ibar)=1.
            Rclass2mom(Zix,Nix,ibar)=1.
            widthc2(Zix,Nix,ibar)=0.2
   50     continue     
          betafiscor(Zix,Nix)=1.
          Z=Zinit-Zix
          N=Ninit-Nix
          oddZ=mod(Z,2)
          oddN=mod(N,2)
          vfiscor(Zix,Nix)=1.
          if (oddZ.eq.0.and.oddN.eq.0) vfiscor(Zix,Nix)=0.86
          if (oddZ.eq.1.and.oddN.eq.0) vfiscor(Zix,Nix)=0.94
          if (oddZ.eq.0.and.oddN.eq.1) vfiscor(Zix,Nix)=0.89
          if (oddZ.eq.1.and.oddN.eq.1) vfiscor(Zix,Nix)=1.02
          fismodelx(Zix,Nix)=fismodel
          if (Ninit-Nix.gt.144.or.fismodel.eq.5) axtype(Zix,Nix,1)=3
          if (fismodel.ne.5) axtype(Zix,Nix,2)=2
          Rtransmom(Zix,Nix,1)=0.6
          hbtransfile(Zix,Nix)='                                       '
          class2file(Zix,Nix)='                                        '
   20   continue     
        levelfile(Zix)='                                               '
        deformfile(Zix)='                                              '
        optmodfileN(Zix)='                                            '
        optmodfileP(Zix)='                                            '
        radialfile(Zix)='                                              '
   10 continue     
      ompenergyfile='                                            '
      radialmodel=2
      do 60 Zix=0,numZ+4
        do 60 Nix=0,numN+4
          massnucleus(Zix,Nix)=0.
          massexcess(Zix,Nix)=0.
          beta2(Zix,Nix,0)=0.
          beta2(Zix,Nix,1)=0.6
          beta2(Zix,Nix,2)=0.8
          beta2(Zix,Nix,3)=1. 
   60 continue     
      msdbins=6
      Emsdmin=0.
      do 70 type=0,6
        Cstrip(type)=1.
        Cknock(type)=1.
   70 continue     
      do 80 type=1,6
        v1adjust(type)=1.
        v2adjust(type)=1.
        v3adjust(type)=1.
        v4adjust(type)=1.
        rvadjust(type)=1.
        avadjust(type)=1.
        w1adjust(type)=1.
        w2adjust(type)=1.
        d1adjust(type)=1.
        d2adjust(type)=1.
        d3adjust(type)=1.
        rvdadjust(type)=1.
        avdadjust(type)=1.
        vso1adjust(type)=1.
        vso2adjust(type)=1.
        rvsoadjust(type)=1.
        avsoadjust(type)=1.
        wso1adjust(type)=1.
        wso2adjust(type)=1.
        rcadjust(type)=1.
        do 90 omptype=1,6
          ompadjustN(type,omptype)=0
          do 92 nr=1,numrange
            ompadjustE1(type,omptype,nr)=0.
            ompadjustE2(type,omptype,nr)=0.
            ompadjustD(type,omptype,nr)=1.
            ompadjusts(type,omptype,nr)=1.
   92     continue     
   90   continue     
   80 continue     
      do 94 omptype=1,6
        ompadjustF(omptype)=.false.
   94 continue     
      flagrescue=.false.
      do 100 mt=1,nummt
        rescuefile(mt)='                                             '
        grescue(mt)=1.
  100 continue     
      jlmmode=0
      lvadjust=1.
      lwadjust=1.
      lv1adjust=1.
      lw1adjust=1.
      lvsoadjust=1.
      lwsoadjust=1.
      if (flaglabddx) then
        nanglerec=numangrec
      else
        nanglerec=1
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
      do 110 i=1,nlines
        call getkeywords(inline(i),word)
        key=word(1)
        value=word(2)
        ch=word(2)(1:1)
c
c Test for keywords
c
c Here, the various model parameters can be set to overrule the default
c values. Most default values will be computed later on, since they
c require more computation (e.g. level density parameters).
c
c iz    : charge number
c ia    : mass number
c val   : help variable 
c ivalue: help variable 
c
        if (key.eq.'elow') then
          read(value,*,end=1000,err=1000) eninclow
          goto 110
        endif
        if (key.eq.'massnucleus') then
          read(word(2),*,end=1000,err=1000) iz
          read(word(3),*,end=1000,err=1000) ia
          read(word(4),*,end=1000,err=1000) val
          Zix=Zinit-iz   
          Nix=Ninit-ia+iz
          if (Zix.lt.0.or.Zix.gt.numZ.or.Nix.lt.0.or.Nix.gt.numN) then
            goto 1040
          else
            massnucleus(Zix,Nix)=val
          endif
          goto 110
        endif
        if (key.eq.'massexcess') then
          read(word(2),*,end=1000,err=1000) iz
          read(word(3),*,end=1000,err=1000) ia
          read(word(4),*,end=1000,err=1000) val
          Zix=Zinit-iz   
          Nix=Ninit-ia+iz
          if (Zix.lt.0.or.Zix.gt.numZ.or.Nix.lt.0.or.Nix.gt.numN) then
            goto 1040
          else
            massexcess(Zix,Nix)=val
          endif
          goto 110
        endif
        if (key.eq.'a') then
          read(word(2),*,end=1000,err=1000) iz
          read(word(3),*,end=1000,err=1000) ia
          read(word(4),*,end=1000,err=1000) val
          Zix=Zinit-iz   
          Nix=Ninit-ia+iz
          if (Zix.lt.0.or.Zix.gt.numZ.or.Nix.lt.0.or.Nix.gt.numN) then
            goto 1040
          else
            alev(Zix,Nix)=val
          endif
          goto 110
        endif
        if (key.eq.'alimit') then
          read(word(2),*,end=1000,err=1000) iz
          read(word(3),*,end=1000,err=1000) ia
          read(word(4),*,end=1000,err=1000) val
          Zix=Zinit-iz   
          Nix=Ninit-ia+iz
          if (Zix.lt.0.or.Zix.gt.numZ.or.Nix.lt.0.or.Nix.gt.numN) then
            goto 1040
          else
            alimit(Zix,Nix)=val
          endif
          goto 110
        endif
        if (key.eq.'gammald') then
          read(word(2),*,end=1000,err=1000) iz
          read(word(3),*,end=1000,err=1000) ia
          read(word(4),*,end=1000,err=1000) val
          Zix=Zinit-iz   
          Nix=Ninit-ia+iz
          if (Zix.lt.0.or.Zix.gt.numZ.or.Nix.lt.0.or.Nix.gt.numN) then
            goto 1040
          else
            gammald(Zix,Nix)=val
          endif
          goto 110
        endif
        if (key.eq.'pair') then
          read(word(2),*,end=1000,err=1000) iz
          read(word(3),*,end=1000,err=1000) ia
          read(word(4),*,end=1000,err=1000) val
          Zix=Zinit-iz   
          Nix=Ninit-ia+iz
          if (Zix.lt.0.or.Zix.gt.numZ.or.Nix.lt.0.or.Nix.gt.numN) then
            goto 1040
          else
            pair(Zix,Nix)=val
          endif
          goto 110
        endif
        if (key.eq.'pshift') then
          ibar=0
          read(word(2),*,end=1000,err=1000) iz
          read(word(3),*,end=1000,err=1000) ia
          read(word(4),*,end=1000,err=1000) val
          Zix=Zinit-iz   
          Nix=Ninit-ia+iz
          if (Zix.lt.0.or.Zix.gt.numZ.or.Nix.lt.0.or.Nix.gt.numN) then
            goto 1040
          else
            read(word(5),*,end=115,err=1000) ibar
            if (ibar.lt.0.or.ibar.gt.numbar) goto 1010
  115       Pshift(Zix,Nix,ibar)=val
          endif
          goto 110
        endif
        if (key.eq.'deltaw') then
          ibar=0
          read(word(2),*,end=1000,err=1000) iz
          read(word(3),*,end=1000,err=1000) ia
          read(word(4),*,end=1000,err=1000) val
          Zix=Zinit-iz   
          Nix=Ninit-ia+iz
          if (Zix.lt.0.or.Zix.gt.numZ.or.Nix.lt.0.or.Nix.gt.numN) then
            goto 1040
          else
            read(word(5),*,end=120,err=1000) ibar
            if (ibar.lt.0.or.ibar.gt.numbar) goto 1010
  120       deltaW(Zix,Nix,ibar)=val
          endif
          goto 110
        endif
        if (key.eq.'exmatch') then
          ibar=0
          read(word(2),*,end=1000,err=1000) iz
          read(word(3),*,end=1000,err=1000) ia
          read(word(4),*,end=1000,err=1000) val
          Zix=Zinit-iz   
          Nix=Ninit-ia+iz
          if (Zix.lt.0.or.Zix.gt.numZ.or.Nix.lt.0.or.Nix.gt.numN) then
            goto 1040
          else
            read(word(5),*,end=130,err=1000) ibar
            if (ibar.lt.0.or.ibar.gt.numbar) goto 1010
  130       Exmatch(Zix,Nix,ibar)=val
          endif
          goto 110
        endif
        if (key.eq.'t') then
          ibar=0
          read(word(2),*,end=1000,err=1000) iz
          read(word(3),*,end=1000,err=1000) ia
          read(word(4),*,end=1000,err=1000) val
          Zix=Zinit-iz   
          Nix=Ninit-ia+iz
          if (Zix.lt.0.or.Zix.gt.numZ.or.Nix.lt.0.or.Nix.gt.numN) then
            goto 1040
          else
            read(word(5),*,end=140,err=1000) ibar
            if (ibar.lt.0.or.ibar.gt.numbar) goto 1010
  140       T(Zix,Nix,ibar)=val
          endif
          goto 110
        endif
        if (key.eq.'e0') then
          ibar=0
          read(word(2),*,end=1000,err=1000) iz
          read(word(3),*,end=1000,err=1000) ia
          read(word(4),*,end=1000,err=1000) val
          Zix=Zinit-iz   
          Nix=Ninit-ia+iz
          if (Zix.lt.0.or.Zix.gt.numZ.or.Nix.lt.0.or.Nix.gt.numN) then
            goto 1040
          else
            read(word(5),*,end=150,err=1000) ibar
            if (ibar.lt.0.or.ibar.gt.numbar) goto 1010
  150       E0(Zix,Nix,ibar)=val
          endif
          goto 110
        endif
        if (key.eq.'nlow') then
          ibar=0
          read(word(2),*,end=1000,err=1000) iz
          read(word(3),*,end=1000,err=1000) ia
          read(word(4),*,end=1000,err=1000) ivalue
          Zix=Zinit-iz   
          Nix=Ninit-ia+iz
          if (Zix.lt.0.or.Zix.gt.numZ.or.Nix.lt.0.or.Nix.gt.numN) then
            goto 1040
          else
            read(word(5),*,end=160,err=1000) ibar
            if (ibar.lt.0.or.ibar.gt.numbar) goto 1010
  160       Nlow(Zix,Nix,ibar)=ivalue
          endif
          goto 110
        endif
        if (key.eq.'ntop') then
          ibar=0
          read(word(2),*,end=1000,err=1000) iz
          read(word(3),*,end=1000,err=1000) ia
          read(word(4),*,end=1000,err=1000) ivalue
          Zix=Zinit-iz   
          Nix=Ninit-ia+iz
          if (Zix.lt.0.or.Zix.gt.numZ.or.Nix.lt.0.or.Nix.gt.numN) then
            goto 1040
          else
            read(word(5),*,end=170,err=1000) ibar
            if (ibar.lt.0.or.ibar.gt.numbar) goto 1010
  170       Ntop(Zix,Nix,ibar)=ivalue
          endif
          goto 110
        endif
        if (key.eq.'s2adjust') then
          ibar=0
          read(word(2),*,end=1000,err=1000) iz
          read(word(3),*,end=1000,err=1000) ia
          read(word(4),*,end=1000,err=1000) val
          Zix=Zinit-iz   
          Nix=Ninit-ia+iz
          if (Zix.lt.0.or.Zix.gt.numZ.or.Nix.lt.0.or.Nix.gt.numN) then
            goto 1040
          else
            read(word(5),*,end=180,err=1000) ibar
            if (ibar.lt.0.or.ibar.gt.numbar) goto 1010
  180       s2adjust(Zix,Nix,ibar)=val
          endif
          goto 110
        endif
        if (key.eq.'krotconstant') then
          ibar=0
          read(word(2),*,end=1000,err=1000) iz
          read(word(3),*,end=1000,err=1000) ia
          read(word(4),*,end=1000,err=1000) val
          Zix=Zinit-iz   
          Nix=Ninit-ia+iz
          if (Zix.lt.0.or.Zix.gt.numZ.or.Nix.lt.0.or.Nix.gt.numN) then
            goto 1040
          else
            read(word(5),*,end=182,err=1000) ibar
            if (ibar.lt.0.or.ibar.gt.numbar) goto 1010
  182       Krotconstant(Zix,Nix,ibar)=val
          endif
          goto 110
        endif
        if (key.eq.'ctable') then
          ibar=0
          read(word(2),*,end=1000,err=1000) iz
          read(word(3),*,end=1000,err=1000) ia
          read(word(4),*,end=1000,err=1000) val
          Zix=Zinit-iz   
          Nix=Ninit-ia+iz
          if (Zix.lt.0.or.Zix.gt.numZ.or.Nix.lt.0.or.Nix.gt.numN) then
            goto 1040
          else
            read(word(5),*,end=185,err=1000) ibar
            if (ibar.lt.0.or.ibar.gt.numbar) goto 1010
  185       ctable(Zix,Nix,ibar)=val
          endif
          goto 110
        endif
        if (key.eq.'ptable') then
          ibar=0
          read(word(2),*,end=1000,err=1000) iz
          read(word(3),*,end=1000,err=1000) ia
          read(word(4),*,end=1000,err=1000) val
          Zix=Zinit-iz   
          Nix=Ninit-ia+iz
          if (Zix.lt.0.or.Zix.gt.numZ.or.Nix.lt.0.or.Nix.gt.numN) then
            goto 1040
          else
            read(word(5),*,end=187,err=1000) ibar
            if (ibar.lt.0.or.ibar.gt.numbar) goto 1010
  187       ptable(Zix,Nix,ibar)=val
          endif
          goto 110
        endif
        if (key.eq.'g') then
          read(word(2),*,end=1000,err=1000) iz
          read(word(3),*,end=1000,err=1000) ia
          read(word(4),*,end=1000,err=1000) val
          Zix=Zinit-iz   
          Nix=Ninit-ia+iz
          if (Zix.lt.0.or.Zix.gt.numZ.or.Nix.lt.0.or.Nix.gt.numN) then
            goto 1040
          else
            g(Zix,Nix)=val
          endif
          goto 110
        endif
        if (key.eq.'gp') then
          read(word(2),*,end=1000,err=1000) iz
          read(word(3),*,end=1000,err=1000) ia
          read(word(4),*,end=1000,err=1000) val
          Zix=Zinit-iz   
          Nix=Ninit-ia+iz
          if (Zix.lt.0.or.Zix.gt.numZ.or.Nix.lt.0.or.Nix.gt.numN) then
            goto 1040
          else
            gp(Zix,Nix)=val
          endif
          goto 110
        endif
        if (key.eq.'gn') then
          read(word(2),*,end=1000,err=1000) iz
          read(word(3),*,end=1000,err=1000) ia
          read(word(4),*,end=1000,err=1000) val
          Zix=Zinit-iz   
          Nix=Ninit-ia+iz
          if (Zix.lt.0.or.Zix.gt.numZ.or.Nix.lt.0.or.Nix.gt.numN) then
            goto 1040
          else
            gn(Zix,Nix)=val
          endif
          goto 110
        endif
        if (key.eq.'egr') then
          igr=1
          read(word(2),*,end=1000,err=1000) iz
          read(word(3),*,end=1000,err=1000) ia
          read(word(4),*,end=1000,err=1000) val
          Zix=Zinit-iz   
          Nix=Ninit-ia+iz
          if (Zix.lt.0.or.Zix.gt.numZ.or.Nix.lt.0.or.Nix.gt.numN) then
            goto 1040
          else
            ch=word(5)(1:1)
            if (ch.eq.'m'.or.ch.eq.'e') then
              if (ch.eq.'m') irad=0
              if (ch.eq.'e') irad=1
              read(word(5)(2:2),*,end=1000,err=1000) lval
              if (lval.lt.1.or.lval.gt.numgam) goto 1030
              read(word(6),*,end=190,err=1000) igr
              if (igr.lt.1.or.igr.gt.2) goto 1020
  190         egr(Zix,Nix,irad,lval,igr)=val
              goto 110
            endif
            goto 1000
          endif
        endif
        if (key.eq.'ggr') then
          igr=1
          read(word(2),*,end=1000,err=1000) iz
          read(word(3),*,end=1000,err=1000) ia
          read(word(4),*,end=1000,err=1000) val
          Zix=Zinit-iz   
          Nix=Ninit-ia+iz
          if (Zix.lt.0.or.Zix.gt.numZ.or.Nix.lt.0.or.Nix.gt.numN) then
            goto 1040
          else
            ch=word(5)(1:1)
            if (ch.eq.'m'.or.ch.eq.'e') then
              if (ch.eq.'m') irad=0
              if (ch.eq.'e') irad=1
              read(word(5)(2:2),*,end=1000,err=1000) lval
              if (lval.lt.1.or.lval.gt.numgam) goto 1030
              read(word(6),*,end=200,err=1000) igr
              if (igr.lt.1.or.igr.gt.2) goto 1020
  200         ggr(Zix,Nix,irad,lval,igr)=val
              goto 110
            endif
            goto 1000
          endif
        endif
        if (key.eq.'sgr') then
          igr=1
          read(word(2),*,end=1000,err=1000) iz
          read(word(3),*,end=1000,err=1000) ia
          read(word(4),*,end=1000,err=1000) val
          Zix=Zinit-iz   
          Nix=Ninit-ia+iz
          if (Zix.lt.0.or.Zix.gt.numZ.or.Nix.lt.0.or.Nix.gt.numN) then
            goto 1040
          else
            ch=word(5)(1:1)
            if (ch.eq.'m'.or.ch.eq.'e') then
              if (ch.eq.'m') irad=0
              if (ch.eq.'e') irad=1
              read(word(5)(2:2),*,end=1000,err=1000) lval
              if (lval.lt.1.or.lval.gt.numgam) goto 1030
              read(word(6),*,end=210,err=1000) igr
              if (igr.lt.1.or.igr.gt.2) goto 1020
  210         sgr(Zix,Nix,irad,lval,igr)=val
              goto 110
            endif
            goto 1000
          endif
        endif
        if (key.eq.'epr') then
          read(word(2),*,end=1000,err=1000) iz
          read(word(3),*,end=1000,err=1000) ia
          read(word(4),*,end=1000,err=1000) val
          Zix=Zinit-iz
          Nix=Ninit-ia+iz
          if (Zix.lt.0.or.Zix.gt.numZ.or.Nix.lt.0.or.Nix.gt.numN) then
            goto 1040
          else
            ch=word(5)(1:1)
            if (ch.eq.'m'.or.ch.eq.'e') then
              if (ch.eq.'m') irad=0
              if (ch.eq.'e') irad=1
              read(word(5)(2:2),*,end=1000,err=1000) lval
              if (lval.lt.1.or.lval.gt.numgam) goto 1030
              epr(Zix,Nix,irad,lval)=val
              goto 110
            endif
            goto 1000
          endif
        endif
        if (key.eq.'gpr') then
          read(word(2),*,end=1000,err=1000) iz
          read(word(3),*,end=1000,err=1000) ia
          read(word(4),*,end=1000,err=1000) val
          Zix=Zinit-iz
          Nix=Ninit-ia+iz
          if (Zix.lt.0.or.Zix.gt.numZ.or.Nix.lt.0.or.Nix.gt.numN) then
            goto 1040
          else
            ch=word(5)(1:1)
            if (ch.eq.'m'.or.ch.eq.'e') then
              if (ch.eq.'m') irad=0
              if (ch.eq.'e') irad=1
              read(word(5)(2:2),*,end=1000,err=1000) lval
              if (lval.lt.1.or.lval.gt.numgam) goto 1030
              gpr(Zix,Nix,irad,lval)=val
              goto 110
            endif
            goto 1000
          endif
        endif
        if (key.eq.'spr') then
          read(word(2),*,end=1000,err=1000) iz
          read(word(3),*,end=1000,err=1000) ia
          read(word(4),*,end=1000,err=1000) val
          Zix=Zinit-iz   
          Nix=Ninit-ia+iz
          if (Zix.lt.0.or.Zix.gt.numZ.or.Nix.lt.0.or.Nix.gt.numN) then
            goto 1040
          else
            ch=word(5)(1:1)
            if (ch.eq.'m'.or.ch.eq.'e') then
              if (ch.eq.'m') irad=0
              if (ch.eq.'e') irad=1
              read(word(5)(2:2),*,end=1000,err=1000) lval
              if (lval.lt.1.or.lval.gt.numgam) goto 1030
              tpr(Zix,Nix,irad,lval)=val
              goto 110
            endif
            goto 1000
          endif
        endif
        if (key.eq.'fiso') then
          do 215 type=0,6
            if (ch.eq.parsym(type)) then
              type2=type
              goto 216
            endif
  215     continue
          goto 1000
  216     continue
          read(word(3),*,end=1000,err=1000) iz
          read(word(4),*,end=1000,err=1000) ia
          read(word(5),*,end=1000,err=1000) val
          Zix=Zinit-iz   
          Nix=Ninit-ia+iz
          if (Zix.lt.0.or.Zix.gt.numZ.or.Nix.lt.0.or.Nix.gt.numN) then
            goto 1040
          else
            fiso(Zix,Nix,type2)=val
          endif
          goto 110
        endif
        if (key.eq.'gamgam') then
          read(word(2),*,end=1000,err=1000) iz
          read(word(3),*,end=1000,err=1000) ia
          read(word(4),*,end=1000,err=1000) val
          Zix=Zinit-iz   
          Nix=Ninit-ia+iz
          if (Zix.lt.0.or.Zix.gt.numZ.or.Nix.lt.0.or.Nix.gt.numN) then
            goto 1040
          else
            gamgam(Zix,Nix)=val
          endif
          goto 110
        endif
        if (key.eq.'d0') then
          read(word(2),*,end=1000,err=1000) iz
          read(word(3),*,end=1000,err=1000) ia
          read(word(4),*,end=1000,err=1000) val
          Zix=Zinit-iz   
          Nix=Ninit-ia+iz
          if (Zix.lt.0.or.Zix.gt.numZ.or.Nix.lt.0.or.Nix.gt.numN) then
            goto 1040
          else
            D0(Zix,Nix)=val*1000.
          endif
          goto 110
        endif
        if (key.eq.'aadjust') then
          read(word(2),*,end=1000,err=1000) iz
          read(word(3),*,end=1000,err=1000) ia
          read(word(4),*,end=1000,err=1000) val
          Zix=Zinit-iz   
          Nix=Ninit-ia+iz
          if (Zix.lt.0.or.Zix.gt.numZ.or.Nix.lt.0.or.Nix.gt.numN) then
            goto 1040
          else
            aadjust(Zix,Nix)=val
          endif
          goto 110
        endif
        if (key.eq.'gnadjust') then
          read(word(2),*,end=1000,err=1000) iz
          read(word(3),*,end=1000,err=1000) ia
          read(word(4),*,end=1000,err=1000) val
          Zix=Zinit-iz   
          Nix=Ninit-ia+iz
          if (Zix.lt.0.or.Zix.gt.numZ.or.Nix.lt.0.or.Nix.gt.numN) then
            goto 1040
          else
            gnadjust(Zix,Nix)=val
          endif
          goto 110
        endif
        if (key.eq.'gpadjust') then
          read(word(2),*,end=1000,err=1000) iz
          read(word(3),*,end=1000,err=1000) ia
          read(word(4),*,end=1000,err=1000) val
          Zix=Zinit-iz   
          Nix=Ninit-ia+iz
          if (Zix.lt.0.or.Zix.gt.numZ.or.Nix.lt.0.or.Nix.gt.numN) then
            goto 1040
          else
            gpadjust(Zix,Nix)=val
          endif
          goto 110
        endif
        if (key.eq.'gamgamadjust') then
          read(word(2),*,end=1000,err=1000) iz
          read(word(3),*,end=1000,err=1000) ia
          read(word(4),*,end=1000,err=1000) val
          Zix=Zinit-iz   
          Nix=Ninit-ia+iz
          if (Zix.lt.0.or.Zix.gt.numZ.or.Nix.lt.0.or.Nix.gt.numN) then
            goto 1040
          else
            gamgamadjust(Zix,Nix)=val
          endif
          goto 110
        endif
        if (key.eq.'s0') then
          read(word(2),*,end=1000,err=1000) iz
          read(word(3),*,end=1000,err=1000) ia
          read(word(4),*,end=1000,err=1000) val
          Zix=Zinit-iz   
          Nix=Ninit-ia+iz
          if (Zix.lt.0.or.Zix.gt.numZ.or.Nix.lt.0.or.Nix.gt.numN) then
            goto 1040
          else
            S0(Zix,Nix)=val
          endif
          goto 110
        endif
        if (key.eq.'etable') then
          read(word(2),*,end=1000,err=1000) iz
          read(word(3),*,end=1000,err=1000) ia
          read(word(4),*,end=1000,err=1000) val
          Zix=Zinit-iz   
          Nix=Ninit-ia+iz
          if (Zix.lt.0.or.Zix.gt.numZ.or.Nix.lt.0.or.Nix.gt.numN) then
            goto 1040
          else
            etable(Zix,Nix)=val
          endif
          goto 110
        endif
        if (key.eq.'ftable') then
          read(word(2),*,end=1000,err=1000) iz
          read(word(3),*,end=1000,err=1000) ia
          read(word(4),*,end=1000,err=1000) val
          Zix=Zinit-iz   
          Nix=Ninit-ia+iz
          if (Zix.lt.0.or.Zix.gt.numZ.or.Nix.lt.0.or.Nix.gt.numN) then
            goto 1040
          else
            ftable(Zix,Nix)=val
          endif
          goto 110
        endif
        if (key.eq.'fisbar') then
          ibar=1
          read(word(2),*,end=1000,err=1000) iz
          read(word(3),*,end=1000,err=1000) ia
          read(word(4),*,end=1000,err=1000) val
          Zix=Zinit-iz
          Nix=Ninit-ia+iz
          if (Zix.lt.0.or.Zix.gt.numZ.or.Nix.lt.0.or.Nix.gt.numN) then
            goto 1040
          else
            read(word(5),*,end=220,err=1000) ibar
            if (ibar.lt.1.or.ibar.gt.numbar) goto 1010
  220       fbarrier(Zix,Nix,ibar)=val
          endif
          goto 110
        endif                   
        if (key.eq.'fishw') then
          ibar=1
          read(word(2),*,end=1000,err=1000) iz
          read(word(3),*,end=1000,err=1000) ia
          read(word(4),*,end=1000,err=1000) val
          Zix=Zinit-iz
          Nix=Ninit-ia+iz
          if (Zix.lt.0.or.Zix.gt.numZ.or.Nix.lt.0.or.Nix.gt.numN) then
            goto 1040
          else
            read(word(5),*,end=230,err=1000) ibar
            if (ibar.lt.1.or.ibar.gt.numbar) goto 1010
  230       fwidth(Zix,Nix,ibar)=val
          endif
          goto 110
        endif                   
        if (key.eq.'rtransmom') then
          ibar=1
          read(word(2),*,end=1000,err=1000) iz
          read(word(3),*,end=1000,err=1000) ia
          read(word(4),*,end=1000,err=1000) val
          Zix=Zinit-iz
          Nix=Ninit-ia+iz
          if (Zix.lt.0.or.Zix.gt.numZ.or.Nix.lt.0.or.Nix.gt.numN) then
            goto 1040
          else
            read(word(5),*,end=240,err=1000) ibar
            if (ibar.lt.1.or.ibar.gt.numbar) goto 1010
  240       Rtransmom(Zix,Nix,ibar)=val
          endif
          goto 110
        endif                   
        if (key.eq.'rclass2mom') then
          ibar=1
          read(word(2),*,end=1000,err=1000) iz
          read(word(3),*,end=1000,err=1000) ia
          read(word(4),*,end=1000,err=1000) val
          Zix=Zinit-iz
          Nix=Ninit-ia+iz
          if (Zix.lt.0.or.Zix.gt.numZ.or.Nix.lt.0.or.Nix.gt.numN) then
            goto 1040
          else
            read(word(5),*,end=250,err=1000) ibar
            if (ibar.lt.1.or.ibar.gt.numbar) goto 1010
  250       Rclass2mom(Zix,Nix,ibar)=val
          endif
          goto 110
        endif                   
        if (key.eq.'class2width') then
          ibar=1
          read(word(2),*,end=1000,err=1000) iz
          read(word(3),*,end=1000,err=1000) ia
          read(word(4),*,end=1000,err=1000) val
          Zix=Zinit-iz
          Nix=Ninit-ia+iz
          if (Zix.lt.0.or.Zix.gt.numZ.or.Nix.lt.0.or.Nix.gt.numN) then
            goto 1040
          else
            read(word(5),*,end=260,err=1000) ibar
            if (ibar.lt.1.or.ibar.gt.numbar) goto 1010
  260       widthc2(Zix,Nix,ibar)=val
          endif
          goto 110
        endif                   
        if (key.eq.'levelfile') then
          read(word(2),*,end=1000,err=1000) iz
          Zix=Zinit-iz
          if (Zix.lt.0.or.Zix.gt.numZ) then
            goto 1040
          else
            levelfile(Zix)=word(3)
            goto 110
          endif
        endif                   
        if (key.eq.'deformfile') then
          read(word(2),*,end=1000,err=1000) iz
          Zix=Zinit-iz
          if (Zix.lt.0.or.Zix.gt.numZ) then
            goto 1040
          else
            deformfile(Zix)=word(3)
            goto 110
          endif
        endif                   
        if (key.eq.'hbtransfile') then
          read(word(2),*,end=1000,err=1000) iz
          read(word(3),*,end=1000,err=1000) ia
          Zix=Zinit-iz
          Nix=Ninit-ia+iz
          if (Zix.lt.0.or.Zix.gt.numZ.or.Nix.lt.0.or.Nix.gt.numN) then
            goto 1040
          else
            hbtransfile(Zix,Nix)=word(4)
            goto 110
          endif
        endif                   
        if (key.eq.'class2file') then
          read(word(2),*,end=1000,err=1000) iz
          read(word(3),*,end=1000,err=1000) ia
          Zix=Zinit-iz
          Nix=Ninit-ia+iz
          if (Zix.lt.0.or.Zix.gt.numZ.or.Nix.lt.0.or.Nix.gt.numN) then
            goto 1040
          else
            class2file(Zix,Nix)=word(4)
            goto 110
          endif
        endif                   
        if (key.eq.'optmodfilen') then
          read(word(2),*,end=1000,err=1000) iz
          Zix=Zinit-iz
          if (Zix.lt.0.or.Zix.gt.numZ) then
            goto 1040
          else
            optmodfileN(Zix)=word(3)
            goto 110
          endif
        endif                   
        if (key.eq.'optmodfilep') then
          read(word(2),*,end=1000,err=1000) iz
          Zix=Zinit-iz
          if (Zix.lt.0.or.Zix.gt.numZ) then
            goto 1040
          else
            optmodfileP(Zix)=word(3)
            goto 110
          endif
        endif                   
        if (key.eq.'radialfile') then
          read(word(2),*,end=1000,err=1000) iz
          Zix=Zinit-iz
          if (Zix.lt.0.or.Zix.gt.numZ) then
            goto 1040
          else
            radialfile(Zix)=word(3)
            goto 110
          endif
        endif                   
        if (key.eq.'ompenergyfile') then
          ompenergyfile=value
          goto 110
        endif                   
        if (key.eq.'beta2') then
          ibar=0
          read(word(2),*,end=1000,err=1000) iz
          read(word(3),*,end=1000,err=1000) ia
          read(word(4),*,end=1000,err=1000) val
          Zix=Zinit-iz   
          Nix=Ninit-ia+iz
          if (Zix.lt.0.or.Zix.gt.numZ.or.Nix.lt.0.or.Nix.gt.numN) then
            goto 1040
          else
            read(word(5),*,end=270,err=1000) ibar
            if (ibar.lt.0.or.ibar.gt.numbar) goto 1010
  270       beta2(Zix,Nix,ibar)=val
          endif
          goto 110
        endif
        if (key.eq.'axtype') then
          ibar=1
          read(word(2),*,end=1000,err=1000) iz
          read(word(3),*,end=1000,err=1000) ia
          read(word(4),*,end=1000,err=1000) ivalue
          Zix=Zinit-iz   
          Nix=Ninit-ia+iz
          if (Zix.lt.0.or.Zix.gt.numZ.or.Nix.lt.0.or.Nix.gt.numN) then
            goto 1040
          else
            read(word(5),*,end=280,err=1000) ibar
            if (ibar.lt.1.or.ibar.gt.numbar) goto 1010
  280       axtype(Zix,Nix,ibar)=ivalue
          endif
          goto 110
        endif
        if (key.eq.'vfiscor') then
          read(word(2),*,end=1000,err=1000) iz
          read(word(3),*,end=1000,err=1000) ia
          read(word(4),*,end=1000,err=1000) val
          Zix=Zinit-iz   
          Nix=Ninit-ia+iz
          if (Zix.lt.0.or.Zix.gt.numZ.or.Nix.lt.0.or.Nix.gt.numN) then
            goto 1040
          else
            vfiscor(Zix,Nix)=val
          endif
          goto 110
        endif
        if (key.eq.'betafiscor') then
          read(word(2),*,end=1000,err=1000) iz
          read(word(3),*,end=1000,err=1000) ia
          read(word(4),*,end=1000,err=1000) val
          Zix=Zinit-iz   
          Nix=Ninit-ia+iz
          if (Zix.lt.0.or.Zix.gt.numZ.or.Nix.lt.0.or.Nix.gt.numN) then
            goto 1040
          else
            betafiscor(Zix,Nix)=val
          endif
          goto 110
        endif
        if (key.eq.'cstrip') then
          do 290 type=0,6
            if (ch.eq.parsym(type)) then
              type2=type
              goto 300
            endif
  290     continue
          goto 1000
  300     read(word(3),*,end=1000,err=1000) Cstrip(type2)
          goto 110
        endif
        if (key.eq.'cknock') then
          do 310 type=0,6
            if (ch.eq.parsym(type)) then
              type2=type
              goto 320
            endif
  310     continue
          goto 1000
  320     read(word(3),*,end=1000,err=1000) Cknock(type2)
          goto 110
        endif
        if (key.eq.'v1adjust') then
          do 330 type=1,6
            if (ch.eq.parsym(type)) then
              type2=type
              goto 340
            endif
  330     continue
          goto 1000
  340     read(word(3),*,end=1000,err=1000) v1adjust(type2)
          goto 110
        endif
        if (key.eq.'v2adjust') then
          do 350 type=1,6
            if (ch.eq.parsym(type)) then
              type2=type
              goto 360
            endif
  350     continue
          goto 1000
  360     read(word(3),*,end=1000,err=1000) v2adjust(type2)
          goto 110
        endif
        if (key.eq.'v3adjust') then
          do 370 type=1,6
            if (ch.eq.parsym(type)) then
              type2=type
              goto 380
            endif
  370     continue
          goto 1000
  380     read(word(3),*,end=1000,err=1000) v3adjust(type2)
          goto 110
        endif
        if (key.eq.'v4adjust') then
          do 390 type=1,6
            if (ch.eq.parsym(type)) then
              type2=type
              goto 400
            endif
  390     continue
          goto 1000
  400     read(word(3),*,end=1000,err=1000) v4adjust(type2)
          goto 110
        endif
        if (key.eq.'rvadjust') then
          do 410 type=1,6
            if (ch.eq.parsym(type)) then
              type2=type
              goto 420
            endif
  410     continue
          goto 1000
  420     read(word(3),*,end=1000,err=1000) rvadjust(type2)
          goto 110
        endif
        if (key.eq.'avadjust') then
          do 430 type=1,6
            if (ch.eq.parsym(type)) then
              type2=type
              goto 440
            endif
  430     continue
          goto 1000
  440     read(word(3),*,end=1000,err=1000) avadjust(type2)
          goto 110
        endif
        if (key.eq.'w1adjust') then
          do 450 type=1,6
            if (ch.eq.parsym(type)) then
              type2=type
              goto 460
            endif
  450     continue
          goto 1000
  460     read(word(3),*,end=1000,err=1000) w1adjust(type2)
          goto 110
        endif
        if (key.eq.'w2adjust') then
          do 470 type=1,6
            if (ch.eq.parsym(type)) then
              type2=type
              goto 480
            endif
  470     continue
          goto 1000
  480     read(word(3),*,end=1000,err=1000) w2adjust(type2)
          goto 110
        endif
        if (key.eq.'d1adjust') then
          do 490 type=1,6
            if (ch.eq.parsym(type)) then
              type2=type
              goto 500
            endif
  490     continue
          goto 1000
  500     read(word(3),*,end=1000,err=1000) d1adjust(type2)
          goto 110
        endif
        if (key.eq.'d2adjust') then
          do 510 type=1,6
            if (ch.eq.parsym(type)) then
              type2=type
              goto 520
            endif
  510     continue
          goto 1000
  520     read(word(3),*,end=1000,err=1000) d2adjust(type2)
          goto 110
        endif
        if (key.eq.'d3adjust') then
          do 530 type=1,6
            if (ch.eq.parsym(type)) then
              type2=type
              goto 540
            endif
  530     continue
          goto 1000
  540     read(word(3),*,end=1000,err=1000) d3adjust(type2)
          goto 110
        endif
        if (key.eq.'rvsoadjust') then
          do 550 type=1,6
            if (ch.eq.parsym(type)) then
              type2=type
              goto 560
            endif
  550     continue
          goto 1000
  560     read(word(3),*,end=1000,err=1000) rvsoadjust(type2)
          goto 110
        endif
        if (key.eq.'avsoadjust') then
          do 570 type=1,6
            if (ch.eq.parsym(type)) then
              type2=type
              goto 580
            endif
  570     continue
          goto 1000
  580     read(word(3),*,end=1000,err=1000) avsoadjust(type2)
          goto 110
        endif
        if (key.eq.'rvdadjust') then
          do 590 type=1,6
            if (ch.eq.parsym(type)) then
              type2=type
              goto  600
            endif
  590     continue
          goto 1000
  600     read(word(3),*,end=1000,err=1000) rvdadjust(type2)
          goto 110
        endif
        if (key.eq.'avdadjust') then
          do 610 type=1,6
            if (ch.eq.parsym(type)) then
              type2=type
              goto 620
            endif
  610     continue
          goto 1000
  620     read(word(3),*,end=1000,err=1000) avdadjust(type2)
          goto 110
        endif
        if (key.eq.'vso1adjust') then
          do 630 type=1,6
            if (ch.eq.parsym(type)) then
              type2=type
              goto 640
            endif
  630     continue
          goto 1000
  640     read(word(3),*,end=1000,err=1000) vso1adjust(type2)
          goto 110
        endif
        if (key.eq.'vso2adjust') then
          do 650 type=1,6
            if (ch.eq.parsym(type)) then
              type2=type
              goto 660
            endif
  650     continue
          goto 1000
  660     read(word(3),*,end=1000,err=1000) vso2adjust(type2)
          goto 110
        endif
        if (key.eq.'wso1adjust') then
          do 670 type=1,6
            if (ch.eq.parsym(type)) then
              type2=type
              goto 680
            endif
  670     continue
          goto 1000
  680     read(word(3),*,end=1000,err=1000) wso1adjust(type2)
          goto 110
        endif
        if (key.eq.'wso2adjust') then
          do 690 type=1,6
            if (ch.eq.parsym(type)) then
              type2=type
              goto 700
            endif
  690     continue
          goto 1000
  700     read(word(3),*,end=1000,err=1000) wso2adjust(type2)
          goto 110
        endif
        if (key.eq.'rcadjust') then
          do 710 type=1,6
            if (ch.eq.parsym(type)) then
              type2=type
              goto 720
            endif
  710     continue
          goto 1000
  720     read(word(3),*,end=1000,err=1000) rcadjust(type2)
          goto 110
        endif
        if (key.eq.'rvadjustf'.or.key.eq.'avadjustf'.or.
     +    key.eq.'rwdadjustf'.or.key.eq.'awdadjustf'.or.
     +    key.eq.'rvsoadjustf'.or.key.eq.'avsoadjustf') then
          if (key.eq.'rvadjustf') omptype=1
          if (key.eq.'avadjustf') omptype=2
          if (key.eq.'rwdadjustf') omptype=3
          if (key.eq.'awdadjustf') omptype=4
          if (key.eq.'rvsoadjustf') omptype=5
          if (key.eq.'avsoadjustf') omptype=6
          do 730 type=1,6
            if (ch.eq.parsym(type)) then
              type2=type
              goto 740
            endif
  730     continue
          goto 1000
  740     ompadjustF(type2)=.true.
          ompadjustN(type2,omptype)=ompadjustN(type2,omptype)+1
          nr=ompadjustN(type2,omptype)
          read(word(3),*,end=1000,err=1000) 
     +      ompadjustE1(type2,omptype,nr)
          read(word(4),*,end=1000,err=1000) 
     +      ompadjustE2(type2,omptype,nr)
          read(word(5),*,end=1000,err=1000) ompadjustD(type2,omptype,nr)
          read(word(6),*,end=110,err=1000) ompadjusts(type2,omptype,nr)
          goto 110
        endif
        if (key.eq.'rescuefile') then
          read(word(2),*,end=1000,err=1000) mt
          rescuefile(mt)=word(3)
          val=1.
          read(word(4),*,end=750,err=1000) val
  750     grescue(mt)=val
          flagrescue=.true.
          goto 110
        endif
        if (key.eq.'jlmmode') then
          read(value,*,end=1000,err=1000) jlmmode
          goto 110
        endif
        if (key.eq.'lvadjust') then
          read(value,*,end=1000,err=1000) lvadjust 
          goto 110
        endif
        if (key.eq.'lwadjust') then
          read(value,*,end=1000,err=1000) lwadjust 
          goto 110
        endif
        if (key.eq.'lv1adjust') then
          read(value,*,end=1000,err=1000) lv1adjust 
          goto 110
        endif
        if (key.eq.'lw1adjust') then
          read(value,*,end=1000,err=1000) lw1adjust 
          goto 110
        endif
        if (key.eq.'lvsoadjust') then
          read(value,*,end=1000,err=1000) lvsoadjust 
          goto 110
        endif
        if (key.eq.'lwsoadjust') then
          read(value,*,end=1000,err=1000) lwsoadjust 
          goto 110
        endif
        if (key.eq.'gnorm') then
          read(value,*,end=1000,err=1000) gnorm 
          goto 110
        endif
        if (key.eq.'spincutmodel') then
          read(value,*,end=1000,err=1000) spincutmodel 
          goto 110
        endif
        if (key.eq.'shellmodel') then
          read(value,*,end=1000,err=1000) shellmodel 
          goto 110
        endif
        if (key.eq.'kvibmodel') then
          read(value,*,end=1000,err=1000) kvibmodel 
          goto 110
        endif
        if (key.eq.'radialmodel') then
          read(value,*,end=1000,err=1000) radialmodel 
          goto 110
        endif
        if (key.eq.'rspincut') then
          read(value,*,end=1000,err=1000) Rspincut 
          goto 110
        endif
        if (key.eq.'alphald') then
          read(value,*,end=1000,err=1000) alphald 
          goto 110
        endif
        if (key.eq.'betald') then
          read(value,*,end=1000,err=1000) betald 
          goto 110
        endif
        if (key.eq.'gammashell1') then
          read(value,*,end=1000,err=1000) gammashell1
          goto 110
        endif
        if (key.eq.'gammashell2') then
          read(value,*,end=1000,err=1000) gammashell2
          goto 110
        endif
        if (key.eq.'pairconstant') then
          read(value,*,end=1000,err=1000) pairconstant
          goto 110
        endif
        if (key.eq.'pshiftconstant') then
          read(value,*,end=1000,err=1000) Pshiftconstant
          goto 110
        endif
        if (key.eq.'ufermi') then
          read(value,*,end=1000,err=1000) Ufermi 
          goto 110
        endif
        if (key.eq.'cfermi') then
          read(value,*,end=1000,err=1000) cfermi 
          goto 110
        endif
        if (key.eq.'ufermibf') then
          read(value,*,end=1000,err=1000) Ufermibf 
          goto 110
        endif
        if (key.eq.'cfermibf') then
          read(value,*,end=1000,err=1000) cfermibf 
          goto 110
        endif
        if (key.eq.'kph') then
          read(value,*,end=1000,err=1000) Kph 
          goto 110
        endif
        if (key.eq.'m2constant') then
          read(value,*,end=1000,err=1000) M2constant 
          goto 110
        endif
        if (key.eq.'m2limit') then
          read(value,*,end=1000,err=1000) M2limit 
          goto 110
        endif
        if (key.eq.'m2shift') then
          read(value,*,end=1000,err=1000) M2shift 
          goto 110
        endif
        if (key.eq.'rpipi') then
          read(value,*,end=1000,err=1000) Rpipi 
          goto 110
        endif
        if (key.eq.'rnunu') then
          read(value,*,end=1000,err=1000) Rnunu 
          goto 110
        endif
        if (key.eq.'rpinu') then
          read(value,*,end=1000,err=1000) Rpinu 
          goto 110
        endif
        if (key.eq.'rnupi') then
          read(value,*,end=1000,err=1000) Rnupi 
          goto 110
        endif
        if (key.eq.'esurf') then
          read(value,*,end=1000,err=1000) Esurf0
          goto 110
        endif
        if (key.eq.'rgamma') then
          read(value,*,end=1000,err=1000) Rgamma 
          goto 110
        endif
        if (key.eq.'msdbins')  then
          read(value,*,end=1000,err=1000) msdbins          
          goto 110
        endif
        if (key.eq.'emsdmin')  then
          read(value,*,end=1000,err=1000) Emsdmin          
          goto 110
        endif
        if (key.eq.'elwidth')  then
          read(value,*,end=1000,err=1000) elwidth          
          goto 110
        endif
        if (key.eq.'xscaptherm')  then
          read(value,*,end=1000,err=1000) xscaptherm          
          goto 110
        endif
        if (key.eq.'xsptherm')  then
          read(value,*,end=1000,err=1000) xsptherm          
          goto 110
        endif
        if (key.eq.'xsalphatherm')  then
          read(value,*,end=1000,err=1000) xsalphatherm          
          goto 110
        endif
        if (key.eq.'anglesrec')  then
          read(value,*,end=1000,err=1000) nanglerec
          goto 110
        endif
        goto 110
 1040   write(*,'(" TALYS-warning: Z,N index out of range,",
     +    " keyword ignored: ",a80)') inline(i)
  110 continue
      return
 1000 write(*,'(" TALYS-error: Wrong input: ",a80)') inline(i)
      stop
 1010 write(*,'(" TALYS-error: 0(1) <= fission barrier <=",i3,
     +  ", ibar index out of range: ",a80)') numbar,inline(i)
      stop
 1020 write(*,'(" TALYS-error: 0 <= resonance number <= 2",
     +  ", igr index out of range: ",a80)') inline(i)
      stop
 1030 write(*,'(" TALYS-error: 0 <= multipole radiation <= ",i1,
     +  ", lval index out of range: ",a80)') numgam,inline(i)
      stop
      end
Copyright (C) 2004  A.J. Koning, S. Hilaire and M.C. Duijvestijn
