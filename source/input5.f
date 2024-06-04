      subroutine input5
c
c +---------------------------------------------------------------------
c | Author: Arjan Koning 
c | Date  : December 8, 2006    
c | Task  : Read input for fifth set of variables
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      character*1  ch
      character*80 word(40),key,value
      integer      Zix,Nix,ibar,irad,lval,type,igr,i,iz,ia,ivalue,type2
      real         val
c
c ************** Defaults for fifth set of input variables *************
c
c k0            : index for incident particle
c flagendf      : flag for information for ENDF-6 file
c eninclow      : minimal incident energy for nuclear model calculations
c Atarget       : mass number of target nucleus
c gnorm         : gamma normalization factor
c Rspincut      : adjustable constant for spin cutoff factor
c spincutmodel  : model for spin cutoff factor for ground state
c shellmodel    : model for shell correction energies
c gammashell2   : gamma-constant for asymptotic level density parameter 
c flagcol       : flag for collective enhancement of level density
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
c Zix           : charge number index for residual nucleus
c numZ          : maximal number of protons away from the initial 
c                 compound nucleus
c Nix           : neutron number index for residual nucleus
c numN          : maximal number of neutrons away from the initial 
c                 compound nucleus
c massnucleus   : mass of nucleus in amu as read from user input file
c massexcess    : mass excess in MeV as read from user input file
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
c beta2         : deformation parameter
c Krotconstant  : normalization constant for rotational enhancement
c c1table       : constant to adjust tabulated level densities
c c2table       : constant to adjust tabulated level densities
c g             : single-particle level density parameter
c gp            : single-particle proton level density parameter
c gn            : single-particle neutron level density parameter
c gamgam        : experimental total radiative width in eV
c D0            : experimental s-wave resonance spacing in eV
c S0            : s-wave strength function
c irad          : variable to indicate M(=0) or E(=1) radiation
c lval          : multipolarity
c numgam        : maximum number of l-values for gamma multipolarity
c igr           : giant resonance
c egr           : energy of GR
c ggr           : width of GR
c sgr           : strength of GR
c axtype        : type of axiality of barrier (1: axial, 2: tri-axial)
c fbarrier      : height of fission barrier
c fwidth        : width of fission barrier 
c Rtransmom     : normalization constant for moment of inertia for 
c                 transition states
c Rclass2mom    : normalization constant for moment of inertia for 
c                 class 2 states
c widthc2       : width of class2 states
c Ntarget       : neutron number of target nucleus
c hbtransfile   : file with head band transition states
c class2file    : file with class 2 transition states
c levelfile     : discrete level file   
c deformfile    : deformation parameter file   
c optmodfileN   : optical model parameter file for neutrons
c optmodfileP   : optical model parameter file for protons
c msdbins       : number of energy points for DWBA calculation for MSD
c Emsdmin       : minimal outgoing energy for MSD calculation
c Cstrip        : adjustable parameter for stripping/pick-up reactions
c Cknock        : adjustable parameter for knockout reactions
c v1adjust....  : adjustable factors for OMP (default 1.)
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
      if (k0.ne.1.or.Atarget.lt.40) then
        gnorm=1.
      else
        gnorm=-1.
      endif
      Rspincut=1.
      spincutmodel=1
      shellmodel=1
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
          alphald=0.057706
          betald=0.081575
          gammashell1=0.57099
          Pshiftconstant=0.
        else
          alphald=0.071601
          betald=0.27159
          gammashell1=0.36510
          Pshiftconstant=0.
        endif                              
      endif                              
      if (ldmodel.eq.2) then
        if (flagcol) then
          alphald=0.066997
          betald=0.055159
          gammashell1=0.61274
          Pshiftconstant=1.2059 
        else
          alphald=0.076417
          betald=0.20856
          gammashell1=0.38829
          Pshiftconstant=0.51642
        endif                              
      endif                              
      if (ldmodel.eq.3) then
        if (flagcol) then
          alphald=0.072524
          betald=0.026177
          gammashell1=0.66344
          Pshiftconstant=-0.36779 
        else
          alphald=0.10073
          betald=0.092402
          gammashell1=0.65768
          Pshiftconstant=0.98198
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
      do 10 Zix=0,numZ
        do 20 Nix=0,numN
          massnucleus(Zix,Nix)=0.
          massexcess(Zix,Nix)=0.
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
            Krotconstant(Zix,Nix,ibar)=1.
   30     continue     
          beta2(Zix,Nix,0)=0.
          beta2(Zix,Nix,1)=0.6
          beta2(Zix,Nix,2)=0.8
          beta2(Zix,Nix,3)=1. 
          c1table(Zix,Nix)=0.
          c2table(Zix,Nix)=0.
          g(Zix,Nix)=0.
          gp(Zix,Nix)=0.
          gn(Zix,Nix)=0.
          gamgam(Zix,Nix)=0.
          D0(Zix,Nix)=0.
          S0(Zix,Nix)=0.
          do 40 irad=0,1
            do 40 lval=1,numgam
              do 40 igr=1,2
                egr(Zix,Nix,irad,lval,igr)=0.
                ggr(Zix,Nix,irad,lval,igr)=0.
                sgr(Zix,Nix,irad,lval,igr)=0.
   40     continue     
          do 50 ibar=1,numbar
            axtype(Zix,Nix,ibar)=1
            fbarrier(Zix,Nix,ibar)=0.
            fwidth(Zix,Nix,ibar)=0.
            Rtransmom(Zix,Nix,ibar)=1.
            Rclass2mom(Zix,Nix,ibar)=1.
            widthc2(Zix,Nix,ibar)=0.2
   50     continue     
          if (Ntarget.gt.144) axtype(Zix,Nix,1)=2
          Rtransmom(Zix,Nix,1)=0.6
          hbtransfile(Zix,Nix)='                                       '
          class2file(Zix,Nix)='                                        '
   20   continue     
        levelfile(Zix)='                                               '
        deformfile(Zix)='                                              '
        optmodfileN(Zix)='                                            '
        optmodfileP(Zix)='                                            '
   10 continue     
      msdbins=6
      Emsdmin=0.
      do 60 type=0,6
        Cstrip(type)=1.
        Cknock(type)=1.
   60 continue     
      do 70 type=1,6
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
   70 continue     
      if (flaglabddx) then
        nanglerec=9
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
c Zinit : charge number of initial compound nucleus
c Ninit : neutron number of initial compound nucleus
c ivalue: help variable 
c Ainit : mass number of initial compound nucleus
c
        if (key.eq.'elow') then
          read(value,*,err=800) eninclow
          goto 110
        endif
        if (key.eq.'massnucleus') then
          read(word(2),*,err=800) iz
          read(word(3),*,err=800) ia
          read(word(4),*,err=800) val
          Zix=Zinit-iz   
          Nix=Ninit-ia+iz
          if (Zix.lt.0.or.Zix.gt.numZ.or.Nix.lt.0.or.Nix.gt.numN) then
            goto 840
          else
            massnucleus(Zix,Nix)=val
          endif
          goto 110
        endif
        if (key.eq.'massexcess') then
          read(word(2),*,err=800) iz
          read(word(3),*,err=800) ia
          read(word(4),*,err=800) val
          Zix=Zinit-iz   
          Nix=Ninit-ia+iz
          if (Zix.lt.0.or.Zix.gt.numZ.or.Nix.lt.0.or.Nix.gt.numN) then
            goto 840
          else
            massexcess(Zix,Nix)=val
          endif
          goto 110
        endif
        if (key.eq.'a') then
          read(word(2),*,err=800) iz
          read(word(3),*,err=800) ia
          read(word(4),*,err=800) val
          Zix=Zinit-iz   
          Nix=Ninit-ia+iz
          if (Zix.lt.0.or.Zix.gt.numZ.or.Nix.lt.0.or.Nix.gt.numN) then
            goto 840
          else
            alev(Zix,Nix)=val
          endif
          goto 110
        endif
        if (key.eq.'alimit') then
          read(word(2),*,err=800) iz
          read(word(3),*,err=800) ia
          read(word(4),*,err=800) val
          Zix=Zinit-iz   
          Nix=Ninit-ia+iz
          if (Zix.lt.0.or.Zix.gt.numZ.or.Nix.lt.0.or.Nix.gt.numN) then
            goto 840
          else
            alimit(Zix,Nix)=val
          endif
          goto 110
        endif
        if (key.eq.'gammald') then
          read(word(2),*,err=800) iz
          read(word(3),*,err=800) ia
          read(word(4),*,err=800) val
          Zix=Zinit-iz   
          Nix=Ninit-ia+iz
          if (Zix.lt.0.or.Zix.gt.numZ.or.Nix.lt.0.or.Nix.gt.numN) then
            goto 840
          else
            gammald(Zix,Nix)=val
          endif
          goto 110
        endif
        if (key.eq.'pair') then
          read(word(2),*,err=800) iz
          read(word(3),*,err=800) ia
          read(word(4),*,err=800) val
          Zix=Zinit-iz   
          Nix=Ninit-ia+iz
          if (Zix.lt.0.or.Zix.gt.numZ.or.Nix.lt.0.or.Nix.gt.numN) then
            goto 840
          else
            pair(Zix,Nix)=val
          endif
          goto 110
        endif
        if (key.eq.'pshift') then
          ibar=0
          read(word(2),*,err=800) iz
          read(word(3),*,err=800) ia
          read(word(4),*,err=800) val
          Zix=Zinit-iz   
          Nix=Ninit-ia+iz
          if (Zix.lt.0.or.Zix.gt.numZ.or.Nix.lt.0.or.Nix.gt.numN) then
            goto 840
          else
            read(word(5),*,end=115,err=800) ibar
            if (ibar.lt.0.or.ibar.gt.numbar) goto 810
  115       Pshift(Zix,Nix,ibar)=val
          endif
          goto 110
        endif
        if (key.eq.'deltaw') then
          ibar=0
          read(word(2),*,err=800) iz
          read(word(3),*,err=800) ia
          read(word(4),*,err=800) val
          Zix=Zinit-iz   
          Nix=Ninit-ia+iz
          if (Zix.lt.0.or.Zix.gt.numZ.or.Nix.lt.0.or.Nix.gt.numN) then
            goto 840
          else
            read(word(5),*,end=120,err=800) ibar
            if (ibar.lt.0.or.ibar.gt.numbar) goto 810
  120       deltaW(Zix,Nix,ibar)=val
          endif
          goto 110
        endif
        if (key.eq.'exmatch') then
          ibar=0
          read(word(2),*,err=800) iz
          read(word(3),*,err=800) ia
          read(word(4),*,err=800) val
          Zix=Zinit-iz   
          Nix=Ninit-ia+iz
          if (Zix.lt.0.or.Zix.gt.numZ.or.Nix.lt.0.or.Nix.gt.numN) then
            goto 840
          else
            read(word(5),*,end=130,err=800) ibar
            if (ibar.lt.0.or.ibar.gt.numbar) goto 810
  130       Exmatch(Zix,Nix,ibar)=val
          endif
          goto 110
        endif
        if (key.eq.'t') then
          ibar=0
          read(word(2),*,err=800) iz
          read(word(3),*,err=800) ia
          read(word(4),*,err=800) val
          Zix=Zinit-iz   
          Nix=Ninit-ia+iz
          if (Zix.lt.0.or.Zix.gt.numZ.or.Nix.lt.0.or.Nix.gt.numN) then
            goto 840
          else
            read(word(5),*,end=140,err=800) ibar
            if (ibar.lt.0.or.ibar.gt.numbar) goto 810
  140       T(Zix,Nix,ibar)=val
          endif
          goto 110
        endif
        if (key.eq.'e0') then
          ibar=0
          read(word(2),*,err=800) iz
          read(word(3),*,err=800) ia
          read(word(4),*,err=800) val
          Zix=Zinit-iz   
          Nix=Ninit-ia+iz
          if (Zix.lt.0.or.Zix.gt.numZ.or.Nix.lt.0.or.Nix.gt.numN) then
            goto 840
          else
            read(word(5),*,end=150,err=800) ibar
            if (ibar.lt.0.or.ibar.gt.numbar) goto 810
  150       E0(Zix,Nix,ibar)=val
          endif
          goto 110
        endif
        if (key.eq.'nlow') then
          ibar=0
          read(word(2),*,err=800) iz
          read(word(3),*,err=800) ia
          read(word(4),*,err=800) ivalue
          Zix=Zinit-iz   
          Nix=Ninit-ia+iz
          if (Zix.lt.0.or.Zix.gt.numZ.or.Nix.lt.0.or.Nix.gt.numN) then
            goto 840
          else
            read(word(5),*,end=160,err=800) ibar
            if (ibar.lt.0.or.ibar.gt.numbar) goto 810
  160       Nlow(Zix,Nix,ibar)=ivalue
          endif
          goto 110
        endif
        if (key.eq.'ntop') then
          ibar=0
          read(word(2),*,err=800) iz
          read(word(3),*,err=800) ia
          read(word(4),*,err=800) ivalue
          Zix=Zinit-iz   
          Nix=Ninit-ia+iz
          if (Zix.lt.0.or.Zix.gt.numZ.or.Nix.lt.0.or.Nix.gt.numN) then
            goto 840
          else
            read(word(5),*,end=170,err=800) ibar
            if (ibar.lt.0.or.ibar.gt.numbar) goto 810
  170       Ntop(Zix,Nix,ibar)=ivalue
          endif
          goto 110
        endif
        if (key.eq.'krotconstant') then
          ibar=0
          read(word(2),*,err=800) iz
          read(word(3),*,err=800) ia
          read(word(4),*,err=800) val
          Zix=Zinit-iz   
          Nix=Ninit-ia+iz
          if (Zix.lt.0.or.Zix.gt.numZ.or.Nix.lt.0.or.Nix.gt.numN) then
            goto 840
          else
            read(word(5),*,end=180,err=800) ibar
            if (ibar.lt.0.or.ibar.gt.numbar) goto 810
  180       Krotconstant(Zix,Nix,ibar)=val
          endif
          goto 110
        endif
        if (key.eq.'c1table') then
          read(word(2),*,err=800) iz
          read(word(3),*,err=800) ia
          read(word(4),*,err=800) val
          Zix=Zinit-iz   
          Nix=Ninit-ia+iz
          if (Zix.lt.0.or.Zix.gt.numZ.or.Nix.lt.0.or.Nix.gt.numN) then
            goto 840
          else
            c1table(Zix,Nix)=val
          endif
          goto 110
        endif
        if (key.eq.'c2table') then
          read(word(2),*,err=800) iz
          read(word(3),*,err=800) ia
          read(word(4),*,err=800) val
          Zix=Zinit-iz   
          Nix=Ninit-ia+iz
          if (Zix.lt.0.or.Zix.gt.numZ.or.Nix.lt.0.or.Nix.gt.numN) then
            goto 840
          else
            c2table(Zix,Nix)=val
          endif
          goto 110
        endif
        if (key.eq.'g') then
          read(word(2),*,err=800) iz
          read(word(3),*,err=800) ia
          read(word(4),*,err=800) val
          Zix=Zinit-iz   
          Nix=Ninit-ia+iz
          if (Zix.lt.0.or.Zix.gt.numZ.or.Nix.lt.0.or.Nix.gt.numN) then
            goto 840
          else
            g(Zix,Nix)=val
          endif
          goto 110
        endif
        if (key.eq.'gp') then
          read(word(2),*,err=800) iz
          read(word(3),*,err=800) ia
          read(word(4),*,err=800) val
          Zix=Zinit-iz   
          Nix=Ninit-ia+iz
          if (Zix.lt.0.or.Zix.gt.numZ.or.Nix.lt.0.or.Nix.gt.numN) then
            goto 840
          else
            gp(Zix,Nix)=val
          endif
          goto 110
        endif
        if (key.eq.'gn') then
          read(word(2),*,err=800) iz
          read(word(3),*,err=800) ia
          read(word(4),*,err=800) val
          Zix=Zinit-iz   
          Nix=Ninit-ia+iz
          if (Zix.lt.0.or.Zix.gt.numZ.or.Nix.lt.0.or.Nix.gt.numN) then
            goto 840
          else
            gn(Zix,Nix)=val
          endif
          goto 110
        endif
        if (key.eq.'egr') then
          igr=1
          read(word(2),*,err=800) iz
          read(word(3),*,err=800) ia
          read(word(4),*,err=800) val
          Zix=Zinit-iz   
          Nix=Ninit-ia+iz
          if (Zix.lt.0.or.Zix.gt.numZ.or.Nix.lt.0.or.Nix.gt.numN) then
            goto 840
          else
            ch=word(5)(1:1)
            if (ch.eq.'m'.or.ch.eq.'e') then
              if (ch.eq.'m') irad=0
              if (ch.eq.'e') irad=1
              read(word(5)(2:2),*,err=800) lval
              if (lval.lt.1.or.lval.gt.numgam) goto 830
              read(word(6),*,end=190,err=800) igr
              if (igr.lt.1.or.igr.gt.2) goto 820
  190         egr(Zix,Nix,irad,lval,igr)=val
              goto 110
            endif
            goto 800
          endif
        endif
        if (key.eq.'ggr') then
          igr=1
          read(word(2),*,err=800) iz
          read(word(3),*,err=800) ia
          read(word(4),*,err=800) val
          Zix=Zinit-iz   
          Nix=Ninit-ia+iz
          if (Zix.lt.0.or.Zix.gt.numZ.or.Nix.lt.0.or.Nix.gt.numN) then
            goto 840
          else
            ch=word(5)(1:1)
            if (ch.eq.'m'.or.ch.eq.'e') then
              if (ch.eq.'m') irad=0
              if (ch.eq.'e') irad=1
              read(word(5)(2:2),*,err=800) lval
              if (lval.lt.1.or.lval.gt.numgam) goto 830
              read(word(6),*,end=200,err=800) igr
              if (igr.lt.1.or.igr.gt.2) goto 820
  200         ggr(Zix,Nix,irad,lval,igr)=val
              goto 110
            endif
            goto 800
          endif
        endif
        if (key.eq.'sgr') then
          igr=1
          read(word(2),*,err=800) iz
          read(word(3),*,err=800) ia
          read(word(4),*,err=800) val
          Zix=Zinit-iz   
          Nix=Ninit-ia+iz
          if (Zix.lt.0.or.Zix.gt.numZ.or.Nix.lt.0.or.Nix.gt.numN) then
            goto 840
          else
            ch=word(5)(1:1)
            if (ch.eq.'m'.or.ch.eq.'e') then
              if (ch.eq.'m') irad=0
              if (ch.eq.'e') irad=1
              read(word(5)(2:2),*,err=800) lval
              if (lval.lt.1.or.lval.gt.numgam) goto 830
              read(word(6),*,end=210,err=800) igr
              if (igr.lt.1.or.igr.gt.2) goto 820
  210         sgr(Zix,Nix,irad,lval,igr)=val
              goto 110
            endif
            goto 800
          endif
        endif
        if (key.eq.'gamgam') then
          read(word(2),*,err=800) iz
          read(word(3),*,err=800) ia
          read(word(4),*,err=800) val
          Zix=Zinit-iz   
          Nix=Ninit-ia+iz
          if (Zix.lt.0.or.Zix.gt.numZ.or.Nix.lt.0.or.Nix.gt.numN) then
            goto 840
          else
            gamgam(Zix,Nix)=val
          endif
          goto 110
        endif
        if (key.eq.'d0') then
          read(word(2),*,err=800) iz
          read(word(3),*,err=800) ia
          read(word(4),*,err=800) val
          Zix=Zinit-iz   
          Nix=Ninit-ia+iz
          if (Zix.lt.0.or.Zix.gt.numZ.or.Nix.lt.0.or.Nix.gt.numN) then
            goto 840
          else
            D0(Zix,Nix)=val
          endif
          goto 110
        endif
        if (key.eq.'s0') then
          read(word(2),*,err=800) iz
          read(word(3),*,err=800) ia
          read(word(4),*,err=800) val
          Zix=Zinit-iz   
          Nix=Ninit-ia+iz
          if (Zix.lt.0.or.Zix.gt.numZ.or.Nix.lt.0.or.Nix.gt.numN) then
            goto 840
          else
            S0(Zix,Nix)=val
          endif
          goto 110
        endif
        if (key.eq.'fisbar') then
          ibar=1
          read(word(2),*,err=800) iz
          read(word(3),*,err=800) ia
          read(word(4),*,err=800) val
          Zix=Zinit-iz
          Nix=Ninit-ia+iz
          if (Zix.lt.0.or.Zix.gt.numZ.or.Nix.lt.0.or.Nix.gt.numN) then
            goto 840
          else
            read(word(5),*,end=220,err=800) ibar
            if (ibar.lt.1.or.ibar.gt.numbar) goto 810
  220       fbarrier(Zix,Nix,ibar)=val
          endif
          goto 110
        endif                   
        if (key.eq.'fishw') then
          ibar=1
          read(word(2),*,err=800) iz
          read(word(3),*,err=800) ia
          read(word(4),*,err=800) val
          Zix=Zinit-iz
          Nix=Ninit-ia+iz
          if (Zix.lt.0.or.Zix.gt.numZ.or.Nix.lt.0.or.Nix.gt.numN) then
            goto 840
          else
            read(word(5),*,end=230,err=800) ibar
            if (ibar.lt.1.or.ibar.gt.numbar) goto 810
  230       fwidth(Zix,Nix,ibar)=val
          endif
          goto 110
        endif                   
        if (key.eq.'rtransmom') then
          ibar=1
          read(word(2),*,err=800) iz
          read(word(3),*,err=800) ia
          read(word(4),*,err=800) val
          Zix=Zinit-iz
          Nix=Ninit-ia+iz
          if (Zix.lt.0.or.Zix.gt.numZ.or.Nix.lt.0.or.Nix.gt.numN) then
            goto 840
          else
            read(word(5),*,end=240,err=800) ibar
            if (ibar.lt.1.or.ibar.gt.numbar) goto 810
  240       Rtransmom(Zix,Nix,ibar)=val
          endif
          goto 110
        endif                   
        if (key.eq.'rclass2mom') then
          ibar=1
          read(word(2),*,err=800) iz
          read(word(3),*,err=800) ia
          read(word(4),*,err=800) val
          Zix=Zinit-iz
          Nix=Ninit-ia+iz
          if (Zix.lt.0.or.Zix.gt.numZ.or.Nix.lt.0.or.Nix.gt.numN) then
            goto 840
          else
            read(word(5),*,end=250,err=800) ibar
            if (ibar.lt.1.or.ibar.gt.numbar) goto 810
  250       Rclass2mom(Zix,Nix,ibar)=val
          endif
          goto 110
        endif                   
        if (key.eq.'class2width') then
          ibar=1
          read(word(2),*,err=800) iz
          read(word(3),*,err=800) ia
          read(word(4),*,err=800) val
          Zix=Zinit-iz
          Nix=Ninit-ia+iz
          if (Zix.lt.0.or.Zix.gt.numZ.or.Nix.lt.0.or.Nix.gt.numN) then
            goto 840
          else
            read(word(5),*,end=260,err=800) ibar
            if (ibar.lt.1.or.ibar.gt.numbar) goto 810
  260       widthc2(Zix,Nix,ibar)=val
          endif
          goto 110
        endif                   
        if (key.eq.'levelfile') then
          read(word(2),*,err=800) iz
          Zix=Zinit-iz
          if (Zix.lt.0.or.Zix.gt.numZ) then
            goto 840
          else
            levelfile(Zix)=word(3)
            goto 110
          endif
        endif                   
        if (key.eq.'deformfile') then
          read(word(2),*,err=800) iz
          Zix=Zinit-iz
          if (Zix.lt.0.or.Zix.gt.numZ) then
            goto 840
          else
            deformfile(Zix)=word(3)
            goto 110
          endif
        endif                   
        if (key.eq.'hbtransfile') then
          read(word(2),*,err=800) iz
          read(word(3),*,err=800) ia
          Zix=Zinit-iz
          Nix=Ninit-ia+iz
          if (Zix.lt.0.or.Zix.gt.numZ.or.Nix.lt.0.or.Nix.gt.numN) then
            goto 840
          else
            hbtransfile(Zix,Nix)=word(4)
            goto 110
          endif
        endif                   
        if (key.eq.'class2file') then
          read(word(2),*,err=800) iz
          read(word(3),*,err=800) ia
          Zix=Zinit-iz
          Nix=Ninit-ia+iz
          if (Zix.lt.0.or.Zix.gt.numZ.or.Nix.lt.0.or.Nix.gt.numN) then
            goto 840
          else
            class2file(Zix,Nix)=word(4)
            goto 110
          endif
        endif                   
        if (key.eq.'optmodfilen') then
          read(word(2),*,err=800) iz
          Zix=Zinit-iz
          if (Zix.lt.0.or.Zix.gt.numZ) then
            goto 840
          else
            optmodfileN(Zix)=word(3)
            goto 110
          endif
        endif                   
        if (key.eq.'optmodfilep') then
          read(word(2),*,err=800) iz
          Zix=Zinit-iz
          if (Zix.lt.0.or.Zix.gt.numZ) then
            goto 840
          else
            optmodfileP(Zix)=word(3)
            goto 110
          endif
        endif                   
        if (key.eq.'beta2') then
          ibar=0
          read(word(2),*,err=800) iz
          read(word(3),*,err=800) ia
          read(word(4),*,err=800) val
          Zix=Zinit-iz   
          Nix=Ninit-ia+iz
          if (Zix.lt.0.or.Zix.gt.numZ.or.Nix.lt.0.or.Nix.gt.numN) then
            goto 840
          else
            read(word(5),*,end=270,err=800) ibar
            if (ibar.lt.0.or.ibar.gt.numbar) goto 810
  270       beta2(Zix,Nix,ibar)=val
          endif
          goto 110
        endif
        if (key.eq.'axtype') then
          ibar=1
          read(word(2),*,err=800) iz
          read(word(3),*,err=800) ia
          read(word(4),*,err=800) ivalue
          Zix=Zinit-iz   
          Nix=Ninit-ia+iz
          if (Zix.lt.0.or.Zix.gt.numZ.or.Nix.lt.0.or.Nix.gt.numN) then
            goto 840
          else
            read(word(5),*,end=280,err=800) ibar
            if (ibar.lt.1.or.ibar.gt.numbar) goto 810
  280       axtype(Zix,Nix,ibar)=ivalue
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
          goto 800
  300     read(word(3),*,err=800) Cstrip(type2)
          goto 110
        endif
        if (key.eq.'cknock') then
          do 310 type=0,6
            if (ch.eq.parsym(type)) then
              type2=type
              goto 320
            endif
  310     continue
          goto 800
  320     read(word(3),*,err=800) Cknock(type2)
          goto 110
        endif
        if (key.eq.'v1adjust') then
          do 330 type=1,6
            if (ch.eq.parsym(type)) then
              type2=type
              goto 340
            endif
  330     continue
          goto 800
  340     read(word(3),*,err=800) v1adjust(type2)
          goto 110
        endif
        if (key.eq.'v2adjust') then
          do 350 type=1,6
            if (ch.eq.parsym(type)) then
              type2=type
              goto 360
            endif
  350     continue
          goto 800
  360     read(word(3),*,err=800) v2adjust(type2)
          goto 110
        endif
        if (key.eq.'v3adjust') then
          do 370 type=1,6
            if (ch.eq.parsym(type)) then
              type2=type
              goto 380
            endif
  370     continue
          goto 800
  380     read(word(3),*,err=800) v3adjust(type2)
          goto 110
        endif
        if (key.eq.'v4adjust') then
          do 390 type=1,6
            if (ch.eq.parsym(type)) then
              type2=type
              goto 400
            endif
  390     continue
          goto 800
  400     read(word(3),*,err=800) v4adjust(type2)
          goto 110
        endif
        if (key.eq.'rvadjust') then
          do 410 type=1,6
            if (ch.eq.parsym(type)) then
              type2=type
              goto 420
            endif
  410     continue
          goto 800
  420     read(word(3),*,err=800) rvadjust(type2)
          goto 110
        endif
        if (key.eq.'avadjust') then
          do 430 type=1,6
            if (ch.eq.parsym(type)) then
              type2=type
              goto 440
            endif
  430     continue
          goto 800
  440     read(word(3),*,err=800) avadjust(type2)
          goto 110
        endif
        if (key.eq.'w1adjust') then
          do 450 type=1,6
            if (ch.eq.parsym(type)) then
              type2=type
              goto 460
            endif
  450     continue
          goto 800
  460     read(word(3),*,err=800) w1adjust(type2)
          goto 110
        endif
        if (key.eq.'w2adjust') then
          do 470 type=1,6
            if (ch.eq.parsym(type)) then
              type2=type
              goto 480
            endif
  470     continue
          goto 800
  480     read(word(3),*,err=800) w2adjust(type2)
          goto 110
        endif
        if (key.eq.'d1adjust') then
          do 490 type=1,6
            if (ch.eq.parsym(type)) then
              type2=type
              goto 500
            endif
  490     continue
          goto 800
  500     read(word(3),*,err=800) d1adjust(type2)
          goto 110
        endif
        if (key.eq.'d2adjust') then
          do 510 type=1,6
            if (ch.eq.parsym(type)) then
              type2=type
              goto 520
            endif
  510     continue
          goto 800
  520     read(word(3),*,err=800) d2adjust(type2)
          goto 110
        endif
        if (key.eq.'d3adjust') then
          do 530 type=1,6
            if (ch.eq.parsym(type)) then
              type2=type
              goto 540
            endif
  530     continue
          goto 800
  540     read(word(3),*,err=800) d3adjust(type2)
          goto 110
        endif
        if (key.eq.'rvsoadjust') then
          do 550 type=1,6
            if (ch.eq.parsym(type)) then
              type2=type
              goto 560
            endif
  550     continue
          goto 800
  560     read(word(3),*,err=800) rvsoadjust(type2)
          goto 110
        endif
        if (key.eq.'avsoadjust') then
          do 570 type=1,6
            if (ch.eq.parsym(type)) then
              type2=type
              goto 580
            endif
  570     continue
          goto 800
  580     read(word(3),*,err=800) avsoadjust(type2)
          goto 110
        endif
        if (key.eq.'rvdadjust') then
          do 590 type=1,6
            if (ch.eq.parsym(type)) then
              type2=type
              goto  600
            endif
  590     continue
          goto 800
  600     read(word(3),*,err=800) rvdadjust(type2)
          goto 110
        endif
        if (key.eq.'avdadjust') then
          do 610 type=1,6
            if (ch.eq.parsym(type)) then
              type2=type
              goto 620
            endif
  610     continue
          goto 800
  620     read(word(3),*,err=800) avdadjust(type2)
          goto 110
        endif
        if (key.eq.'vso1adjust') then
          do 630 type=1,6
            if (ch.eq.parsym(type)) then
              type2=type
              goto 640
            endif
  630     continue
          goto 800
  640     read(word(3),*,err=800) vso1adjust(type2)
          goto 110
        endif
        if (key.eq.'vso2adjust') then
          do 650 type=1,6
            if (ch.eq.parsym(type)) then
              type2=type
              goto 660
            endif
  650     continue
          goto 800
  660     read(word(3),*,err=800) vso2adjust(type2)
          goto 110
        endif
        if (key.eq.'wso1adjust') then
          do 670 type=1,6
            if (ch.eq.parsym(type)) then
              type2=type
              goto 680
            endif
  670     continue
          goto 800
  680     read(word(3),*,err=800) wso1adjust(type2)
          goto 110
        endif
        if (key.eq.'wso2adjust') then
          do 690 type=1,6
            if (ch.eq.parsym(type)) then
              type2=type
              goto 700
            endif
  690     continue
          goto 800
  700     read(word(3),*,err=800) wso2adjust(type2)
          goto 110
        endif
        if (key.eq.'rcadjust') then
          do 710 type=1,6
            if (ch.eq.parsym(type)) then
              type2=type
              goto 720
            endif
  710     continue
          goto 800
  720     read(word(3),*,err=800) rcadjust(type2)
          goto 110
        endif
        if (key.eq.'gnorm') then
          read(value,*,err=800) gnorm 
          goto 110
        endif
        if (key.eq.'spincutmodel') then
          read(value,*,err=800) spincutmodel 
          goto 110
        endif
        if (key.eq.'shellmodel') then
          read(value,*,err=800) shellmodel 
          goto 110
        endif
        if (key.eq.'rspincut') then
          read(value,*,err=800) Rspincut 
          goto 110
        endif
        if (key.eq.'alphald') then
          read(value,*,err=800) alphald 
          goto 110
        endif
        if (key.eq.'betald') then
          read(value,*,err=800) betald 
          goto 110
        endif
        if (key.eq.'gammashell1') then
          read(value,*,err=800) gammashell1
          goto 110
        endif
        if (key.eq.'gammashell2') then
          read(value,*,err=800) gammashell2
          goto 110
        endif
        if (key.eq.'pairconstant') then
          read(value,*,err=800) pairconstant
          goto 110
        endif
        if (key.eq.'pshiftconstant') then
          read(value,*,err=800) Pshiftconstant
          goto 110
        endif
        if (key.eq.'ufermi') then
          read(value,*,err=800) Ufermi 
          goto 110
        endif
        if (key.eq.'cfermi') then
          read(value,*,err=800) cfermi 
          goto 110
        endif
        if (key.eq.'ufermibf') then
          read(value,*,err=800) Ufermibf 
          goto 110
        endif
        if (key.eq.'cfermibf') then
          read(value,*,err=800) cfermibf 
          goto 110
        endif
        if (key.eq.'kph') then
          read(value,*,err=800) Kph 
          goto 110
        endif
        if (key.eq.'m2constant') then
          read(value,*,err=800) M2constant 
          goto 110
        endif
        if (key.eq.'m2limit') then
          read(value,*,err=800) M2limit 
          goto 110
        endif
        if (key.eq.'m2shift') then
          read(value,*,err=800) M2shift 
          goto 110
        endif
        if (key.eq.'rpipi') then
          read(value,*,err=800) Rpipi 
          goto 110
        endif
        if (key.eq.'rnunu') then
          read(value,*,err=800) Rnunu 
          goto 110
        endif
        if (key.eq.'rpinu') then
          read(value,*,err=800) Rpinu 
          goto 110
        endif
        if (key.eq.'rnupi') then
          read(value,*,err=800) Rnupi 
          goto 110
        endif
        if (key.eq.'esurf') then
          read(value,*,err=800) Esurf0
          goto 110
        endif
        if (key.eq.'rgamma') then
          read(value,*,err=800) Rgamma 
          goto 110
        endif
        if (key.eq.'msdbins')  then
          read(value,*,err=800) msdbins          
          goto 110
        endif
        if (key.eq.'emsdmin')  then
          read(value,*,err=800) Emsdmin          
          goto 110
        endif
        if (key.eq.'elwidth')  then
          read(value,*,err=800) elwidth          
          goto 110
        endif
        if (key.eq.'anglesrec')  then
          read(value,*,err=800) nanglerec
          goto 110
        endif
        goto 110
  840   write(*,'(" TALYS-warning: Z,N index out of range,",
     +    " keyword ignored: ",a80)') inline(i)
  110 continue
      return
  800 write(*,'(" TALYS-error: Wrong input: ",a80)') inline(i)
      stop
  810 write(*,'(" TALYS-error: 0(1) <= fission barrier <=",i3,
     +  ", ibar index out of range: ",a80)') numbar,inline(i)
      stop
  820 write(*,'(" TALYS-error: 0 <= resonance number <= 2",
     +  ", igr index out of range: ",a80)') inline(i)
      stop
  830 write(*,'(" TALYS-error: 0 <= multipole radiation <= ",i1,
     +  ", lval index out of range: ",a80)') numgam,inline(i)
      stop
      end
Copyright (C) 2004  A.J. Koning, S. Hilaire and M.C. Duijvestijn
