      subroutine input5
c
c +---------------------------------------------------------------------
c | Author: Arjan Koning 
c | Date  : December 1, 2004
c | Task  : Read input for fifth set of variables
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      character*1 ch
      integer     Zix,Nix,ibar,irad,lval,igr,i,iz,ia,i2,i3,ivalue
      real        value
c
c ************** Defaults for fifth set of input variables *************
c
c Atarget     : mass number of target nucleus
c gnorm       : gamma normalization factor
c Rspincut    : adjustable constant for spin cutoff factor
c ldmodel     : level density model 
c alphad      : alpha-constant for asymptotic level density parameter
c betald      : beta-constant for asymptotic level density parameter    
c gammashell1 : gamma-constant for asymptotic level density parameter 
c gammashell2 : gamma-constant for asymptotic level density parameter 
c Ufermi      : energy of Fermi distribution for damping of ground-state
c             : rotational effects
c cfermi      : width of Fermi distribution for damping of ground-state
c             : rotational effects
c Ufermibf    : energy of Fermi distribution for damping of barrier
c             : rotational effects
c cfermibf    : width of Fermi distribution for damping of barrier
c             : rotational effects              
c Kph         : constant for single-particle level density parameter 
c               (g=A/Kph)
c M2constant  : overall constant for matrix element in exciton model
c M2limit     : constant for asymptotical value for matrix element
c M2shift     : constant for energy shift for matrix element
c Rpinu       : ratio for two-component matrix element
c Esurf0      : well depth for surface interaction
c Rgamma      : adjustable parameter for pre-equilibrium gamma decay    
c elwidth     : width of elastic peak in MeV
c Zix         : charge number index for residual nucleus
c numZ        : maximal number of protons away from the initial compound
c               nucleus
c Nix         : neutron number index for residual nucleus
c numN        : maximal number of neutrons away from the initial 
c               compound nucleus
c alev        : level density parameter
c alimit      : asymptotic level density parameter
c gammald     : gamma-constant for asymptotic level density parameter
c pair        : total pairing correction
c ibar        : fission barrier
c numbar      : number of fission barriers
c deltaW      : shell correction in nuclear mass
c Exmatch     : matching point for Ex
c T           : nuclear temperature
c E0          : constant of temperature formula
c Nlow        : lowest discrete level for temperature matching
c Ntop        : highest discrete level for temperature matching
c beta2       : deformation parameter
c Krotconstant: normalization constant for rotational enhancement
c c1table     : constant to adjust tabulated level densities
c c2table     : constant to adjust tabulated level densities
c g           : single-particle level density parameter
c gp          : single-particle proton level density parameter
c gn          : single-particle neutron level density parameter
c gamgam      : experimental total radiative width in eV
c D0          : experimental s-wave resonance spacing in eV
c S0          : s-wave strength function
c irad        : variable to indicate M(=0) or E(=1) radiation
c lval        : multipolarity
c numgam      : maximum number of l-values for gamma multipolarity
c igr         : giant resonance
c egr         : energy of GR
c ggr         : width of GR
c sgr         : strength of GR
c axtype      : type of axiality of barrier (1: axial, 2: tri-axial)
c fbarrier    : height of fission barrier
c fwidth      : width of fission barrier 
c Rtransmom   : normalization constant for moment of inertia for 
c               transition states
c Rclass2mom  : normalization constant for moment of inertia for 
c               class 2 states
c widthc2     : width of class2 states
c Ntarget     : neutron number of target nucleus
c hbtransfile : file with head band transition states
c class2file  : file with class 2 transition states
c levelfile   : discrete level file   
c deformfile  : deformation parameter file   
c optmodfileN : optical model parameter file for neutrons
c optmodfileP : optical model parameter file for protons
c msdbins     : number of energy points for DWBA calculation for MSD
c Emsdmin     : minimal outgoing energy for MSD calculation
c flaglabddx  : flag for calculation of DDX in LAB system
c nanglerec   : number of recoil angles
c
c Advice of S. Goriely: no gamma normalization for A < 40.
c
      if (Atarget.lt.40) then
        gnorm=1.
      else
        gnorm=-1.
      endif
      Rspincut=1.
c
c Energy dependent level density systematics from M. Duijvestijn.
c 
      if (ldmodel.eq.1) then
        alphald=0.0666
        betald=0.258
        gammashell1=0.459
        gammashell2=0.
      else
        alphald=0.0564
        betald=0.083
        gammashell1=0.620
        gammashell2=0.
      endif                              
      Ufermi=30.
      cfermi=10.
      Ufermibf=45.
      cfermibf=10.
      Kph=15.
      M2constant=1.0
      M2limit=1.0
      M2shift=1.0
      Rpinu=1.5
      Esurf0=-1.
      Rgamma=2.
      elwidth=0.5
      do 10 Zix=0,numZ
        do 20 Nix=0,numN
          alev(Zix,Nix)=0.
          alimit(Zix,Nix)=0.
          gammald(Zix,Nix)=-1.
          pair(Zix,Nix)=0.
          do 30 ibar=0,numbar
            deltaW(Zix,Nix,ibar)=0.
            Exmatch(Zix,Nix,ibar)=0.
            T(Zix,Nix,ibar)=0.
            E0(Zix,Nix,ibar)=0.
            Nlow(Zix,Nix,ibar)=-1
            Ntop(Zix,Nix,ibar)=-1
            Krotconstant(Zix,Nix,ibar)=1.
   30     continue     
          beta2(Zix,Nix,0)=0.
          beta2(Zix,Nix,1)=0.6
          beta2(Zix,Nix,2)=0.8
          beta2(Zix,Nix,3)=1. 
          c1table(Zix,Nix)=1.
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
      if (flaglabddx) then
        nanglerec=9
      else
        nanglerec=1
      endif
c
c **************** Read fifth set of input variables *******************
c
c Here, the various model parameters can be set to overrule the default
c values. Most default values will be computed later on, since they
c require more computation (e.g. level density parameters).
c
c nlines: number of input lines
c inline: input line                 
c iz    : charge number
c ia    : mass number
c value : help variable 
c Zinit : charge number of initial compound nucleus
c Ninit : neutron number of initial compound nucleus
c ivalue: help variable 
c Ainit : mass number of initial compound nucleus
c
c First the keyword is identified. Next the corresponding value is read.
c Erroneous input is immediately checked.
c
      do 110 i=1,nlines
        if (inline(i)(1:2).eq.'a ') then
          read(inline(i)(3:80),*,err=600) iz,ia,value
          Zix=Zinit-iz   
          Nix=Ninit-ia+iz
          if (Zix.lt.0.or.Zix.gt.numZ.or.Nix.lt.0.or.Nix.gt.numN) then
            goto 640
          else
            alev(Zix,Nix)=value
          endif
          goto 110
        endif
        if (inline(i)(1:7).eq.'alimit ') then
          read(inline(i)(8:80),*,err=600) iz,ia,value
          Zix=Zinit-iz   
          Nix=Ninit-ia+iz
          if (Zix.lt.0.or.Zix.gt.numZ.or.Nix.lt.0.or.Nix.gt.numN) then
            goto 640
          else
            alimit(Zix,Nix)=value
          endif
          goto 110
        endif
        if (inline(i)(1:8).eq.'gammald ') then
          read(inline(i)(9:80),*,err=600) iz,ia,value
          Zix=Zinit-iz   
          Nix=Ninit-ia+iz
          if (Zix.lt.0.or.Zix.gt.numZ.or.Nix.lt.0.or.Nix.gt.numN) then
            goto 640
          else
            gammald(Zix,Nix)=value
          endif
          goto 110
        endif
        if (inline(i)(1:5).eq.'pair ') then
          read(inline(i)(6:80),*,err=600) iz,ia,value
          Zix=Zinit-iz   
          Nix=Ninit-ia+iz
          if (Zix.lt.0.or.Zix.gt.numZ.or.Nix.lt.0.or.Nix.gt.numN) then
            goto 640
          else
            pair(Zix,Nix)=value
          endif
          goto 110
        endif
        if (inline(i)(1:7).eq.'deltaw ') then
          ibar=0
          read(inline(i)(8:80),*,err=600) iz,ia,value
          Zix=Zinit-iz   
          Nix=Ninit-ia+iz
          if (Zix.lt.0.or.Zix.gt.numZ.or.Nix.lt.0.or.Nix.gt.numN) then
            goto 640
          else
            read(inline(i)(8:80),*,end=120,err=600) iz,ia,value,ibar
            if (ibar.lt.0.or.ibar.gt.numbar) goto 610
  120       deltaW(Zix,Nix,ibar)=value
          endif
          goto 110
        endif
        if (inline(i)(1:8).eq.'exmatch ') then
          ibar=0
          read(inline(i)(9:80),*,err=600) iz,ia,value
          Zix=Zinit-iz   
          Nix=Ninit-ia+iz
          if (Zix.lt.0.or.Zix.gt.numZ.or.Nix.lt.0.or.Nix.gt.numN) then
            goto 640
          else
            read(inline(i)(9:80),*,end=130,err=600) iz,ia,value,ibar
            if (ibar.lt.0.or.ibar.gt.numbar) goto 610
  130       Exmatch(Zix,Nix,ibar)=value
          endif
          goto 110
        endif
        if (inline(i)(1:2).eq.'t ') then
          ibar=0
          read(inline(i)(3:80),*,err=600) iz,ia,value
          Zix=Zinit-iz   
          Nix=Ninit-ia+iz
          if (Zix.lt.0.or.Zix.gt.numZ.or.Nix.lt.0.or.Nix.gt.numN) then
            goto 640
          else
            read(inline(i)(3:80),*,end=140,err=600) iz,ia,value,ibar
            if (ibar.lt.0.or.ibar.gt.numbar) goto 610
  140       T(Zix,Nix,ibar)=value
          endif
          goto 110
        endif
        if (inline(i)(1:3).eq.'e0 ') then
          ibar=0
          read(inline(i)(4:80),*,err=600) iz,ia,value
          Zix=Zinit-iz   
          Nix=Ninit-ia+iz
          if (Zix.lt.0.or.Zix.gt.numZ.or.Nix.lt.0.or.Nix.gt.numN) then
            goto 640
          else
            read(inline(i)(4:80),*,end=150,err=600) iz,ia,value,ibar
            if (ibar.lt.0.or.ibar.gt.numbar) goto 610
  150       E0(Zix,Nix,ibar)=value
          endif
          goto 110
        endif
        if (inline(i)(1:5).eq.'nlow ') then
          ibar=0
          read(inline(i)(6:80),*,err=600) iz,ia,ivalue
          Zix=Zinit-iz   
          Nix=Ninit-ia+iz
          if (Zix.lt.0.or.Zix.gt.numZ.or.Nix.lt.0.or.Nix.gt.numN) then
            goto 640
          else
            read(inline(i)(6:80),*,end=160,err=600) iz,ia,ivalue,ibar
            if (ibar.lt.0.or.ibar.gt.numbar) goto 610
  160       Nlow(Zix,Nix,ibar)=ivalue
          endif
          goto 110
        endif
        if (inline(i)(1:5).eq.'ntop ') then
          ibar=0
          read(inline(i)(6:80),*,err=600) iz,ia,ivalue
          Zix=Zinit-iz   
          Nix=Ninit-ia+iz
          if (Zix.lt.0.or.Zix.gt.numZ.or.Nix.lt.0.or.Nix.gt.numN) then
            goto 640
          else
            read(inline(i)(6:80),*,end=170,err=600) iz,ia,ivalue,ibar
            if (ibar.lt.0.or.ibar.gt.numbar) goto 610
  170       Ntop(Zix,Nix,ibar)=ivalue
          endif
          goto 110
        endif
        if (inline(i)(1:13).eq.'krotconstant ') then
          ibar=0
          read(inline(i)(14:80),*,err=600) iz,ia,value
          Zix=Zinit-iz   
          Nix=Ninit-ia+iz
          if (Zix.lt.0.or.Zix.gt.numZ.or.Nix.lt.0.or.Nix.gt.numN) then
            goto 640
          else
            read(inline(i)(14:80),*,end=180,err=600) iz,ia,value,ibar
            if (ibar.lt.0.or.ibar.gt.numbar) goto 610
  180       Krotconstant(Zix,Nix,ibar)=value
          endif
          goto 110
        endif
        if (inline(i)(1:8).eq.'c1table ') then
          read(inline(i)(9:80),*,err=600) iz,ia,value
          Zix=Zinit-iz   
          Nix=Ninit-ia+iz
          if (Zix.lt.0.or.Zix.gt.numZ.or.Nix.lt.0.or.Nix.gt.numN) then
            goto 640
          else
            c1table(Zix,Nix)=value
          endif
          goto 110
        endif
        if (inline(i)(1:8).eq.'c2table ') then
          read(inline(i)(9:80),*,err=600) iz,ia,value
          Zix=Zinit-iz   
          Nix=Ninit-ia+iz
          if (Zix.lt.0.or.Zix.gt.numZ.or.Nix.lt.0.or.Nix.gt.numN) then
            goto 640
          else
            c2table(Zix,Nix)=value
          endif
          goto 110
        endif
        if (inline(i)(1:2).eq.'g ') then
          read(inline(i)(3:80),*,err=600) iz,ia,value
          Zix=Zinit-iz   
          Nix=Ninit-ia+iz
          if (Zix.lt.0.or.Zix.gt.numZ.or.Nix.lt.0.or.Nix.gt.numN) then
            goto 640
          else
            g(Zix,Nix)=value
          endif
          goto 110
        endif
        if (inline(i)(1:3).eq.'gp ') then
          read(inline(i)(4:80),*,err=600) iz,ia,value
          Zix=Zinit-iz   
          Nix=Ninit-ia+iz
          if (Zix.lt.0.or.Zix.gt.numZ.or.Nix.lt.0.or.Nix.gt.numN) then
            goto 640
          else
            gp(Zix,Nix)=value
          endif
          goto 110
        endif
        if (inline(i)(1:3).eq.'gn ') then
          read(inline(i)(4:80),*,err=600) iz,ia,value
          Zix=Zinit-iz   
          Nix=Ninit-ia+iz
          if (Zix.lt.0.or.Zix.gt.numZ.or.Nix.lt.0.or.Nix.gt.numN) then
            goto 640
          else
            gn(Zix,Nix)=value
          endif
          goto 110
        endif
        if (inline(i)(1:4).eq.'egr ') then
          igr=1
          read(inline(i)(5:80),*,err=600) iz,ia,value
          Zix=Zinit-iz   
          Nix=Ninit-ia+iz
          if (Zix.lt.0.or.Zix.gt.numZ.or.Nix.lt.0.or.Nix.gt.numN) then
            goto 640
          else
            do 190 i2=5,80
              ch=inline(i)(i2:i2)
              if (ch.eq.'m'.or.ch.eq.'e') then
                if (inline(i)(i2+1:i2+1).eq.'0') goto 190
                if (inline(i)(i2+1:i2+1).eq.'-') goto 190
                if (inline(i)(i2+1:i2+1).eq.'+') goto 190
                if (ch.eq.'m') irad=0
                if (ch.eq.'e') irad=1
                read(inline(i)(i2+1:i2+1),*,err=600) lval
                if (lval.lt.1.or.lval.gt.numgam) goto 630
                read(inline(i)(i2+2:80),*,end=200,err=600) igr
                if (igr.lt.1.or.igr.gt.2) goto 620
  200           egr(Zix,Nix,irad,lval,igr)=value
                goto 110
              endif
  190       continue
            goto 600
          endif
        endif
        if (inline(i)(1:4).eq.'ggr ') then
          igr=1
          read(inline(i)(5:80),*,err=600) iz,ia,value
          Zix=Zinit-iz   
          Nix=Ninit-ia+iz
          if (Zix.lt.0.or.Zix.gt.numZ.or.Nix.lt.0.or.Nix.gt.numN) then
            goto 640
          else
            do 210 i2=5,80
              ch=inline(i)(i2:i2)
              if (ch.eq.'m'.or.ch.eq.'e') then
                if (inline(i)(i2+1:i2+1).eq.'0') goto 210
                if (inline(i)(i2+1:i2+1).eq.'-') goto 210
                if (inline(i)(i2+1:i2+1).eq.'+') goto 210
                if (ch.eq.'m') irad=0
                if (ch.eq.'e') irad=1
                read(inline(i)(i2+1:i2+1),*,err=600) lval
                if (lval.lt.1.or.lval.gt.numgam) goto 630
                read(inline(i)(i2+2:80),*,end=220,err=600) igr
                if (igr.lt.1.or.igr.gt.2) goto 620
  220           ggr(Zix,Nix,irad,lval,igr)=value
                goto 110
              endif
  210       continue
            goto 600
          endif
        endif
        if (inline(i)(1:4).eq.'sgr ') then
          igr=1
          read(inline(i)(5:80),*,err=600) iz,ia,value
          Zix=Zinit-iz   
          Nix=Ninit-ia+iz
          if (Zix.lt.0.or.Zix.gt.numZ.or.Nix.lt.0.or.Nix.gt.numN) then
            goto 640
          else
            do 230 i2=5,80
              ch=inline(i)(i2:i2)
              if (ch.eq.'m'.or.ch.eq.'e') then
                if (inline(i)(i2+1:i2+1).eq.'0') goto 230
                if (inline(i)(i2+1:i2+1).eq.'-') goto 230
                if (inline(i)(i2+1:i2+1).eq.'+') goto 230
                if (ch.eq.'m') irad=0
                if (ch.eq.'e') irad=1
                read(inline(i)(i2+1:i2+1),*,err=600) lval
                if (lval.lt.1.or.lval.gt.numgam) goto 630
                read(inline(i)(i2+2:80),*,end=240,err=600) igr
                if (igr.lt.1.or.igr.gt.2) goto 620
  240           sgr(Zix,Nix,irad,lval,igr)=value
                goto 110
              endif
  230       continue
            goto 600
          endif
        endif
        if (inline(i)(1:7).eq.'gamgam ') then
          read(inline(i)(8:80),*,err=600) iz,ia,value
          Zix=Zinit-iz   
          Nix=Ninit-ia+iz
          if (Zix.lt.0.or.Zix.gt.numZ.or.Nix.lt.0.or.Nix.gt.numN) then
            goto 640
          else
            gamgam(Zix,Nix)=value
          endif
          goto 110
        endif
        if (inline(i)(1:3).eq.'d0 ') then
          read(inline(i)(4:80),*,err=600) iz,ia,value
          Zix=Zinit-iz   
          Nix=Ninit-ia+iz
          if (Zix.lt.0.or.Zix.gt.numZ.or.Nix.lt.0.or.Nix.gt.numN) then
            goto 640
          else
            D0(Zix,Nix)=value
          endif
          goto 110
        endif
        if (inline(i)(1:3).eq.'s0 ') then
          read(inline(i)(4:80),*,err=600) iz,ia,value
          Zix=Zinit-iz   
          Nix=Ninit-ia+iz
          if (Zix.lt.0.or.Zix.gt.numZ.or.Nix.lt.0.or.Nix.gt.numN) then
            goto 640
          else
            S0(Zix,Nix)=value
          endif
          goto 110
        endif
        if (inline(i)(1:7).eq.'fisbar ') then
          ibar=1
          read(inline(i)(8:80),*,err=600) iz,ia,value
          Zix=Zinit-iz
          Nix=Ninit-ia+iz
          if (Zix.lt.0.or.Zix.gt.numZ.or.Nix.lt.0.or.Nix.gt.numN) then
            goto 640
          else
            read(inline(i)(8:80),*,end=250,err=600) iz,ia,value,ibar
            if (ibar.lt.1.or.ibar.gt.numbar) goto 610
  250       fbarrier(Zix,Nix,ibar)=value
          endif
          goto 110
        endif                   
        if (inline(i)(1:6).eq.'fishw ') then
          ibar=1
          read(inline(i)(7:80),*,err=600) iz,ia,value
          Zix=Zinit-iz
          Nix=Ninit-ia+iz
          if (Zix.lt.0.or.Zix.gt.numZ.or.Nix.lt.0.or.Nix.gt.numN) then
            goto 640
          else
            read(inline(i)(7:80),*,end=260,err=600) iz,ia,value,ibar
            if (ibar.lt.1.or.ibar.gt.numbar) goto 610
  260       fwidth(Zix,Nix,ibar)=value
          endif
          goto 110
        endif                   
        if (inline(i)(1:10).eq.'rtransmom ') then
          ibar=1
          read(inline(i)(11:80),*,err=600) iz,ia,value
          Zix=Zinit-iz
          Nix=Ninit-ia+iz
          if (Zix.lt.0.or.Zix.gt.numZ.or.Nix.lt.0.or.Nix.gt.numN) then
            goto 640
          else
            read(inline(i)(11:80),*,end=270,err=600) iz,ia,value,ibar
            if (ibar.lt.1.or.ibar.gt.numbar) goto 610
  270       Rtransmom(Zix,Nix,ibar)=value
          endif
          goto 110
        endif                   
        if (inline(i)(1:11).eq.'rclass2mom ') then
          ibar=1
          read(inline(i)(12:80),*,err=600) iz,ia,value
          Zix=Zinit-iz
          Nix=Ninit-ia+iz
          if (Zix.lt.0.or.Zix.gt.numZ.or.Nix.lt.0.or.Nix.gt.numN) then
            goto 640
          else
            read(inline(i)(12:80),*,end=280,err=600) iz,ia,value,ibar
            if (ibar.lt.1.or.ibar.gt.numbar) goto 610
  280       Rclass2mom(Zix,Nix,ibar)=value
          endif
          goto 110
        endif                   
        if (inline(i)(1:12).eq.'class2width ') then
          ibar=1
          read(inline(i)(13:80),*,err=600) iz,ia,value
          Zix=Zinit-iz
          Nix=Ninit-ia+iz
          if (Zix.lt.0.or.Zix.gt.numZ.or.Nix.lt.0.or.Nix.gt.numN) then
            goto 640
          else
            read(inline(i)(13:80),*,end=290,err=600) iz,ia,value,ibar
            if (ibar.lt.1.or.ibar.gt.numbar) goto 610
  290       widthc2(Zix,Nix,ibar)=value
          endif
          goto 110
        endif                   
        if (inline(i)(1:10).eq.'levelfile ') then
          do 300 i2=11,80
            ch=inline(i)(i2:i2)
            if (ch.ne.' ') then
              read(inline(i)(i2:80),*,err=600) iz
              Zix=Zinit-iz
              if (Zix.lt.0.or.Zix.gt.numZ) then
                goto 640
              else
                goto 310
              endif
            endif
  300     continue                
  310     do 320 i2=11,80
            ch=inline(i)(i2:i2)
            if (ch.ne.' '.and.(ch.lt.'0'.or.ch.gt.'9')) then
              do 330 i3=i2,80
                ch=inline(i)(i3:i3)
                if (ch.eq.' ') then
                  levelfile(Zix)=inline(i)(i2:i3-1)
                  goto 110
                endif
  330         continue                
            endif
  320     continue                
        endif                   
        if (inline(i)(1:11).eq.'deformfile ') then
          do 340 i2=12,80
            ch=inline(i)(i2:i2)
            if (ch.ne.' ') then
              read(inline(i)(i2:80),*,err=600) iz
              Zix=Zinit-iz
              if (Zix.lt.0.or.Zix.gt.numZ) then
                goto 640
              else
                goto 350
              endif
            endif
  340     continue                
  350     do 360 i2=12,80
            ch=inline(i)(i2:i2)
            if (ch.ne.' '.and.(ch.lt.'0'.or.ch.gt.'9')) then
              do 370 i3=i2,80
                ch=inline(i)(i3:i3)
                if (ch.eq.' ') then
                  deformfile(Zix)=inline(i)(i2:i3-1)
                  goto 110
                endif
  370         continue                
            endif
  360     continue                
        endif                   
        if (inline(i)(1:12).eq.'hbtransfile ') then
          do 380 i2=13,80
            ch=inline(i)(i2:i2)
            if (ch.ne.' ') then
              read(inline(i)(i2:80),*,err=600) iz,ia
              Zix=Zinit-iz
              Nix=Ninit-ia+iz
              if (Zix.lt.0.or.Zix.gt.numZ.or.Nix.lt.0.or.Nix.gt.numN)
     +          then
                goto 640
              else
                goto 390
              endif
            endif
  380     continue                
  390     do 400 i2=14,80
            ch=inline(i)(i2:i2)
            if (ch.ne.' '.and.(ch.lt.'0'.or.ch.gt.'9')) then
              do 410 i3=i2,80
                ch=inline(i)(i3:i3)
                if (ch.eq.' ') then
                  hbtransfile(Zix,Nix)=inline(i)(i2:i3-1)
                  goto 110
                endif
  410         continue                
            endif
  400     continue                
        endif                   
        if (inline(i)(1:11).eq.'class2file ') then
          do 420 i2=12,80
            ch=inline(i)(i2:i2)
            if (ch.ne.' ') then
              read(inline(i)(i2:80),*,err=600) iz,ia
              Zix=Zinit-iz
              Nix=Ninit-ia+iz
              if (Zix.lt.0.or.Zix.gt.numZ.or.Nix.lt.0.or.Nix.gt.numN)
     +          then
                goto 640
              else
                goto 430
              endif
            endif
  420     continue                
  430     do 440 i2=13,80
            ch=inline(i)(i2:i2)
            if (ch.ne.' '.and.(ch.lt.'0'.or.ch.gt.'9')) then
              do 450 i3=i2,80
                ch=inline(i)(i3:i3)
                if (ch.eq.' ') then
                  class2file(Zix,Nix)=inline(i)(i2:i3-1)
                  goto 110
                endif
  450         continue                
            endif
  440     continue                
        endif                   
        if (inline(i)(1:12).eq.'optmodfilen ') then
          do 460 i2=13,80
            ch=inline(i)(i2:i2)
            if (ch.ne.' ') then
              read(inline(i)(i2:80),*,err=600) iz
              Zix=Zinit-iz
              if (Zix.lt.0.or.Zix.gt.numZ) then
                goto 640
              else
                goto 470
              endif
            endif
  460     continue                
  470     do 480 i2=14,80
            ch=inline(i)(i2:i2)
            if (ch.ne.' '.and.(ch.lt.'0'.or.ch.gt.'9')) then
              do 490 i3=i2,80
                ch=inline(i)(i3:i3)
                if (ch.eq.' ') then
                  optmodfileN(Zix)=inline(i)(i2:i3-1)
                  goto 110
                endif
  490         continue                
            endif
  480     continue                
        endif                   
        if (inline(i)(1:12).eq.'optmodfilep ') then
          do 500 i3=14,80
            ch=inline(i)(i2:i2)
            if (ch.ne.' ') then
              read(inline(i)(i2:80),*,err=600) iz
              Zix=Zinit-iz
              if (Zix.lt.0.or.Zix.gt.numZ) then
                goto 640
              else
                goto 510
              endif
            endif
  500     continue                
  510     do 520 i2=14,80
            ch=inline(i)(i2:i2)
            if (ch.ne.' '.and.(ch.lt.'0'.or.ch.gt.'9')) then
              do 530 i3=i2,80
                ch=inline(i)(i3:i3)
                if (ch.eq.' ') then
                  optmodfileP(Zix)=inline(i)(i2:i3-1)
                  goto 110
                endif
  530         continue                
            endif
  520     continue                
        endif                   
        if (inline(i)(1:6).eq.'beta2 ') then
          ibar=0
          read(inline(i)(7:80),*,err=600) iz,ia,value
          Zix=Zinit-iz   
          Nix=Ninit-ia+iz
          if (Zix.lt.0.or.Zix.gt.numZ.or.Nix.lt.0.or.Nix.gt.numN) then
            goto 640
          else
            read(inline(i)(7:80),*,end=540,err=600) iz,ia,value,ibar
            if (ibar.lt.0.or.ibar.gt.numbar) goto 610
  540       beta2(Zix,Nix,ibar)=value
          endif
          goto 110
        endif
        if (inline(i)(1:7).eq.'axtype ') then
          ibar=1
          read(inline(i)(8:80),*,err=600) iz,ia,ivalue
          Zix=Zinit-iz   
          Nix=Ninit-ia+iz
          if (Zix.lt.0.or.Zix.gt.numZ.or.Nix.lt.0.or.Nix.gt.numN) then
            goto 640
          else
            read(inline(i)(8:80),*,end=550,err=600) iz,ia,ivalue,ibar
            if (ibar.lt.1.or.ibar.gt.numbar) goto 610
  550       axtype(Zix,Nix,ibar)=ivalue
          endif
          goto 110
        endif
        if (inline(i)(1:6).eq.'gnorm ')
     +    read(inline(i)(7:80),*,err=600) gnorm 
        if (inline(i)(1:9).eq.'rspincut ')
     +    read(inline(i)(10:80),*,err=600) Rspincut 
        if (inline(i)(1:8).eq.'alphald ')
     +    read(inline(i)(9:80),*,err=600) alphald 
        if (inline(i)(1:7).eq.'betald ')
     +    read(inline(i)(8:80),*,err=600) betald 
        if (inline(i)(1:12).eq.'gammashell1 ')
     +    read(inline(i)(13:80),*,err=600) gammashell1
        if (inline(i)(1:12).eq.'gammashell2 ')
     +    read(inline(i)(13:80),*,err=600) gammashell2
        if (inline(i)(1:7).eq.'ufermi ')
     +    read(inline(i)(8:80),*,err=600) Ufermi 
        if (inline(i)(1:7).eq.'cfermi ')
     +    read(inline(i)(8:80),*,err=600) cfermi 
        if (inline(i)(1:9).eq.'ufermibf ')
     +    read(inline(i)(10:80),*,err=600) Ufermibf 
        if (inline(i)(1:9).eq.'cfermibf ')
     +    read(inline(i)(10:80),*,err=600) cfermibf 
        if (inline(i)(1:4).eq.'kph ')
     +    read(inline(i)(5:80),*,err=600) Kph 
        if (inline(i)(1:11).eq.'m2constant ')
     +    read(inline(i)(12:80),*,err=600) M2constant 
        if (inline(i)(1:8).eq.'m2limit ')
     +    read(inline(i)(9:80),*,err=600) M2limit 
        if (inline(i)(1:8).eq.'m2shift ')
     +    read(inline(i)(9:80),*,err=600) M2shift 
        if (inline(i)(1:6).eq.'rpinu ')
     +    read(inline(i)(7:80),*,err=600) Rpinu 
        if (inline(i)(1:6).eq.'esurf ')
     +    read(inline(i)(7:80),*,err=600) Esurf0
        if (inline(i)(1:7).eq.'rgamma ')
     +    read(inline(i)(8:80),*,err=600) Rgamma 
        if (inline(i)(1:8).eq.'msdbins ') 
     +    read(inline(i)(9:80),*,err=600) msdbins          
        if (inline(i)(1:8).eq.'emsdmin ') 
     +    read(inline(i)(9:80),*,err=600) Emsdmin          
        if (inline(i)(1:8).eq.'elwidth ') 
     +    read(inline(i)(9:80),*,err=600) elwidth          
        if (inline(i)(1:10).eq.'anglesrec ') 
     +    read(inline(i)(11:80),*,err=600) nanglerec
        goto 110
  600   write(*,'("TALYS-error: Wrong input: ",a80)') inline(i)
        stop
  610   write(*,'("TALYS-error: 0(1) <= fission barrier <=",i3,$)') 
     +    numbar
        write(*,'(", ibar index out of range: ",a80)') inline(i)
        stop
  620   write(*,'("TALYS-error: 0 <= resonance number <= 2",$)') 
        write(*,'(", igr index out of range: ",a80)') inline(i)
        stop
  630   write(*,'("TALYS-error: 0 <= multipole radiation <= ",i1,$)') 
     +    numgam
        write(*,'(", lval index out of range: ",a80)') inline(i)
        stop
  640   write(*,'("TALYS-warning: Z,N index out of range,",$)')
        write(*,'(" keyword ignored: ",a80)') inline(i)
  110 continue
      return
      end
Copyright (C) 2004  A.J. Koning, S. Hilaire and M.C. Duijvestijn
