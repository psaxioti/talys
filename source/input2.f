      subroutine input2
c
c +---------------------------------------------------------------------
c | Author: Arjan Koning
c | Date  : January 24, 2023
c | Task  : Read input for second set of variables
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      logical      fcol,flagassign
      character*1  ch
      character*132 word(40),key,value,cval,line
      integer      type,Zix,Nix,col(0:numZ,0:numN),i,ip,i2,iz,ia,ldmod,
     +             nex,class,ival,ibar,irad,lval,igr
      real         val,sfthall,sfexpall,sfth,sfexp
c
c ************* Defaults for second set of input variables *************
c
c outtype     : type of outgoing particles
c maxZ,numZ   : maximal number of protons away from the initial
c               compound nucleus
c maxN,numN   : maximal number of neutrons away from the initial
c               compound nucleus
c nbins0      : number of continuum excitation energy bins
c segment     : number of segments to divide emission energy grid
c nlevmax     : maximum number of included discrete levels for target
c Ltarget     : excited level of target
c nlevmaxres  : maximum number of included discrete levels for residual
c               nucleus
c nlevbin     : number of excited levels for binary nucleus
c k0          : index of incident particle
c Lisoinp     : user assignment of target isomer number
c isomer      : definition of isomer in seconds
c core        : even-even core for weakcoupling (-1 or 1)
c gammax      : number of l-values for gamma multipolarity
c nangle      : number of angles
c numang      : maximum number of angles
c nanglecont  : number of angles for continuum
c maxenrec    : number of recoil energies
c massmodel   : model for theoretical nuclear mass
c disctable   : table with discrete levels
c flagmicro   : flag for completely microscopic Talys calculation
c ldmodelall  : level density model for all nuclides
c ldmodelCN   : level density model for compound nucleus
c flagcolall  : flag for collective enhancement of level density for all
c               nuclides
c fislim      : mass above which nuclide fissions
c wmode       : designator for width fluctuation model
c WFCfactor   : enhancement factor for WFC: 1: Moldauer, 
c               2: Ernebjerg and Herman
c preeqmode   : designator for pre-equilibrium model
c mpreeqmode  : designator for multiple pre-equilibrium model
c phmodel     : particle-hole state density model
c nlev        : number of excited levels for nucleus
c ldmodel     : level density model
c skipCN      : flag to skip compound nucleus in evaporation chain
c col         : help variable
c flagcol     : flag for collective enhancement of level density
c numlev      : maximum number of included discrete levels
c strength    : model for E1 gamma-ray strength function
c flagomponly : flag to execute ONLY an optical model calculation
c flagjlm     : flag for using semi-microscopic JLM OMP
c flagequi    : flag to use equidistant excitation bins instead of
c               logarithmic bins
c flagequispec: flag to use equidistant bins for emission spectra
c flagpopMeV  : flag to use initial population per MeV instead of
c               histograms
c flagmassdis : flag for calculation of fission fragment mass yields
c flagracap   : flag for radiative capture model
c ldmodelracap: level density model for direct radiative capture
c spectfacexp : experimental spectroscopic factor
c spectfacth  : theoretical spectroscopic factor
c maxZrp      : maximal number of protons away from the initial
c               compound nucleus before new residual evaporation
c maxNrp      : maximal number of neutrons away from the initial
c               compound nucleus before new residual evaporation
c
      outtype=' '
      maxZ=numZ-2
      maxN=numN-2
      nbins0=40
      flagequispec=.false.
      segment=1
      nlevmax=max(30,Ltarget)
      nlevmaxres=10
      do type=0,6
        if (type.le.2.or.type.eq.6) then
          nlevbin(type)=10
        else
          nlevbin(type)=5
        endif
      enddo
      nlevbin(k0)=nlevmax
      Lisoinp=-1
      isomer=1.
      core=-1
      gammax=2
      nangle=numang
      nanglecont=18
      maxenrec=numenrec
      massmodel=2
      disctable=1
      phmodel=1
      if (flagmicro) then
        ldmodelall=5
        strength=8
        flagjlm=.true.
      else
        ldmodelall=1
        strength=9
        flagjlm=.false.
      endif
      if (k0.le.1 .and. Atarget.gt.fislim) ldmodelall = 5
      ldmodelCN=0
      if (Atarget.gt.fislim) then
        flagcolall=.true.
      else
        flagcolall=.false.
      endif
      preeqmode=2
      if (k0.eq.1) then
        wmode=1
      else
        wmode=2
      endif
      WFCfactor=1
      mpreeqmode=2
      sfthall=1.
      sfexpall=0.347
      if (mod(Atarget,2).ne.0) sfexpall=1.
      nlev=0
      ldmodel=0
      skipCN=0
      col=0
      flagcol=flagcolall
      spectfacth=0.
      spectfacexp=0.
      flagomponly=.false.
      flagequi=.true.
      flagpopMeV=.false.
      flagmassdis=.false.
      flagracap=.false.
      ldmodelracap=3
      maxZrp=numZ-2
      maxNrp=numN-2
c
c **************** Read second set of input variables ******************
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
c
c Test for keywords
c
c parsym   : symbol of particle
c Zinit    : charge number of initial compound nucleus
c Ninit    : neutron number of initial compound nucleus
c getvalues: subroutine to assign values to keywords
c fcol     : flag for collective enhancement
c sfexp    : variable for spectrocopic factor
c sfexpall : variable for spectrocopic factor
c sfth     : variable for spectrocopic factor
c sfthall  : variable for spectrocopic factor
c
        if (key.eq.'ejectiles') then
          ip=-1
          do 210 i2=2,40
            ch=word(i2)(1:1)
            do type=0,6
              if (ch.eq.parsym(type)) then
                ip=ip+1
                if (ip.le.6) outtype(ip)=ch
                goto 210
              endif
            enddo
            if (ip.eq.-1) goto 300
  210     continue
          cycle
        endif
        if (key.eq.'maxz') then
          read(value,*,end=300,err=300) maxZ
          cycle
        endif
        if (key.eq.'maxn') then
          read(value,*,end=300,err=300) maxN
          cycle
        endif
        if (key.eq.'bins') then
          read(value,*,end=300,err=300) nbins0
          cycle
        endif
        if (key.eq.'segment') then
          read(value,*,end=300,err=300) segment
          cycle
        endif
        if (key.eq.'maxlevelstar') then
          read(value,*,end=300,err=300) nlevmax
          nlevmax=max(nlevmax,Ltarget)
          nlevbin(k0)=nlevmax
          cycle
        endif
        if (key.eq.'maxlevelsres') then
          read(value,*,end=300,err=300) nlevmaxres
          cycle
        endif
        if (key.eq.'liso') then
          read(value,*,end=300,err=300) Lisoinp
          cycle
        endif
        if (key.eq.'isomer') then
          read(value,*,end=300,err=300) isomer
          cycle
        endif
        if (key.eq.'core') then
          read(value,*,end=300,err=300) core
          cycle
        endif
        if (key.eq.'gammax') then
          read(value,*,end=300,err=300) gammax
          cycle
        endif
        if (key.eq.'angles') then
          read(value,*,end=300,err=300) nangle
          cycle
        endif
        if (key.eq.'anglescont') then
          read(value,*,end=300,err=300) nanglecont
          cycle
        endif
        if (key.eq.'maxenrec') then
          read(value,*,end=300,err=300) maxenrec
          cycle
        endif
        if (key.eq.'massmodel') then
          read(value,*,end=300,err=300) massmodel
          cycle
        endif
        if (key.eq.'disctable') then
          read(value,*,end=300,err=300) disctable
          cycle
        endif
        if (key.eq.'ldmodel') then
          read(value,*,end=300,err=300) ldmod
          read(word(3),*,end=230,err=300) iz
          read(word(4),*,end=300,err=300) ia
  230     if (word(3).eq.' ') then
            ldmodelall=ldmod
          else
            Zix=Zinit-iz
            Nix=Ninit-ia+iz
            if (Zix.lt.0.or.Zix.gt.numZ.or.Nix.lt.0.or.Nix.gt.numN) then
              goto 1000
            else
              ldmodel(Zix,Nix)=ldmod
            endif
          endif
          cycle
        endif
        if (key.eq.'ldmodelcn') then
          read(value,*,end=300,err=300) ldmodelCN
          cycle
        endif
        if (key.eq.'colenhance') then
          if (ch.eq.'n') fcol=.false.
          if (ch.eq.'y') fcol=.true.
          if (ch.ne.'y'.and.ch.ne.'n') goto 300
          read(word(3),*,end=240,err=300) iz
          read(word(4),*,end=300,err=300) ia
  240     if (word(3).eq.' ') then
            flagcolall=fcol
          else
            Zix=Zinit-iz
            Nix=Ninit-ia+iz
            if (Zix.lt.0.or.Zix.gt.numZ.or.Nix.lt.0.or.Nix.gt.numN) then
              goto 1000
            else
              flagcol(Zix,Nix)=fcol
              col(Zix,Nix)=1
            endif
          endif
          cycle
        endif
        if (key.eq.'skipcn') then
          read(word(2),*,end=300,err=300) iz
          read(word(3),*,end=300,err=300) ia
          Zix=Zinit-iz
          Nix=Ninit-ia+iz
          if (Zix.lt.0.or.Zix.gt.numZ.or.Nix.lt.0.or.Nix.gt.numN) then
            goto 1000
          else
            skipCN(Zix,Nix)=1
          endif
          cycle
        endif
        if (key.eq.'widthmode') then
          read(value,*,end=300,err=300) wmode
          cycle
        endif
        if (key.eq.'wfcfactor') then
          read(value,*,end=300,err=300) WFCfactor
          cycle
        endif
        if (key.eq.'preeqmode') then
          read(value,*,end=300,err=300) preeqmode
          cycle
        endif
        if (key.eq.'mpreeqmode') then
          read(value,*,end=300,err=300) mpreeqmode
          cycle
        endif
        if (key.eq.'phmodel') then
          read(value,*,end=300,err=300) phmodel
          cycle
        endif
        if (key.eq.'nlevels') then
          class=2
          call getvalues(class,word,Zix,Nix,type,
     +      ibar,irad,lval,igr,val,ival,cval,flagassign)
          if (flagassign) nlev(Zix,Nix)=ival
          cycle
        endif
        if (key.eq.'maxlevelsbin') then
          class=7
          call getvalues(class,word,Zix,Nix,type,
     +      ibar,irad,lval,igr,val,ival,cval,flagassign)
          if (flagassign) nlevbin(type)=ival
          cycle
        endif
        if (key.eq.'strength') then
          read(value,*,end=300,err=300) strength
          cycle
        endif
        if (key.eq.'omponly') then
          if (ch.eq.'n') flagomponly=.false.
          if (ch.eq.'y') flagomponly=.true.
          if (ch.ne.'y'.and.ch.ne.'n') goto 300
          cycle
        endif
        if (key.eq.'jlmomp') then
          if (ch.eq.'n') flagjlm=.false.
          if (ch.eq.'y') flagjlm=.true.
          if (ch.ne.'y'.and.ch.ne.'n') goto 300
          cycle
        endif
        if (key.eq.'equidistant') then
          if (ch.eq.'n') flagequi=.false.
          if (ch.eq.'y') flagequi=.true.
          if (ch.ne.'y'.and.ch.ne.'n') goto 300
          cycle
        endif
        if (key.eq.'equispec') then
          if (ch.eq.'n') flagequispec=.false.
          if (ch.eq.'y') flagequispec=.true.
          if (ch.ne.'y'.and.ch.ne.'n') goto 300
          cycle
        endif
        if (key.eq.'popmev') then
          if (ch.eq.'n') flagpopMeV=.false.
          if (ch.eq.'y') flagpopMeV=.true.
          if (ch.ne.'y'.and.ch.ne.'n') goto 300
          cycle
        endif
        if (key.eq.'massdis') then
          if (ch.eq.'n') flagmassdis=.false.
          if (ch.eq.'y') flagmassdis=.true.
          if (ch.ne.'y'.and.ch.ne.'n') goto 300
          cycle
        endif
        if (key.eq.'racap') then
          if (ch.eq.'n') flagracap=.false.
          if (ch.eq.'y') flagracap=.true.
          if (ch.ne.'y'.and.ch.ne.'n') goto 300
          cycle
        endif
        if (key.eq.'ldmodelracap') then
          read(value,*,end=300,err=300) ldmodelracap
          cycle
        endif
        if (key.eq.'maxzrp') then
          read(value,*,end=300,err=300) maxZrp
          cycle
        endif
        if (key.eq.'maxnrp') then
          read(value,*,end=300,err=300) maxNrp
          cycle
        endif
        if (key.eq.'sfth') then
          read(value,*,end=300,err=300) sfth
          read(word(3),*,end=250,err=300) iz
          read(word(4),*,end=300,err=300) ia
  250     if (word(3).eq.' ') then
            sfthall=sfth
          else
            Zix=Zinit-iz
            Nix=Ninit-ia+iz
            if (Zix.lt.0.or.Zix.gt.numZ.or.Nix.lt.0.or.Nix.gt.numN) then
              goto 1000
            else
              spectfacth(Zix,Nix)=sfth
            endif
          endif
          cycle
        endif
        if (key.eq.'sfexp') then
          nex=-1
          Zix=0
          Nix=0
          read(value,*,end=300,err=300) sfexp
          read(word(3),*,end=260,err=300) iz
          read(word(4),*,end=260,err=300) ia
          read(word(5),*,end=260,err=300) nex
  260     if (word(3).eq.' ') then
            sfexpall=sfexp
          else
            Zix=Zinit-iz
            Nix=Ninit-ia+iz
            if (Zix.lt.0.or.Zix.gt.numZ.or.Nix.lt.0.or.Nix.gt.numN)
     +        goto 1000
            if (word(3).eq.' '.or.word(4).eq.' ') nex=iz
            if (nex.eq.-1) then
              sfexpall=sfexp
            else
              if (nex.lt.0.or.nex.gt.numlev) then
                write(*,'(" TALYS-error: 0 <= nex <= numlev ",a)')
     +            trim(line)
                stop
              else
                spectfacexp(Zix,Nix,nex)=sfexp
              endif
            endif
          endif
          cycle
        endif
        cycle
 1000   write(*,'(" TALYS-warning: Z,N index out of range,",
     +    " keyword ignored: ",a)') trim(line)
      enddo
c
c Set level density models and spectroscopic factors per nucleus
c
      if (ldmodelCN.gt.0) then
        ldmodel(0,0)=ldmodelCN
      else
        ldmodelCN=ldmodelall
      endif
      do Nix=0,numN
        do Zix=0,numZ
          if (ldmodel(Zix,Nix).eq.0) ldmodel(Zix,Nix)=ldmodelall
          if (col(Zix,Nix).eq.0) flagcol(Zix,Nix)=flagcolall
          if (spectfacth(Zix,Nix).eq.0.) spectfacth(Zix,Nix)=sfthall
          do nex=0,numlev
            if (spectfacexp(Zix,Nix,nex).eq.0.)
     +        spectfacexp(Zix,Nix,nex)=sfexpall
          enddo
        enddo
      enddo
      return
  300 write(*,'(" TALYS-error: Wrong input: ",a)') trim(line)
      stop
      end
Copyright (C)  2023 A.J. Koning, S. Hilaire and S. Goriely
