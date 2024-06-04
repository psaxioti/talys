      subroutine input2
c
c +---------------------------------------------------------------------
c | Author: Arjan Koning
c | Date  : December 21, 2011
c | Task  : Read input for second set of variables
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      logical      fcol
      character*1  ch
      character*80 word(40),key,value
      integer      type,Zix,Nix,col(0:numZ,0:numN),i,ip,i2,iz,ia,ldmod,
     +             ivalue,type2
c
c ************* Defaults for second set of input variables *************
c
c outtype    : type of outgoing particles
c maxZ,numZ  : maximal number of protons away from the initial
c              compound nucleus
c maxN,numN  : maximal number of neutrons away from the initial
c              compound nucleus
c nbins0     : number of continuum excitation energy bins
c segment    : number of segments to divide emission energy grid
c nlevmax    : maximum number of included discrete levels for target
c nlevmaxres : maximum number of included discrete levels for residual
c              nucleus
c nlevbin    : number of excited levels for binary nucleus
c k0         : index of incident particle
c isomer     : definition of isomer in seconds
c core       : even-even core for weakcoupling (-1 or 1)
c gammax     : number of l-values for gamma multipolarity
c nangle     : number of angles
c numang     : maximum number of angles
c nanglecont : number of angles for continuum
c maxenrec   : number of recoil energies
c massmodel  : model for theoretical nuclear mass
c flagmicro  : flag for completely microscopic Talys calculation
c ldmodelall : level density model for all nuclides
c flagcolall : flag for collective enhancement of level density for all
c              nuclides
c wmode      : designator for width fluctuation model
c preeqmode  : designator for pre-equilibrium model
c mpreeqmode : designator for multiple pre-equilibrium model
c phmodel    : particle-hole state density model
c nlev       : number of excited levels for nucleus
c ldmodel    : level density model
c col        : help variable
c flagcol    : flag for collective enhancement of level density 
c flagomponly: flag to execute ONLY an optical model calculation
c
      do 10 type=0,6
        outtype(type)=' '
   10 continue
      maxZ=numZ-2
      maxN=numN-2
      nbins0=40
      segment=1
      nlevmax=20
      nlevmaxres=10
      do 20 type=0,6
        if (type.le.2.or.type.eq.6) then
          nlevbin(type)=10
        else
          nlevbin(type)=5
        endif
   20 continue
      nlevbin(k0)=nlevmax
      isomer=1.
      core=-1
      gammax=2
      nangle=numang
      nanglecont=18
      maxenrec=numenrec
      massmodel=2
      phmodel=1
      if (flagmicro) then
        ldmodelall=5
      else
        ldmodelall=1
      endif
      flagfission=.false.
      flagcolall=.false.
      if (Atarget.gt.209) then
        flagfission=.true.
        flagcolall=.true.
      endif
      preeqmode=2
      wmode=1
      mpreeqmode=2
      do 30 Nix=0,numN
        do 30 Zix=0,numZ
          nlev(Zix,Nix)=0
          ldmodel(Zix,Nix)=0
          col(Zix,Nix)=0
          flagcol(Zix,Nix)=flagcolall
   30 continue
      flagomponly=.false.
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
      do 110 i=1,nlines
        call getkeywords(inline(i),word)
        key=word(1)
        value=word(2)
        ch=word(2)(1:1)
c
c Test for keywords
c
c parsym: symbol of particle
c Zinit : charge number of initial compound nucleus
c Ninit : neutron number of initial compound nucleus
c
        if (key.eq.'ejectiles') then
          ip=-1
          do 210 i2=2,40
            ch=word(i2)(1:1)
            do 220 type=0,6
              if (ch.eq.parsym(type)) then
                ip=ip+1
                if (ip.le.6) outtype(ip)=ch
                goto 210
              endif
  220       continue
            if (ip.eq.-1) goto 300
  210     continue
          goto 110
        endif
        if (key.eq.'maxz') then
          read(value,*,end=300,err=300) maxZ
          goto 110
        endif
        if (key.eq.'maxn') then
          read(value,*,end=300,err=300) maxN
          goto 110
        endif
        if (key.eq.'bins') then
          read(value,*,end=300,err=300) nbins0
          goto 110
        endif
        if (key.eq.'segment') then
          read(value,*,end=300,err=300) segment
          goto 110
        endif
        if (key.eq.'maxlevelstar') then
          read(value,*,end=300,err=300) nlevmax
          goto 110
        endif
        if (key.eq.'maxlevelsres') then
          read(value,*,end=300,err=300) nlevmaxres
          goto 110
        endif
        if (key.eq.'isomer') then
          read(value,*,end=300,err=300) isomer
          goto 110
        endif
        if (key.eq.'core') then
          read(value,*,end=300,err=300) core
          goto 110
        endif
        if (key.eq.'gammax') then
          read(value,*,end=300,err=300) gammax
          goto 110
        endif
        if (key.eq.'angles') then
          read(value,*,end=300,err=300) nangle
          goto 110
        endif
        if (key.eq.'anglescont') then
          read(value,*,end=300,err=300) nanglecont
          goto 110
        endif
        if (key.eq.'maxenrec') then
          read(value,*,end=300,err=300) maxenrec
          goto 110
        endif
        if (key.eq.'massmodel') then
          read(value,*,end=300,err=300) massmodel
          goto 110
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
              write(*,'(" TALYS-warning: Z,N index out of range,",
     +          " keyword ignored: ",a80)') inline(i)
            else
              ldmodel(Zix,Nix)=ldmod
            endif
          endif
          goto 110
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
              write(*,'(" TALYS-warning: Z,N index out of range,",
     +          " keyword ignored: ",a80)') inline(i)
            else
              flagcol(Zix,Nix)=fcol
              col(Zix,Nix)=1
            endif
          endif
          goto 110
        endif
        if (key.eq.'widthmode') then
          read(value,*,end=300,err=300) wmode
          goto 110
        endif
        if (key.eq.'preeqmode') then
          read(value,*,end=300,err=300) preeqmode
          goto 110
        endif
        if (key.eq.'mpreeqmode') then
          read(value,*,end=300,err=300) mpreeqmode
          goto 110
        endif
        if (key.eq.'phmodel') then
          read(value,*,end=300,err=300) phmodel
          goto 110
        endif
        if (key.eq.'nlevels') then
          read(word(2),*,end=300,err=300) iz
          read(word(3),*,end=300,err=300) ia
          read(word(4),*,end=300,err=300) ivalue
          Zix=Zinit-iz
          Nix=Ninit-ia+iz
          if (Zix.lt.0.or.Zix.gt.numZ.or.Nix.lt.0.or.Nix.gt.numN) then
            write(*,'(" TALYS-warning: Z,N index out of range,",
     +        " keyword ignored: ",a80)') inline(i)
          else
            nlev(Zix,Nix)=ivalue
          endif
          goto 110
        endif
        if (key.eq.'maxlevelsbin') then
          do 250 type=0,6
            if (ch.eq.parsym(type)) then
              type2=type
              goto 260
            endif
  250     continue
          goto 300
  260     read(word(3),*,end=300,err=300) nlevbin(type2)
          goto 110
        endif
        if (key.eq.'omponly') then
          if (ch.eq.'n') flagomponly=.false.
          if (ch.eq.'y') flagomponly=.true.
          if (ch.ne.'y'.and.ch.ne.'n') goto 300
          goto 110
        endif
  110 continue
c
c Set level density models per nucleus
c
      do 310 Nix=0,numN
        do 310 Zix=0,numZ
          if (ldmodel(Zix,Nix).eq.0) ldmodel(Zix,Nix)=ldmodelall
          if (col(Zix,Nix).eq.0) flagcol(Zix,Nix)=flagcolall
          if (ldmodel(Zix,Nix).ge.4) flagcol(Zix,Nix)=.false.
  310 continue
      return
  300 write(*,'(" TALYS-error: Wrong input: ",a80)') inline(i)
      stop
      end
Copyright (C) 2004  A.J. Koning, S. Hilaire and M.C. Duijvestijn
