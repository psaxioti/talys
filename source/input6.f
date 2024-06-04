      subroutine input6
c
c +---------------------------------------------------------------------
c | Author: Arjan Koning
c | Date  : January 29, 2023
c | Task  : Read input for sixth set of variables
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      character*1   ch
      character*132 word(40),key,value,line
      integer       type,i,i2,ivalue,istat
      real          val
c
c ************* Defaults for sixth set of input variables **************
c
c flagastro   : flag for calculation of astrophysics reaction rate
c transpower  : power for transmission coefficient limit
c transeps    : absolute limit for transmission coefficient
c xseps       : limit for cross sections
c popeps      : limit for population cross section per nucleus
c Rfiseps     : ratio for limit for fission cross section per nucleus
c fileelastic : flag for elastic angular distribution on separate file
c filespectrum: designator for spectrum on separate file
c ddxecount   : counter for double-differential cross section files
c ddxacount   : counter for double-differential cross section files
c numfile     : maximum number of separate output files
c fileddxe    : designator for double-differential cross sections on
c               separate file: angular distribution
c fileddxa    : designator for double-differential cross sections on
c               separate file: spectrum per angle
c numlev      : maximum number of included discrete levels
c fileangle   : designator for angular distributions on separate file
c filediscrete: flag for discrete level cross sections on separate
c               file
c filetotal   : flag for total cross sections on separate file
c fileresidual: flag for residual production cross sections on
c               separate file
c flagcompo   : flag for output of cross section components
c filechannels: flag for exclusive channel cross sections on
c               separate file
c filefission : flag for fission cross sections on separate file
c filegamdis  : flag for gamma-ray intensities on separate file
c flagblock   : flag to block spectra, angle and gamma files
c flagexc     : flag for output of excitation functions
c flagendf    : flag for information for ENDF-6 file
c flagendfdet : flag for detailed ENDF-6 information per channel
c filerecoil  : flag for recoil spectra on separate file
c flagfission : flag for fission
c flagchannels: flag for exclusive channels calculation
c flagrecoil  : flag for calculation of recoils
c flaggamdis  : flag for output of discrete gamma-ray intensities
c flagdisc    : flag for output of discrete state cross sections
c filedensity : flag for level densities on separate files
c filepsf     : flag for photon strength functions on separate files
c flagastro   : flag for calculation of astrophysics reaction rate
c flagintegral: flag for calculation of effective cross section using
c               integral spectrum
c flagsacs    : flag for statistical analysis of cross sections
c Nflux       : number of reactions with integral data
c xsfluxfile  : TALYS cross section file for integral data
c fluxname    : name of integral spectrum
c integralexp : experimental effective cross section
c
      transpower=5
      transeps=1.e-8
      xseps=1.e-7
      popeps=1.e-3
      Rfiseps=1.e-3
      if (flagmassdis) then
        transpower=10
        transeps=1.e-12
        xseps=1.e-12
        popeps=1.e-10
        Rfiseps=1.e-9
      endif
      if (flagastro) then
        transpower=15
        transeps=1.e-18
        xseps=1.e-17
        popeps=1.e-13
        Rfiseps=1.e-6
      endif
      fileelastic=.false.
      filespectrum=.false.
      ddxecount=0
      ddxacount=0
      fileddxe=0.
      fileddxa=0.
c
c Explicit double-differential cross sections for deuteron ENDF files
c
      if (flagendf.and.k0.eq.3) then
        do type=1,2
          filespectrum(type)=.true.
          ddxacount(type)=18
          do i=1,7
            fileddxa(type,i)=5.*(i-1)
          enddo
          do i=8,14
            fileddxa(type,i)=30.+10.*(i-7)
          enddo
          do i=15,18
            fileddxa(type,i)=100.+20.*(i-14)
          enddo
        enddo
      endif
      fileangle=.false.
      filediscrete=.false.
      flagblock=.false.
      filetotal=.false.
      fileresidual=.false.
      flagcompo=.false.
      filechannels=.false.
      filerecoil=.false.
      filefission=.false.
      filegamdis=.false.
c
c If the results of TALYS are used to create ENDF-6 data files,
c several output flags are automatically set.
c
      if (flagendf) then
        flagblock=.true.
        fileelastic=.true.
        filetotal=.true.
        fileresidual=.true.
        if (flagendfdet) filechannels=.true.
        if (flagfission) filefission=.true.
        if (flagrecoil) filerecoil=.true.
        if (flagendfdet) filegamdis=.true.
        filespectrum=.true.
        if (flagendfdet) then
          fileangle=.true.
          filediscrete=.true.
        endif
      endif
c
c If the results of TALYS are written as excitation functions in the
c output file, several output flags are automatically set.
c
      if (flagexc) then
        filetotal=.true.
        fileresidual=.true.
        if (flagchannels) filechannels=.true.
        if (flagfission) filefission=.true.
        if (flaggamdis) filegamdis=.true.
        if (flagdisc) filediscrete=.true.
      endif
      filedensity=.false.
      filepsf=.false.
      if (flagastro) fileresidual=.true.
      flagintegral=.false.
      flagsacs=.false.
      Nflux=0
      xsfluxfile='                                                '
      fluxname='                                                  '
      integralexp=0.
c
c ***************** Read sixth set of input variables ******************
c
c nlines     : number of input lines
c getkeywords: subroutine to retrieve keywords and values from input
c              line
c inline     : input line
c word       : words on input line
c key        : keyword
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
c parsym: symbol of particle
c ivalue: counter
c
       if (key.eq.'transpower') then
          read(value,*,end=300,err=300) transpower
          cycle
        endif
        if (key.eq.'transeps') then
          read(value,*,end=300,err=300) transeps
          cycle
        endif
        if (key.eq.'xseps') then
          read(value,*,end=300,err=300) xseps
          cycle
        endif
        if (key.eq.'popeps') then
          read(value,*,end=300,err=300) popeps
          cycle
        endif
        if (key.eq.'rfiseps') then
          read(value,*,end=300,err=300) Rfiseps
          cycle
        endif
        if (key.eq.'filespectrum') then
          do 210 i2=2,40
            ch=word(i2)(1:1)
            do 220 type=0,6
              if (ch.eq.parsym(type)) then
                filespectrum(type)=.true.
                goto 210
              endif
  220       continue
  210     continue
          cycle
        endif
        if (key.eq.'fileelastic') then
          if (ch.eq.'n') fileelastic=.false.
          if (ch.eq.'y') fileelastic=.true.
          if (ch.ne.'y'.and.ch.ne.'n') goto 300
          cycle
        endif
        if (key.eq.'filetotal') then
          if (ch.eq.'n') filetotal=.false.
          if (ch.eq.'y') filetotal=.true.
          if (ch.ne.'y'.and.ch.ne.'n') goto 300
          cycle
        endif
        if (key.eq.'fileresidual') then
          if (ch.eq.'n') fileresidual=.false.
          if (ch.eq.'y') fileresidual=.true.
          if (ch.ne.'y'.and.ch.ne.'n') goto 300
          cycle
        endif
        if (key.eq.'components') then
          if (ch.eq.'n') flagcompo=.false.
          if (ch.eq.'y') flagcompo=.true.
          if (ch.ne.'y'.and.ch.ne.'n') goto 300
          cycle
        endif
        if (key.eq.'filechannels') then
          if (ch.eq.'n') filechannels=.false.
          if (ch.eq.'y') filechannels=.true.
          if (ch.ne.'y'.and.ch.ne.'n') goto 300
          cycle
        endif
        if (key.eq.'filerecoil') then
          if (ch.eq.'n') filerecoil=.false.
          if (ch.eq.'y') filerecoil=.true.
          if (ch.ne.'y'.and.ch.ne.'n') goto 300
          cycle
        endif
        if (key.eq.'filefission') then
          if (ch.eq.'n') filefission=.false.
          if (ch.eq.'y') filefission=.true.
          if (ch.ne.'y'.and.ch.ne.'n') goto 300
          cycle
        endif
        if (key.eq.'fileddxe') then
          do 230 type=0,6
            if (ch.eq.parsym(type)) then
              ddxecount(type)=ddxecount(type)+1
              if (ddxecount(type).gt.numfile) goto 330
              goto 240
            endif
  230     continue
          goto 300
  240     read(word(3),*,err=300,end=300) val
          fileddxe(type,ddxecount(type))=val
          cycle
        endif
        if (key.eq.'fileddxa') then
          do 250 type=0,6
            if (ch.eq.parsym(type)) then
              ddxacount(type)=ddxacount(type)+1
              if (ddxacount(type).gt.numfile) goto 340
              goto 260
            endif
  250     continue
          goto 300
  260     read(word(3),*,err=300,end=300) val
          fileddxa(type,ddxacount(type))=val
          cycle
        endif
        if (key.eq.'integral') then
          if (ch.eq.'n') flagintegral=.false.
          if (ch.eq.'y') flagintegral=.true.
          if (ch.ne.'y'.and.ch.ne.'n') then
            Nflux=Nflux+1
            if (k0.gt.1) goto 350
            if (Nflux.gt.numflux) goto 360
            xsfluxfile(Nflux)=value
            fluxname(Nflux)=word(3)
            flagintegral=.true.
            read(word(4),*,iostat=istat) integralexp(Nflux)
            if (istat.ne.0) cycle
          endif
          cycle
        endif
        if (key.eq.'sacs') then
          if (ch.eq.'n') flagsacs=.false.
          if (ch.eq.'y') flagsacs=.true.
          if (ch.ne.'y'.and.ch.ne.'n') goto 300
          cycle
        endif
        if (key.eq.'block') then
          if (ch.eq.'n') flagblock=.false.
          if (ch.eq.'y') flagblock=.true.
          if (ch.ne.'y'.and.ch.ne.'n') goto 300
          cycle
        endif
        if (key.eq.'fileangle') then
          read(value,*,end=300,err=300) ivalue
          if (ivalue.lt.0.or.ivalue.gt.numlev) goto 310
          fileangle(ivalue)=.true.
          cycle
        endif
        if (key.eq.'filediscrete') then
          read(value,*,end=300,err=300) ivalue
          if (ivalue.lt.0.or.ivalue.gt.numlev) goto 320
          filediscrete(ivalue)=.true.
          cycle
        endif
        if (key.eq.'filegamdis') then
          if (ch.eq.'n') filegamdis=.false.
          if (ch.eq.'y') filegamdis=.true.
          if (ch.ne.'y'.and.ch.ne.'n') goto 300
          cycle
        endif
        if (key.eq.'filedensity') then
          if (ch.eq.'n') filedensity=.false.
          if (ch.eq.'y') filedensity=.true.
          if (ch.ne.'y'.and.ch.ne.'n') goto 300
          cycle
        endif
        if (key.eq.'filepsf') then
          if (ch.eq.'n') filepsf=.false.
          if (ch.eq.'y') filepsf=.true.
          if (ch.ne.'y'.and.ch.ne.'n') goto 300
          cycle
        endif
      enddo
      return
  300 write(*,'(" TALYS-error: Wrong input: ",a)') trim(line)
      stop
  310 write(*,'(" TALYS-error: 0 <= fileangle <=",i3,
     +  ", fileangle index out of range: ",a)') numlev,trim(line)
      stop
  320 write(*,'(" TALYS-error: 0 <= filediscrete <=",i3,
     +  ", filediscrete index out of range: ",a)') 
     +  numlev,trim(line)
      stop
  330 write(*,'(" TALYS-error: number of fileddxe <=",i3,
     +  ", index out of range: ",a)') numfile,trim(line)
      stop
  340 write(*,'(" TALYS-error: number of fileddxa <=",i3,
     +  ", index out of range: ",a)') numfile,trim(line)
  350 write(*,'(" TALYS-error: effective cross section can only be",
     +  " calculated for incident photons and neutrons")')
  360 write(*,'(" TALYS-error: number of integral data sets <=",i3,
     +  ", index out of range: ",a)') numflux,trim(line)
      stop
      end
Copyright (C)  2023 A.J. Koning, S. Hilaire and S. Goriely
