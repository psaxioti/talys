      subroutine input6
c
c +---------------------------------------------------------------------
c | Author: Arjan Koning 
c | Date  : May 11, 2007
c | Task  : Read input for sixth set of variables
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      character*1  ch
      character*80 word(40),key,value
      integer      type,i,i2,ivalue
      real         val
c
c ************* Defaults for sixth set of input variables **************
c
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
c filechannels: flag for exclusive channel cross sections on
c               separate file      
c filefission : flag for fission cross sections on separate file
c filegamdis  : flag for gamma-ray intensities on separate file 
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
c flagastro   : flag for calculation of astrophysics reaction rate
c
      fileelastic=.false.
      do 10 type=0,6
        filespectrum(type)=.false.
        ddxecount(type)=0
        ddxacount(type)=0
        do 20 i=1,numfile
          fileddxe(type,i)=0.
          fileddxa(type,i)=0.
   20   continue
   10 continue
      do 30 i=0,numlev
        fileangle(i)=.false.
        filediscrete(i)=.false.
   30 continue
      filetotal=.false.
      fileresidual=.false.
      filechannels=.false.
      filerecoil=.false.
      filefission=.false.
      filegamdis=.false.
c
c If the results of TALYS are used to create ENDF-6 data files,
c several output flags are automatically set.
c 
      if (flagendf) then
        fileelastic=.true.
        filetotal=.true.
        fileresidual=.true.
        if (flagendfdet) filechannels=.true.
        if (flagfission) filefission=.true.
        if (flagrecoil) filerecoil=.true.
        if (flagendfdet) filegamdis=.true.
        do 40 type=0,6
          filespectrum(type)=.true.
   40   continue
        if (flagendfdet) then
          do 50 i=0,numlev
            fileangle(i)=.true.
            filediscrete(i)=.true.
   50     continue
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
        if (flagdisc) then
          do 60 i=0,numlev
            filediscrete(i)=.true.
   60     continue
        endif
      endif
      filedensity=.true.
      if (flagastro) fileresidual=.true.
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
      do 110 i=1,nlines
        call getkeywords(inline(i),word)
        key=word(1)
        value=word(2)
        ch=word(2)(1:1)
c
c Test for keywords
c  
c parsym: symbol of particle
c  
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
          goto 110
        endif
        if (key.eq.'fileelastic') then
          if (ch.eq.'n') fileelastic=.false.
          if (ch.eq.'y') fileelastic=.true.
          if (ch.ne.'y'.and.ch.ne.'n') goto 300
          goto 110
        endif
        if (key.eq.'filetotal') then
          if (ch.eq.'n') filetotal=.false.
          if (ch.eq.'y') filetotal=.true.
          if (ch.ne.'y'.and.ch.ne.'n') goto 300
          goto 110
        endif
        if (key.eq.'fileresidual') then
          if (ch.eq.'n') fileresidual=.false.
          if (ch.eq.'y') fileresidual=.true.
          if (ch.ne.'y'.and.ch.ne.'n') goto 300
          goto 110
        endif
        if (key.eq.'filechannels') then
          if (ch.eq.'n') filechannels=.false.
          if (ch.eq.'y') filechannels=.true.
          if (ch.ne.'y'.and.ch.ne.'n') goto 300
          goto 110
        endif
        if (key.eq.'filerecoil') then
          if (ch.eq.'n') filerecoil=.false.
          if (ch.eq.'y') filerecoil=.true.
          if (ch.ne.'y'.and.ch.ne.'n') goto 300
          goto 110
        endif
        if (key.eq.'filefission') then
          if (ch.eq.'n') filefission=.false.
          if (ch.eq.'y') filefission=.true.
          if (ch.ne.'y'.and.ch.ne.'n') goto 300
          goto 110
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
          goto 110
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
          goto 110
        endif
        if (key.eq.'fileangle') then
          read(value,*,end=300,err=300) ivalue
          if (ivalue.lt.0.or.ivalue.gt.numlev) goto 310
          fileangle(ivalue)=.true.
          goto 110
        endif
        if (key.eq.'filediscrete') then
          read(value,*,end=300,err=300) ivalue
          if (ivalue.lt.0.or.ivalue.gt.numlev) goto 320
          filediscrete(ivalue)=.true.
          goto 110
        endif
        if (key.eq.'filegamdis') then
          if (ch.eq.'n') filegamdis=.false.
          if (ch.eq.'y') filegamdis=.true.
          if (ch.ne.'y'.and.ch.ne.'n') goto 300
          goto 110
        endif
        if (key.eq.'filedensity') then
          if (ch.eq.'n') filedensity=.false.
          if (ch.eq.'y') filedensity=.true.
          if (ch.ne.'y'.and.ch.ne.'n') goto 300
          goto 110
        endif
  110 continue
      return
  300 write(*,'(" TALYS-error: Wrong input: ",a80)') inline(i)
      stop                                                     
  310 write(*,'(" TALYS-error: 0 <= fileangle <=",i3, 
     +  ", fileangle index out of range: ",a80)') numlev,inline(i)
      stop         
  320 write(*,'(" TALYS-error: 0 <= filediscrete <=",i3,
     +  ", filediscrete index out of range: ",a80)') numlev,inline(i)
      stop         
  330 write(*,'(" TALYS-error: number of fileddxe <=",i3,
     +  ", index out of range: ",a80)') numfile,inline(i)
      stop         
  340 write(*,'(" TALYS-error: number of fileddxa <=",i3,
     +  ", index out of range: ",a80)') numfile,inline(i)
      stop         
      end
Copyright (C) 2004  A.J. Koning, S. Hilaire and M.C. Duijvestijn
