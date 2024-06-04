      subroutine input6
c
c +---------------------------------------------------------------------
c | Author: Arjan Koning 
c | Date  : June 28, 2004
c | Task  : Read input for sixth set of variables
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      character*1 ch
      integer     type,i,i2,ivalue
      real        value
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
c filerecoil  : flag for recoil spectra on separate file 
c flagfission : flag for fission
c k0          : index of incident particle 
c flagchannels: flag for exclusive channels calculation
c flagrecoil  : flag for calculation of recoils 
c flaggamdis  : flag for output of discrete gamma-ray intensities
c flagdisc    : flag for output of discrete state cross sections 
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
        if (k0.eq.1) filechannels=.true.
        if (flagfission) filefission=.true.
        if (flagrecoil) filerecoil=.true.
        if (k0.eq.1) filegamdis=.true.
        do 40 type=0,6
          filespectrum(type)=.true.
   40   continue
        if (k0.eq.1) then
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
c
c ***************** Read sixth set of input variables ******************
c
c nlines: number of input lines
c ch    : character
c inline: input line                 
c parsym: symbol of particle
c
c First the keyword is identified. Next the corresponding value is read.
c Erroneous input is immediately checked.
c  
      do 110 i=1,nlines
        if (inline(i)(1:13).eq.'filespectrum ') then
          do 120 i2=14,80
            ch=inline(i)(i2:i2)
            if (ch.ne.' ') then 
              do 130 type=0,6
                if (ch.eq.parsym(type)) then
                  filespectrum(type)=.true.
                  goto 120
                endif
  130         continue
              goto 300
            endif
  120     continue
          goto 110
        endif
        if (inline(i)(1:12).eq.'fileelastic ') then
          do 140 i2=13,80
            ch=inline(i)(i2:i2)
            if (ch.ne.' ') then
              if (ch.eq.'n') fileelastic=.false.
              if (ch.eq.'y') fileelastic=.true.
              if (ch.ne.'y'.and.ch.ne.'n') goto 300
              goto 110
            endif
  140     continue
        endif
        if (inline(i)(1:10).eq.'filetotal ') then
          do 150 i2=11,80
            ch=inline(i)(i2:i2)
            if (ch.ne.' ') then
              if (ch.eq.'n') filetotal=.false.
              if (ch.eq.'y') filetotal=.true.
              if (ch.ne.'y'.and.ch.ne.'n') goto 300
              goto 110
            endif
  150     continue
        endif
        if (inline(i)(1:13).eq.'fileresidual ') then
          do 160 i2=14,80
            ch=inline(i)(i2:i2)
            if (ch.ne.' ') then
              if (ch.eq.'n') fileresidual=.false.
              if (ch.eq.'y') fileresidual=.true.
              if (ch.ne.'y'.and.ch.ne.'n') goto 300
              goto 110
            endif
  160     continue
        endif
        if (inline(i)(1:13).eq.'filechannels ') then
          do 170 i2=14,80
            ch=inline(i)(i2:i2)
            if (ch.ne.' ') then
              if (ch.eq.'n') filechannels=.false.
              if (ch.eq.'y') filechannels=.true.
              if (ch.ne.'y'.and.ch.ne.'n') goto 300
              goto 110
            endif
  170     continue
        endif
        if (inline(i)(1:11).eq.'filerecoil ') then
          do 180 i2=12,80
            ch=inline(i)(i2:i2)
            if (ch.ne.' ') then
              if (ch.eq.'n') filerecoil=.false.
              if (ch.eq.'y') filerecoil=.true.
              if (ch.ne.'y'.and.ch.ne.'n') goto 300
              goto 110
            endif
  180     continue
        endif
        if (inline(i)(1:12).eq.'filefission ') then
          do 190 i2=13,80
            ch=inline(i)(i2:i2)
            if (ch.ne.' ') then
              if (ch.eq.'n') filefission=.false.
              if (ch.eq.'y') filefission=.true.
              if (ch.ne.'y'.and.ch.ne.'n') goto 300
              goto 110
            endif
  190     continue
        endif
        if (inline(i)(1:9).eq.'fileddxe ') then
          do 200 i2=10,80
            ch=inline(i)(i2:i2)
            if (ch.ne.' ') then
              do 210 type=0,6
                if (ch.eq.parsym(type)) then
                  ddxecount(type)=ddxecount(type)+1
                  if (ddxecount(type).gt.numfile) goto 330
                  goto 220
                endif
  210         continue           
              goto 300
  220         read(inline(i)(i2+1:80),*,err=300,end=300) value
              fileddxe(type,ddxecount(type))=value
              goto 110
            endif
  200     continue
        endif
        if (inline(i)(1:9).eq.'fileddxa ') then
          do 230 i2=10,80
            ch=inline(i)(i2:i2)
            if (ch.ne.' ') then
              do 240 type=0,6
                if (ch.eq.parsym(type)) then
                  ddxacount(type)=ddxacount(type)+1
                  if (ddxacount(type).gt.numfile) goto 340
                  goto 250
                endif
  240         continue           
              goto 300
  250         read(inline(i)(i2+1:80),*,err=300,end=300) value
              fileddxa(type,ddxacount(type))=value
              goto 110
            endif
  230     continue
        endif
        if (inline(i)(1:10).eq.'fileangle ') then
          read(inline(i)(11:80),*,err=300) ivalue
          if (ivalue.lt.0.or.ivalue.gt.numlev) goto 310
          fileangle(ivalue)=.true.
          goto 110
        endif
        if (inline(i)(1:13).eq.'filediscrete ') then
          read(inline(i)(14:80),*,err=300) ivalue
          if (ivalue.lt.0.or.ivalue.gt.numlev) goto 320
          filediscrete(ivalue)=.true.
          goto 110
        endif
        if (inline(i)(1:11).eq.'filegamdis ') then
          do 260 i2=12,80
            ch=inline(i)(i2:i2)
            if (ch.ne.' ') then
              if (ch.eq.'n') filegamdis=.false.
              if (ch.eq.'y') filegamdis=.true.
              if (ch.ne.'y'.and.ch.ne.'n') goto 300
              goto 110
            endif
  260     continue
        endif
  110 continue
      return
  300 write(*,'("TALYS-error: Wrong input: ",a80)') inline(i)
      stop                                                     
  310 write(*,'("TALYS-error: 0 <= fileangle <=",i3,$)') numlev
      write(*,'(", fileangle index out of range: ",a80)') inline(i)
      stop         
  320 write(*,'("TALYS-error: 0 <= filediscrete <=",i3,$)') numlev
      write(*,'(", filediscrete index out of range: ",a80)') inline(i)
      stop         
  330 write(*,'("TALYS-error: number of fileddxe <=",i3,$)') numfile
      write(*,'(", index out of range: ",a80)') inline(i)
      stop         
  340 write(*,'("TALYS-error: number of fileddxa <=",i3,$)') numfile
      write(*,'(", index out of range: ",a80)') inline(i)
      stop         
      end
Copyright (C) 2004  A.J. Koning, S. Hilaire and M.C. Duijvestijn
