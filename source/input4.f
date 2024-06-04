      subroutine input4
c
c +---------------------------------------------------------------------
c | Author: Arjan Koning 
c | Date  : September 14, 2004
c | Task  : Read input for fourth set of variables
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      character*1 ch
      integer     i,i2,type
c
c ************** Defaults for fourth set of input variables ************
c
c flagmain    : flag for main output 
c flagbasic   : flag for output of basic information and results
c flagpop     : flag for output of population 
c flagcheck   : flag for output of numerical checks 
c flagoutomp  : flag for output of optical model parameters
c flagdirect  : flag for output of direct reaction cross sections
c flaginverse : flag for output of transmission coefficients and inverse
c               reaction cross sections
c flaggamma   : flag for output of gamma-ray information
c flaglevels  : flag for output of discrete level information
c flagdensity : flag for output of level densities       
c flagdisc    : flag for output of discrete state cross sections 
c flagfisout  : flag for output of fission information   
c flagfission : flag for fission
c flagtransen : flag for output of transmission coefficients per energy
c flagpeout   : flag for output of pre-equilibrium results 
c flagang     : flag for output of angular distributions 
c ddxmode     : mode for double-differential cross sections: 0: None,
c               1: Angular distributions, 2: Spectra per angle, 3: Both
c flaglegendre: flag for output of Legendre coefficients
c flagspec    : flag for output of spectra 
c flagrecoil  : flag for calculation of recoils
c flagddx     : flag for output of double-differential cross sections
c flagoutdwba : flag for output of DWBA cross sections for MSD
c flaggamdis  : flag for output of discrete gamma-ray intensities
c flageciscomp: flag for compound nucleus calculation by ECIS
c flagoutecis : flag for output of ECIS results
c numinc      : number of incident energies
c flagexc     : flag for output of excitation functions
c flagnatural : flag for calculation of natural element
c flagadd     : flag for addition of discrete states to spectra 
c flagaddel   : flag for addition of elastic peak to spectra
c flagelectron: flag for application of electron conversion coefficient
c flagspher   : flag to force spherical optical model
c flagpartable: flag for output of model parameters on separate file
c maxchannel, : maximal number of outgoing particles in individual
c numchannel    channel description (e.g. this is 3 for (n,2np))
c Ztarget     : charge number of target nucleus
c fismodel    : fission model
c fismodelalt : alternative fission model for default barriers 
c flagendf    : flag for information for ENDF-6 file
c flagchannels: flag for exclusive channels calculation 
c k0          : index of incident particle 
c
      flagmain=.true.
      flagpop=flagbasic
      flagcheck=flagbasic
      flagoutomp=flagbasic
      flagdirect=flagbasic
      flaginverse=flagbasic
      flaggamma=flagbasic
      flaglevels=flagbasic
      flagdensity=flagbasic
      flagdisc=flagbasic
      flagfisout=flagbasic
      if (.not.flagfission) flagfisout=.false.
      flagtransen=.true.
      flagpeout=.false.
      flagang=.false.
      ddxmode=0
      flaglegendre=.false.
      flagspec=.false.
      if (flagrecoil) flagspec=.true.
      flagddx=.false.
      flagoutdwba=.false.
      flaggamdis=.false.
      flagoutecis=flageciscomp
c
c By default, we assume that with more than one incident energy output 
c of excitation functions (e.g. residual production cross sections as a 
c function of incident energy) are wanted.
c
      if (numinc.eq.1) then
        flagexc=.false.
      else
        flagexc=.true.
      endif
      if (flagnatural) flagexc=.true.
      flagadd=.true.
      flagaddel=flagadd
      flagelectron=.false.
      flagspher=.false.      
      flagpartable=.false.
      maxchannel=4
      fismodel=1
      fismodelalt=4     
c
c If the results of TALYS are used to create ENDF-6 data files,
c several output flags are automatically set.
c
      if (flagendf) then
        flagcheck=.true.
        flagdisc=.true.
        flagadd=.false.
        flagaddel=.false.
        if (flagfission) flagfisout=.true.
        flagang=.true.
        flaglegendre=.true.
        flagspec=.true.
        flagexc=.true.
        flagelectron=.true.
        if (k0.eq.1) then
          flaggamdis=.true.
          flagchannels=.true.
        endif
      endif
c
c **************** Read fourth set of input variables ******************
c
c nlines : number of input lines
c ch     : character
c inline : input line                 
c flagrot: flag for use of rotational optical model per
c          outgoing particle, if available       
c
c First the keyword is identified. Next the corresponding value is read.
c Erroneous input is immediately checked.
c  
      do 10 i=1,nlines
        if (inline(i)(1:8).eq.'outmain ') then
          do 20 i2=9,80
            ch=inline(i)(i2:i2)
            if (ch.ne.' ') then
              if (ch.eq.'n') flagmain=.false.
              if (ch.eq.'y') flagmain=.true.
              if (ch.ne.'y'.and.ch.ne.'n') goto 400
              goto 10
            endif
   20     continue
        endif
        if (inline(i)(1:14).eq.'outpopulation ') then
          do 30 i2=15,80
            ch=inline(i)(i2:i2)
            if (ch.ne.' ') then
              if (ch.eq.'n') flagpop=.false.
              if (ch.eq.'y') flagpop=.true.
              if (ch.ne.'y'.and.ch.ne.'n') goto 400
              goto 10
            endif
   30     continue
        endif
        if (inline(i)(1:9).eq.'outcheck ') then
          do 40 i2=10,80
            ch=inline(i)(i2:i2)
            if (ch.ne.' ') then
              if (ch.eq.'n') flagcheck=.false.
              if (ch.eq.'y') flagcheck=.true.
              if (ch.ne.'y'.and.ch.ne.'n') goto 400
              goto 10
            endif
   40     continue
        endif
        if (inline(i)(1:7).eq.'outomp ') then
          do 50 i2=8,80
            ch=inline(i)(i2:i2)
            if (ch.ne.' ') then
              if (ch.eq.'n') flagoutomp=.false.
              if (ch.eq.'y') flagoutomp=.true.
              if (ch.ne.'y'.and.ch.ne.'n') goto 400
              goto 10
            endif
   50     continue
        endif
        if (inline(i)(1:11).eq.'outinverse ') then
          do 60 i2=12,80
            ch=inline(i)(i2:i2)
            if (ch.ne.' ') then
              if (ch.eq.'n') flaginverse=.false.
              if (ch.eq.'y') flaginverse=.true.
              if (ch.ne.'y'.and.ch.ne.'n') goto 400
              goto 10
            endif
   60     continue
        endif
        if (inline(i)(1:15).eq.'outtransenergy ') then
          do 70 i2=16,80
            ch=inline(i)(i2:i2)
            if (ch.ne.' ') then
              if (ch.eq.'n') flagtransen=.false.
              if (ch.eq.'y') flagtransen=.true.
              if (ch.ne.'y'.and.ch.ne.'n') goto 400
              goto 10
            endif
   70     continue
        endif
        if (inline(i)(1:9).eq.'outgamma ') then
          do 80 i2=10,80
            ch=inline(i)(i2:i2)
            if (ch.ne.' ') then
              if (ch.eq.'n') flaggamma=.false.
              if (ch.eq.'y') flaggamma=.true.
              if (ch.ne.'y'.and.ch.ne.'n') goto 400
              goto 10
            endif
   80     continue
        endif
        if (inline(i)(1:10).eq.'outlevels ') then
          do 90 i2=11,80
            ch=inline(i)(i2:i2)
            if (ch.ne.' ') then
              if (ch.eq.'n') flaglevels=.false.
              if (ch.eq.'y') flaglevels=.true.
              if (ch.ne.'y'.and.ch.ne.'n') goto 400
              goto 10
            endif
   90     continue
        endif
        if (inline(i)(1:11).eq.'outdensity ') then
          do 100 i2=12,80
            ch=inline(i)(i2:i2)
            if (ch.ne.' ') then
              if (ch.eq.'n') flagdensity=.false.
              if (ch.eq.'y') flagdensity=.true.
              if (ch.ne.'y'.and.ch.ne.'n') goto 400
              goto 10
            endif
  100     continue
        endif
        if (inline(i)(1:11).eq.'outfission ') then
          do 110 i2=12,80
            ch=inline(i)(i2:i2)
            if (ch.ne.' ') then
              if (ch.eq.'n') flagfisout=.false.
              if (ch.eq.'y') flagfisout=.true.
              if (ch.ne.'y'.and.ch.ne.'n') goto 400
              goto 10
            endif
  110     continue
        endif
        if (inline(i)(1:18).eq.'outpreequilibrium ') then
          do 120 i2=19,80
            ch=inline(i)(i2:i2)
            if (ch.ne.' ') then
              if (ch.eq.'n') flagpeout=.false.
              if (ch.eq.'y') flagpeout=.true.
              if (ch.ne.'y'.and.ch.ne.'n') goto 400
              goto 10
            endif
  120     continue
        endif
        if (inline(i)(1:12).eq.'outdiscrete ') then
          do 130 i2=13,80
            ch=inline(i)(i2:i2)
            if (ch.ne.' ') then
              if (ch.eq.'n') flagdisc=.false.
              if (ch.eq.'y') flagdisc=.true.
              if (ch.ne.'y'.and.ch.ne.'n') goto 400
              goto 10
            endif
  130     continue
        endif
        if (inline(i)(1:11).eq.'outspectra ') then
          do 140 i2=12,80
            ch=inline(i)(i2:i2)
            if (ch.ne.' ') then
              if (ch.eq.'n') flagspec=.false.
              if (ch.eq.'y') flagspec=.true.
              if (ch.ne.'y'.and.ch.ne.'n') goto 400
              goto 10
            endif
  140     continue
        endif
        if (inline(i)(1:9).eq.'outangle ') then
          do 150 i2=10,80
            ch=inline(i)(i2:i2)
            if (ch.ne.' ') then
              if (ch.eq.'n') flagang=.false.
              if (ch.eq.'y') flagang=.true.
              if (ch.ne.'y'.and.ch.ne.'n') goto 400
              goto 10
            endif
  150     continue
        endif
        if (inline(i)(1:12).eq.'outlegendre ') then
          do 160 i2=13,80
            ch=inline(i)(i2:i2)
            if (ch.ne.' ') then
              if (ch.eq.'n') flaglegendre=.false.
              if (ch.eq.'y') flaglegendre=.true.
              if (ch.ne.'y'.and.ch.ne.'n') goto 400
              goto 10
            endif
  160     continue
        endif
        if (inline(i)(1:8).eq.'ddxmode ') then
          read(inline(i)(9:80),*,err=400) ddxmode
          if (ddxmode.eq.0) flagddx=.false.
          if (ddxmode.gt.0) then
            flagddx=.true.
            flagspec=.true.
          endif
          goto 10
        endif
        if (inline(i)(1:8).eq.'outdwba ') then
          do 170 i2=9,80
            ch=inline(i)(i2:i2)
            if (ch.ne.' ') then
              if (ch.eq.'n') flagoutdwba=.false.
              if (ch.eq.'y') flagoutdwba=.true.
              if (ch.ne.'y'.and.ch.ne.'n') goto 400
              goto 10
            endif
  170     continue
        endif
        if (inline(i)(1:10).eq.'outgamdis ') then
          do 180 i2=11,80
            ch=inline(i)(i2:i2)
            if (ch.ne.' ') then
              if (ch.eq.'n') flaggamdis=.false.
              if (ch.eq.'y') flaggamdis=.true.
              if (ch.ne.'y'.and.ch.ne.'n') goto 400
              goto 10
            endif
  180     continue
        endif
        if (inline(i)(1:14).eq.'outexcitation ') then
          do 190 i2=15,80
            ch=inline(i)(i2:i2)
            if (ch.ne.' ') then
              if (ch.eq.'n') flagexc=.false.
              if (ch.eq.'y') flagexc=.true.
              if (ch.ne.'y'.and.ch.ne.'n') goto 400
              goto 10
            endif
  190     continue
        endif
        if (inline(i)(1:8).eq.'outecis ') then
          do 200 i2=9,80
            ch=inline(i)(i2:i2)
            if (ch.ne.' ') then
              if (ch.eq.'n') flagoutecis=.false.
              if (ch.eq.'y') flagoutecis=.true.
              if (ch.ne.'y'.and.ch.ne.'n') goto 400
              goto 10
            endif
  200     continue
        endif
        if (inline(i)(1:10).eq.'outdirect ') then
          do 210 i2=11,80
            ch=inline(i)(i2:i2)
            if (ch.ne.' ') then
              if (ch.eq.'n') flagdirect=.false.
              if (ch.eq.'y') flagdirect=.true.
              if (ch.ne.'y'.and.ch.ne.'n') goto 400
              goto 10
            endif
  210     continue
        endif
        if (inline(i)(1:12).eq.'adddiscrete ') then
          do 220 i2=13,80
            ch=inline(i)(i2:i2)
            if (ch.ne.' ') then
              if (ch.eq.'n') flagadd=.false.
              if (ch.eq.'y') flagadd=.true.
              if (ch.ne.'y'.and.ch.ne.'n') goto 400
              goto 10
            endif
  220     continue
        endif
        if (inline(i)(1:11).eq.'addelastic ') then
          do 230 i2=12,80
            ch=inline(i)(i2:i2)
            if (ch.ne.' ') then
              if (ch.eq.'n') flagaddel=.false.
              if (ch.eq.'y') flagaddel=.true.
              if (ch.ne.'y'.and.ch.ne.'n') goto 400
              goto 10
            endif
  230     continue
        endif
        if (inline(i)(1:13).eq.'electronconv ') then
          do 240 i2=14,80
            ch=inline(i)(i2:i2)
            if (ch.ne.' ') then
              if (ch.eq.'n') flagelectron=.false.
              if (ch.eq.'y') flagelectron=.true.
              if (ch.ne.'y'.and.ch.ne.'n') goto 400
              goto 10
            endif
  240     continue
        endif                    
        if (inline(i)(1:10).eq.'spherical ') then
          do 250 i2=11,80
            ch=inline(i)(i2:i2)
            if (ch.ne.' ') then
              if (ch.eq.'n') flagspher=.false.
              if (ch.eq.'y') then
                flagspher=.true.
                do 260 type=1,6
                  flagrot(type)=.false.
  260           continue
              endif
              if (ch.ne.'y'.and.ch.ne.'n') goto 400
              goto 10
            endif
  250     continue
        endif               
        if (inline(i)(1:9).eq.'partable ') then
          do 270 i2=10,80
            ch=inline(i)(i2:i2)
            if (ch.ne.' ') then
              if (ch.eq.'n') flagpartable=.false.
              if (ch.eq.'y') flagpartable=.true.
              if (ch.ne.'y'.and.ch.ne.'n') goto 400
              goto 10
            endif
  270     continue
        endif                    
        if (inline(i)(1:11).eq.'maxchannel ') 
     +    read(inline(i)(12:80),*,err=400) maxchannel
        if (inline(i)(1:9).eq.'fismodel ')
     +    read(inline(i)(10:80),*,err=400) fismodel
        if (inline(i)(1:12).eq.'fismodelalt ')
     +    read(inline(i)(13:80),*,err=400) fismodelalt   
        goto 10
  400   write(*,'("TALYS-error: Wrong input: ",a80)') inline(i)
        stop
   10 continue
      return
      end
Copyright (C) 2004  A.J. Koning, S. Hilaire and M.C. Duijvestijn
