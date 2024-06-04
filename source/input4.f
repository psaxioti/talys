      subroutine input4
c
c +---------------------------------------------------------------------
c | Author: Arjan Koning 
c | Date  : October 19, 2006   
c | Task  : Read input for fourth set of variables
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      character*1  ch
      character*80 word(40),key,value
      integer      i,type
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
c flagcol     : flag for collective enhancement of level density
c flagpartable: flag for output of model parameters on separate file
c maxchannel, : maximal number of outgoing particles in individual
c numchannel    channel description (e.g. this is 3 for (n,2np))
c Ztarget     : charge number of target nucleus
c pairmodel   : model for preequilibrium pairing energy
c fismodel    : fission model
c fismodelalt : alternative fission model for default barriers 
c flagendf    : flag for information for ENDF-6 file
c flagchannels: flag for exclusive channels calculation 
c flagendfdet : flag for detailed ENDF-6 information per channel
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
      pairmodel=2
      fismodel=1
      fismodelalt=4     
      if (flagfission) then
        flagcol=.true.
      else
        flagcol=.false.
      endif
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
        if (flagendfdet) then
          flaggamdis=.true.
          flagchannels=.true.
        endif
      endif
c
c **************** Read fourth set of input variables ******************
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
      do 10 i=1,nlines
        call getkeywords(inline(i),word)
        key=word(1)
        value=word(2)
        ch=word(2)(1:1)
c
c Test for keywords
c
c flagrot: flag for use of rotational optical model per
c          outgoing particle, if available       
c
        if (key.eq.'outmain') then
          if (ch.eq.'n') flagmain=.false.
          if (ch.eq.'y') flagmain=.true.
          if (ch.ne.'y'.and.ch.ne.'n') goto 200
          goto 10
        endif
        if (key.eq.'outpopulation') then
          if (ch.eq.'n') flagpop=.false.
          if (ch.eq.'y') flagpop=.true.
          if (ch.ne.'y'.and.ch.ne.'n') goto 200
          goto 10
        endif
        if (key.eq.'outcheck') then
          if (ch.eq.'n') flagcheck=.false.
          if (ch.eq.'y') flagcheck=.true.
          if (ch.ne.'y'.and.ch.ne.'n') goto 200
          goto 10
        endif
        if (key.eq.'outomp') then
          if (ch.eq.'n') flagoutomp=.false.
          if (ch.eq.'y') flagoutomp=.true.
          if (ch.ne.'y'.and.ch.ne.'n') goto 200
          goto 10
        endif
        if (key.eq.'outinverse') then
          if (ch.eq.'n') flaginverse=.false.
          if (ch.eq.'y') flaginverse=.true.
          if (ch.ne.'y'.and.ch.ne.'n') goto 200
          goto 10
        endif
        if (key.eq.'outtransenergy') then
          if (ch.eq.'n') flagtransen=.false.
          if (ch.eq.'y') flagtransen=.true.
          if (ch.ne.'y'.and.ch.ne.'n') goto 200
          goto 10
        endif
        if (key.eq.'outgamma') then
          if (ch.eq.'n') flaggamma=.false.
          if (ch.eq.'y') flaggamma=.true.
          if (ch.ne.'y'.and.ch.ne.'n') goto 200
          goto 10
        endif
        if (key.eq.'outlevels') then
          if (ch.eq.'n') flaglevels=.false.
          if (ch.eq.'y') flaglevels=.true.
          if (ch.ne.'y'.and.ch.ne.'n') goto 200
          goto 10
        endif
        if (key.eq.'outdensity') then
          if (ch.eq.'n') flagdensity=.false.
          if (ch.eq.'y') flagdensity=.true.
          if (ch.ne.'y'.and.ch.ne.'n') goto 200
          goto 10
        endif
        if (key.eq.'outfission') then
          if (ch.eq.'n') flagfisout=.false.
          if (ch.eq.'y') flagfisout=.true.
          if (ch.ne.'y'.and.ch.ne.'n') goto 200
          goto 10
        endif
        if (key.eq.'outpreequilibrium') then
          if (ch.eq.'n') flagpeout=.false.
          if (ch.eq.'y') flagpeout=.true.
          if (ch.ne.'y'.and.ch.ne.'n') goto 200
          goto 10
        endif
        if (key.eq.'outdiscrete') then
          if (ch.eq.'n') flagdisc=.false.
          if (ch.eq.'y') flagdisc=.true.
          if (ch.ne.'y'.and.ch.ne.'n') goto 200
          goto 10
        endif
        if (key.eq.'outspectra') then
          if (ch.eq.'n') flagspec=.false.
          if (ch.eq.'y') flagspec=.true.
          if (ch.ne.'y'.and.ch.ne.'n') goto 200
          goto 10
        endif
        if (key.eq.'outangle') then
          if (ch.eq.'n') flagang=.false.
          if (ch.eq.'y') flagang=.true.
          if (ch.ne.'y'.and.ch.ne.'n') goto 200
          goto 10
        endif
        if (key.eq.'outlegendre') then
          if (ch.eq.'n') flaglegendre=.false.
          if (ch.eq.'y') flaglegendre=.true.
          if (ch.ne.'y'.and.ch.ne.'n') goto 200
          goto 10
        endif
        if (key.eq.'ddxmode') then
          read(value,*,err=200) ddxmode
          if (ddxmode.eq.0) flagddx=.false.
          if (ddxmode.gt.0) then
            flagddx=.true.
            flagspec=.true.
          endif
          goto 10
        endif
        if (key.eq.'outdwba') then
          if (ch.eq.'n') flagoutdwba=.false.
          if (ch.eq.'y') flagoutdwba=.true.
          if (ch.ne.'y'.and.ch.ne.'n') goto 200
          goto 10
        endif
        if (key.eq.'outgamdis') then
          if (ch.eq.'n') flaggamdis=.false.
          if (ch.eq.'y') flaggamdis=.true.
          if (ch.ne.'y'.and.ch.ne.'n') goto 200
          goto 10
        endif
        if (key.eq.'outexcitation') then
          if (ch.eq.'n') flagexc=.false.
          if (ch.eq.'y') flagexc=.true.
          if (ch.ne.'y'.and.ch.ne.'n') goto 200
          goto 10
        endif
        if (key.eq.'outecis') then
          if (ch.eq.'n') flagoutecis=.false.
          if (ch.eq.'y') flagoutecis=.true.
          if (ch.ne.'y'.and.ch.ne.'n') goto 200
          goto 10
        endif
        if (key.eq.'outdirect') then
          if (ch.eq.'n') flagdirect=.false.
          if (ch.eq.'y') flagdirect=.true.
          if (ch.ne.'y'.and.ch.ne.'n') goto 200
          goto 10
        endif
        if (key.eq.'adddiscrete') then
          if (ch.eq.'n') flagadd=.false.
          if (ch.eq.'y') flagadd=.true.
          if (ch.ne.'y'.and.ch.ne.'n') goto 200
          goto 10
        endif
        if (key.eq.'addelastic') then
          if (ch.eq.'n') flagaddel=.false.
          if (ch.eq.'y') flagaddel=.true.
          if (ch.ne.'y'.and.ch.ne.'n') goto 200
          goto 10
        endif
        if (key.eq.'electronconv') then
          if (ch.eq.'n') flagelectron=.false.
          if (ch.eq.'y') flagelectron=.true.
          if (ch.ne.'y'.and.ch.ne.'n') goto 200
          goto 10
        endif                    
        if (key.eq.'spherical') then
          if (ch.eq.'n') flagspher=.false.
          if (ch.eq.'y') then
            flagspher=.true.
            do 110 type=1,6
              flagrot(type)=.false.
  110       continue
          endif
          if (ch.ne.'y'.and.ch.ne.'n') goto 200
          goto 10
        endif               
        if (key.eq.'colenhance') then
          if (ch.eq.'n') flagcol=.false.
          if (ch.eq.'y') flagcol=.true.
          if (ch.ne.'y'.and.ch.ne.'n') goto 200
          goto 10
        endif
        if (key.eq.'partable') then
          if (ch.eq.'n') flagpartable=.false.
          if (ch.eq.'y') flagpartable=.true.
          if (ch.ne.'y'.and.ch.ne.'n') goto 200
          goto 10
        endif                    
        if (key.eq.'maxchannel') then
          read(value,*,err=200) maxchannel
          goto 10
        endif                    
        if (key.eq.'pairmodel') then
          read(value,*,err=200) pairmodel
          goto 10
        endif                    
        if (key.eq.'fismodel') then
          read(value,*,err=200) fismodel
          goto 10
        endif                    
        if (key.eq.'fismodelalt') then
          read(value,*,err=200) fismodelalt   
          goto 10
        endif                    
   10 continue
      return
  200 write(*,'(" TALYS-error: Wrong input: ",a80)') inline(i)
      stop
      end
Copyright (C) 2004  A.J. Koning, S. Hilaire and M.C. Duijvestijn
