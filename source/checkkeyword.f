      subroutine checkkeyword
c
c +---------------------------------------------------------------------
c | Author: Arjan Koning
c | Date  : November 22, 2006
c | Task  : Check for errors in keywords
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      integer      numkey,i,j
      parameter    (numkey=211)
      character*80 keyword(numkey),word(40),key
c
c Although it is difficult to prevent the user from all possible input
c errors, we can check for the use of wrong keywords and for unphysical
c values for most of the input variables.
c
c *********************** Check for wrong keywords *********************
c
c keyword: keyword
c numkey : number of keywords
c
c TALYS will stop if a keyword is incorrect
c
      data (keyword(i),i=1,numkey) /
     +  ' ', 'a',
     +  'abundance', 'adddiscrete',
     +  'addelastic', 'alimit', 'alphald', 'angles',
     +  'anglescont', 'anglesrec', 'asys',
     +  'autorot', 'avadjust', 'avdadjust', 'avsoadjust', 'axtype', 
     +  'beta2', 'betald', 'bins', 'c1table', 
     +  'c2table', 'cfermi', 'cfermibf', 'channelenergy',
     +  'channels', 'cknock', 'class2', 'class2file',
     +  'class2width', 'colenhance', 'compound',
     +  'core', 'cstrip', 'd0', 'd1adjust', 'd2adjust', 
     +  'd3adjust', 'ddxmode', 'deformfile', 'deltaw',
     +  'e0', 'eciscalc',
     +  'eciscompound', 'ecisdwba', 'ecissave',
     +  'egr', 'ejectiles',
     +  'electronconv', 'element',
     +  'elow', 'elwidth',
     +  'emsdmin', 'endf', 'endfdetail',
     +  'energy', 'esurf', 'exmatch', 'expmass',
     +  'ffevaporation', 'fileangle',
     +  'filechannels', 'fileddxa',
     +  'fileddxe', 'filedensity', 'filediscrete',
     +  'fileelastic', 'filefission',
     +  'filegamdis', 'filerecoil', 'fileresidual',
     +  'filespectrum', 'filetotal',
     +  'fisbar', 'fishw',
     +  'fismodel', 'fismodelalt',
     +  'fission', 'fullhf', 'g',
     +  'gamgam',
     +  'gammald', 'gammashell1', 'gammashell2', 'gammax',
     +  'ggr',
     +  'giantresonance', 'gn', 'gnorm',
     +  'gp', 'gshell', 'hbtransfile',
     +  'inccalc', 'isomer',
     +  'kph', 'krotconstant', 'labddx', 'ldmodel', 'levelfile',
     +  'localomp', 'ltarget',
     +  'm2constant', 'm2limit', 'm2shift',
     +  'mass', 'massdis', 'massexcess', 'massmodel', 'massnucleus',
     +  'maxband', 'maxchannel', 'maxenrec',
     +  'maxlevelsbin', 'maxlevelsres',
     +  'maxlevelstar', 'maxn', 'maxrot',
     +  'maxz', 'mpreeqmode',
     +  'msdbins',
     +  'multipreeq', 'nlevels',
     +  'nlow', 'ntop', 'nulldev',
     +  'onestep', 'optmod',
     +  'optmodall', 'optmodfilen', 'optmodfilep',
     +  'outangle', 'outbasic', 'outcheck',
     +  'outdensity', 'outdirect',
     +  'outdiscrete', 'outdwba',
     +  'outecis', 'outexcitation',
     +  'outfission', 'outgamdis',
     +  'outgamma', 'outinverse',
     +  'outlegendre', 'outlevels',
     +  'outmain', 'outomp',
     +  'outpopulation', 'outpreequilibrium',
     +  'outspectra', 'outtransenergy',
     +  'pair', 'pairconstant', 'pairmodel', 'parity', 'partable', 
     +  'popeps', 'preeqcomplex', 'preeqmode',
     +  'preeqspin', 'preeqsurface',
     +  'preequilibrium', 'projectile', 'pshift', 'pshiftconstant',
     +  'rcadjust', 'rclass2mom', 'reaction', 'recoil',
     +  'recoilaverage', 'relativistic', 'rfiseps', 'rgamma',
     +  'rnunu', 'rnupi', 'rotational', 'rpinu', 'rpipi',
     +  'rspincut', 'rtransmom', 'rvadjust', 'rvdadjust', 
     +  'rvsoadjust', 's0', 'segment', 'sgr', 'shellmodel',
     +  'spherical', 'spincutmodel', 'statepot',
     +  'strength', 'strucpath', 'sysreaction',
     +  't', 'transeps',
     +  'transpower', 'twocomponent',
     +  'ufermi', 'ufermibf',
     +  'v1adjust', 'v2adjust', 'v3adjust', 'v4adjust',
     +  'vso1adjust', 'vso2adjust',
     +  'w1adjust', 'w2adjust', 'wso1adjust', 'wso2adjust',
     +  'widthfluc', 'widthmode', 'xseps'/
c
c A keyword can be de-activated by putting a # in front of it.
c All first words of the input lines are checked against the list
c of keywords.
c
c nlines     : number of input lines
c getkeywords: subroutine to retrieve keywords and values from input 
c              line
c inline     : input line                 
c word       : words on input line
c key        : keyword
c
c The keyword is identified.
c
      do 10 i=1,nlines
        call getkeywords(inline(i),word)
        key=word(1)
        if (key(1:1).eq.'#') goto 10
        do 20 j=1,numkey
          if (keyword(j).eq.key) goto 10
   20   continue
        write(*,'(/" TALYS-error: Wrong keyword: ",a20)') key
        stop
   10 continue
      return
      end
Copyright (C) 2004  A.J. Koning, S. Hilaire and M.C. Duijvestijn
