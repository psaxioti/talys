      subroutine checkkeyword
c
c +---------------------------------------------------------------------
c | Author: Arjan Koning
c | Date  : December 11, 2009
c | Task  : Check for errors in keywords
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      integer      numkey,i,j
      parameter    (numkey=265)
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
     +  ' ', 'a', 'aadjust',
     +  'abundance', 'adddiscrete',
     +  'addelastic', 'alimit', 'alphald', 'alphaomp', 'angles',
     +  'anglescont', 'anglesrec', 'astro', 'astrogs', 'asys',
     +  'autorot', 'avadjust', 'avadjustf', 'avdadjust', 'avdadjustf', 
     +  'avsoadjust', 'avsoadjustf', 'axtype', 'best', 'bestpath', 
     +  'beta2', 'betafiscor', 'betald', 'bins', 'cfermi',
     +  'cfermibf', 'cglobal', 'channelenergy',
     +  'channels', 'cknock', 'class2', 'class2file',
     +  'class2width', 'colenhance', 'colldamp', 'compound', 'core',
     +  'coulomb', 'cstrip', 'ctable', 'ctmglobal', 'd0', 'd1adjust',
     +  'd2adjust', 'd3adjust', 'ddxmode', 'deformfile', 'deltaw', 
     +  'dispersion', 'e0', 'eciscalc',
     +  'eciscompound', 'ecisdwba', 'ecissave',
     +  'egr', 'ejectiles', 'electronconv', 'element',
     +  'elow', 'elwidth', 'emsdmin', 'endf', 'endfdetail', 'endfecis',
     +  'energy', 'epr', 'esurf', 'etable', 'exmatch', 'expmass',
     +  'ffevaporation', 'fileangle', 'filechannels', 'fileddxa',
     +  'fileddxe', 'filedensity', 'filediscrete',
     +  'fileelastic', 'filefission', 'filegamdis', 'filerecoil', 
     +  'fileresidual', 'filespectrum', 'filetotal', 'fisbar', 'fishw',
     +  'fismodel', 'fismodelalt', 'fiso', 'fission', 'ftable', 
     +  'fullhf', 'g', 'gamgam', 'gamgamadjust',
     +  'gammald', 'gammashell1', 'gammashell2', 'gammax',
     +  'ggr', 'giantresonance', 'gn', 'gnadjust', 'gnorm',
     +  'gp', 'gpadjust', 'gpr', 'gshell', 'hbstate', 'hbtransfile',
     +  'inccalc', 'isomer', 'jlmmode', 'jlmomp', 'kph', 'krotconstant',
     +  'kvibmodel', 'labddx', 'ldmodel', 'levelfile',
     +  'localomp', 'ltarget', 'lvadjust', 'lv1adjust', 'lwadjust',
     +  'lw1adjust', 'lvsoadjust', 'lwsoadjust', 'm2constant', 
     +  'm2limit', 'm2shift', 'mass', 'massdis', 'massexcess', 
     +  'massmodel', 'massnucleus', 'maxband', 'maxchannel', 
     +  'maxenrec', 'maxlevelsbin', 'maxlevelsres',
     +  'maxlevelstar', 'maxn', 'maxrot', 'maxz', 'micro', 
     +  'mpreeqmode', 'msdbins', 'multipreeq', 'nlevels',
     +  'nlow', 'ntop', 'nulldev', 'ompenergyfile', 'omponly', 
     +  'onestep', 'optmod', 'optmodall', 'optmodfilen', 'optmodfilep',
     +  'outangle', 'outbasic', 'outcheck', 'outdensity', 'outdirect',
     +  'outdiscrete', 'outdwba', 'outecis', 'outexcitation',
     +  'outfission', 'outgamdis', 'outgamma', 'outinverse',
     +  'outlegendre', 'outlevels', 'outmain', 'outomp',
     +  'outpopulation', 'outpreequilibrium',
     +  'outspectra', 'outtransenergy',
     +  'pair', 'pairconstant', 'pairmodel', 'parity', 'partable', 
     +  'pglobal', 'phmodel', 'popeps', 'preeqcomplex', 'preeqmode',
     +  'preeqspin', 'preeqsurface', 'preequilibrium', 
     +  'projectile', 'pshift', 'pshiftconstant', 'ptable',
     +  'radialfile', 'radialmodel', 'rcadjust', 'rclass2mom',
     +  'reaction', 'recoil', 'recoilaverage', 'relativistic', 
     +  'rescuefile', 'rfiseps', 'rgamma', 'rnunu', 'rnupi', 
     +  'rotational', 'rpinu', 'rpipi', 'rspincut', 'rtransmom', 
     +  'rvadjust', 'rvadjustf', 'rvdadjust', 'rvdadjustf', 
     +  'rvsoadjust', 'rvsoadjustf', 's0', 's2adjust',
     +  'segment', 'sgr', 'shellmodel', 'soswitch', 'spherical', 
     +  'spincutmodel', 'spr', 'statepot', 'strength', 'strengthm1', 
     +  'strucpath', 'sysreaction', 't', 'transeps', 'transpower', 
     +  'twocomponent', 'ufermi', 'ufermibf', 'urr', 
     +  'v1adjust', 'v2adjust', 'v3adjust', 'v4adjust',
     +  'vfiscor', 'vso1adjust', 'vso2adjust',
     +  'w1adjust', 'w2adjust', 'wso1adjust', 'wso2adjust',
     +  'widthfluc', 'widthmode', 'xsalphatherm', 'xscaptherm',
     +  'xseps', 'xsptherm'/
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
