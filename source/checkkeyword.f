      subroutine checkkeyword
c
c +---------------------------------------------------------------------
c | Author: Arjan Koning
c | Date  : December 2, 2004
c | Task  : Check for errors in keywords
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      integer      numkey,i,ii,jj,kk
      parameter    (numkey=170)
      character*80 keyword(numkey),inkey
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
     +  ' ', 'a ',
     +  'abundance ', 'adddiscrete ',
     +  'addelastic ', 'alimit ', 'alphald ',
     +  'angles ', 'anglescont ', 'anglesrec ',
     +  'asys ', 'autorot ', 'axtype ', 'beta2 ', 'betald', 'bins ',
     +  'c1table ', 'c2table ',
     +  'cfermi ', 'cfermibf', 'channelenergy ',
     +  'channels ','class2 ', 'class2file ', 'class2width ', 
     +  'colldamp ', 'compound ', 'core ',
     +  'd0 ', 'ddxmode ', 'deformfile ',
     +  'deltaw ',
     +  'e0 ', 'eciscalc ',
     +  'eciscompound ', 'ecisdwba ', 'ecissave ',
     +  'egr ', 'ejectiles ',
     +  'electronconv ', 'element ',
     +  'elow ', 'elwidth ',
     +  'emsdmin ', 'endf ',
     +  'energy ', 'esurf ', 'exmatch ',
     +  'ffevaporation ', 'fileangle ',
     +  'filechannels ', 'fileddxa ',
     +  'fileddxe ', 'filediscrete ',
     +  'fileelastic ', 'filefission ',
     +  'filegamdis ', 'filerecoil ', 'fileresidual ',
     +  'filespectrum ', 'filetotal ',
     +  'fisbar ', 'fishw ',
     +  'fismodel ', 'fismodelalt ',
     +  'fission ', 'fullhf ', 'g ',
     +  'gamgam ',
     +  'gammald ', 'gammashell1 ', 'gammashell2 ', 'gammax ',
     +  'gnorm ', 'ggr ',
     +  'giantresonance ', 'gn ',
     +  'gp ', 'gshell ', 'hbtransfile ',
     +  'inccalc ', 'isomer ',
     +  'kph ', 'krotconstant ', 'labddx ', 'ldmodel ', 'levelfile ',
     +  'localomp ', 'ltarget ',
     +  'm2constant ', 'm2limit ', 'm2shift ',
     +  'mass ', 'massdis ', 'maxband ', 'maxchannel ',
     +  'maxenrec ',
     +  'maxlevelsbin ', 'maxlevelsres ',
     +  'maxlevelstar ', 'maxn ', 'maxrot ',
     +  'maxz ', 'mpreeqmode ',
     +  'msdbins ',
     +  'multipreeq ', 'nlevels ',
     +  'nlow ', 'ntop ',
     +  'onestep ', 'optmod ',
     +  'optmodall ', 'optmodfilen ', 'optmodfilep ',
     +  'outangle ', 'outbasic ', 'outcheck ',
     +  'outdensity ', 'outdirect ',
     +  'outdiscrete ', 'outdwba ',
     +  'outecis ', 'outexcitation ',
     +  'outfission ', 'outgamdis ',
     +  'outgamma ', 'outinverse ',
     +  'outlegendre ', 'outlevels ',
     +  'outmain ', 'outomp ',
     +  'outpopulation ', 'outpreequilibrium ',
     +  'outspectra ', 'outtransenergy ',
     +  'pair ', 'partable ', 'popeps ',
     +  'preeqcomplex ', 'preeqmode ',
     +  'preeqspin ', 'preeqsurface ',
     +  'preequilibrium ', 'projectile ', 'rclass2mom ', 
     +  'recoil ', 'recoilaverage ', 
     +  'relativistic ', 'rfiseps ', 'rgamma ',
     +  'rotational ', 'rpinu ',
     +  'rspincut ', 'rtransmom ', 's0 ',
     +  'segment ', 'sgr ',
     +  'spherical ', 'statepot ',
     +  'strength ', 'sysreaction ',
     +  't ', 'transeps ',
     +  'transpower ', 'twocomponent ',
     +  'ufermi ', 'ufermibf', 'widthfluc ',
     +  'widthmode ', 'xseps '/
c
c A keyword can be de-activated by putting a # in the first column.
c All first words of the input lines are checked against the list
c of keywords.
c
c nlines: number of input lines
c inkey : keyword from input
c
      do 10 i=1,nlines
        inkey='                                                        '
        if (inline(i)(1:1).eq.'#') goto 10
        do 20 ii=1,80
          if (inline(i)(ii:ii).eq.' ') then
            inkey(1:ii)=inline(i)(1:ii)
            goto 30
          endif
   20   continue
   30   do 40 jj=1,numkey
          if (keyword(jj)(1:ii).eq.inkey(1:ii)) then
            if (jj.eq.1) then
              do 50 kk=2,80
                if (inline(i)(kk:kk).ne.' ') then
                  write(*,'(/"Remove leading blanks: ",a20)') inline(i)
                  stop
                endif
   50         continue
              goto 10
            else
              goto 10
            endif
          endif
   40   continue
        write(*,'(/"Wrong keyword: ",a20)') inkey
        stop
   10 continue
      return
      end
Copyright (C) 2004  A.J. Koning, S. Hilaire and M.C. Duijvestijn
