      subroutine input3
c
c +---------------------------------------------------------------------
c | Author: Arjan Koning and Marieke Duijvestijn
c | Date  : December 2, 2004
c | Task  : Read input for third set of variables
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      character*1  ch
      character*72 ompf
      integer      Zix,Nix,type,i,i2,i3,i4,iz,ia
c
c ************** Defaults for third set of input variables *************
c
c Some defaults are set on the basis of the input variables specified
c before. This can of course be overruled in the input file.
c
c flagecissave: flag for saving ECIS input and output files
c flageciscalc: flag for new ECIS calculation for outgoing particles
c               and energy grid
c flaginccalc : flag for new ECIS calculation for incident channel
c flagrel     : flag for relativistic kinematics
c flagcomp    : flag for compound nucleus calculation
c flagfullhf  : flag for full spin dependence of transmission
c               coefficients         
c ewfc        : off-set incident energy for width fluctuation 
c               calculation
c epreeq      : on-set incident energy for preequilibrium calculation
c emulpre     : on-set incident energy for multiple preequilibrium 
c               calculation
c flagpespin  : flag for pre-equilibrium spin distribution or compound
c               spin distribution for pre-equilibrium cross section
c maxband     : highest vibrational level added to rotational model
c maxrot      : number of included excited rotational levels
c k0          : index for incident particle
c strength    : strength function of Kopecky-Uhl (1) or Brink-Axel (2)
c flagpecomp  : flag for Kalbach complex particle emission model 
c flagsurface : flag for surface effects in exciton model
c flaggiant0  : flag for collective contribution from giant resonances
c flag2comp   : flag for two-component pre-equilibrium model
c flagchannels: flag for exclusive channels calculation 
c flagfission : flag for fission
c Atarget     : mass number of target nucleus
c flagcolldamp: flag for damping of collective effects
c flagclass2  : flag for class2 states in fission
c flagbasic   : flag for output of basic information and results
c flageciscomp: flag for compound nucleus calculation by ECIS
c flagecisdwba: flag for new ECIS calculation for DWBA for MSD
c flagonestep : flag for continuum one-step direct only
c flaglocalomp: flag for local (y) or global (n) optical model 
c flagompall  : flag for new optical model calculation for all 
c               residual nuclei
c flagautorot : flag for automatic rotational coupled channels
c               calculations for A > 150
c flagstate   : flag for optical model potential for each excited state 
c Zix         : charge number index for residual nucleus
c numZ        : maximal number of protons away from the initial 
c               compound nucleus
c Nix         : neutron number index for residual nucleus
c numN        : maximal number of neutrons away from the initial 
c               compound nucleus
c optmod      : file with optical model parameters
c flagsys     : flag for reaction cross section from systematics
c flagrot     : flag for use of rotational optical model per 
c               outgoing particle, if available
c flagasys    : flag for all level density parameters a from 
c               systematics 
c flaggshell  : flag for energy dependence of single particle level 
c               density parameter g
c flagmassdis : flag for calculation of fission fragment mass yields
c flagffevap  : flag for calculation of particle evaporation from 
c               fission fragment mass yields
c flagendf    : flag for information for ENDF-6 file  
c flagrecoil  : flag for calculation of recoils 
c flaglabddx  : flag for calculation of DDX in LAB system 
c flagrecoilav: flag for average velocity in recoil calculation
c flagEchannel: flag for channel energy for emission spectrum
c
      flagecissave=.false.
      flageciscalc=.true.
      flaginccalc=.true.
      flagrel=.true.
      flagcomp=.true.
      flagfullhf=.false.
      ewfc=-1.
      epreeq=-1.
      emulpre=20.
      flagpespin=.false.
      maxband=0
      maxrot=2
      if (k0.ge.1) then
        strength=1
        flagpecomp=.true.
        flagsurface=.true.
      else
        strength=2
        flagpecomp=.false.
        flagsurface=.false.
      endif
      if (k0.eq.1.or.k0.eq.2) then
        flaggiant0=.true.
      else
        flaggiant0=.false.
      endif
      flag2comp=.true.
      flagchannels=.false.
      flagfission=.false.
      if (Atarget.gt.209) flagfission=.true.
      flagcolldamp=.false.
      flagclass2=.false.
      flagbasic=.false.
      flageciscomp=.false.
      flagecisdwba=.true.
      flagonestep=.false.
      flaglocalomp=.true.
      flagompall=.false.
      flagautorot=.false.
      flagstate=.false.
      do 10 Zix=0,numZ
        do 10 Nix=0,numN
          do 10 type=1,6
            optmod(Zix,Nix,type)=
     +        '                                                        '
   10 continue
      do 20 type=0,2
        flagsys(type)=.false.
   20 continue   
      do 30 type=3,6
        flagsys(type)=.true.
   30 continue   
      flagrot(0)=.false.
      do 40 type=1,2
        flagrot(type)=.true.
   40 continue   
      do 50 type=3,6
        flagrot(type)=.false.
   50 continue   
      flagasys=.true.
      flaggshell=.false.
      flagmassdis=.false.
      flagffevap=.false.
      flagendf=.false.
      flagrecoil=.false.
      flaglabddx=.false.
      flagrecoilav=.false.
      flagEchannel=.false.
c
c **************** Read third set of input variables *******************
c
c ch    : character
c nlines: number of input lines
c inline: input line                 
c Zinit : charge number of initial compound nucleus
c Ninit : neutron number of initial compound nucleus
c parsym: symbol of particle 
c
c First the keyword is identified. Next the corresponding value is read.
c Erroneous input is immediately checked.    
c
      do 110 i=1,nlines
        if (inline(i)(1:8).eq.'maxband ') 
     +    read(inline(i)(9:80),*,err=700) maxband
        if (inline(i)(1:7).eq.'maxrot ') 
     +    read(inline(i)(8:80),*,err=700) maxrot
        if (inline(i)(1:9).eq.'strength ') 
     +    read(inline(i)(10:80),*,err=700) strength
        if (inline(i)(1:9).eq.'ecissave ') then
          do 120 i2=10,80
            ch=inline(i)(i2:i2)
            if (ch.ne.' ') then
              if (ch.eq.'n') flagecissave=.false.
              if (ch.eq.'y') flagecissave=.true.
              if (ch.ne.'y'.and.ch.ne.'n') goto 700
              goto 110
            endif
  120     continue
        endif
        if (inline(i)(1:9).eq.'eciscalc ') then
          do 130 i2=10,80
            ch=inline(i)(i2:i2)
            if (ch.ne.' ') then
              if (ch.eq.'n') flageciscalc=.false.
              if (ch.eq.'y') flageciscalc=.true.
              if (ch.ne.'y'.and.ch.ne.'n') goto 700
              goto 110
            endif
  130     continue
        endif
        if (inline(i)(1:8).eq.'inccalc ') then
          do 140 i2=9,80
            ch=inline(i)(i2:i2)
            if (ch.ne.' ') then
              if (ch.eq.'n') flaginccalc=.false.
              if (ch.eq.'y') flaginccalc=.true.
              if (ch.ne.'y'.and.ch.ne.'n') goto 700
              goto 110
            endif
  140     continue
        endif
        if (inline(i)(1:13).eq.'relativistic ') then
          do 150 i2=14,80
            ch=inline(i)(i2:i2)
            if (ch.ne.' ') then
              if (ch.eq.'n') flagrel=.false.
              if (ch.eq.'y') flagrel=.true.
              if (ch.ne.'y'.and.ch.ne.'n') goto 700
              goto 110
            endif
  150     continue
        endif
        if (inline(i)(1:9).eq.'compound ') then
          do 160 i2=10,80
            ch=inline(i)(i2:i2)
            if (ch.ne.' ') then
              if (ch.eq.'n') flagcomp=.false.
              if (ch.eq.'y') flagcomp=.true.
              if (ch.ne.'y'.and.ch.ne.'n') goto 700
              goto 110
            endif
  160     continue
        endif
        if (inline(i)(1:10).eq.'widthfluc ') then
          do 170 i2=11,80
            ch=inline(i)(i2:i2)
            if (ch.ne.' ') then
              if (ch.eq.'y') then
                ewfc=20.
                goto 110
              endif
              if (ch.eq.'n') then
                ewfc=0.
                goto 110
              endif
              goto 180
            endif
  170     continue
  180     read(inline(i)(11:80),*,err=700) ewfc
        endif
        if (inline(i)(1:15).eq.'preequilibrium ') then
          do 190 i2=16,80
            ch=inline(i)(i2:i2)
            if (ch.ne.' ') then
              if (ch.eq.'y') then 
                epreeq=0.
                goto 110
              endif
              if (ch.eq.'n') then 
                epreeq=250.
                goto 110
              endif
              goto 200
            endif
  190     continue
  200     read(inline(i)(16:80),*,err=700) epreeq
        endif
        if (inline(i)(1:11).eq.'multipreeq ') then
          do 210 i2=12,80
            ch=inline(i)(i2:i2)
            if (ch.ne.' ') then
              if (ch.eq.'y') then 
                emulpre=0.
                goto 110
              endif
              if (ch.eq.'n') then 
                emulpre=250.
                goto 110
              endif
              goto 220
            endif
  210     continue
  220     read(inline(i)(12:80),*,err=700) emulpre
        endif
        if (inline(i)(1:10).eq.'preeqspin ') then
          do 230 i2=11,80
            ch=inline(i)(i2:i2)
            if (ch.ne.' ') then
              if (ch.eq.'n') flagpespin=.false.
              if (ch.eq.'y') flagpespin=.true.
              if (ch.ne.'y'.and.ch.ne.'n') goto 700
              goto 110
            endif
  230     continue
        endif
        if (inline(i)(1:15).eq.'giantresonance ') then
          do 240 i2=16,80
            ch=inline(i)(i2:i2)
            if (ch.ne.' ') then
              if (ch.eq.'n') flaggiant0=.false.
              if (ch.eq.'y') flaggiant0=.true.
              if (ch.ne.'y'.and.ch.ne.'n') goto 700
              goto 110
            endif
  240     continue
        endif
        if (inline(i)(1:13).eq.'preeqsurface ') then
          do 250 i2=14,80
            ch=inline(i)(i2:i2)
            if (ch.ne.' ') then
              if (ch.eq.'n') flagsurface=.false.
              if (ch.eq.'y') flagsurface=.true.
              if (ch.ne.'y'.and.ch.ne.'n') goto 700
              goto 110
            endif
  250     continue
        endif
        if (inline(i)(1:13).eq.'preeqcomplex ') then
          do 260 i2=14,80
            ch=inline(i)(i2:i2)
            if (ch.ne.' ') then
              if (ch.eq.'n') flagpecomp=.false.
              if (ch.eq.'y') flagpecomp=.true.
              if (ch.ne.'y'.and.ch.ne.'n') goto 700
              goto 110
            endif
  260     continue
        endif
        if (inline(i)(1:13).eq.'twocomponent ') then
          do 270 i2=14,80
            ch=inline(i)(i2:i2)
            if (ch.ne.' ') then
              if (ch.eq.'n') flag2comp=.false.
              if (ch.eq.'y') flag2comp=.true.
              if (ch.ne.'y'.and.ch.ne.'n') goto 700
              goto 110
            endif
  270     continue
        endif
        if (inline(i)(1:9).eq.'channels ') then
          do 280 i2=10,80
            ch=inline(i)(i2:i2)
            if (ch.ne.' ') then
              if (ch.eq.'n') flagchannels=.false.
              if (ch.eq.'y') flagchannels=.true.
              if (ch.ne.'y'.and.ch.ne.'n') goto 700
              goto 110
            endif
  280     continue
        endif
        if (inline(i)(1:8).eq.'fission ') then
          do 290 i2=9,80
            ch=inline(i)(i2:i2)
            if (ch.ne.' ') then
              if (ch.eq.'n') flagfission=.false.
              if (ch.eq.'y') flagfission=.true.
              if (ch.ne.'y'.and.ch.ne.'n') goto 700
              goto 110
            endif
  290     continue
        endif
        if (inline(i)(1:9).eq.'colldamp ') then
          do 300 i2=10,80
            ch=inline(i)(i2:i2)
            if (ch.ne.' ') then
              if (ch.eq.'n') flagcolldamp=.false.
              if (ch.eq.'y') flagcolldamp=.true.
              if (ch.ne.'y'.and.ch.ne.'n') goto 700
              goto 110
            endif
  300     continue
        endif
        if (inline(i)(1:7).eq.'class2 ') then
          do 310 i2=8,80
            ch=inline(i)(i2:i2)
            if (ch.ne.' ') then
              if (ch.eq.'n') flagclass2=.false.
              if (ch.eq.'y') flagclass2=.true.
              if (ch.ne.'y'.and.ch.ne.'n') goto 700
              goto 110
            endif
  310     continue
        endif
        if (inline(i)(1:9).eq.'outbasic ') then
          do 320 i2=10,80
            ch=inline(i)(i2:i2)
            if (ch.ne.' ') then
              if (ch.eq.'n') flagbasic=.false.
              if (ch.eq.'y') flagbasic=.true.
              if (ch.ne.'y'.and.ch.ne.'n') goto 700
              goto 110
            endif
  320     continue
        endif
        if (inline(i)(1:13).eq.'eciscompound ') then
          do 330 i2=14,80
            ch=inline(i)(i2:i2)
            if (ch.ne.' ') then
              if (ch.eq.'n') flageciscomp=.false.
              if (ch.eq.'y') flageciscomp=.true.
              if (ch.ne.'y'.and.ch.ne.'n') goto 700
              goto 110
            endif
  330     continue
        endif
        if (inline(i)(1:9).eq.'ecisdwba ') then
          do 340 i2=10,80
            ch=inline(i)(i2:i2)
            if (ch.ne.' ') then
              if (ch.eq.'n') flagecisdwba=.false.
              if (ch.eq.'y') flagecisdwba=.true.
              if (ch.ne.'y'.and.ch.ne.'n') goto 700
              goto 110
            endif
  340     continue
        endif
        if (inline(i)(1:8).eq.'onestep ') then
          do 350 i2=9,80
            ch=inline(i)(i2:i2)
            if (ch.ne.' ') then
              if (ch.eq.'n') flagonestep=.false.
              if (ch.eq.'y') flagonestep=.true.
              if (ch.ne.'y'.and.ch.ne.'n') goto 700
              goto 110
            endif
  350     continue
        endif
        if (inline(i)(1:9).eq.'localomp ') then
          do 360 i2=10,80
            ch=inline(i)(i2:i2)
            if (ch.ne.' ') then
              if (ch.eq.'n') flaglocalomp=.false.
              if (ch.eq.'y') flaglocalomp=.true.
              if (ch.ne.'y'.and.ch.ne.'n') goto 700
              goto 110
            endif
  360     continue
        endif
        if (inline(i)(1:10).eq.'optmodall ') then
          do 370 i2=11,80
            ch=inline(i)(i2:i2)
            if (ch.ne.' ') then
              if (ch.eq.'n') flagompall=.false.
              if (ch.eq.'y') flagompall=.true.
              if (ch.ne.'y'.and.ch.ne.'n') goto 700
              goto 110
            endif
  370     continue
        endif
        if (inline(i)(1:9).eq.'statepot ') then
          do 380 i2=10,80
            ch=inline(i)(i2:i2)
            if (ch.ne.' ') then
              if (ch.eq.'n') flagstate=.false.
              if (ch.eq.'y') flagstate=.true.
              if (ch.ne.'y'.and.ch.ne.'n') goto 700
              goto 110
            endif
  380     continue
        endif
        if (inline(i)(1:7).eq.'optmod ') then 
          do 390 i2=8,80
            ch=inline(i)(i2:i2)
            if (ch.ne.' ') then
              read(inline(i)(i2:80),*,err=700) iz,ia
              Zix=Zinit-iz
              Nix=Ninit-ia+iz
              if (Zix.lt.0.or.Zix.gt.numZ.or.Nix.lt.0.or.Nix.gt.numN) 
     +          then
                write(*,'("TALYS-warning: Z,N index out of range,",$)')
                write(*,'(" keyword ignored: ",a80)') inline(i)
              else    
                goto 400
              endif   
            endif
  390     continue
  400     do 410 i2=8,80
            ch=inline(i)(i2:i2)
            if (ch.ne.' '.and.(ch.lt.'0'.or.ch.gt.'9')) then
              do 420 i3=i2,80
                ch=inline(i)(i3:i3)
                if (ch.eq.' ') then
                  ompf=inline(i)(i2:i3-1)
                  do 430 i4=i3,80
                    ch=inline(i)(i4:i4)
                    if (ch.ne.' ') then
                      do 440 type=1,6
                        if (ch.eq.parsym(type)) then
                          optmod(Zix,Nix,type)=ompf
                          goto 110
                        endif
  440                 continue
                    endif
  430             continue
                  optmod(Zix,Nix,1)=ompf
                  goto 110
                endif
  420         continue
            endif
  410     continue
        endif
        if (inline(i)(1:12).eq.'sysreaction ') then
          do 450 type=0,6
            flagsys(type)=.false.
  450     continue
          do 460 i2=13,80
            ch=inline(i)(i2:i2)
            if (ch.ne.' ') then
              do 470 type=0,6
                if (ch.eq.parsym(type)) then
                  flagsys(type)=.true.
                  goto 460
                endif
  470         continue
              goto 700
            endif
  460     continue
          goto 110
        endif
        if (inline(i)(1:11).eq.'rotational ') then
          do 480 type=1,6
            flagrot(type)=.false.
  480     continue
          do 490 i2=12,80
            ch=inline(i)(i2:i2)
            if (ch.ne.' ') then
              do 500 type=1,6
                if (ch.eq.parsym(type)) then
                  flagrot(type)=.true.
                  goto 490
                endif
  500         continue
              goto 700
            endif
  490     continue
          goto 110
        endif
        if (inline(i)(1:5).eq.'asys ') then
          do 510 i2=6,80
            ch=inline(i)(i2:i2)
            if (ch.ne.' ') then
              if (ch.eq.'n') flagasys=.false.
              if (ch.eq.'y') flagasys=.true.
              if (ch.ne.'y'.and.ch.ne.'n') goto 700
              goto 110
            endif
  510     continue
        endif
        if (inline(i)(1:7).eq.'gshell ') then
          do 520 i2=8,80
            ch=inline(i)(i2:i2)
            if (ch.ne.' ') then
              if (ch.eq.'n') flaggshell=.false.
              if (ch.eq.'y') flaggshell=.true.
              if (ch.ne.'y'.and.ch.ne.'n') goto 700
              goto 110
            endif
  520     continue
        endif
        if (inline(i)(1:8).eq.'massdis ') then
          do 530 i2=9,80
            ch=inline(i)(i2:i2)
            if (ch.ne.' ') then
              if (ch.eq.'n') flagmassdis=.false.
              if (ch.eq.'y') flagmassdis=.true.
              if (ch.ne.'y'.and.ch.ne.'n') goto 700
              goto 110
            endif
  530     continue
        endif
        if (inline(i)(1:14).eq.'ffevaporation ') then
          do 540 i2=15,80
            ch=inline(i)(i2:i2)
            if (ch.ne.' ') then
              if (ch.eq.'n') flagffevap=.false.
              if (ch.eq.'y') flagffevap=.true.
              if (ch.ne.'y'.and.ch.ne.'n') goto 700
              goto 110
            endif
  540     continue
        endif
        if (inline(i)(1:5).eq.'endf ') then
          do 550 i2=6,80
            ch=inline(i)(i2:i2)
            if (ch.ne.' ') then
              if (ch.eq.'n') flagendf=.false.
              if (ch.eq.'y') flagendf=.true.
              if (ch.ne.'y'.and.ch.ne.'n') goto 700
              goto 110
            endif
  550     continue
        endif            
        if (inline(i)(1:7).eq.'recoil ') then
          do 560 i2=8,80
            ch=inline(i)(i2:i2)
            if (ch.ne.' ') then
              if (ch.eq.'n') flagrecoil=.false.
              if (ch.eq.'y') flagrecoil=.true.
              if (ch.ne.'y'.and.ch.ne.'n') goto 700
              goto 110
            endif
  560     continue
        endif
        if (inline(i)(1:7).eq.'labddx ') then
          do 580 i2=8,80
            ch=inline(i)(i2:i2)
            if (ch.ne.' ') then
              if (ch.eq.'n') flaglabddx=.false.
              if (ch.eq.'y') flaglabddx=.true.
              if (ch.ne.'y'.and.ch.ne.'n') goto 700
              goto 110
            endif
  580     continue
        endif
        if (inline(i)(1:14).eq.'recoilaverage ') then
          do 590 i2=15,80
            ch=inline(i)(i2:i2)
            if (ch.ne.' ') then
              if (ch.eq.'n') flagrecoilav=.false.
              if (ch.eq.'y') flagrecoilav=.true.
              if (ch.ne.'y'.and.ch.ne.'n') goto 700
              goto 110
            endif
  590     continue
        endif
        if (inline(i)(1:14).eq.'channelenergy ') then
          do 600 i2=15,80
            ch=inline(i)(i2:i2)
            if (ch.ne.' ') then
              if (ch.eq.'n') flagEchannel=.false.
              if (ch.eq.'y') flagEchannel=.true.
              if (ch.ne.'y'.and.ch.ne.'n') goto 700
              goto 110
            endif
  600     continue
        endif
        if (inline(i)(1:7).eq.'fullhf ') then
          do 610 i2=8,80
            ch=inline(i)(i2:i2)
            if (ch.ne.' ') then
              if (ch.eq.'n') flagfullhf=.false.
              if (ch.eq.'y') flagfullhf=.true.
              if (ch.ne.'y'.and.ch.ne.'n') goto 700
              goto 110
            endif
  610     continue
        endif
        if (inline(i)(1:8).eq.'autorot ') then
          do 620 i2=9,80
            ch=inline(i)(i2:i2)
            if (ch.ne.' ') then
              if (ch.eq.'n') flagautorot=.false.
              if (ch.eq.'y') flagautorot=.true.
              if (ch.ne.'y'.and.ch.ne.'n') goto 700
              goto 110
            endif
  620     continue
        endif
        goto 110
  700   write(*,'("TALYS-error: Wrong input: ",a80)') inline(i)
        stop
  110 continue
      return
      end
Copyright (C) 2004  A.J. Koning, S. Hilaire and M.C. Duijvestijn
