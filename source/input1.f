      subroutine input1
c
c +---------------------------------------------------------------------
c | Author: Arjan Koning 
c | Date  : June 26, 2007
c | Task  : Read input for first set of variables
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      logical      projexist,massexist,elemexist,enerexist,lexist
      character*1  ch
      character*80 word(40),key,value
      integer      i,i2,inull,type,iz,J
      real         Ein
c
c ************ Read first set of variables from input lines ************
c
c energyfile : file with incident energies
c projexist  : logical for existence of projectile
c massexist  : logical for existence of mass
c elemexist  : logical for existence of element
c enerexist  : logical for existence of energy
c flagnatural: flag for calculation of natural element
c eninc      : incident energy in MeV
c Ztarget    : charge number of target nucleus
c Starget    : symbol of target nucleus
c
c 1. Initializations
c
      energyfile='                                                     '
      projexist=.false.
      massexist=.false.
      elemexist=.false.
      enerexist=.false.
      flagnatural=.false.
      eninc(1)=0.
      Ztarget=0
      Starget='  '
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
c 2. The projectile is read
c
c ptype0: type of incident particle
c
        if (key.eq.'projectile') then
          projexist=.true.
          ptype0=ch
          goto 10
        endif
c
c 3. The target mass is read
c
c Atarget: mass number of target nucleus
c
        if (key.eq.'mass') then
          massexist=.true.
          read(value,*,end=400,err=400) Atarget
          goto 10
        endif
c
c 4. The nuclear symbol or charge number is read
c
c numelem: number of elements   
c
        if (key.eq.'element') then
          elemexist=.true.
          if (ch.ge.'0'.and.ch.le.'9') then 
            read(value,*,end=400,err=400) Ztarget 
            if (Ztarget.lt.1.or.Ztarget.gt.numelem) goto 400
            goto 10
          else
            read(value,'(a2)',end=400,err=400) Starget
            Starget(1:1)=char(ichar(Starget(1:1))-32)
            goto 10
          endif
        endif
c
c 5. The incident energy or file with incident energies is read
c
        if (key.eq.'energy') then
          enerexist=.true.
          if ((ch.ge.'0'.and.ch.le.'9').or.ch.eq.'.') then 
            read(value,*,end=400,err=400) eninc(1)
            goto 10
          else
            eninc(1)=0.
            energyfile=value
            goto 10
          endif
        endif
   10 continue
c
c The four main keywords MUST be present in the input file.
c
      if (.not.projexist) then
        write(*,'(" TALYS-error: projectile must be given")') 
        stop
      endif
      if (.not.massexist) then
        write(*,'(" TALYS-error: mass must be given")') 
        stop
      endif
      if (.not.elemexist) then
        write(*,'(" TALYS-error: element must be given")') 
        stop
      endif
      if (.not.enerexist) then
        write(*,'(" TALYS-error: energy must be given")') 
        stop
      endif
c
c Manual input of structure path and null device.
c
c path   : directory containing structure files to be read
c lenpath: length of pathname
c nulldev: null device
c
      do 50 i=1,nlines
        call getkeywords(inline(i),word)
        key=word(1)
        value=word(2)
        ch=word(2)(1:1)
        if (key.eq.'strucpath') then
          lenpath=0
          do 60 i2=11,80
            ch=inline(i)(i2:i2)
            if (ch.ne.' ') then
              lenpath=lenpath+1
              path(lenpath:lenpath)=ch
            endif
   60     continue
        endif
        inull=0
        if (key.eq.'nulldev') then
          nulldev='             '
          do 70 i2=9,80
            ch=inline(i)(i2:i2)
            if (ch.ne.' ') then
              inull=inull+1
              nulldev(inull:inull)=ch
            endif
   70     continue
        endif
   50 continue
c
c Test to check accessibility of structure files and null device
c
      if (path(lenpath:lenpath).ne.'/') then
        lenpath=lenpath+1
        path(lenpath:lenpath)='/'
      endif
      if (lenpath.gt.60) then
        write(*,'(" TALYS-warning: path name should contain 60",
     +    " characters or less")')
      endif
      inquire (file=path(1:lenpath)//'abundance/z001',exist=lexist)
      if (.not.lexist) then
        write(*,'(" TALYS-error: Structure database not installed:",
     +    " change path in machine.f or strucpath keyword",
     +    " in input file")')
        stop
      endif
      if (inull.gt.13) then
        write(*,'(" TALYS-warning: null device should contain 13",
     +    " characters or less")')
      endif
c
c ************* Process first set of input variables *******************
c
c 1. Assignment of index k0 to incident particle 
c
c flaginitpop: flag for initial population distribution
c parsym     : symbol of particle
c k0         : index of incident particle
c
c Throughout TALYS, the initial particle can always be identified
c by the index k0.
c
      flaginitpop=.false.
      do 110 type=0,6
        if (ptype0.eq.parsym(type)) then
          k0=type
          goto 200
        endif
  110 continue
c
c It is also possible to define a population distribution as the
c initial state, through the keyword projectile 0. In that case,
c we assign a photon projectile.
c
      if (ptype0.eq.'0') then
        flaginitpop=.true.
        k0=0
      endif
c
c 2. Identification of target and initial compound nucleus
c
c nuc: symbol of nucleus
c
  200 if (Ztarget.eq.0) then
        do 210 iz=1,numelem
          if (nuc(iz).eq.Starget) then
            Ztarget=iz
            goto 220
          endif
  210   continue
      else
        Starget=nuc(Ztarget)
      endif
c
c A calculation for a natural element is specified by target mass 0
c
c abundance: subroutine for natural abundances
c iso      : counter for isotope
c isotope  : isotope number of residual nucleus
c Ntarget  : neutron number of target nucleus
c Zinit    : charge number of initial compound nucleus
c parZ     : charge number of particle
c Ninit    : neutron number of initial compound nucleus
c parN     : neutron number of particle
c Ainit    : mass number of initial compound nucleus
c
  220 if (Atarget.eq.0) then
        flagnatural=.true.
        if (iso.eq.1) call abundance
        Atarget=isotope(iso)
      endif
      Ntarget=Atarget-Ztarget
      Zinit=Ztarget+parZ(k0)
      Ninit=Ntarget+parN(k0)
      Ainit=Zinit+Ninit
c
c 3. Determine incident energies
c
c Ein: incident energy
c nin: counter for incident energies
c
c 1. Normal case of a projectile with a target nucleus
c
      Ein=0.
      if (.not.flaginitpop) then
c
c A. If no incident energy is given in the input file, incident energies
c    should be read from a file.
c
        eninc(0)=0.
        if (eninc(1).eq.0.) then
          inquire (file=energyfile,exist=lexist)
          if (.not.lexist) then   
            write(*,'(" TALYS-error: give a single incident energy ",
     +        "in the input file using the energy keyword ")')
            write(*,'(14x,"or specify a range of incident energies ",
     +        "in a file ",a73)') energyfile
            stop
          endif
          nin=0
          open (unit=2,status='old',file=energyfile)
 310      read(2,*,end=320,err=410) Ein
          if (Ein.ne.0.) then
            nin=nin+1
c
c There is a maximum number of incident energies
c
c numenin : maximal number of incident energies 
c
            if (nin.gt.numenin) then
              write(*,'(" TALYS-error: there are more than",i3,
     +          " incident energies in file ",a73)') numenin,energyfile
              write(*,'(" numenin in talys.cmb should be increased")')
              stop
            endif
            eninc(nin)=Ein
c
c Incident energies must be given in ascending order
c
            if (eninc(nin).le.eninc(nin-1)) then
              write(*,'(" TALYS-error: incident energies must",
     +          " be given in ascending order")')
              stop
            endif
          endif
          goto 310
 320      close (unit=2)
c
c The number of incident energies is numinc
c
c numinc: number of incident energies 
c
          numinc=nin
          if (numinc.eq.0) then
            write(*,'(" TALYS-error: there are no",
     +        " incident energies in file ",a73)') energyfile
            stop
          endif
c
c The minimum and maximum value of all the incident energies is 
c determined.
c
c enincmin: minimum incident energy 
c enincmax: maximum incident energy 
c
          enincmin=eninc(1)
          enincmax=eninc(1)
          do 330 nin=2,numinc
            enincmin=min(enincmin,eninc(nin))
            enincmax=max(enincmax,eninc(nin))
 330      continue
        else
c
c B. Single value given in the user input file
c
          numinc=1
          enincmin=eninc(1)
          enincmax=eninc(1)
        endif
      else
c
c 2. Population distribution as the initial state
c
c npopbins: number of excitation energy bins for population distribution
c npopJ   : number of spins for population distribution
c numbins : maximal number of continuum excitation energy bins
c numJ    : maximal J-value
c Exdist  : excitation energy of population distribution
c Pdistex : population distribution, spin-independent
c Pdist   : population distribution per spin
c
        inquire (file=energyfile,exist=lexist)
        if (.not.lexist) then   
          write(*,'(" TALYS-error: if projectile 0, specify a range",
     +        " of excitation energies in a file ",a73)') energyfile
          stop
        endif
        open (unit=2,status='old',file=energyfile)
        read(2,*,end=410,err=410) npopbins,npopJ,eninc(1)
        if (npopbins.lt.2.or.npopbins.gt.numbins) then
          write(*,'(" TALYS-error: 2 <= bins <=",i3," in population ",
     +      "distribution file")') numbins
          stop
        endif
        if (npopJ.lt.0.or.npopJ.gt.numJ+1) then
          write(*,'(" TALYS-error: 0 <= number of spins <=",i3,
     +      " + 1 in population distribution file")') numJ
          stop
        endif
        Exdist(0)=0.
        Pdistex(0)=0.
        do 340 J=0,numJ
          Pdist(0,J)=0.
 340    continue
c
c Only excitation energy distribution (no spins)
c 
        if (npopJ.eq.0) then
          do 350 nin=1,npopbins
            read(2,*,end=410,err=410) Exdist(nin),Pdistex(nin)
 350      continue
        else
c
c Spin-dependent excitation energy distribution (no total)
c 
          do 360 nin=1,npopbins
            read(2,*,end=410,err=410) Exdist(nin),
     +        (Pdist(nin,J),J=0,npopJ-1)
 360      continue
        endif
        close (unit=2)
        numinc=1
        enincmin=eninc(1)
        enincmax=eninc(1)
      endif
      return
  400 write(*,'(" TALYS-error: Wrong input: ",a80)') inline(i)
      stop
  410 write(*,'(" TALYS-error: Problem in file ",a73)') energyfile
      write(*,'(" after E= ",e12.5)') Ein
      stop
      end
Copyright (C) 2004  A.J. Koning, S. Hilaire and M.C. Duijvestijn
