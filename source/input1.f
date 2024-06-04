      subroutine input1
c
c +---------------------------------------------------------------------
c | Author: Arjan Koning 
c | Date  : June 28, 2004
c | Task  : Read input for first set of variables
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      logical     projexist,massexist,elemexist,enerexist,lexist
      character*1 ch
      integer     i,i2,type,iz
      real        Ein
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
c
c 2. The projectile is read
c
c nlines: number of input lines
c ch    : character
c inline: input line                 
c ptype0: type of incident particle
c
      do 10 i=1,nlines
        if (inline(i)(1:11).eq.'projectile ') then
          projexist=.true.
          do 20 i2=12,80
            ch=inline(i)(i2:i2)
            if (ch.ne.' ') then 
              ptype0=inline(i)(i2:i2)
              goto 10
            endif
   20     continue
        endif
c
c 3. The target mass is read
c
c Atarget: mass number of target nucleus
c
        if (inline(i)(1:5).eq.'mass ') then
          massexist=.true.
          read(inline(i)(6:80),*,err=400) Atarget
          goto 10
        endif
c
c 4. The nuclear symbol or charge number is read
c
c numelem: number of elements   
c Starget: symbol of target nucleus
c
        if (inline(i)(1:8).eq.'element ') then
          elemexist=.true.
          do 30 i2=9,80
            ch=inline(i)(i2:i2)
            if (ch.ne.' ') then 
              if (ch.ge.'0'.and.ch.le.'9') then 
                read(inline(i)(i2:80),*,err=400) Ztarget 
                if (Ztarget.lt.1.or.Ztarget.gt.numelem) goto 400
                goto 10
              else
                read(inline(i)(i2:i2+1),'(a2)',err=400) Starget
                Starget(1:1)=char(ichar(Starget(1:1))-32)
                goto 10
              endif
            endif
   30     continue
        endif
c
c 5. The incident energy or file with incident energies is read
c
        if (inline(i)(1:7).eq.'energy ') then
          enerexist=.true.
          do 40 i2=8,80
            ch=inline(i)(i2:i2)
            if (ch.ne.' ') then 
              if ((ch.ge.'0'.and.ch.le.'9').or.ch.eq.'.') then 
                read(inline(i)(i2:80),*,err=400) eninc(1)
                goto 10
              else
                eninc(1)=0.
                read(inline(i)(i2:80),'(a73)') energyfile
                goto 10
              endif
            endif
   40     continue
        endif
   10 continue
c
c The four main keywords MUST be present in the input file.
c
      if (.not.projexist) then
        write(*,'("TALYS-error: projectile must be given")') 
        stop
      endif
      if (.not.massexist) then
        write(*,'("TALYS-error: mass must be given")') 
        stop
      endif
      if (.not.elemexist) then
        write(*,'("TALYS-error: element must be given")') 
        stop
      endif
      if (.not.enerexist) then
        write(*,'("TALYS-error: energy must be given")') 
        stop
      endif
c
c ************* Process first set of input variables *******************
c
c 1. Assignment of index k0 to incident particle 
c
c parsym: symbol of particle
c k0    : index of incident particle
c
c Throughout TALYS, the initial particle can always be identified
c by the index k0.
c
      do 110 type=0,6
        if (ptype0.eq.parsym(type)) then
          k0=type
          goto 200
        endif
  110 continue
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
c nin: counter for incident energies
c Ein: incident energy
c
c A. If no incident energy is given in the input file, incident energies
c    should be read from a file.
c
      if (eninc(1).eq.0.) then
        inquire (file=energyfile,exist=lexist)
        if (.not.lexist) then   
          write(*,'("TALYS-error: give a single incident energy ",$)')
          write(*,'("in the input file using the energy keyword ")')
          write(*,'(13x,"or specify a range of incident energies ",$)')
          write(*,'("in a file ",a73)') energyfile
          stop
        endif
        nin=0
        eninc(0)=0.
        open (unit=2,status='old',file=energyfile)
 310    read(2,*,end=320,err=410) Ein
        if (Ein.ne.0.) then
          nin=nin+1
c
c There is a maximum number of incident energies
c
c numenin : maximal number of incident energies 
c
          if (nin.gt.numenin) then
            write(*,'("TALYS-error: there are more than",i3,$)') 
     +        numenin
            write(*,'(" incident energies in file ",a73)') energyfile
            write(*,'(" numenin in talys.cmb should be increased")')
            stop
          endif
          eninc(nin)=Ein
c
c Incident energies must be given in ascending order
c
          if (eninc(nin).le.eninc(nin-1)) then
            write(*,'("TALYS-error: incident energies must",$)') 
            write(*,'(" be given in ascending order")')
            stop
          endif
        endif
        goto 310
 320    close (unit=2)
c
c The number of incident energies is numinc
c
c numinc: number of incident energies 
c
        numinc=nin
        if (numinc.eq.0) then
          write(*,'("TALYS-error: there are no",$)') 
          write(*,'(" incident energies in file ",a73)') energyfile
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
 330    continue
      else
c
c B. Single value given in the user input file
c
        numinc=1
        enincmin=eninc(1)
        enincmax=eninc(1)
      endif
      return
  400 write(*,'("TALYS-error: Wrong input: ",a80)') inline(i)
      stop
  410 write(*,'("TALYS-error: Wrong energy in file ",a73)') energyfile
      write(*,'(" after E= ",e12.5)') Ein
      stop
      end
Copyright (C) 2004  A.J. Koning, S. Hilaire and M.C. Duijvestijn
