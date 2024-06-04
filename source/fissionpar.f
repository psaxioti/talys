      subroutine fissionpar(Zix,Nix)
c
c +---------------------------------------------------------------------
c | Author: Stephane Hilaire, Marieke Duijvestijn and Arjan Koning
c | Date  : October 14, 2004
c | Task  : Fission parameters
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      logical      lexist
      character*4  fischar
      character*90 fisfile,hbsfile,c2file
      integer      Zix,Nix,fislocal,Z,N,A,ia,n1,i,j,n2,il,modz,modn
      real         bar1,bar2,hw1,hw2,egs,lbar0,esp
c
c ****************** Read fission barrier parameters *******************
c
c Note that next to the chosen fission model (fismodel), there is always
c an alternative fission model (fismodelalt) which comes into play if 
c fission parameters for the first choice model are not available.

c Determine whether barrier parameters have been provided in input
c
c Zix     : charge number index for residual nucleus
c Nix     : neutron number index for residual nucleus
c nfisbar : number of fission barrier parameters
c fbarrier: height of fission barrier
c fwidth  : width of fission barrier
c fislocal: fission model
c fismodel: fission model
c
      nfisbar(Zix,Nix)=0
      do 10 i=1,numbar
        if (fbarrier(Zix,Nix,i).ne.0..or.fwidth(Zix,Nix,i).ne.0.) 
     +    nfisbar(Zix,Nix)=nfisbar(Zix,Nix)+1
   10 continue
      fislocal=fismodel
c
c Fission parameters from database
c
c ZZ,Z   : charge number of residual nucleus
c NN,N   : neutron number of residual nucleus
c AA,A   : mass number of residual nucleus
c fischar: help variable
c
      Z=ZZ(Zix,Nix,0)
      N=NN(Zix,Nix,0)
      A=AA(Zix,Nix,0)
      fischar='z   '
      write(fischar(2:4),'(i3.3)') Z
c
c Fismodel 1: Experimental parameters
c
c fisfile  : fission file
c path     : directory containing structure files to be read
c lenpath  : length of pathname
c axtype   : type of axiality of barrier (1: axial, 2: tri-axial)
c bar1,bar2: inner and outer barrier heights
c hw1,hw2  : inner and outer barrier curvatures
c deltaW   : shell correction in nuclear mass  
c
c Fission barriers from database may have been overruled by user input.
c As starting point, we take the RIPL values as compiled by Maslov.
c
      if (fislocal.eq.1) then
        fisfile=path(1:lenpath)//'fission/barrier/'//fischar
        inquire (file=fisfile,exist=lexist)
        if (.not.lexist) goto 100
        open (unit=2,status='old',file=fisfile)
  110   read(2,'(4x,i4,4x,1x,2(f8.2),5x,2(f8.2))',end=100)
     +    ia,bar1,hw1,bar2,hw2
        if (A.ne.ia) goto 110
        if (axtype(Zix,Nix,1).eq.2.and.deltaW(Zix,Nix,1).eq.0.) 
     +    deltaW(Zix,Nix,1)=2.5
        if (axtype(Zix,Nix,1).eq.1.and.deltaW(Zix,Nix,1).eq.0.) 
     +    deltaW(Zix,Nix,1)=1.5
        if (fbarrier(Zix,Nix,1).eq.0.) fbarrier(Zix,Nix,1)=bar1
        if (fwidth(Zix,Nix,1).eq.0.) fwidth(Zix,Nix,1)=hw1
        if (fbarrier(Zix,Nix,2).eq.0.) fbarrier(Zix,Nix,2)=bar2
        if (fwidth(Zix,Nix,2).eq.0.) fwidth(Zix,Nix,2)=hw2
        nfisbar(Zix,Nix)=2
        close (unit=2)
c
c Read fission states
c
c modz,modn: help variables
c hbsfile  : file with head band transition states
c c2file   : file with class 2 states
c
        modz=mod(Z,2)
        modn=mod(N,2)
        if (modz.eq.0) then
          if (modn.eq.0) then
            hbsfile=path(1:lenpath)//'fission/states/hbstates.ee'
            c2file=path(1:lenpath)//'fission/states/class2states.ee'
          else
            hbsfile=path(1:lenpath)//'fission/states/hbstates.eo'
            c2file=path(1:lenpath)//'fission/states/class2states.eo'
          endif
        else
          if (modn.eq.0) then
            hbsfile=path(1:lenpath)//'fission/states/hbstates.oe'
            c2file=path(1:lenpath)//'fission/states/class2states.oe'
          else
            hbsfile=path(1:lenpath)//'fission/states/hbstates.oo'
            c2file=path(1:lenpath)//'fission/states/class2states.oo'
          endif
        endif
c
c Use user-defined files for head band and class 2 transition states
c
        if (hbtransfile(Zix,Nix)(1:1).ne.' ') 
     +    hbsfile=hbtransfile(Zix,Nix)
        if (class2file(Zix,Nix)(1:1).ne.' ') 
     +    c2file=class2file(Zix,Nix)
c
c Read head band transition states
c
c n1,n2   : help variables
c nfistrhb: number of head band transition states for barrier
c fecont  : start of continuum energy
c efistrhb: energy of head band transition states 
c jfistrhb: spin of head band transition states 
c pfistrhb: parity of head band transition states 
c
        open (unit=2,status='old',file=hbsfile)
        n1=nfisbar(Zix,Nix)
        do 210 i=1,n1
          read(2,'(4x,i4,f8.3)') nfistrhb(Zix,Nix,i),fecont(Zix,Nix,i)
          do 220 j=1,nfistrhb(Zix,Nix,i)
            read(2,'(4x,f11.6,f6.1,i5)') efistrhb(Zix,Nix,i,j),
     +        jfistrhb(Zix,Nix,i,j),pfistrhb(Zix,Nix,i,j)
 220      continue
 210    continue
        close (unit=2)                       
c
c Class2 states
c
c nclass2 : number of sets of class2 states   
c nfisc2hb: number of class2 states for barrier
c efisc2hb: energy of class2 states 
c jfisc2hb: spin of class2 states 
c pfisc2hb: parity of class2 states 
c
        open (unit=2,status='old',file=c2file)
        n2=1
        nclass2(Zix,Nix)=n2
        do 230 i=1,n2
          read(2,'(4x,i4)') nfisc2hb(Zix,Nix,i)
          do 240 j=1,nfisc2hb(Zix,Nix,i)
            read(2,'(4x,f11.6,f6.1,i5)') efisc2hb(Zix,Nix,i,j),
     +        jfisc2hb(Zix,Nix,i,j),pfisc2hb(Zix,Nix,i,j)
  240     continue
  230   continue
        close (unit=2)
      endif
c
c Fismodel 2: Mamdouh parameters
c
      if (fislocal.eq.2) then
        fisfile=path(1:lenpath)//'fission/mamdouh/'//fischar
        inquire (file=fisfile,exist=lexist)
        if (.not.lexist) goto 100
        open (unit=2,status='old',file=fisfile)
  310   read(2,'(4x,i4,2(24x,f8.2))',end=100) ia,bar1,bar2
        if (A.ne.ia) goto 310
        if (fbarrier(Zix,Nix,1).eq.0.) fbarrier(Zix,Nix,1)=bar1
        if (fbarrier(Zix,Nix,2).eq.0.) fbarrier(Zix,Nix,2)=bar2
        if (fbarrier(Zix,Nix,1).eq.0..or.fbarrier(Zix,Nix,2).eq.0.) then
          nfisbar(Zix,Nix)=1
        else
          nfisbar(Zix,Nix)=2
        endif
        close (unit=2)
      endif
  100 if (fismodel.le.2.and.nfisbar(Zix,Nix).eq.0) fislocal=fismodelalt
c
c Fismodel 3: Sierk
c
c barsierk: subroutine for fission barrier heights, rotating gs energy 
c           and lbar0    
c egs     : rotating ground state energy
c lbar0   : l-value for which barrier height becomes zero 
c il      : angular momentum
c
      if (fislocal.eq.3) then
        il=0
        nfisbar(Zix,Nix)=1
        call barsierk(Z,A,il,bar1,egs,lbar0)
        if (fbarrier(Zix,Nix,1).eq.0.) fbarrier(Zix,Nix,1)=bar1
        if (fwidth(Zix,Nix,1).eq.0.) fwidth(Zix,Nix,1)=0.24
      endif
c
c Fismodel 4: Rotating Liquid Drop Model
c
c rldm: subroutine for saddle point energies, rotating gs energy
c esp : saddle point energy
c
      if (fislocal.eq.4) then
        il=0
        nfisbar(Zix,Nix)=1
        call rldm(Z,A,il,egs,esp)
        if (fbarrier(Zix,Nix,1).eq.0.) fbarrier(Zix,Nix,1)=esp-egs
        if (fwidth(Zix,Nix,1).eq.0.) fwidth(Zix,Nix,1)=0.24
      endif
c
c ************************* Default parameters *************************
c
c minertia  : moment of inertia of fission barrier deformation
c Rtransmom : normalization constant for moment of inertia for
c             transition states
c Irigid    : rigid body value of moment of inertia
c minertc2  : moment of inertia for class2 states
c Rclass2mom: normalization constant for moment of inertia for
c             class 2 states
c
      if (fwidth(Zix,Nix,1).eq.0.) fwidth(Zix,Nix,1)=1.
      if (deltaW(Zix,Nix,1).eq.0..and.Ztarget.gt.84) 
     +  deltaW(Zix,Nix,1)=2.5    
      if (fwidth(Zix,Nix,2).eq.0.) fwidth(Zix,Nix,2)=0.6
      if (deltaW(Zix,Nix,2).eq.0.) deltaW(Zix,Nix,2)=0.6
      if (nfisbar(Zix,Nix).eq.1.and.fbarrier(Zix,Nix,1).eq.0.) then
        fbarrier(Zix,Nix,1)=fbarrier(Zix,Nix,2)
        fwidth(Zix,Nix,1)=fwidth(Zix,Nix,2)
        deltaW(Zix,Nix,1)=deltaW(Zix,Nix,2)
      endif
      do 410 i=1,numbar
        minertia(Zix,Nix,i)=Rtransmom(Zix,Nix,i)*Irigid(Zix,Nix,i)
        minertc2(Zix,Nix,i)=Rclass2mom(Zix,Nix,i)*Irigid(Zix,Nix,i)
  410 continue
c
c ********** Rotational bands on transition and class2 states **********
c
c rotband   : subroutine to build rotational bands on transition states
c flagclass2: flag for class2 states in fission
c rotclass2 : subroutine to build rotational bands on class2 states 
c
      call rotband(Zix,Nix)
      if (flagclass2) call rotclass2(Zix,Nix) 
      return
      end  
Copyright (C) 2004  A.J. Koning, S. Hilaire and M.C. Duijvestijn
