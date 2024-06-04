      subroutine massdisout
c
c +---------------------------------------------------------------------
c | Author: Arjan Koning
c | Date  : October 18, 2015
c | Task  : Output of fission fragment yields
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      character*1  ftype(2)
      character*90 yieldfile,fpfile
      integer      i,iz,ia,in,nen
      real         fiseps
c
c ****************** Output of fission yields **************************
c
      write(*,'(/" ++++++++++ FISSION YIELDS ++++++++++"/)')
      ftype(1)='A'
      ftype(2)='N'
      fiseps=Rfiseps*xsfistot
      do i=1,2
c
c Write results to separate files
c
c yieldfile : file with fission yields
c natstring : string extension for file names
c iso       : counter for isotope
c Einc0     : incident energy in MeV
c parsym    : symbol of particle
c k0        : index of incident particle
c Atarget   : mass number of target nucleus
c Starget   : symbol of target nucleus
c Ztarget   : charge number of target nucleus
c yieldApre : pre-neutron emission mass yield
c yieldApost: post-neutron emission corrected mass yield
c
        yieldfile='yield'//ftype(i)//'000.000.fis'//natstring(iso)
        if (Einc0.lt.0.001) then
          write(yieldfile(7:13),'(1p,e7.1)') Einc0
        else
          write(yieldfile(7:13),'(f7.3)') Einc0
          write(yieldfile(7:9),'(i3.3)') int(Einc0)
        endif
        open (unit=1,status='unknown',file=yieldfile)
        write(1,'("# ",a1," + ",i3,a2,": fission yields")')
     +    parsym(k0),Atarget,Starget
        write(1,'("# E-incident = ",1p,e12.5)') Einc0
        write(1,'("# ")')
        write(1,'("# ")')
        if (i.eq.1) then
          write(*,'(" Fission yields as function of A"/)')
          write(*,'("  ",a1,"    Pre-n yield   Post-n yield ",
     +      "      Pre-n xs      Post-n xs")') ftype(i)
          write(1,'("# ",a1,"    Pre-n yield   Post-n yield ",
     +      "      Pre-n xs      Post-n xs")') ftype(i)
          do 10 ia=1,Atarget
            write(*,'(i3,1p,4e15.4)') ia,yieldApre(ia),yieldApost(ia),
     +        xsfpApre(ia),xsfpApost(ia)
            write(1,'(i3,1p,4e15.4)') ia,yieldApre(ia),yieldApost(ia),
     +        xsfpApre(ia),xsfpApost(ia)
  10      continue
          write(*,'(/"Tot",1p,4e15.4)') yieldtotpre,yieldtotpost,
     +      xsfptotpre,xsfptotpost
        else
          write(*,'(/" Fission yields as function of N")')
         write(*,'(/"  ",a1,"    Pre-n yield   Post-n yield")') ftype(i)
         write(1,'("# ",a1,"    Pre-n yield   Post-n yield")') ftype(i)
          do 20 in=1,Ntarget
            write(*,'(i3,1p,2e15.4)') in,yieldnpre(in),yieldnpost(in)
            write(1,'(i3,1p,2e15.4)') in,yieldnpre(in),yieldnpost(in)
  20      continue
          write(*,'(/"Tot",1p,2e15.4)') yieldtotpre,yieldtotpost
        endif
        close (unit=1)
      enddo
c
c Write ff/fp production
c
c fpexist    : flag for existence of fission product
c fpfile     : file with fission product
c yieldZApre : pre-neutron emission isotopic yield
c yieldZApost: post-neutron emission corrected isotopic yield
c
      write(*,'(/" FF/FP production"/)')
      write(*,'("    Z    A    Pre-n yield   Post-n yield ",
     +  "      Pre-n xs      Post-n xs")')
      do 210 ia=1,Atarget
        do 220 iz=1,Ztarget
          in=ia-iz
          if (in.lt.1.or.in.gt.numneu) goto 220
          if (xsfpZApre(iz,in).lt.fiseps.and.
     +      xsfpZApost(iz,in).lt.fiseps.and..not.fpexist(iz,in))
     +      goto 220
          fpfile='fp000000.tot'//natstring(iso)
          write(fpfile(3:8),'(2i3.3)') iz,ia
          if (.not.fpexist(iz,in)) then
            fpexist(iz,in)=.true.
            open (unit=1,status='unknown',file=fpfile)
            write(1,'("# ",a1," + ",i3,a2,": Fission product yield of ",
     +        i3,a2)') parsym(k0),Atarget,Starget,ia,nuc(iz)
            write(1,'("# ")')
            write(1,'("# # energies =",i6)') numinc
            write(1,'("# ")')
            write(1,'("# E-incident Pre-n yield Post-n yield "
     +        "  Pre-n xs   Post-n xs")')
            do 230 nen=1,nin0-1
              write(1,'(1p,5e12.5)') eninc(nen),0.,0.,0.,0.
  230       continue
          else
            open (unit=1,status='old',file=fpfile)
            do 240 nen=1,nin0+4
              read(1,*,end=250,err=250)
  240       continue
          endif
          write(*,'(2i5,1p,4e15.4)') iz,ia,yieldZApre(iz,in),
     +      yieldZApost(iz,in),xsfpZApre(iz,in),xsfpZApost(iz,in)
          write(1,'(1p,5e12.5)') Einc0,yieldZApre(iz,in),
     +      yieldZApost(iz,in),xsfpZApre(iz,in),xsfpZApost(iz,in)
  250     close (unit=1)
  220   continue
c
c Write cumulative ff/fp production
c
        if (xsfpApre(ia).lt.fiseps.and.
     +    xsfpApost(ia).lt.fiseps.and..not.fpaexist(ia)) goto 210
        fpfile='fp000000.tot'//natstring(iso)
        write(fpfile(6:8),'(i3.3)') ia
        if (.not.fpaexist(ia)) then
          fpaexist(ia)=.true.
          open (unit=1,status='unknown',file=fpfile)
          write(1,'("# ",a1," + ",i3,a2,": Fission product yield of A=",
     +      i3)') parsym(k0),Atarget,Starget,ia
          write(1,'("# ")')
          write(1,'("# # energies =",i6)') numinc
          write(1,'("# ")')
          write(1,'("# E-incident Pre-n yield Post-n yield",
     +      "   Pre-n xs   Post-n xs")')
          do 310 nen=1,nin0-1
            write(1,'(1p,5e12.5)') eninc(nen),0.,0.,0.,0.
  310     continue
        else
          open (unit=1,status='old',file=fpfile)
          do 320 nen=1,nin0+4
            read(1,*,end=330,err=330)
  320     continue
        endif
        write(1,'(1p,5e12.5)') Einc0,yieldApre(ia),yieldApost(ia),
     +    xsfpApre(ia),xsfpApost(ia)
  330   close (unit=1)
  210 continue
      write(*,'(/"Total     ",1p,4e15.4)') yieldtotpre,yieldtotpost,
     +  xsfptotpre,xsfptotpost
      write(*,'(/" Cumulative FF/FP production"/)')
      write(*,'("    A    Pre-n yield   Post-n yield ",
     +  "      Pre-n xs      Post-n xs")')
      do 410 ia=1,Atarget
        write(*,'(i5,1p,4e15.4)') ia,yieldApre(ia),yieldApost(ia),
     +    xsfpApre(ia),xsfpApost(ia)
  410 continue
      write(*,'(/"Total",1p,4e15.4)') yieldtotpre,yieldtotpost,
     +  xsfptotpre,xsfptotpost
      return
      end
Copyright (C)  2013 A.J. Koning, S. Hilaire and S. Goriely
