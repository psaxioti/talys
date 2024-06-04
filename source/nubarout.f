      subroutine nubarout
c
c +---------------------------------------------------------------------
c | Author: Arjan Koning
c | Date  : October 10, 2015
c | Task  : Output of number of fission neutrons and spectra
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      character*90 nufile
      integer      i,ia,type,nen
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
c
c nu per number, P(nu) and nubar as function of Z and A
c
      if (nin0.eq.numinclow+1.and.numinclow.gt.0) then
        do nen=1,numinclow
          do type=0,6
            fnubar(nen,type)=nubar(type)
          enddo
        enddo
      endif
      write(*,'(/" +++ AVERAGE NUMBER OF PROMPT FISSION NEUTRONS +++")')
      do 10 type=0,6
        if (nubar(type).eq.0.) goto 10
        write(*,'(/" nubar for ",a8,f10.5)') parname(type),nubar(type)
        nufile='Pnux000.000.fis'//natstring(iso)
        nufile(4:4)=parsym(type)   
        if (Einc0.lt.0.001) then
          write(nufile(5:11),'(1p,e7.1)') Einc0
        else
          write(nufile(5:11),'(f7.3)') Einc0
          write(nufile(5:7),'(i3.3)') int(Einc0)
        endif
        open (unit=1,status='unknown',file=nufile)
        write(1,'("# ",a1," + ",i3,a2,": Prompt ",a8," multiplicity ",
     +    "distribution ")') parsym(k0),Atarget,Starget,parname(type)
        write(1,'("# E-incident = ",f7.3," MeV")') Einc0
        write(1,'("# Mean value (nubar-prompt) = ",f10.5)') nubar(type)
        write(1,'("# ")')
        write(*,'(/" P(nu) for ",a8/)') parname(type)
        write(*,'(" number    nu ")')
        write(1,'("# number   nu ")')
        do 20 i=1,numnu
          if (Pdisnu(type,i).gt.0.) then
            write(*,'(i3,1p,e15.4)') i,Pdisnu(type,i)
            write(1,'(i3,1p,e15.4)') i,Pdisnu(type,i)
          endif
  20    continue
        close (unit=1)
        nufile='nuxA000.000.fis'//natstring(iso)
        nufile(3:3)=parsym(type)   
        if (Einc0.lt.0.001) then
          write(nufile(5:11),'(1p,e7.1)') Einc0
        else
          write(nufile(5:11),'(f7.3)') Einc0
          write(nufile(5:7),'(i3.3)') int(Einc0)
        endif
        open (unit=1,status='unknown',file=nufile)
        write(1,'("# ",a1," + ",i3,a2,": Average prompt ",a8,
     +    " multiplicity as function of mass")') parsym(k0),Atarget,
     +    Starget,parname(type)
        write(1,'("# E-incident = ",f7.3," MeV")') Einc0
        write(1,'("# Mean value (nubar-prompt) = ",f10.5)') nubar(type)
        write(1,'("# ")')
        write(*,'(/"  nu(A)"/)')
        write(*,'("  A      Pre-n nu       Post-n nu")')
        write(1,'("# A      Pre-n nu       Post-n nu")')
        do 30 ia=1,Atarget
          write(*,'(i3,1p,2e15.4)') ia,nupre(type,ia),nupost(type,ia)
          write(1,'(i3,1p,2e15.4)') ia,nupre(type,ia),nupost(type,ia)
  30    continue
        close (unit=1)
c
c Write nubar
c
c nubarexist : flag for existence of nubar file
c numinc     : number of incident energies
c eninc      : incident energy in MeV
c
        nufile='nubarx.tot'//natstring(iso)//'       '
        nufile(6:6)=parsym(type)   
        if (.not.nubarexist(type)) then
          nubarexist(type)=.true.
          open (unit=1,status='unknown',file=nufile)
          write(1,'("# ",a1," + ",i3,a2,
     +      ": Average prompt ",a8," multiplicity (nubar-prompt)")')
     +    parsym(k0),Atarget,Starget,parname(type)
          write(1,'("# ")')
          write(1,'("# # energies =",i6)') numinc
          write(1,'("# ")')
          write(1,'("# E-in           nubar")')
          do 40 nen=1,numinclow
            write(1,'(1p,2e12.5)') eninc(nen),fnubar(nen,type)
   40     continue
          do 50 nen=numinclow+1,nin0-1
            write(1,'(1p,2e12.5)') eninc(nen),0.
   50     continue
        else
          open (unit=1,status='old',file=nufile)
          do 60 nen=1,nin0+4
            read(1,*,end=70,err=70)
   60     continue
        endif
        write(1,'(1p,2e12.5)') Einc0,nubar(type)
   70   close (unit=1)
c
c Write PFNS, PFGS, etc.
c
        if (fymodel.le.2) goto 10
        write(*,'(/" +++ Compound fission ",a8," spectrum +++")')
     +    parname(type)
        nufile='cfxs000.000.fis'//natstring(iso)
        nufile(3:3)=parsym(type)   
        if (Einc0.lt.0.001) then
          write(nufile(5:11),'(1p,e7.1)') Einc0
        else
          write(nufile(5:11),'(f7.3)') Einc0
          write(nufile(5:7),'(i3.3)') int(Einc0)
        endif
        open (unit=1,status='unknown',file=nufile)
        write(1,'("# ",a1," + ",i3,a2,": Compound fission ",a8,
     +    " spectrum ")') parsym(k0),Atarget,Starget,parname(type)
        write(1,'("# E-incident = ",f7.3," MeV")') Einc0
        write(1,'("# ")') 
        write(1,'("# E-average  = ",f7.3," MeV")') Eavpfns(type)
ctest   write(*,'(/"       E-out         spectrum    Maxwell ratio"/)')
ctest   write(1,'("#      E-out          spectrum    Maxwell ratio")')
        write(*,'(/"       E-out         spectrum"/)')
        write(1,'("#      E-out          spectrum")')
        do 110 nen=1,numen2
          if (pfns(type,nen).gt.0.) then
ctest       write(*,'(1p,3e15.4)') espec(type,nen),pfns(type,nen),
ctest+        maxpfns(type,nen)
ctest       write(1,'(1p,3e15.4)') espec(type,nen),pfns(type,nen),
ctest+        maxpfns(type,nen)
            write(*,'(1p,2e15.4)') espec(type,nen),pfns(type,nen)
            write(1,'(1p,2e15.4)') espec(type,nen),pfns(type,nen)
          endif
 110    continue
        close (unit=1)
   10 continue
      return
      end
Copyright (C)  2013 A.J. Koning, S. Hilaire and S. Goriely