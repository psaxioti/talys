      subroutine nudisout
c
c +---------------------------------------------------------------------
c | Author: Arjan Koning
c | Date  : February 28, 2022
c | Task  : Output of number of fission neutrons
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      character*8   Estring
      character*132 nufile,Eline
      integer       i,ia,iz,in,type
c
c Write results to separate files
c
c nufile    : file for nubar
c yieldfile : file with fission yields
c natstring : string extension for file names
c iso       : counter for isotope
c Einc0     : incident energy in MeV
c parsym    : symbol of particle
c k0        : index of incident particle
c Atarget   : mass number of target nucleus
c Ztarget   : charge number of target nucleus
c
c nu per number, P(nu) and nubar as function of Z and A
c
      Estring='0000.000'
      if (Einc0.lt.0.001) then
        write(Estring(1:8),'(es8.2)') Einc0
      else
        write(Estring(1:8),'(f8.3)') Einc0
        write(Estring(1:4),'(i4.4)') int(Einc0)
      endif
      Eline='# E-incident =          MeV'
      if (Einc0.lt.0.001) then
        write(Eline(16:23),'(es8.2)') Einc0
      else
        write(Eline(16:23),'(f8.3)') Einc0
      endif
      write(*,'(/" +++ NUMBER OF PROMPT FISSION NEUTRONS ",
     +  "AND GAMMAS +++")')
      do 10 type=0,6
        if (parskip(type)) goto 10
        if (nubar(type).eq.0.) goto 10
c
c P(nu)
c
        write(*,'(/" nubar for ",a8,f10.5)') parname(type),nubar(type)
        nufile='Pnux'//Estring//'.fis'//natstring(iso)
        nufile(4:4)=parsym(type)
        open (unit=1,file=nufile,status='replace')
        write(1,'("# ",a1," + ",a,": Prompt ",a8," multiplicity ",
     +    "distribution ")') parsym(k0),trim(targetnuclide),
     +    parname(type)
        write(1,'(a)') trim(Eline)
        write(1,'("# Mean value (nubar-prompt) = ",f10.5)') nubar(type)
        write(1,'("# Average P(nu)             = ",f10.5)') 
     +    Pdisnuav(type)
        write(*,'(/" P(nu) for ",a8/)') parname(type)
        write(*,'(" number    nu ")')
        write(1,'("# number   nu ")')
        do 20 i=0,numnu
          if (Pdisnu(type,i).gt.0..or.i.le.4) then
            write(*,'(i3,es15.4)') i,Pdisnu(type,i)
            write(1,'(i3,es15.4)') i,Pdisnu(type,i)
          endif
  20    continue
        close (unit=1)
        write(*,'(/" Average P(nu) for ",a8,f10.5/)') parname(type),
     +    Pdisnuav(type)
c
c nu(A)
c
        nufile='nuxA'//Estring//'.fis'//natstring(iso)
        nufile(3:3)=parsym(type)
        open (unit=1,file=nufile,status='replace')
        write(1,'("# ",a1," + ",a,": Average prompt ",a8,
     +    " multiplicity as function of mass")') parsym(k0),
     +    trim(targetnuclide),parname(type)
        write(1,'(a)') trim(Eline)
        write(1,'("# Mean value (nubar-prompt) = ",f10.5)') nubar(type)
        write(1,'("# ")')
        write(*,'(/"  nu(A) for ",a8/)') parname(type)
        write(*,'("  A        nu   ")')
        write(1,'("# A        nu   ")')
        do 30 ia=1,Atarget
          write(*,'(i3,es17.6)') ia,nuA(type,ia)
          write(1,'(i3,es17.6)') ia,nuA(type,ia)
  30    continue
        close (unit=1)
c
c nu(Z,A)
c
        nufile='nuxZA'//Estring//'.fis'//natstring(iso)
        nufile(3:3)=parsym(type)
        open (unit=1,file=nufile,status='replace')
        write(1,'("# ",a1," + ",a,": Average prompt ",a8,
     +    " multiplicity as function of nucleus")') parsym(k0),
     +    trim(targetnuclide),parname(type)
        write(1,'(a)') trim(Eline)
        write(1,'("# Mean value (nubar-prompt) = ",f10.5)') nubar(type)
        write(1,'("# ")')
        write(*,'(/"  nu(Z,A) for ",a8/)') parname(type)
        write(*,'("  Z   A       nu")')
        write(1,'("# Z   A       nu")')
        do 35 iz=1,Ztarget
          do 35 ia=iz+1,Atarget
            in=ia-iz
            if (in.le.numneu) then
              if (nuZA(type,iz,in).gt.0.) then
                write(*,'(i3,i4,es17.6)') iz,ia,nuZA(type,iz,in)
                write(1,'(i3,i4,es17.6)') iz,ia,nuZA(type,iz,in)
              endif
            endif
   35   continue
        close (unit=1)
   10 continue
      if (fymodel <=2) return
c
c E-average(Z,A) and E-average(A)
c
      nufile='EavZA'//Estring//'.fis'//natstring(iso)
      open (unit=1,file=nufile,status='replace')
      write(1,'("# ",a1," + ",a,": Average emission ",
     +  "energy per (Z,A)")') parsym(k0),trim(targetnuclide)
      write(1,'(a)') trim(Eline)
      write(1,'("# ")')
      write(1,'("# ")')
      write(1,'("# Z  A      gamma    neutron")')
      write(*,'(/"  Average emission energy per (Z,A)")')
      write(*,'(/"  Z  A      gamma    neutron"/)')
      do 110 iz=1,Ztarget
        do 110 ia=iz+1,Atarget
          in=ia-iz
          if (in.le.numneu) then
            if (EaverageZA(0,iz,in).gt.0..or.EaverageZA(1,iz,in).gt.0.)
     +        then
                write(*,'(i3,i4,2es12.5)') iz,ia,
     +            (EaverageZA(type,iz,in),type=0,1)
                write(1,'(i3,i4,2es12.5)') iz,ia,
     +            (EaverageZA(type,iz,in),type=0,1)
            endif
          endif
  110 continue
      close (unit=1)
      nufile='EavA'//Estring//'.fis'//natstring(iso)
      open (unit=1,file=nufile,status='replace')
      write(1,'("# ",a1," + ",a,": Average emission ",
     +  "energy per A")') parsym(k0),trim(targetnuclide)
      write(1,'(a)') trim(Eline)
      write(1,'("# ")')
      write(1,'("# ")')
      write(1,'("# A     gamma    neutron")')
      write(*,'(/"  Average emission energy per A")')
      write(*,'(/"  A     gamma    neutron"/)')
      do 120 ia=1,Atarget
        write(*,'(i3,2es12.5)') ia,(EaverageA(type,ia),type=0,1)
        write(1,'(i3,2es12.5)') ia,(EaverageA(type,ia),type=0,1)
  120 continue
      close (unit=1)
      return
      end
Copyright (C)  2021 A.J. Koning, S. Hilaire and S. Goriely
