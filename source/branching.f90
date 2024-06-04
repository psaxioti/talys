      subroutine branching
c
c +---------------------------------------------------------------------
c | Author: Arjan Koning
c | Date  : January 29, 2023
c | Task  : Best set of branching ratios
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      logical       lexist
      character*10  branchchar
      character*132 line,word(40)
      character*132 branchfile
      integer       Zix,Nix,iz,ia,iword,ilev0,ilev1,nbr,k
      real          sum,bra
c
c ******************** Read best branching ratios **********************
c
c branchfile : branching ratio file
c branchchar : part of branching ratio file
c Ztarget    : charge number of target nucleus
c Atarget    : mass number of target nucleus
c Zix        : charge number index for residual nucleus
c Nix        : neutron number index for residual nucleus
c Zinit      : charge number of initial compound nucleus
c Ninit      : neutron number of initial compound nucleus
c branchratio: gamma-ray branching ratio to level
c bra        : gamma-ray branching ratio to level
c branchlevel: level to which branching takes place
c nbranch    : number of branching levels
c nbr        : number of branching levels
c
      branchchar='000.branch'
      write(branchchar(1:3),'(i3.3)') Atarget
      branchfile=trim(path)//'levels/branch/'//
     +  trim(nuc(Ztarget))//branchchar
      inquire (file=branchfile,exist=lexist)
      if (lexist) then
        open (unit=2,file=branchfile,status='old')
   10   read(2,'(a)',end=100) line
        call getkeywords(line,word)
        read(word(2),*,end=1000,err=1000) iz
        read(word(3),*,end=1000,err=1000) ia
        read(word(4),*,end=1000,err=1000) ilev0
        read(word(5),*,end=1000,err=1000) nbr
        Zix=Zinit-iz
        Nix=Ninit-ia+iz
        if (Zix.lt.0.or.Zix.gt.numZ.or.Nix.lt.0.or.Nix.gt.numN)
     +    goto 10
        if (ilev0.lt.0..or.ilev0.gt.numlev) goto 1060
        if (nbr.lt.0..or.nbr.gt.numlev) goto 1060
        if (nbranch(Zix,Nix,ilev0).ne.0) goto 10
        iword=5
        sum=0.
        do k=1,nbr
          iword=iword+1
          read(word(iword),*,end=1000,err=1000) ilev1
          if (ilev1.lt.0..or.ilev1.gt.numlev) goto 1060
          iword=iword+1
          read(word(iword),*,end=1000,err=1000) bra
          if (bra.lt.0.) goto 1070
          branchlevel(Zix,Nix,ilev0,k)=ilev1
          branchratio(Zix,Nix,ilev0,k)=bra
          sum=sum+bra
        enddo
        if (sum.gt.0.) then
          do k=1,nbr
            branchratio(Zix,Nix,ilev0,k)=branchratio(Zix,Nix,ilev0,k)
     +        /sum
          enddo
        endif
        nbranch(Zix,Nix,ilev0)=nbr
        goto 10
  100   close (unit=2)
      endif
      return
 1000 write(*,'(" TALYS-error: Wrong branching ratio: ",a)') trim(line)
      stop
 1060 write(*,'(" TALYS-error: 0 <= level number <= ",i4,
     +  ", ilev0, ilev1 or nbr index out of range: ",a)')
     +  numlev,trim(line)
      stop
 1070 write(*,'(" TALYS-error: 0 <= branching ratio",
     +  ", br index out of range: ",a)') trim(line)
      stop
      end
Copyright (C)  2023 A.J. Koning, S. Hilaire and S. Goriely
