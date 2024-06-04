      subroutine xsfit(Z,A)
c
c +---------------------------------------------------------------------
c | Author: Arjan Koning
c | Date  : October 25, 2023
c | Task  : Adjusted parameters to fit cross sections
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      integer numfit
      parameter (numfit=10)
      logical       lexist,first
      logical       flagassign
      character*1   yesno,ch
      character*2   ele
      character*132 word(40),key,value,cval,line
      character*132 colf,proj,elf
      character*132 ffile(numfit),parfile,keystring(20)
      integer       Z,A,fitpar(numfit),iz,ia,i,k,psf,fis,ldm,aomp,Npar,
     +              Lt,Zix,Nix,igr,type,irad,istat,class,ibar,lval,ival,
     +              keyix
      real          val
c
c ************************ Initialization ******************************
c
      fitpar = 0
      ffile = 'x'
      if (flagnnfit) fitpar(1) = 1
      if (flagngfit) fitpar(2) = 1
      if (flagnafit) fitpar(3) = 1
      if (flagnffit) fitpar(4) = 1
      if (flaganfit) fitpar(5) = 1
      if (flagdnfit) fitpar(6) = 1
      if (flagpnfit) fitpar(7) = 1
      if (flaggnfit) fitpar(8) = 1
      if (flaggamgamfit) fitpar(9) = 1
      if (flagmacsfit) fitpar(10) = 1
      ffile(1) = 'nn.par'
      ffile(2) = 'ng.par'
      ffile(3) = 'na.par'
      ffile(4) = 'nf.par'
      ffile(5) = 'an.par'
      ffile(6) = 'dn.par'
      ffile(7) = 'pn.par'
      ffile(8) = 'gn.par'
      ffile(9) = 'gamgam.par'
      ffile(10) = 'macs.par'
c
c ************************* Read parameters ****************************
c
      first = .true.
      do k = 1, numfit
        if (fitpar(k) == 0) cycle
        parfile = trim(path)//'best/fits/'//trim(ffile(k))
        inquire (file = parfile, exist = lexist)
        if ( .not. lexist) cycle
        Lt = 0
        open (unit = 2, file = parfile, status = 'old')
        Loop1: do
          do
            read(2, '(a)', iostat = istat) line
            if (istat == - 1) exit Loop1
            key = 'projectile'
            keyix=index(line,trim(key))
            if (keyix > 0) read(line(keyix+len_trim(key):132),'(a)',
     +        iostat = istat) proj
            key = 'element'
            keyix=index(line,trim(key))
            if (keyix > 0) read(line(keyix+len_trim(key):132),'(a)',
     +        iostat = istat) elf
            key = 'mass'
            keyix=index(line,trim(key))
            if (keyix > 0) read(line(keyix+len_trim(key):132),*,
     +        iostat = istat) ia
            key = 'ldmodel'
            keyix=index(line,trim(key))
            if (keyix > 0) read(line(keyix+len_trim(key):132),*,
     +        iostat = istat) ldm
            key = 'colenhance'
            keyix=index(line,trim(key))
            if (keyix > 0) read(line(keyix+len_trim(key):132),'(a)',
     +        iostat = istat) colf
            key = 'strength'
            keyix=index(line,trim(key))
            if (keyix > 0) read(line(keyix+len_trim(key):132),*,
     +        iostat = istat) psf
            key = 'fismodel'
            keyix=index(line,trim(key))
            if (keyix > 0) read(line(keyix+len_trim(key):132),*,
     +        iostat = istat) fis
            key = 'alphaomp'
            keyix=index(line,trim(key))
            if (keyix > 0) read(line(keyix+len_trim(key):132),*,
     +        iostat = istat) aomp
            key = '##parameters'
            keyix=index(line,trim(key))
            if (keyix > 0) then
              read(line(keyix+len_trim(key):132),*, iostat = istat) Npar
              do i = 1, Npar
                read(2, '(a)') keystring(i)
              enddo
            endif
            key = '#####'
            keyix=index(line,trim(key))
            if (keyix > 0) exit
          enddo
          ele = trim(adjustl(elf))
          ch = ele(1:1)
          if (ch  >= 'a'  .and. ch <= 'z') ele(1:1)=achar(iachar(ch)-32)
          if (trim(ele) /= trim(nuc(Z))) cycle
          if (ia /= A) cycle
          if (Lt /= Ltarget) cycle
          if (ldm /= ldmodel(0, 0)) cycle
          if (k/=4.and.trim(adjustl(colf))/=yesno(flagcol(0,0))) cycle
          if ((k==2.or.k==8.or.k==9.or.k==10).and.psf/=strength) cycle
          if ((k == 3 .or. k == 5) .and. aomp /= alphaomp) cycle
          if (k == 4 .and. fis /= fismodel) cycle
          if (first) then
            open (unit = 1, file = 'adjust.dat', status = 'replace')
            first = .false.
            write(1, '("# Adjusted parameters from fitting")')
          else
            open (unit=1,file='adjust.dat',status='old',
     +        position='append')
          endif
          do i = 1, Npar
            write(1, '(a)') trim(keystring(i))
          enddo
          close(1)
          do i = 1, Npar
            line = keystring(i)
            call getkeywords(line, word)
            key = word(1)
            value = word(2)
            ch = word(2)(1:1)
            Zix = 0
            Nix = 0
            type = 0
            lval = 0
            ibar = 0
            igr = 1
            if (key == 'rvadjust') then
              class = 6
              call getvalues(class, word, Zix, Nix, type, ibar, irad,
     +          lval, igr, val, ival, cval, flagassign)
              if (flagassign) then
                rvadjust(type) = val
                ompadjustp(type) = .true.
              endif
              cycle
            endif
            if (key == 'rwdadjust') then
              class = 6
              call getvalues(class, word, Zix, Nix, type, ibar, irad,
     +          lval, igr, val, ival, cval, flagassign)
              if (flagassign) then
                rwdadjust(type) = val
                ompadjustp(type) = .true.
              endif
              cycle
            endif
            if (key == 'awdadjust') then
              class = 6
              call getvalues(class, word, Zix, Nix, type, ibar, irad,
     +          lval, igr, val, ival, cval, flagassign)
              if (flagassign) then
                awdadjust(type) = val
                ompadjustp(type) = .true.
              endif
              cycle
            endif
            if (key == 'gadjust') then
              class = 1
              call getvalues(class, word, Zix, Nix, type, ibar, irad,
     +          lval, igr, val, ival, cval, flagassign)
              if (flagassign) gadjust(Zix, Nix) = val
              cycle
            endif
            if (key == 'wtable') then
              class = 5
              call getvalues(class, word, Zix, Nix, type, ibar, irad,
     +          lval, igr, val, ival, cval, flagassign)
              if (flagassign) then
                wtable(Zix, Nix, irad, lval) = val
                gamadjust(Zix, Nix) = .true.
              endif
              cycle
            endif
            if (key == 'fisbaradjust') then
              ibar = 1
              class = 3
              call getvalues(class, word, Zix, Nix, type, ibar, irad,
     +          lval, igr, val, ival, cval, flagassign)
              if (flagassign) then
                fbaradjust(Zix, Nix, ibar) = val
                fisadjust(Zix, Nix) = .true.
              endif
              cycle
            endif
            if (key == 'fishwadjust') then
              ibar = 1
              class = 3
              call getvalues(class, word, Zix, Nix, type, ibar, irad,
     +          lval, igr, val, ival, cval, flagassign)
              if (flagassign) then
                fwidthadjust(Zix, Nix, ibar) = val
                fisadjust(Zix, Nix) = .true.
              endif
              cycle
            endif
            if (key == 'vfiscor') then
              class = 1
              call getvalues(class, word, Zix, Nix, type, ibar, irad,
     +          lval, igr, val, ival, cval, flagassign)
              if (flagassign) vfiscor(Zix, Nix) = val
              cycle
            endif
            if (key == 'ctable') then
              class = 3
              call getvalues(class, word, Zix, Nix, type, ibar, irad,
     +          lval, igr, val, ival, cval, flagassign)
              if (flagassign) then
                ctable(Zix, Nix, ibar) = val
                ldadjust(Zix, Nix) = .true.
              endif
              cycle
            endif
            if (key == 'ptable') then
              class = 3
              call getvalues(class, word, Zix, Nix, type, ibar, irad,
     +          lval, igr, val, ival, cval, flagassign)
              if (flagassign) then
                ptable(Zix, Nix, ibar) = val
                ldadjust(Zix, Nix) = .true.
              endif
              cycle
            endif
            if (key == 'ctableadjust') then
              class = 3
              call getvalues(class, word, Zix, Nix, type, ibar, irad,
     +          lval, igr, val, ival, cval, flagassign)
              if (flagassign) then
                ctableadjust(Zix, Nix, ibar) = val
                ldadjust(Zix, Nix) = .true.
              endif
              cycle
            endif
            if (key == 'ptableadjust') then
              class = 3
              call getvalues(class, word, Zix, Nix, type, ibar, irad,
     +          lval, igr, val, ival, cval, flagassign)
              if (flagassign) then
                ptableadjust(Zix, Nix, ibar) = val
                ldadjust(Zix, Nix) = .true.
              endif
              cycle
            endif
            if (key == 's2adjust') then
              class = 3
              call getvalues(class, word, Zix, Nix, type, ibar, irad,
     +          lval, igr, val, ival, cval, flagassign)
              if (flagassign) s2adjust(Zix, Nix, ibar) = val
              cycle
            endif
            if (key == 'risomer') then
              class = 1
              call getvalues(class, word, Zix, Nix, type, ibar, irad,
     +          lval, igr, val, ival, cval, flagassign)
              if (flagassign) risomer(Zix, Nix) = val
              cycle
            endif
            if (key == 'cstrip') then
              class = 6
              call getvalues(class, word, Zix, Nix, type, ibar, irad,
     +          lval, igr, val, ival, cval, flagassign)
              if (flagassign) Cstrip(type) = val
              cycle
            endif
          enddo
          exit
        enddo Loop1
      enddo
      return
      end
Copyright (C)  2023 A.J. Koning, S. Hilaire and S. Goriely
