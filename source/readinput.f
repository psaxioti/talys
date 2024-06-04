      subroutine readinput
c
c +---------------------------------------------------------------------
c | Author: Arjan Koning
c | Date  : July 2, 2011
c | Task  : Read input
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      character*80 str
      integer      i,k
c
c ************************** User Input ********************************
c
c iso     : counter for isotope
c inline  : input line
c numlines: maximum number of input lines
c nlines  : number of input lines
c
c We read the complete input file first as a set of character strings.
c The actual keywords will be read from these later on. For natural
c elements, the input file only needs to be read once.
c
      if (iso.ne.1) return
      i=1
   10 read(*,'(a80)',end=100) inline(i)
      i=i+1
      if (i.gt.numlines) then
        write(*,'(" TALYS-error: Number of input lines exceeds ",i5)')
     +    numlines
        write(*,'(" numlines in talys.cmb should be increased")')
        stop
      endif
      goto 10
  100 nlines=i-1
c
c ************** Convert uppercase to lowercase characters *************
c
c For easy handling of all the input parameters, the whole input is
c converted to lowercase characters, with the exception of filenames or
c other character strings.
c
      do 110 i=1,nlines
        str(1:80)=inline(i)(1:80)
        do 120 k=1,80
          if (inline(i)(k:k).ge.'A'.and.inline(i)(k:k).le.'Z')
     +      inline(i)(k:k)=char(ichar(inline(i)(k:k))+32)
  120   continue
        do 130 k=0,60
          if (inline(i)(k+1:k+7).eq.'energy ') then
            inline(i)(k+8:80)=str(k+8:80)
            goto 110
          endif
          if (inline(i)(k+1:k+7).eq.'optmod ') then
            inline(i)(k+8:80)=str(k+8:80)
            goto 110
          endif
          if (inline(i)(k+1:k+7).eq.'nulldev') then
            inline(i)(k+8:80)=str(k+8:80)
            goto 110
          endif
          if (inline(i)(k+1:k+8).eq.'bestpath') then
            inline(i)(k+9:80)=str(k+9:80)
            goto 110
          endif
          if (inline(i)(k+1:k+8).eq.'integral') then
            inline(i)(k+9:80)=str(k+9:80)
            goto 110
          endif
          if (inline(i)(k+1:k+9).eq.'abundance') then
            inline(i)(k+10:80)=str(k+10:80)
            goto 110
          endif
          if (inline(i)(k+1:k+9).eq.'levelfile') then
            inline(i)(k+10:80)=str(k+10:80)
            goto 110
          endif
          if (inline(i)(k+1:k+9).eq.'strucpath') then
            inline(i)(k+10:80)=str(k+10:80)
            goto 110
          endif
          if (inline(i)(k+1:k+10).eq.'class2file') then
            inline(i)(k+11:80)=str(k+11:80)
            goto 110
          endif
          if (inline(i)(k+1:k+10).eq.'deformfile') then
            inline(i)(k+11:80)=str(k+11:80)
            goto 110
          endif
          if (inline(i)(k+1:k+10).eq.'radialfile') then
            inline(i)(k+11:80)=str(k+11:80)
            goto 110
          endif
          if (inline(i)(k+1:k+10).eq.'rescuefile') then
            inline(i)(k+11:80)=str(k+11:80)
            goto 110
          endif
          if (inline(i)(k+1:k+11).eq.'hbtransfile') then
            inline(i)(k+12:80)=str(k+12:80)
            goto 110
          endif
          if (inline(i)(k+1:k+11).eq.'optmodfilen') then
            inline(i)(k+12:80)=str(k+12:80)
            goto 110
          endif
          if (inline(i)(k+1:k+11).eq.'optmodfilep') then
            inline(i)(k+12:80)=str(k+12:80)
            goto 110
          endif
          if (inline(i)(k+1:k+13).eq.'ompenergyfile') then
            inline(i)(k+14:80)=str(k+14:80)
            goto 110
          endif
  130   continue
  110 continue
      return
      end
Copyright (C) 2004  A.J. Koning, S. Hilaire and M.C. Duijvestijn
