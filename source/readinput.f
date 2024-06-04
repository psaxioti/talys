      subroutine readinput
c
c +---------------------------------------------------------------------
c | Author: Arjan Koning 
c | Date  : June 28, 2004
c | Task  : Read input
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      integer i,k
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
        write(*,'(" Number of input lines exceeds ",i5)') numlines
        write(*,'(" numlines in talys.cmb should be increased")')
        stop
      endif
      goto 10
  100 nlines=i-1
c
c ************** Convert uppercase to lowercase characters *************
c
c For easy handling of all the input parameters, everything is converted
c to lowercase characters.
c
      do 110 i=1,nlines
        do 120 k=1,80
          if (inline(i)(k:k).ge.'A'.and.inline(i)(k:k).le.'Z') 
     +      inline(i)(k:k)=char(ichar(inline(i)(k:k))+32)
  120   continue
  110 continue
      return
      end
Copyright (C) 2004  A.J. Koning, S. Hilaire and M.C. Duijvestijn
