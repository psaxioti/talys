      subroutine machine
c
c +---------------------------------------------------------------------
c | Author: Arjan Koning 
c | Date  : December 5, 2004
c | Task  : Machine dependent statements
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      logical lexist
      integer i,istat
c
c ********************* Set directory for structure data ***************
c
c path   : directory containing structure files to be read
c lenpath: length of pathname
c
      path='/home/l2634/akoning/import/talys/structure/'
      lenpath=0
      do 10 i=1,60
        if (path(i:i).eq.' ') goto 110
        lenpath=lenpath+1
   10 continue                       
c
c Test to check accessibility of structure files
c
  110 inquire (file=path(1:lenpath)//'abundance/z001',exist=lexist)
      if (lexist) return
      write(*,'("TALYS-error: Structure database not installed:",$)')
      write(*,'(" change path in machine.f")')
      stop
      end
Copyright (C) 2004  A.J. Koning, S. Hilaire and M.C. Duijvestijn
