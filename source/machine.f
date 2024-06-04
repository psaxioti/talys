      subroutine machine
c
c +---------------------------------------------------------------------
c | Author: Arjan Koning 
c | Date  : October 5, 2006
c | Task  : Machine dependent statements
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      logical lexist
      integer i
c
c ********************* Set null device ********************************
c
c The null device is a "black hole" for output that is produced, but not
c of interest to the user. Some ECIS output files fall in this category.
c To ensure compatibility with Unix, Linux, Windows and other systems
c a null device string is used, of which the default setting is given  
c here. The input file may also be used to alter this setting, through 
c the nulldev keyword (suggestion of Michael Borchard).
c
c nulldev: null device
c
      nulldev='/dev/null    '
cWindows
c     nulldev='NUL          '
c
c ********************* Set directory for structure data ***************
c
c path   : directory containing structure files to be read
c lenpath: length of pathname
c
c The maximum length of the path is 60 characters
c
      path='/home/finux01b/akoning/talys/structure/'
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
      write(*,'(" TALYS-error: Structure database not installed:",
     +  " change path in machine.f")')
      stop
      end
Copyright (C) 2004  A.J. Koning, S. Hilaire and M.C. Duijvestijn
