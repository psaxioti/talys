      subroutine dwbaread(nen1,nen2)
c
c +---------------------------------------------------------------------
c | Author: Arjan Koning and Eric Bauge
c | Date  : August 11, 2004
c | Task  : Read ECIS results for DWBA for MSD
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      character*72     line
      integer          nen1,nen2,J,iang
      double precision xs
c
c ********************** Read DWBA cross sections **********************
c
c maxJmsd   : maximal spin for MSD calculation
c xsdwin    : DWBA cross section as a function of incident energy,
c             outgoing energy and angular momentum
c flagddx   : flag for output of double-differential cross sections
c nanglecont: number of angles for continuum
c xsdw      : DWBA angular distribution as a function of incident
c             energy, outgoing energy, angular momentum and angle
c
      do 10 J=0,maxJmsd
  20    read(3,'(a72)',end=10) line
        if(line(1:1).eq.'<') goto 20
        read(line,*) xs    
        xsdwin(nen1,nen2,J,0)=real(xs)
   10 continue
      if (flagddx) then
        do 30 iang=0,nanglecont
   40      read(7,'(a72)') line
           if(line(1:1).eq.'<') goto 40
   30   continue
        do 50 J=0,maxJmsd
          do 50 iang=0,nanglecont
   60       read(7,'(a72)',end=50) line
            if(line(1:1).eq.'<') goto 60
            read(line,'(15x,e12.5)') xs
            xsdw(nen1,nen2,J,iang,0)=real(xs)
   50   continue
      endif
      return
      end
Copyright (C) 2004  A.J. Koning, S. Hilaire and M.C. Duijvestijn
