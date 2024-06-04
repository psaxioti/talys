      subroutine particles
c
c +---------------------------------------------------------------------
c | Author: Arjan Koning
c | Date  : June 30, 2004
c | Task  : Determine included light particles
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      logical default
      integer type,type2
c
c ************ Determine outgoing particles to be included *************
c
c parinclude : logical to include outgoing particle
c parskip    : logical to skip outgoing particle
c flagfission: flag for fission
c outtype    : type of outgoing particles
c default    : logical to determine default outgoing particles
c parsym     : symbol of particle
c
c The default is to include all particles from photons to alpha's as 
c competing channels. If specific outgoing particles in the input are
c given, only those will be included as competing channels. The logicals
c parinclude and parskip will be referred to throughout TALYS. They are 
c always each others' opposite.
c
c 1. Default: All competing channels included
c
      do 10 type=-1,6
        parinclude(type)=.true.
        parskip(type)=.false.
   10 continue
      if (.not.flagfission) then
        parinclude(-1)=.false.
        parskip(-1)=.true.
      endif
c
c Check if default is to be used
c
      default=.true.
      do 20 type=0,6
        if (outtype(type).ne.' ') then
          default=.false.    
          goto 40
        endif
   20 continue
      do 30 type=0,6
        outtype(type)=parsym(type)
   30 continue
c
c 2. No default, but specific competing outgoing particles are 
c    requested. Now, we first turn all competing channels off, and then 
c    determine which are to be included and which are to be skipped.
c
   40 if (.not.default) then
        do 50 type=0,6
          parinclude(type)=.false.
   50   continue
        do 60 type=0,6
          do 60 type2=0,6
            if (outtype(type).eq.parsym(type2)) parinclude(type2)=.true.
   60   continue
        do 70 type=0,6
          parskip(type)=.not.parinclude(type)
   70   continue
      endif
c
c The incident particle is always included as outgoing particle
c
c k0: index of incident particle
c
      parinclude(k0)=.true.
      parskip(k0)=.false.
      return
      end
Copyright (C) 2004  A.J. Koning, S. Hilaire and M.C. Duijvestijn
