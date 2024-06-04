      function omega(p,h,gs,Eex,rJ)
c
c +---------------------------------------------------------------------
c | Author: Arjan Koning
c | Date  : August 11, 2004
c | Task  : Particle-hole state density
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      logical surfwell
      integer p,h,Ji
      real    omega,gs,Eex,rJ,phdens
c
c ******************** Particle-hole state density *********************
c
c omega   : particle-hole state density
c p       : particle number
c h       : hole number
c gs      : single-particle level density parameter
c Eex     : excitation energy
c rJ      : spin
c surfwell: flag for surface effects in finite well 
c RnJ     : spin distribution for particle-hole states
c phdens  : particle-hole state density
c Efermi  : depth of Fermi well
c
      Ji=int(rJ)
      surfwell=.false.
      omega=(2.*rJ+1.)*RnJ(2,Ji)*phdens(p,h,gs,Eex,Efermi,surfwell)
      return
      end
Copyright (C) 2004  A.J. Koning, S. Hilaire and M.C. Duijvestijn
