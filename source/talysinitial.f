      subroutine talysinitial
c
c +---------------------------------------------------------------------
c | Author: Arjan Koning 
c | Date  : July 28, 2005
c | Task  : Initialization of nuclear structure
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
c
c ********** Initialization of constants and nuclear structure *********
c
c particles   : subroutine to determine included light particles
c nuclides    : subroutine for properties of nuclides
c flagreaction: flag for calculation of nuclear reactions
c grid        : subroutine for energy and angle grid
c flagmain    : flag for main output
c mainout     : subroutine for main output
c timer       : subroutine for output of execution time
c
      call particles
      call nuclides
      if (flagreaction) call grid
      if (flagmain) call mainout
      if (.not.flagreaction) then
        if (flagmain) call timer
        stop
      endif
      return
      end
Copyright (C) 2004  A.J. Koning, S. Hilaire and M.C. Duijvestijn
