      function cmsd(Ein)
c
c +---------------------------------------------------------------------
c | Author: Arjan Koning    
c | Date  : August 11, 2004
c | Task  : Normalization factor for MSD
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      real cmsd,Ein
c
c ***************** Normalization factor for MSD ***********************
c
c cmsd      : normalization factor for MSD
c Ein       : incident energy
c M2constant: constant for matrix element in exciton model (here used
c             for MSD model)
c Atarget   : mass number of target nucleus
c
      cmsd=M2constant*3.3e6/(Atarget**3)/Ein
      return
      end
Copyright (C) 2004  A.J. Koning, S. Hilaire and M.C. Duijvestijn
