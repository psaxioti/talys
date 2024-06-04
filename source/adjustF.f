      subroutine adjustF(type,E,factor)
c
c +---------------------------------------------------------------------
c | Author: Arjan Koning
c | Date  : December 7, 2009
c | Task  : Local optical model geometry adjustment
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      integer type,omptype,nr
      real    E,factor(6),elow,eup,emid,D,sigma,offset
c
c ************************* Woods-Saxon function ***********************
c
c factor     : Woods-Saxon multiplication factor
c E          : incident energy
c ompadjustN : number of energy ranges for local OMP adjustment
c ompadjustE1: start energy of local OMP adjustment
c ompadjustE2: end energy of local OMP adjustment
c ompadjustD : depth of local OMP adjustment
c ompadjusts : variance of local OMP adjustment
c offset     : offset to guarantee continuity
c elow,eup,..: help variables
c
      do 10 omptype=1,6
        factor(omptype)=1.
   10 continue
      do 20 omptype=1,6
        do 30 nr=1,ompadjustN(type,omptype)
          elow=ompadjustE1(type,omptype,nr)
          eup=ompadjustE2(type,omptype,nr)
          if (E.gt.elow.and.E.lt.eup) then
            emid=0.5*(elow+eup)
            D=0.01*ompadjustD(type,omptype,nr)
            sigma=ompadjusts(type,omptype,nr)
            if (sigma.eq.0.) sigma=(eup-emid)/2.
            offset=-D*exp(-(eup-emid)**2/(2.*sigma**2))
            factor(omptype)=1.+D*exp(-(E-emid)**2/(2.*sigma**2))+offset
            goto 20
          endif
   30   continue
   20 continue
      return
      end
Copyright (C) 2009  A.J. Koning, S. Hilaire and M.C. Duijvestijn
