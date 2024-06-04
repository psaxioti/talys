      function match(Eex)
c
c +---------------------------------------------------------------------
c | Author: Arjan Koning
c | Date  : August 17, 2004
c | Task  : Matching function
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      integer          i
      real             match,Eex,dEx,temp,logrhof
      double precision rhof,factor1,factor2
c
c *********************** Matching function ****************************
c
c match       : matching function
c Eex         : excitation energy
c pol1        : subroutine for interpolation of first order
c temp,temprho: nuclear temperature
c rhof        : value for total level density 
c NLo,EL,NP,EP: matching level numbers and energies
c factor1,2   : help variables
c
      dEx=0.1
      i=max(int(Eex/dEx),1)
      call pol1(i*dEx,(i+1)*dEx,temprho(i),temprho(i+1),Eex,temp)
      if (temp.gt.0.) then
        call pol1(i*dEx,(i+1)*dEx,logrho(i),logrho(i+1),Eex,logrhof)
        rhof=exp(logrhof)
        factor1=exp(-Eex/temp)
        factor2=exp(EP/temp)
        if (EL.ne.0.) factor2=factor2-exp(EL/temp)
        match=real(temp*rhof*factor1*factor2)+NLo-NP
      endif
      return
      end
Copyright (C) 2004  A.J. Koning, S. Hilaire and M.C. Duijvestijn
