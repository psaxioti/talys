      function spindis(sc,Rspin)
c
c +---------------------------------------------------------------------
c | Author: Arjan Koning
c | Date  : July 6, 2022
c | Task  : Wigner spin distribution
c +---------------------------------------------------------------------
c
c ******************* Declarations and common blocks *******************
c
      real    sc,Rspin,sigma22,spindis
c
c *********************** Wigner formula ******************************
c
c sc     : spin cutoff factor
c Rspin  : spin
c sigma22: 2 * spin cutoff factor
c spindis: Wigner spin distribution
c
      sigma22=2.*sc
      spindis=(2.*Rspin+1.)/sigma22*exp(-(Rspin+0.5)**2/sigma22)
      return
      end
Copyright (C)  2022 A.J. Koning, S. Hilaire and S. Goriely
