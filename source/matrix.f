      subroutine matrix(A,n)
c
c +---------------------------------------------------------------------
c | Author: Arjan Koning and Marieke Duijvestijn
c | Date  : August 9, 2004
c | Task  : Matrix element for exciton model 
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      integer A,n
c
c ****************** Energy dependent matrix element *******************
c
c A         : mass number of compound nucleus
c n         : exciton number
c M2        : square of matrix element
c M2constant: overall constant for matrix element in exciton model
c M2limit   : constant for asymptotical value for matrix element
c Ecomp     : total energy of composite system
c M2shift   : constant for energy shift for matrix element         
c flag2comp : flag for two-component pre-equilibrium model
c M2pipi    : square of proton-proton matrix element
c M2nunu    : square of neutron-neutron matrix element
c M2pinu    : square of proton-neutron matrix element
c Rpinu     : ratio for two-component matrix element
c M2nupi    : square of neutron-proton matrix element
c
c We use a parameterization for the matrix element that is reliable 
c between 7 and 200 MeV.
c For the one component exciton model the constant for the matrix 
c element needs to be corrected for two-component effects. Empirically,
c we find that we need to multiply the matrix element by a ratio of
c 0.50.
c We generally use the exact numerical solution (preeqmode 2) for the 
c transition rates. If the analytical solution is used (preeqmode 1), 
c we need to correct the squared matrix element by 20%.
c
      M2=M2constant*(M2limit*6.8+4.2e5/((Ecomp/n+M2shift*10.7)**3))/
     +  (A**3)
      if (preeqmode.eq.1) M2=1.20*M2
      if (flag2comp) then
        M2pipi=M2
        M2nunu=M2
        M2pinu=Rpinu*M2
        M2nupi=Rpinu*M2
      else
        M2=M2*0.50
      endif
c
c ************* Constants for optical model transition rates ***********
c
c Wompfac: adjustable constant for OMP based transition rates
c
      if (preeqmode.eq.3) then
        Wompfac(1)=M2constant*0.55/(1.+2.*Rpinu)
        Wompfac(2)=M2constant*0.55*2.*Rpinu/(1.+2.*Rpinu)
        Wompfac(0)=0.5*(Wompfac(1)+Wompfac(2))
      endif
      return
      end
Copyright (C) 2004  A.J. Koning, S. Hilaire and M.C. Duijvestijn
