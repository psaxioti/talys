      double precision function mliquid(Z,A)
c
c +---------------------------------------------------------------------
c | Author: Arjan Koning 
c | Date  : July 2, 2004
c | Task  : Myers-Swiatecki liquid drop mass
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      integer          Z,A,N,oddZ,oddN
      double precision a1,a2,kappa,c3,c4,Mn,Mh,rA,rZ,rN,factor,Ev,Esur,
     +                 Ecoul,delta,Eldm
c
c ************** Spherical Myers-Swiatecki parameters ******************
c
c mliquid  : function for liquid drop mass
c AA,rA,A  : mass number of residual nucleus
c ZZ,rZ,Z  : charge number of residual nucleus
c a1,a2,.. : Myers-Swiatecki parameters
c Mn,Mh    : mass excess of neutron and hydrogen
c NN,rN,N  : neutron number of residual nucleus
c factor   : help variable
c Ev       : volume energy
c Esur     : surface energy
c twothird : 2/3
c Ecoul    : Coulomb energy
c onethird : 1/3
c oddZ,oddN: help variables
c delta    : pairing energy
c Eldm     : liquid drop energy
c amu      : atomic mass unit in MeV
c
c Myers-Swiatecki model: Nucl. Phys. 81 (1966) 1. 
c We use the original M-S parameters, see Mengoni and Nakajima, 
c J. Nucl. Sci. Technol.  31[2], p 151 (1994).
c
      data a1,a2,kappa,c3,c4 /15.677,18.56,1.79,0.717,1.21129/
      data Mn,Mh/ 8.07144,7.28899/
      N=A-Z
      rA=real(A)
      rZ=real(Z)
      rN=real(N)
      factor=1.-kappa*((rN-rZ)/rA)**2
      Ev=-a1*factor*rA
      Esur=a2*factor*rA**twothird
      Ecoul=c3*rZ**2/(rA**onethird)-c4*rZ**2/rA
      oddZ=mod(Z,2)
      oddN=mod(N,2)
      if (oddZ.ne.oddN) delta=0.
      if (oddZ.eq.0.and.oddN.eq.0) delta=-11./sqrt(rA)
      if (oddZ.eq.1.and.oddN.eq.1) delta=11./sqrt(rA)
      Eldm=Z*Mh+N*mn+Ev+Esur+Ecoul+delta
      mliquid=A+Eldm/amu
      return
      end
Copyright (C) 2004  A.J. Koning, S. Hilaire and M.C. Duijvestijn
