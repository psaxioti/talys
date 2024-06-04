      subroutine dtheory(Zix,Nix,l,E)
c
c +---------------------------------------------------------------------
c | Author: Arjan Koning 
c | Date  : January 8, 2009
c | Task  : Theoretical calculation of average neutron spacings
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      integer          l,i,k,Nres,Zix,Nix,jj2b,jj2e,J2b,J2e,parity,J2,J
      real             E
      double precision rho,rhosum,density
c
c ********************** Level density parameters **********************
c
c Zix    : charge number index for residual nucleus
c Nix    : neutron number index for residual nucleus
c l      : orbital angular momentum
c E      : incident energy
c Nres   : maximal neutron number index for residual nucleus
c numN   : maximal number of neutrons away from the initial compound
c          nucleus
c J2b    : 2 * start of J summation 
c J2e    : 2 * end of J summation 
c jdis   : spin of level
c rhosum : help variable
c parity : parity
c sgn    : +1 for even argument, -1 for odd argument
c density: level density 
c S      : separation energy per particle 
c J      : J-value
c ldmodel: level density model 
c Djltheo: theoretical resonance spacing per J,l value
c Dltheo : theoretical resonance spacing per l value
c
      do 10 i=0,numl
        Dltheo(i)=0.
        do 10 k=0,numJ
          Djltheo(k,i)=0.
   10 continue
      Nres=min(numN,Nix+1)
      call levels(Zix,Nres)
      jj2b=int(abs(2.*(jdis(Zix,Nres,0)-0.5)))
      jj2e=int(2.*(jdis(Zix,Nres,0)+0.5))
      J2b=jj2b-2*l
      J2e=jj2e+2*l
      rhosum=0.
      do 20 J2=J2b,J2e,2
        if (J2.lt.0) goto 20
        parity=sgn((J2e-J2)/2)
        rho=density(Zix,Nix,max(0.,S(Zix,Nix,1)+E),0.5*J2,parity,0,
     +    ldmodel)
        J=J2/2
        Djltheo(J,l)=real(1.e6/rho)
        rhosum=rhosum+rho
   20 continue
      if (rhosum.ge.1.e-10) Dltheo(l)=real(1.e6/rhosum)
      return 
      end
Copyright (C) 2004  A.J. Koning, S. Hilaire and M.C. Duijvestijn
