      subroutine d0theory(Zix,Nix)
c
c +---------------------------------------------------------------------
c | Author: Arjan Koning 
c | Date  : September 13, 2005
c | Task  : Theoretical calculation of D0
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      integer          Zix,Nix,Nres,J2b,J2e,parity,J2
      double precision rhosum,density
c
c ********************** Level density parameters **********************
c
c Zix    : charge number index for residual nucleus
c Nix    : neutron number index for residual nucleus
c rhosum : help variable
c Nres   : maximal neutron number index for residual nucleus
c numN   : maximal number of neutrons away from the initial compound
c          nucleus
c J2b    : 2 * start of J summation 
c J2e    : 2 * end of J summation 
c jdis   : spin of level
c parity : parity
c sgn    : +1 for even argument, -1 for odd argument
c density: level density 
c S      : separation energy per particle 
c ldmodel: level density model 
c Dtheo  : theoretical s-wave resonance spacing 
c
      rhosum=0.
      Nres=min(numN,Nix+1)
      call levels(Zix,Nres)
      J2b=int(abs(2.*(jdis(Zix,Nres,0)-0.5)))
      J2e=int(2.*(jdis(Zix,Nres,0)+0.5))
      do 10 J2=J2b,J2e,2
        parity=sgn((J2e-J2)/2)
        rhosum=rhosum+density(Zix,Nix,max(0.,S(Zix,Nix,1)),
     +    0.5*J2,parity,0,ldmodel)
   10 continue
      if (rhosum.ne.0.) Dtheo(Zix,Nix)=real(1.e6/rhosum)
      return 
      end
Copyright (C) 2004  A.J. Koning, S. Hilaire and M.C. Duijvestijn
