      subroutine D0theory(Zix,Nix)
c
c +---------------------------------------------------------------------
c | Author: Arjan Koning 
c | Date  : July 7, 2004
c | Task  : Theoretical calculation of D0
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      integer          Zix,Nix,J2b,J2e,J2
      double precision rhosum,density
c
c ********************** Level density parameters **********************
c
c Zix    : charge number index for residual nucleus
c Nix    : neutron number index for residual nucleus
c J2b    : 2 * start of J summation 
c J2e    : 2 * end of J summation 
c jdis   : spin of level]
c rhosum : help variable
c density: level density 
c S      : separation energy per particle 
c ldmodel: level density model 
c Dtheo  : theoretical s-wave resonance spacing 
c
      rhosum=0.
      J2b=int(abs(2.*(jdis(Zix,Nix,0)-0.5)))
      J2e=int(2.*(jdis(Zix,Nix,0)+0.5))
      do 10 J2=J2b,J2e,2
         rhosum=rhosum+
     +     density(Zix,Nix,max(0.,S(Zix,Nix,1)),0.5*J2,0,ldmodel)
   10 continue
      if (rhosum.ne.0.) Dtheo(Zix,Nix)=real(1.e6/rhosum)
      return 
      end
Copyright (C) 2004  A.J. Koning, S. Hilaire and M.C. Duijvestijn
