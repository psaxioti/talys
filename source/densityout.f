      subroutine densityout(Zix,Nix)
c
c +---------------------------------------------------------------------
c | Author: Arjan Koning 
c | Date  : October 16, 2004
c | Task  : Output of level density   
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      integer          Zix,Nix,Z,N,A,ibar,odd,J,nex,nen1,nen2,nen,i0,i
      real             aldmatch,SS,Eex,ignatyuk,spincut,ald,Krot,Kvib,C,
     +                 Ncum
      double precision dens,gilcam,density
c
c ********************** Level density parameters **********************
c
c Zix        : charge number index for residual nucleus
c Nix        : neutron number index for residual nucleus
c ZZ,Z       : charge number of residual nucleus
c NN,N       : neutron number of residual nucleus
c AA,A       : mass number of residual nucleus
c S,SS       : separation energy per particle 
c nuc        : symbol of nucleus
c flagfission: flag for fission
c ibar       : fission barrier
c nfisbar    : number of fission barrier parameters
c alev       : level density parameter
c aldmatch   : function to determine effective level density parameter
c
      Z=ZZ(Zix,Nix,0)
      N=NN(Zix,Nix,0)
      A=AA(Zix,Nix,0)
      SS=S(Zix,Nix,1)
      write(*,'(/"Level density parameters for Z=",i3," N=",i3,$)') Z,N
      write(*,'(" (",i3,a2,") "/)') A,nuc(Z)
      if (flagfission) then
        write(*,'(20x," g.s.     Fission barriers ")')
        write(*,'(28x,3(5x,i1,4x))') (ibar,ibar=1,nfisbar(Zix,Nix))
        write(*,'()')
      endif
      write(*,'("a(Sn)           :",f10.5)') alev(Zix,Nix)
      if (flagfission) write(*,'("a-effective     :",10x,3f10.5)') 
     +  (aldmatch(Zix,Nix,SS,ibar),ibar=1,nfisbar(Zix,Nix))
c
c Theoretical D0
c
c D0      : experimental s-wave resonance spacing in eV  
c dD0     : uncertainty in D0         
c D0theory: subroutine for theoretical calculation of D0   
c
      write(*,'("Experimental D0 :",f10.2," eV +- ",f12.5)')
     +  D0(Zix,Nix),dD0(Zix,Nix)
      write(*,'("Theoretical D0  :",f10.2," eV")') Dtheo(Zix,Nix)
c
c Other parameters
c
c alimit     : asymptotic level density parameter
c gammald    : gamma-constant for asymptotic level density parameter
c pair       : total pairing correction
c deltaW     : shell correction in nuclear mass 
c Exmatch    : matching point for Ex
c Nlast      : last discrete level
c Nlow       : lowest discrete level for temperature matching
c Ntop       : highest discrete level for temperature matching
c T          : nuclear temperature
c E0         : constant of temperature formula
c scutoffdisc: spin cutoff factor for discrete level region
c spincut    : spin cutoff factor 
c ignatyuk   : function for energy dependent level density parameter a
c Ufermi     : energy of Fermi distribution for damping of ground-state
c            : rotational effects
c cfermi     : width of Fermi distribution for damping of ground-state
c            : rotational effects
c Ufermibf   : energy of Fermi distribution for damping of barrier
c            : rotational effects
c cfermibf   : width of Fermi distribution for damping of barrier
c            : rotational effects     
c
      write(*,'("Asymptotic a    :",f10.5)') alimit(Zix,Nix)
      write(*,'("Damping gamma   :",f10.5)') gammald(Zix,Nix)
      write(*,'("Pairing energy  :",f10.5)') pair(Zix,Nix)
      write(*,'("Shell correction:",4f10.5)') (deltaW(Zix,Nix,ibar),
     +  ibar=0,nfisbar(Zix,Nix))
      write(*,'("Matching Ex     :",4f10.5)') (Exmatch(Zix,Nix,ibar),
     +  ibar=0,nfisbar(Zix,Nix))
      write(*,'("Last disc. level:",4(7x,i3))') (Nlast(Zix,Nix,ibar),
     +  ibar=0,nfisbar(Zix,Nix))
      write(*,'("Nlow            :",4(7x,i3))') (Nlow(Zix,Nix,ibar),
     +  ibar=0,nfisbar(Zix,Nix))
      write(*,'("Ntop            :",4(7x,i3))') (Ntop(Zix,Nix,ibar),
     +  ibar=0,nfisbar(Zix,Nix))
      write(*,'("Temperature     :",4f10.5)') (T(Zix,Nix,ibar),
     +  ibar=0,nfisbar(Zix,Nix))
      write(*,'("E0              :",4f10.5)') (E0(Zix,Nix,ibar),
     +  ibar=0,nfisbar(Zix,Nix))
      write(*,'("disc. spin cut  :",4f10.5)') 
     +  (sqrt(scutoffdisc(Zix,Nix,ibar)),ibar=0,nfisbar(Zix,Nix))
      write(*,'("spin cut (Sn)   :",4f10.5)') 
     +  (sqrt(spincut(Zix,Nix,ignatyuk(Zix,Nix,SS,ibar),SS,ibar)),
     +  ibar=0,nfisbar(Zix,Nix))
      if ((ldmodel.eq.1.and.flagfission).or.ldmodel.eq.2) then
        write(*,'("Krotconstant    :",4f10.5)') 
     +    (Krotconstant(Zix,Nix,ibar),ibar=0,nfisbar(Zix,Nix))
        write(*,'("Ufermi          :",f10.5)') Ufermi
        write(*,'("cfermi          :",f10.5)') cfermi
        if (ibar.gt.0) then
          write(*,'("Ufermibf        :",f10.5)') Ufermibf
          write(*,'("cfermibf        :",f10.5)') cfermibf
        endif
      endif
c
c ********************** Total level density ***************************
c
c odd: odd (1) or even (0) nucleus
c
      odd=mod(A,2)
      do 110 ibar=0,nfisbar(Zix,Nix)
        write(*,'(/"Total level density for ",$)')
        if (ibar.eq.0) then
          write(*,'(" ground state"/)')
        else
          write(*,'(" fission barrier",i3/)') ibar
        endif
        write(*,'("   Ex     a    sp cut  total ",$)')
        write(*,'(9("  JP=",f4.1,"+"),$)') (real(J+0.5*odd),J=0,8)
        if ((ldmodel.eq.1.and.ibar.gt.0).or.ldmodel.eq.2) then
          write(*,'("     Krot      Kvib"/)')
        else
          write(*,'(/)')
        endif
c
c Tabulated level densities
c
c ldmodel   : level density model 
c ldexist   : flag for existence of level density table
c edens     : energy grid for tabulated level densities
c ldtottable: total level density from table  
c ldtable   : level density from table
c pardis    : parity distribution         
c c1table   : constant to adjust tabulated level densities
c c2table   : constant to adjust tabulated level densities
c
        if (ldmodel.eq.3.and.ldexist(Zix,Nix,ibar)) then
          do 120 nex=1,55    
            Eex=edens(nex)
            write(*,'(f6.2,14x,1p,11e10.3)') Eex,
     +        ldtottable(Zix,Nix,nex,ibar),
     +        (ldtable(Zix,Nix,nex,J,ibar)*pardis,J=0,9)
  120     continue
          if (c1table(Zix,Nix).ne.1..or.c2table(Zix,Nix).ne.0.) then
            write(*,'(/"Normalization: f(Ex)= c1*exp(c2*sqrt(Ex))")')
            write(*,'("           c1=",f10.5)') c1table(Zix,Nix)
            write(*,'("           c2=",f10.5)') c2table(Zix,Nix)
          endif
        else
c
c Analytical level densities
c
c enincmax   : maximum incident energy
c k0         : index of incident particle 
c Eex        : excitation energy
c ald        : level density parameter
c dens,gilcam: Gilbert-Cameron level density
c colenhance : subroutine for collective enhancement  
c Krot       : rotational enhancement factor
c Kvib       : vibrational enhancement factor
c density    : level density 
c
          nen1=5*int(Exmatch(Zix,Nix,ibar)+1.)
          nen2=int(enincmax+S(0,0,k0)-Exmatch(Zix,Nix,ibar)-1.)
          do 130 nen=1,nen1+nen2   
            if (nen.le.nen1) then
              Eex=0.2*real(nen)
            else
              Eex=real(0.2*nen1+(nen-nen1))
            endif
            ald=ignatyuk(Zix,Nix,Eex,ibar)
            call colenhance(Zix,Nix,Eex,ald,ibar,Krot,Kvib)
            dens=Krot*Kvib*gilcam(Zix,Nix,Eex,ald,ibar)
            write(*,'(f6.2,2f7.3,1p,10e10.3,$)') Eex,ald,
     +        sqrt(spincut(Zix,Nix,ald,Eex,ibar)),dens,
     +        (density(Zix,Nix,Eex,real(J+0.5*odd),ibar,ldmodel),J=0,8)
            if ((ldmodel.eq.1.and.ibar.gt.0).or.ldmodel.eq.2) then
              write(*,'(2f9.3)') Krot,Kvib
            else
              write(*,'()')
            endif
  130     continue
        endif
  110 continue
c
c Cumulative number of discrete levels vs. constant temperature formula.
c
c Ncum    : number of cumulative levels (integral of constant 
c           temperature level density)
c C       : integration constant
c nlevmax2: maximum number of levels
c
      if (ldmodel.le.2) then
        write(*,'(/"Discrete levels versus constant temperature",$)')
        write(*,'(" formula"/)')
        write(*,'(" Energy  Level   N_cumulative"/)')
        i0=Nlow(Zix,Nix,0)
        C=i0-exp((edis(Zix,Nix,i0)- E0(Zix,Nix,0))/T(Zix,Nix,0))
        do 210 i=1,nlevmax2(Zix,Nix)
          Eex=edis(Zix,Nix,i)
          Ncum=exp((Eex-E0(Zix,Nix,0))/T(Zix,Nix,0))+C
          if (Ncum.lt.1.e7) write(*,'(f8.4,i4,f12.3)') Eex,i,Ncum
  210   continue    
      endif
      return 
      end
Copyright (C) 2004  A.J. Koning, S. Hilaire and M.C. Duijvestijn
