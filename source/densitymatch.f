      subroutine densitymatch(Zix,Nix)
c
c +---------------------------------------------------------------------
c | Author: Arjan Koning and Stephane Hilaire
c | Date  : July 10, 2006
c | Task  : Level density matching solution
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      integer          Zix,Nix,ibar,i,A,nEx,j
      real             ald,Exm,ignatyuk,tmp,Um,logrholoc(-1:1),Exend,
     +                 dEx,Eex,U,Krot,Kvib,Kcoll,P,logrhomatch,rhomatch
      double precision fermi
c
c ************* Determine level density matching parameters ************
c
c For non-fissile nuclei, there are no fission barriers and only the
c loop for ibar=0 is performed, i.e. for level densities relative to the
c ground state.
c
c Zix       : charge number index for residual nucleus
c Nix       : neutron number index for residual nucleus
c nfisbar   : number of fission barrier parameters
c Ntop,NP   : highest discrete level for temperature matching   
c Nlow,NLo  : lowest discrete level for temperature matching
c EL,EP     : matching level energies
c edis      : energy of level
c efistrrot : energy of rotational transition states
c
      do 10 ibar=0,nfisbar(Zix,Nix)
c
c 1. Determine matching levels and energies. The matching level is set
c    to the last discrete level. This obviously needs to be worked on.
c
        NP=Ntop(Zix,Nix,ibar)
        NLo=Nlow(Zix,Nix,ibar)
c
c A. Parameters for the ground state level density
c
        if (ibar.eq.0) then
          EL=edis(Zix,Nix,NLo)
          EP=edis(Zix,Nix,NP)
        else
c
c B. Parameters for the fission level density
c
          EL=efistrrot(Zix,Nix,ibar,NLo)
          EP=efistrrot(Zix,Nix,ibar,NP)
        endif
        if (ldmodel.eq.2.or.ldmodel.eq.3.or.ldexist(Zix,Nix,ibar)) 
     +    goto 10
c
c 2. Solve level density matching problem for Gilbert-Cameron model
c
c AA,A      : mass number of residual nucleus
c Exend     : end of possible energy region
c nEx       : number of energy points
c logrho    : logarithm of level density
c temprho   : temperature    
c Eex       : excitation energy
c pair      : total pairing correction
c ald       : level density parameter
c ignatyuk  : function for energy dependent level density parameter a
c colenhance: subroutine for collective enhancement  
c Kvib      : vibrational enhancement factor
c Krot      : rotational enhancement factor  
c Kcoll     : total collective enhancement
c fermi     : function for Fermi gas level density formula
c
c Calculate logarithm of level density and its derivative (temperature)
c on excitation energy grid.
c
        A=AA(Zix,Nix,0)
        Exend=20.+300./A
        dEx=0.1
        nEx=int(Exend/dEx)
        do 20 i=nEx,1,-1
          logrho(i)=0.
          temprho(i)=0.
          do 30 j=-1,1
            Eex=dEx*(i+0.5*j)
            U=Eex-pair(Zix,Nix)
            if (U.gt.0.) then
              ald=ignatyuk(Zix,Nix,Eex,ibar)
              call colenhance(Zix,Nix,Eex,ald,ibar,Krot,Kvib,Kcoll)
              logrholoc(j)=real(log(Kcoll*
     +          fermi(Zix,Nix,ald,Eex,pair(Zix,Nix),ibar)))
            else
              logrholoc(j)=0.
            endif
   30     continue
          logrho(i)=logrholoc(0)
          if (logrholoc(1).ne.logrholoc(-1)) 
     +      temprho(i)=dEx/(logrholoc(1)-logrholoc(-1))
          if (temprho(i).le.0.1) temprho(i)=temprho(i+1)
   20   continue
c
c Exmatch,Exm: matching point for Ex 
c matching   : subroutine to determine matching between temperature and
c              Fermi-gas region
c
        if (Exmatch(Zix,Nix,ibar).eq.0.) then
          call matching(Zix,Nix,Exm,ibar)
          Exmatch(Zix,Nix,ibar)=Exm
        else
c
c If the Exmatch given in the input is unphysical, we choose an 
c empirical value.
c
          ald=ignatyuk(Zix,Nix,Exmatch(Zix,Nix,ibar),ibar)
          P=max(pair(Zix,Nix),0.)
          if (Exmatch(Zix,Nix,ibar).gt.2.25/ald+P) then
            Exm=Exmatch(Zix,Nix,ibar)
          else
            Exm=max(2.8+266./A+P,0.1)
          endif
        endif
c
c 3. Determine parameters for constant temperature region
c
c T       : nuclear temperature
c pol1    : subroutine for interpolation of first order
c Um,tmp  : help variables
c E0      : constant of temperature formula
c rhomatch: level density at matching point
c
        if (T(Zix,Nix,ibar).eq.0.) then
          i=int(Exm/dEx)
          if (i.gt.0) then
            call pol1(i*dEx,(i+1)*dEx,temprho(i),temprho(i+1),Exm,tmp) 
            T(Zix,Nix,ibar)=tmp
          endif
        endif
        if (T(Zix,Nix,ibar).eq.0.) then
          ald=ignatyuk(Zix,Nix,Exm,ibar)
          Um=Exm-pair(Zix,Nix)
          T(Zix,Nix,ibar)=1./(sqrt(ald/Um)-1.5/Um)
        endif
        if (E0(Zix,Nix,ibar).eq.1.e-20) then
          i=int(Exm/dEx)
          call pol1(i*dEx,(i+1)*dEx,logrho(i),logrho(i+1),Exm,
     +      logrhomatch) 
          rhomatch=exp(logrhomatch)
          tmp=T(Zix,Nix,ibar)
          E0(Zix,Nix,ibar)=Exm-tmp*log(tmp*rhomatch)
        endif
   10 continue
      return
      end
Copyright (C) 2004  A.J. Koning, S. Hilaire and M.C. Duijvestijn
