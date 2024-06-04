      subroutine densitymatch(Zix,Nix)
c
c +---------------------------------------------------------------------
c | Author: Arjan Koning 
c | Date  : October 15, 2004
c | Task  : Level density matching solution
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      integer          Zix,Nix,ibar,imin,i,A,nEx,j
      real             ald,sigsum,denom,rj,Exm,ignatyuk,tmp,Um,
     +                 logrholoc(-1:1),Exend,dEx,Eex,U,Krot,Kvib,
     +                 logrhomatch,rhomatch
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
c imin      : starting index (=0 for ground states)
c                            (=1 for transition states)
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
          imin=0
          EL=edis(Zix,Nix,NLo)
          EP=edis(Zix,Nix,NP)
        else
c
c B. Parameters for the fission level density
c
          imin=1
          EL=efistrrot(Zix,Nix,ibar,NLo)
          EP=efistrrot(Zix,Nix,ibar,NP)
        endif
c
c 2. Determine spin cutoff parameter for discrete level region
c
c sigsum,denom: help variables
c jdis        : spin of level
c jfistrrot   : spin of rotational transition states
c scutoffdisc : spin cutoff factor for discrete level region
c
        sigsum=0.
        denom=0.
        do 20 i=imin,NP
          if (ibar.eq.0) then
            rj=jdis(Zix,Nix,i)
          else
            rj=jfistrrot(Zix,Nix,ibar,i)
          endif
          sigsum=sigsum+rj*(rj+1.)*(2.*rj+1.)
          denom=denom+2.*rj+1.
   20   continue
        scutoffdisc(Zix,Nix,ibar)=sigsum/(3.*denom)
c
c 3. Solve level density matching problem
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
c fermi     : function for Fermi gas level density formula
c
c Calculate logarithm of level density and its derivative (temperature)
c on excitation energy grid.
c
        A=AA(Zix,Nix,0)
        Exend=20.+300./A
        dEx=0.1
        nEx=int(Exend/dEx)
        do 30 i=nEx,1,-1
          logrho(i)=0.
          temprho(i)=0.
          do 40 j=-1,1
            Eex=dEx*(i+0.5*j)
            U=Eex-pair(Zix,Nix)
            if (U.gt.0.) then
              ald=ignatyuk(Zix,Nix,Eex,ibar)
              call colenhance(Zix,Nix,Eex,ald,ibar,Krot,Kvib)       
              logrholoc(j)=real(log(Krot*Kvib*
     +          fermi(Zix,Nix,ald,Eex,pair(Zix,Nix),ibar)))
            else
              logrholoc(j)=0.
            endif
   40     continue
          logrho(i)=logrholoc(0)
          if (logrholoc(1).ne.logrholoc(-1)) 
     +      temprho(i)=dEx/(logrholoc(1)-logrholoc(-1))
          if (temprho(i).le.0.1) temprho(i)=temprho(i+1)
   30   continue
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
          if (Exmatch(Zix,Nix,ibar).gt.2.25/ald+pair(Zix,Nix)) then
            Exm=Exmatch(Zix,Nix,ibar)
          else
            Exm=2.8+266./A+pair(Zix,Nix)
          endif
        endif
c
c 4. Determine parameters for constant temperature region
c
c T       : nuclear temperature
c pol1    : subroutine for interpolation of first order
c Um,tmp  : help variables
c E0      : constant of temperature formula
c rhomatch: level density at matching point
c
        ald=ignatyuk(Zix,Nix,Exm,ibar)
        if (T(Zix,Nix,ibar).eq.0.) then
          i=int(Exm/dEx)
          if (i.gt.0) then
            call pol1(i*dEx,(i+1)*dEx,temprho(i),temprho(i+1),Exm,tmp) 
            T(Zix,Nix,ibar)=tmp
          else
            Um=Exm-pair(Zix,Nix)
            T(Zix,Nix,ibar)=1./(sqrt(ald/Um)-1.5/Um)
          endif
        endif
        if (E0(Zix,Nix,ibar).eq.0.) then
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
