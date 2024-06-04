      subroutine strucinitial
c
c +---------------------------------------------------------------------
c | Author: Arjan Koning
c | Date  : October 14, 2004
c | Task  : Initialization of arrays for various structure parameters
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      integer Zix,Nix,type,i,k,l,irad,nen,nex
c
c *********** Initialization of nuclear structure arrays ***************
c
c Masses and separation energies
c
c Nix     : neutron number index for residual nucleus
c numN    : maximal number of neutrons away from the initial 
c           compound nucleus
c Zix     : charge number index for residual nucleus
c numZ    : maximal number of protons away from the initial 
c           compound nucleus
c nucmass : mass of nucleus
c expmexc : experimental mass excess
c thmexc  : theoretical mass excess           
c numpar  : number of particles
c specmass: specific mass for target nucleus
c redumass: reduced mass
c S       : separation energy per particle
c
      do 10 Nix=0,numN+4
        do 10 Zix=0,numZ+4
          nucmass(Zix,Nix)=0.
          expmexc(Zix,Nix)=0.    
          thmexc(Zix,Nix)=0.    
   10 continue
      do 20 Nix=0,numN
        do 20 Zix=0,numZ
          do 20 type=0,numpar
            specmass(Zix,Nix,type)=0.    
            redumass(Zix,Nix,type)=0.    
            S(Zix,Nix,type)=0.    
   20 continue
c
c Level and deformation parameters
c
c numlev2    : maximum number of levels    
c jassign    : flag for assignment of spin
c passign    : flag for assignment of parity  
c parlev     : parity of level 
c edis       : energy of level
c jdis       : spin of level
c leveltype  : type of level (rotational (R) or vibrational (V))
c indexlevel : level index   
c indexcc    : level index for coupled channel
c vibband    : band number of level  
c lband      : angular momentum 
c Kmag,Kband : magnetic quantum number
c iph,iphonon: phonon (1 or 2)
c deform     : deformation parameter       
c defpar     : deformation parameter  
c numlev     : maximum number of included discrete levels
c ENSDF      : string from original ENSDF discrete level file
c tau        : lifetime of state in seconds     
c bassign    : flag for assignment of branching ratio 
c branchratio: gamma-ray branching ratio to level
c conv       : conversion coefficient                
c colltype   : type of collectivity (D, V or R)
c deftype    : deformation length (D) or parameter (B)
c nlevmax2   : maximum number of levels  
c ndef       : number of collective levels           
c nrot       : number of deformation parameters
c numrotcc   : number of rotational deformation parameters
c rotpar     : deformation parameters for rotational nucleus
c
      do 110 i=0,numlev2
        do 110 Nix=0,numN
          do 110 Zix=0,numZ
            jassign(Zix,Nix,i)=' ' 
            passign(Zix,Nix,i)=' ' 
            parlev(Zix,Nix,i)=1 
            edis(Zix,Nix,i)=0.
            jdis(Zix,Nix,i)=0.
            leveltype(Zix,Nix,i)='D'
            indexlevel(Zix,Nix,i)=0
            indexcc(Zix,Nix,i)=0
            vibband(Zix,Nix,i)=0
            lband(Zix,Nix,i)=0
            Kband(Zix,Nix,i)=0
            iphonon(Zix,Nix,i)=0
            deform(Zix,Nix,i)=0.
            defpar(Zix,Nix,i)=0.
  110 continue
      do 120 i=0,numlev
        do 120 Nix=0,numN
          do 120 Zix=0,numZ
            ENSDF(Zix,Nix,i)='                  ' 
            tau(Zix,Nix,i)=0.
  120 continue
      do 130 k=0,numlev
        do 130 i=0,numlev
          do 130 Nix=0,numN
            do 130 Zix=0,numZ
              bassign(Zix,Nix,i,k)=' '
              branchratio(Zix,Nix,i,k)=0.
              conv(Zix,Nix,i,k)=0.
  130 continue
      do 140 Nix=0,numN
        do 140 Zix=0,numZ
          colltype(Zix,Nix)=' '
          deftype(Zix,Nix)='B'
          nlevmax2(Zix,Nix)=0
          ndef(Zix,Nix)=0
          nrot(Zix,Nix)=0.
  140 continue
      do 150 i=0,numrotcc
        do 150 Nix=0,numN
          do 150 Zix=0,numZ
            rotpar(Zix,Nix,i)=0.
  150 continue
c
c Resonance parameters
c
c dD0    : uncertainty in D0
c dS0    : uncertainty in S0
c dgamgam: uncertainty in gamgam            
c
      do 210 Nix=0,numN
        do 210 Zix=0,numZ
          dD0(Zix,Nix)=0.
          dS0(Zix,Nix)=0.
          dgamgam(Zix,Nix)=0.
  210 continue
c
c Gamma parameters
c
c gammax: number of l-values for gamma multipolarity
c ngr   : number of GR  
c kgr   : constant for gamma-ray strength function         
c
      do 310 l=1,gammax
        do 310 irad=0,1
          do 310 Nix=0,numN
            do 310 Zix=0,numZ
              ngr(Zix,Nix,irad,l)=1
              kgr(Zix,Nix,irad,l)=0.
  310 continue                         
c
c Optical model parameters
c
c ompglobal  : flag for use of global optical model
c rc0,rv0,...: optical model parameters
c numomp     : number of energies on optical model file
c eomp       : energies on optical model file
c vomp       : optical model parameters from file
c
      do 410 k=1,numpar
        do 410 Nix=0,numN
          do 410 Zix=0,numZ
            ompglobal(Zix,Nix,k)=.false.
            rc0(Zix,Nix,k)=0.
            rv0(Zix,Nix,k)=0.
            av0(Zix,Nix,k)=0.
            v1(Zix,Nix,k)=0.
            v2(Zix,Nix,k)=0.
            v3(Zix,Nix,k)=0.
            w1(Zix,Nix,k)=0.
            w2(Zix,Nix,k)=0.
            rvd0(Zix,Nix,k)=0.
            avd0(Zix,Nix,k)=0.
            d1(Zix,Nix,k)=0.
            d2(Zix,Nix,k)=0.
            d3(Zix,Nix,k)=0.
            rvso0(Zix,Nix,k)=0.
            avso0(Zix,Nix,k)=0.
            vso1(Zix,Nix,k)=0.
            vso2(Zix,Nix,k)=0.
            wso1(Zix,Nix,k)=0.
            wso2(Zix,Nix,k)=0.
  410 continue
      do 420 nen=0,numomp
        do 420 k=1,numpar
          eomp(k,nen)=0.
  420 continue
      do 430 i=1,19
        do 430 nen=0,numomp
          do 430 k=1,numpar
            vomp(k,nen,i) =0.
  430 continue
c
c Fission parameters
c
c nfisbar   : number of fission barrier parameters
c nclass2   : number of sets of class2 states             
c numbar    : number of fission barriers 
c nfistrhb  : number of head band transition states for barrier
c nfisc2hb  : number of class2 states for barrier
c minertia  : moment of inertia of fission barrier deformation
c fecont    : start of continuum energy
c minertc2  : moment of inertia for class2 states
c nfistrrot : number of rotational transition states for barrier
c nfisc2rot : number of rotational class2 states per set
c Emaxclass2: maximum energy for class2 states
c pfistrhb  : parity of head band transition states                
c pfisc2hb  : parity of class2 states     
c efistrhb  : energy of head band transition states
c jfistrhb  : spin of head band transition states
c efisc2hb  : energy of class2 states
c jfisc2hb  : spin of class2 states
c pfistrrot : parity of rotational transition states     
c efistrrot : energy of rotational transition states
c jfistrrot : spin of rotational transition states
c pfisc2rot : parity of rotational class2 states                       
c efisc2rot : energy of rotational class2 states
c jfisc2rot : spin of rotational class2 states
c
      do 510 Nix=0,numN
        do 510 Zix=0,numZ
          nfisbar(Zix,Nix)=0.
          nclass2(Zix,Nix)=0.
  510 continue
      do 520 i=1,numbar
        do 520 Nix=0,numN
          do 520 Zix=0,numZ
            nfistrhb(Zix,Nix,i)=0
            nfisc2hb(Zix,Nix,i)=0
            minertia(Zix,Nix,i)=0.
            fecont(Zix,Nix,i)=0.
            minertc2(Zix,Nix,i)=0.
            nfistrrot(Zix,Nix,i)=0.
            nfisc2rot(Zix,Nix,i)=0.
            Emaxclass2(Zix,Nix,i)=0.
  520 continue
      do 530 k=0,numlev
        do 530 i=1,numbar
          do 530 Nix=0,numN
            do 530 Zix=0,numZ
              pfistrhb(Zix,Nix,i,k)=1
              pfisc2hb(Zix,Nix,i,k)=1
              efistrhb(Zix,Nix,i,k)=0.
              jfistrhb(Zix,Nix,i,k)=0.
              efisc2hb(Zix,Nix,i,k)=0.
              jfisc2hb(Zix,Nix,i,k)=0.
  530 continue
      do 540 k=0,numrot
        do 540 i=1,numbar
          do 540 Nix=0,numN
            do 540 Zix=0,numZ
              pfistrrot(Zix,Nix,i,k)=1
              efistrrot(Zix,Nix,i,k)=0.
              jfistrrot(Zix,Nix,i,k)=0.
              pfisc2rot(Zix,Nix,i,k)=1
              efisc2rot(Zix,Nix,i,k)=0.
              jfisc2rot(Zix,Nix,i,k)=0.
  540 continue
c
c Level density parameters
c
c Nlast      : last discrete level
c scutoffdisc: spin cutoff factor for discrete level region 
c ldexist    : flag for existence of level density table  
c ldmodel    : level density model       
c edens      : energy grid for tabulated level densities
c Dtheo      : theoretical s-wave resonance spacing 
c
      do 610 i=0,numbar
        do 610 Nix=0,numN
          do 610 Zix=0,numZ
            Nlast(Zix,Nix,i)=0.
            scutoffdisc(Zix,Nix,i)=0.
            ldexist(Zix,Nix,i)=.false.
  610 continue
c
c Set energy grid for tabulated level densities 
c
      if (ldmodel.eq.3) then
        edens(0)=0.
        do 620 nex=1,20
          edens(nex)=0.25*nex
  620   continue
        do 630 nex=21,30
          edens(nex)=5.+0.5*(nex-20)
  630   continue
        do 640 nex=31,40
          edens(nex)=10.+nex-30
  640   continue
        edens(41)=22.5
        edens(42)=25.
        edens(43)=30.
        do 650 nex=44,55
          edens(nex)=30.+10.*(nex-43)
  650   continue                                
      endif
      do 660 Nix=0,numN
        do 660 Zix=0,numZ
          Dtheo(Zix,Nix)=0.
  660 continue
c
c Weak coupling parameters
c
c jcore: spin of level of core nucleus
c pcore: parity of level of core nucleus      
c
      do 710 i=0,numlev2
        do 710 Nix=0,numN
          do 710 Zix=0,numZ
            jcore(Zix,Nix,i)=0.
            pcore(Zix,Nix,i)=1
  710 continue
c
c Giant resonance sum rules
c
c betagr : deformation parameter for giant resonance
c Egrcoll: energy of giant resonance
c Ggrcoll: width of giant resonance        
c
      do 810 l=0,3
        do 810 i=1,2
          betagr(l,i)=0.
          Egrcoll(l,i)=0.
          Ggrcoll(l,i)=0.
  810 continue    
c
c Q-values and threshold energies 
c
c Qres   : Q-value for residual nucleus
c Ethresh: threshold incident energy for residual nucleus 
c
      do 910 i=0,numlev
        do 910 Nix=0,numN
          do 910 Zix=0,numZ
            Qres(Zix,Nix,i)=0.
            Ethresh(Zix,Nix,i)=0.
  910 continue    
      return
      end
Copyright (C) 2004  A.J. Koning, S. Hilaire and M.C. Duijvestijn
