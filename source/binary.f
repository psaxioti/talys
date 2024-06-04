      subroutine binary
c
c +---------------------------------------------------------------------
c | Author: Arjan Koning and Stephane Hilaire
c | Date  : October 14, 2004
c | Task  : Binary reaction results
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      integer type,Zix,Nix,NL,nex,J,parity,nen,Z,N,A,odd
      real    factor,Eex,ald,ignatyuk,spindis,xscompall,Eaveragesum,
     +        frac,Eaverage(0:numpar)
c
c * Add direct and pre-equilibrium cross sections to population arrays *
c
c parskip     : logical to skip outgoing particle
c Zindex,Zix  : charge number index for residual nucleus
c Nindex,Nix  : neutron number index for residual nucleus
c Nlast,NL    : last discrete level
c xsdirdisc   : direct cross section for discrete state
c jdis        : spin of level
c parity      : parity
c parlev      : parity of level
c xspop       : population cross section
c xspopex     : population cross section summed over spin and parity
c xspopex0    : binary population cross section for discrete states
c xspopnuc    : population cross section per nucleus
c xsdirdisctot: direct cross section summed over discrete states
c xsbinary    : cross section from initial compound to residual nucleus 
c
c Depending on whether the nucleus is even or odd, the quantum number J
c appearing in the array xspop represents either J or J+0.5.
c
      do 10 type=0,6
        if (parskip(type)) goto 10
        Zix=Zindex(0,0,type)
        Nix=Nindex(0,0,type)
        NL=Nlast(Zix,Nix,0)
        do 20 nex=0,NL
          if (xsdirdisc(type,nex).ne.0.) then
            J=int(jdis(Zix,Nix,nex))
            parity=parlev(Zix,Nix,nex)
            xspop(Zix,Nix,nex,J,parity)=xspop(Zix,Nix,nex,J,parity)+
     +        xsdirdisc(type,nex)
            xspopex(Zix,Nix,nex)=xspopex(Zix,Nix,nex)+
     +        xsdirdisc(type,nex)
          endif
          xspopex0(type,nex)=xspopex(Zix,Nix,nex)
   20   continue
        xspopnuc(Zix,Nix)=xspopnuc(Zix,Nix)+xsdirdisctot(type)
        xsbinary(type)=xsbinary(type)+xsdirdisctot(type)
c
c
c ******* Assign spin distribution to pre-equilibrium population *******
c
c flagpreeq : flag for pre-equilibrium calculation
c flagpespin: flag for pre-equilibrium spin distribution or compound
c             spin distribution for pre-equilibrium cross section   
c maxex     : maximum excitation energy bin for compound nucleus
c maxJph    : maximal spin for particle-hole states 
c factor    : help variable
c preeqpop  : pre-equilibrium population cross section
c preeqpopex: pre-equilibrium population cross section summed over
c             spin and parity
c Eex,Ex    : excitation energy
c ald       : level density parameter
c ignatyuk  : function for energy dependent level density parameter a 
c spindis   : Wigner spin distribution
c pardis    : parity distribution               
c xspreeqtot: preequilibrium cross section per particle type
c xsgrtot   : total smoothed giant resonance cross section
c             incident channel 
c
        if (flagpreeq) then
          if (.not.flagpespin) then
            do 30 nex=NL+1,maxex(Zix,Nix)
              if (xspopex(Zix,Nix,nex).ne.0.) then
                do 40 parity=-1,1,2
                  do 40 J=0,maxJph
                    factor=xspop(Zix,Nix,nex,J,parity)/
     +                xspopex(Zix,Nix,nex)
                    preeqpop(type,nex,J,parity)=factor*
     +                preeqpopex(type,nex)
   40           continue
              else
                Eex=Ex(Zix,Nix,nex)
                ald=ignatyuk(Zix,Nix,Eex,0)
                do 50 parity=-1,1,2
                  do 50 J=0,maxJph
                    factor=spindis(Zix,Nix,Eex,ald,real(J),0)*pardis
                    preeqpop(type,nex,J,parity)=factor*
     +                preeqpopex(type,nex)
   50           continue
              endif
   30       continue
          endif
          do 60 nex=NL+1,maxex(Zix,Nix)
            xspopex(Zix,Nix,nex)=xspopex(Zix,Nix,nex)+
     +        preeqpopex(type,nex)
            do 70 parity=-1,1,2
              do 70 J=0,maxJph                                         
                xspop(Zix,Nix,nex,J,parity)=
     +            xspop(Zix,Nix,nex,J,parity)+
     +            preeqpop(type,nex,J,parity)
   70       continue
   60     continue
          xspopnuc(Zix,Nix)=xspopnuc(Zix,Nix)+xspreeqtot(type)+
     +      xsgrtot(type)
          xsbinary(type)=xsbinary(type)+xspreeqtot(type)+xsgrtot(type)
        endif
c
c ************* Other total binary cross sections **********************
c
c xscompdisctot: compound cross section summed over discrete states 
c k0           : index of incident particle
c Ltarget      : excited level of target
c xsdisc       : total cross section for discrete state 
c xscompdisc   : compound cross section for discrete state      
c xsdisctot    : total cross section summed over discrete states     
c xsdircont    : direct cross section for continuum  
c xsdirect     : total direct cross section             
c xscompcont   : compound cross section for continuum
c xseps        : limit for cross sections                 
c xsconttot    : total cross section for continuum 
c xscompound   : total compound cross section    
c xscompel     : compound elastic cross section
c xselastot    : total elastic cross section (shape + compound)
c xselasinc    : total elastic cross section (neutrons only) 
c xsnonel      : non-elastic cross section 
c xsreacinc    : reaction cross section for incident channel
c xscompall    : total compound cross section summed over particles
c xsdirdiscsum : total direct cross section
c xspreeqsum   : total preequilibrium cross section summed over
c                particles
c xsgrsum      : sum over giant resonance cross sections   
c xscompnonel  : total compound non-elastic cross section        
c
        xscompdisctot(type)=0.
        do 80 nex=0,NL
          if (type.eq.k0.and.nex.eq.Ltarget) goto 80
          xsdisc(type,nex)=xspopex0(type,nex)
          xscompdisc(type,nex)=xsdisc(type,nex)-xsdirdisc(type,nex)
          xscompdisctot(type)=xscompdisctot(type)+xscompdisc(type,nex)
   80   continue
        xsdisctot(type)=xsdirdisctot(type)+xscompdisctot(type)
        xsdircont(type)=xspreeqtot(type)+xsgrtot(type)
        xsdirect(type)=xsdirdisctot(type)+xsdircont(type)
        if (xscompcont(type).lt.xseps) xscompcont(type)=0.
        xsconttot(type)=xscompcont(type)+xsdircont(type)
        xscompound(type)=xscompdisctot(type)+xscompcont(type)           
   10 continue
      xscompel=xspopex0(k0,Ltarget)
      xselastot=xselasinc+xscompel
      xsnonel=xsreacinc-xscompel
      xscompall=max(xsreacinc-xsdirdiscsum-xspreeqsum-xsgrsum,0.)
      xscompnonel=xscompall-xscompel    
c
c ***************** Create binary feeding channels *********************
c
c flagchannels: flag for exclusive channels calculation
c feedbinary  : feeding from first compound nucleus
c
c This is necessary for exclusive cross sections
c
      if (flagchannels) then
        do 110 type=0,6
          if (parskip(type)) goto 110
          Zix=Zindex(0,0,type)
          Nix=Nindex(0,0,type)
          do 120 nex=0,maxex(Zix,Nix)
            feedbinary(type,nex)=xspopex(Zix,Nix,nex)
  120     continue
  110   continue
        feedbinary(k0,Ltarget)=0.
      endif               
c
c *************** Interpolate decay on emission spectrum ***************
c
c flagspec     : flag for output of spectra
c flagrecoil   : flag for calculation of recoils
c binaryspectra: subroutine for creation of binary spectra
c
      if (flagspec.or.flagrecoil) call binaryspectra
c
c ********************** Average emission energy ***********************
c
c binemissum : integrated binary emission spectrum
c Eaveragesum: help variable
c ebegin     : first energy point of energy grid
c nendisc    : last discrete bin
c frac       : help variable
c Etop       : top of outgoing energy bin
c eoutdis    : outgoing energy of discrete state reaction
c xsbinemis  : cross section for emission from first compound nucleus
c deltaE     : energy bin around outgoing energies
c egrid      : outgoing energy grid      
c Eaverage   : average outgoing energy
c
      if (flagrecoil.or.flagspec) then
        do 210 type=0,6
          if (parskip(type)) goto 210
          binemissum(type)=0.
          Eaveragesum=0.
          do 220 nen=ebegin(type),nendisc(type)
            binemissum(type)=binemissum(type)+xsbinemis(type,nen)*
     +        deltaE(nen)
            Eaveragesum=Eaveragesum+egrid(nen)*xsbinemis(type,nen)*
     +        deltaE(nen)
  220     continue
          Zix=Zindex(0,0,type)
          Nix=Nindex(0,0,type)
          NL=Nlast(Zix,Nix,0)
          if (eoutdis(type,NL).gt.0.) then
            frac=Etop(nendisc(type))-eoutdis(type,NL)
            binemissum(type)=binemissum(type)-
     +        xsbinemis(type,nendisc(type))*frac
            Eaveragesum=Eaveragesum-egrid(nendisc(type))*
     +        xsbinemis(type,nendisc(type))*frac
          endif
          if (binemissum(type).gt.0.) then
            Eaverage(type)=Eaveragesum/binemissum(type)
          else
            Eaverage(type)=0.
          endif
  210   continue
      endif
c
c ************ Output of population after binary emission **************
c
c flagpop    : flag for output of population
c flagfission: flag for fission
c ZZ,Z       : charge number of residual nucleus
c NN,N       : neutron number of residual nucleus
c AA,A       : mass number of residual nucleus
c parname    : name of particle
c nuc        : symbol of nucleus
c eendhigh   : last energy point for energy grid for any particle
c flagcheck  : flag for output of numerical checks
c binnorm    : normalization factor for binary spectra
c odd        : odd (1) or even (0) nucleus
c Exmax      : maximum excitation energy for residual nucleus 
c deltaEx    : excitation energy bin for population arrays
c
      if (flagpop) then
        write(*,'(/"########## BINARY CHANNELS ###########")')
        write(*,'(/"++++++++++ BINARY CROSS SECTIONS ++++++++++"/)')
        if (flagfission) 
     +    write(*,'("fission  channel",23x,":",1p,e12.5)') 
     +    xsbinary(-1)
        do 310 type=0,6
          if (parskip(type)) goto 310
          Z=ZZ(0,0,type)
          N=NN(0,0,type)
          A=AA(0,0,type)
          write(*,'(a8," channel to Z=",i3," N=",i3," (",i3,a2,$)')
     +      parname(type),Z,N,A,nuc(Z)
          write(*,'("):",1p,e12.5)') xsbinary(type)
  310   continue
        if (flagspec) then
          write(*,'(/"Binary emission spectra"/)')
          write(*,'(" Energy ",7(2x,a8,2x)/)') (parname(type),type=0,6)
          do 320 nen=ebegin(0),eendhigh
            write(*,'(f8.3,1p,7e12.5)') egrid(nen),
     +        (xsbinemis(type,nen),type=0,6)
  320     continue
        endif
        if (flagspec.and.flagcheck) then
          write(*,'(/"++++++++++ CHECK OF INTEGRATED ",$)')
          write(*,'("BINARY EMISSION SPECTRA ++++++++++"/)')
          write(*,'(12x,"Continuum cross section  Integrated",$)')
          write(*,'(" spectrum  Compound normalization",$)')
          write(*,'(" Average emission energy"/)')
          do 330 type=0,6
            if (parskip(type)) goto 330
            write(*,'(a8,1p,3(10x,e12.5),0p,10x,f8.3)') parname(type),
     +        xscompcont(type)+xspreeqtot(type)+xsgrtot(type),
     +        binemissum(type),binnorm(type),Eaverage(type)
  330     continue
        endif
        write(*,'(/"++++++++++ POPULATION AFTER BINARY EMISSION",$)')
        write(*,'(" ++++++++++")')
        do 340 type=0,6
          if (parskip(type)) goto 340
          Zix=Zindex(0,0,type)
          Nix=Nindex(0,0,type)
          NL=Nlast(Zix,Nix,0)
          Z=ZZ(0,0,type)
          N=NN(0,0,type)
          A=AA(0,0,type)
          if (xspopnuc(Zix,Nix).eq.0.) goto 340
          odd=mod(A,2)
          write(*,'(/"Population of Z=",i3," N=",i3,$)') Z,N
          write(*,'(" (",i3,a2,") after binary ",$)') A,nuc(Z)
          write(*,'(a8," emission:",1p,e12.5)') parname(type),
     +      xspopnuc(Zix,Nix)
          write(*,'("Maximum excitation energy:",f8.3,$)')
     +      Exmax(Zix,Nix)
          write(*,'(" Discrete levels:",i3,$)') NL
          if (maxex(Zix,Nix).gt.NL) then
            write(*,'(" Continuum bins:",i3,$)') maxex(Zix,Nix)-NL
            write(*,'(" Continuum bin size:",f8.3/)') deltaEx(Zix,Nix)
          else
            write(*,'(/)')
          endif
          write(*,'("bin    Ex    Popul. ",$)')
          write(*,'(5("   J=",f4.1,"-   J=",f4.1,"+")/)') 
     +      (J+0.5*odd,J+0.5*odd,J=0,4)
          do 350 nex=0,maxex(Zix,Nix)
            write(*,'(i3,f8.3,1p,11e10.3)') nex,Ex(Zix,Nix,nex),
     +        xspopex(Zix,Nix,nex),((xspop(Zix,Nix,nex,J,parity),
     +        parity=-1,1,2),J=0,4)
  350     continue
  340   continue
      endif
c
c Remove compound elastic scattering from population of target state.
c
c parZ      : charge number of particle
c parN      : neutron number of particle           
c targetspin: spin of target        
c targetP   : parity of target 
c
      xspopex(parZ(k0),parN(k0),Ltarget)=0.
      xspop(parZ(k0),parN(k0),Ltarget,int(targetspin),targetP)=0.
c
c **************************** Recoils *********************************
c
c binaryrecoil: subroutine for recoil for binary reaction
c
      if (flagrecoil) call binaryrecoil
      return
      end
Copyright (C) 2004  A.J. Koning, S. Hilaire and M.C. Duijvestijn
