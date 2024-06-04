      subroutine strucinitial
c
c +---------------------------------------------------------------------
c | Author: Arjan Koning
c | Date  : January 31, 2023
c | Task  : Initialization of arrays for various structure parameters
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      logical       lexist
      character*2   phstring1(14)
      character*4   phstring2(72)
      character*132 denfile
      integer       Zix,Nix,Z,N,A,l,i,k,nex,nen,ipar,j
      real          Eout,degrid,Eeps
c
c *********** Initialization of nuclear structure arrays ***************
c
c Masses and separation energies
c
c Nix         : neutron number index for residual nucleus
c numN        : maximal number of neutrons away from the initial
c               compound nucleus
c Zix         : charge number index for residual nucleus
c numZ        : maximal number of protons away from the initial
c               compound nucleus
c nucmass     : mass of nucleus
c expmass     : experimental mass
c thmass      : theoretical mass
c expmexc     : experimental mass excess
c thmexc      : theoretical mass excess
c dumexc      : theoretical mass excess from Duflo-Zuker formula
c beta4       : deformation parameters
c gsspin      : ground state spin
c gsparity    : ground state parity
c ldparexist  : flag for existence of tabulated level density parameters
c numpar      : number of particles
c specmass    : specific mass for target nucleus
c redumass    : reduced mass
c S           : separation energy per particle
c numen       : maximum number of outgoing energies
c egrid       : energies of basic energy grid in MeV
c Q           : Q value
c coullimit   : energy limit for charged particle OMP calculation
c numelem     : number of elements
c nummass     : number of masses
c lmaxinc     : maximal l-value for transmission coefficients for
c               incident channel
c flagreaction: flag for calculation of nuclear reactions
c
      nucmass=0.
      expmass=0.
      thmass=0.
      expmexc=0.
      thmexc=0.
      dumexc=0.
      beta4=0.
      gsparity=1
      gsspin=0.
      do Nix=0,numN+4
        do Zix=0,numZ+4
          Z=Zinit-Zix
          N=Ninit-Nix
          A=Z+N
          if (mod(A,2).eq.0) then
            gsspin(Zix,Nix)=0.
          else
            gsspin(Zix,Nix)=0.5
          endif
        enddo
      enddo
      ldparexist=.false.
      specmass=0.
      redumass=0.
      S=0.
      egrid=0.
      Q=0.
      coullimit=0.
      lmaxinc=0
c
c Level and deformation parameters
c
c flaginvecis: logical for calculating inverse channel OMP
c Ltarget    : excited level of target
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
c tau        : lifetime of state in seconds
c bassign    : flag for assignment of branching ratio
c conv       : conversion coefficient
c colltype   : type of collectivity (D, V or R)
c deftype    : deformation length (D) or parameter (B)
c nlevmax2   : maximum number of levels
c ndef       : number of collective levels
c nrot       : number of deformation parameters
c Lisomer    : level number of isomer
c Nisomer    : number of isomers for this nuclide
c deltaEx    : excitation energy bin for population arrays
c sfactor    : spin factor
c numrotcc   : number of rotational deformation parameters
c rotpar     : deformation parameters for rotational nucleus
c
      flaginvecis=.true.
      Ltarget0=Ltarget
      jassign=' '
      passign=' '
      parlev=1
      levnum=0
      edis=0.
      jdis=0.
      tau=0.
      leveltype='D'
      indexlevel=0
      indexcc=0
      vibband=0
      lband=0
      Kband=0
      iphonon=0
      deform=0.
      defpar=0.
      bassign=' '
      conv=0.
      colltype=' '
      deftype='B'
      nlevmax2=0
      ndef=0
      nrot=0
      deltaEx=0.
      sfactor=0.
      Nisomer=0
      Lisomer=0
      rotpar=0.
c
c Resonance parameters
c
c swaveth : theoretical strength function for s-wave
c dD0     : uncertainty in D0
c dgamgam : uncertainty in gamgam
c gamgamth: theoretical total radiative width
c
      swaveth=0.
      dD0=0.
      dgamgam=0.
      gamgamth=0.
c
c Decay data parameters
c
c lambda: decay rate per isotope
c prate : production rate per isotope
c rtyp  : type of beta decay, beta-: 1 , beta+: 2 (from ENDF format)
c Thalf : half life of nuclide in sec.
c Td    : half life per time unit
c
      lambda=0.
      prate=0.
      rtyp=0
      Thalf=1.e30
      Td=0
c
c Gamma parameters
c
c numgam    : maximum number of l-values for gamma multipolarity
c ngr       : number of GR
c kgr       : constant for gamma-ray strength function
c qrpaexist : flag for existence of tabulated QRPA strength functions
c numgamqrpa: number of energies for QRPA strength function
c eqrpa     : energy grid for QRPA strength function
c fqrpa     : tabulated QRPA strength function
c xsracap   : direct-semidirect radiative capture cross section
c xsracapEM : direct-semidirect radiative capture cross section as
c             function of type
c
      do l=1,numgam
        kgr(l)=pi2h2c2/(2*l+1.)
      enddo
      ngr=1
      qrpaexist=.false.
      eqrpa=0.
      fqrpa=0.
      xsracap=0.
      xsracapEM=0.
c
c If no nuclear reaction calculation is requested, we skip a large
c part of this subroutine to speed up the calculation.
c
      if (.not.flagreaction) goto 600
c
c Optical model parameters
c
c ompglobal  : flag for use of global optical model
c rc0,rv0,...: optical model parameters
c disp       : flag for dispersive optical model
c jlmexist   : flag for existence of tabulated radial matter density
c normjlm    : JLM potential normalization factors
c numomp     : number of energies on optical model file
c eomp       : energies on optical model file
c vomp       : optical model parameters from file
c numNph     : maximal number of neutrons away from the initial
c              compound nucleus for multiple pre-equilibrium emission
c numZph     : maximal number of protons away from the initial
c              compound nucleus for multiple pre-equilibrium emission
c wvol       : absorption part of the optical potential averaged over
c              the volume
c flagjlm    : flag for using semi-microscopic JLM OMP
c flagracap  : flag for radiative capture model
c alphaomp   : alpha optical model (1=normal, 2= McFadden-Satchler,
c              3-5= folding potential, 6,8= Avrigeanu, 7=Nolte)
c rhojlmn    : density for neutrons
c rhojlmp    : density for protons
c potjlm     : JLM potential depth values
c radjlm     : radial points for JLM potential
c xstotadjust: total cross section adjustment
c xseladjust : elastic cross section adjustment
c xsnonadjust: nonelastic cross section adjustment
c threshnorm : normalization factor at trheshold
c Rprime     : potential scattering radius
c Emaxtalys  : maximum acceptable energy for TALYS
c nubarexist : flag for existence of nubar file
c
      ompglobal=.false.
      ef=0.
      rc0=0.
      rv0=0.
      av0=0.
      v1=0.
      v2=0.
      v3=0.
      w1=0.
      w2=0.
      w3=0.
      w4=0.
      rvd0=0.
      avd0=0.
      d1=0.
      d2=0.
      d3=0.
      rvso0=0.
      avso0=0.
      vso1=0.
      vso2=0.
      wso1=0.
      wso2=0.
      disp=.false.
      jlmexist=.false.
      normjlm=1.
      omplines=0
      eomp=0.
      vomp=0.
      V0=0.
      Vjoin=0.
      Wjoin=0.
      wvol=0.
      if (flagjlm.or.flagracap.or.(alphaomp.ge.3.and.alphaomp.le.5))
     +  then
        rhojlmp=0.
        rhojlmn=0.
        potjlm=0.
        radjlm=0.
      endif
      xstotadjust=0.
      xseladjust=0.
      xsnonadjust=0.
      threshnorm=1.
      Rprime=0.
      Eompbeg0=0.
      Eompbeg1=0.
      Eompend1=Emaxtalys
      Eompend0=Emaxtalys
c
c Fission parameters
c
c flagfission: flag for fission
c nfisbar    : number of fission barrier parameters
c nclass2    : number of sets of class2 states
c numbar     : number of fission barriers
c nfistrhb   : number of head band transition states for barrier
c nfisc2hb   : number of class2 states for barrier
c minertia   : moment of inertia of fission barrier deformation
c fecont     : start of continuum energy
c minertc2   : moment of inertia for class2 states
c nfistrrot  : number of rotational transition states for barrier
c nfisc2rot  : number of rotational class2 states per set
c Emaxclass2 : maximum energy for class2 states
c pfistrhb   : parity of head band transition states
c pfisc2hb   : parity of class2 states
c efistrhb   : energy of head band transition states
c jfistrhb   : spin of head band transition states
c efisc2hb   : energy of class2 states
c jfisc2hb   : spin of class2 states
c pfistrrot  : parity of rotational transition states
c efistrrot  : energy of rotational transition states
c jfistrrot  : spin of rotational transition states
c pfisc2rot  : parity of rotational class2 states
c efisc2rot  : energy of rotational class2 states
c jfisc2rot  : spin of rotational class2 states
c
      nfisbar=0
      nclass2=0
      if (flagfission) then
        nfistrhb=0
        nfisc2hb=0
        minertia=0.
        fecont=0.
        minertc2=0.
        nfistrrot=0
        nfisc2rot=0
        Emaxclass2=0.
        pfistrhb=1
        pfisc2hb=1
        efistrhb=0.
        jfistrhb=0.
        efisc2hb=0.
        jfisc2hb=0.
        pfistrrot=1
        efistrrot=0.
        jfistrrot=0.
        pfisc2rot=1
        efisc2rot=0.
        jfisc2rot=0.
        eintfis=0.
        rhofis=0.
      endif
      betafis=0.
      vfis=0.
      Vpos=0.
      Vheight=0.
      Vwidth=0.
c
c Level density parameters
c
c Nlast      : last discrete level
c Ediscrete  : energy of middle of discrete level region
c scutoffdisc: spin cutoff factor for discrete level region
c delta      : energy shift
c ldexist    : flag for existence of level density table
c edens      : energy grid for tabulated level densities
c ldmodel    : level density model
c nendens    : number of energies for level density grid
c nenphdens  : number of energies for particle-hole state density grid
c Edensmax   : maximum energy on level density table
c Ephdensmax : maximum energy on particle-hole state density table
c ENSDF      : string from original ENSDF discrete level file
c D0theo     : theoretical s-wave resonance spacing
c Econd      : condensation energy
c Ucrit      : critical U
c Scrit      : critical entropy
c Dcrit      : critical determinant
c aldcrit    : critical level density parameter
c ldtable    : level density from table
c ldtottableP: total level density per parity from table
c ldtottable : total level density from table
c rhogrid    : integrated level density
c
c
c Set energy grid for tabulated level densities
c
  600 edens(0)=0.
      do nex=1,20
        edens(nex)=0.25*nex
      enddo
      do nex=21,30
        edens(nex)=5.+0.5*(nex-20)
      enddo
      do nex=31,40
        edens(nex)=10.+nex-30
      enddo
      edens(41)=22.5
      edens(42)=25.
      edens(43)=30.
      do nex=44,60
        edens(nex)=30.+10.*(nex-43)
      enddo
      nenphdens=60
      Ephdensmax=200.
      ENSDF='                  '
      do Nix=0,numN
        do Zix=0,numZ
          if (ldmodel(Zix,Nix).le.3) nendens(Zix,Nix)=60
          if (ldmodel(Zix,Nix).eq.4) then
            nendens(Zix,Nix)=55
            Edensmax(Zix,Nix)=150.
          endif
          if (ldmodel(Zix,Nix).ge.5) then
            nendens(Zix,Nix)=60
            Edensmax(Zix,Nix)=200.
          endif
        enddo
      enddo
      Nlast=0
      Ediscrete=0.
      scutoffdisc=1.
      delta=0.
      ldexist=.false.
      D0theo=0.
      D1theo=0.
      Tcrit=0.
      Ucrit=0.
      Econd=0.
      Scrit=0.
      aldcrit=0.
      do Nix=0,numN
        do Zix=0,numZ
          if (ldmodel(Zix,Nix).ge.4) then
            do i=0,numbar
              do ipar=-1,1,2
                do j=0,numJ
                  do nex=0,numdens
                    ldtottable(Zix,Nix,nex,i)=0.
                    ldtottableP(Zix,Nix,nex,ipar,i)=0.
                    ldtable(Zix,Nix,nex,j,ipar,i)=0.
                  enddo
                enddo
              enddo
            enddo
          endif
        enddo
      enddo
      rhogrid=0.
      if (.not.flagreaction) return
c
c Weak coupling parameters
c
c jcore: spin of level of core nucleus
c pcore: parity of level of core nucleus
c
      jcore=0.
      pcore=1
c
c Particle-hole state densities
c
c phmodel : particle-hole state density model
c phexist2: flag for existence of particle-hole state density table
c phtable2: particle-hole state density from table
c
      if (phmodel.eq.2) then
        RnJ=0.
        phexist2=.false.
        phexist1=.false.
        phtable2=0.
        phtable1=0.

c Configurations for microscopic particle-hole state densities
c
c Nphconf2 : number of 2-component particle-hole configurations
c Nphconf1 : number of 1-component particle-hole configurations
c phstring1: help variable
c phstring2: help variable
c ppitable : proton particle number from table
c hpitable : proton hole number from table
c pnutable : neutron particle number from table
c hnutable : neutron hole number from table
c pptable  : particle number from table
c hhtable  : hole number from table
c
        denfile=trim(path)//'density/ph/Fe.ph'
        inquire (file=denfile,exist=lexist)
        if (lexist) then
          Nphconf2=72
          Nphconf1=14
          open (unit=2,file=denfile,status='old')
          read(2,'(/////,9x,72(a4,5x),1x,14(a2,7x))')
     +      (phstring2(i),i=1,72),(phstring1(k),k=1,14)
          do 770 i=1,Nphconf2
            read(phstring2(i),'(4i1)') ppitable(i),hpitable(i),
     +        pnutable(i),hnutable(i)
  770     enddo
          do 780 i=1,Nphconf1
            read(phstring1(i),'(2i1)') pptable(i),hhtable(i)
  780     enddo
          close(2)
        else
          Nphconf2=0
          Nphconf1=0
        endif
      endif
c
c Giant resonance sum rules
c
c betagr : deformation parameter for giant resonance
c Egrcoll: energy of giant resonance
c Ggrcoll: width of giant resonance
c
      betagr=0.
      Egrcoll=0.
      Ggrcoll=0.
c
c Q-values and threshold energies
c
c Qres   : Q-value for residual nucleus
c Ethresh: threshold incident energy for residual nucleus
c
      Qres=0.
      Ethresh=0.
c
c Reaction flags
c
c flagwidth  : flag for width fluctuation calculation
c flagpreeq  : flag for pre-equilibrium calculation
c flagcompang: flag for compound angular distribution calculation
c flaggiant  : flag for collective contribution from giant resonances
c flagmulpre : flag for multiple pre-equilibrium calculation
c
c These flags will be reset later in subroutine energies.f, depending
c on the incident energy.
c
      flagwidth=.false.
      flagpreeq=.false.
      flagcompang=.false.
      flaggiant=.false.
      flagmulpre=.false.
c
c Flags for existence of files
c
c rpexist     : flag for existence of residual production cross section
c fisexist    : flag for existence of fission cross section
c rpisoexist  : flag for existence of isomeric residual production cross
c               section
c gamexist    : flag for existence of gamma production cross section
c numin,....  : maximal number of ejectile in channel description
c chanexist   : flag for existence of exclusive cross section
c spchanexist : flag for existence of exclusive spectra
c gamchanexist: flag for existence of exclusive discrete gamma-rays
c chanfisexist: flag for existence of exclusive fission cross section
c chanisoexist: flag for existence of exclusive isomeric cross section
c fpexist     : flag for existence of fission product
c Nrescue     : number of energies for adjustment factors
c urrexist    : flag for existence of URR
c lminU,lmaxU : minimal and maximal orbital angular momentum
c JminU,JmaxU : minimal and maximal total angular momentum
c flagurrendf : flag for URR info to ENDF
c
      rpexist=.false.
      fisexist=.false.
      recexist=.false.
      rpisoexist=.false.
      gamexist=.false.
      chanexist=.false.
      spchanexist=.false.
      recchanexist=.false.
      spfischanexist=.false.
      gamchanexist=.false.
      chanfisexist=.false.
      chanisoexist=.false.
      spexist1=.false.
      spexist2=.false.
      ddxexist1=.false.
      ddxexist2=.false.
      ddxexist3=.false.
      ddxexist4=.false.
      legexist=.false.
      angexist=.false.
      breakupexist=.false.
      if (nin0.eq.0) then
        fpexist=.false.
        fpaexist=.false.
        nubarexist=.false.
      endif
      lminU=numl
      lmaxU=0
      flagurrendf=.false.
      JminU=numJ
      JmaxU=0
      urrexist=.false.
      Nrescue=0
c
c PFNS energy grid
c
      Eout=0.
      degrid=0.001
      Epfns=0.
      NEpfns=0
      nen=0
 1110 Eout=Eout+degrid
      Eeps=Eout+1.e-4
      if (Eeps.gt.50.) goto 1120
      if (nen.eq.numpfns) goto 1120
      nen=nen+1
      Epfns(nen)=Eout
      if (Eeps.gt.0.01) degrid=0.01
      if (Eeps.gt.0.2) degrid=0.02
      if (Eeps.gt.0.5) degrid=0.05
      if (Eeps.gt.3.) degrid=0.1
      if (Eeps.gt.10.) degrid=0.5
      goto 1110
 1120 NEpfns=nen
      dEpfns=0.
      dEpfns(1)=Epfns(1)
      do nen=2,NEpfns-1
        dEpfns(nen)=0.5*(Epfns(nen+1)-Epfns(nen-1))
      enddo
      dEpfns(NEpfns)=0.5*(Epfns(NEpfns)-Epfns(NEpfns-1))
      return
      end
Copyright (C)  2023 A.J. Koning, S. Hilaire and S. Goriely
