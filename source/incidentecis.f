      subroutine incidentecis
c
c +---------------------------------------------------------------------
c | Author: Arjan Koning
c | Date  : November 13, 2004
c | Task  : ECIS calculation for incident energy
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      logical rotational,vibrational
      integer Zix,Nix,i1,i,ii,Z,A
      real    Ein
c
c ********************** Set ECIS input parameters *********************
c
c flaginccalc     : flag for new ECIS calculation for incident channel
c legendre        : logical for output of Legendre coefficients     
c Einc,Ein        : incident energy in MeV
c rmatch          : matching radius
c projmass,parmass: mass of projectile
c k0              : index of incident particle
c njmax           : maximal number of j-values in ECIS
c Atarget         : mass number of target nucleus
c onethird        : 1/3
c numl            : maximum l-value (set in talys.cmb)
c spin,parspin    : spin of incident particle
c resmass,tarmass : mass of target nucleus
c prodZ           : product of charges of projectile and target nucleus
c Ztarget         : charge number of target nucleus
c parZ            : charge number of particle
c Nband           : number of vibrational bands
c angbeg          : first angle
c anginc          : angle increment
c nangle          : number of angles
c angend          : last angle
c Zindex          : charge number index for residual nucleus
c Nindex          : neutron number index for residual nucleus
c
c Specific ECIS flags:
c ecis2(13)=T : output transmission coefficients 
c ecis2(14)=T : output of elastic angular distribution 
c ecis2(15)=T : output of Legendre coefficients        
c Extra ECIS-flags added by A. Koning:
c ecis2(10)=F : output of polarization or Rutherford cross section
c ecis2(20)=T : output of reaction cross section 
c ecis2(30)=T : output of direct inelastic cross section
c ecis2(40)=T : output of total and total elastic cross section 
c
      if (flaginccalc) 
     +  open (unit=9,status='unknown',file='ecisinc.inp')
      legendre=.true.
      Ein=Einc
      rmatch=0.
c
c We use a simple formula to estimate the required number of j-values:
c    njmax=2.4*k*R;  R=1.25*A**1/3 ; k=0.22*sqrt(m(in amu)E(in MeV))
c and we always take a minimum of njmax=20.
c
      projmass=parmass(k0)
      njmax=int(2.4*1.25*(real(Atarget)**onethird)*0.22*
     +  sqrt(projmass*Ein))
      njmax=max(njmax,20)
      njmax=min(njmax,numl)
      spin=parspin(k0)
      resmass=tarmass
      prodZ=real(Ztarget*parZ(k0))
      Nband=0
      if (k0.eq.1) then
        angbeg=0.
      else
        angbeg=0.00001
      endif
      anginc=180./nangle
      angend=180.
      Zix=Zindex(0,0,k0)
      Nix=Nindex(0,0,k0)
c
c Standard ECIS inputs for phenomenological optical potentials
c
c ecis1,ecis2: 100 input flags ('T' or 'F') for ECIS
c colltype   : type of collectivity (D, V or R)
c flagrot    : flag for use of rotational optical model per
c              outgoing particle, if available 
c rotational : flag for rotational input
c vibrational: flag for vibrational input      
c title      : title of ECIS input file
c ncoll      : number of nuclear states
c npp        : number of optical potentials
c iterm      : number of iterations
c idvib      : identifier for existence of vibrational state inside
c              rotational model        
c Elevel     : energy of level   
c tarspin    : spin of target nucleus
c tarparity  : parity of target nucleus     
c
      ecis1='FFFFFTFFFFFFFFFFFFFFFFFFTFFTFFFFFFFFFFFFFFFFFFFFFF'
      ecis2='FFFFFFFFFTFFTTTFTTTTFTTFTFFFFTFFFFFFFFFTFFFFFFFFFF'
c
c 1. Spherical nucleus
c
      if (colltype(Zix,Nix).eq.'S'.or..not.flagrot(k0)) then
        rotational=.false.
        vibrational=.false.
        title='Spherical optical model                           '
        ncoll=1
        npp=1
        iterm=1
        idvib(1)=1
        Elevel(1)=0.
        tarspin=0.
        tarparity='+'
      else
c
c 2. Deformed nucleus
c
c ndef       : number of collective levels           
c indexlevel : level index
c leveltype  : type of level (rotational (R) or vibrational (V))  
c vibband    : band number of level
c maxband    : highest vibrational level added to rotational model
c jdis       : spin of level
c cparity    : parity (character)
c parlev     : parity of level   
c Jlevel     : spin of level
c Plevel     : parity of level
c iph,iphonon: phonon (1 or 2)
c iband      : band number of level
c Jband,lband: angular momentum 
c Kmag       : magnetic quantum number  
c vibbeta    : vibrational deformation parameter 
c defpar     : deformation parameter
c flagstate  : flag for optical model potential for each excited state 
c
        iterm=0
        tarspin=jdis(Zix,Nix,0)
        tarparity=cparity(parlev(Zix,Nix,0))
        i1=0
        do 10 i=1,ndef(Zix,Nix)
          ii=indexlevel(Zix,Nix,i)
          if (leveltype(Zix,Nix,ii).ne.'V'.and.
     +      leveltype(Zix,Nix,ii).ne.'R') goto 10
          if (colltype(Zix,Nix).eq.'R'.and.
     +      vibband(Zix,Nix,i).gt.maxband) goto 10    
          i1=i1+1
          idvib(i1)=vibband(Zix,Nix,i)
          Elevel(i1)=edis(Zix,Nix,ii)
          Jlevel(i1)=jdis(Zix,Nix,ii)
          Plevel(i1)=cparity(parlev(Zix,Nix,ii))
          iph(i1)=iphonon(Zix,Nix,i)
          if (iph(i1).eq.2) ecis1(2:2)='T'
          Jband(i1)=lband(Zix,Nix,i)
          iband(i1)=vibband(Zix,Nix,i)
          Kmag(i1)=Kband(Zix,Nix,i)
          vibbeta(i1)=defpar(Zix,Nix,i)
          Nband=max(Nband,iband(i1))
   10   continue
        ncoll=i1
        if (flagstate) then
          npp=ncoll
        else
          npp=1
        endif
        ecis1(12:12)='T'
c
c 2a. Vibrational model
c
        if (colltype(Zix,Nix).eq.'V') then
          rotational=.false.
          vibrational=.true.
          title='Vibrational optical model                         '
          do 20 i=1,ncoll
            idvib(i)=0
   20     continue
        else
c
c 2b. Rotational model
c
c coulbar : Coulomb barrier  
c Nrotbeta: number of deformation parameters for rotational nucleus
c rotbeta : deformation parameters for rotational nucleus  
c nrot    : number of deformation parameters
c rotpar  : deformation parameters for rotational nucleus 
c iqm     : largest order of deformation
c iqmax   : maximum l-value of multipole expansion
c deftype : deformation length (D) or parameter (B)
c
          rotational=.true.
          vibrational=.false.
          ecis1(1:1)='T'
c
c Some input flags for ECIS are energy dependent for the rotational
c model. For rotational nuclei, the switch at 3 MeV needs to be made 
c according to Pascal Romain. 
c
          if (Ein.le.3.) then
            ecis1(21:21)='T'
            ecis1(42:42)='T'
          else
            ecis1(21:21)='F'
            ecis1(42:42)='F'
          endif
          if (k0.gt.1.and.Ein.le.0.05*coulbar(k0).and.
     +      Ein.le.2.*Elevel(ncoll)) Ein=0.1*Elevel(ncoll)
          Nrotbeta=nrot(Zix,Nix)
          do 30 i=1,Nrotbeta
            rotbeta(i)=rotpar(Zix,Nix,i)
   30     continue
          if (colltype(Zix,Nix).eq.'R') then
            title='Symmetric rotational optical model                '
            iqm=2*Nrotbeta
          else
            title='Asymmetric rotational optical model               '
            ecis1(3:3)='T'
            iqm=2*(Nrotbeta-1)
          endif
          iqmax=8
        endif
      endif
c
c ************** Calculate optical potential parameters ****************
c
c flagrel   : flag for relativistic kinematics
c optical   : subroutine for determination of optical potential
c flagoutomp: flag for output of optical model parameters
c ZZ        : charge number of residual nucleus
c AA        : mass number of residual nucleus
c parname   : name of particle
c nuc       : symbol of nucleus
c v,rv,.....: optical model parameters
c
      if (deftype(Zix,Nix).eq.'B') ecis1(6:6)='F'
      if (flagrel) ecis1(8:8)='T'
      call optical(Zix,Nix,k0,Ein)
c
c Output of optical model parameters, if requested.
c 
      if (flagoutomp) then
        Z=ZZ(0,0,k0)
        A=AA(0,0,k0)
        write(*,'(/"+++++++++ OPTICAL MODEL PARAMETERS FOR ",$)')
        write(*,'("INCIDENT CHANNEL ++++++++++")')
        write(*,'(/10x,a8," on ",i3,a2/)') parname(k0),A,nuc(Z)
        write(*,'(" Energy",4x,"V",5x,"rv",4x,"av",4x,"W",5x,"rw",$)')
        write(*,'(4x,"aw",4x,"Vd",3x,"rvd",3x,"avd",4x,"Wd",$)')
        write(*,'(3x,"rwd",3x,"awd",3x,"Vso",3x,"rvso",2x,"avso",$)')
        write(*,'(2x,"Wso",3x,"rwso",2x,"awso",2x,"rc",/)')
        write(*,'(f7.3,1x,6(f6.2,f6.3,f6.3),f6.3)')
     +    Ein,v,rv,av,w,rw,aw,vd,rvd,avd,wd,rwd,awd,vso,rvso,
     +    avso,wso,rwso,awso,rc
      endif
      if (.not.flaginccalc) return
c
c ******************* Write ECIS input file ****************************
c
c ecisinput: subroutine to create ECIS input file
c
      call ecisinput(Zix,Nix,k0,Ein,rotational,vibrational)
      write(9,'("fin")') 
      close (unit=9)
      legendre=.false.
c
c ************ ECIS calculation for incident energy ********************
c
c flagoutecis: flag for output of ECIS results
c ecis97t    : subroutine ecis97, adapted for TALYS
c ecisstatus : status of ECIS file  
c
      if (flagoutecis) then
        call ecis97t('ecisinc.inp  ','ecisinc.out  ','ecis97.inccs ',
     +    'ecis97.incres')  
      else
        call ecis97t('ecisinc.inp  ','/dev/null    ','ecis97.inccs ',
     +    'ecis97.incres')  
      endif
      open (unit=9,status='unknown',file='ecisinc.inp')
      close (unit=9,status=ecisstatus)
      return
      end
Copyright (C) 2004  A.J. Koning, S. Hilaire and M.C. Duijvestijn
