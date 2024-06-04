      subroutine omppar(Zix,Nix)
c
c +---------------------------------------------------------------------
c | Author: Arjan Koning 
c | Date  : December 2, 2004
c | Task  : Optical model parameters
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      logical      lexist
      character*1  omptype
      character*4  ompchar
      character*72 optmodfile
      character*90 ompfile
      integer      Zix,Nix,Z,N,A,k,ia,i,nomp,ii,nen
c
c ************** Read optical model parameters from database ***********
c
c Zix         : charge number index for residual nucleus
c Nix         : neutron number index for residual nucleus
c ZZ,Z        : charge number of residual nucleus
c NN,N        : neutron number of residual nucleus
c AA,A        : mass number of residual nucleus
c omptype     : type of optical model (spherical or coupled)
c flaglocalomp: flag for local (y) or global (n) optical model
c ompchar     : help variable
c ompfile     : optical model parameter file
c optmodfileN : optical model parameter file for neutrons
c optmodfileP : optical model parameter file for protons
c path        : directory containing structure files to be read
c lenpath     : length of pathname
c ia          : mass number from level file
c nomp        : number of particles for optical model parameters
c ef          : Fermi energy
c rc0         : Coulomb radius
c rv0,av0     : real volume radius, diffuseness
c v1,v2,v3    : components for V
c w1,w2       : components for W
c rvd0,avd0   : real surface radius, diffuseness
c d1,d2,d3    : components for Wd
c rvso0,avso0 : real spin-orbit radius, diffuseness
c vso1,vso2   : components for Vso
c wso1,wso2   : components for Wso
c colltype    : type of collectivity (D, V or R)  
c
      Z=ZZ(Zix,Nix,0)
      N=NN(Zix,Nix,0)
      A=AA(Zix,Nix,0)
      omptype=' '
      if (.not.flaglocalomp) goto 100
      ompchar='z   '
      write(ompchar(2:4),'(i3.3)') Z
      do 10 k=1,2
        if (k.eq.1) then
          if (optmodfileN(Zix)(1:1).ne.' ') then
            ompfile=optmodfileN(Zix)
          else   
            ompfile=path(1:lenpath)//'optical/neutron/'//ompchar
          endif   
        else
          if (optmodfileP(Zix)(1:1).ne.' ') then
            ompfile=optmodfileP(Zix)
          else   
            ompfile=path(1:lenpath)//'optical/proton/'//ompchar
          endif   
        endif
        omptype=' '
        inquire (file=ompfile,exist=lexist)
        if (.not.lexist) goto 10
        open (unit=2,status='old',file=ompfile)
   20   read(2,'(4x,2i4,3x,a1)',end=60) ia,nomp,omptype
        if (A.ne.ia) then
          do 30 i=1,nomp
            do 30 ii=1,4
              read(2,'()')
   30     continue
          goto 20
        endif
        do 40 i=1,nomp
          read(2,'(4x,f7.2,f8.3)') ef(Zix,Nix,k),rc0(Zix,Nix,k)
          read(2,'(2f8.3,f6.1,f10.4,f9.6,f6.1,f7.1)') rv0(Zix,Nix,k),
     +      av0(Zix,Nix,k),v1(Zix,Nix,k),v2(Zix,Nix,k),v3(Zix,Nix,k),
     +      w1(Zix,Nix,k),w2(Zix,Nix,k)
          read(2,'(2f8.3,f6.1,f10.4,f7.2)') rvd0(Zix,Nix,k),
     +      avd0(Zix,Nix,k),d1(Zix,Nix,k),d2(Zix,Nix,k),d3(Zix,Nix,k)
          read(2,'(2f8.3,f6.1,f10.4,f6.1,f7.1)') rvso0(Zix,Nix,k),
     +      avso0(Zix,Nix,k),vso1(Zix,Nix,k),vso2(Zix,Nix,k),
     +      wso1(Zix,Nix,k),wso2(Zix,Nix,k)
c
c Only non-dispersive potentials are included in this version of TALYS.
c
          goto 60
   40   continue
   60   close (unit=2)
c
c Reduce d1 parameter (of Wd) in case of coupled-channels, unless
c already specified in OMP parameterization.
c
        if (colltype(Zix,Nix).ne.'S'.and.omptype.ne.'C')
     +    d1(Zix,Nix,k)=0.85*d1(Zix,Nix,k)
   10 continue
c
c ************************* Global optical model ***********************
c
c 1. Neutrons
c
c ompglobal: flag for use of global optical model
c onethird : 1/3
c twothird : 2/3
c
c Test if local OMP has been assigned.
c
  100 if (rv0(Zix,Nix,1).eq.0.) then
        ompglobal(Zix,Nix,1)=.true.
        ef(Zix,Nix,1)=-11.2814+0.02646*A
        rv0(Zix,Nix,1)=1.3039-0.4054*A**(-onethird)
        av0(Zix,Nix,1)=0.6778-1.487e-4*A
        v1(Zix,Nix,1)=59.30-21.0*real(N-Z)/A-0.024*A
        v2(Zix,Nix,1)=7.228e-3-1.48e-6*A
        v3(Zix,Nix,1)=1.994e-5-2.0e-8*A
        w1(Zix,Nix,1)=12.195+0.0167*A
        w2(Zix,Nix,1)=73.55+0.0795*A
        rvd0(Zix,Nix,1)=1.3424-0.01585*A**onethird
        avd0(Zix,Nix,1)=0.5446-1.656e-4*A
        d1(Zix,Nix,1)=16.0-16.0*real(N-Z)/A
        d2(Zix,Nix,1)=0.0180+3.802e-3/(1.+exp((A-156.)/8.0))
        d3(Zix,Nix,1)=11.5
        vso1(Zix,Nix,1)=5.922+0.0030*A
        vso2(Zix,Nix,1)=0.0040
        rvso0(Zix,Nix,1)=1.1854-0.647*A**(-onethird)
        avso0(Zix,Nix,1)=0.59    
        wso1(Zix,Nix,1)=-3.1
        wso2(Zix,Nix,1)=160.
        rc0(Zix,Nix,1)=0.
c
c Reduce d1 parameter (of Wd) in case of coupled-channels, unless
c already specified in OMP parameterization.
c
        if (colltype(Zix,Nix).ne.'S'.and.omptype.ne.'C')
     +    d1(Zix,Nix,1)=0.85*d1(Zix,Nix,1)
      endif
c
c 2. Protons
c
      if (rv0(Zix,Nix,2).eq.0.) then
        ompglobal(Zix,Nix,2)=.true.
        ef(Zix,Nix,2)=-8.4075+0.01378*A
        rv0(Zix,Nix,2)=1.3039-0.4054*A**(-onethird)
        av0(Zix,Nix,2)=0.6778-1.487e-4*A
        v1(Zix,Nix,2)=59.30+21.0*real(N-Z)/A-0.024*A
        v2(Zix,Nix,2)=7.067e-3+4.23e-6*A
        v3(Zix,Nix,2)=1.729e-5+1.136e-8*A
        w1(Zix,Nix,2)=14.667+0.009629*A
        w2(Zix,Nix,2)=73.55+0.0795*A
        rvd0(Zix,Nix,2)=1.3424-0.01585*A**onethird
        avd0(Zix,Nix,2)=0.5187+5.205e-4*A
        d1(Zix,Nix,2)=16.0+16.0*real(N-Z)/A
        d2(Zix,Nix,2)=0.0180+3.802e-3/(1.+exp((A-156.)/8.0))
        d3(Zix,Nix,2)=11.5
        vso1(Zix,Nix,2)=5.922+0.0030*A
        vso2(Zix,Nix,2)=0.0040
        rvso0(Zix,Nix,2)=1.1854-0.647*A**(-onethird)
        avso0(Zix,Nix,2)=0.59    
        wso1(Zix,Nix,2)=-3.1
        wso2(Zix,Nix,2)=160.
        rc0(Zix,Nix,2)=1.198+0.697*A**(-twothird)+12.994*A**(-5./3.)
c
c Reduce d1 parameter (of Wd) in case of coupled-channels, unless
c already specified in OMP parameterization.
c
        if (colltype(Zix,Nix).ne.'S'.and.omptype.ne.'C')
     +    d1(Zix,Nix,2)=0.85*d1(Zix,Nix,2)
      endif
c
c ************** Optical model file from user input file ***************
c
c optmod,optmodfile: file with optical model parameters
c omplines         : number of lines on optical model file
c eomp             : energies on optical model file
c vomp             : optical model parameters from file
c
      do 210 k=1,6
        optmodfile=optmod(Zix,Nix,k)
        if (optmodfile(1:1).ne.' ') then
          open (unit=2,status='old',file=optmodfile)   
          read(2,'(8x,i4)') omplines(k)
          eomp(k,0)=0.
          do 220 nen=1,omplines(k)
            read(2,'(f7.3,6(f7.3,2f6.3),f6.3)',err=300) eomp(k,nen),
     +        (vomp(k,nen,ii),ii=1,19)    
  220     continue
        endif
  210 continue
      return
  300 write(*,'("TALYS-error: Format error in ",a72)') optmodfile
      stop
      end
Copyright (C) 2004  A.J. Koning, S. Hilaire and M.C. Duijvestijn
