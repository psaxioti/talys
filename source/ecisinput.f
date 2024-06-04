      subroutine ecisinput(Zix,Nix,kopt,e,rotational,vibrational)
c
c +---------------------------------------------------------------------
c | Author: Arjan Koning
c | Date  : November 13, 2004
c | Task  : Create ECIS input file
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      logical      rotational,vibrational
      character*32 Eformat
      integer      Zix,Nix,kopt,i
      real         e,eopt
c
c *********************** Write standard input *************************
c
c Zix        : charge number index for residual nucleus
c Nix        : neutron number index for residual nucleus
c kopt       : optical model index for fast particle
c e          : incident energy in MeV
c rotational : flag for rotational input
c vibrational: flag for vibrational input
c title      : title of ECIS input file
c ecis1,ecis2: 100 input flags ('T' or 'F') for ECIS
c ncoll      : number of nuclear states
c njmax      : maximal number of j-values in ECIS
c iterm      : number of iterations   
c npp        : number of optical potentials    
c rmatch     : matching radius
c legendre   : logical for output of Legendre coefficients
c Eformat    : format string
c tarspin    : spin of target nucleus
c idvib      : identifier for existence of vibrational state inside 
c              rotational model
c tarparity  : parity of target nucleus
c spin       : spin of incident particle
c projmass   : mass of projectile      
c resmass    : mass of target nucleus
c prodZ      : product of charges of projectile and target nucleus
c Jlevel     : spin of level
c Plevel     : parity of level
c Elevel     : energy of level
c iph        : phonon (1 or 2)
c iband      : band number of level
c Nband      : number of vibrational bands  
c Kmag       : magnetic quantum number
c vibbeta    : vibrational deformation parameter 
c iqm        : largest order of deformation
c iqmax      : maximum l-value of multipole expansion
c rotbeta    : deformation parameters for rotational nucleus
c Nrotbeta   : number of deformation parameters for rotational nucleus
c
      write(9,'(a72)') title
      write(9,'(a50)') ecis1
      write(9,'(a50)') ecis2
      write(9,'(4i5)') ncoll,njmax,iterm,npp
      write(9,'(10x,f10.5,10x,"    1.e-10    1.e-10    1.e-30")') rmatch
      if (legendre) write(9,'()')
      if (e.ge.0.01) then
        Eformat='(f5.2,2i2,a1,5f10.5)'
      else
        Eformat='(f5.2,2i2,a1,1p,e10.3,0p,4f10.5)'
      endif
      write(9,fmt=Eformat) tarspin,idvib(1),1,tarparity,e,spin,
     +  projmass,resmass,prodZ
      if (vibrational) write(9,'()')
      do 10 i=2,ncoll
        write(9,'(f5.2,2i2,a1,5f10.5)') Jlevel(i),idvib(i),min(i,npp),
     +    Plevel(i),Elevel(i)
        if (vibrational) then
          write(9,'(3i5)') iph(i),iband(i),1
        else
          if (idvib(i).ge.1) write(9,'(2i5)') iph(i),iband(i)
          if (ecis1(3:3).eq.'T') write(9,'()')
        endif
   10 continue
      do 20 i=1,Nband
        write(9,'(2i5,f10.5)') Jband(i),Kmag(i),vibbeta(i)
   20 continue
      if (rotational) then
        write(9,'(2i5,f10.1)') iqm,iqmax,tarspin
        write(9,'(7f10.5)') (rotbeta(i),i=1,Nrotbeta)
      endif
c
c The optical model parameters can be calculated at the ground state
c and every excited state.
c
c eopt         : incident energy
c optical      : subroutine for determination of optical potential
c v,rv,av      : real volume potential, radius, diffuseness
c vd,rvd,avd   : real surface potential, radius, diffuseness
c w,rw,aw      : imaginary volume potential, radius, diffuseness
c wd,rwd,awd   : imaginary surface potential, radius, diffuseness
c vso,rvso,avso: real spin-orbit potential, radius, diffuseness
c wso,rwso,awso: imaginary spin-orbit potential, radius, diffuseness
c rc           : Coulomb radius
c angbeg       : first angle 
c anginc       : angle increment
c angend       : last angle 
c
      do 30 i=1,npp
        eopt=e-real(Elevel(i)*(resmass+projmass)/resmass)
        call optical(Zix,Nix,kopt,eopt)
        write(9,'(3f10.5)') v,rv,av
        write(9,'(3f10.5)') w,rw,aw
        write(9,'(3f10.5)') vd,rvd,avd
        write(9,'(3f10.5)') wd,rwd,awd
        write(9,'(3f10.5)') vso,rvso,avso
        write(9,'(3f10.5)') wso,rwso,awso
        write(9,'(3f10.5)') rc,0.,0.
        write(9,'(3f10.5)') 0.,0.,0.
   30 continue
      write(9,'(3f10.5)') angbeg,anginc,angend
      return
      end
Copyright (C) 2004  A.J. Koning, S. Hilaire and M.C. Duijvestijn
