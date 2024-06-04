      subroutine densitypar(Zix,Nix)
c
c +---------------------------------------------------------------------
c | Author: Arjan Koning 
c | Date  : December 2, 2004
c | Task  : Level density parameters
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      logical          lexist
      character*4      denchar
      character*90     denfile
      logical          inpdeltaW,inpalimit
      integer          Zix,Nix,Z,N,A,ia,ibar,nmax,oddZ,oddN
      real             gammaphen,ald,Spair,fU,factor,argum
      double precision mliquid
c
c ******* Total level density parameter a for Fermi Gas model **********
c
c Zix        : charge number index for residual nucleus
c Nix        : neutron number index for residual nucleus
c ZZ,Z       : charge number of residual nucleus
c NN,N       : neutron number of residual nucleus
c AA,A       : mass number of residual nucleus
c gammaphen  : gamma-constant for asymptotic level density parameter
c gammashell1: gamma-constant for asymptotic level density parameter
c gammashell2: gamma-constant for asymptotic level density parameter
c onethird   : 1/3
c
c Energy dependent level density systematics from M. Duijvestijn.
c
      Z=ZZ(Zix,Nix,0)
      N=NN(Zix,Nix,0)
      A=AA(Zix,Nix,0)
      gammaphen=gammashell1/(A**onethird)+gammashell2
c
c Read experimental values from level density parameter file
c
c alev,ald: level density parameter
c flagasys: flag for all level density parameters a from systematics
c denfile : level density parameter file
c ia      : mass number from level density table
c
c Level density parameters from the table can always be overruled
c by values given in the input file. With flagasys, all experimental
c level density parameters a from the table can be overruled by the 
c systematics.
c
      if (alev(Zix,Nix).eq.0..and..not.flagasys) then
        denchar='z   '
        write(denchar(2:4),'(i3.3)') Z
        denfile=path(1:lenpath)//'density/ground/fermi/'//denchar
        inquire (file=denfile,exist=lexist)
        if (.not.lexist) goto 20
        open (unit=2,status='old',file=denfile)
   10   read(2,'(4x,i4,f8.4)',end=20) ia,ald
        if (A.ne.ia) goto 10
        alev(Zix,Nix)=ald
   20   close (unit=2)
      endif
c
c ****** Initialization of other default level density parameters ******
c
c ibar     : fission barrier
c nfisbar  : number of fission barrier parameters
c Nlast    : last discrete level
c nlev     : number of excited levels for nucleus
c nfistrrot: number of rotational transition states for barrier 
c Ntop     : highest discrete level for temperature matching
c nmax     : help variable
c Nlow     : lowest discrete level for temperature matching
c
      do 110 ibar=0,nfisbar(Zix,Nix)
        if (ibar.eq.0) then
          Nlast(Zix,Nix,ibar)=nlev(Zix,Nix)
        else
          Nlast(Zix,Nix,ibar)=max(nfistrrot(Zix,Nix,ibar),1)
        endif
        if (ibar.eq.0.and.Ntop(Zix,Nix,ibar).eq.-1) then
          denchar='z   '
          write(denchar(2:4),'(i3.3)') Z
          denfile=path(1:lenpath)//'density/ground/nmax/'//denchar
          inquire (file=denfile,exist=lexist)
          if (.not.lexist) goto 120
c
c We allow a maximum of Ntop=30
c
          open (unit=2,status='old',file=denfile)
  130     read(2,'(4x,2i4)',end=120) ia,nmax
          if (A.ne.ia) goto 130
          Ntop(Zix,Nix,0)=min(nmax,30)
  120     close (unit=2)
        endif
        if (Ntop(Zix,Nix,ibar).eq.-1) Ntop(Zix,Nix,ibar)=
     +    Nlast(Zix,Nix,ibar)
        if (Nlow(Zix,Nix,ibar).eq.-1) Nlow(Zix,Nix,ibar)=2
        if (Ntop(Zix,Nix,ibar).le.2) Nlow(Zix,Nix,ibar)=0
  110 continue
c
c There are many input possibilities for the energy dependent level
c density parameter of the Ignatyuk formula. The required parameters
c are alev, alimit, gammald and deltaW. The Ignatyuk formula implies
c that they can not all be given at the same time in the input file,
c or, for alev, in a table. We re-determine alev if deltaW, gammald and 
c alimit are given by input.
c
c deltaW   : shell correction in nuclear mass 
c alimit   : asymptotic level density parameter
c gammald  : gamma-constant for asymptotic level density parameter
c inpdeltaW: logical to determine existence of input value for deltaW
c nucmass  : mass of nucleus
c mliquid  : function for liquid drop mass
c amu      : atomic mass unit in MeV
c oddZ,oddN: help variables
c inpalimit: logical to determine existence of input value for alimit
c alphald  : alpha-constant for asymptotic level density parameter
c betald   : beta-constant for asymptotic level density parameter
c twothird : 2/3
c pair     : total pairing correction
c Spair    : help variable
c S        : separation energy per particle
c
      if (alev(Zix,Nix).ne.0.and.deltaW(Zix,Nix,0).ne.0..and.
     +  alimit(Zix,Nix).ne.0.and.gammald(Zix,Nix).ne.-1.) 
     +  alev(Zix,Nix)=0.
      inpdeltaW=.true.
      if (deltaW(Zix,Nix,0).eq.0.) then
        inpdeltaW=.false.
        deltaW(Zix,Nix,0)=real((nucmass(Zix,Nix)-mliquid(Z,A))*amu)
      endif
      oddZ=mod(Z,2)
      oddN=mod(N,2)
      inpalimit=.true.
      if (alimit(Zix,Nix).eq.0.) then
        inpalimit=.false.
        alimit(Zix,Nix)=alphald*A+betald*(A**twothird)
      endif
      if (pair(Zix,Nix).eq.0.)
     +  pair(Zix,Nix)=(2.-oddZ-oddN)*12./sqrt(real(A))
      Spair=S(Zix,Nix,1)-pair(Zix,Nix)
      Spair=max(Spair,1.)
c
c 1. If no experimental level density parameter is available, 
c    i.e. as determined from the neutron resonance spacing, use the
c    Ignatyuk formula to derive the level density parameter at the 
c    separation energy.
c
c fU,factor: help variables
c
      if (alev(Zix,Nix).eq.0.) then
        if (gammald(Zix,Nix).eq.-1.) gammald(Zix,Nix)=gammaphen
        fU=1.-exp(-gammald(Zix,Nix)*Spair)
        factor=1.+fU*deltaW(Zix,Nix,0)/Spair
        alev(Zix,Nix)=alimit(Zix,Nix)*factor
        alev(Zix,Nix)=max(alev(Zix,Nix),1.)
      else
c
c 2. If an experimental level density parameter is available, then we
c    impose the extra boundary boundary condition that it should be 
c    equal to the energy dependent level density parameter at the
c    neutron separation energy. There are various possibilities:
c    If no further input is given, this eliminates gammald as a free
c    parameter. If all level density parameters are given, we first 
c    readjust gammald.
c
c argum: help variable
c
        if (gammald(Zix,Nix).eq.-1.) then
          argum=1.-Spair/deltaW(Zix,Nix,0)*
     +      (alev(Zix,Nix)/alimit(Zix,Nix)-1.)
          if (argum.gt.0..and.argum.lt.1.) then
            gammald(Zix,Nix)=-1./Spair*log(argum)
          else
c
c If gammald can not be solved or is unphysical (this may happen when 
c both deltaW and alev come from experiment), we re-adjust the shell 
c correction.
c
            gammald(Zix,Nix)=gammaphen
            fU=1.-exp(-gammald(Zix,Nix)*Spair)
            factor=alev(Zix,Nix)/alimit(Zix,Nix)-1.
            deltaW(Zix,Nix,0)=Spair*factor/fU
          endif
        else
c
c Determine alimit if alev, gammald and deltaW are given by input.
c
          if (inpdeltaW.or..not.inpalimit) then
            fU=1.-exp(-gammald(Zix,Nix)*Spair)
            factor=1.+fU*deltaW(Zix,Nix,0)/Spair
            alimit(Zix,Nix)=alev(Zix,Nix)/factor
          endif
c
c Determine deltaW if alev, gammald and alimit are given by input.
c
          if (inpalimit) then
            fU=1.-exp(-gammald(Zix,Nix)*Spair)
            factor=alev(Zix,Nix)/alimit(Zix,Nix)-1.
            deltaW(Zix,Nix,0)=Spair*factor/fU
          endif
        endif
      endif
c
c ************** Single-particle level density parameter g *************
c
c g           : single-particle level density parameter
c pi2         : pi**2
c Kph         : constant for single-particle level density parameter
c               (g=A/Kph)          
c gp          : single-particle proton level density parameter
c gn          : single-particle neutron level density parameter
c
c One component
c
      if (g(Zix,Nix).eq.0.) then
        g(Zix,Nix)=A/Kph
      endif
c
c Two component
c
      if (gp(Zix,Nix).eq.0.) then
        gp(Zix,Nix)=Z/Kph
      endif
      if (gn(Zix,Nix).eq.0.) then
        gn(Zix,Nix)=N/Kph
      endif
      return
      end
Copyright (C) 2004  A.J. Koning, S. Hilaire and M.C. Duijvestijn
