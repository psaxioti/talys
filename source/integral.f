        subroutine integral
c
c +---------------------------------------------------------------------
c | Author: Arjan Koning
c | Date  : July 30, 2011
c | Task  : Calculate effective cross section for integral spectrum
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      character*80  fluxfile
      integer       i,istat,nen,Nspec,Nxs,nen0
      real          Efluxup(0:numenin),Eflux(0:numenin),fspec(numenin),
     +              fluxsum,Exs(0:numenin),xs(0:numenin),xseff,Efl,Ea1,
     +              Eb1,xsa,xsb,xsf
c
c ********* Read integral spectrum from experimental database **********
c
c parsym  : symbol of particle
c k0      : index of incident particle
c Atarget : mass number of target nucleus
c nuc     : symbol of nucleus
c Ztarget : charge number of target nucleus
c Nflux   : number of reactions with integral data
c path    : directory containing structure files to be read
c lenpath : length of pathname
c fluxfile: file with experimental integral spectrum
c fluxname: name of experimental flux
c Nspec   : number of spectral energies
c fspec   : spectrum values
c
      open (unit=1,status='unknown',file='integral.dat')
      write(1,'("# ",a1," + ",i3,a2,
     +  ": Effective cross sections from integral data")')
     +  parsym(k0),Atarget,nuc(Ztarget)
      write(1,'("# Channel      Flux    Effective cross section (b)")')
      do 10 i=1,Nflux
        fluxfile=path(1:lenpath)//'flux/spectrum.'//fluxname(i)
        open (unit=2,status='old',file=fluxfile,iostat=istat)
        if (istat.eq.0) then
          read(2,'(1x,i4)') Nspec
          do 110 nen=Nspec,1,-1
            read(2,*) nen0,Efluxup(nen),fspec(nen)
  110     continue
          close (unit=2)
c
c Determine middle of histograms and energy bins. Calculate the
c integral of the flux for normalization.
c
          fluxsum=0.
          Efluxup(0)=Efluxup(1)
          do 120 nen=1,Nspec
            Eflux(nen)=0.5*(Efluxup(nen-1)+Efluxup(nen))*1.e-6
            fluxsum=fluxsum+fspec(nen)
  120     continue
        else
          write(*,'(" TALYS-error: integral spectrum file ",a80,
     +      "does not exist")') fluxfile
          stop
        endif
c
c *************** Read cross sections from TALYS output files **********
c
c Nxs: number of incident energies
c
        open (unit=3,status='old',file=xsfluxfile(i),iostat=istat)
        if (istat.eq.0) then
          Exs(0)=0.
          xs(0)=0.
          read(3,'(///,14x,i3,/)') Nxs
          do 210 nen=1,Nxs
            read(3,*) Exs(nen),xs(nen)
  210     continue
          close (unit=3)
        else
          write(*,'(" TALYS-error: cross section file ",a80,
     +      "does not exist")') xsfluxfile(i)
          stop
        endif
c
c *************** Calculate effective cross section by folding *********
c
c xseff; effective cross section
c
c Interpolate cross section grid on the grid of the integral spectrum
c
        xseff=0.
        do 310 nen=1,Nspec
          Efl=Eflux(nen)
          if (Efl.ge.Exs(0).and.Efl.le.Exs(Nxs)) then
            call locate(Exs,0,Nxs,Efl,nen0)
            Ea1=Exs(nen0)
            Eb1=Exs(nen0+1)
            xsa=xs(nen0)
            xsb=xs(nen0+1)
            call pol1(Ea1,Eb1,xsa,xsb,Efl,xsf)
            xseff=xseff+xsf*fspec(nen)
          endif
  310   continue
        xseff=0.001*xseff/fluxsum
        write(1,'(2a15,1p,e12.5)') xsfluxfile(i),fluxname(i),xseff
   10 continue
      close (unit=1)
      return
      end
Copyright (C) 2004  A.J. Koning, S. Hilaire and M.C. Duijvestijn
