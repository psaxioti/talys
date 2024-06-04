        subroutine integral
c
c +---------------------------------------------------------------------
c | Author: Arjan Koning
c | Date  : December 5, 2020
c | Task  : Calculate effective cross section for integral spectrum
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      integer numP
      parameter (numP=1000000)
      logical        lexist
      character*1    isostring(numflux)
      character*3    Astring,ext
      character*8    reac,reacstr(numflux)
      character*132  word(40)
      character*132  xsfile,fluxfile
      character*200  line
      integer        i,istat,nen,Nspec,Nxs,nen0,is,k,L
      real           Eflux(0:numenin),fspec(numenin),
     +               fluxsum,Exs(0:numP),xs(0:numP),xseff,Efl,Ea1,
     +               Eb1,xsa,xsb,xsf,xsexp,ratio,igr,Eup,elow,dum,sp
c
c ********* Read reaction channels with integral cross sections ********
c
c This is for the case where a simple 'integral y' is given in the input 
c file, i.e. no explicit information per case
c
      if (Nflux.eq.0.and.Liso.eq.0) then
        Astring='   '
        write(Astring(1:3),'(i3.3)') Atarget
        xsfile=trim(path)//'integral/sacs/'//trim(Starget)//Astring//
     +    '.sacs'
        open (unit=2,file=xsfile,status='old',iostat=istat)
        if (istat.eq.0) then
          i=1
          do 
            read(2,'(a)',iostat=istat) line
            if (istat.ne.0) exit
            call getkeywords(line(1:132),word)
            reacstr(i)=trim(word(3))
            isostring(i)=trim(word(4))
            fluxname(i)=trim(word(5))
            read(line(41:49),*) integralexp(i)
c           read(line(51:59),*) dintegralexp(i)
            i=i+1
          enddo 
          close (unit=2)
          Nflux=i-1
        endif
      endif
c
c Determine corresponding TALYS output file
c
      do i=1,Nflux
        reac=''
        ext='tot'
        if (trim(reacstr(i)).eq.'n,g') reac='xs000000'
        if (trim(reacstr(i)).eq.'n,p') reac='xs010000'
        if (trim(reacstr(i)).eq.'n,d') reac='xs001000'
        if (trim(reacstr(i)).eq.'n,t') reac='xs000100'
        if (trim(reacstr(i)).eq.'n,h') reac='xs000010'
        if (trim(reacstr(i)).eq.'n,a') reac='xs000001'
        if (trim(reacstr(i)).eq.'n,2n') reac='xs200000'
        if (trim(reacstr(i)).eq.'n,np') reac='xs110000'
        if (trim(reacstr(i)).eq.'n,na') reac='xs100001'
        if (trim(reacstr(i)).eq.'n,2a') reac='xs000002'
        if (trim(reacstr(i)).eq.'n,2na') reac='xs200001'
        if (trim(reacstr(i)).eq.'n,3n') reac='xs300000'
        if (trim(reacstr(i)).eq.'n,4n') reac='xs400000'
        if (trim(reacstr(i)).eq.'n,pa') reac='xs010001'
        if (trim(reacstr(i)).eq.'n,da') reac='xs001001'
        if (trim(reacstr(i)).eq.'n,nt') reac='xs100100'
        if (trim(reacstr(i)).eq.'n,n2a') reac='xs100002'
        if (trim(reacstr(i)).eq.'n,3na') reac='xs300001'
        if (isostring(i).eq.'g') ext='L00'
        if (isostring(i).eq.'m') then
          ext='L  '
          write(ext(2:3),'(i2.2)') L
          do L=1,numlev
            write(ext(2:3),'(i2.2)') L
            inquire(file=reac//'.'//ext,exist=lexist)
            if (lexist) exit
          enddo
        endif
        xsfluxfile(i)=reac//'.'//ext
      enddo
c
c ********* Read integral spectrum from experimental database **********
c
c parsym  : symbol of particle
c k0      : index of incident particle
c Atarget : mass number of target nucleus
c Starget : symbol of target nucleus
c Ztarget : charge number of target nucleus
c Nflux   : number of reactions with integral data
c path    : directory containing structure files to be read
c Eflux   : energy of bin for flux
c Efl     : energy of bin for flux
c fluxfile: file with experimental integral spectrum
c fluxname: name of experimental flux
c Nspec   : number of spectral energies
c fspec   : spectrum values
c
      open (unit=1,file='integral.dat',status='replace')
      write(1,'("# ",a1," + ",a,
     +  ": Effective cross sections from integral data")')
     +  parsym(k0),trim(targetnuclide)
      write(1,'("# Channel      Spectrum        Eff. xs (b)",
     +  " Exp. xs (b)         Ratio")')
      do 10 i=1,Nflux
        fluxfile=trim(path)//'integral/spectra/'//trim(fluxname(i))//
     +    '.txt'
        open (unit=2,file=fluxfile,status='old',iostat=istat)
        if (istat.eq.0) then
          read(2,'()')
          read(2,'()')
          nen=0
          fluxsum=0
          do
            read(2,*,iostat=istat) igr,Eup,elow,dum,sp
            if (istat.ne.0) exit
            nen=nen+1
            fspec(nen)=sp
            Eflux(nen)=0.5*(Eup+Elow)*1.e-6
            fluxsum=fluxsum+sp
          enddo
          close (unit=2)
          Nspec=nen
        else
          write(*,'(" TALYS-warning: integral spectrum file ",a,
     +      " does not exist")') trim(fluxfile)
          goto 10
        endif
c
c *************** Read cross sections from TALYS output files **********
c
c Exs: energy of cross section file
c Nxs: number of incident energies
c
        is=1
  200   open (unit=3,file=xsfluxfile(i),status='old',iostat=istat)
        if (istat.eq.0) then
          Exs(0)=0.
          xs(0)=0.
          read(3,'(///,14x,i6,/)') Nxs
          do 210 nen=1,Nxs
            read(3, *, iostat=istat ) Exs(nen), xs(nen)
            if (istat == -1) exit
  210     continue
          close (unit=3)
        else
          if (is.lt.numlev) then
            is=is+1
            do 220 k=1,20
              if (xsfluxfile(i)(k:k).eq.'L') then
                write(xsfluxfile(i)(k+1:k+2),'(i2.2)') is
                goto 200
              endif
  220       continue
          endif
          write(*,'(" TALYS-warning: cross section file ",a,
     +      " does not exist")') trim(xsfluxfile(i))
          goto 10
        endif
c
c *************** Calculate effective cross section by folding *********
c
c xseff: effective cross section
c xsexp: experimental effective cross section
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
        xsexp=integralexp(i)
        if (xsexp.gt.0.) then
          ratio=xseff/xsexp
        else
          ratio=1.
        endif
        write(1,'(2a15,2es12.5,f15.5)') xsfluxfile(i),fluxname(i),
     +    xseff,xsexp,ratio
   10 continue
      close (unit=1)
      return
      end
Copyright (C)  2016 A.J. Koning, S. Hilaire and S. Goriely
