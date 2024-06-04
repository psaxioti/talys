      subroutine astroout
c
c +---------------------------------------------------------------------
c | Author: Stephane Goriely, Stephane Hilaire, Arjan Koning
c | Date  : April 7, 2022
c | Task  : Output of astrophysical reaction rates
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      integer       numZ1,numN1
      parameter     (numZ1=numZ+1,numN1=numN+1)
      logical       lexist
      character*7   machar
      character*132 mafile
      character*13  astrofile
      character*1   targetparity,yesno
      integer       Acomp,Zcomp,Ncomp,i,Z,A,ires,iresprod,iwriterp,ia,
     +              zrespro(numZ1*numN1),arespro(numZ1*numN1),nex,
     +              maxAastro
      real          rateastrorp(numT,numZ1*numN1),rateastroeps,
     +              ratio,macs,dmacs,xs,dxs,branch
c
c ************************ Output of reaction rates ********************
c
c numZ1      : numZ + 1
c numN1      : numN + 1
c Acomp      : mass number index for compound nucleus
c maxAastro  : maximal number of nucleons away from initial compound
c              nucleus for astrophysical calculations
c Zcomp      : proton number index for compound nucleus
c maxZastro  : maximal number of protons away from initial compound
c              nucleus for astrophysical calculations
c Ncomp      : neutron number index for compound nucleus
c maxNastro  : maximal number of neutrons away from initial compound
c              nucleus for astrophysical calculations
c nTmax      : effective number of temperatures for Maxwellian
c rateastro  : thermonuclear reaction rate factor
c macsastro  : thermonuclear reaction cross section
c flagastrogs: flag for calculation of astrophysics reaction rate with
c              target in ground state only
c ZZ,Z       : charge number of residual nucleus
c AA,A       : mass number of residual nucleus
c nuc        : symbol of nucleus
c T9         : Temperature grid in 10**9 K
c partf      : integrated partition function
c edis       : energy of level
c astrofile: file with astro results
c targetparity: parity of target
c branch     : branching ratio to a given excited state
c
      write(*,'(/" 8. Thermonuclear reaction rates")')
      maxAastro=maxZastro+maxNastro
      do Acomp=0,maxAastro
        do 20 Zcomp=0,maxZastro
          Ncomp=Acomp-Zcomp
          if (Ncomp.lt.0.or.Ncomp.gt.maxNastro) goto 20
          do i=1,nTmax
            if (rateastro(Zcomp,Ncomp,i).gt.0.) goto 40
          enddo
          goto 20
   40     Z=ZZ(Zcomp,Ncomp,0)
          A=AA(Zcomp,Ncomp,0)
         if (nTmax.eq.1) then
            write(*,'(/" Reaction rate for Z=",i3," A=",i3," (",i3,a2,
     +        ") at <E>=",f8.5," MeV (Excited States Contribution : ",
     +        a1,")"/)') Z,A,A,nuc(Z),astroE,yesno(flagastrogs)
          else
            write(*,'(/" Reaction rate for Z=",
     +        i3," A=",i3," (",i3,a2,")"/)') Z,A,A,nuc(Z)
          endif
          write(*,'("    T        G(T)        Rate       MACS "/)')
          do i=1,nTmax
            write(*,'(1x,f8.4,3es12.5)') T9(i),partf(i),
     +        rateastro(Zcomp,Ncomp,i),macsastro(Zcomp,Ncomp,i)
          enddo
          if (flagastroex) then
            do nex=0,Nlast(Zcomp,Ncomp,0)
              if (nex.gt.0.and.tau(Zcomp,Ncomp,nex).eq.0.) cycle
              write(*,'(/" Reaction rate for Z=",
     +          i3," A=",i3," (",i3,a2,") to the excited states L",
     +          i2.2," at E=",f12.5," MeV",/)') Z,A,A,nuc(Z),nex,
     +          edis(Zcomp,Ncomp,nex)
              write(*,'("    T       Rate         MACS      ",
     +          "Branching"/)')
              do i=1,nTmax
                branch=0.
                if (rateastro(Zcomp,Ncomp,i).gt.0.) branch=
     +            rateastroex(Zcomp,Ncomp,i,nex)/
     +            rateastro(Zcomp,Ncomp,i)
                write(*,'(1x,f8.4,3es12.5)') T9(i),
     +            rateastroex(Zcomp,Ncomp,i,nex),
     +            macsastroex(Zcomp,Ncomp,i,nex),branch
              enddo
            enddo
          endif
          if (flagracap.and.Zcomp.eq.0.and.Acomp.eq.0) then
            write(*,'(/"    T      Rate(Eq)    Rate(DC)  ",
     +        "  MACS(Eq)    MACS(DC)  "/)')
            do i=1,nTmax
              write(*,'(1x,f8.4,4es12.5)') T9(i),
     +          rateastro(Zcomp,Ncomp,i)-rateastroracap(i),
     +          rateastroracap(i),
     +          macsastro(Zcomp,Ncomp,i)-macsastroracap(i),
     +          macsastroracap(i)
            enddo
          endif
   20   continue
      enddo
c
c Comparison with experimental MACS at 30 keV
c
c machar : part of filename for MACS
c mafile : file with MACS
c dmacs: uncertainty of MACS
c rateastrorp: rate for astro residual production
c
      if (astroE.ge.0.029.and.astroE.le.0.031) then
        Z=Ztarget
        A=Atarget
        ratio=0.
        macs=0.
        dmacs=0.
        machar=trim(nuc(Z))//'.macs'
        mafile=trim(path)//'gamma/macs/'//machar
        inquire (file=mafile,exist=lexist)
        if (lexist) then
          open (unit=2,file=mafile,status='old')
   70     read(2,'(4x,i4,2(es12.4))',end=80) ia,xs,dxs
          if (A.ne.ia) goto 70
          macs=xs
          dmacs=dxs
   80     close (unit=2)
        endif
        if (macs.gt.0.) ratio=macsastro(0,0,1)/macs
        astrofile='macs.g'
        open (unit=1,file=astrofile,status='replace')
        write(1,'("# Z   A     MACS(mb)    Exp(mb)     dEXP",
     +    "(mb)   MACS/Exp")')
        write(1,'(2i4,4es12.5)') Z,A,macsastro(0,0,1),macs,dmacs,
     +    ratio
        close (unit=1)
      endif
c
c Write results to separate files
c
c rateastroeps: cutoff value
c zresprod    : help variable
c Atarget     : mass number of target nucleus
c Ztarget     : charge number of target nucleus
c parsym      : symbol of particle
c k0          : index of incident particle
c flagfission : flag for fission
c
      astrofile='astrorate.g'
      open (unit=1,file=astrofile,status='replace')
      if (nTmax.eq.1) then
        write(1,'("# Reaction rate for ",a,"(",a1,",g) at <E>=",
     +    f8.5," MeV (Excited States Contribution : ",a1,")")')
     +    trim(targetnuclide),parsym(k0),astroE,yesno(flagastrogs)
      else
        write(1,'("# Reaction rate for ",a,"(",a1,",g)")')
     +    trim(targetnuclide),parsym(k0)
      endif
      write(1,'("#    T       Rate       MACS       G(T)")')
      do i=1,nTmax
        write(1,'(f8.4,3es12.5)') T9(i),rateastro(0,0,i),
     +    macsastro(0,0,i),partf(i)
      enddo
      close (unit=1)
c
c output partial rates(n,g) to given excited states in a specific file
c
      if (.not.flagastroex) goto 240
      do nex=1,Nlast(0,0,0)
        if (tau(0,0,nex).ne.0.) goto 230
      enddo
      goto 240
  230 do nex=0,Nlast(0,0,0)
        if (nex.eq.0.or.tau(0,0,nex).ne.0.) then
          astrofile='astrorate.g'
          write(astrofile(12:13),'(i2.2)') nex
          open (unit=1,file=astrofile,status='replace')
          if (nTmax.eq.1) then
            write(1,'("# Reaction rate for ",a,"(",a1,",g) at <E>=",
     +        f8.5," MeV (Exc. States Cont. : ",a1,
     +        ") to final level L",i2.2," at E=",f12.5," MeV")')
     +        trim(targetnuclide),parsym(k0),astroE,
     +        yesno(flagastrogs),nex,edis(0,0,nex)
          else
            write(1,'("# Reaction rate for Z=",i3," A=",i3,2x,a,"(",
     +        a1,",g) to the excited state L",i2.2," at E=",f12.5,
     +        " MeV")') Ztarget,Atarget,trim(targetnuclide),parsym(k0),
     +        nex,edis(0,0,nex)
          endif
          write(1,'("#   T9      Rate         MACS      Branching")')
          do i=1,nTmax
            write(1,'(1x,f8.4,3es12.5)') T9(i),
     +        rateastroex(0,0,i,nex),
     +        macsastroex(0,0,i,nex),
     +        rateastroex(0,0,i,nex)/rateastro(0,0,i)
          enddo
          close (unit=1)
        endif
      enddo
  240 continue
      astrofile='astrorate.p'
      open (unit=1,file=astrofile,status='replace')
      if (nTmax.eq.1) then
        write(1,'("# Reaction rate for ",a,"(",a1,",p) at <E>=",
     +    f8.5," MeV (Excited States Contribution : ",a1,")")')
     +    trim(targetnuclide),parsym(k0),astroE,yesno(flagastrogs)
      else
        write(1,'("# Reaction rate for ",a,"(",a1,",p)")')
     +    trim(targetnuclide),parsym(k0)
      endif
      write(1,'("#    T       Rate       MACS       G(T)")')
      do i=1,nTmax
        write(1,'(f8.4,3es12.5)') T9(i),rateastro(1,0,i),
     +    macsastro(1,0,i),partf(i)
      enddo
      close (unit=1)
c
c output partial rates(p,g) to given excited states in a specific file
c
      if (.not.flagastroex) goto 340
      do nex=1,Nlast(1,0,0)
        if (tau(1,0,nex).ne.0.) goto 330
      enddo
      goto 340
  330 do nex=0,Nlast(1,0,0)
        if (nex.eq.0.or.tau(1,0,nex).ne.0.) then
          astrofile='astrorate.p'
          write(astrofile(12:13),'(i2.2)') nex
          open (unit=1,file=astrofile,status='replace')
          if (nTmax.eq.1) then
            write(1,'("# Reaction rate for ",a,"(",a1,",p) at <E>=",
     +        f8.5," MeV (Exc. States Cont. : ",a1,
     +        ") to final level L",i2.2," at E=",f12.5," MeV")')
     +        trim(targetnuclide),parsym(k0),astroE,
     +        yesno(flagastrogs),nex,edis(1,0,nex)
          else
            write(1,'("# Reaction rate for Z=",i3," A=",i3,2x,a,"(",
     +        a1,",p) to the excited state L",i2.2," at E=",f12.5,
     +        " MeV")') Ztarget,Atarget,trim(targetnuclide),parsym(k0),
     +        nex,edis(1,0,nex)
          endif
          write(1,'("#   T       Rate         MACS      Branching")')
          do i=1,nTmax
            write(1,'(1x,f8.4,3es12.5)') T9(i),
     +        rateastroex(1,0,i,nex),
     +        macsastroex(1,0,i,nex),
     +        rateastroex(1,0,i,nex)/rateastro(1,0,i)
          enddo
          close (unit=1)
        endif
      enddo
  340 continue
      astrofile='astrorate.a'
      open (unit=1,file=astrofile,status='replace')
      if (nTmax.eq.1) then
        write(1,'("# Reaction rate for ",a,"(",a1,",a) at <E>=",
     +    f8.5," MeV (Excited States Contribution : ",a1,")")')
     +    trim(targetnuclide),parsym(k0),astroE,yesno(flagastrogs)
      else
        write(1,'("# Reaction rate for ",a,"(",a1,",a)")')
     +  trim(targetnuclide),parsym(k0)
      endif
      write(1,'("#    T       Rate       MACS       G(T)")')
      do i=1,nTmax
        write(1,'(f8.4,3es12.5)') T9(i),rateastro(2,2,i),
     +    macsastro(2,2,i),partf(i)
      enddo
      close (unit=1)
c
c output partial rates(p,g) to given excited states in a specific file
c
c ires   : counter
c iresprod   : counter
c iwriterp   : counter
c zrespro: Z of residual product
c arespro: A of residual product
c
      if (.not.flagastroex) goto 440
      do nex=1,Nlast(2,2,0)
        if (tau(2,2,nex).ne.0.) goto 430
      enddo
      goto 440
  430 do nex=0,Nlast(2,2,0)
        if (nex.eq.0.or.tau(2,2,nex).ne.0.) then
          astrofile='astrorate.a'
          write(astrofile(12:13),'(i2.2)') nex
          open (unit=1,file=astrofile,status='replace')
          if (nTmax.eq.1) then
            write(1,'("# Reaction rate for ",a,"(",a1,",a) at <E>=",
     +        f8.5," MeV (Exc. States Cont. : ",a1,
     +        ") to final level L",i2.2," at E=",f12.5," MeV")')
     +        trim(targetnuclide),parsym(k0),astroE,
     +        yesno(flagastrogs),nex,edis(2,2,nex)
          else
            write(1,'("# Reaction rate for Z=",i3," A=",i3,2x,a,"(",
     +        a1,",a) to the excited state L",i2.2," at E=",f12.5,
     +        " MeV")') Ztarget,Atarget,trim(targetnuclide),parsym(k0),
     +        nex,edis(2,2,nex)
          endif
          write(1,'("#   T       Rate         MACS      Branching")')
          do i=1,nTmax
            write(1,'(1x,f8.4,3es12.5)') T9(i),
     +        rateastroex(2,2,i,nex),
     +        macsastroex(2,2,i,nex),
     +        rateastroex(2,2,i,nex)/rateastro(2,2,i)
          enddo
          close (unit=1)
        endif
      enddo
  440 continue
      if (flagfission) then
        astrofile='astrorate.f'
        open (unit=1,file=astrofile,status='replace')
        if (nTmax.eq.1) then
          write(1,'("# Reaction rate for ",a,"(",a1,",f) at <E>=",
     +      f8.5," MeV (Excited States Contribution : ",a1,")")')
     +      trim(targetnuclide),parsym(k0),astroE,
     +      yesno(flagastrogs)
        else
          write(1,'("# Reaction rate for ",a,"(",a1,",a)")')
     +    trim(targetnuclide),parsym(k0)
        endif
        write(1,'("#    T       Rate       MACS       G(T)")')
        do i=1,nTmax
          write(1,'(f8.4,3es12.5)') T9(i),rateastrofis(i),
     +      macsastrofis(i),partf(i)
        enddo
        close (unit=1)
      endif
      rateastroeps=1.e-10
      iresprod=0
      do Acomp=0,maxAastro
        do Zcomp=0,maxZastro
          Ncomp=Acomp-Zcomp
          if (Ncomp.lt.0.or.Ncomp.gt.maxNastro) cycle
          if (Zcomp.eq.0.and.Ncomp.eq.0) cycle
          if (Zcomp.eq.0.and.Ncomp.eq.1) cycle
          if (Zcomp.eq.1.and.Ncomp.eq.0) cycle
          if (Zcomp.eq.2.and.Ncomp.eq.2) cycle
          iwriterp=0
          do i=1,nTmax
            if (rateastro(Zcomp,Ncomp,i).gt.rateastroeps) iwriterp=1
          enddo
          if (iwriterp.eq.1) then
            iresprod=iresprod+1
            do i=1,nTmax
              rateastrorp(i,iresprod)=rateastro(Zcomp,Ncomp,i)
            enddo
            zrespro(iresprod)=ZZ(Zcomp,Ncomp,0)
            arespro(iresprod)=AA(Zcomp,Ncomp,0)
          endif
        enddo
      enddo
      astrofile='astrorate.tot'
      targetparity='+'
      if (targetP.eq.-1) targetparity='-'
      open (unit=1,file=astrofile,status='replace')
      if (nTmax.eq.1) then
        write(1,'("# Reaction rate for ",2i4,a2,"+",a1,3x,i4,
     +    " reactions  Jp(GS)=",f5.1,a1," at <E>=",f8.5,
     +    " MeV (Excited States Contribution : ",a1,")")')
     +    Ztarget,Atarget,Starget,parsym(k0),iresprod+5,targetspin,
     +    targetparity,astroE,yesno(flagastrogs)
      else
        if (.not.flagastrogs) then
          if (nonthermlev.eq.-1) then
            write(1,'("# Reaction rate for ",2i4,a2,"+",a1,3x,i4,
     +        " reactions  Jp(GS)=",f5.1,a1,
     +        " Tot: fully thermalized target")')
     +        Ztarget,Atarget,Starget,parsym(k0),iresprod+5,targetspin,
     +        targetparity
          endif
          if (nonthermlev.gt.0) then
            write(1,'("# Reaction rate for ",2i4,a2,"+",a1,3x,i4,
     +        " reactions  Jp(GS)=",f5.1,a1,
     +        " GS: thermalized target on the GS",
     +        " excluding Level nb=",i2)')
     +        Ztarget,Atarget,Starget,parsym(k0),iresprod+5,targetspin,
     +        targetparity,nonthermlev
          endif
          if (nonthermlev.eq.0) then
            write(1,'("# Reaction rate for ",2i4,a2,"+",a1,3x,i4,
     +        " reactions  Jp(GS)=",f5.1,a1,
     +        " Isom: thermalized target in Level nb=",i2,
     +        " and excluding Level nb=",i2)')
     +        Ztarget,Atarget,Starget,parsym(k0),iresprod+5,targetspin,
     +        targetparity,Ltarget,nonthermlev
          endif
        else
          if (Ltarget.eq.0) then
            write(1,'("# Reaction rate for ",2i4,a2,"+",a1,3x,i4,
     +        " reactions  Jp(GS)=",f5.1,a1,
     +        " GS: non-thermalized target in its ground state")')
     +        Ztarget,Atarget,Starget,parsym(k0),iresprod+5,targetspin,
     +        targetparity
          endif
          if (Ltarget.gt.0) then
            write(1,'("# Reaction rate for ",2i4,a2,"+",a1,3x,i4,
     +        " reactions  Jp(GS)=",f5.1,a1,
     +        " Isom: non-thermalized target in its excited state = L",
     +        i2.2)')
     +        Ztarget,Atarget,Starget,parsym(k0),iresprod+5,targetspin,
     +        targetparity,Ltarget
          endif
        endif
      endif
      write(1,'("     T9      G(T)    (",a1,",g)",i3,a2,"  (",a1,
     +  ",n)",i3,a2,"  (",a1,",p)",i3,a2,"  (",a1,",a)",i3,a2,1x,
     +  "   Fission  ",
     +  200(6x,i3,a2,1x))')
     +  parsym(k0),AA(0,0,0),nuc(ZZ(0,0,0)),
     +  parsym(k0),AA(0,1,0),nuc(ZZ(0,1,0)),
     +  parsym(k0),AA(1,0,0),nuc(ZZ(1,0,0)),
     +  parsym(k0),AA(2,2,0),nuc(ZZ(2,2,0)),
     +  (arespro(ires),nuc(zrespro(ires)),ires=1,iresprod)
      do i=1,nTmax
        write(1,'(f8.4,200es12.5)') T9(i),partf(i),rateastro(0,0,i),
     +    rateastro(0,1,i),rateastro(1,0,i),rateastro(2,2,i),
     +    rateastrofis(i),(rateastrorp(i,ires),ires=1,iresprod)
      enddo
      close (unit=1)
      return
      end
Copyright (C)  2016 A.J. Koning, S. Hilaire and S. Goriely
