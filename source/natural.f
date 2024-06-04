      subroutine natural
c
c +---------------------------------------------------------------------
c | Author: Arjan Koning 
c | Date  : June 11, 2007
c | Task  : Calculation for natural elements
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      logical      lexist,elexist,specexist,resexist,xsexist,
     +             fissionexist,fyexist
      character*13 totfile,prodfile
      character*15 fisfile
      character*16 resfile,xsfile
      character*20 discfile,specfile,fyfile
      integer      i,k,k2,type,iang,lenfile,zbeg,zend,abeg,aend,iz,ia,
     +             npart,ih,it,id,ip,in,nex,nen
      real         en(numen2),xst(9),xstotnat(9,numen2),
     +             xsprodnat(numen2),xsyieldnat(numen2),xsnat(0:numen2),
     +             xs1nat(0:numen2),xs2nat(0:numen2),xs,y,xs1,xs2,
     +             xs3nat(numen2),xs4nat(numen2),xs3,xs4
c
c **************** Create runs and directories per isotope *************
c
c isonum       : number of isotopes
c iso          : counter for isotope
c talysinput   : subroutine for user input and defaults
c talysinitial : subroutine for initialization of nuclear structure
c talysreaction: subroutine with reaction models      
c
c For mono-isotopic nuclides we are done.
c
      if (isonum.eq.1) return
c
c Do a full TALYS calculation for each isotope
c
      do 10 iso=2,isonum
        call talysinput
        call talysinitial
        call talysreaction
   10 continue
c
c ******** Merge output files in results for natural elements **********
c
c 1. Total cross sections
c
c numen2    : maximum number of outgoing energies
c en        : incident energy
c xsprodnat : production cross sections for natural element
c xsyieldnat: yields for natural element
c xst       : help variable
c xstotnat  : total cross sections for natural element
c filetotal : flag for total cross sections on separate file
c totfile   : file with total cross sections
c numinc    : number of incident energies
c natstring : string extension for file names
c abun      : isotopic abundance
c parsym    : symbol of particle   
c nuc       : symbol of nucleus 
c
      do 110 k=1,numen2
        en(k)=0.
        xsprodnat(k)=0.
        xsyieldnat(k)=0.
        do 110 k2=1,9
          xst(k2)=0.
          xstotnat(k2,k)=0.
  110 continue
      if (filetotal) then
        do 120 i=1,isonum
          totfile='total.tot'//natstring(i)
          inquire (file=totfile,exist=lexist)
          if (lexist) then
            open (2,status='old',file=totfile)
            read(2,'(////)',end=150,err=150)
            do 130 k=1,numinc
              read(2,'(f10.3,2x,1p,9e11.4)',end=150,err=130) en(k),
     +          (xst(k2),k2=1,9)
              do 140 k2=1,9
                xstotnat(k2,k)=xstotnat(k2,k)+abun(i)*xst(k2)
  140         continue
  130       continue
  150       close (unit=2)
          endif
  120   continue
        open (3,status='unknown',file='total.tot')
        write(3,'("# ",a1," + nat-",a2," Total cross sections")')
     +    parsym(k0),nuc(Ztarget)
        write(3,'("# ")')
        write(3,'("# ")')
        write(3,'("# # energies =",i3)') numinc
        write(3,'("#    E      Non-elastic  Elastic     Total",
     +    "     Comp. el.  Shape el.  Reaction",
     +    " Comp. nonel   Direct   Pre-equil.")')              
        do 160 k=1,numinc
          write(3,'(f10.3,2x,1p,9e11.4)') en(k),(xstotnat(k2,k),k2=1,9)
  160   continue
        close (unit=3)
c
c 2. Particle production cross sections
c
c prodfile: file with total particle production cross sections
c parname : name of particle   
c
        do 210 type=1,6
          do 220 k=1,numinc
            xsprodnat(k)=0.
            xsyieldnat(k)=0.
  220     continue
          do 230 i=1,isonum
            prodfile=' prod.tot'//natstring(i)
            write(prodfile(1:1),'(a1)') parsym(type)
            inquire (file=prodfile,exist=lexist)
            if (lexist) then
              open (2,status='old',file=prodfile)
              read(2,'(////)',end=250,err=250)
              do 240 k=1,numinc
                read(2,'(f10.3,2e12.5)',end=250,err=240) en(k),xs,y
                xsprodnat(k)=xsprodnat(k)+abun(i)*xs
                xsyieldnat(k)=xsyieldnat(k)+abun(i)*y
  240         continue
  250         close (unit=2)
            endif
  230     continue
          open (3,status='unknown',file=prodfile(1:9))
          write(3,'("# ",a1," + nat-",a2," Total ",a8," production")')
     +      parsym(k0),nuc(Ztarget),parname(type)
          write(3,'("# ")')
          write(3,'("# ")')
          write(3,'("# # energies =",i3)') numinc
          write(3,'("#    E         xs         Yield")')       
          do 260 k=1,numinc
            write(3,'(1p,e10.3,2e12.5)') en(k),xsprodnat(k),
     +        xsyieldnat(k)
  260     continue
          close (unit=3)
  210   continue
      endif
c
c 3. Elastic scattering angular distribution
c
c nangle   : number of angles 
c xsnat,...: cross section for natural element
c elexist  : logical to determine existence of elastic scattering file
c discfile : file with elastic scattering angular distribution
c angle,ang: angle
c eninc    : incident energy in MeV
c
      do 310 k=1,numinc
        do 320 iang=0,nangle
          xsnat(iang)=0.
          xs1nat(iang)=0.
          xs2nat(iang)=0.
  320   continue
        elexist=.false.
        do 330 i=1,isonum
          discfile='nn       ang.L00'//natstring(i)
          write(discfile(3:9),'(f7.3)') eninc(k)
          write(discfile(3:5),'(i3.3)') int(eninc(k))
          inquire (file=discfile,exist=lexist)
          if (lexist) then
            elexist=.true.
            open (2,status='old',file=discfile)
            read(2,'(////)',end=350,err=350)
            do 340 iang=0,nangle
              read(2,'(f5.1,3e16.5)',end=350,err=340) angle(iang),xs,
     +          xs1,xs2
              xsnat(iang)=xsnat(iang)+abun(i)*xs
              xs1nat(iang)=xs1nat(iang)+abun(i)*xs1
              xs2nat(iang)=xs2nat(iang)+abun(i)*xs2
  340       continue
  350       close (unit=2)
          endif
  330   continue
        if (elexist) then
          open (3,status='unknown',file=discfile(1:16))
          write(3,'("# ",a1," + nat-",a2," Elastic scattering", 
     +      " angular distribution")') parsym(k0),nuc(Ztarget)
          write(3,'("# E-incident = ",f7.3)') eninc(k)
          write(3,'("# ")')
          write(3,'("# # angles   =",i3)') nangle+1
          write(3,'("#   E         xs            Direct",
     +      "         Compound")')     
          do 360 iang=0,nangle
            write(3,'(f5.1,1p,3e16.5)') angle(iang),xsnat(iang),
     +        xs1nat(iang),xs2nat(iang)
  360     continue
          close (unit=3)
        endif
  310 continue
c
c 4. Composite particle spectra
c
c specexist: logical to determine existence of spectrum file
c specfile : file with composite particle spectra
c lenfile  : length of file
c xs,xs1,..: help variables
c ebegin   : first energy point of energy grid  
c eendout  : last energy point of energy grid 
c
      do 410 k=1,numinc
        do 420 type=1,6
          do 430 k2=1,numen2
            xsnat(k2)=0.
            xs1nat(k2)=0.
            xs2nat(k2)=0.
            xs3nat(k2)=0.
            xs4nat(k2)=0.
  430     continue
          specexist=.false.
          do 440 i=1,isonum
            specfile=parsym(type)//'spec000.000.tot'//natstring(i)
            write(specfile(6:12),'(f7.3)') eninc(k)
            write(specfile(6:8),'(i3.3)') int(eninc(k))
            inquire (file=specfile,exist=lexist)
            if (lexist) then
              specexist=.true.
              open (2,status='old',file=specfile)
              read(2,'(////)',end=460,err=460)
              do 450 k2=1,numen2
                lenfile=k2-1
                read(2,'(f7.3,5e12.5)',end=460,err=450) en(k2),xs,xs1,
     +            xs2,xs3,xs4
                xsnat(k2)=xsnat(k2)+abun(i)*xs
                xs1nat(k2)=xs1nat(k2)+abun(i)*xs1
                xs2nat(k2)=xs2nat(k2)+abun(i)*xs2
                xs3nat(k2)=xs3nat(k2)+abun(i)*xs3
                xs4nat(k2)=xs4nat(k2)+abun(i)*xs4
  450         continue
  460         close (unit=2)
            endif
  440     continue
          if (specexist) then
            open (3,status='unknown',file=specfile(1:16))
            write(3,'("# ",a1," + nat-",a2,": ",a8," spectrum")')
     +        parsym(k0),nuc(Ztarget),parname(type)
            write(3,'("# E-incident = ",f7.3)') eninc(k)
            write(3,'("# ")')
            write(3,'("# # energies =",i3)') 
     +        eendout(type)-ebegin(type)+1
            write(3,'("# E-out    Total       Direct    Pre-equil.",
     +        "  Mult. preeq  Compound")')             
            do 470 k2=1,lenfile
              write(3,'(f7.3,1p,5e12.5)') en(k2),xsnat(k2),xs1nat(k2),
     +          xs2nat(k2),xs3nat(k2),xs4nat(k2)
  470       continue
            close (unit=3)
          endif
  420   continue
  410 continue
c
c 5. Residual production cross sections
c
c zbeg,..: help variables
c Ztarget: charge number of target nucleus
c isotope: isotope number of residual nucleus
c numZ   : maximal number of protons away from the initial compound 
c          nucleus
c numN   : maximal number of neutrons away from the initial compound 
c          nucleus
c
      zbeg=Ztarget-numZ-2
      zend=Ztarget+2
      abeg=isotope(1)-numN-2
      aend=isotope(isonum)+4
      do 510 iz=zbeg,zend
        do 510 ia=abeg,aend
c
c Total
c
c resexist: logical to determine existence of residual production file
c resfile : file with residual production cross sections
c
          do 520 k=1,numinc
            xsnat(k)=0.
  520     continue
          resexist=.false.
          do 530 i=1,isonum
            resfile='rp000000.tot'//natstring(i)
            write(resfile(3:8),'(2i3.3)') iz,ia
            inquire (file=resfile,exist=lexist)
            if (lexist) then
              resexist=.true.
              open (2,status='old',file=resfile)
              read(2,'(////)',end=550,err=550)
              do 540 k=1,numinc
                read(2,'(f10.3,e12.5)',end=550,err=540) en(k),xs
                xsnat(k)=xsnat(k)+abun(i)*xs
  540         continue
  550         close (unit=2)
            endif
  530     continue
          if (resexist) then
            open (3,status='unknown',file=resfile(1:12))
            write(3,'("# ",a1," + nat-",a2,": Production of ",i3,a2,
     +        " - Total")') parsym(k0),nuc(Ztarget),ia,nuc(iz)
            write(3,'("# ")')
            write(3,'("# ")')
            write(3,'("# # energies =",i3)') numinc
            write(3,'("#    E         xs")')                
            do 560 k=1,numinc
              write(3,'(1p,e10.3,e12.5)') en(k),xsnat(k)
  560       continue
            close (unit=3)
          endif
c
c Per ground state and isomer
c
          write(resfile(10:12),'("L00")')
          inquire (file=resfile,exist=lexist)
          if (.not.lexist) goto 510
          do 570 nex=0,numlev
            do 580 k=1,numinc
              xsnat(k)=0.
  580       continue
            resexist=.false.
            do 590 i=1,isonum
              resfile='rp000000.L00'//natstring(i)
              write(resfile(3:8),'(2i3.3)') iz,ia
              write(resfile(11:12),'(i2.2)') nex
              inquire (file=resfile,exist=lexist)
              if (lexist) then
                resexist=.true.
                open (2,status='old',file=resfile)
                read(2,'(////)',end=610,err=610)
                do 600 k=1,numinc
                  read(2,'(f10.3,e12.5)',end=610,err=600) en(k),xs
                  xsnat(k)=xsnat(k)+abun(i)*xs
  600           continue
  610           close (unit=2)
              endif
  590       continue
            if (resexist) then
              open (3,status='unknown',file=resfile(1:12))
              write(3,'("# ",a1," + nat-",a2,": Production of ",i3,a2,
     +          " - Level",i3)') parsym(k0),nuc(Ztarget),ia,nuc(iz),nex
              write(3,'("# ")')
              write(3,'("# ")')
              write(3,'("# # energies =",i3)') numinc
              write(3,'("#    E         xs")')                
              do 620 k=1,numinc
                write(3,'(1p,e10.3,e12.5)') en(k),xsnat(k)
  620         continue
              close (unit=3)
            endif
  570     continue
  510 continue
c
c 6. Exclusive channel cross sections
c
c npart     : number of particles in outgoing channel
c maxchannel: maximal number of outgoing particles in individual
c             channel description (e.g. this is 3 for (n,2np))     
c numin,....: maximal number of ejectile in channel description  
c xsfile    : file with channel cross sections
c
      do 710 npart=0,maxchannel
      do 710 ia=0,numia
      do 710 ih=0,numih
      do 710 it=0,numit
      do 710 id=0,numid
      do 710 ip=0,numip
      do 710 in=0,numin    
        if (in+ip+id+it+ih+ia.ne.npart) goto 710
        do 720 k=1,numinc
          xsnat(k)=0.
  720   continue
        xsexist=.false.
        do 730 i=1,isonum
          xsfile='xs000000.tot'//natstring(i)
          write(xsfile(3:8),'(6i1)') in,ip,id,it,ih,ia
          inquire (file=xsfile,exist=lexist)
          if (lexist) then
            xsexist=.true.
            open (2,status='old',file=xsfile)
            read(2,'(////)',end=750,err=750)
            do 740 k=1,numinc
              read(2,'(f10.3,e12.5)',end=750,err=740) en(k),xs
              xsnat(k)=xsnat(k)+abun(i)*xs
  740       continue
  750       close (unit=2)
          endif
  730   continue
        if (xsexist) then
          open (3,status='unknown',file=xsfile(1:12))
          write(3,'("# ",a1," + nat-",a2)') parsym(k0),nuc(Ztarget)
          write(3,'("# ")') 
          write(3,'("# ")') 
          write(3,'("# # energies =",i3)') numinc
          write(3,'("#    E         xs")')
          do 760 k=1,numinc
            write(3,'(1p,e10.3,e12.5)') en(k),xsnat(k)
  760     continue
          close (unit=3)
        endif
  710 continue
c
c 7. Fission cross sections
c
c fissionexist: logical to determine existence of fission file
c fisfile     : file with fission cross sections
c
      do 810 k=1,numinc
        xsnat(k)=0.
  810 continue
      fissionexist=.false.
      do 820 i=1,isonum
        fisfile='fission.tot'//natstring(i)
        inquire (file=fisfile,exist=lexist)
        if (lexist) then
          fissionexist=.true.
          open (2,status='old',file=fisfile)
          read(2,'(////)',end=840,err=840)
          do 830 k=1,numinc
            read(2,'(f10.3,e12.5)',end=840,err=830) en(k),xs
            xsnat(k)=xsnat(k)+abun(i)*xs
  830     continue
  840     close (unit=2)
        endif
  820 continue
      if (fissionexist) then
        open (3,status='unknown',file='fission.tot')
        write(3,'("# ",a1," + nat-",a2," Total fission cross ",
     +    "section")') parsym(k0),nuc(Ztarget)
        write(3,'("# ")')
        write(3,'("# ")')
        write(3,'("# # energies =",i3)') numinc
        write(3,'("#    E         xs")')
        do 850 k=1,numinc
          write(3,'(f7.3,1p,e12.5)') en(k),xsnat(k)
  850   continue
        close (unit=3)
      endif
c
c 8. Fission yields                   
c
c Atarget   : mass number of target nucleus
c flagffevap: flag for calculation of particle evaporation from
c             fission fragment mass yields    
c
c Excitation function per fission product
c
      do 910 iz=1,Ztarget   
        do 910 ia=1,Atarget   
          do 920 k=1,numinc
            xs1nat(k)=0.
            xs2nat(k)=0.
  920     continue
          resexist=.false.
          do 930 i=1,isonum
            resfile='fp000000.tot'//natstring(i)
            write(resfile(3:8),'(2i3.3)') iz,ia
            inquire (file=resfile,exist=lexist)
            if (lexist) then
              resexist=.true.
              open (2,status='old',file=resfile)
              read(2,'(////)',end=950,err=950)
              do 940 k=1,numinc
                read(2,'(e10.3,e12.4,3x,e12.4)',end=950,err=940) 
     +            en(k),xs1,xs2
                xs1nat(k)=xs1nat(k)+abun(i)*xs1
                xs2nat(k)=xs2nat(k)+abun(i)*xs2
  940         continue
  950         close (unit=2)
            endif
  930     continue
          if (resexist) then
            open (3,status='unknown',file=resfile(1:12))
            write(3,'("# ",a1," + nat-",a2,": ff yield of ",i3,a2)')
     +        parsym(k0),nuc(Ztarget),ia,nuc(iz)
            write(3,'("# ")')
            write(3,'("# ")')
            write(3,'("# # energies =",i3)') numinc
            write(3,'("#    E         xs")')                
            if (flagffevap) then
              write(3,'("# E-incident   FF Yield   FP yield")')
              do 960 nen=1,numinc
                write(3,'(1p,e10.3,e12.4,3x,e12.4)') eninc(nen),
     +            xs1nat(nen),xs2nat(nen)
  960         continue
            else
              write(3,'("# E-incident   FF Yield")')
              do 970 nen=1,numinc
                write(3,'(1p,e10.3,e12.4)') eninc(nen),xs1nat(nen)
  970         continue
            endif
            close (unit=3)
          endif
  910 continue
c
c Mass distribution per incident energy  
c
c fyfile: file with fission yields
c
      do 1010 k=1,numinc
        do 1020 ia=1,Atarget
          xs1nat(ia)=0.
          xs2nat(ia)=0.
 1020   continue
        fyexist=.false.
        do 1030 i=1,isonum
          fyfile='yield000.000.fis'//natstring(i)
          write(fyfile(6:12),'(f7.3)') eninc(k)
          write(fyfile(6:8),'(i3.3)') int(eninc(k))
          inquire (file=fyfile,exist=lexist)
          if (lexist) then
            fyexist=.true.
            open (2,status='old',file=fyfile)
            read(2,'(////)',end=1050,err=1050)
            do 1040 ia=1,isotope(i)
              read(2,'(3x,2e15.4)',end=1050,err=1040) xs1,xs2
              xs1nat(ia)=xs1nat(ia)+abun(i)*xs1
              xs2nat(ia)=xs2nat(ia)+abun(i)*xs2
 1040       continue
 1050       close (unit=2)
          endif
 1030   continue
        if (fyexist) then
          open (3,status='unknown',file=fyfile(1:16))
          write(3,'("# ",a1," +  nat-",a2,": mass yields")')
     +      parsym(k0),nuc(Ztarget)
          write(3,'("# E-incident = ",f7.3)') eninc(k)
          write(3,'("# ")')
          write(3,'("# ")')
          write(3,'("# Mass    Yield   Corrected yield")')   
          do 1060 ia=1,Atarget
            write(3,'(i3,3x,1p,e12.4,3x,e12.4)') ia,xs1nat(ia),
     +        xs2nat(ia)
 1060     continue         
          close (unit=3)
        endif
 1010 continue
      return
      end
Copyright (C) 2004  A.J. Koning, S. Hilaire and M.C. Duijvestijn
