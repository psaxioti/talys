      subroutine natural
c
c +---------------------------------------------------------------------
c | Author: Arjan Koning 
c | Date  : December 6, 2006
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
      real         en(numen2),xst(9),xstotnat(9,0:numen2),
     +             xsprodnat(0:numen2),xsyieldnat(0:numen2),
     +             xsnat(0:numen2),xs1nat(0:numen2),xs2nat(0:numen2),xs,
     +             y,xs1,xs2,xs3nat(0:numen2),xs4nat(0:numen2),xs3,xs4
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
c numinc    : number of incident energies
c xsprodnat : production cross sections for natural element
c xsyieldnat: yields for natural element
c xstotnat  : total cross sections for natural element
c filetotal : flag for total cross sections on separate file
c totfile   : file with total cross sections
c natstring : string extension for file names
c en        : incident energy
c xst       : help variable
c abun      : isotopic abundance
c parsym    : symbol of particle   
c nuc       : symbol of nucleus 
c
      do 110 k=1,numinc
        xsprodnat(k)=0.
        xsyieldnat(k)=0.
        do 110 k2=1,9
          xstotnat(k2,k)=0.
  110 continue
      if (filetotal) then
        do 120 i=1,isonum
          totfile='total.tot'//natstring(i)
          inquire (file=totfile,exist=lexist)
          if (lexist) then
            open (2,status='old',file=totfile)
            read(2,'(////)')
            do 130 k=1,numinc
              read(2,'(f10.3,2x,1p,9e11.4)',err=130) en(k),
     +          (xst(k2),k2=1,9)
              do 140 k2=1,9
                xstotnat(k2,k)=xstotnat(k2,k)+abun(i)*xst(k2)
  140         continue
  130       continue
            close (unit=2)
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
        do 150 k=1,numinc
          write(3,'(f10.3,2x,1p,9e11.4)') en(k),(xstotnat(k2,k),k2=1,9)
  150   continue
        close (unit=3)
c
c 2. Particle production cross sections
c
c prodfile: file with total particle production cross sections
c parname : name of particle   
c
        do 160 type=1,6
          do 170 k=1,numinc
            xsprodnat(k)=0.
            xsyieldnat(k)=0.
  170     continue
          do 180 i=1,isonum
            prodfile=' prod.tot'//natstring(i)
            write(prodfile(1:1),'(a1)') parsym(type)
            inquire (file=prodfile,exist=lexist)
            if (lexist) then
              open (2,status='old',file=prodfile)
              read(2,'(////)')
              do 190 k=1,numinc
                read(2,'(f10.3,2e12.5)',err=190) en(k),xs,y
                xsprodnat(k)=xsprodnat(k)+abun(i)*xs
                xsyieldnat(k)=xsyieldnat(k)+abun(i)*y
  190         continue
              close (unit=2)
            endif
  180     continue
          open (3,status='unknown',file=prodfile(1:9))
          write(3,'("# ",a1," + nat-",a2," Total ",a8," production")')
     +      parsym(k0),nuc(Ztarget),parname(type)
          write(3,'("# ")')
          write(3,'("# ")')
          write(3,'("# # energies =",i3)') numinc
          write(3,'("#    E         xs         Yield")')       
          do 200 k=1,numinc
            write(3,'(1p,e10.3,2e12.5)') en(k),xsprodnat(k),
     +        xsyieldnat(k)
  200     continue
          close (unit=3)
  160   continue
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
      do 210 k=1,numinc
        do 220 iang=0,nangle
          xsnat(iang)=0.
          xs1nat(iang)=0.
          xs2nat(iang)=0.
  220   continue
        elexist=.false.
        do 230 i=1,isonum
          discfile='nn       ang.L00'//natstring(i)
          write(discfile(3:9),'(f7.3)') eninc(k)
          write(discfile(3:5),'(i3.3)') int(eninc(k))
          inquire (file=discfile,exist=lexist)
          if (lexist) then
            elexist=.true.
            open (2,status='old',file=discfile)
            read(2,'(////)')
            do 240 iang=0,nangle
              read(2,'(f5.1,3e16.5)',err=240) angle(iang),xs,xs1,xs2
              xsnat(iang)=xsnat(iang)+abun(i)*xs
              xs1nat(iang)=xs1nat(iang)+abun(i)*xs1
              xs2nat(iang)=xs2nat(iang)+abun(i)*xs2
  240       continue
            close (unit=2)
          endif
  230   continue
        if (elexist) then
          open (3,status='unknown',file=discfile(1:16))
          write(3,'("# ",a1," + nat-",a2," Elastic scattering", 
     +      " angular distribution")') parsym(k0),nuc(Ztarget)
          write(3,'("# E-incident = ",f7.3)') eninc(k)
          write(3,'("# ")')
          write(3,'("# # angles   =",i3)') nangle+1
          write(3,'("#   E         xs            Direct",
     +      "         Compound")')     
          do 250 iang=0,nangle
            write(3,'(f5.1,1p,3e16.5)') angle(iang),xsnat(iang),
     +        xs1nat(iang),xs2nat(iang)
  250     continue
          close (unit=3)
        endif
  210 continue
c
c 4. Composite particle spectra
c
c numen2   : maximum number of outgoing energies
c specexist: logical to determine existence of spectrum file
c specfile : file with composite particle spectra
c lenfile  : length of file
c xs,xs1,..: help variables
c ebegin   : first energy point of energy grid  
c eendout  : last energy point of energy grid 
c
      do 310 k=1,numinc
        do 320 type=1,6
          do 330 k2=1,numen2
            xsnat(k2)=0.
            xs1nat(k2)=0.
            xs2nat(k2)=0.
            xs3nat(k2)=0.
            xs4nat(k2)=0.
  330     continue
          specexist=.false.
          do 340 i=1,isonum
            specfile=parsym(type)//'spec000.000.tot'//natstring(i)
            write(specfile(6:12),'(f7.3)') eninc(k)
            write(specfile(6:8),'(i3.3)') int(eninc(k))
            inquire (file=specfile,exist=lexist)
            if (lexist) then
              specexist=.true.
              open (2,status='old',file=specfile)
              read(2,'(////)')
              do 350 k2=1,numen2
                lenfile=k2-1
                read(2,'(f7.3,5e12.5)',err=350,end=360) en(k2),xs,xs1,
     +            xs2,xs3,xs4
                xsnat(k2)=xsnat(k2)+abun(i)*xs
                xs1nat(k2)=xs1nat(k2)+abun(i)*xs1
                xs2nat(k2)=xs2nat(k2)+abun(i)*xs2
                xs3nat(k2)=xs3nat(k2)+abun(i)*xs3
                xs4nat(k2)=xs4nat(k2)+abun(i)*xs4
  350         continue
  360         close (unit=2)
            endif
  340     continue
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
            do 370 k2=1,lenfile
              write(3,'(f7.3,1p,5e12.5)') en(k2),xsnat(k2),xs1nat(k2),
     +          xs2nat(k2),xs3nat(k2),xs4nat(k2)
  370       continue
            close (unit=3)
          endif
  320   continue
  310 continue
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
      do 410 iz=zbeg,zend
        do 410 ia=abeg,aend
c
c Total
c
c resexist: logical to determine existence of residual production file
c resfile : file with residual production cross sections
c
          do 420 k=1,numinc
            xsnat(k)=0.
  420     continue
          resexist=.false.
          do 430 i=1,isonum
            resfile='rp000000.tot'//natstring(i)
            write(resfile(3:8),'(2i3.3)') iz,ia
            inquire (file=resfile,exist=lexist)
            if (lexist) then
              resexist=.true.
              open (2,status='old',file=resfile)
              read(2,'(////)')
              do 440 k=1,numinc
                read(2,'(f10.3,e12.5)',err=440) en(k),xs
                xsnat(k)=xsnat(k)+abun(i)*xs
  440         continue
              close (unit=2)
            endif
  430     continue
          if (resexist) then
            open (3,status='unknown',file=resfile(1:12))
            write(3,'("# ",a1," + nat-",a2,": Production of ",i3,a2,
     +        " - Total")') parsym(k0),nuc(Ztarget),ia,nuc(iz)
            write(3,'("# ")')
            write(3,'("# ")')
            write(3,'("# # energies =",i3)') numinc
            write(3,'("#    E         xs")')                
            do 450 k=1,numinc
              write(3,'(1p,e10.3,e12.5)') en(k),xsnat(k)
  450       continue
            close (unit=3)
          endif
c
c Per ground state and isomer
c
          write(resfile(10:12),'("L00")')
          inquire (file=resfile,exist=lexist)
          if (.not.lexist) goto 410
          do 460 nex=0,numlev
            do 470 k=1,numinc
              xsnat(k)=0.
  470       continue
            resexist=.false.
            do 480 i=1,isonum
              resfile='rp000000.L00'//natstring(i)
              write(resfile(3:8),'(2i3.3)') iz,ia
              write(resfile(11:12),'(i2.2)') nex
              inquire (file=resfile,exist=lexist)
              if (lexist) then
                resexist=.true.
                open (2,status='old',file=resfile)
                read(2,'(////)')
                do 490 k=1,numinc
                  read(2,'(f10.3,e12.5)',err=490) en(k),xs
                  xsnat(k)=xsnat(k)+abun(i)*xs
  490           continue
                close (unit=2)
              endif
  480       continue
            if (resexist) then
              open (3,status='unknown',file=resfile(1:12))
              write(3,'("# ",a1," + nat-",a2,": Production of ",i3,a2,
     +          " - Level",i3)') parsym(k0),nuc(Ztarget),ia,nuc(iz),nex
              write(3,'("# ")')
              write(3,'("# ")')
              write(3,'("# # energies =",i3)') numinc
              write(3,'("#    E         xs")')                
              do 500 k=1,numinc
                write(3,'(1p,e10.3,e12.5)') en(k),xsnat(k)
  500         continue
              close (unit=3)
            endif
  460     continue
  410 continue
c
c 6. Exclusive channel cross sections
c
c npart     : number of particles in outgoing channel
c maxchannel: maximal number of outgoing particles in individual
c             channel description (e.g. this is 3 for (n,2np))     
c numin,....: maximal number of ejectile in channel description  
c xsfile    : file with channel cross sections
c
      do 510 npart=0,maxchannel
      do 510 ia=0,numia
      do 510 ih=0,numih
      do 510 it=0,numit
      do 510 id=0,numid
      do 510 ip=0,numip
      do 510 in=0,numin    
        if (in+ip+id+it+ih+ia.ne.npart) goto 510
        do 520 k=1,numinc
          xsnat(k)=0.
  520   continue
        xsexist=.false.
        do 530 i=1,isonum
          xsfile='xs000000.tot'//natstring(i)
          write(xsfile(3:8),'(6i1)') in,ip,id,it,ih,ia
          inquire (file=xsfile,exist=lexist)
          if (lexist) then
            xsexist=.true.
            open (2,status='old',file=xsfile)
            read(2,'(////)')
            do 540 k=1,numinc
              read(2,'(f10.3,e12.5)',err=540) en(k),xs
              xsnat(k)=xsnat(k)+abun(i)*xs
  540       continue
            close (unit=2)
          endif
  530   continue
        if (xsexist) then
          open (3,status='unknown',file=xsfile(1:12))
          write(3,'("# ",a1," + nat-",a2)') parsym(k0),nuc(Ztarget)
          write(3,'("# ")') 
          write(3,'("# ")') 
          write(3,'("# # energies =",i3)') numinc
          write(3,'("#    E         xs")')
          do 550 k=1,numinc
            write(3,'(1p,e10.3,e12.5)') en(k),xsnat(k)
  550     continue
          close (unit=3)
        endif
  510 continue
c
c 7. Fission cross sections
c
c fissionexist: logical to determine existence of fission file
c fisfile     : file with fission cross sections
c
      do 610 k=1,numinc
        xsnat(k)=0.
  610 continue
      fissionexist=.false.
      do 620 i=1,isonum
        fisfile='fission.tot'//natstring(i)
        inquire (file=fisfile,exist=lexist)
        if (lexist) then
          fissionexist=.true.
          open (2,status='old',file=fisfile)
          read(2,'(////)')
          do 630 k=1,numinc
            read(2,'(f10.3,e12.5)',err=630) en(k),xs
            xsnat(k)=xsnat(k)+abun(i)*xs
  630     continue
          close (unit=2)
        endif
  620 continue
      if (fissionexist) then
        open (3,status='unknown',file='fission.tot')
        write(3,'("# ",a1," + nat-",a2," Total fission cross ",
     +    "section")') parsym(k0),nuc(Ztarget)
        write(3,'("# ")')
        write(3,'("# ")')
        write(3,'("# # energies =",i3)') numinc
        write(3,'("#    E         xs")')
        do 640 k=1,numinc
          write(3,'(f7.3,1p,e12.5)') en(k),xsnat(k)
  640   continue
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
      do 710 iz=1,Ztarget   
        do 710 ia=1,Atarget   
          do 720 k=1,numinc
            xs1nat(k)=0.
            xs2nat(k)=0.
  720     continue
          resexist=.false.
          do 730 i=1,isonum
            resfile='fp000000.tot'//natstring(i)
            write(resfile(3:8),'(2i3.3)') iz,ia
            inquire (file=resfile,exist=lexist)
            if (lexist) then
              resexist=.true.
              open (2,status='old',file=resfile)
              read(2,'(////)')
              do 740 k=1,numinc
                read(2,'(e10.3,e12.4,3x,e12.4)',err=740) en(k),xs1,xs2
                xs1nat(k)=xs1nat(k)+abun(i)*xs1
                xs2nat(k)=xs2nat(k)+abun(i)*xs2
  740         continue
              close (unit=2)
            endif
  730     continue
          if (resexist) then
            open (3,status='unknown',file=resfile(1:12))
            write(3,'("# ",a1," + nat-",a2,": ff yield of ",i3,a2)')
     +        ia,nuc(iz),parsym(k0),nuc(Ztarget)
            write(3,'("# ")')
            write(3,'("# ")')
            write(3,'("# # energies =",i3)') numinc
            write(3,'("#    E         xs")')                
            if (flagffevap) then
              write(3,'("# E-incident   FF Yield   FP yield")')
              do 750 nen=1,numinc
                write(3,'(1p,e10.3,e12.4,3x,e12.4)') eninc(nen),
     +            xs1nat(nen),xs2nat(nen)
  750         continue
            else
              write(3,'("# E-incident   FF Yield")')
              do 760 nen=1,numinc
                write(3,'(1p,e10.3,e12.4)') eninc(nen),xs1nat(nen)
  760         continue
            endif
            close (unit=3)
          endif
  710 continue
c
c Mass distribution per incident energy  
c
c fyfile: file with fission yields
c
      do 810 k=1,numinc
        do 820 ia=1,Atarget
          xs1nat(ia)=0.
          xs2nat(ia)=0.
  820   continue
        fyexist=.false.
        do 830 i=1,isonum
          fyfile='yield000.000.fis'//natstring(i)
          write(fyfile(6:12),'(f7.3)') eninc(k)
          write(fyfile(6:8),'(i3.3)') int(eninc(k))
          inquire (file=fyfile,exist=lexist)
          if (lexist) then
            fyexist=.true.
            open (2,status='old',file=fyfile)
            read(2,'(////)')
            do 840 ia=1,isotope(i)
              read(2,'(3x,2e15.4)',err=840) xs1,xs2
              xs1nat(ia)=xs1nat(ia)+abun(i)*xs1
              xs2nat(ia)=xs2nat(ia)+abun(i)*xs2
  840       continue
            close (unit=2)
          endif
  830   continue
        if (fyexist) then
          open (3,status='unknown',file=fyfile(1:16))
          write(3,'("# ",a1," +  nat-",a2,": mass yields")')
     +      parsym(k0),nuc(Ztarget)
          write(3,'("# E-incident = ",f7.3)') eninc(k)
          write(3,'("# ")')
          write(3,'("# ")')
          write(3,'("# Mass    Yield   Corrected yield")')   
          do 850 ia=1,Atarget
            write(3,'(i3,3x,1p,e12.4,3x,e12.4)') ia,xs1nat(ia),
     +        xs2nat(ia)
  850     continue         
          close (unit=3)
        endif
  810 continue
      return
      end
Copyright (C) 2004  A.J. Koning, S. Hilaire and M.C. Duijvestijn
