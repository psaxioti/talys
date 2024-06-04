      subroutine angleout
c
c +---------------------------------------------------------------------
c | Author: Arjan Koning 
c | Date  : October 16, 2004
c | Task  : Output of discrete angular distributions
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      character*16 discfile
      integer      LL,iang,Zix,Nix,i,type
c
c **************** Elastic scattering angular distribution *************
c
c 1. Legendre coefficients
c
c flaglegendre: flag for output of Legendre coefficients
c J2end       : 2 * end of J summation
c tleg        : total Legendre coefficient
c k0          : index of incident particle
c Ltarget     : excited level of target
c dleg        : direct reaction Legendre coefficient
c cleg        : compound nucleus Legendre coefficient
c tlegnor     : total Legendre coefficient normalized to 1
c fileelastic : flag for elastic angular distribution on separate file
c parsym      : symbol of particle
c Atarget     : mass number of target nucleus
c nuc         : symbol of nucleus
c Ztarget     : charge number of target nucleus  
c Einc        : incident energy in MeV 
c
      write(*,'(/"8. Discrete state angular distributions")')
      if (flaglegendre) then
        write(*,'(/"8a1. Legendre coefficients for elastic",$)')
        write(*,'(" scattering"/)')
        write(*,'("  L       Total           Direct",$)')
        write(*,'("         Compound       Normalized"/)')
        do 10 LL=0,J2end
          write(*,'(i3,1p,4e16.5)') LL,tleg(k0,Ltarget,LL),
     +      dleg(k0,Ltarget,LL),cleg(k0,Ltarget,LL),
     +      tlegnor(k0,Ltarget,LL)
   10   continue
c
c Write results to separate file
c
        if (fileelastic) then
          discfile='         leg.L00'
          write(discfile(1:2),'(2a1)') parsym(k0),parsym(k0)
          write(discfile(3:9),'(f7.3)') Einc
          write(discfile(3:5),'(i3.3)') int(Einc)
          open (unit=1,status='unknown',file=discfile)  
          write(1,'("# ",a1," + ",i3,a2,$)') parsym(k0),Atarget,
     +      nuc(Ztarget)
          write(1,'(" Elastic scattering Legendre coefficients")')
          write(1,'("# E-incident = ",f7.3)') Einc
          write(1,'("# ")')
          write(1,'("# # coeff.   =",i3)') J2end+1
          write(1,'("#  L       Total           Direct",$)')
          write(1,'("        Compound       Normalized")')
          do 20 LL=0,J2end
            write(1,'(i3,1p,4e16.5)') LL,tleg(k0,Ltarget,LL),
     +        dleg(k0,Ltarget,LL),cleg(k0,Ltarget,LL),
     +        tlegnor(k0,Ltarget,LL)
   20     continue
          close (unit=1)
        endif    
      endif
c
c 2. Angular distributions
c
c nangle  : number of angles
c angle   : angle
c discad  : discrete state angular distribution
c directad: direct angular distribution
c compad  : compound angular distribution
c ruth    : elastic/Rutherford ratio
c
      write(*,'(/"8a2. Elastic scattering angular distribution"/)')
      if (k0.eq.1) then
        write(*,'("Angle        Total          Direct",$)')
        write(*,'("         Compound"/)')
        do 30 iang=0,nangle
          write(*,'(f5.1,1p,3e16.5)') angle(iang),
     +      discad(k0,Ltarget,iang),directad(k0,Ltarget,iang),
     +      compad(k0,Ltarget,iang)
   30   continue
c
c Write results to separate file
c
        if (fileelastic) then
          discfile='nn       ang.L00'
          write(discfile(3:9),'(f7.3)') Einc
          write(discfile(3:5),'(i3.3)') int(Einc)
          open (unit=1,status='unknown',file=discfile)  
          write(1,'("# ",a1," + ",i3,a2,$)') parsym(k0),Atarget,
     +      nuc(Ztarget)
          write(1,'(" Elastic scattering angular distribution")')
          write(1,'("# E-incident = ",f7.3)') Einc
          write(1,'("# ")')
          write(1,'("# # angles   =",i3)') nangle+1
          write(1,'("#   E         xs            Direct",$)')
          write(1,'("         Compound")')     
          do 40 iang=0,nangle
            write(1,'(f5.1,1p,3e16.5)') angle(iang),
     +        discad(k0,Ltarget,iang),directad(k0,Ltarget,iang),
     +        compad(k0,Ltarget,iang)
   40     continue
          close (unit=1)
        endif    
      else
        write(*,'("Angle        Total            Direct",$)')
        write(*,'("       Compound       c.s/Rutherford"/)')
        do 50 iang=0,nangle
          write(*,'(f5.1,1p,4e16.5)') angle(iang),
     +      discad(k0,Ltarget,iang),directad(k0,Ltarget,iang),
     +      compad(k0,Ltarget,iang),ruth(iang)
   50   continue
c
c Write results to separate file
c
        if (fileelastic) then
          discfile='         ang.L00'
          write(discfile(1:2),'(2a1)') parsym(k0),parsym(k0)
          write(discfile(3:9),'(f7.3)') Einc
          write(discfile(3:5),'(i3.3)') int(Einc)
          open (unit=1,status='unknown',file=discfile)  
          write(1,'("# ",a1," + ",i3,a2,$)') parsym(k0),Atarget,
     +      nuc(Ztarget)
          write(1,'(" Elastic scattering angular distribution")')
          write(1,'("# E-incident = ",f7.3)') Einc
          write(1,'("# ")')
          write(1,'("# # angles   =",i3)') nangle+1
          write(1,'("#   E         xs            Direct",$)')
          write(1,'("         Compound    c.s./Rutherford")')     
          do 60 iang=0,nangle
            write(1,'(f5.1,1p,4e16.5)') angle(iang),
     +        discad(k0,Ltarget,iang),directad(k0,Ltarget,iang),
     +        compad(k0,Ltarget,iang),ruth(iang)
   60     continue
          close (unit=1)
        endif    
      endif
c
c ************** Inelastic scattering angular distributions ************
c
c Zindex,Zix: charge number index for residual nucleus
c Nindex,Nix: neutron number index for residual nucleus
c nlev      : number of levels for nucleus
c xsdisc    : total cross section for discrete state
c fileangle : designator for angular distributions on separate file 
c
      Zix=Zindex(0,0,k0)
      Nix=Nindex(0,0,k0)
c
c 1. Legendre coefficients
c
      if (flaglegendre) then
        write(*,'(/"8b1. Legendre coefficients for inelastic",$)')
        write(*,'(" scattering")')
        do 110 i=0,nlev(Zix,Nix)
          if (i.eq.Ltarget) goto 110
          if (xsdisc(k0,i).eq.0.) goto 110
          write(*,'(/"   Level ",i2/)') i
          write(*,'("  L       Total           Direct",$)')
          write(*,'("         Compound       Normalized"/)')
          do 120 LL=0,J2end
            write(*,'(i3,1p,4e16.5)') LL,tleg(k0,i,LL),dleg(k0,i,LL),
     +        cleg(k0,i,LL),tlegnor(k0,i,LL)
  120     continue
c
c Write results to separate file
c
          if (fileangle(i)) then
            discfile='         leg.L00'
            write(discfile(1:2),'(2a1)') parsym(k0),parsym(k0)
            write(discfile(3:9),'(f7.3)') Einc
            write(discfile(3:5),'(i3.3)') int(Einc)
            write(discfile(15:16),'(i2.2)') i
            open (unit=1,status='unknown',file=discfile)  
            write(1,'("# ",a1," + ",i3,a2,$)') parsym(k0),Atarget,
     +        nuc(Ztarget)
            write(1,'(" Inelastic scattering Legendre coefficients",$)')
            write(1,'(" - Level",i3)') i
            write(1,'("# E-incident = ",f7.3)') Einc
            write(1,'("# ")')
            write(1,'("# # coeff.   =",i3)') J2end+1
            write(1,'("#  L       Total           Direct",$)')
            write(1,'("        Compound       Normalized")')
            do 130 LL=0,J2end
              write(1,'(i3,1p,4e16.5)') LL,tleg(k0,i,LL),dleg(k0,i,LL),
     +          cleg(k0,i,LL),tlegnor(k0,i,LL)
  130       continue
            close (unit=1)
          endif    
  110   continue
      endif
c
c 2. Angular distributions
c
      write(*,'(/"8b2. Inelastic angular distributions")')
      do 140 i=1,nlev(Zix,Nix)
        if (i.eq.Ltarget) goto 140
        if (xsdisc(k0,i).eq.0.) goto 140
        write(*,'(/"   Level ",i2/)') i
        write(*,'("Angle       Total         Direct        Compound"/)')
        do 150 iang=0,nangle
          write(*,'(f5.1,1p,3e15.5)') angle(iang),discad(k0,i,iang),
     +      directad(k0,i,iang),compad(k0,i,iang)
  150   continue
c
c Write results to separate file
c
        if (fileangle(i)) then
          discfile='         ang.L00'
          write(discfile(1:2),'(2a1)') parsym(k0),parsym(k0)
          write(discfile(3:9),'(f7.3)') Einc
          write(discfile(3:5),'(i3.3)') int(Einc)
          write(discfile(15:16),'(i2.2)') i
          open (unit=1,status='unknown',file=discfile)  
          write(1,'("# ",a1," + ",i3,a2,$)') parsym(k0),Atarget,
     +      nuc(Ztarget)
          write(1,'(" Inelastic scattering angular distribution",$)')
          write(1,'(" - Level",i3)') i
          write(1,'("# E-incident = ",f7.3)') Einc
          write(1,'("# ")')
          write(1,'("# # angles   =",i3)') nangle+1
          write(1,'("#  E         xs           Direct       Compound")')
          do 160 iang=0,nangle
            write(1,'(f5.1,1p,3e15.5)') angle(iang),discad(k0,i,iang),
     +        directad(k0,i,iang),compad(k0,i,iang)
  160     continue
          close (unit=1)
        endif    
  140 continue
c
c ********** Non-inelastic scattering angular distributions ************
c
c parskip  : logical to skip outgoing particle
c xsdisctot: total cross section summed over discrete states 
c
      write(*,'(/"8c. Angular distributions for other reactions")')
      do 210 type=0,6
        if (parskip(type)) goto 210
        if (type.eq.k0) goto 210
        if (xsdisctot(type).eq.0.) goto 210
        Zix=Zindex(0,0,type)
        Nix=Nindex(0,0,type)
c
c 1. Legendre coefficients
c
        if (flaglegendre) then
          write(*,'(/"8c1. Legendre coefficients for (",a1,",",a1,")")')
     +      parsym(k0),parsym(type)
          do 220 i=0,nlev(Zix,Nix)
            if (xsdisc(type,i).eq.0.) goto 220
            write(*,'(/"   Level ",i2/)') i
            write(*,'("  L       Total           Direct         ",$)')
            write(*,'("Compound       Normalized"/)')
            do 230 LL=0,J2end
              write(*,'(i3,1p,4e16.5)') LL,tleg(type,i,LL),
     +          dleg(type,i,LL),cleg(type,i,LL),tlegnor(type,i,LL)
  230       continue
c
c Write results to separate file
c
            if (fileangle(i)) then
              discfile='         leg.L00'
              write(discfile(1:2),'(2a1)') parsym(k0),parsym(type)
              write(discfile(3:9),'(f7.3)') Einc
              write(discfile(3:5),'(i3.3)') int(Einc)
              write(discfile(15:16),'(i2.2)') i
              open (unit=1,status='unknown',file=discfile)  
              write(1,'("# ",a1," + ",i3,a2,$)') parsym(k0),Atarget,
     +          nuc(Ztarget)
              write(1,'(" (",a1,",",a1,") Legendre coefficients",$)') 
     +          parsym(k0),parsym(type)
              write(1,'(" - Level",i3)') i
              write(1,'("# E-incident = ",f7.3)') Einc
              write(1,'("# ")')
              write(1,'("# # coeff.   =",i3)') J2end+1
              write(1,'("#  L       Total           Direct",$)')
              write(1,'("        Compound       Normalized")')
              do 240 LL=0,J2end
                write(1,'(i3,1p,4e16.5)') LL,tleg(type,i,LL),
     +            dleg(type,i,LL),cleg(type,i,LL),tlegnor(type,i,LL)
  240       continue
            close (unit=1)
          endif    
  220     continue
        endif
c
c 2. Angular distributions
c
        write(*,'(/"8c2. (",a1,",",a1,") angular distributions")') 
     +    parsym(k0),parsym(type)
        do 250 i=0,nlev(Zix,Nix)
          if (xsdisc(type,i).eq.0.) goto 250
          write(*,'(/"   Level ",i2/)') i
          write(*,'("Angle      Total          Direct",$)')
          write(*,'("        Compound"/)')
          do 260 iang=0,nangle
            write(*,'(f5.1,1p,3e15.5)') angle(iang),discad(type,i,iang),
     +        directad(type,i,iang),compad(type,i,iang)
  260     continue
c
c Write results to separate file
c
          if (fileangle(i)) then
            discfile='         ang.L00'
            write(discfile(1:2),'(2a1)') parsym(k0),parsym(type)
            write(discfile(3:9),'(f7.3)') Einc
            write(discfile(3:5),'(i3.3)') int(Einc)
            write(discfile(15:16),'(i2.2)') i
            open (unit=1,status='unknown',file=discfile)  
            write(1,'("# ",a1," + ",i3,a2,$)') parsym(k0),Atarget,
     +        nuc(Ztarget)
            write(1,'(" (",a1,",",a1,") angular distributions",$)') 
     +        parsym(k0),parsym(type)
            write(1,'(" - Level",i3)') i
            write(1,'("# E-incident = ",f7.3)') Einc
            write(1,'("# ")')
            write(1,'("# # angles   =",i3)') nangle+1
            write(1,'("#   E        xs           Direct",$)')
            write(1,'("        Compound")')     
            do 270 iang=0,nangle
              write(1,'(f5.1,1p,3e15.5)') angle(iang),
     +          discad(type,i,iang),directad(type,i,iang),
     +          compad(type,i,iang)
  270       continue
            close (unit=1)
          endif    
  250   continue
  210 continue
      return
      end
Copyright (C) 2004  A.J. Koning, S. Hilaire and M.C. Duijvestijn
