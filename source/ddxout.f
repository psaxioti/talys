      subroutine ddxout
c
c +---------------------------------------------------------------------
c | Author: Arjan Koning 
c | Date  : September 1, 2004
c | Task  : Output of double-differential cross sections
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      character*13 ddxfile
      integer      type,nen,iang,i
      real         Eo(0:numen2),enf,fac,xsa,xsb,xs1,xs2,xs3,xs4,xs5,
     +             Eout,angf
c
c ******************* Double-differential cross sections ***************
c
c ddxmode      : mode for double-differential cross sections: 0: None,
c                1: Angular distributions, 2: Spectra per angle, 3: Both
c parskip      : logical to skip outgoing particle
c xsparticle   : total particle production cross section
c ebegin       : first energy point of energy grid
c eendout      : last energy point of energy grid
c espec,Eo     : outgoing energy grid
c xssumout     : cross section summed over mechanisms
c parname      : name of particle
c nanglecont   : number of angles for continuum
c anglecont    : angle in degrees for continuum
c xssumoutad   : angular distribution summed over mechanisms
c xsdiscoutad  : smoothed angular distribution for discrete state
c xspreeqoutad : preequilibrium angular distribution per particle type
c xsmpreeqoutad: multiple preequilibrium angular distribution
c xscompoutad  : compound emission angular distribution
c flagrecoil   : flag for calculation of recoils              
c flaglabddx   : flag for calculation of DDX in LAB system
c iejlab       : number of ejectile lab bins  
c Eejlab       : center of ejectile lab bin                
c ddxejlab     : LAB ddx array
c ddxecount    : counter for double-differential cross section files
c fileddxe     : designator for double-differential cross sections on
c                separate file: angular distribution
c locate       : subroutine to find value in ordered table
c deltaE       : energy bin around outgoing energies  
c parsym       : symbol of particle 
c Einc         : incident energy in MeV
c xsa,....     : help variables
c
c 1. Angular distributions per outgoing energy
c
      if (ddxmode.eq.1.or.ddxmode.eq.3) then
        write(*,'(/"9. Double-differential cross sections per",$)')
        write(*,'(" outgoing energy")') 
        do 10 type=1,6
          if (parskip(type)) goto 10
          if (xsparticle(type).eq.0.) goto 10
          do 20 nen=ebegin(type),eendout(type)
            Eo(nen)=espec(type,nen)
            if (xssumout(type,nen).eq.0.) goto 20
            write(*,'(/"DDX for outgoing ",a8," at ",f7.3," MeV"/)') 
     +        parname(type),Eo(nen)
            write(*,'("Angle   Total      Direct  ",$)')
            write(*,'("   Pre-equil.  Mult. preeq   Compound"/)')
            do 30 iang=0,nanglecont
              write(*,'(f5.1,1p,5e12.5)') anglecont(iang),
     +          xssumoutad(type,nen,iang),xsdiscoutad(type,nen,iang),
     +          xspreeqoutad(type,nen,iang),
     +          xsmpreeqoutad(type,nen,iang),xscompoutad(type,nen,iang)
   30       continue
   20     continue
          if (flagrecoil.and.flaglabddx) then
            do 40 nen=1,iejlab(type)
              write(*,'(/"DDX for outgoing ",a8," at ",f7.3," MeV",$)') 
     +          parname(type),Eejlab(type,nen)
              write(*,'(" in LAB frame")')
              write(*,'("Angle   Cross section"/)')
              do 50 iang=0,nanglecont   
                write(*,'(f5.1,1p,e12.5)') anglecont(iang),
     +            ddxejlab(type,nen,iang)
   50         continue
   40       continue
          endif
          do 60 i=1,ddxecount(type)
            enf=fileddxe(type,i)
            call locate(Eo,ebegin(type),eendout(type),enf,nen)
            fac=(enf-Eo(nen))/deltaE(nen)
            ddxfile=' ddx000_0.mev'
            write(ddxfile(1:1),'(a1)') parsym(type)
            write(ddxfile(5:9),'(f5.1)') enf
            write(ddxfile(5:7),'(i3.3)') int(enf)
            open (unit=1,status='unknown',file=ddxfile)     
            write(1,'("# ",a1," + ",i3,a2,": ",a8," DDX spectrum")')
     +        parsym(k0),Atarget,nuc(Ztarget),parname(type)
            write(1,'("# E-incident = ",f7.3)') Einc
            write(1,'("# ")')
            write(1,'("# # angles =",i3)') nanglecont+1
            write(1,'("# Angle    Total       Direct    Pre-equil.",$)')
            write(1,'("  Mult. preeq  Compound")')      
            do 70 iang=0,nanglecont
              xsa=xssumoutad(type,nen,iang)
              xsb=xssumoutad(type,nen+1,iang)
              xs1=xsa+fac*(xsb-xsa)
              xsa=xsdiscoutad(type,nen,iang)
              xsb=xsdiscoutad(type,nen+1,iang)
              xs2=xsa+fac*(xsb-xsa)
              xsa=xspreeqoutad(type,nen,iang)
              xsb=xspreeqoutad(type,nen+1,iang)
              xs3=xsa+fac*(xsb-xsa)
              xsa=xsmpreeqoutad(type,nen,iang)
              xsb=xsmpreeqoutad(type,nen+1,iang)
              xs4=xsa+fac*(xsb-xsa)
              xsa=xscompoutad(type,nen,iang)
              xsb=xscompoutad(type,nen+1,iang)
              xs5=xsa+fac*(xsb-xsa)
              write(1,'(f5.1,1p,5e12.5)') anglecont(iang),xs1,xs2,xs3,
     +          xs4,xs5
   70       continue
            close (unit=1)
   60     continue
   10   continue
      endif
c
c 2. Emission spectra per outgoing angle
c
c anginc   : angle increment
c ddxacount: counter for double-differential cross section files
c fileddxa : designator for double-differential cross sections on
c            separate file: spectrum per angle
c
      if (ddxmode.eq.2.or.ddxmode.eq.3) then
        anginc=180./nanglecont
        write(*,'(/"9. Double-differential cross sections per",$)')
        write(*,'(" outgoing angle")') 
        do 110 type=1,6
          if (parskip(type)) goto 110
          if (xsparticle(type).eq.0.) goto 110
          do 120 iang=0,nanglecont
            write(*,'(/"DDX for outgoing ",a8," at ",f7.3," degrees")') 
     +        parname(type),anglecont(iang)
            write(*,'(/"  E-out    Total      Direct  ",$)')
            write(*,'("   Pre-equil. Mult. preeq   Compound"/)')
            do 130 nen=ebegin(type),eendout(type)
              Eout=espec(type,nen)
              write(*,'(f7.3,1p,5e12.5)') Eout,
     +          xssumoutad(type,nen,iang),xsdiscoutad(type,nen,iang),
     +          xspreeqoutad(type,nen,iang),
     +          xsmpreeqoutad(type,nen,iang),xscompoutad(type,nen,iang)
  130       continue
  120     continue
          if (flagrecoil.and.flaglabddx) then
            do 140 iang=0,nanglecont
              write(*,'(/"DDX for outgoing ",a8," at ",f7.3,$)') 
     +          parname(type),anglecont(iang)
              write(*,'(" degrees in LAB frame"/)')
              write(*,'("Energy  Cross section"/)')
              do 150 nen=1,iejlab(type)
                write(*,'(f7.3,1p,e12.5)') Eejlab(type,nen),
     +            ddxejlab(type,nen,iang)
  150         continue
  140       continue
          endif
          do 160 i=1,ddxacount(type)
            angf=fileddxa(type,i)
            call locate(anglecont,0,nanglecont,angf,iang)
            fac=(angf-anglecont(iang))/anginc
            ddxfile=' ddx000_0.deg'
            write(ddxfile(1:1),'(a1)') parsym(type)
            write(ddxfile(5:9),'(f5.1)') angf
            write(ddxfile(5:7),'(i3.3)') int(angf)
            open (unit=1,status='unknown',file=ddxfile)     
            write(1,'("# ",a1," + ",i3,a2,": ",a8," DDX spectrum")')
     +        parsym(k0),Atarget,nuc(Ztarget),parname(type)
            write(1,'("# E-incident = ",f7.3)') Einc
            write(1,'("# ")')
            write(1,'("# # energies =",i3)')
     +        eendout(type)-ebegin(type)+1
            write(1,'("# E-out    Total       Direct    Pre-equil.",$)')
            write(1,'("  Mult. preeq  Compound")')      
            do 170 nen=ebegin(type),eendout(type)
              Eout=espec(type,nen)
              xsa=xssumoutad(type,nen,iang)
              xsb=xssumoutad(type,nen,iang+1)
              xs1=xsa+fac*(xsb-xsa)
              xsa=xsdiscoutad(type,nen,iang)
              xsb=xsdiscoutad(type,nen,iang+1)
              xs2=xsa+fac*(xsb-xsa)
              xsa=xspreeqoutad(type,nen,iang)
              xsb=xspreeqoutad(type,nen,iang+1)
              xs3=xsa+fac*(xsb-xsa)
              xsa=xsmpreeqoutad(type,nen,iang)
              xsb=xsmpreeqoutad(type,nen,iang+1)
              xs4=xsa+fac*(xsb-xsa)
              xsa=xscompoutad(type,nen,iang)
              xsb=xscompoutad(type,nen,iang+1)
              xs5=xsa+fac*(xsb-xsa)
              write(1,'(f7.3,1p,5e12.5)') Eout,xs1,xs2,xs3,xs4,xs5
  170       continue
            close (unit=1)
  160     continue
  110   continue
      endif
      return
      end
Copyright (C) 2004  A.J. Koning, S. Hilaire and M.C. Duijvestijn
