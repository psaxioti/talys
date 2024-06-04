      subroutine directread
c
c +---------------------------------------------------------------------
c | Author: Arjan Koning and Eric Bauge
c | Date  : August 4, 2004
c | Task  : Read ECIS results for direct cross section
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      character*72     line
      integer          type,Zix,Nix,NL,i,l,iang,itype
      real             levelenergy,xsdwbatot
      double precision xs,dl
c
c **************** Read direct, inelastic cross sections ***************
c
c k0              : index of incident particle
c parskip         : logical to skip outgoing particle  
c Zindex,Zix      : charge number index for residual nucleus
c Nindex,Nix      : neutron number index for residual nucleus
c Nlast,NL        : last discrete level
c numlev2         : maximum number of levels
c deform          : deformation parameter 
c eoutdis         : outgoing energy of discrete state reaction
c edis,levelenergy: energy of level
c eninccm         : center-of-mass incident energy in MeV
c parA            : mass number of particle  
c xs              : help variable
c xsdirdisc       : direct cross section for discrete state    
c dorigin         : origin of direct cross section (Direct or Preeq)
c
      open (unit=3,status='unknown',file='ecis97.dircs')
      open (unit=7,status='unknown',file='ecis97.dirres')
      do 10 type=k0,k0
        if (parskip(type)) goto 10
        Zix=Zindex(0,0,type)
        Nix=Nindex(0,0,type)   
        NL=Nlast(Zix,Nix,0)   
c
c 1. Direct collective states
c 
        do 20 i=0,numlev2
          if (i.eq.0.and.type.eq.k0) goto 20               
          if (deform(Zix,Nix,i).eq.0.) goto 20
          if (eoutdis(type,i).le.0.) goto 20 
          levelenergy=edis(Zix,Nix,i)
          if (eninccm.le.levelenergy+0.1*parA(type)) goto 20
  30      read(3,'(a72)',end=120) line
          if (line(1:1).eq.'<') goto 30
          read(line,*) xs       
          xsdirdisc(type,i)=real(xs)
          if (i.le.NL) dorigin(type,i)='Direct'
c
c ******************* Direct reaction Legendre coefficients ************
c
c We read the Legendre coefficients for the direct component of the 
c reaction only, the compound nucleus coefficients are calculated by 
c TALYS later on.
c
c i      : level
c l      : l-value
c dleg,dl: direct reaction Legendre coefficient
c
        
  110     read(7,'(a72)',end=120) line
          if (line(1:1).eq.'<') goto 110
          if (line(5:5).eq.'0') goto 120
          if (line(5:5).eq.'2') then
            read(line,'(5x,i5,e20.10)') l,dl
            if (i.le.NL) dleg(type,i,l)=real(dl)
          endif
          goto 110
  120     backspace 7
c
c ************************ Angular distributions ***********************
c
c nangle  : number of angles
c itype,xs: help variables
c directad: direct angular distribution
c
c We first skip the elastic angular distribution
c
          iang=-1
  210     read(7,'(a72)') line
          if (line(1:1).eq.'<') goto 210
          read(line,'(i3,12x,e12.5)') itype,xs         
          if (itype.eq.0) iang=iang+1
          if (iang.lt.nangle) goto 210
          iang=-1
  220     read(7,'(a72)') line
          if (line(1:1).eq.'<') goto 220
          read(line,'(i3,12x,e12.5)') itype,xs    
          if (itype.eq.0) then
            iang=iang+1
            directad(type,i,iang)=real(xs)
          endif
          if (iang.lt.nangle) goto 220
  230     read(7,'(a72)',end=240) line
          if (line(1:1).eq.'<') goto 230
          if (line(3:3).ne.'0'.and.line(3:3).ne.' ') goto 230
          backspace 7
   20   continue
c
c 2. Giant resonance states
c 
c Egrcoll: energy of giant resonance  
c
  240   if (type.eq.k0) then
          do 310 l=0,3
            do 310 i=1,2
              if (betagr(l,i).eq.0.) goto 310
              levelenergy=Egrcoll(l,i)
              if (eninccm.le.levelenergy+0.1*parA(type)) goto 310
c
c Giant resonance cross section
c
c xsgrcoll: giant resonance cross section
c
  320         read(3,'(a72)') line
              if (line(1:1).eq.'<') goto 320
              read(line,*) xs      
              xsgrcoll(k0,l,i)=real(xs)
c
c Skip the Legendre coefficients
c
  330         read(7,'(a72)',end=340) line
              if (line(1:1).eq.'<') goto 330
              if (line(5:5).eq.'0') goto 340
              goto 330
  340         backspace 7
c
c Giant resonance angular distribution
c
c nanglecont: number of angles for continuum
c grcollad  : giant resonance angular distribution
c
              iang=-1
  350         read(7,'(a72)') line
              if (line(1:1).eq.'<') goto 350
              read(line,'(i3,12x,e12.5)') itype,xs       
              if (itype.eq.0) iang=iang+1
              if (iang.lt.nanglecont) goto 350
              iang=-1
  360         read(7,'(a72)')line
              if (line(1:1).eq.'<') goto 360
              read(line,'(i3,12x,e12.5)') itype,xs   
              if (itype.eq.0) then
                iang=iang+1
                grcollad(k0,l,i,iang)=real(xs)
              endif
              if (iang.lt.nanglecont) goto 360
  370         read(7,'(a72)',end=310) line
              if (line(1:1).eq.'<') goto 370
              if (line(3:3).ne.'0'.and.line(3:3).ne.' ') goto 370
              backspace 7
  310     continue
        endif
c
c ************* Create total direct inelastic cross section ************
c
c xsdwbatot    : direct DWBA cross section summed over discrete states
c xsdirdisctot : direct cross section summed over discrete states
c xscollconttot: total collective cross section in the continuum
c xsdirdiscsum : total direct cross section
c ecisstatus   : status of ECIS file  
c
        xsdwbatot=0.
        do 410 i=0,numlev2
          if (i.eq.0.and.type.eq.k0) goto 410               
          if (deform(Zix,Nix,i).ne.0.) then
            if (i.le.NL) then
              xsdwbatot=xsdwbatot+xsdirdisc(type,i)
            else
              xscollconttot=xscollconttot+xsdirdisc(type,i)
            endif
          endif
  410   continue
        xsdirdisctot(type)=xsdirdisctot(type)+xsdwbatot
        xsdirdiscsum=xsdirdiscsum+xsdwbatot
   10 continue
      close (unit=3,status=ecisstatus)
      close (unit=7,status=ecisstatus)
      return
      end
Copyright (C) 2004  A.J. Koning, S. Hilaire and M.C. Duijvestijn
