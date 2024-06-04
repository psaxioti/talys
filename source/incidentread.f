      subroutine incidentread
c
c +---------------------------------------------------------------------
c | Author: Arjan Koning and Eric Bauge
c | Date  : July 11, 2004
c | Task  : Read ECIS results for incident energy 
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      logical          contlinep,contline
      character*72     line
      character*126    linecc
      integer          infilecs,infiletr,A,Zix,Nix,iread,lparold,J,
     +                 levelcc(4),lcc(4),i,lpar,ispin,l,lev,l1,l2,
     +                 ispin1,ispin2,ii,iang
      real             oddJ,rjcc(4),jres,rj,factor,jres1,jres2
      double precision xs,tjlcc(4),Tcoef,Tcoef1,Tcoef2,dl,teps
c
c ************ Read total, reaction and elastic cross section **********
c
c flaginccalc: flag for new ECIS calculation for incident channel
c infilecs   : file with total cross sections
c infiletr   : file with transmission coefficients, Legendre 
c              coefficients and angular distributions
c line,linecc: line of ECIS-output
c k0         : index of incident particle
c xs         : help variable
c xstotinc   : total cross section (neutrons only) for incident channel
c xsreacinc  : reaction cross section for incident channel
c xsoptinc   : optical model reaction cross section for incident channel
c xselasinc  : total elastic cross section (neutrons only) for incident 
c              channel
c
c If the ECIS calculation has already been done in a previous run, we 
c can read from existing files.
c
      if (flaginccalc) then
        open (unit=3,status='unknown',file='ecis97.inccs')
        open (unit=7,status='unknown',file='ecis97.incres')
        infilecs=3
        infiletr=7
        open (unit=13,status='unknown',file='incident.cs')
        open (unit=17,status='unknown',file='incident.res')
   10   read(3,'(a72)',end=20) line
        write(13,'(a72)') line
        goto 10
   20   rewind 3
        read(7,'(a126)',end=30) linecc
        if (linecc(3:3).ne.' '.and.linecc(8:8).eq.' ') then
          write(17,'(a126)') linecc
        else
          write(17,'(a80)') linecc
        endif
        goto 20
   30   rewind 7
      else
        infilecs=13
        infiletr=17
      endif
      if (k0.eq.1) then
   40   read(infilecs,'(a126)') linecc
        if(linecc(1:1).eq.'<') goto 40
        read(linecc,'(e12.5)') xs
        xstotinc=real(xs)
      endif
   50 read(infilecs,'(a126)') linecc
      if(linecc(1:1).eq.'<') goto 50
      read(linecc,'(e12.5)') xs
      xsreacinc=real(xs)
      xsoptinc=real(xs)
      if (k0.eq.1) then
   60   read(infilecs,'(a126)')linecc
        if(linecc(1:1).eq.'<') goto 60
        read(linecc,'(e12.5)') xs      
        xselasinc=real(xs)
      endif
c
c ******************* Read transmission coefficients *******************
c
c AA,A        : mass number of residual nucleus      
c Zindex,Zix  : charge number index for residual nucleus
c Nindex,Nix  : neutron number index for residual nucleus
c parA        : mass number of particle
c oddJ,iread  : help variables
c lpar,lparold: variable to determine other parity block of transmission
c               coefficients
c
      A=AA(0,0,k0)
      Zix=Zindex(0,0,k0)
      Nix=Nindex(0,0,k0)
      oddJ=0.5*mod(A+parA(k0),2)
      iread=0
      lparold=-1
c
c 1. Spin 1/2 particles: Neutrons, protons, tritons and Helium-3
c
      if (k0.ne.3.and.k0.ne.6) then
c
c A. Transform rotational transmission coefficients into spherical
c    equivalents
c
c colltype    : type of collectivity (D, V or R)
c flagrot     : flag for use of rotational optical model per
c               outgoing particle, if available  
c J           : compound nucleus spin
c contlinep   : help variable to determine new block of data
c levelcc     : level number
c lcc,l       : orbital angular momentum
c rjcc,rj     : j-value
c Tjlinc,tjlcc: transmission coefficients as a function of spin 
c               and l for the incident channel
c factor      : help variable
c targetspin2 : 2 * spin of target
c ispin       : spin index
c
        if (colltype(Zix,Nix).ne.'S'.and.flagrot(k0)) then
          J=-1
          contlinep=.false.
  110     read(infiletr,'(a126)',end=130) linecc
          if (linecc(1:1).eq.'<') goto 110
          if (linecc(5:5).ne.' ') goto 130
          read(linecc,'(4(1x,i2,i4,f6.1,2x,e14.7,3x))')
     +      (levelcc(i),lcc(i),rjcc(i),tjlcc(i),i=1,4)
          if (levelcc(1).ne.1) then
            contlinep=.false.
            goto 110
          endif
          contline=.false.     
          if (levelcc(1).eq.1.and.levelcc(2).eq.1.and.levelcc(3).eq.1.
     +      and.levelcc(4).eq.1) contline=.true.      
          lpar=mod(lcc(1),2)
          if (iread.eq.0) then
            lparold=lpar
            iread=1
          endif
          if (lpar.ne.lparold) then
            J=-1
            lparold=lpar
            contlinep=.false.
          endif                    
          if (.not.contlinep) J=J+1
          contlinep=contline    
          do 120 i=1,4
            if (levelcc(i).eq.1) then
              l=lcc(i)
              rj=rjcc(i)
              factor=(2.*(J+oddJ)+1.)/(2.*rj+1.)/(targetspin2+1.)
              ispin=int(2.*(rj-real(l)))
              Tjlinc(ispin,l)=Tjlinc(ispin,l)+
     +          factor*max(real(tjlcc(i)),0.)
            endif
  120     continue
          goto 110                               
  130     backspace infiletr
        else
c
c B Spherical transmission coefficients
c                                                   
c lev                : level number
c l,l1,l2            : l-value
c jres,jres1,jres2   : j-value
c Tcoef,Tcoef1,Tcoef2: transmission coefficients
c ispin1,ispin2      : spin index
c
  140     read(infiletr,'(a72)',end=150) line
          if (line(1:1).eq.'<') goto 140
          if (line(3:3).ne.'1') goto 150
          read(line,'(i3,i4,f6.1,e16.7)') lev,l,jres,Tcoef
          ispin=int(2.*(jres-real(l)))
          Tjlinc(ispin,l)=max(real(Tcoef),0.)
          goto 140
  150     backspace infiletr
        endif
      endif
c
c 2. Spin 1 particles: Deuterons
c
      if (k0.eq.3) then
c
c A. Transform rotational transmission coefficients into spherical
c    equivalents
c
        if (colltype(Zix,Nix).ne.'S'.and.flagrot(k0)) then
          J=-1
          contlinep=.false.
  160     read(infiletr,'(a126)',end=180) linecc
          if (linecc(1:1).eq.'<') goto 160
          if (linecc(5:5).ne.' ') goto 180
          read(linecc,'(4(1x,i2,i4,f6.1,2x,e14.7,3x))')
     +      (levelcc(i),lcc(i),rjcc(i),tjlcc(i),i=1,4)
          if (levelcc(1).ne.1) then
            contlinep=.false.
            goto 160
          endif
          contline=.false.     
          if (levelcc(1).eq.1.and.levelcc(2).eq.1.and.levelcc(3).eq.1.
     +      and.levelcc(4).eq.1) contline=.true.      
          lpar=mod(lcc(1),2)
          if (iread.eq.0) then
            lparold=lpar
            iread=1
          endif
          if (lpar.ne.lparold) then
            J=-1
            lparold=lpar
            contlinep=.false.
          endif                    
          if (.not.contlinep) J=J+1
          contlinep=contline    
          do 170 i=1,4
            if (levelcc(i).eq.1) then
              l=lcc(i)
              rj=rjcc(i)
              factor=(2.*(J+oddJ)+1.)/(2.*rj+1.)/(targetspin2+1.)
              ispin=int(rj-real(l))
              Tjlinc(ispin,l)=Tjlinc(ispin,l)+
     +          factor*max(real(tjlcc(i)),0.)
            endif
  170     continue
          goto 160                               
  180     backspace infiletr
        else
c
c B Spherical transmission coefficients
c                                                   
  190     read(infiletr,'(a72)',end=200) line
          if (line(1:1).eq.'<') goto 190
          if (line(3:3).ne.'1') goto 200
          read(line,'(2(3x,i4,f6.1,e16.7,3x))') l1,jres1,Tcoef1,l2,
     +      jres2,Tcoef2
          ispin1=int(jres1-real(l1))
          Tjlinc(ispin1,l1)=max(real(Tcoef1),0.)
          if (Tcoef2.ne.0.) then
            ispin2=int(jres2-real(l2))
            Tjlinc(ispin2,l2)=max(real(Tcoef2),0.)
          endif
          goto 190
  200     backspace infiletr
        endif
      endif
c
c 3. Spin 0 particles: Alpha-particles
c
      if (k0.eq.6) then
  210   read(infiletr,'(a72)',end=220) line
        if (line(1:1).eq.'<') goto 210
        if (line(5:5).ne.' ') goto 220
        read(line,'(i3,i4,f6.1,e16.7)') lev,l,jres,Tcoef
        if (lev.ne.1) goto 210
        Tjlinc(0,l)=max(real(Tcoef),0.)
        goto 210
  220   backspace infiletr
      endif
c
c *************** Processing of transmission coefficients **************
c
c Transmission coefficients averaged over spin and determination of
c maximal l-value. ECIS stops its output of transmission coefficients 
c somewhat too early. For the highest l values the transmission 
c coefficient for (l+spin) is not written in the output. Since these 
c are small numbers we put them equal to the value for (l-spin).
c
c numl      : maximal number of l-values
c Tlinc     : transmission coefficients as a function of l for the
c             incident channel, averaged over spin
c translimit: limit for transmission coefficient 
c teps      : help variable
c transeps  : absolute limit for transmission coefficient
c lmaxinc   : maximal l-value for transmission coefficients for incident
c             channel
c
c 1. Spin 1/2 particles: Neutrons, protons, tritons and Helium-3
c
      if (k0.ne.3.and.k0.ne.6) then
        do 310 l=0,numl
          if (Tjlinc(-1,l).ne.0.and.Tjlinc(1,l).eq.0) 
     +      Tjlinc(1,l)=Tjlinc(-1,l)
          if (Tjlinc(-1,l).eq.0.and.Tjlinc(1,l).ne.0.and.l.gt.0) 
     +      Tjlinc(-1,l)=Tjlinc(1,l)
          Tlinc(l)=((l+1)*Tjlinc(1,l)+l*Tjlinc(-1,l))/(2*l+1)
          teps=Tlinc(0)*translimit/(2*l+1)
          teps=max(teps,transeps)
          if (Tjlinc(-1,l).lt.teps.and.Tjlinc(1,l).lt.teps) then
            lmaxinc=l-1
            goto 400
          endif
          lmaxinc=l
  310   continue
      endif
c
c 2. Spin 1 particles: Deuterons
c
      if (k0.eq.3) then
        do 320 l=0,numl
          if (Tjlinc(-1,l).ne.0.and.Tjlinc(0,l).eq.0) 
     +      Tjlinc(0,l)=Tjlinc(-1,l)
          if (Tjlinc(-1,l).ne.0.and.Tjlinc(1,l).eq.0) 
     +      Tjlinc(1,l)=Tjlinc(-1,l)
          if (Tjlinc(-1,l).eq.0.and.Tjlinc(1,l).ne.0.and.l.gt.0) 
     +      Tjlinc(-1,l)=Tjlinc(1,l)
          Tlinc(l)=((2*l+3)*Tjlinc(1,l)+(2*l+1)*Tjlinc(0,l)+
     +      (2*l-1)*Tjlinc(-1,l))/(3*(2*l+1))
          teps=Tlinc(0)*translimit/(2*l+1)
          teps=max(teps,transeps)
          if (Tjlinc(-1,l).lt.teps.and.Tjlinc(0,l).lt.teps.and.
     +      Tjlinc(1,l).lt.teps) then
            lmaxinc=l-1
            goto 400
          endif
          lmaxinc=l
  320   continue
      endif
c
c 3. Spin 0 particles: Alpha-particles
c
      if (k0.eq.6) then
        do 330 l=0,numl
          Tlinc(l)=Tjlinc(0,l)
          teps=Tlinc(0)*translimit/(2*l+1)
          teps=max(teps,transeps)
          if (Tlinc(l).lt.teps) then
            lmaxinc=l-1
            goto 400
          endif
          lmaxinc=l
  330   continue
      endif
c
c ******************* Direct reaction Legendre coefficients ************
c
c We read the Legendre coefficients for the direct component of the 
c reaction only. The compound nucleus coefficients are calculated by 
c TALYS later on. For coupled-channels reactions, the inelastic
c Legendre coefficients are also read.
c
c i      : level number
C l      : l-value
c indexcc: level index for coupled channel
c dleg,dl: direct reaction Legendre coefficient
c
  400 read(infiletr,'(a72)',end=410) line
      if (line(3:3).ne.' ') goto 410
      read(line,'(2i5,e20.10)') i,l,dl
      ii=i-1
      if (i.ne.1) ii=indexcc(Zix,Nix,i)
      dleg(k0,ii,l)=real(dl)
      goto 400
  410 backspace infiletr
c
c ******************* Read elastic angular distributions ***************
c
c iang    : running variable for angle
c nangle  : number of angles
c xs      : help variable
c directad: direct angular distribution
c ruth    : elastic/rutherford ratio
c
c For charged particles, we also read the elastic/rutherford ratio.
c
      do 510 iang=0,nangle
        read(infiletr,'(15x,e12.5)') xs
        directad(k0,0,iang)=real(xs)
        if (k0.gt.1) then
          read(infiletr,'(15x,e12.5)') ruth(iang)
        else
          read(infiletr,'()') 
        endif
  510 continue
c
c For coupled-channels calculations, we also read the discrete inelastic
c angular distributions and cross sections.
c
c xscoupled    : inelastic cross section from coupled channels
c xsdirdisc    : direct cross section for discrete state
c Nlast        : last discrete level   
c xsdirdisctot : direct cross section summed over discrete states
c dorigin      : origin of direct cross section (Direct or Preeq)
c xscollconttot: total collective cross section in the continuum    
c xsdirdiscsum : total direct cross section
c
      if (colltype(Zix,Nix).ne.'S'.and.flagrot(k0)) then
        xscoupled=0.
  520   read(infiletr,'(a72)',end=540) line
        if (line(3:3).ne.' '.and.line(8:8).eq.' ') then
          backspace infiletr
          goto 540
        endif
        backspace infiletr
        do 530 iang=0,nangle
          read(infiletr,'(a72)',end=540) line
          if (line(1:1).eq.'<') then
            backspace infiletr
            goto 540
          endif
          read(line,'(15x,e12.5,i3)',end=540) xs,i
          read(infiletr,'()',end=540) 
          ii=indexcc(Zix,Nix,i+1)
          directad(k0,ii,iang)=real(xs)
  530   continue
        goto 520
  540   read(infilecs,'(a72)',end=560) line
        if (line(1:1).eq.'<') then
          backspace infilecs
          goto 560
        endif
        read(line,'(e12.5,i3)') xs,i
        ii=indexcc(Zix,Nix,i+1)
        xsdirdisc(k0,ii)=real(xs)
        if (ii.le.Nlast(Zix,Nix,0)) then
          xsdirdisctot(k0)=xsdirdisctot(k0)+xsdirdisc(k0,ii)
        else
          xscollconttot=xscollconttot+xsdirdisc(k0,ii)
        endif       
        xscoupled=xscoupled+xsdirdisc(k0,ii)
        dorigin(k0,ii)='Direct'
        goto 540
  560   xsdirdiscsum=xsdirdisctot(k0)
      endif
c
c Close files
c
c ecisstatus: status of ECIS file  
c nin       : counter for incident energy
c numinc    : number of incident energies
c
      close (unit=3,status=ecisstatus)
      close (unit=7,status=ecisstatus)
      if (nin.eq.numinc) then
        close (unit=13,status=ecisstatus)
        close (unit=17,status=ecisstatus)
      endif
      return
      end
Copyright (C) 2004  A.J. Koning, S. Hilaire and M.C. Duijvestijn
