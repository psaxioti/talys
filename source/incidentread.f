      subroutine incidentread
c
c +---------------------------------------------------------------------
c | Author: Arjan Koning, Eric Bauge and Pascal Romain
c | Date  : September 14, 2006
c | Task  : Read ECIS results for incident energy 
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      character*72     line
      integer          infilecs,infiletr,infileang,infileleg,infilein,
     +                 Zix,Nix,i,nJ,k,nS,lev,l,ispin,nL,ii,nSt,iSt,iang,
     +                 itype
      real             groundspin2,rj,jres,factor,xsr
      double precision xs,Tcoef,dl,teps
c
c ************ Read total, reaction and elastic cross section **********
c
c flaginccalc: flag for new ECIS calculation for incident channel
c infilecs   : file with total cross sections
c infiletr   : file with transmission coefficients
c infileang  : file with angular distributions
c infileleg  : file with Legendre coefficients
c line       : line of ECIS-output
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
        open (unit=3,status='unknown',file='ecis03.inccs')
        open (unit=7,status='unknown',file='ecis03.inctr')
        open (unit=8,status='unknown',file='ecis03.incang')
        open (unit=9,status='unknown',file='ecis03.incleg')
        open (unit=10,status='unknown',file='ecis03.incin')
        infilecs=3
        infiletr=7
        infileang=8
        infileleg=9
        infilein=10
        open (unit=13,status='unknown',file='incident.cs')
        open (unit=17,status='unknown',file='incident.tr')
        open (unit=18,status='unknown',file='incident.ang')
        open (unit=19,status='unknown',file='incident.leg')
        open (unit=20,status='unknown',file='incident.in')
   10   read(3,'(a72)',end=20) line
        write(13,'(a72)') line
        goto 10
   20   rewind 3
        read(7,'(a72)',end=30) line
        write(17,'(a72)') line
        goto 20
   30   rewind 7
   40   read(8,'(a72)',end=50) line
        write(18,'(a72)') line
        goto 40
   50   rewind 8
   60   read(9,'(a72)',end=70) line
        write(19,'(a72)') line
        goto 60
   70   rewind 9
   80   read(10,'(a72)',end=90) line
        write(20,'(a72)') line
        goto 80
   90   rewind 10
      else
        infilecs=13
        infiletr=17
        infileang=18
        infileleg=19
        infilein=20
      endif
      read(infilecs,'()') 
      if (k0.eq.1) then
        read(infilecs,*) xs
        xstotinc=real(xs)
      endif
      read(infilecs,*) xs
      xsreacinc=real(xs)
      xsoptinc=real(xs)
      if (k0.eq.1) then
        read(infilecs,*) xs
        xselasinc=real(xs)
      endif
c
c ******************* Read transmission coefficients *******************
c
c Zindex,Zix : charge number index for residual nucleus
c Nindex,Nix : neutron number index for residual nucleus
c groundspin2: 2 * spin of ground state
c jdis       : spin of level
c nJ         : number of total J values for transmission coeffients
c rJ         : compound nucleus spin
c nS         : number of transmission coeffients per J-value
c lev        : level number
c l          : orbital angular momentum
c jres       : j-value
c Tjl,Tcoef  : transmission coefficients as a function of spin and 
c              l-value
c colltype   : type of collectivity (D, V or R)
c flagrot    : flag for use of rotational optical model per
c              outgoing particle, if available  
c factor     : help variable
c ispin      : spin index
c parspin    : spin of particle
c
c For rotational nuclei, the rotational transmission 
c coefficients are transformed into into spherical equivalents.
c
      Zix=Zindex(0,0,k0)
      Nix=Nindex(0,0,k0)
      groundspin2=int(2.*jdis(Zix,Nix,0))
      read(infiletr,'(45x,i5)') nJ
      do 110 i=1,nJ
        read(infiletr,'(f5.1,2x,i5)') rj,nS
        do 120 k=1,nS
          read(infiletr,'(i3,i4,f6.1,e16.7)') lev,l,jres,Tcoef
          if (lev.eq.1) then
            if (colltype(Zix,Nix).ne.'S'.and.flagrot(k0)) then
              factor=(2.*rj+1.)/(2.*jres+1.)/(groundspin2+1.)
            else
              factor=1.
            endif
            if (parspin(k0).eq.0.5) then
              ispin=int(2.*(jres-real(l)))
            else
              ispin=int(jres-real(l))
            endif
            Tjlinc(ispin,l)=Tjlinc(ispin,l)+factor*max(real(Tcoef),0.)
          endif
  120   continue
  110 continue
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
        do 210 l=0,numl
          if (Tjlinc(-1,l).ne.0.and.Tjlinc(1,l).eq.0) 
     +      Tjlinc(1,l)=Tjlinc(-1,l)
          if (Tjlinc(-1,l).eq.0.and.Tjlinc(1,l).ne.0.and.l.gt.0) 
     +      Tjlinc(-1,l)=Tjlinc(1,l)
          Tlinc(l)=((l+1)*Tjlinc(1,l)+l*Tjlinc(-1,l))/(2*l+1)
          teps=Tlinc(0)*translimit/(2*l+1)
          teps=max(teps,transeps)
          if (Tjlinc(-1,l).lt.teps.and.Tjlinc(1,l).lt.teps) then
            lmaxinc=l-1
            goto 300
          endif
          lmaxinc=l
  210   continue
      endif
c
c 2. Spin 1 particles: Deuterons
c
      if (k0.eq.3) then
        do 220 l=0,numl
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
            goto 300
          endif
          lmaxinc=l
  220   continue
      endif
c
c 3. Spin 0 particles: Alpha-particles
c
      if (k0.eq.6) then
        do 230 l=0,numl
          Tlinc(l)=Tjlinc(0,l)
          teps=Tlinc(0)*translimit/(2*l+1)
          teps=max(teps,transeps)
          if (Tlinc(l).lt.teps) then
            lmaxinc=l-1
            goto 300
          endif
          lmaxinc=l
  230   continue
      endif
c
c ******************* Direct reaction Legendre coefficients ************
c
c We read the Legendre coefficients for the direct component of the 
c reaction only. The compound nucleus coefficients are calculated by 
c TALYS later on. For coupled-channels reactions, the inelastic
c Legendre coefficients are also read.
c
c nL     : number of Legendre coefficients
c i      : level number
c l      : l-value
c indexcc: level index for coupled channel
c dleg,dl: direct reaction Legendre coefficient
c
  300 read(infileleg,'()') 
      read(infileleg,'(5x,i5)')  nL
  310 do 320 k=1,nL
        read(infileleg,'(2i5,e20.10)') i,l,dl
        ii=i-1
        if (i.ne.1) ii=indexcc(Zix,Nix,i)
        dleg(k0,ii,l)=real(dl)
  320 continue
      read(infileleg,'(a72)',end=400) line
      if (line(2:9).eq.'legendre') then
        backspace infileleg
        goto 400
      endif
      read(line,'(5x,i5)') nL
      goto 310
        
c
c ******************* Read elastic angular distributions ***************
c
c nSt     : number of states
c iang    : running variable for angle
c nangle  : number of angles
c xs      : help variable
c directad: direct angular distribution
c ruth    : elastic/rutherford ratio
c
c For charged particles, we also read the elastic/rutherford ratio.
c
  400 read(infileang,'(45x,i5)') nSt
      read(infileang,'(12x,i3)') nS
      do 410 iang=0,nangle
        do 410 k=1,nS
          read(infileang,'(i3,12x,e12.5)') itype,xs
          xsr=1.e38
          if (xs.le.1.e38) xsr=real(xs)
          if (itype.eq.0) directad(k0,0,iang)=xsr
          if (itype.eq.1) ruth(iang)=xsr
  410 continue
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
        do 510 ist=2,nSt
          read(infileang,'(i5,7x,i3)') i,nS
          do 520 iang=0,nangle
            do 530 k=1,nS
              read(infileang,'(i3,12x,e12.5)') itype,xs
              ii=indexcc(Zix,Nix,i)
              if (itype.eq.0) directad(k0,ii,iang)=real(xs)
  530       continue
  520     continue
  510   continue
        read(infilein,'(a72)',end=560) line
  540   read(infilein,'(a72)',end=560) line
        if (line(1:1).eq.'<') then
          backspace infilein
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
      close (unit=8,status=ecisstatus)
      close (unit=9,status=ecisstatus)
      close (unit=10,status=ecisstatus)
      if (nin.eq.numinc) then
        close (unit=13,status=ecisstatus)
        close (unit=17,status=ecisstatus)
        close (unit=18,status=ecisstatus)
        close (unit=19,status=ecisstatus)
        close (unit=20,status=ecisstatus)
      endif
      return
      end
Copyright (C) 2004  A.J. Koning, S. Hilaire and M.C. Duijvestijn
