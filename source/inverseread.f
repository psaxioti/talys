      subroutine inverseread(Zcomp,Ncomp)
c
c +---------------------------------------------------------------------
c | Author: Stephane Hilaire, Arjan Koning and Eric Bauge
c | Date  : July 7, 2004
c | Task  : Read ECIS results for outgoing particles and energy grid
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      logical          lexist,contline,contlinep
      character*126    linein
      integer          Zcomp,Ncomp,type,nen,A,Zix,Nix,iread,J,
     +                 levelcc(4),lcc(4),lpar,lparold,lpar0,i,l,l1,l2,
     +                 ispin,lev,ispin1,ispin2
      real             oddJ,rjcc(4),rj,factor,jres,jres1,jres2
      double precision xs,tjlcc(4),Tcoef,Tcoef1,Tcoef2,teps
c
c ************ Read total, reaction and elastic cross section **********
c
c Zcomp     : charge number index for compound nucleus
c Ncomp     : neutron number index for compound nucleus
c csfile    : file with inverse reaction cross sections      
c parskip   : logical to skip outgoing particle
c ebegin    : first energy point of energy grid
c eendmax   : last energy point of energy grid for maximum incident 
c             energy
c linein    : input line
c xs        : help variable
c xstot     : total cross section (neutrons only) for
c xsreac    : reaction cross section
c xsopt     : optical model reaction cross section
c xselas    : total elastic cross section (neutrons only)
c ecisstatus: status of ECIS file  
c
      inquire (file=csfile,exist=lexist)
      if (.not.lexist) then
        write(*,'("TALYS-error: The first calculation of a run",$)')
        write(*,'(" should always be done with ecissave y and",$)')
        write(*,'(" eciscalc y")')
        stop
      endif    
      open (unit=3,status='unknown',file=csfile)
      do 10 type=1,6
        if (parskip(type)) goto 10
        do 20 nen=ebegin(type),eendmax(type)
          if (type.eq.1) then
   30       read(3,'(a126)') linein
            if(linein(1:1).eq.'<') goto 30
            read(linein,*) xs         
            xstot(type,nen)=real(xs)
          endif
   40     read(3,'(a126)') linein
          if(linein(1:1).eq.'<') goto 40
          read(linein,*) xs         
          xsreac(type,nen)=real(xs)
          xsopt(type,nen)=real(xs)
          if (type.eq.1) then
   50       read(3,'(a126)') linein
            if(linein(1:1).eq.'<') goto 50
            read(linein,*) xs         
            xselas(type,nen)=real(xs)
          endif
   20   continue
   10 continue
      close (unit=3,status=ecisstatus)
c
c ******************* Read transmission coefficients *******************
c
c transfile : file with transmission coefficients
c AA,A      : mass number of residual nucleus
c Zindex    : charge number index for residual nucleus
c Nindex    : neutron number index for residual nucleus    
c oddJ,iread: help variable
c
c The transmission coefficient Tjl has four indices: particle type,
c energy, spin and  l-value. For spin-1/2 particles, we use the
c indices -1 and 1 for the two spin values. For spin-1 particles,
c we use -1, 0 and 1 and for spin-0 particles we use 0 only.
c
      open(unit=7,status='unknown',file=transfile)
      A=AA(Zcomp,Ncomp,0)
      oddJ=0.5*mod(A,2)
      do 110 type=1,6
        if (parskip(type)) goto 110
        Zix=Zindex(Zcomp,Ncomp,type)
        Nix=Nindex(Zcomp,Ncomp,type)
        nen=ebegin(type)-1
        if (nen.gt.eendmax(type)) goto 110
        iread=0
c
c 1. Spin 1/2 particles: Neutrons, protons, tritons and Helium-3
c
        if (type.ne.3.and.type.ne.6) then
c
c A. Only for rotational nuclei: Transform rotational transmission 
c    coefficients into spherical equivalents.
c
c colltype     : type of collectivity (D, V or R)
c flagrot      : flag for use of rotational optical model per
c                outgoing particle, if available  
c J            : compound nucleus spin
c contlinep    : help variable
c contline     : help variable
c levelcc      : level number
c lcc,l        : orbital angular momentum
c rjcc,rj      : j-value
c Tjl,tjlcc    : transmission coefficients as a function of
c                particle type, energy, spin and l-value
c lpar         : variable to determine other parity block of 
c                transmission coefficients
c lparold,lpar0: help variables
c factor       : help variable
c targetspin2  : 2 * spin of target
c ispin        : spin index
c 
          if (colltype(Zix,Nix).ne.'S'.and.flagrot(type)) then
            J=-1
            contlinep=.false.
  120       read(7,'(a126)',end=110) linein
            if (linein(1:1).eq.'<') goto 120
            read(linein,'(4(1x,i2,i4,f6.1,2x,e14.7,3x))',end=110)
     +        (levelcc(i),lcc(i),rjcc(i),tjlcc(i),i=1,4)
            if (levelcc(1).ne.1) then
              contlinep=.false.
              goto 120
            endif
            contline=.false.
            if (levelcc(1).eq.1.and.levelcc(2).eq.1.and.levelcc(3).eq.1.
     +        and.levelcc(4).eq.1) contline=.true.  
            lpar=mod(lcc(1),2)
            if (iread.eq.0) then
              lparold=lpar
              lpar0=lpar
              iread=1
              nen=nen+1
            endif
            if (lpar.ne.lparold.and.lpar.eq.lpar0) nen=nen+1
            if (nen.gt.eendmax(type)) then
              backspace 7
              goto 110
            endif
            if (lpar.ne.lparold) then
              J=-1
              lparold=lpar 
              contlinep=.false.
            endif
            if (.not.contlinep) J=J+1
            contlinep=contline
            do 130 i=1,4
              if (levelcc(i).eq.1) then
                l=lcc(i)
                rj=rjcc(i)
                factor=(2.*(J+oddJ)+1.)/(2.*rj+1.)/(targetspin2+1.)
                ispin=int(2.*(rj-real(l)))
                Tjl(type,nen,ispin,l)=Tjl(type,nen,ispin,l)+
     +            factor*max(real(tjlcc(i)),0.)
              endif
  130       continue
            goto 120
          else
c
c B Spherical transmission coefficients
c
c lev                : level number
c l,l1,l2            : l-value
c jres,jres1,jres2   : j-value
c Tcoef,Tcoef1,Tcoef2: transmission coefficients
c
  140       read(7,'(a126)',end=110) linein
            if (linein(1:1).eq.'<') goto 140
            read(linein,'(i3,i4,f6.1,e16.7)',end=110) lev,l,jres,Tcoef
            if (l.eq.0) nen=nen+1
            if (nen.gt.eendmax(type)) then
              backspace 7
              goto 110
            endif
            ispin=int(2.*(jres-real(l)))
            Tjl(type,nen,ispin,l)=max(real(Tcoef),0.)
            goto 140
          endif                     
        endif
c
c 2. Spin 1 particles: Deuterons
c
c A. Only for rotational nuclei: Transform rotational transmission 
c    coefficients into spherical equivalents.
c
        if (type.eq.3) then
          if (colltype(Zix,Nix).ne.'S'.and.flagrot(type)) then
            J=-1
            contlinep=.false.
  150       read(7,'(a126)',end=110) linein
            if (linein(1:1).eq.'<') goto 150
            read(linein,'(4(1x,i2,i4,f6.1,2x,e14.7,3x))',end=110) 
     +        (levelcc(i),lcc(i),rjcc(i),tjlcc(i),i=1,4)
            if (levelcc(1).ne.1) then
              contlinep=.false.
              goto 150
            endif
            contline=.false.
            if (levelcc(1).eq.1.and.levelcc(2).eq.1.and.levelcc(3).eq.1.
     +        and.levelcc(4).eq.1) contline=.true.  
            lpar=mod(lcc(1),2)
            if (iread.eq.0) then
              lparold=lpar
              lpar0=lpar
              iread=1
              nen=nen+1
            endif
            if (lpar.ne.lparold.and.lpar.eq.lpar0) nen=nen+1
            if (nen.gt.eendmax(type)) then
              backspace 7
              goto 110
            endif
            if (lpar.ne.lparold) then
              J=-1
              lparold=lpar 
              contlinep=.false.
            endif
            if (.not.contlinep) J=J+1
            contlinep=contline
            do 160 i=1,4
              if (levelcc(i).eq.1) then
                l=lcc(i)
                rj=rjcc(i)
                factor=(2.*(J+oddJ)+1.)/(2.*rj+1.)/(targetspin2+1.)
                ispin=int(rj-real(l))
                Tjl(type,nen,ispin,l)=Tjl(type,nen,ispin,l)+
     +            factor*max(real(tjlcc(i)),0.)
              endif
  160       continue
            goto 150
          else
c
c B Spherical transmission coefficients
c
  170       read(7,'(a126)',end=110) linein
            if (linein(1:1).eq.'<') goto 170
            read(linein,'(2(3x,i4,f6.1,e16.7,3x))',end=110) l1,jres1,
     +        Tcoef1,l2,jres2,Tcoef2
            if (l1.eq.0) nen=nen+1
            if (nen.gt.eendmax(type)) then
              backspace 7
              goto 110
            endif
            ispin1=int(jres1-real(l1))
            Tjl(type,nen,ispin1,l1)=max(real(Tcoef1),0.)
            if (Tcoef2.ne.0.) then
              ispin2=int(jres2-real(l2))
              Tjl(type,nen,ispin2,l2)=max(real(Tcoef2),0.)
            endif
            goto 170
          endif                  
        endif                  
c
c 3. Spin 0 particles: Alpha-particles
c
        if (type.eq.6) then
  180     read(7,'(a126)',end=110) linein
          if (linein(1:1).eq.'<') goto 180
          read(linein,'(i3,i4,f6.1,e16.7)',end=110) lev,l,jres,Tcoef
          if (lev.ne.1) goto 180
          if (l.eq.0) nen=nen+1
          if (nen.gt.eendmax(type)) then
            backspace 7
            goto 110
          endif
          Tjl(type,nen,0,l)=max(real(Tcoef),0.)
          goto 180
        endif
  110 continue
c
c ************** Processing of transmission coefficients ***************
c
c Transmission coefficients averaged over spin and determination of
c maximal l-value. ECIS stops its output of transmission coefficients 
c somewhat too early. For the highest l values the transmission 
c coefficient for (l+spin) is not written in the output. Since these 
c are small numbers we put them equal to the value for (l-spin).
c
c eend      : last energy point of energy grid 
c numl      : maximal number of l-values in TALYS
c Tl        : transmission coefficients as a function of
c             particle type, energy and l-value (averaged over spin)
c translimit: limit for transmission coefficient 
c teps      : help variable
c transeps  : limit for transmission coefficient
c lmax      : maximal l-value for transmission coefficients
c
      do 210 type=1,6
        if (parskip(type)) goto 210
c
c 1. Spin 1/2 particles: Neutrons, protons, tritons and Helium-3
c
        if (type.ne.3.and.type.ne.6) then
          do 220 nen=ebegin(type),eend(type)
            do 230 l=0,numl
              if (Tjl(type,nen,-1,l).ne.0.and.Tjl(type,nen,1,l).eq.0) 
     +          Tjl(type,nen,1,l)=Tjl(type,nen,-1,l)
              if (Tjl(type,nen,-1,l).eq.0.and.Tjl(type,nen,1,l).ne.0.
     +          and.l.gt.0) Tjl(type,nen,-1,l)=Tjl(type,nen,1,l)
              Tl(type,nen,l)=((l+1)*Tjl(type,nen,1,l)+
     +          l*Tjl(type,nen,-1,l))/(2*l+1)
              teps=Tl(type,nen,0)*translimit/(2*l+1)
              teps=max(teps,transeps)
              if (Tjl(type,nen,-1,l).lt.teps.and.Tjl(type,nen,1,l).
     +          lt.teps) then
                lmax(type,nen)=l-1
                goto 220
              endif
              lmax(type,nen)=l
  230       continue
  220     continue
        endif
c
c 2. Spin 1 particles: Deuterons
c
        if (type.eq.3) then
          do 240 nen=ebegin(type),eend(type)
            do 250 l=0,numl
              if (Tjl(type,nen,-1,l).ne.0.and.Tjl(type,nen,0,l).eq.0) 
     +          Tjl(type,nen,0,l)=Tjl(type,nen,-1,l)
              if (Tjl(type,nen,-1,l).ne.0.and.Tjl(type,nen,1,l).eq.0) 
     +          Tjl(type,nen,1,l)=Tjl(type,nen,-1,l)
              if (Tjl(type,nen,-1,l).eq.0.and.Tjl(type,nen,1,l).ne.0.
     +          and.l.gt.0) 
     +          Tjl(type,nen,-1,l)=Tjl(type,nen,1,l)
              Tl(type,nen,l)=((2*l+3)*Tjl(type,nen,1,l)+
     +          (2*l+1)*Tjl(type,nen,0,l)+
     +          (2*l-1)*Tjl(type,nen,-1,l))/(3*(2*l+1))
              teps=Tl(type,nen,0)*translimit/(2*l+1)
              teps=max(teps,transeps)
              if (Tjl(type,nen,-1,l).lt.teps.and.Tjl(type,nen,0,l).
     +          lt.teps.and.Tjl(type,nen,1,l).lt.teps) then
                lmax(type,nen)=l-1
                goto 240
              endif
              lmax(type,nen)=l
  250       continue
  240     continue
        endif
c
c 3. Spin 0 particles: Alpha-particles
c
        if (type.eq.6) then
          do 260 nen=ebegin(type),eend(type)
            do 270 l=0,numl
              Tl(type,nen,l)=Tjl(type,nen,0,l)
              teps=Tl(type,nen,0)*translimit/(2*l+1)
              teps=max(teps,transeps)
              if (Tl(type,nen,l).lt.teps) then
                lmax(type,nen)=l-1
                goto 260
              endif
              lmax(type,nen)=l
  270       continue
  260     continue
        endif
  210 continue
      close (unit=7,status=ecisstatus)
      return
      end
Copyright (C) 2004  A.J. Koning, S. Hilaire and M.C. Duijvestijn
