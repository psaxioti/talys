      subroutine input2
c
c +---------------------------------------------------------------------
c | Author: Arjan Koning
c | Date  : December 2, 2004
c | Task  : Read input for second set of variables
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      character*1 ch
      integer     type,Zix,Nix,i,ip,i2,iz,ia,ivalue,type2
c
c ************* Defaults for second set of input variables *************
c
c outtype   : type of outgoing particles
c maxZ,numZ : maximal number of protons away from the initial 
c             compound nucleus
c maxN,numN : maximal number of neutrons away from the initial 
c             compound nucleus
c nbins     : number of continuum excitation energy bins
c segment   : number of segments to divide emission energy grid
c nlevmax   : maximum number of included discrete levels for target
c nlevmaxres: maximum number of included discrete levels for residual
c             nucleus
c nlevbin   : number of excited levels for binary nucleus  
c k0        : index of incident particle 
c Ltarget   : excited level of target
c isomer    : definition of isomer in seconds
c core      : even-even core for weakcoupling (-1 or 1)
c gammax    : number of l-values for gamma multipolarity
c transpower: power for transmission coefficient limit
c transeps  : absolute limit for transmission coefficient
c xseps     : limit for cross sections 
c popeps    : limit for population cross section per nucleus
c Rfiseps   : ratio for limit for fission cross section per nucleus
c eninclow  : minimal incident energy for nuclear model calculations
c nangle    : number of angles
c numang    : maximum number of angles
c nanglecont: number of angles for continuum
c maxenrec  : number of recoil energies
c ldmodel   : level density model
c wmode     : designator for width fluctuation model
c preeqmode : designator for pre-equilibrium model
c mpreeqmode: designator for multiple pre-equilibrium model
c nlev      : number of excited levels for nucleus  
c
      do 10 type=0,6
        outtype(type)=' '
   10 continue
      maxZ=numZ-2
      maxN=numN-2
      nbins=40
      segment=1
      nlevmax=20
      nlevmaxres=10
      do 20 type=0,6
        if (type.le.2.or.type.eq.6) then
          nlevbin(type)=10
        else
          nlevbin(type)=5
        endif
   20 continue
      nlevbin(k0)=nlevmax
      Ltarget=0
      isomer=1.
      core=-1
      gammax=2
      transpower=5
      transeps=1.e-8
      xseps=1.e-7
      popeps=1.e-3
      Rfiseps=1.e-3
      eninclow=0.
      nangle=numang
      nanglecont=18
      maxenrec=numenrec
      ldmodel=1
      wmode=1
      preeqmode=2
      mpreeqmode=1
      do 30 Nix=0,numN
        do 30 Zix=0,numZ
          nlev(Zix,Nix)=0
   30 continue
c
c **************** Read second set of input variables ******************
c
c nlines: number of input lines
c ch    : character
c inline: input line                 
c parsym: symbol of particle
c Zinit : charge number of initial compound nucleus 
c Ninit : neutron number of initial compound nucleus 
c
c First the keyword is identified. Next the corresponding value is read.
c Erroneous input is immediately checked.
c
      do 110 i=1,nlines
        if (inline(i)(1:10).eq.'ejectiles ') then
          ip=-1
          do 120 i2=11,80
            ch=inline(i)(i2:i2)
            if (ch.ne.' ') then 
              do 130 type=0,6
                if (ch.eq.parsym(type)) then
                  ip=ip+1
                  if (ip.gt.6) goto 200
                  outtype(ip)=ch
                  goto 120
                endif
  130         continue
              goto 200
            endif
  120     continue
          goto 110
        endif
        if (inline(i)(1:5).eq.'maxz ') 
     +    read(inline(i)(6:80),*,err=200) maxZ
        if (inline(i)(1:5).eq.'maxn ') 
     +    read(inline(i)(6:80),*,err=200) maxN
        if (inline(i)(1:5).eq.'bins ') 
     +    read(inline(i)(6:80),*,err=200) nbins
        if (inline(i)(1:8).eq.'segment ') 
     +    read(inline(i)(9:80),*,err=200) segment
        if (inline(i)(1:13).eq.'maxlevelstar ') 
     +    read(inline(i)(14:80),*,err=200) nlevmax
        if (inline(i)(1:13).eq.'maxlevelsres ') 
     +    read(inline(i)(14:80),*,err=200) nlevmaxres
        if (inline(i)(1:8).eq.'ltarget ') 
     +    read(inline(i)(9:80),*,err=200) Ltarget
        if (inline(i)(1:7).eq.'isomer ') 
     +    read(inline(i)(8:80),*,err=200) isomer
        if (inline(i)(1:5).eq.'core ') 
     +    read(inline(i)(6:80),*,err=200) core
        if (inline(i)(1:7).eq.'gammax ') 
     +    read(inline(i)(8:80),*,err=200) gammax
        if (inline(i)(1:11).eq.'transpower ') 
     +    read(inline(i)(12:80),*,err=200) transpower
        if (inline(i)(1:9).eq.'transeps ') 
     +    read(inline(i)(10:80),*,err=200) transeps
        if (inline(i)(1:6).eq.'xseps ') 
     +    read(inline(i)(7:80),*,err=200) xseps
        if (inline(i)(1:7).eq.'popeps ') 
     +    read(inline(i)(8:80),*,err=200) popeps
        if (inline(i)(1:8).eq.'Rfiseps ') 
     +    read(inline(i)(9:80),*,err=200) Rfiseps
        if (inline(i)(1:5).eq.'elow ') 
     +    read(inline(i)(6:80),*,err=200) eninclow
        if (inline(i)(1:7).eq.'angles ') 
     +    read(inline(i)(8:80),*,err=200) nangle
        if (inline(i)(1:11).eq.'anglescont ') 
     +    read(inline(i)(12:80),*,err=200) nanglecont
        if (inline(i)(1:9).eq.'maxenrec ') 
     +    read(inline(i)(10:80),*,err=200) maxenrec
        if (inline(i)(1:8).eq.'ldmodel ') 
     +    read(inline(i)(9:80),*,err=200) ldmodel
        if (inline(i)(1:10).eq.'widthmode ') 
     +    read(inline(i)(11:80),*,err=200) wmode
        if (inline(i)(1:10).eq.'preeqmode ') 
     +    read(inline(i)(11:80),*,err=200) preeqmode
        if (inline(i)(1:11).eq.'mpreeqmode ') 
     +    read(inline(i)(12:80),*,err=200) mpreeqmode
        if (inline(i)(1:8).eq.'nlevels ') then
          read(inline(i)(9:80),*,err=200) iz,ia,ivalue
          Zix=Zinit-iz
          Nix=Ninit-ia+iz
          if (Zix.lt.0.or.Zix.gt.numZ.or.Nix.lt.0.or.Nix.gt.numN) then
            write(*,'("TALYS-warning: Z,N index out of range,",$)')
            write(*,'(" keyword ignored: ",a80)') inline(i)
          else
            nlev(Zix,Nix)=ivalue
          endif
        endif           
        if (inline(i)(1:13).eq.'maxlevelsbin ') then
          do 140 i2=14,80
            ch=inline(i)(i2:i2)
            if (ch.ne.' ') then       
              do 150 type=0,6
                if (ch.eq.parsym(type)) then
                  type2=type
                  goto 160
                endif
  150         continue
              goto 200
  160         read(inline(i)(i2+1:80),*,err=200) nlevbin(type2)
              goto 110
            endif
  140     continue
        endif           
  110 continue
      return
  200 write(*,'("TALYS-error: Wrong input: ",a80)') inline(i)
      stop                                                     
      end
Copyright (C) 2004  A.J. Koning, S. Hilaire and M.C. Duijvestijn
