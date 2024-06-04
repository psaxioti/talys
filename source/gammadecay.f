      subroutine gammadecay(Zix,Nix)
c
c +---------------------------------------------------------------------
c | Author: Arjan Koning 
c | Date  : July 1, 2004
c | Task  : Scheme for discrete gamma decay
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      integer     numgamdisc
      parameter   (numgamdisc=numlev*numlev/2)
      character*7 decayfile
      integer     Zix,Nix,i,j,k,ng,Ngam(numlev),Z,A,type
      real        flux(0:numlev),egam(numlev,numgamdisc),
     +            br(numlev,numgamdisc),total,yieldg(numlev),egamtmp,
     +            brtmp
c
c ********************* Construct gamma decay scheme *******************
c
c Zix        : charge number index for residual nucleus
c Nix        : neutron number index for residual nucleus
c nlev       : number of excited levels for nucleus
c flux       : flux in level (normalized to 1)
c tau        : lifetime of state in seconds
c branchratio: gamma-ray branching ratio to level
c egam       : gamma energy
c edis       : energy of level
c br         : branching ratio multiplied by initial flux
c total      : help variable
c yieldg     : total discrete gamma yield per level
c Ngam       : number of gamma ray lines per level
c
      do 10 i=1,nlev(Zix,Nix)
        ng=0
        do 20 j=0,i-1
          flux(j)=0.
   20   continue
        flux(i)=1.
        total=0.
        do 30 j=i,1,-1
          if (tau(Zix,Nix,j).ne.0.) goto 30
          do 40 k=0,j-1
            if (flux(j).eq.0.) goto 40
            if (branchratio(Zix,Nix,j,k).eq.0.) goto 40
            ng=ng+1
            egam(i,ng)=edis(Zix,Nix,j)-edis(Zix,Nix,k)
            br(i,ng)=branchratio(Zix,Nix,j,k)*flux(j)
            flux(k)=flux(k)+br(i,ng)
            total=total+br(i,ng)
   40     continue
   30   continue
        Ngam(i)=ng
        yieldg(i)=total
        if (total.eq.0.) goto 10
c
c Normalize intensities to 1.
c
        do 50 j=1,ng
            br(i,j)=br(i,j)/total
   50   continue
c
c Sort discrete gamma ray lines in descending order
c
        do 60 j=1,ng
          do 70 k=1,j
            if (egam(i,j).gt.egam(i,k)) then
              egamtmp=egam(i,k)
              brtmp=br(i,k)
              egam(i,k)=egam(i,j)
              br(i,k)=br(i,j)
              egam(i,j)=egamtmp
              br(i,j)=brtmp
            endif
   70     continue
   60   continue
   10 continue
c
c Write decay information to file
c
c ZZ       : charge number of residual nucleus   
c AA       : mass number of residual nucleus   
c decayfile: decay file
c parsym   : symbol of particle
c nuc      : symbol of nucleus
c
      Z=ZZ(Zix,Nix,0)
      A=AA(Zix,Nix,0)
      type=2*Zix+Nix
      decayfile='decay. '     
      write(decayfile(7:7),'(a1)') parsym(type)
      open (unit=1,status='unknown',file=decayfile)   
      write(1,'("# ",i3,a2," Discrete gamma decay")') A,nuc(Z)
      write(1,'("# ")')
      write(1,'("# ")')
      write(1,'("# # levels   =",i3)') nlev(Zix,Nix)
      write(1,'("#   E        fraction              ")')    
      do 110 i=1,nlev(Zix,Nix)
        write(1,'(2i4,f10.6)') i,Ngam(i),yieldg(i)
        do 120 j=1,Ngam(i)
          write(1,'(f11.6,1p,e12.5)') egam(i,j),br(i,j)
  120   continue
  110 continue
      close (unit=1)
      return
      end
Copyright (C) 2004  A.J. Koning, S. Hilaire and M.C. Duijvestijn
