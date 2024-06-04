      subroutine massdis
c
c +---------------------------------------------------------------------
c | Author: Marieke Duijvestijn and Arjan Koning
c | Date  : June 19, 2007
c | Task  : Fission fragment yields
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      character*90 yieldfile,fpfile
      real         fiseps,yielda(nummass),yieldacor(nummass),
     +             yieldaz(nummass,numelem),yieldazcor(nummass,numelem),
     +             fisepsA,partfisxs
      integer      k,i,Zcomp,Ncomp,Z,A,Zix,Nix,nexend,iskip,istep,nex,
     +             nen
c
c ************************** Mass yields *******************************
c
c fiseps    : limit for fission cross section per nucleus
c Rfiseps   : ratio for limit for fission cross section per nucleus
c xsfistot  : total fission cross section 
c nummass   : number of masses 
c yielda    : mass yield
c yieldacor : corrected mass yield
c numelem   : number of elements  
c yieldaz   : isotopic yield
c yieldazcor: corrected isotopic yield
c
c Initialization
c
      fiseps=Rfiseps*xsfistot
      do 10 k=1,nummass
        yielda(k)=0.
        yieldacor(k)=0.
        do 10 i=1,numelem
          yieldaz(k,i)=0.
          yieldazcor(k,i)=0.
   10 continue
c
c Loop over nuclides
c
c Zcomp      : charge number index for compound nucleus
c maxZ       : maximal number of protons away from the initial
c              compound nucleus
c Ncomp      : neutron number index for compound nucleus
c maxN       : maximal number of neutrons away from the initial
c              compound nucleus 
c ZZ,Z       : charge number of residual nucleus  
c AA,A       : mass number of residual nucleus  
c Zindex,Zix : charge number index for residual nucleus
c Nindex,Nix : neutron number index for residual nucleus
c maxex      : maximum excitation energy bin for compound nucleus
c fiseps     : limit for fission cross section per excitation energy bin
c iskip,istep: help variables
c xsbinary   : cross section from initial compound to residual nucleus
c Ex         : excitation energy
c partfisxs  : partial fission cross section
c fisfeedex  : fission contribution from excitation energy bin  
c brosafy    : subroutine for fission fragment yields based on Brosa 
c              model 
c disa       : normalised fission fragment mass yield per excitation 
c              energy bin 
c disacor    : normalised fission product mass yield per excitation 
c              energy bin 
c disaz      : normalised fission fragment isotope yield 
c              per excitation energy bin 
c disazcor   : normalised fission product isotope yield 
c              per excitation energy bin 
c
      do 20 Zcomp=0,maxZ
        do 20 Ncomp=0,maxN
c
c Brosa parameters available between Z=72 and Z=96, outside this range 
c we adopt the values of the boundary nuclides
c
          Z=ZZ(Zcomp,Ncomp,0)
          A=AA(Zcomp,Ncomp,0)
          Zix=Zindex(Zcomp,Ncomp,0)
          Nix=Nindex(Zcomp,Ncomp,0)
          if (xsfeed(Zcomp,Ncomp,-1).lt.fiseps) goto 20
          if (Zcomp.eq.0.and.Ncomp.eq.0) then
            nexend=maxex(Zcomp,Ncomp)+1
          else
            nexend=maxex(Zcomp,Ncomp)
          endif
          fisepsA=fiseps/max(3*maxex(Zcomp,Ncomp),1)
          iskip=0
          istep=4
          do 30 nex=nexend,0,-1
            if (nex.eq.maxex(Zcomp,Ncomp)+1) then
              excfis=Etotal
              partfisxs=xsbinary(-1)
            else
              if (mod(iskip,istep).ne.0) then
                iskip=iskip+1
                goto 30
              endif
              if (nex-istep+1.lt.0) goto 30
              if (Ex(Zcomp,Ncomp,nex-istep+1).ge.30.) then
                partfisxs=0.
                do 40 i=0,istep-1
                  partfisxs=partfisxs+fisfeedex(Zcomp,Ncomp,nex-i)
   40           continue
                if (partfisxs.ne.0) then
                  excfis=0.
                  do 50 i=0,istep-1
                    excfis=excfis+fisfeedex(Zcomp,Ncomp,nex-i)*
     +                Ex(Zcomp,Ncomp,nex-i)
   50             continue
                  excfis=excfis/partfisxs
                endif
                iskip=1
              else
                excfis=Ex(Zcomp,Ncomp,nex)
                partfisxs=fisfeedex(Zcomp,Ncomp,nex)
              endif
            endif
            if (partfisxs.gt.fisepsA) then
              call brosafy(Zix,Nix)
              do 60 k=1,A
                yielda(k)=yielda(k)+disa(k)*partfisxs
                yieldacor(k)=yieldacor(k)+disacor(k)*partfisxs
                do 60 i=1,Z
                  yieldaz(k,i)=yieldaz(k,i)+disaz(k,i)*partfisxs
                  yieldazcor(k,i)=yieldazcor(k,i)+disazcor(k,i)
     +              *partfisxs
   60           continue
            endif
   30     continue
   20 continue
c
c Write results to separate files
c
c yieldfile: file with fission yields
c natstring: string extension for file names
c iso      : counter for isotope
c Einc     : incident energy in MeV
c parsym   : symbol of particle
c k0       : index of incident particle
c Atarget  : mass number of target nucleus
c nuc      : symbol of nucleus
c Ztarget  : charge number of target nucleus    
c
      yieldfile='yield000.000.fis'//natstring(iso)
      write(yieldfile(6:12),'(f7.3)') Einc
      write(yieldfile(6:8),'(i3.3)') int(Einc)
      open (unit=1,status='unknown',file=yieldfile)
      write(1,'("# ",a1," + ",i3,a2,": mass yields")')
     +  parsym(k0),Atarget,nuc(Ztarget)
      write(1,'("# E-incident = ",f7.3)') Einc
      write(1,'("# ")')
      write(1,'("# ")') 
      write(1,'("# Mass    Yield   Corrected yield")')
      do 70 k=1,Atarget
        write(1,'(i3,3x,1p,e12.4,3x,e12.4)') k,yielda(k),yieldacor(k)
 70   continue
      close (unit=1)
c
c Write ff/fp residual production
c
c fpexist   : flag for existence of fission product
c fpfile    : file with fission product
c numinc    : number of incident energies   
c flagffevap: flag for calculation of particle evaporation from
c             fission fragment mass yields    
c eninc     : incident energy in MeV
c
      do 80 k=1,Atarget
        do 90 i=1,Ztarget
          if (yieldaz(k,i).lt.1.e-3.and..not.fpexist(i,k)) goto 90
          fpfile='fp000000.tot'//natstring(iso)
          write(fpfile(3:8),'(2i3.3)') i,k
          if (.not.fpexist(i,k)) then
            fpexist(i,k)=.true.
            open (unit=1,status='unknown',file=fpfile)
            write(1,'("# ",a1," + ",i3,a2,": ff yield of ",i3,a2)')
     +        parsym(k0),Atarget,nuc(Ztarget),k,nuc(i)
            write(1,'("# ")')
            write(1,'("# # energies =",i3)') numinc  
            write(1,'("# ")')
            if (flagffevap) then
              write(1,'("# E-incident   FF Yield   FP yield")')
              do 100 nen=1,nin-1
                write(1,'(1p,e10.3,e12.4,3x,e12.4)') eninc(nen),0.,0.
  100         continue      
            else
              write(1,'("# E-incident   FF Yield")')
              do 110 nen=1,nin-1
                write(1,'(1p,e10.3,e12.4)') eninc(nen),0.
  110         continue      
            endif
          else
            open (unit=1,status='old',file=fpfile)
            do 120 nen=1,nin+4
              read(1,*,end=130,err=130)
  120       continue       
          endif
          if (flagffevap) then
            write(1,'(1p,e10.3,e12.4,3x,e12.4)')
     +        Einc,yieldaz(k,i),yieldazcor(k,i)
          else
            write(1,'(1p,e10.3,e12.4)') Einc,yieldaz(k,i)
          endif
  130     close (unit=1)
   90   continue
   80 continue              
      end
Copyright (C) 2004  A.J. Koning, S. Hilaire and M.C. Duijvestijn
