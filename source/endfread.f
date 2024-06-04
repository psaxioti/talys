      subroutine endfread
c
c +---------------------------------------------------------------------
c | Author: Arjan Koning
c | Date  : October 31, 2007
c | Task  : Read ECIS results for incident particle on ENDF-6 energy 
c |         grid
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      integer          nen,Z,A
      real             tripathi,e,enuc
      double precision xs
c
c ************ Read total, reaction and elastic cross section **********
c
c nen6     : total number of energies 
c e6,e     : energies of ENDF-6 energy grid in MeV
c k0       : index of incident particle
c coullimit: energy limit for charged particle OMP calculation
c xs       : help variable
c xstot6   : total cross section (neutrons only) for ENDF-6 file
c xsreac6  : reaction cross section for ENDF-6 file
c xsopt6   : optical model reaction cross section for ENDF-6 file
c xselas6  : total elastic cross section (neutrons only) for ENDF-6 file
c
      open (unit=3,status='unknown',file='ecis06.endfcs')
      do 10 nen=1,nen6
        e=real(e6(nen))
        if (k0.gt.1.and.e.lt.coullimit(k0)) goto 10
        read(3,'()')
        if (k0.eq.1) then
          read(3,*) xs
          xstot6(nen)=real(xs)
        endif
        read(3,*) xs
        xsreac6(nen)=real(xs)
        xsopt6(nen)=real(xs)
        if (k0.eq.1) then
          read(3,*) xs
          xselas6(nen)=real(xs)
        endif
   10 continue
      close (unit=3)
      open (unit=10,status='unknown',file='ecis06.endfin')
      close (unit=10,status=ecisstatus)
c
c ************ Normalization with semi-empirical results ***************
c  
c flagsys   : flag for reaction cross section from systematics      
c ZZ,Z      : charge number of residual nucleus
c AA,A      : mass number of residual nucleus 
c enuc      : incident energy in MeV per nucleon    
c parA      : mass number of particle     
c tripathi  : function for semi-empirical reaction cross section of
c             Tripathi et al.
c parZ      : charge number of particle
c threshnorm: normalization factor at trheshold
c
      if (flagsys(k0)) then
        Z=ZZ(0,0,k0)
        A=AA(0,0,k0)
        do 110 nen=1,nen6
          if (xsopt6(nen).eq.0.) goto 110
          e=real(e6(nen))
          enuc=e/parA(k0)
          xs=tripathi(parZ(k0),parA(k0),Z,A,enuc)
          if (xs.eq.0.) xs=xsopt(k0,nen)*threshnorm(k0)
          xsreac6(nen)=xs
          if (k0.eq.1) xselas6(nen)=xselas6(nen)+xsopt6(nen)-xs   
  110   continue
      endif
c
c **************** Write total cross sections to file ******************
c  
c parsym   : symbol of particle  
c Atarget  : mass number of target nucleus   
c nuc      : symbol of nucleus 
c Ztarget  : charge number of target nucleus     
c numinclow: number of incident energies below Elow
c
      open (unit=1,status='unknown',file='endf.tot')
      write(1,'("# ",a1," + ",i3,a2," Total cross sections")')
     +  parsym(k0),Atarget,nuc(Ztarget)
      write(1,'("# ")')
      write(1,'("# ")')
      write(1,'("# # energies =",i3)') nen6+numinclow
      write(1,'("#    E         Reaction  Elastic     Total")')
      do 210 nen=1,numinclow
        write(1,'(1p,e12.5,2x,3e11.4)') eninc(nen),fxsreacinc(nen),
     +    fxselasinc(nen),fxstotinc(nen)
  210 continue
      do 220 nen=1,nen6
        write(1,'(1p,e12.5,2x,3e11.4)') e6(nen),xsreac6(nen),
     +    xselas6(nen),xstot6(nen)
  220 continue
      close (unit=1)
      return
      end
Copyright (C) 2004  A.J. Koning, S. Hilaire and M.C. Duijvestijn
