      subroutine masses
c
c +---------------------------------------------------------------------
c | Author: Arjan Koning 
c | Date  : August 18, 2004
c | Task  : Nuclear masses 
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      logical          lexist
      character*4      masschar
      character*90     massfile
      integer          Zix,Z,Nbegin,Nend,Abegin,Aend,ia,N,Nix,A,type
      real             ebin
      double precision expmass,thmass,expmexc1,thmexc1
c
c ************************ Read nuclear masses *************************
c
c Zix     : charge number index for residual nucleus
c maxZ    : maximal number of protons away from initial compound nucleus
c Z       : charge number of residual nucleus
c Zinit   : charge number of initial compound nucleus
c Nbegin  : first N to be included
c Ninit   : neutron number of initial compound nucleus
c maxN    : maximal number of neutrons away from initial compound 
c           nucleus 
c Nend    : last N to be included
c Abegin  : first A to be included
c Aend    : last A to be included
c masschar: help variable
c massfile: mass file
c path    : directory containing structure files to be read
c lenpath : length of pathname
c ia      : mass number from mass table
c expmass : experimental mass
c thmass  : theoretical mass
c expmexc : experimental mass excess
c thmexc  : theoretical mass excess
c N       : neutron number of residual nucleus
c Nix     : neutron number index for residual nucleus
c nucmass : mass of nucleus 
c
c We read both the experimental masses, from Audi-Wapstra (1995), and 
c the theoretical masses, from Moller (1995), from the masstable.
c The experimental nuclear mass is adopted, when available. We also read
c the experimental and theoretical mass excess, to enable a more precise
c calculation of separation energies. If we need separation energies and
c specific masses for nuclides up to (maxZ,maxN), we need nuclear masses
c up to (maxZ+4,maxN+4). 
c
      do 10 Zix=0,maxZ+4
        Z=Zinit-Zix
        Nbegin=Ninit-maxN-4
        Nend=Ninit
        Abegin=Z+Nbegin
        Aend=Z+Nend
        masschar='z   '
        write(masschar(2:4),'(i3.3)') Z
        massfile=path(1:lenpath)//'masses/'//masschar
        inquire (file=massfile,exist=lexist)
        if (.not.lexist) goto 10
        open (unit=1,status='old',file=massfile)
   20   read(1,'(4x,i4,4f12.6)',end=30) ia,expmass,thmass,expmexc1,
     +    thmexc1
        if (ia.lt.Abegin) goto 20
        if (ia.gt.Aend) goto 30
        N=ia-Z
        Nix=Ninit-N
        if (expmass.ne.0.) then
          nucmass(Zix,Nix)=expmass
        else
          nucmass(Zix,Nix)=thmass
        endif
        expmexc(Zix,Nix)=expmexc1
        thmexc(Zix,Nix)=thmexc1
        goto 20
   30   close (unit=1)
   10 continue
c
c The target nucleus MUST be present in the masstable. This is to avoid
c unbound nuclei.
c
c parZ: charge number of particle
c parN: neutron number of particle
c k0  : index of incident particle
c
      if (nucmass(parZ(k0),parN(k0)).eq.0.) then
        write(*,'("TALYS-error: Target nucleus not in masstable")')
        stop                                                    
      endif
c
c ********* Use analytical mass formula for remaining nuclei ***********
c
c duflo  : Analytical mass formula of Duflo-Zuker 
c ebin   : total binding energy
c parmass: mass of particle in a.m.u.
c amu    : atomic mass unit in MeV  
c A      : mass number of residual nucleus
c
c If a residual nucleus is not in the experimental/theoretical mass 
c table, we use the analytical formula of Duflo-Zuker.
c
      do 110 Zix=0,maxZ+4
        Z=Zinit-Zix
        if (Z.le.0) goto 110
        do 120 Nix=0,maxN+4
          if (thmexc(Zix,Nix).eq.0.) then
            N=Ninit-Nix
            if (N.le.0) goto 110
            call duflo(N,Z,ebin)
            nucmass(Zix,Nix)=Z*parmass(2)+N*parmass(1)-ebin/amu
            A=Z+N
            thmexc(Zix,Nix)=(nucmass(Zix,Nix)-A)*amu
          endif
  120   continue
  110 continue
c
c ********************* Reduced and specific masses ********************
c
c specmass: specific mass for target nucleus
c redumass: reduced mass
c
      do 210 Zix=0,maxZ+2
        do 210 Nix=0,maxN+2
          do 220 type=0,6
            if (parskip(type)) goto 220
            specmass(Zix,Nix,type)=nucmass(Zix,Nix)/(nucmass(Zix,Nix)+
     +        parmass(type))
            redumass(Zix,Nix,type)=specmass(Zix,Nix,type)*parmass(type)
  220     continue
  210 continue
      return
      end
Copyright (C) 2004  A.J. Koning, S. Hilaire and M.C. Duijvestijn
