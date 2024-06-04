      subroutine endfenergies
c
c +---------------------------------------------------------------------
c | Author: Arjan Koning
c | Date  : November 13, 2007
c | Task  : Energy grid for ENDF-6 file
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      integer          nen,type,Zix,Nix,nex,idc,nen2
      real             Eout,degrid,Eeps
      double precision emax,ee,e6tmp
c
c ************************ Basic ENDF-6 energy grid ********************
c
c Eout,e6,Eeps: energies of ENDF-6 energy grid in MeV
c eninclow    : minimal incident energy for nuclear model calculations
c degrid      : energy increment
c enincmax    : maximum incident energy
c 
c The basic ENDF-6 energy grid we use is:
c
c         0.001 -   0.01 MeV  : dE= 0.001 MeV
c         0.01  -   0.1 MeV   : dE= 0.01 MeV
c         0.1   -   1 MeV     : dE= 0.05 MeV
c         1     -   4 MeV     : dE= 0.1 MeV
c         4     -   7 MeV     : dE= 0.2 MeV
c         7     -  20 MeV     : dE= 0.5 MeV
c        20     - 100 MeV     : dE= 1.0 MeV
c       100     - 200 MeV     : dE= 2.0 MeV
c
c This grid ensures that the total, elastic and reaction cross section  
c are calculated on a sufficiently precise energy grid.
c
      e6(1)=eninclow
      Eout=0.
      degrid=0.001
      nen=1
   10 Eout=Eout+degrid
      if (Eout.gt.e6(1)) then
        nen=nen+1
        e6(nen)=Eout
      endif
      Eeps=Eout+1.e-4
      if (Eeps.gt.enincmax) goto 100
      if (Eeps.gt.0.01) degrid=0.01
      if (Eeps.gt.0.1) degrid=0.05
      if (Eeps.gt.1.) degrid=0.1
      if (Eeps.gt.4.) degrid=0.2
      if (Eeps.gt.8.) degrid=0.5
      if (Eeps.gt.20.) degrid=1. 
      if (Eeps.gt.100.) degrid=2. 
      if (Eeps.gt.200.) goto 100
      goto 10
c
c *************** Add partial thresholds to energy grid ****************
c
c k0        : index of incident particle
c emax      : help variable
c Zindex,Zix: charge number index for residual nucleus
c Nindex,Nix: neutron number index for residual nucleus
c Nlast     : last discrete level  
c Ltarget   : excited level of target 
c Ethresh   : threshold incident energy for residual nucleus
c idnum     : counter for exclusive channel 
c idchannel : identifier for exclusive channel
c Ethrexcl  : threshold incident energy for exclusive channel
c nen6      : total number of energies
c
  100 if (k0.eq.1) then
        emax=min(enincmax,20.)
        do 110 type=1,6
          Zix=Zindex(0,0,type)
          Nix=Nindex(0,0,type)
          do 120 nex=0,Nlast(Zix,Nix,0)
            if (type.eq.k0.and.nex.eq.Ltarget) goto 120
            ee=Ethresh(Zix,Nix,nex)
            if (ee.gt.eninclow.and.ee.le.emax) then
              nen=nen+1
              e6(nen)=ee
            endif
  120     continue
  110   continue
        do 130 idc=0,idnum
          if (idchannel(idc).eq.100000) goto 130
          if (idchannel(idc).eq.10000) goto 130
          if (idchannel(idc).eq.1000) goto 130
          if (idchannel(idc).eq.100) goto 130
          if (idchannel(idc).eq.10) goto 130
          if (idchannel(idc).eq.1) goto 130
          ee=Ethrexcl(idc,0)
          if (ee.gt.eninclow.and.ee.le.emax) then
            nen=nen+1
            e6(nen)=ee
          endif
  130   continue
      endif
      nen6=nen
c
c *************************** Sort energies ****************************
c
      do 210 nen=1,nen6
        do 210 nen2=nen,nen6
          if (e6(nen).le.e6(nen2)) goto 210
          e6tmp=e6(nen)
          e6(nen)=e6(nen2)
          e6(nen2)=e6tmp
  210 continue
      return
      end
Copyright (C) 2004  A.J. Koning, S. Hilaire and M.C. Duijvestijn
