      function spincut(Zix,Nix,ald,Eex,ibar)
c
c +---------------------------------------------------------------------
c | Author: Arjan Koning
c | Date  : October 14, 2004
c | Task  : Spin cutoff factor
c +---------------------------------------------------------------------
c
c ******************* Declarations and common blocks *******************
c
      include "talys.cmb"
      integer Zix,Nix,ibar,A,Ns2d
      real    spincut,ald,Eex,scutconst,b2,aldm,ignatyuk,Umatch,s2m,s2d,
     +        Es2m,Es2d,ascutoff,bscutoff,U
c
c *********************** Spin cutoff parameter ************************
c
c spincut     : spin cutoff factor
c Zix         : charge number index for residual nucleus
c Nix         : neutron number index for residual nucleus
c ald         : level density parameter
c Eex         : excitation energy 
c ibar        : fission barrier
c
c We interpolate between the discrete level spin cutoff factor and 
c the value at the matching energy.
c
c 1. Below the matching energy
c
c AA,A        : mass number of nucleus
c ldmodel     : level density model
c scutconst   : constant for spin cutoff factor
c Rspincut    : adjustable constant for spin cutoff factor 
c b2,beta2    : deformation parameter
c Exmatch     : matching point for Ex
c Umatch      : Exmatch - pairing       
c aldm        : level density parameter at matching energy
c ignatyuk    : function for energy dependent level density parameter a 
c s2m,Es2d,...: help variables
c twothird    : 2/3
c pair        : pairing energy
c
      A=AA(Zix,Nix,0)
      if (ldmodel.eq.1) then
        scutconst=0.0888*Rspincut
      else
        scutconst=0.01389*Rspincut
      endif
      b2=beta2(Zix,Nix,ibar)
      if (Eex.le.Exmatch(Zix,Nix,ibar)) then
        aldm=ignatyuk(Zix,Nix,Exmatch(Zix,Nix,ibar),ibar)
        Umatch=Exmatch(Zix,Nix,ibar)-pair(Zix,Nix)
        if (Umatch.gt.0.) then
          if (ldmodel.eq.1) then
            s2m=scutconst*A**twothird*sqrt(aldm*Umatch)
          else
            s2m=scutconst*A**twothird*A*sqrt(Umatch/aldm)
          endif
        else
c
c von Egidy et al systematics, Nuc. Phys. A481, 189 (1988)
c

          s2m=(0.98*(A**0.29))**2
        endif
        if (ibar.ne.0) s2m=(1.+b2/3.)*s2m
c
c Discrete spin cutoff-factor may not be zero (this is the case if e.g.
c only the 0+ ground state is known).
c
c scutoffdisc : spin cutoff factor for discrete level region 
c edis        : energy of level
c Ntop        : highest discrete level for temperature matching
c ascutoff,...: help variables
c efistrrot   : energy of rotational transition states
c
        if (scutoffdisc(Zix,Nix,ibar).eq.0.) 
     +    scutoffdisc(Zix,Nix,ibar)=s2m
        s2d=scutoffdisc(Zix,Nix,ibar)
        if (s2d.ge.s2m) s2d=s2m
        Es2m=Exmatch(Zix,Nix,ibar)
        Es2d=0.5*edis(Zix,Nix,max(Ntop(Zix,Nix,0),1))
        if (ibar.ne.0) then
          Ns2d=max(Ntop(Zix,Nix,ibar),1)
          Es2d=0.5*efistrrot(Zix,Nix,ibar,Ns2d)
        endif
        ascutoff=(s2m-s2d)/(Es2m-Es2d)
        bscutoff=s2d-ascutoff*Es2d
        spincut=ascutoff*Eex+bscutoff
        if (spincut.le.0) spincut=s2m
      else
c
c 2. Above the matching energy
c
c U: excitation energy minus pairing energy
c
        U=Eex-pair(Zix,Nix)
        if (U.gt.0.) then
          if (ldmodel.eq.1) then
            spincut=scutconst*A**twothird*sqrt(ald*U)
          else
            spincut=scutconst*A**twothird*A*sqrt(U/ald)
          endif
        else
          spincut=(0.98*(A**0.29))**2
        endif
        if (ibar.ne.0) spincut=(1.+b2/3.)*spincut
      endif
      return
      end
Copyright (C) 2004  A.J. Koning, S. Hilaire and M.C. Duijvestijn
