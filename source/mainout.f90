      subroutine mainout
c
c +---------------------------------------------------------------------
c | Author: Arjan Koning
c | Date  : December 29, 2023
c | Task  : Main output
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      integer Zcomp,Ncomp,Zix,Nix,type,i,J,parity
c
c *************************** Code and version *************************
c
      write(*,'(/"    TALYS-1.97 (Version: December 29, 2023)"/)')
      write(*,'(" Copyright (C) 2023  A.J. Koning, S. Hilaire ",
     +  "and S. Goriely"/)')
      write(*,'(" Dimensions - Cross sections: mb, Energies: MeV, ",
     +  "Angles: degrees")')
c
c ***************** Write input file and default parameters ************
c
c inputout: subroutine to write the input parameters
c
      call inputout
c
c ********************** Main nuclear parameters ***********************
c
c Zcomp      : charge number index for compound nucleus
c Ncomp      : neutron number index for compound nucleus
c Zindex     : charge number index for residual nucleus
c Nindex     : neutron number index for residual nucleus
c k0         : index of incident particle
c parname    : name of particle
c parmass    : mass of particle in a.m.u.
c Atarget    : mass number of target nucleus
c tarmass    : mass of target nucleus
c Ltarget    : excited level of target
c edis       : energy of level
c jdis       : spin of level
c cparity    : parity of level (character)
c parlev     : parity of level
c tau        : lifetime of state in seconds
c parskip    : logical to skip outgoing particle
c numinc     : number of incident energies
c flaginitpop: flag for initial population distribution
c eninc      : incident energy in MeV
c parsym     : symbol of particle
c Q          : Q-value for target nucleus
c flagomponly: flag to execute ONLY an optical model calculation
c flagcomp   : flag for compound nucleus calculation
c
      Zcomp=0
      Ncomp=0
      Zix=Zindex(Zcomp,Ncomp,k0)
      Nix=Nindex(Zcomp,Ncomp,k0)
      write(*,'(/" ########## BASIC REACTION PARAMETERS ##########"/)')
      write(*,'(" Projectile           : ",a8,t37,
     +  "Mass in a.m.u.      : ",f10.6)') parname(k0),parmass(k0)
      write(*,'(" Target               : ",a,t37,
     +  "Mass in a.m.u.      : ",f10.6)') trim(targetnuclide),tarmass
      if (Ltarget.ne.0) then
        write(*,'(/" Excited target level : Number  Energy  ",
     +    "Spin Parity Lifetime(sec)")')
        write(*,'(24x,i3,4x,f7.4,2x,f4.1,3x,a1,4x,es10.3)')
     +    Ltarget,edis(Zix,Nix,Ltarget),jdis(Zix,Nix,Ltarget),
     +    cparity(parlev(Zix,Nix,Ltarget)),tau(Zix,Nix,Ltarget)
      endif
      write(*,'(/" Included channels:")')
      do type=-1,6
        if (parskip(type)) cycle
        write(*,'(21x,a8)') parname(type)
      enddo
      if (flagomponly.and..not.flagcomp) return
c
c Projectile
c
      if (.not.flaginitpop) then
        if (numinc.eq.1) then
          write(*,'(/,"     1 incident energy (LAB):"/)')
        else
          write(*,'(/,i6," incident energies (LAB):"/)') numinc
        endif
        do i=1,numinc
          if (eninc(i).lt.0.001) then
            write(*,'(1x,es10.3)') eninc(i)
          else
            write(*,'(1x,f10.3)') eninc(i)
          endif
        enddo
      else
c
c Initial population distribution
c
c npopE   : number of energies for population distribution
c npopJ   : number of spins for population distribution
c Exdist  : excitation energy of population distribution
c Pdistex : population distribution, spin-independent
c Pdist   : population distribution per spin and parity
c

        write(*,'(/," Initial population distribution - Bins: ",i4,
     +    " Spins: ",i3," Maximum excitation energy:",f12.5/)')
     +    npopE,npopJ,eninc(1)
        if (npopJ.eq.0) then
          write(*,'("    Ex     Population "/)')
          do i=1,npopE
            write(*,'(2es10.3)') EdistE(i),PdistE(i)
          enddo
        else
          write(*,'(" Parity   Ex ",11("      J=",i2)/)')
     +      (J,J=0,10)
          do parity=-1,1,2
            do i=1,npopE
              write(*,'(i6,12es10.3)') parity,EdistE(i),
     +          (PdistJP(i,J,parity),J=0,10)
            enddo
          enddo
        endif
      endif
      write(*,'(/" Q-values for binary reactions:"/)')
      do type=0,6
        if (parskip(type)) cycle
        write(*,'(" Q(",a1,",",a1,"):",f9.5)') parsym(k0),parsym(type),
     +    Q(type)
      enddo
      if (flagcheck) call arraysize
c
c * Write nuclear structure parameters for target and compound nucleus *
c
c flaglevels   : flag for output of discrete level information
c levelsout    : subroutine for output of discrete levels
c flagdensity  : flag for output of level densities
c filedensity  : flag for level densities on separate files
c densityout   : subroutine for output of level density parameters
c flagfisout   : flag for output of fission information
c fissionparout: subroutine for output for fission parameters
c strucwrite   : flag for output of nuclear structure info
c
      if (flaglevels) call levelsout(Zix,Nix)
      if (flagdensity) call densityout(Zix,Nix)
      if (flagfisout) call fissionparout(Zix,Nix)
      strucwrite(Zix,Nix)=.true.
      if (parskip(0)) return
      if (k0.ne.0) then
        if (flaglevels) call levelsout(Zcomp,Ncomp)
        if (flagdensity.or.filedensity) call densityout(Zcomp,Ncomp)
        if (flagfisout) call fissionparout(Zcomp,Ncomp)
        strucwrite(Zcomp,Ncomp)=.true.
      endif
      return
      end
Copyright (C)  2023 A.J. Koning, S. Hilaire and S. Goriely
