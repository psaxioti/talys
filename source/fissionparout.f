      subroutine fissionparout(Zix,Nix)
c
c +---------------------------------------------------------------------
c | Author: Arjan Koning    
c | Date  : September 8, 2004
c | Task  : Output for fission parameters
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      integer Zix,Nix,Z,N,A,i,j
c
c ****************** Output of fission barrier parameters **************
c
c Zix       : charge number index for residual nucleus
c Nix       : neutron number index for residual nucleus
c ZZ        : charge number of residual nucleus
c NN        : neutron number of residual nucleus
c AA        : mass number of residual nucleus
c nuc       : symbol of nucleus
c nfisbar   : number of fission barrier parameters
c flagclass2: flag for class2 states in fission
c nclass2   : number of sets of class2 states    
c axtype    : type of axiality of barrier (1: axial, 2: tri-axial)
c fbarrier  : height of fission barrier
c fwidth    : width of fission barrier
c Rtransmom : normalization constant for moment of inertia for
c             transition states    
c minertia  : moment of intertia of fission barrier deformation
c nfistrhb  : number of head band transition states for barrier
c fecont    : start of continuum energy
c
c 1. Main fission parameters
c
      Z=ZZ(Zix,Nix,0)
      N=NN(Zix,Nix,0)
      A=AA(Zix,Nix,0)  
      write(*,'(/"Fission information for Z=",i3," N=",i3,$)') Z,N
      write(*,'(" (",i3,a2,") "/)') A,nuc(Z)
      write(*,'("Number of fission barriers           :",i3)') 
     +  nfisbar(Zix,Nix)
      if (flagclass2) 
     +  write(*,'("Number of sets of class2 states      :",i3)') 
     +  nclass2(Zix,Nix)
      do 10 i=1,nfisbar(Zix,Nix)
        write(*,'(/"Parameters for fission barrier",i3/)') i
        write(*,'("Type of axiality                     :",i3,$)') 
     +    axtype(Zix,Nix,i)
        if (axtype(Zix,Nix,i).eq.1) then
          write(*,'(" (axial)")') 
        else
          write(*,'(" (tri-axial)")') 
        endif
        write(*,'("Height of fission barrier ",i1,"          :",f8.3)') 
     +    i,fbarrier(Zix,Nix,i)
        write(*,'("Width of fission barrier ",i1,"           :",f8.3)') 
     +    i,fwidth(Zix,Nix,i)
        write(*,'("Rtransmom                            :",f8.3)') 
     +    Rtransmom(Zix,Nix,i)
        write(*,'("Moment of inertia                    :",f8.3)') 
     +    minertia(Zix,Nix,i)
        write(*,'("Number of head band transition states:",i3)') 
     +    nfistrhb(Zix,Nix,i)
        write(*,'("Start of continuum energy            :",f8.3)') 
     +    fecont(Zix,Nix,i)
c
c 2. Head band transition states
c
c efistrhb: energy of head band transition states 
c jfistrhb: spin of head band transition states 
c cparity : parity of level (character)
c pfistrhb: parity of head band transition states 
c
        write(*,'(/"Head band transition states"/)') 
        write(*,'(" no.    E    spin    parity"/)') 
        do 20 j=1,nfistrhb(Zix,Nix,i)
          write(*,'(i4,f8.3,f6.1,3x,a1)') j,efistrhb(Zix,Nix,i,j),
     +      jfistrhb(Zix,Nix,i,j),cparity(pfistrhb(Zix,Nix,i,j))
   20   continue
c
c 3. Rotational bands transition states
c
c nfistrrot: number of rotational transition states for barrier
c efistrrot: energy of rotational transition states
c jfistrrot: spin of rotational transition states
c pfistrrot: parity of rotational transition states
c
        write(*,'(/"Rotational bands"/)') 
        write(*,'(" no.    E    spin    parity"/)') 
        do 30 j=1,nfistrrot(Zix,Nix,i)
          write(*,'(i4,f8.3,f6.1,3x,a1)') j,efistrrot(Zix,Nix,i,j),
     +      jfistrrot(Zix,Nix,i,j),cparity(pfistrrot(Zix,Nix,i,j))
   30   continue
   10 continue
c
c 4. Class2 states               
c
c Rclass2mom: normalization constant for moment of inertia for
c             class 2 states             
c minertc2  : moment of inertia for class2 states
c nfisc2hb  : number of class2 states for barrier
c widthc2   : width of class2 states 
c efisc2hb  : energy of class2 states 
c jfisc2hb  : spin of class2 states 
c pfisc2hb  : parity of class2 states 
c
      if (flagclass2) then
        do 40 i=1,nclass2(Zix,Nix)
          write(*,'(/"Parameters for set",i3," of class2 states"/)') i
          write(*,'("Rclass2mom                         :",f8.3)') 
     +      Rclass2mom(Zix,Nix,i)
          write(*,'("Moment of inertia                  :",f8.3)') 
     +      minertc2(Zix,Nix,i)
          write(*,'("Number of class2 states            :",i3)') 
     +      nfisc2hb(Zix,Nix,i)
          write(*,'("Width of class2 states (MeV)       :",f8.3)')
     +      widthc2(Zix,Nix,i)
          write(*,'(/"Class 2 states"/)') 
          write(*,'(" no.    E    spin    parity"/)') 
          do 50 j=1,nfisc2hb(Zix,Nix,i)
            write(*,'(i4,f8.3,f6.1,3x,a1)') j,efisc2hb(Zix,Nix,i,j),
     +        jfisc2hb(Zix,Nix,i,j),cparity(pfisc2hb(Zix,Nix,i,j))
   50     continue
c
c 5. Rotational bands
c
c nfisc2rot: number of class2 rotational transition states for barrier
c efisc2rot: energy of class2 rotational transition states
c jfisc2rot: spin of class2 rotational transition states
c pfisc2rot: parity of class2 rotational transition states
c
          write(*,'(/"Rotational bands"/)')
          write(*,'(" no.    E    spin    parity"/)')
          do 60 j=1,nfisc2rot(Zix,Nix,i)
            write(*,'(i4,f8.3,f6.1,3x,a1)') j,efisc2rot(Zix,Nix,i,j),
     +        jfisc2rot(Zix,Nix,i,j),cparity(pfisc2rot(Zix,Nix,i,j))
   60     continue
   40   continue
      endif
      return
      end  
Copyright (C) 2004  A.J. Koning, S. Hilaire and M.C. Duijvestijn
