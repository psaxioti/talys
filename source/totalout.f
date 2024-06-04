      subroutine totalout
c
c +---------------------------------------------------------------------
c | Author: Arjan Koning 
c | Date  : September 1, 2004
c | Task  : Output of total cross sections
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      character*9 totfile
      integer     nen,i
c
c ********************* Total cross sections ***************************
c
c Einc        : incident energy in MeV
c nin         : counter for incident energy
c eninccm     : center-of-mass incident energy in MeV
c k0          : index of incident particle
c xstotinc    : total cross section (neutrons only) for incident channel
c xselasinc   : total elastic cross section (neutrons only) for incident
c               channel
c xsreacinc   : reaction cross section for incident channel
c xscompel    : compound elastic cross section
c xsnonel     : non-elastic cross section
c xsdirdiscsum: total direct cross section
c xspreeqsum  : total preequilibrium cross section summed over particles
c flaggiant   : flag for collective contribution from giant resonances
c xsgrsum     : sum over giant resonance cross sections
c xscompnonel : total compound non-elastic cross section
c xselastot   : total elastic cross section (shape + compound)
c
      write(*,'(/"########### REACTION SUMMARY FOR E=",f7.3,$)') Einc
      write(*,'(" ###########"/)') 
      write(*,'("Center-of-mass energy: ",f7.3/)') eninccm
      write(*,'("1. Total (binary) cross sections"/)') 
      if (k0.eq.1) write(*,'("Total           =",1p,e12.5)') xstotinc
      if (k0.eq.1) write(*,'("  Shape elastic   =",1p,e12.5)') xselasinc
      write(*,'("  Reaction        =",1p,e12.5)') xsreacinc
      write(*,'("    Compound elastic=",1p,e12.5)') xscompel
      write(*,'("    Non-elastic     =",1p,e12.5)') xsnonel
      write(*,'("      Direct          =",1p,e12.5)') xsdirdiscsum
      write(*,'("      Pre-equilibrium =",1p,e12.5)') xspreeqsum
      if (flaggiant) write(*,'("      Giant resonance =",1p,e12.5)') 
     +  xsgrsum
      write(*,'("      Compound non-el =",1p,e12.5)') xscompnonel
      if (k0.eq.1) write(*,'("    Total elastic   =",1p,e12.5)') 
     +  xselastot
c
c Write results to separate file
c
c filetotal: flag for total cross sections on separate file
c numinclow: number of incident energies below Elow 
c parsym   : symbol of particle
c k0       : index of incident particle
c Atarget  : mass number of target nucleus
c nuc      : symbol of nucleus
c Ztarget  : charge number of target nucleus                   
c numinc   : number of incident energies 
c
      if (filetotal) then
        totfile='total.tot'
        if (nin.eq.numinclow+1) then
          open (unit=1,status='unknown',file=totfile)
          write(1,'("# ",a1," + ",i3,a2," Total cross sections")')
     +      parsym(k0),Atarget,nuc(Ztarget)
          write(1,'("# ")')
          write(1,'("# ")')
          write(1,'("# # energies =",i3)') numinc
          write(1,'("#    E      Non-elastic  Elastic     Total",$)')
          write(1,'("     Comp. el.  Shape el.  Reaction",$)')
          write(1,'(" Comp. nonel   Direct   Pre-equil.")')
          do 10 nen=1,numinclow
            write(1,'(1p,e10.3,2x,9e11.4)') eninc(nen),(0.,i=1,9)
   10     continue
        else
          open (unit=1,status='old',file=totfile)
          do 20 nen=1,nin+4
            read(1,*) 
   20     continue
        endif
        write(1,'(1p,e10.3,2x,9e11.4)') Einc,xsnonel,xselastot,
     +    xstotinc,xscompel,xselasinc,xsreacinc,xscompnonel,
     +    xsdirdiscsum,xspreeqsum
        close (unit=1)
      endif
      return
      end
Copyright (C) 2004  A.J. Koning, S. Hilaire and M.C. Duijvestijn
