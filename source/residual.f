      subroutine residual
c
c +---------------------------------------------------------------------
c | Author: Arjan Koning 
c | Date  : September 2, 2004
c | Task  : Residual production cross sections
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      integer Acomp,Zcomp,Ncomp,nex
c
c ************************ Cross sections ******************************
c
c maxA      : maximal number of nucleons away from the initial compound
c             nucleus
c maxZ      : maximal number of protons away from the initial compound 
c             nucleus
c maxN      : maximal number of neutrons away from the initial compound
c             nucleus
c xsresprod : total residual production (= reaction) cross section 
c Acomp     : mass number index for compound nucleus
c xsmassprod: residual production cross section per mass unit
c Zcomp     : charge number index for compound nucleus
c Ncomp     : neutron number index for compound nucleus
c xspopnuc  : population cross section per nucleus
c Nlast     : last discrete level
c tau       : lifetime of state in seconds
c xsbranch  : branching ratio for isomeric cross section
c xspopex   : population cross section summed over spin and parity
c
      maxA=maxZ+maxN
      xsresprod=0.
      do 10 Acomp=0,maxA
        xsmassprod(Acomp)=0.
        do 20 Zcomp=0,maxZ
          Ncomp=Acomp-Zcomp
          if (Ncomp.lt.0.or.Ncomp.gt.maxN) goto 20
          if (xspopnuc(Zcomp,Ncomp).eq.0.) goto 40
          xsresprod=xsresprod+xspopnuc(Zcomp,Ncomp)
          xsmassprod(Acomp)=xsmassprod(Acomp)+xspopnuc(Zcomp,Ncomp)
          do 30 nex=0,Nlast(Zcomp,Ncomp,0)
            if (nex.eq.0.or.tau(Zcomp,Ncomp,nex).ne.0.) then
              xsbranch(Zcomp,Ncomp,nex)=xspopex(Zcomp,Ncomp,nex)/
     +          xspopnuc(Zcomp,Ncomp)
            endif
   30     continue
c
c For non-threshold reactions (positive Q-value) we always assign a
c minimum value to the exclusive cross section. (The transmission
c coefficients for these reactions might have been zero (from ECIS),
c but non-threshold reactions theoretically have a non-zero cross
c section.)
c
c Qres : Q-value for residual nucleus
c xseps: limit for cross sections
c
   40     if (Qres(Zcomp,Ncomp,0).gt.0..and.
     +      xspopnuc(Zcomp,Ncomp).le.xseps) xspopnuc(Zcomp,Ncomp)=xseps
   20   continue
   10 continue
      return
      end
Copyright (C) 2004  A.J. Koning, S. Hilaire and M.C. Duijvestijn
