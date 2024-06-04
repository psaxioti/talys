      subroutine compnorm
c
c +---------------------------------------------------------------------
c | Author: Arjan Koning 
c | Date  : December 8, 2006
c | Task  : Normalization of compound nucleus cross section
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      integer Zix,Nix,parity,J2,jj2beg,jj2end,jj2,lbeg,lend,l,updown
      real    xsreacsum,pik2,cfratio,norm
c
c There is a small difference between the reaction cross section as 
c calculated by ECIS and the sum over transmission coefficients.
c We therefore normalize the level population accordingly.
c
c ** Set (kinematical) factors for compound cross section calculation **
c
c xsreacsum  : reaction cross section constructed from transmission
c              coefficients
c pik2       : pi/(k**2) in mb
c pi         : pi
c wavenum    : wave number
c Zix        : charge number index for residual nucleus
c parZ       : charge number of particle
c k0         : index of incident particle
c Nix        : neutron number index for residual nucleus
c parN       : neutron number of particle
c CNfactor   : factor for compound nucleus cross section: 
c              pi/[ k**2 (2s+1)(2I+1) ]
c parspin    : spin of incident particle
c targetspin2: 2 * spin of target
c J2beg      : 2 * start of J summation
c targetspin : spin of target
c J2end      : 2 * end of J summation
c lmaxinc    : maximal l-value for transmission coefficients for
c              incident channel                                        
c numJ       : maximal J-value 
c flagwidth  : flag for width fluctuation calculation
c
      xsreacsum=0.
      pik2=10.*pi/(wavenum*wavenum)
      Zix=parZ(k0)
      Nix=parN(k0)
      CNfactor=pik2/(2.*parspin(k0)+1.)/(targetspin2+1.)
      if (k0.eq.0) CNfactor=0.5*CNfactor
      J2beg=mod(int(2.*(targetspin+parspin(k0))),2)
      J2end=int(2*(lmaxinc+parspin(k0)+targetspin))
      J2end=min(J2end,2*numJ)
c
c ******************* Loop over incoming channels **********************
c
c In order to get do loops running over integer values, certain quantum 
c numbers are multiplied by 2, which can be seen from a 2 present in the
c corresponding variable names. For each loop, the begin and end point 
c is determined from the triangular rule.
c
c 10: Sum over compound nucleus parity
c
c parity : parity
c pardif : difference between target and compound nucleus parity
c targetP: parity of target
c
      do 10 parity=-1,1,2
        pardif=abs(targetP-parity)/2
c
c 20: Sum over total angular momentum J (J2) of compound nucleus
c
c J2: 2 * J
c
        do 20 J2=J2beg,J2end,2
c
c 30: Sum over j (jj2) of incident channel
c
c jj2beg: 2 * start of j summation
c jj2end: 2 * end of j summation
c jj2   : 2 * j
c
          jj2beg=int(abs(J2-targetspin2))
          jj2end=int(J2+targetspin2)
          do 30 jj2=jj2beg,jj2end,2
c
c 40: Sum over l of incident channel
c
c lbeg: start of l summation
c lend: end of l summation
c
            lbeg=int(abs(0.5*jj2-parspin(k0)))
            lend=int(0.5*jj2+parspin(k0))
            do 40 l=lbeg,lend
c
c Check parity conservation and make index for transmission 
c coefficient. Sum transmission coefficients for incident channel.
c Add partial contribution to the reaction cross section.
c
c updown : spin index for transmission coefficient
c spin2  : 2 * spin of particle (usually)
c Tjlinc : transmission coefficients as a function of j and l 
c          for the incident channel
c
              if (l.gt.lmaxinc.or.mod(l,2).ne.pardif) goto 40
              updown=int(jj2-2.*l)/spin2(k0)
              xsreacsum=xsreacsum+CNfactor*(J2+1.)*Tjlinc(updown,l)
   40       continue
   30     continue
   20   continue
   10 continue
c
c Create the compound nucleus formation cross section xsflux and the 
c associated normalization factors.
c
c colltype    : type of collectivity (D, V or R)
c flagrot     : flag for use of rotational optical model per
c               outgoing particle, if available  
c norm        : normalization factor    
c xsflux      : cross section flux
c xsdirdiscsum: total direct cross section
c xscoupled   : inelastic cross section from coupled channels
c xspreeqsum  : total preequilibrium cross section summed over particles
c xsgrsum     : sum over giant resonance cross sections   
c xsreacinc   : reaction cross section for incident channel
c cfratio     : compound formation ratio
c flagcheck   : flag for output of numerical checks
c
c For coupled-channels calculations, the transmission coefficients are 
c already depleted by the contribution going into the strongly coupled
c inelastic channels.
c
      if (colltype(Zix,Nix).ne.'S'.and.flagrot(k0)) then
        norm=1.
        xsflux=xsreacsum-xsdirdiscsum+xscoupled-xspreeqsum-xsgrsum
        cfratio=xsflux/xsreacsum
      else
        norm=xsreacsum/xsreacinc
        xsflux=xsreacinc-xsdirdiscsum-xspreeqsum-xsgrsum
        cfratio=xsflux/xsreacinc
      endif
      CNfactor=CNfactor*cfratio/norm
      if (flagcheck) then
        write(*,'(/" ++++++++++ Normalization of reaction cross",
     +    " section ++++++++++"/)')
        write(*,'(" Reaction cross section          :",f11.5," (A)")') 
     +    xsreacinc
        write(*,'(" Sum over T(j,l)                 :",f11.5," (B)")') 
     +    xsreacsum
        write(*,'(" Compound nucleus formation c.s. :",f11.5," (C)")') 
     +    xsflux
        write(*,'(" Ratio C/B                       :",f11.5)') 
     +    cfratio
      endif
      return
      end
Copyright (C) 2004  A.J. Koning, S. Hilaire and M.C. Duijvestijn
