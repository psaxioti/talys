      subroutine preeqcomplex
c
c +---------------------------------------------------------------------
c | Author: Arjan Koning and Vivian Demetriou
c | Date  : December 14, 2006
c | Task  : Pre-equilibrium complex particle emission
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      logical flagknock,flaginel,surfwell
      integer type,ndelta,ndeltapi,ndeltanu,ppi,hpi,pnu,hnu,A,i,j,type2,
     +        nen,k,l,parity
      real    proj2sp1,ejecmass,ejec2sp1,term1,Kap,Va,term2,term3,term4,
     +        base,term5,termps,V1well,surface,XNT,gsn,gsp,Ewell,
     +        termki,AKO,Ccl,phi,denom,Pn(6),gscomp(6),emax,dE,total,
     +        Eout,sigav,denomki(6),termk0,termin0,factor1,Eres,omegaNT,
     +        omegaph,phdens2,U,factor
c
c ************************** Kalbach model *****************************
c
c The complex particle emission model is described in 
c C. Kalbach, "Preequilibrium reactions with complex channels", 
c Phys. Rev. C71, 034606 (2005).
c
c projmass: mass of projectile
c parmass : mass of particle in a.m.u.
c k0      : index of incident particle
c proj2sp1: 2*spin +1 of projectile
c parspin : spin of particle 
c parskip : logical to skip outgoing particle
c
c Factors for pickup and stripping processes. For reactions involving
c only neutrons and protons this is thus not considered.
c 
      projmass=parmass(k0)
      proj2sp1=2.*parspin(k0)+1.
      do 10 type=0,6
        if (parskip(type)) goto 10
        if (k0.le.2.and.type.le.2) goto 10
c
c Calculation of terms independent of emission energy.
c
c ejecmass: mass of ejectile
c ejec2sp1: 2*spin +1 of ejectile
c term1,..: help variables
c Kap     : alpha enhancement factor
c Einc    : incident energy in MeV
c ndelta  : number of transferred particles
c parA    : mass number of particle                    
c ndeltapi: number of transferred protons
c parZ    : charge number of particle
c ndeltanu: number of transferred neutrons
c parN    : neutron number of particle
c Va      : potential drop
c
        ejecmass=parmass(type)
        ejec2sp1=2.*parspin(type)+1.
        term1=ejec2sp1*ejecmass/(proj2sp1*projmass)
        Kap=1.
        if ((k0.eq.1.or.k0.eq.2).and.type.eq.6) Kap=12.
        if (k0.eq.6.and.(type.eq.1.or.type.eq.2)) 
     +    Kap=12.-11.*max(Einc-20.,0.)/Einc
        ndelta=abs(parA(k0)-parA(type))
        ndeltapi=parZ(k0)-parZ(type)
        ndeltanu=parN(k0)-parN(type)
        Va=12.5*projmass
        term2=Kap/projmass*(projmass/(Einc+Va))**(2*ndelta)
c
c Initial configuration for pickup, stripping or t-h charge exchange
c
c ppi: proton particle number
c hpi: proton hole number
c pnu: neutron particle number
c hnu: neutron hole number         
c
        ppi=max(ndeltapi,0)
        hpi=max(-ndeltapi,0)
        pnu=max(ndeltanu,0)
        hnu=max(-ndeltanu,0)
c
c Further terms
c
c AA,A   : mass number of residual nucleus 
c base   : help variable
c Ztarget: charge number of target nucleus 
c Atarget: mass number of target nucleus 
c termps : term for pickup and stripping
c Cstrip : adjustable parameter for stripping/pick-up reactions
c
        A=AA(0,0,type)
        if (k0.eq.1) then
          term3=(5500./real(A))**ndelta
        else
          term3=(3800./real(A))**ndelta
        endif
        if (parA(k0).lt.parA(type)) term4=1./(80.*Einc)
        if (parA(k0).gt.parA(type)) term4=1./(580.*sqrt(Einc))
        if (parA(k0).eq.parA(type)) term4=1./(1160.*sqrt(Einc))
        base=real(2.*Ztarget)/real(Atarget)
        term5=base**(2*(parZ(k0)+2)*hpi+2*pnu)
        termps=Cstrip(type)*term1*term2*term3*term4*term5
c
c XNT function 
c
c V1well,Ewell: depth of potential well 
c XNT         : probability of exciting particle-hole pair
c gsn,gsp,gsa : single-particle level density parameter
c Ntarget     : neutron number of target nucleus 
c Kph         : constant for single-particle level density parameter
c               (g=A/Kph)          
c
        V1well=17.
        if (k0.eq.1) V1well=surface(1,Einc)
        if (k0.eq.2) V1well=surface(2,Einc)
        if (k0.eq.5.or.k0.eq.6) V1well=25.
        XNT=sqrt(Einc/projmass)*7./(V1well*Atarget*Atarget)*
     +    (pnu**2+ppi**2+hnu**2+1.5*hpi**2)
        gsn=Ntarget/Kph
        gsp=Ztarget/Kph
        if (ndeltapi.eq.0) then
          Ewell=V1well*base
        else
          Ewell=V1well
        endif
c
c Factors for knockout and inelastic processes.
c
c termki       : term for knockout and inelastic
c AKO          : Pauli correction factor
c flagknock....: help variables
c xsreacinc    : reaction cross section for incident channel    
c Ccl          : adjustment factor
c phi          : time fraction for alpha cluster
c denom,Pn     : help variables
c gscomp       : single-particle level density parameter
c
c Calculation of terms independent of emission energy.
c
        termki=0.
        AKO=0.
        flagknock=((k0.eq.1.or.k0.eq.2).and.type.eq.6)
        flaginel=(k0.eq.type)
        if (.not.(flagknock.or.flaginel)) goto 100
        Ccl=1./14.
        phi=0.08
        if (Ntarget.gt.116.and.Ntarget.lt.126) 
     +    phi=0.02+0.06*(126-Ntarget)/10.
        if (Ntarget.ge.126.and.Ntarget.lt.129) 
     +    phi=0.02+0.06*(Ntarget-126)/3.
        denom=Atarget-2.*phi*Ztarget+0.5*phi*Ztarget
        Pn(1)=(Ntarget-phi*Ztarget)/denom
        Pn(2)=(Ztarget-phi*Ztarget)/denom
        Pn(6)=0.5*phi*Ztarget/denom
        gscomp(1)=Atarget/13.
        gscomp(2)=Atarget/13.
        gscomp(3)=Atarget/52.
        gscomp(4)=Atarget/156.
        gscomp(5)=Atarget/156.
        gscomp(6)=Atarget/208.
        if (flagknock) 
     +    AKO=1./(2.*(gscomp(k0)**2))+1./(2.*(gscomp(6)**2))
c
c Denominator of cluster emission formula     
c
c emax      : maximal emission energy for particle channel
c eninccm   : center-of-mass incident energy in MeV
c Q         : Q-value for target nucleus
c dE,total  : help variables
c coulbar   : Coulomb barrier  
c ebegin    : first energy point of energy grid 
c eend      : last energy point of energy grid 
c egrid,Eout: energies of basic energy grid in MeV
c sigav     : average cross section for emission channel
c xsreac    : reaction cross section
c deltaE    : energy bin around outgoing energies
c
        do 20 type2=1,6
          denomki(type2)=0.
          if (parskip(type2)) goto 20
          emax=eninccm+Q(type2)
          dE=emax-coulbar(type2)
          if (dE.gt.2.) then
            total=0.
            do 30 nen=ebegin(type2),eend(type2)
              Eout=egrid(nen)
              if (Eout.lt.coulbar(type2)) goto 30
              total=total+xsreac(type2,nen)*deltaE(nen)
   30       continue
            sigav=total/dE
          else
            sigav=xsreac(type2,eend(type2))
          endif
          denomki(type2)=(2.*parspin(type2)+1.)*sigav*
     +      (emax+2.*coulbar(type2))*(emax-coulbar(type2))**2
   20   continue
c
c Knockout and inelastic terms
c
c termk0,termin0: help variable
c Cknock        : adjustable parameter for knockout reactions
c
        if (flagknock) then
          termk0=projmass*denomki(k0)*
     +      gscomp(k0)*gscomp(6)**2/(6.*gscomp(k0))+
     +      ejecmass*denomki(6)*
     +      gscomp(k0)*gscomp(6)**2/(6.*gscomp(6))
          if (termk0.ne.0.) then
            term1=Ccl*xsreacinc*ejec2sp1*ejecmass
            termki=Cknock(type)*term1*Pn(6)*gscomp(k0)*gscomp(6)/termk0
          endif
        else
          do 50 type2=1,6
            if (type2.ge.3.and.type2.le.5) goto 50
            termin0=projmass*denomki(k0)*
     +        gscomp(k0)*gscomp(type2)**2/(6.*gscomp(k0))+
     +        projmass*denomki(type2)*
     +        gscomp(k0)*gscomp(type2)**2/(6.*gscomp(type2))
            if (termin0.ne.0.) then
              term1=Ccl*xsreacinc*proj2sp1*projmass
              termki=termki+Cknock(type)*term1*Pn(type2)*gscomp(type2)*
     +          gscomp(type2)/termin0
            endif
   50     continue
        endif
c
c Calculation of pre-equilibrium complex particle emission.
c
c factor1: help variable
c Eres   : total energy of residual system
c Etotal : total energy of compound system (target + projectile)
c S      : separation energy per particle                       
c
  100   do 110 nen=ebegin(type),eend(type)
          Eout=egrid(nen)
          factor1=xsreac(type,nen)*Eout
          Eres=Etotal-S(0,0,type)-Eout
c
c Check if outgoing energy exceeds maximal possible energy
c    
          if (Eres.lt.0.) goto 110
c
c Stripping/pick-up terms that depend on emission energy.
c For reactions in which only one hole is left, e.g. (p,d),
c the finite well function leads to a discontinuity at the well 
c depth. Therefore, for this case the well depth is set equal to
c the Fermi energy.
c
c omegaNT  : state density function for pickup and stripping
c surfwell : flag for surface effects in finite well  
c omegaph  : particle-hole state density
c phdens2  : function for two-component particle-hole state density
c xspreeqps: preequilibrium cross section per particle type and
c            outgoing energy for pickup and stripping
c U        : help variable
c xspreeqki: preequilibrium cross section per particle type and
c            outgoing energy for knockout and inelastic
c
          omegaNT=0.
          surfwell=.false.
          do 120 i=0,3
            do 120 j=0,3-i
              omegaph=phdens2(ppi+i,hpi+i,pnu+j,hnu+j,gsp,gsn,Eres,
     +          Ewell,surfwell)
              omegaNT=omegaNT+XNT**(i+j)*omegaph
  120     continue
          do 130 i=0,ppi
            do 130 j=0,hpi
              do 130 k=0,pnu
                do 130 l=0,hnu
                  if (i+j+k+l.ne.0) omegaNT=omegaNT+phdens2(ppi-i,
     +              hpi-j,pnu-k,hnu-l,gsp,gsn,Eres,Ewell,surfwell) 
  130     continue
          xspreeqps(type,nen)=termps*omegaNT*factor1

c
c Knockout term that depends on emission energy.
c
          if (flagknock.or.flaginel) then
            U=max(Eres-AKO,0.)
            xspreeqki(type,nen)=termki*factor1*U
          endif
c
c If the pre-equilibrium spin distribution is chosen, we assume that the
c spin distribution for pickup, stripping and knockout is the same as in
c the exciton model.
c
c flagpespin: flag for pre-equilibrium spin distribution or compound
c             spin distribution for pre-equilibrium cross section
c xspreeq   : preequilibrium cross section per particle type and
c             outgoing energy       
c parity    : parity
c maxJph    : maximal spin for particle-hole states
c xspreeqJP : preequilibrium cross section per particle type,
c             outgoing energy, spin and parity                    
c
          if (flagpespin.and.xspreeq(type,nen).ne.0.) then
            do 140 parity=-1,1,2
              do 140 J=0,maxJph
                factor=xspreeqJP(type,nen,J,parity)/xspreeq(type,nen)
                xspreeqJP(type,nen,J,parity)=
     +            xspreeqJP(type,nen,J,parity)+
     +            factor*(xspreeqps(type,nen)+xspreeqki(type,nen))
  140       continue
          endif
          xspreeq(type,nen)=xspreeq(type,nen)+xspreeqps(type,nen)+
     +      xspreeqki(type,nen)
  110   continue
   10 continue
      return
      end
Copyright (C) 2004  A.J. Koning, S. Hilaire and M.C. Duijvestijn
