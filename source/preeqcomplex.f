      subroutine preeqcomplex
c
c +---------------------------------------------------------------------
c | Author: Arjan Koning
c | Date  : October 16, 2004
c | Task  : Pre-equilibrium complex particle emission
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      logical nalpha,alphan,surfwell
      integer type,ndelta,ndeltapi,ndeltanu,ppi,hpi,pnu,hnu,A,i,j,nen,k,
     +        l,parity
      real    term1,term2,term3,Kap,Va,term4,term5,term6,term7,base,
     +        termps,V1well,XNT(0:3,0:3),gsn,gsp,Ewell,termki,AKO,phi,
     +        denom,gsa,ga,gb,Pb,emaxa,emaxb,dE,total,Eout,sigava,
     +        sigavb,term8a,term8b,term8,factor1,Eres,omegaps,omegaph,
     +        phdens2,UKO,factor
c
c ************************** Kalbach model *****************************
c
c The complex particle emission model is described in the PRECO-2000
c manual by C. Kalbach, Chapter 4.1, March 2001.
c
c parskip : logical to skip outgoing particle
c k0      : index of incident particle
c term1,..: help variables
c parspin : spin of particle 
c parmass : mass of particle in a.m.u.
c Einc    : incident energy in MeV
c Kap     : alpha enhancement factor
c ndelta  : number of transferred particles
c parA    : mass number of particle                    
c ndeltapi: number of transferred protons
c parZ    : charge number of particle
c ndeltanu: number of transferred neutrons
c parN    : neutron number of particle
c Va      : potential drop
c ppi     : proton particle number
c hpi     : proton hole number
c pnu     : neutron particle number
c hnu     : neutron hole number         
c AA,A    : mass number of residual nucleus 
c
c Factors for pickup and stripping processes. For reactions involving
c only neutrons and protons this is thus not considered.
c 
      do 10 type=0,6
        if (parskip(type)) goto 10
        if (k0.le.2.and.type.le.2) goto 10
c
c Calculation of terms independent of emission energy.
c
        term1=(2.*parspin(type)+1.)*parmass(type)
        term2=(2.*parspin(k0)+1.)*parmass(k0)
        term3=parmass(k0)*Einc
        if ((k0.eq.6.and.type.ne.6).or.(k0.le.2.and.type.eq.6)) then
          Kap=12.
        else
          Kap=1.
        endif
        ndelta=abs(parA(k0)-parA(type))
        ndeltapi=abs(parZ(k0)-parZ(type))
        ndeltanu=abs(parN(k0)-parN(type))
        Va=12.5*parmass(k0)
        term4=Kap*(parmass(k0)/(Einc+Va))**(2*ndelta)
        ppi=0
        hpi=ndeltapi
        pnu=0
        hnu=ndeltanu
        A=AA(0,0,type)
c
c Revised constant for incident neutrons, C. Kalbach, private 
c communication, January 2002.
c
c base   : help variable
c Ztarget: charge number of target nucleus 
c Atarget: mass number of target nucleus 
c termps : term for pickup and stripping
c
        if (k0.eq.1) then
          term5=(6500./real(A))**ndelta
        else
          term5=(3800./real(A))**ndelta
        endif
        term6=0.0127
        base=real(2.*Ztarget)/real(Atarget)
c
c Purely empirical trick to prevent huge (n,alpha) cross sections
c for nuclides with Z=N (like ca40).
c
        base=min(base,0.83)
        term7=base**(6*hpi)
        termps=term1/term2/term3*term4*term5*term6*term7
c
c Revised function for XNT, C. Kalbach, private communication, 
c January 2002.
c
c V1well,Ewell: depth of potential well 
c XNT         : probability of exciting particle-hole pair
c gsn,gsp,gsa : single-particle level density parameter
c Ntarget     : neutron number of target nucleus 
c Kph         : constant for single-particle level density parameter
c               (g=A/Kph)          
c
        if (k0.eq.1) V1well=7.
        if (k0.ge.2.and.k0.le.4) V1well=17.
        if (k0.eq.5.or.k0.eq.6) V1well=25.
c
c Problems with Kalbach's parameterization above 100 MeV.
c An energy dependent well depth is implemented to prevent 
c divergence of the production cross sections.
c
        V1well=V1well+(38.-V1well)*Einc**4/(Einc**4+150.**4)
        do 20 i=0,3
          do 20 j=0,3-i
            XNT(i,j)=sqrt(Einc/parmass(k0))*7./
     +        (V1well*Atarget*Atarget)*(pnu**2+ppi**2+hnu**2+1.5*hpi**2)
   20   continue
        gsn=Ntarget/Kph
        gsp=Ztarget/Kph
        Ewell=12.5
c
c Factors for knockout and inelastic processes. This term only 
c contributes when alpha-particles are involved.
c
c termki       : term for knockout and inelastic
c AKO          : Pauli correction factor
c nalpha,alphan: help variables
c xsreacinc    : reaction cross section for incident channel    
c phi          : time fraction for alpha cluster
c denom,Pb,dE  : help variables
c ga,gb        : single-particle level density parameter
c eninccm      : center-of-mass incident energy in MeV
c emaxa,emaxb  : maximal emission energy for particle channel
c Q            : Q-value for target nucleus
c coulbar      : Coulomb barrier  
c ebegin       : first energy point of energy grid 
c eend         : last energy point of energy grid 
c egrid,Eout   : energies of basic energy grid in MeV
c sigava,sigavb: average cross section for emission channel
c xsreac       : reaction cross section
c deltaE       : energy bin around outgoing energies
c
c Calculation of terms independent of emission energy.
c
        termki=0.
        AKO=0.
        nalpha=((k0.eq.1.or.k0.eq.2).and.type.eq.6)
        alphan=((type.eq.1.or.type.eq.2).and.k0.eq.6)
        if (.not.(nalpha.or.alphan)) goto 100
        term1=xsreacinc/12.*(2.*parspin(type)+1.)*parmass(type)
        if (Ntarget.le.116) phi=0.08
        if (Ntarget.gt.116.and.Ntarget.lt.126) 
     +    phi=0.02+0.06*(126-Ntarget)/10.
        if (Ntarget.ge.126.and.Ntarget.lt.129) 
     +    phi=0.02+0.06*(Ntarget-126)/3.
        if (Ntarget.ge.129) phi=0.08
        denom=Atarget-2.*phi*Ztarget+0.5*phi*Ztarget
        gsa=Atarget/208.
        if (k0.eq.1.and.type.eq.6) then
          ga=gsn
          gb=gsa
          Pb=0.5*phi*Ztarget/denom
        endif
        if (k0.eq.6.and.type.eq.1) then
          ga=gsa
          gb=gsn
          Pb=(Ntarget-phi*Ztarget)/denom
        endif
        if (k0.eq.2.and.type.eq.6) then
          ga=gsp
          gb=gsa
          Pb=0.5*phi*Ztarget/denom
        endif
        if (k0.eq.6.and.type.eq.2) then
          ga=gsa
          gb=gsp
          Pb=(Ztarget-phi*Ztarget)/denom
        endif
        AKO=0.5/(ga*ga)+0.5/(gb*gb)
        term2=Pb*ga*gb
        emaxa=eninccm
        emaxb=eninccm+Q(type)
        dE=emaxa-coulbar(k0)
        if (dE.gt.2.) then
          total=0.
          do 30 nen=ebegin(k0),eend(k0)
            Eout=egrid(nen)
            if (Eout.lt.coulbar(k0)) goto 30
            total=total+xsreac(k0,nen)*deltaE(nen)
   30     continue
          sigava=total/dE
        else
          sigava=xsreac(k0,eend(k0))
        endif
        dE=emaxb-coulbar(type)
        if (dE.gt.2.) then
          total=0.
          do 40 nen=ebegin(type),eend(type)
            Eout=egrid(nen)
            if (Eout.lt.coulbar(type)) goto 40
            total=total+xsreac(type,nen)*deltaE(nen)
   40     continue
          sigavb=total/dE
        else
          sigavb=xsreac(type,eend(type))
        endif
        term8a=(2.*parspin(k0)+1.)*parmass(k0)*sigava*
     +    (emaxa+2.*coulbar(k0))*(emaxa-coulbar(k0))**2*ga*gb**2/(6.*ga)
        term8b=(2.*parspin(type)+1.)*parmass(type)*sigavb*
     +    (emaxb+2.*coulbar(type))*(emaxb-coulbar(type))**2*
     +    ga*gb**2/(6.*gb)
        term8=term8a+term8b
        if (term8.ne.0.) termki=term1*term2/term8
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
          if (Eres.lt.0.) goto 10
c
c Stripping/pick-up terms that depend on emission energy.
c For reactions in which only one hole is left, e.g. (p,d),
c the finite well function leads to a discontinuity at the well 
c depth. Therefore, for this case the well depth is set equal to
c the Fermi energy.
c
c omegaps  : state density function for pickup and stripping
c surfwell : flag for surface effects in finite well  
c omegaph  : particle-hole state density
c phdens2  : function for two-component particle-hole state density
c xspreeqps: preequilibrium cross section per particle type and
c            outgoing energy for pickup and stripping
c UKO      : help variable
c xspreeqki: preequilibrium cross section per particle type and
c            outgoing energy for knockout and inelastic
c
          omegaps=0.
          surfwell=.false.
          do 120 i=0,3
            do 120 j=0,3-i
              omegaph=phdens2(ppi+i,hpi+i,pnu+j,hnu+j,gsp,gsn,Eres,
     +          Ewell,surfwell)
              omegaps=omegaps+XNT(i,j)**(i+j)*omegaph
  120     continue
          do 130 i=0,ppi
            do 130 j=0,hpi
              do 130 k=0,pnu
                do 130 l=0,hnu
                  if (i+j+k+l.ne.0) omegaps=omegaps+phdens2(ppi-i,
     +              hpi-j,pnu-k,hnu-l,gsp,gsn,Eres,Ewell,surfwell) 
  130     continue
          xspreeqps(type,nen)=termps*omegaps*factor1
c
c Knockout term that depends on emission energy.
c
          UKO=max(Eres-AKO,0.)
          xspreeqki(type,nen)=termki*factor1*UKO
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
