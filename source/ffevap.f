      subroutine ffevap
c
c +---------------------------------------------------------------------
c | Author: Arjan Koning
c | Date  : November 1, 2015
c | Task  : Evaporation of fission fragments
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      integer iz,ia,Zix,Nix,nex,type,nen,npar,in,ih,it,id,ip,
     +        ident,idc,iaa,inn
      real    xsexcpart(0:numpar,numin),sumpost,sumpfns,Esumpfns,E,Eav,
     +        fiseps,maxwell,summax,dE,xsc,yZA,yA
c
c ********************** Loop over fission fragments *******************
c
c Do a full TALYS calculation for each fission fragment and incident
c energy
c
      flagffruns=.true.
      do type=0,6
        do nen=0,numen2
          pfns(type,nen)=0.
          maxpfns(type,nen)=0.
        enddo
      enddo
      fiseps=Rfiseps*xsfistot
      write(*,'(/" ########## Start of loop over fission fragments"/)')
      do 20 ia=1,Atarget0
        do 30 iz=1,Ztarget0
          in=ia-iz
          if (in.lt.1.or.in.gt.numneu) goto 30
          if (xsfpZApre(iz,in).lt.fiseps.and..not.fpexist(iz,in))
     +      goto 30
          Aff=ia
          Zff=iz
          call evaptalys
c
c Add fission product cross sections
c
          do 110 Zix=0,maxZ
            do 120 Nix=0,maxN
              xsfpZApost(iz,in)=xsfpZApost(iz,in)+xspopnuc(Zix,Nix)
              xsfpApost(ia)=xsfpApost(ia)+xspopnuc(Zix,Nix)
              do 130 nex=0,Nlast(Zix,Nix,0)
                if (nex.eq.0.or.tau(Zix,Nix,nex).ne.0.) 
     +            xsfpex(iz,in,nex)=xsfpex(iz,in,nex)+
     +            xspopex(Zix,Nix,nex)
  130         continue
  120       continue
  110     continue
c
c Add prompt fission particle and gamma production and spectra
c
          if (xsinitpop.gt.0.) then
            do 202 npar=1,numin
              do 204 type=0,6
                xsexcpart(type,npar)=0.
  204         continue
              do 206 iaa=0,numia
              do 206 ih=0,numih
              do 206 it=0,numit
              do 206 id=0,numid
              do 206 ip=0,numip
              do 206 inn=0,numin
                if (inn+ip+id+it+ih+iaa.ne.npar) goto 206
                ident=100000*inn+10000*ip+1000*id+100*it+10*ih+iaa
                do 208 idc=0,idnum
                  if (idchannel(idc).eq.ident) then
                    xsc=xschannel(idc)
                    if (inn.gt.0) xsexcpart(1,inn)=xsexcpart(1,inn)+xsc
                    if (ip.gt.0) xsexcpart(2,ip)=xsexcpart(2,ip)+xsc
                    if (id.gt.0) xsexcpart(3,id)=xsexcpart(3,id)+xsc
                    if (it.gt.0) xsexcpart(4,it)=xsexcpart(4,it)+xsc
                    if (ih.gt.0) xsexcpart(5,ih)=xsexcpart(5,ih)+xsc
                    if (iaa.gt.0) xsexcpart(6,iaa)=xsexcpart(6,iaa)+xsc
                  endif
  208           continue
  206         continue
  202       continue
            yA=yieldApre(ia)
            yZA=yieldZApre(iz,in)
            do 210 type=0,6
              nubar(type)=nubar(type)+yZA*multiplicity(type)
              if (yA.gt.0.) then
                nupre(type,ia)=nupre(type,ia)+yZA/yA*multiplicity(type)
                do 220 npar=1,numin
                  Pdisnu(type,npar)=Pdisnu(type,npar)+
     +              xsexcpart(type,npar)/xsinitpop
  220           continue
                do 230 nen=0,numen2
                  pfns(type,nen)=pfns(type,nen)+yZA*xssumout(type,nen)
  230           continue
              endif
  210       continue
          endif
   30   continue
   20 continue
      write(*,'(/" ########## End of loop over fission fragments"/)')
c
c Average energy and relation to Maxwellian
c
      do 290 type=0,6
        sumpfns=0.
        Esumpfns=0.
        do 292 nen=1,numen2
          dE=espec(type,nen)-espec(type,nen-1)
          sumpfns=sumpfns+pfns(type,nen)*dE
          Esumpfns=Esumpfns+espec(type,nen)*pfns(type,nen)*dE
          maxpfns(type,nen)=0.
  292   continue
        if (sumpfns.gt.0.) then
          Eavpfns(type)=Esumpfns/sumpfns
        else
          Eavpfns(type)=0.
        endif
        summax=0.
        do 294 nen=1,numen2
          dE=espec(type,nen)-espec(type,nen-1)
          E=espec(type,nen)
          Eav=Eavpfns(type)
          maxwell=sqrt(E)*exp(-E/Eav)
          summax=summax+maxwell*dE
  294   continue
        do 296 nen=1,numen2
          E=espec(type,nen)
          Eav=Eavpfns(type)
          maxwell=sqrt(E)*exp(-E/Eav)
          if (maxwell.gt.0..and.sumpfns.gt.0.) 
     +      maxpfns(type,nen)=pfns(type,nen)/maxwell*summax/sumpfns
  296   continue
  290 continue
      sumpost=0.
      do 310 ia=1,Atarget0
        sumpost=sumpost+xsfpApost(ia)
  310 continue
      sumpost=0.5*sumpost
      if (sumpost.gt.0.) then
        xsfptotpost=0.
        yieldtotpost=0.
        do 320 iz=1,Ztarget0
          do 330 ia=iz+1,Atarget0
            in=ia-iz
            if (in.gt.numneu) goto 330
            if (xsfpZApost(iz,in).eq.0.) goto 330
            yieldZApost(iz,in)=xsfpZApost(iz,in)/sumpost
            yieldApost(ia)=yieldApost(ia)+yieldZApost(iz,in)
            yieldnpost(in)=yieldnpost(in)+yieldZApost(iz,in)
            xsfptotpost=xsfptotpost+xsfpZApost(iz,in)
            yieldtotpost=yieldtotpost+yieldZApost(iz,in)
  330     continue
  320   continue
      endif
c
c Reset variables to those of original target.
c
c flagffruns: flag to denote that run is for fission fragment
c
      flagffruns=.false.
      call talysinput
      flagmain=.false.
      call talysinitial
      flagmain=.true.
      if (.not.flagompall) call basicxs(0,0)
      if (parinclude(0)) call gamma(0,0)
      if (enincmax.ge.epreeq.or.flagracap) then
        call preeqinit
        call excitoninit
      endif
      if (flagracap) call racapinit
      if (flagcomp) call compoundinit
      if (flagastro) call astroinit
c 
c Output
c
      call massdisout
      call nubarout
      return
      end
Copyright (C)  2013 A.J. Koning, S. Hilaire and S. Goriely
