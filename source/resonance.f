        subroutine resonance
c
c +---------------------------------------------------------------------
c | Author: Arjan Koning
c | Date  : March 6, 2015
c | Task  : Reconstruction and broadening of resonance information
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      logical      lexist
      character*12 xsfile,rpfile,pfile
      character*14 resfile
      character*72 infile,outfile
      character*80 string,headstring(5)
      integer      MF,MT,i,nlin,j,k,NR,NPr,NP0,NT,Ntot,Nh,ifile
      real         x(3),y(3),E(numP),xs(numP),Et(numP),xst(numP)
c
c prepro   : PREPRO routines for making cross section pointwise
c Liso     : isomeric number of target
c path     : directory containing structure files to be read
c lenpath  : length of pathname
c reslib   : library with resonance parameters
c lenreslib: length of library name with resonance parameters
c Tres     : temperature for broadening low energy cross sections
c flaggroup: flag for output of low energy groupwise cross sections
c
      resfile='z000a000i0.res'
      write(resfile(2:4),'(i3.3)') Ztarget
      write(resfile(6:8),'(i3.3)') Atarget
      write(resfile(10:10),'(i1)') Liso
      infile=path(1:lenpath)//'resfiles/'//reslib(1:lenreslib)//'/'
     +  //resfile
      inquire (file=infile,exist=lexist)
      if (.not.lexist) then
        write(*,'(" TALYS-warning: ",a72," does not exist")') infile
        return
      endif
      if (flaggroup) then
        outfile=resfile(1:11)//'group'
      else
        if (Tres.gt.0.) then
          outfile=resfile(1:11)//'point'
        else
          outfile=resfile(1:11)//'point0'
        endif
      endif
      call prepro(infile,outfile,Tres,flaggroup)
      open (unit=1,status='unknown',file=outfile)
   10 read(1,'(a80)',err=100,end=100) string
      read(string(71:72),'(i2)') MF
      if (MF.ne.3) goto 10
      read(string(73:75),'(i3)',end=100,err=100) MT
      read(1,'(a80)',err=100,end=100) string
      read(string(45:55),'(i11)',err=100) NR
      read(string(56:66),'(i11)',err=100) NPr
      nlin=1+(NR-1)/3
      do 110 i=1,nlin
        read(1,'(a80)',err=100) string
  110 continue
      nlin=1+(NPr-1)/3
      k=0
      do 120 i=1,nlin
        read(1,'(6e11.6)',err=100) (x(j),y(j),j=1,3)
        do 130 j=1,3
          k=k+1
          E(k)=x(j)*1.e-6
          xs(k)=y(j)*1.e3
  130   continue
  120 continue
      NPr=NPr-1
      NP0=NPR
      do 140 k=NP0,1,-1
        if (xs(k).gt.0.) then
          NPR=k
          goto 150
        endif
  140 continue
  150 read(1,'(a80)',err=100,end=100) string
      if (MT.eq.1) xsfile='totalxs.tot'
      if (MT.eq.2) xsfile='elastic.tot'
      if (MT.eq.18) xsfile='fission.tot'
      rpfile='rp000000.tot'
      if (MT.eq.102) then 
        xsfile='xs000000.tot'
        write(rpfile(3:5),'(i3.3)') Ztarget
        write(rpfile(6:8),'(i3.3)') Atarget+1
      endif
      if (MT.eq.103)  then
        xsfile='xs010000.tot'
        write(rpfile(3:5),'(i3.3)') Ztarget-1
        write(rpfile(6:8),'(i3.3)') Atarget
      endif
      if (MT.eq.107) then
        xsfile='xs000001.tot'
        write(rpfile(3:5),'(i3.3)') Ztarget-2
        write(rpfile(6:8),'(i3.3)') Atarget-3
      endif
c
c Read TALYS output files
c
      do 210 ifile=1,2
        if (ifile.eq.1) then
          pfile=xsfile
        else
          if (MT.lt.102) goto 210
          rpfile=rpfile
        endif
        inquire (file=pfile,exist=lexist)
        if (.not.lexist) goto 10
        open (unit=2,status='unknown',file=pfile)
        do i=1,5
          read(2,'(a80)') headstring(i)
        enddo
        read(headstring(4)(15:80),*) NT
        do k=1,NT
          read(2,*) Et(k),xst(k)
        enddo
        close (2)
        Nh=NT+1
        do k=1,NT
          if (Et(k).gt.E(NPr)) then
            Nh=k
            goto 160
          endif
        enddo
  160   Ntot=NPr+NT-Nh+1
        write(headstring(4)(15:20),'(i6)') Ntot
c
c Write final output files
c
        open (unit=2,status='unknown',file=pfile)
        do i=1,5
          write(2,'(a80)') headstring(i)
        enddo
        do k=1,NPr
          write(2,'(1p,2e12.5)') E(k),xs(k)
        enddo
        do k=Nh,NT
          write(2,'(1p,2e12.5)') Et(k),xst(k)
        enddo
        close(2)
  210 continue
      goto 10
  100 close(1)
      return
      end
Copyright (C)  2013 A.J. Koning, S. Hilaire and S. Goriely
