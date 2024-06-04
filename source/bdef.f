      subroutine bdef(amass,zchar,epscloc,coul,sym,b)
c
c +---------------------------------------------------------------------
c | Author: Marieke Duijvestijn
c | Date  : September 9, 2004
c | Task  : Computes the binding energy of a deformed nucleus with 
c |         eccentricity by the droplet model without shell corrections.
c +---------------------------------------------------------------------
c
c *************************** Comments *********************************
c
c This subroutine is based on the subroutine bdef originally developed 
c by U. Brosa.
c
c ****************** Declarations and common blocks ********************
c
      external fcoul,fsurf                                             
      real     b,sym,coul,epscloc,zchar,amass,a3w,fcoul,hsc,fsurf
c
c **********************************************************************
c
      a3w=amass**(1./3.)
      coul=0.7053/a3w*fcoul(epscloc)-1.153/amass
      hsc=15.4941*amass-17.9439*a3w**2*fsurf(epscloc)
      sym=1.7826/amass**2*hsc
      b=-hsc+sym*(amass-zchar-zchar)**2+coul*zchar**2
      return                                                            
      end                                                               
