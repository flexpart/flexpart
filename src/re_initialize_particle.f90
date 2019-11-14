! SPDX-FileCopyrightText: FLEXPART 1998-2019, see flexpart_license.txt
! SPDX-License-Identifier: GPL-3.0-or-later

subroutine re_initialize_particle(zp,ust,wst,h,sigmaw,wp,nrand,ol)
!                                      i   i  i   i  i    io  io    i 
!=============== CBL skewed vertical profiles and formulation of LHH 1996 with profile of w3 from lHB 2000                                       ======
!=============== LHH formulation has been modified to account for variable density profiles and backward in time or forward in time simulations  ======            
!=============== this routine re-initiaalize particle velocity if a numerical instability in the cbl scheme generated a NaN value                ======
!=============== the particle velocity is extracted from the updraft and downdraft distribution as required                                      ======
!=============== the re-initialization si not perfect see e.g. Cassiani et al(2015) BLM                                                          ====== 
!======================================================================================================================================================   
!======================================================================================================================================================   
  use par_mod, only:pi
  use com_mod, only:ldirect,rannumb

  implicit none


  real :: usurad2,usurad2p,C0,costluar4,eps 
  parameter  (usurad2=0.7071067812,usurad2p=0.3989422804,C0=2,costluar4=0.66667,eps=0.000001)

  integer idum,nrand
  real :: wp,zp,ust,wst,h,dens,ddens,sigmaw,dsigmawdz,tlw,dcas,dcas1 !,ran3,gasdev
  real :: w3,w2
  real ::  z, &    
       skew, &
       skew2, &
       radw2, &
       fluarw,fluarw2, &
       rluarw, &
       xluarw, &
       aluarw, &
       bluarw, &
       sigmawa, &
       sigmawb, &  
       ath, &
       bth, &
       wb,wa 
  real timedir
  real ol,transition

!---------------------------------------------------------------------------  
!timedir direction of time forward (1) or backward(-1)
  nrand=nrand+1
  dcas1=rannumb(nrand)
  timedir=ldirect 
  z=zp/h
  transition=1.

  if (-h/ol.lt.15) transition=((sin((((-h/ol)+10.)/10.)*pi)))/2.+0.5

  w2=sigmaw*sigmaw  
  w3=(((1.2*z*((1.-z)**(3./2.)))+eps)*wst**3)*transition 
  skew=w3/(w2**1.5)
  skew2=skew*skew
  radw2=sqrt(w2) !sigmaw

  fluarw=costluar4*skew**0.333333333333333
  fluarw2=fluarw*fluarw
  rluarw=(1.+fluarw2)**3.*skew2/((3.+fluarw2)**2.*fluarw2)  !-> r
  xluarw=rluarw**0.5 !(1.+fluarw2)**1.5*skew/((3.+fluarw2)*fluarw)    !----> r^1/2

  aluarw=0.5*(1.-xluarw/(4.+rluarw)**0.5)
  bluarw=1.-aluarw

  sigmawa=radw2*(bluarw/(aluarw*(1.+fluarw2)))**0.5
  sigmawb=radw2*(aluarw/(bluarw*(1.+fluarw2)))**0.5

  wa=(fluarw*sigmawa)
  wb=(fluarw*sigmawb)



  if ((sign(1.,wp)*timedir).gt.0) then !updraft	  
100 wp=(dcas1*sigmawa+wa)
    if (wp.lt.0)  then
      nrand=nrand+1
      dcas1=rannumb(nrand)
      goto 100
    end if
    wp=wp*timedir
  else if ((sign(1.,wp)*timedir).lt.0) then !downdraft	  
101 wp=(dcas1*sigmawb-wb)
    if (wp.gt.0)  then 
      nrand=nrand+1
      dcas1=rannumb(nrand)
      goto 101
    end if
    wp=wp*timedir
  end if

  return
end subroutine re_initialize_particle
