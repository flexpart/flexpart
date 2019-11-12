subroutine initialize_cbl_vel(idum,zp,ust,wst,h,sigmaw,wp, ol)
!                              i/o   i  i   i  i     i  o   i  

  use par_mod, only:pi
  use com_mod, only:ldirect
  use random_mod, only: gasdev, ran3

  implicit none

!===============================================================================
! CBL skewed vertical profiles and formulation of LHH 1996 with profile of w3
! from LHB 2000
! LHH formulation has been modified to account for variable density profiles and
! backward in time or forward in time simulations
! see Cassiani et al. BLM 2014 doi  for explanations and references
!===============================================================================


  real :: usurad2,usurad2p,C0,costluar4,eps 
  parameter  (usurad2=0.7071067812,usurad2p=0.3989422804,C0=2,costluar4=0.66667,eps=0.000001)

  integer idum
  real :: wp,zp,ust,wst,h,dens,ddens,sigmaw,dsigmawdz,tlw,dcas,dcas1!,ran3,gasdev
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
  real ol, transition

!---------------------------------------------------------------------------  
  timedir=ldirect !direction of time forward (1) or backward(-1)
  z=zp/h


  transition=1.
  if (-h/ol.lt.15) transition=((sin((((-h/ol)+10.)/10.)*pi)))/2.+0.5  !see also in cbl.f90

  w2=sigmaw*sigmaw
  w3=(((1.2*z*((1.-z)**(3./2.)))+eps)*wst**3) *transition

  skew=w3/(w2**1.5)
  skew2=skew*skew

  radw2=sqrt(w2) !sigmaw

  fluarw=costluar4*skew**0.333333333333333
  fluarw2=fluarw*fluarw
  rluarw=(1.+fluarw2)**3.*skew2/((3.+fluarw2)**2.*fluarw2)  !-> r
  xluarw=rluarw**0.5 !----> r^1/2

  aluarw=0.5*(1.-xluarw/(4.+rluarw)**0.5)
  bluarw=1.-aluarw

  sigmawa=radw2*(bluarw/(aluarw*(1.+fluarw2)))**0.5
  sigmawb=radw2*(aluarw/(bluarw*(1.+fluarw2)))**0.5

  wa=(fluarw*sigmawa)
  wb=(fluarw*sigmawb)

  dcas=ran3(idum) 

  if (dcas.le.aluarw) then
    dcas1=gasdev(idum)
    wp=timedir*(dcas1*sigmawa+wa)
  else
    dcas1=gasdev(idum)
    wp=timedir*(dcas1*sigmawb-wb)
  end if

  return
end subroutine initialize_cbl_vel

