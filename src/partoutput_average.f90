subroutine partoutput_average(itime)
  !                             i
  !*****************************************************************************
  !                                                                            *
  !     Dump all particle positions                                            *
  !                                                                            *
  !     Author: A. Stohl                                                       *
  !                                                                            *
  !     12 March 1999                                                          *
  !                                                                            *
  !*****************************************************************************
  !                                                                            *
  ! Variables:                                                                 *
  !                                                                            *
  !*****************************************************************************

  use par_mod
  use com_mod

  implicit none

  real(kind=dp) :: jul
  integer :: itime,i,j,jjjjmmdd,ihmmss
  integer :: ix,jy,ixp,jyp,indexh,m,il,ind,indz,indzp
  real :: xlon,ylat
  real :: dt1,dt2,dtt,ddx,ddy,rddx,rddy,p1,p2,p3,p4,dz1,dz2,dz
  real :: topo,hm(2),hmixi,pv1(2),pvprof(2),pvi,qv1(2),qvprof(2),qvi
  real :: tt1(2),ttprof(2),tti,rho1(2),rhoprof(2),rhoi
  real :: tr(2),tri,zlim
  character :: adate*8,atime*6

  integer(kind=2) :: ishort_xlon(maxpart),ishort_ylat(maxpart),ishort_z(maxpart)
  integer(kind=2) :: ishort_topo(maxpart),ishort_tro(maxpart),ishort_hmix(maxpart)
  integer(kind=2) :: ishort_pv(maxpart),ishort_rho(maxpart),ishort_qv(maxpart)
  integer(kind=2) :: ishort_tt(maxpart),ishort_uu(maxpart),ishort_vv(maxpart)
  integer(kind=2) :: ishort_energy(maxpart)


  ! Determine current calendar date, needed for the file name
  !**********************************************************

  jul=bdate+real(itime,kind=dp)/86400._dp
  call caldate(jul,jjjjmmdd,ihmmss)
  write(adate,'(i8.8)') jjjjmmdd
  write(atime,'(i6.6)') ihmmss


  ! Some variables needed for temporal interpolation
  !*************************************************

  dt1=real(itime-memtime(1))
  dt2=real(memtime(2)-itime)
  dtt=1./(dt1+dt2)

  ! Open output file and write the output
  !**************************************

    open(unitpartout_average,file=path(2)(1:length(2))//'partposit_average_'//adate// &
         atime,form='unformatted',access='direct',status='replace',recl=24)


  ! Write current time to file
  !***************************

!  write(unitpartout_average) itime,numpart
  do i=1,numpart

  ! Take only valid particles
  !**************************

    if (itra1(i).eq.itime) then
      part_av_topo(i)=part_av_topo(i)/float(npart_av(i))
      part_av_z(i)=part_av_z(i)/float(npart_av(i))
      part_av_pv(i)=part_av_pv(i)/float(npart_av(i))
      part_av_qv(i)=part_av_qv(i)/float(npart_av(i))
      part_av_tt(i)=part_av_tt(i)/float(npart_av(i))
      part_av_uu(i)=part_av_uu(i)/float(npart_av(i))
      part_av_vv(i)=part_av_vv(i)/float(npart_av(i))
      part_av_rho(i)=part_av_rho(i)/float(npart_av(i))
      part_av_tro(i)=part_av_tro(i)/float(npart_av(i))
      part_av_hmix(i)=part_av_hmix(i)/float(npart_av(i))
      part_av_energy(i)=part_av_energy(i)/float(npart_av(i))

      part_av_cartx(i)=part_av_cartx(i)/float(npart_av(i))
      part_av_carty(i)=part_av_carty(i)/float(npart_av(i))
      part_av_cartz(i)=part_av_cartz(i)/float(npart_av(i))

! Project Cartesian coordinates back onto Earth's surface
! #######################################################
      xlon=atan2(part_av_cartx(i),-1.*part_av_carty(i))
      ylat=atan2(part_av_cartz(i),sqrt(part_av_cartx(i)*part_av_cartx(i)+ &
           part_av_carty(i)*part_av_carty(i)))
      xlon=xlon/pi180
      ylat=ylat/pi180


  ! Convert all data to integer*2 variables (from -32768 to 32767) for compressed output
  !*************************************************************************************

      if (xlon.gt.180.) xlon=xlon-360.
      if (xlon.lt.-180.) xlon=xlon+360.
      ishort_xlon(i)=nint(xlon*180.)
      ishort_ylat(i)=nint(ylat*360.)

      zlim=(part_av_z(i)*2.-32000.)
      zlim=min(zlim,32766.)
      zlim=max(zlim,-32766.)
      ishort_z(i)=nint(zlim)

      zlim=(part_av_topo(i)*2.-32000.)
      zlim=min(zlim,32766.)
      zlim=max(zlim,-32766.)
      ishort_topo(i)=nint(zlim)

      zlim=(part_av_tro(i)*2.-32000.)
      zlim=min(zlim,32766.)
      zlim=max(zlim,-32766.)
      ishort_tro(i)=nint(zlim)

      zlim=(part_av_hmix(i)*2.-32000.)
      zlim=min(zlim,32766.)
      zlim=max(zlim,-32766.)
      ishort_hmix(i)=nint(zlim)

      zlim=(part_av_rho(i)*20000.-32000.)
      zlim=min(zlim,32766.)
      zlim=max(zlim,-32766.)
      ishort_rho(i)=nint(zlim)

      zlim=(part_av_qv(i)*1000000.-32000.)
      zlim=min(zlim,32766.)
      zlim=max(zlim,-32766.)
      ishort_qv(i)=nint(zlim)

      zlim=(part_av_pv(i)*100.)
      zlim=min(zlim,32766.)
      zlim=max(zlim,-32766.)
      ishort_pv(i)=nint(zlim)

      zlim=((part_av_tt(i)-273.15))*300.
      zlim=min(zlim,32766.)
      zlim=max(zlim,-32766.)
      ishort_tt(i)=nint(zlim)

      zlim=(part_av_uu(i)*200.)
      zlim=min(zlim,32766.)
      zlim=max(zlim,-32766.)
      ishort_uu(i)=nint(zlim)

      zlim=(part_av_vv(i)*200.)
      zlim=min(zlim,32766.)
      zlim=max(zlim,-32766.)
      ishort_vv(i)=nint(zlim)

      zlim=(part_av_energy(i)-300000.)/30.
      zlim=min(zlim,32766.)
      zlim=max(zlim,-32766.)
      ishort_energy(i)=nint(zlim)

  ! Turn on for readable test output
  !*********************************

!        write(119,*) itime,i,xlon,ylat,part_av_z(i),part_av_topo(i),part_av_tro(i), &
!        part_av_hmix(i),part_av_rho(i),part_av_qv(i),part_av_pv(i),part_av_tt(i), &
!        ishort_uu(i),ishort_vv(i)
    endif

! Re-initialize averages, and set number of steps to zero
    npart_av(i)=0
    part_av_topo(i)=0.
    part_av_z(i)=0.
    part_av_cartx(i)=0.
    part_av_carty(i)=0.
    part_av_cartz(i)=0.
    part_av_pv(i)=0.
    part_av_qv(i)=0.
    part_av_tt(i)=0.
    part_av_uu(i)=0.
    part_av_vv(i)=0.
    part_av_rho(i)=0.
    part_av_tro(i)=0.
    part_av_hmix(i)=0.
    part_av_energy(i)=0.
    
  end do

  ! Write the output
  !*****************
  do i=1,numpart
    if (itra1(i).eq.itime) then
      write(unitpartout_average,rec=i) ishort_xlon(i),ishort_ylat(i),ishort_z(i), &
           ishort_topo(i),ishort_tro(i),ishort_hmix(i),ishort_rho(i),ishort_qv(i),ishort_pv(i), &
           ishort_tt(i),ishort_uu(i),ishort_vv(i)  ! ,ishort_energy(i)
    endif
  enddo


  close(unitpartout_average)

end subroutine partoutput_average
