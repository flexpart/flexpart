!**********************************************************************
! Copyright 1998,1999,2000,2001,2002,2005,2007,2008,2009,2010         *
! Andreas Stohl, Petra Seibert, A. Frank, Gerhard Wotawa,             *
! Caroline Forster, Sabine Eckhardt, John Burkhart, Harald Sodemann   *
!                                                                     *
! This file is part of FLEXPART.                                      *
!                                                                     *
! FLEXPART is free software: you can redistribute it and/or modify    *
! it under the terms of the GNU General Public License as published by*
! the Free Software Foundation, either version 3 of the License, or   *
! (at your option) any later version.                                 *
!                                                                     *
! FLEXPART is distributed in the hope that it will be useful,         *
! but WITHOUT ANY WARRANTY; without even the implied warranty of      *
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the       *
! GNU General Public License for more details.                        *
!                                                                     *
! You should have received a copy of the GNU General Public License   *
! along with FLEXPART.  If not, see <http://www.gnu.org/licenses/>.   *
!**********************************************************************

subroutine releaseparticles(itime)
  !                              o
  !*****************************************************************************
  !                                                                            *
  !     This subroutine releases particles from the release locations.         *
  !                                                                            *
  !     It searches for a "vacant" storage space and assigns all particle      *
  !     information to that space. A space is vacant either when no particle   *
  !     is yet assigned to it, or when it's particle is expired and, thus,     *
  !     the storage space is made available to a new particle.                 *
  !                                                                            *
  !     Author: A. Stohl                                                       *
  !                                                                            *
  !     29 June 2002                                                           *
  !                                                                            *
  !*****************************************************************************
  !                                                                            *
  ! Variables:                                                                 *
  ! itime [s]            current time                                          *
  ! ireleasestart, ireleaseend          start and end times of all releases    *
  ! npart(maxpoint)      number of particles to be released in total           *
  ! numrel               number of particles to be released during this time   *
  !                      step                                                  *
  !                                                                            *
  !*****************************************************************************

  use point_mod
  use xmass_mod
  use par_mod
  use com_mod
  use random_mod, only: ran1

  implicit none

  !real xaux,yaux,zaux,ran1,rfraction,xmasssave(maxpoint)
  real :: xaux,yaux,zaux,rfraction
  real :: topo,rhoaux(2),r,t,rhoout,ddx,ddy,rddx,rddy,p1,p2,p3,p4
  real :: dz1,dz2,dz,xtn,ytn,xlonav,timecorrect(maxspec),press,pressold
  real :: presspart,average_timecorrect
  integer :: itime,numrel,i,j,k,n,ix,jy,ixp,jyp,ipart,minpart,ii
  integer :: indz,indzp,kz,ngrid
  integer :: nweeks,ndayofweek,nhour,jjjjmmdd,ihmmss,mm
  real(kind=dp) :: juldate,julmonday,jul,jullocal,juldiff
  real,parameter :: eps=nxmax/3.e5,eps2=1.e-6

  integer :: idummy = -7
  !save idummy,xmasssave
  !data idummy/-7/,xmasssave/maxpoint*0./



  ! Determine the actual date and time in Greenwich (i.e., UTC + correction for daylight savings time)
  !*****************************************************************************

  julmonday=juldate(19000101,0)          ! this is a Monday
  jul=bdate+real(itime,kind=dp)/86400._dp    ! this is the current day
  call caldate(jul,jjjjmmdd,ihmmss)
  mm=(jjjjmmdd-10000*(jjjjmmdd/10000))/100
  if ((mm.ge.4).and.(mm.le.9)) jul=jul+1._dp/24._dp   ! daylight savings time in summer


  ! For every release point, check whether we are in the release time interval
  !***************************************************************************

  minpart=1
  do i=1,numpoint
    if ((itime.ge.ireleasestart(i)).and. &! are we within release interval?
         (itime.le.ireleaseend(i))) then

  ! Determine the local day and time
  !*********************************

      xlonav=xlon0+(xpoint2(i)+xpoint1(i))/2.*dx  ! longitude needed to determine local time
      if (xlonav.lt.-180.) xlonav=xlonav+360.
      if (xlonav.gt.180.) xlonav=xlonav-360.
      jullocal=jul+real(xlonav,kind=dp)/360._dp   ! correct approximately for time zone to obtain local time

      juldiff=jullocal-julmonday
      nweeks=int(juldiff/7._dp)
      juldiff=juldiff-real(nweeks,kind=dp)*7._dp
      ndayofweek=int(juldiff)+1              ! this is the current day of week, starting with Monday
      nhour=nint((juldiff-real(ndayofweek-1,kind=dp))*24._dp)    ! this is the current hour
      if (nhour.eq.0) then
        nhour=24
        ndayofweek=ndayofweek-1
        if (ndayofweek.eq.0) ndayofweek=7
      endif

  ! Calculate a species- and time-dependent correction factor, distinguishing between
  ! area (those with release starting at surface) and point (release starting above surface) sources
  ! Also, calculate an average time correction factor (species independent)
  !*****************************************************************************
      average_timecorrect=0.
      do k=1,nspec
        if (zpoint1(i).gt.0.5) then      ! point source
          timecorrect(k)=point_hour(k,nhour)*point_dow(k,ndayofweek)
        else                             ! area source
          timecorrect(k)=area_hour(k,nhour)*area_dow(k,ndayofweek)
        endif
        average_timecorrect=average_timecorrect+timecorrect(k)
      end do
      average_timecorrect=average_timecorrect/real(nspec)

  ! Determine number of particles to be released this time; at start and at end of release,
  ! only half the particles are released
  !*****************************************************************************

      if (ireleasestart(i).ne.ireleaseend(i)) then
        rfraction=abs(real(npart(i))*real(lsynctime)/ &
             real(ireleaseend(i)-ireleasestart(i)))
        if ((itime.eq.ireleasestart(i)).or. &
             (itime.eq.ireleaseend(i))) rfraction=rfraction/2.

  ! Take the species-average time correction factor in order to scale the
  ! number of particles released this time
  !**********************************************************************
        rfraction=rfraction*average_timecorrect

        rfraction=rfraction+xmasssave(i)  ! number to be released at this time
        numrel=int(rfraction)
        xmasssave(i)=rfraction-real(numrel)
      else
        numrel=npart(i)
      endif

      xaux=xpoint2(i)-xpoint1(i)
      yaux=ypoint2(i)-ypoint1(i)
      zaux=zpoint2(i)-zpoint1(i)
      do j=1,numrel                       ! loop over particles to be released this time
        do ipart=minpart,maxpart          ! search for free storage space

  ! If a free storage space is found, attribute everything to this array element
  !*****************************************************************************

          if (itra1(ipart).ne.itime) then

  ! Particle coordinates are determined by using a random position within the release volume
  !*****************************************************************************

  ! Determine horizontal particle position
  !***************************************

            xtra1(ipart)=xpoint1(i)+ran1(idummy)*xaux
            if (xglobal) then
              if (xtra1(ipart).gt.real(nxmin1)) xtra1(ipart)= &
                   xtra1(ipart)-real(nxmin1)
              if (xtra1(ipart).lt.0.) xtra1(ipart)= &
                   xtra1(ipart)+real(nxmin1)
            endif
            ytra1(ipart)=ypoint1(i)+ran1(idummy)*yaux

  ! Assign mass to particle: Total mass divided by total number of particles.
  ! Time variation has partly been taken into account already by a species-average
  ! correction factor, by which the number of particles released this time has been
  ! scaled. Adjust the mass per particle by the species-dependent time correction factor
  ! divided by the species-average one
  !*****************************************************************************
            do k=1,nspec
               xmass1(ipart,k)=xmass(i,k)/real(npart(i)) &
                    *timecorrect(k)/average_timecorrect
  !            write (*,*) 'xmass1: ',xmass1(ipart,k),ipart,k
  ! Assign certain properties to particle
  !**************************************
            end do
            nclass(ipart)=min(int(ran1(idummy)*real(nclassunc))+1, &
                 nclassunc)
            numparticlecount=numparticlecount+1
            if (mquasilag.eq.0) then
              npoint(ipart)=i
            else
              npoint(ipart)=numparticlecount
            endif
            idt(ipart)=mintime               ! first time step
            itra1(ipart)=itime
            itramem(ipart)=itra1(ipart)
            itrasplit(ipart)=itra1(ipart)+ldirect*itsplit


  ! Determine vertical particle position
  !*************************************

            ztra1(ipart)=zpoint1(i)+ran1(idummy)*zaux

  ! Interpolation of topography and density
  !****************************************

  ! Determine the nest we are in
  !*****************************

            ngrid=0
            do k=numbnests,1,-1
              if ((xtra1(ipart).gt.xln(k)+eps).and. &
                   (xtra1(ipart).lt.xrn(k)-eps).and. &
                   (ytra1(ipart).gt.yln(k)+eps).and. &
                   (ytra1(ipart).lt.yrn(k)-eps)) then
                ngrid=k
                goto 43
              endif
            end do
43          continue

  ! Determine (nested) grid coordinates and auxiliary parameters used for interpolation
  !*****************************************************************************

            if (ngrid.gt.0) then
              xtn=(xtra1(ipart)-xln(ngrid))*xresoln(ngrid)
              ytn=(ytra1(ipart)-yln(ngrid))*yresoln(ngrid)
              ix=int(xtn)
              jy=int(ytn)
              ddy=ytn-real(jy)
              ddx=xtn-real(ix)
            else
              ix=int(xtra1(ipart))
              jy=int(ytra1(ipart))
              ddy=ytra1(ipart)-real(jy)
              ddx=xtra1(ipart)-real(ix)
            endif
            ixp=ix+1
            jyp=jy+1
            rddx=1.-ddx
            rddy=1.-ddy
            p1=rddx*rddy
            p2=ddx*rddy
            p3=rddx*ddy
            p4=ddx*ddy

            if (ngrid.gt.0) then
              topo=p1*oron(ix ,jy ,ngrid) &
                   + p2*oron(ixp,jy ,ngrid) &
                   + p3*oron(ix ,jyp,ngrid) &
                   + p4*oron(ixp,jyp,ngrid)
            else
              topo=p1*oro(ix ,jy) &
                   + p2*oro(ixp,jy) &
                   + p3*oro(ix ,jyp) &
                   + p4*oro(ixp,jyp)
            endif

  ! If starting height is in pressure coordinates, retrieve pressure profile and convert zpart1 to meters
  !*****************************************************************************
            if (kindz(i).eq.3) then
              presspart=ztra1(ipart)
              do kz=1,nz
                if (ngrid.gt.0) then
                  r=p1*rhon(ix ,jy ,kz,2,ngrid) &
                       +p2*rhon(ixp,jy ,kz,2,ngrid) &
                       +p3*rhon(ix ,jyp,kz,2,ngrid) &
                       +p4*rhon(ixp,jyp,kz,2,ngrid)
                  t=p1*ttn(ix ,jy ,kz,2,ngrid) &
                       +p2*ttn(ixp,jy ,kz,2,ngrid) &
                       +p3*ttn(ix ,jyp,kz,2,ngrid) &
                       +p4*ttn(ixp,jyp,kz,2,ngrid)
                else
                  r=p1*rho(ix ,jy ,kz,2) &
                       +p2*rho(ixp,jy ,kz,2) &
                       +p3*rho(ix ,jyp,kz,2) &
                       +p4*rho(ixp,jyp,kz,2)
                  t=p1*tt(ix ,jy ,kz,2) &
                       +p2*tt(ixp,jy ,kz,2) &
                       +p3*tt(ix ,jyp,kz,2) &
                       +p4*tt(ixp,jyp,kz,2)
                endif
                press=r*r_air*t/100.
                if (kz.eq.1) pressold=press

                if (press.lt.presspart) then
                  if (kz.eq.1) then
                    ztra1(ipart)=height(1)/2.
                  else
                    dz1=pressold-presspart
                    dz2=presspart-press
                    ztra1(ipart)=(height(kz-1)*dz2+height(kz)*dz1) &
                         /(dz1+dz2)
                  endif
                  goto 71
                endif
                pressold=press
              end do
71            continue
            endif

  ! If release positions are given in meters above sea level, subtract the
  ! topography from the starting height
  !***********************************************************************

            if (kindz(i).eq.2) ztra1(ipart)=ztra1(ipart)-topo
            if (ztra1(ipart).lt.eps2) ztra1(ipart)=eps2   ! Minimum starting height is eps2
            if (ztra1(ipart).gt.height(nz)-0.5) ztra1(ipart)= &
                 height(nz)-0.5 ! Maximum starting height is uppermost level - 0.5 meters



  ! For special simulations, multiply particle concentration air density;
  ! Simply take the 2nd field in memory to do this (accurate enough)
  !***********************************************************************
  !AF IND_SOURCE switches between different units for concentrations at the source
  !Af    NOTE that in backward simulations the release of particles takes place at the
  !Af         receptor and the sampling at the source.
  !Af          1="mass"
  !Af          2="mass mixing ratio"
  !Af IND_RECEPTOR switches between different units for concentrations at the receptor
  !Af          1="mass"
  !Af          2="mass mixing ratio"

  !Af switches for the releasefile:
  !Af IND_REL =  1 : xmass * rho
  !Af IND_REL =  0 : xmass * 1

  !Af ind_rel is defined in readcommand.f

            if (ind_rel .eq. 1) then

  ! Interpolate the air density
  !****************************

              do ii=2,nz
                if (height(ii).gt.ztra1(ipart)) then
                  indz=ii-1
                  indzp=ii
                  goto 6
                endif
              end do
6             continue

              dz1=ztra1(ipart)-height(indz)
              dz2=height(indzp)-ztra1(ipart)
              dz=1./(dz1+dz2)

              if (ngrid.gt.0) then
                do n=1,2
                  rhoaux(n)=p1*rhon(ix ,jy ,indz+n-1,2,ngrid) &
                       +p2*rhon(ixp,jy ,indz+n-1,2,ngrid) &
                       +p3*rhon(ix ,jyp,indz+n-1,2,ngrid) &
                       +p4*rhon(ixp,jyp,indz+n-1,2,ngrid)
                end do
              else
                do n=1,2
                  rhoaux(n)=p1*rho(ix ,jy ,indz+n-1,2) &
                       +p2*rho(ixp,jy ,indz+n-1,2) &
                       +p3*rho(ix ,jyp,indz+n-1,2) &
                       +p4*rho(ixp,jyp,indz+n-1,2)
                end do
              endif
              rhoout=(dz2*rhoaux(1)+dz1*rhoaux(2))*dz
              rho_rel(i)=rhoout


  ! Multiply "mass" (i.e., mass mixing ratio in forward runs) with density
  !********************************************************************

              do k=1,nspec
                xmass1(ipart,k)=xmass1(ipart,k)*rhoout
              end do
            endif


            numpart=max(numpart,ipart)
            goto 34      ! Storage space has been found, stop searching
          endif
        end do
        if (ipart.gt.maxpart) goto 996

34      minpart=ipart+1
      end do
      endif
  end do


  return

996   continue
  write(*,*) '#####################################################'
  write(*,*) '#### FLEXPART MODEL SUBROUTINE RELEASEPARTICLES: ####'
  write(*,*) '####                                             ####'
  write(*,*) '#### ERROR - TOTAL NUMBER OF PARTICLES REQUIRED  ####'
  write(*,*) '#### EXCEEDS THE MAXIMUM ALLOWED NUMBER. REDUCE  ####'
  write(*,*) '#### EITHER NUMBER OF PARTICLES PER RELEASE POINT####'
  write(*,*) '#### OR REDUCE NUMBER OF RELEASE POINTS.         ####'
  write(*,*) '#####################################################'
  stop

end subroutine releaseparticles
