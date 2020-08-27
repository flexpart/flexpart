! SPDX-FileCopyrightText: FLEXPART 1998-2019, see flexpart_license.txt
! SPDX-License-Identifier: GPL-3.0-or-later

subroutine coordtrafo

  !**********************************************************************
  !                                                                     *
  !             FLEXPART MODEL SUBROUTINE COORDTRAFO                    *
  !                                                                     *
  !**********************************************************************
  !                                                                     *
  !             AUTHOR:      G. WOTAWA                                  *
  !             DATE:        1994-02-07                                 *
  !             LAST UPDATE: 1996-05-18   A. STOHL                      *
  !                                                                     *
  !**********************************************************************
  !                                                                     *
  ! DESCRIPTION: This subroutine transforms x and y coordinates of      *
  ! particle release points to grid coordinates.                        *
  !                                                                     *
  !**********************************************************************

  use point_mod
  use par_mod
  use com_mod

  implicit none

  integer :: i,j,k
  real :: yrspc ! small real number relative to x

  if (numpoint.eq.0) goto 30

  ! TRANSFORM X- AND Y- COORDINATES OF STARTING POINTS TO GRID COORDINATES
  !***********************************************************************

  do i=1,numpoint
    xpoint1(i)=(xpoint1(i)-xlon0)/dx
    xpoint2(i)=(xpoint2(i)-xlon0)/dx
    ypoint1(i)=(ypoint1(i)-ylat0)/dy
    ypoint2(i)=(ypoint2(i)-ylat0)/dy
  end do

15   continue


  ! CHECK IF RELEASE POINTS ARE WITHIN DOMAIN
  !******************************************
  
  yrspc = spacing(real(nymin1,kind=sp))
  
  do i=1,numpoint
    if (sglobal.and.(ypoint1(i).lt.1.e-6)) ypoint1(i)=1.e-6
    if (nglobal.and.(ypoint2(i).gt.real(nymin1,kind=dp)-1.e-5)) &
         ypoint2(i)=real(nymin1,kind=dp)-10*yrspc
    if ((ypoint1(i).lt.1.e-6).or.(ypoint1(i).ge.real(nymin1,kind=dp)-1.e-6) &
       .or.(ypoint2(i).lt.1.e-6).or.(ypoint2(i).ge.real(nymin1,kind=dp)-yrspc) &
       .or.((.not.xglobal).and.((xpoint1(i).lt.1.e-6).or. &
       (xpoint1(i).ge.real(nxmin1,kind=dp)-1.e-6).or.(xpoint2(i).lt.1.e-6).or. &
       (xpoint2(i).ge.real(nxmin1,kind=dp)-1.e-6)))) then
      write(*,*) ' NOTICE: RELEASE POINT OUT OF DOMAIN DETECTED.'
      write(*,*) ' IT IS REMOVED NOW ... '
      if (i.le.1000) then
         write(*,*) ' COMMENT: ',compoint(i)
      else
         write(*,*) ' COMMENT: ',compoint(1001)
      endif
      if (i.lt.numpoint) then
        do j=i+1,numpoint
          xpoint1(j-1)=xpoint1(j)
          ypoint1(j-1)=ypoint1(j)
          xpoint2(j-1)=xpoint2(j)
          ypoint2(j-1)=ypoint2(j)
          zpoint1(j-1)=zpoint1(j)
          zpoint2(j-1)=zpoint2(j)
          npart(j-1)=npart(j)
          kindz(j-1)=kindz(j)
          ireleasestart(j-1)=ireleasestart(j)
          ireleaseend(j-1)=ireleaseend(j)
          if (j.le.1000) compoint(j-1)=compoint(j)
          do k=1,nspec
            xmass(j-1,k)=xmass(j,k)
          end do
        end do
      endif

      numpoint=numpoint-1
      if (numpoint.gt.0) goto 15
    endif
  end do

30   if(numpoint.eq.0) then
    write(*,*) ' FLEXPART MODEL SUBROUTINE COORDTRAFO: ERROR ! '
    write(*,*) ' NO PARTICLE RELEASES ARE DEFINED!'
    write(*,*) ' CHECK FILE RELEASES...'
    stop
  endif

end subroutine coordtrafo
