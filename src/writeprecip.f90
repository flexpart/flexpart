subroutine writeprecip(itime,imem)

  !*****************************************************************************
  !                                                                            *
  !  This routine produces a file containing total precipitation for each      *
  !  releases point.                                                           *
  !                                                                            *
  !     Author: S. Eckhardt                                                    * 
  !     7 Mai 2017                                                             *
  !*****************************************************************************

  use point_mod
  use par_mod
  use com_mod

  implicit none

  integer :: jjjjmmdd,ihmmss,itime,i
  real(kind=dp) :: jul
  character :: adate*8,atime*6

  integer :: ix,jy,imem
  real :: xp1,yp1

  
  if (itime.eq.0) then
      open(unitprecip,file=path(2)(1:length(2))//'wetscav_precip.txt', &
       form='formatted',err=998)
  else
      open(unitprecip,file=path(2)(1:length(2))//'wetscav_precip.txt', &
       ACCESS='APPEND',form='formatted',err=998)
  endif

  jul=bdate+real(itime,kind=dp)/86400._dp
  call caldate(jul,jjjjmmdd,ihmmss)
  write(adate,'(i8.8)') jjjjmmdd
  write(atime,'(i6.6)') ihmmss

  do i=1,numpoint
    xp1=xpoint1(i)*dx+xlon0 !lat, long (real) coord
    yp1=ypoint1(i)*dy+ylat0 !lat, long (real) coord
    ix=int((xpoint1(i)+xpoint2(i))/2.)
    jy=int((ypoint1(i)+ypoint2(i))/2.)
    write(unitprecip,*)  jjjjmmdd, ihmmss, & 
           xp1,yp1,lsprec(ix,jy,1,imem),convprec(ix,jy,1,imem) !time is the same as in the ECMWF windfield
! units mm/h, valid for the time given in the windfield
  end do

  close(unitprecip)

  return


998   write(*,*) ' #### FLEXPART MODEL ERROR!   THE FILE         #### '
  write(*,*) ' #### '//path(2)(1:length(2))//'header_txt'//' #### '
  write(*,*) ' #### CANNOT BE OPENED. IF A FILE WITH THIS    #### '
  write(*,*) ' #### NAME ALREADY EXISTS, DELETE IT AND START #### '
  write(*,*) ' #### THE PROGRAM AGAIN.                       #### '
  stop

end subroutine writeprecip
