subroutine caldate(juldate,yyyymmdd,hhmiss)
  !                      i       o       o
  !*****************************************************************************
  !                                                                            *
  !     Calculates the Gregorian date from the Julian date                     *
  !                                                                            *
  !     AUTHOR: Andreas Stohl (21 January 1994), adapted from Numerical Recipes*
  !                                                                            *
  !     Variables:                                                             *
  !     dd             Day                                                     *
  !     hh             Hour                                                    *
  !     hhmiss         Hour, Minute, Second                                    *
  !     ja,jb,jc,jd,je help variables                                          *
  !     jalpha         help variable                                           *
  !     juldate        Julian Date                                             *
  !     julday         help variable                                           *
  !     mi             Minute                                                  *
  !     mm             Month                                                   *
  !     ss             Seconds                                                 *
  !     yyyy           Year                                                    *
  !     yyyymmdd       Year, Month, Day                                        *
  !                                                                            *
  !     Constants:                                                             *
  !     igreg          help constant                                           *
  !                                                                            *
  !*****************************************************************************

  use par_mod, only: dp

  implicit none

  integer           :: yyyymmdd,yyyy,mm,dd,hhmiss,hh,mi,ss
  integer           :: julday,ja,jb,jc,jd,je,jalpha
  real(kind=dp)     :: juldate
  integer,parameter :: igreg=2299161

  julday=int(juldate)
  if(julday.ge.igreg)then
    jalpha=int(((julday-1867216)-0.25)/36524.25)
    ja=julday+1+jalpha-int(0.25*jalpha)
  else
    ja=julday
  endif
  jb=ja+1524
  jc=int(6680.+((jb-2439870)-122.1)/365.25)
  jd=365*jc+int(0.25*jc)
  je=int((jb-jd)/30.6001)
  dd=jb-jd-int(30.6001*je)
  mm=je-1
  if (mm.gt.12) mm=mm-12
  yyyy=jc-4715
  if (mm.gt.2) yyyy=yyyy-1
  if (yyyy.le.0) yyyy=yyyy-1

  yyyymmdd=10000*yyyy+100*mm+dd
  hh=int(24._dp*(juldate-real(julday,kind=dp)))
  mi=int(1440._dp*(juldate-real(julday,kind=dp))-60._dp*real(hh,kind=dp))
  ss=nint(86400._dp*(juldate-real(julday,kind=dp))-3600._dp*real(hh,kind=dp)- &
       60._dp*real(mi,kind=dp))
  if (ss.eq.60) then  ! 60 seconds = 1 minute
    ss=0
    mi=mi+1
  endif
  if (mi.eq.60) then
    mi=0
    hh=hh+1
  endif
  hhmiss=10000*hh+100*mi+ss

end subroutine caldate
