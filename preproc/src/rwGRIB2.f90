MODULE RWGRIB2

CONTAINS

  SUBROUTINE READLATLON(FILENAME,FELD,MAXL,MAXB,MLEVEL,MPAR)

    USE GRIB_API

    IMPLICIT NONE

    integer                            ::  ifile
    integer                            ::  iret
    integer                            ::  n,mk,parid,nm
    integer                            ::  i,k
    integer,dimension(:),allocatable   ::  igrib
    integer                            ::  numberOfPointsAlongAParallel
    integer                            ::  numberOfPointsAlongAMeridian
    real, dimension(:), allocatable    ::  values
    integer                            ::  numberOfValues
    real,dimension(maxl,maxb,mlevel)   ::  feld  
    integer:: maxl,maxb,mlevel,mstride,mpar(:),irest,div
    integer :: l(size(mpar))
    character*(*):: filename                             

    call grib_open_file(ifile, TRIM(FILENAME),'r')

    ! count the messages in the file
    call grib_count_in_file(ifile,n)
    allocate(igrib(n))
    igrib=-1

    ! Load the messages from the file.
    DO i=1,n
       call grib_new_from_file(ifile,igrib(i), iret)
    END DO

    ! we can close the file
    call grib_close_file(ifile)

    nm=size(mpar)
    div=mlevel/nm
    l=0

    ! Loop on all the messages in memory
    iloop:  DO i=1,n
       !      write(*,*) 'processing message number ',i
       !     get as a integer
       call grib_get(igrib(i),'numberOfPointsAlongAParallel', &
            numberOfPointsAlongAParallel)

       !     get as a integer
       call grib_get(igrib(i),'numberOfPointsAlongAMeridian', &
            numberOfPointsAlongAMeridian)

       call grib_get(igrib(i),'numberOfVerticalCoordinateValues',mk)

       call grib_get_size(igrib(i),'values',numberOfValues)
       !      write(*,*) 'numberOfValues=',numberOfValues

       allocate(values(numberOfValues), stat=iret)
       !     get data values
       call grib_get(igrib(i),'values',values)

       call grib_get(igrib(i),'paramId',parid)

       kloop:  do k=1,nm
          if(parid .eq. mpar(k)) then
             l(k)=l(k)+1
             feld(:,:,(k-1)*div+l(k))=reshape(values,(/maxl,maxb/))
             !         print*,(k-1)*div+l(k),parid
             exit kloop
          endif
       enddo kloop
       if(k .gt. nm .and. parid .ne. mpar(nm)) then
          write(*,*) k,nm,parid,mpar(nm)
          write(*,*) 'ERROR readlatlon: parameter ',parid,'is not',mpar
          stop
       endif

       !      print*,i
      deallocate(values)
    END DO iloop
    write(*,*) 'readlatlon: ',i-1,' records read'

    DO i=1,n
       call grib_release(igrib(i))
    END DO

    deallocate(igrib)

  END SUBROUTINE READLATLON

  SUBROUTINE WRITELATLON(iunit,igrib,FELD,MAXL,MAXB,MLEVEL,&
       MLEVELIST,MPAR)

    USE GRIB_API

    IMPLICIT NONE

    INTEGER IFIELD,MLEVEL,MNAUF,I,J,K,L,IERR,JOUT
    INTEGER MPAR,MAXL,MAXB,LEVMIN,LEVMAX
    INTEGER IUNIT,igrib
    REAL ZSEC4(MAXL*MAXB)
    REAL    FELD(MAXL,MAXB,MLEVEL)
    CHARACTER*(*) MLEVELIST
    INTEGER ILEVEL(MLEVEL),MLINDEX(MLEVEL+1),LLEN

    ! parse MLEVELIST

    LLEN=len(trim(MLEVELIST))
    if(index(MLEVELIST,'to') .ne. 0 .or. index(MLEVELIST,'TO') .ne. 0) THEN
       i=index(MLEVELIST,'/')
       read(MLEVELIST(1:i-1),*) LEVMIN
       i=index(MLEVELIST,'/',.true.)
       read(MLEVELIST(i+1:LLEN),*) LEVMAX
       l=0
       do i=LEVMIN,LEVMAX
          l=l+1
          ILEVEL(l)=i
       enddo
    else
       l=1
       MLINDEX(1)=0
       do i=1,LLEN
          if(MLEVELIST(i:i) .eq. '/') THEN
             l=l+1
             MLINDEX(l)=i
	  endif
       enddo
       MLINDEX(l+1)=LLEN+1
       do i=1,l
	  read(MLEVELIST(MLINDEX(i)+1:MLINDEX(i+1)-1),*) ILEVEL(i)
       enddo
    endif

    DO k=1,l
       call grib_set(igrib,"level",ILEVEL(k))
       call grib_set(igrib,"paramId",MPAR)
       if(l .ne. mlevel) then
          zsec4(1:maxl*maxb)=RESHAPE(FELD(:,:,ILEVEL(k)),(/maxl*maxb/))
       else
          zsec4(1:maxl*maxb)=RESHAPE(FELD(:,:,k),(/maxl*maxb/))
       endif
       call grib_set(igrib,"values",zsec4)

       call grib_write(igrib,iunit)

    ENDDO

  END SUBROUTINE WRITELATLON

  SUBROUTINE READSPECTRAL(FILENAME,CXMN,MNAUF,MLEVEL,&
       MAXLEV,MPAR)

    USE GRIB_API

    IMPLICIT NONE

    character*(*),               intent(in)    :: filename      

    integer                            ::  ifile
    integer                            ::  iret
    integer                            ::  n,mk,div,nm,k
    integer                            ::  i,j,parid
    integer,dimension(:),allocatable   ::  igrib
    real, dimension(:), allocatable    ::  values
    integer                            ::  numberOfValues,maxlev  
    REAL:: CXMN(0:(MNAUF+1)*(MNAUF+2)-1,MLEVEL)
    integer:: maxl,maxb,mlevel,mstride,mpar(:),mnauf,ioffset,ipar,ilev,l(size(mpar))

    call grib_open_file(ifile, TRIM(FILENAME),'r')

    ! count the messages in the file
    call grib_count_in_file(ifile,n)
    allocate(igrib(n))
    igrib=-1

    ! Load the messages from the file.
    DO i=1,n
       call grib_new_from_file(ifile,igrib(i), iret)
    END DO

    ! we can close the file
    call grib_close_file(ifile)

    l=0
    ! Loop on all the messages in memory
    iloop: DO i=1,n
       ! write(*,*) 'processing message number ',i
       !     get as a integer
       call grib_get(igrib(i),'pentagonalResolutionParameterJ', j)

       call grib_get_size(igrib(i),'values',numberOfValues)
       !   write(*,*) 'numberOfValues=',numberOfValues

       call grib_get(igrib(i),'numberOfVerticalCoordinateValues',mk)

       call grib_get(igrib(i),'level',ilev)

       allocate(values(numberOfValues), stat=iret)
       !     get data values
       call grib_get(igrib(i),'values',values)

       call grib_get(igrib(i),'paramId',parid)
       nm=size(mpar)
       div=mlevel/nm
       kloop:  do k=1,nm
          if(parid .eq. mpar(k)) then
             l(k)=l(k)+1
             cxmn(:,(k-1)*div+l(k))=values(1:(MNAUF+1)*(MNAUF+2))
             !         print*,(k-1)*div+l(k),parid
             exit kloop
          endif

       enddo kloop
       if(k .gt. nm .and. parid .ne. mpar(nm)) then
          write(*,*) k,nm,parid,mpar(nm)
          write(*,*) 'ERROR readspectral: parameter ',parid,'is not',mpar
          stop
       endif

       !      print*,i
      deallocate(values)

    END DO iloop

    write(*,*) 'readspectral: ',i-1,' records read'

    DO i=1,n
       call grib_release(igrib(i))
    END DO

    deallocate(igrib)

  END SUBROUTINE READSPECTRAL

END MODULE RWGRIB2
