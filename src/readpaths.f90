! SPDX-FileCopyrightText: FLEXPART 1998-2019, see flexpart_license.txt
! SPDX-License-Identifier: GPL-3.0-or-later

subroutine readpaths

  !*****************************************************************************
  !                                                                            *
  !     Reads the pathnames, where input/output files are expected to be.      *
  !     The file pathnames must be available in the current working directory. *
  !                                                                            *
  !     Author: A. Stohl                                                       *
  !                                                                            *
  !     1 February 1994                                                        *
  !     last modified                                                          *
  !     HS, 7.9.2012                                                           *
  !     option to give pathnames file as command line option                   *
  !                                                                            *
  !*****************************************************************************
  !                                                                            *
  ! Variables:                                                                 *
  ! length(numpath)    lengths of the path names                               *
  ! path(numpath)      pathnames of input/output files                         *
  !                                                                            *
  ! Constants:                                                                 *
  ! numpath            number of pathnames to be read in                       *
  !                                                                            *
  !*****************************************************************************

  use par_mod
  use com_mod

  implicit none

  integer   :: i
  character(256) :: string_test 
  character(1) :: character_test 

  ! Read the pathname information stored in unitpath
  !*************************************************

  open(unitpath,file=trim(pathfile),status='old',err=999)

  do i=1,numpath
    read(unitpath,'(a)',err=998) path(i)
    length(i)=index(path(i),' ')-1

    
    string_test = path(i)
    character_test = string_test(length(i):length(i))
    !print*, 'character_test,  string_test ', character_test,  string_test 
      if ((character_test .NE. '/') .AND. (i .LT. 4))  then
         print*, 'WARNING: path not ending in /' 
         print*, path(i)
         path(i) = string_test(1:length(i)) // '/'
         length(i)=length(i)+1
         print*, 'fix: padded with /' 
         print*, path(i)
         print*, 'length(i) increased 1' 
      endif
  end do

  ! Check whether any nested subdomains are to be used
  !***************************************************

  do i=1,maxnests
  ! ESO 2016 Added 'end'/'err' in case user forgot '====' at end of file and
  ! maxnests > numbnests
    read(unitpath,'(a)', end=30, err=30) path(numpath+2*(i-1)+1)
    read(unitpath,'(a)', end=30, err=30) path(numpath+2*(i-1)+2)
    if (path(numpath+2*(i-1)+1)(1:5).eq.'=====') goto 30
    length(numpath+2*(i-1)+1)=index(path(numpath+2*(i-1)+1),' ')-1
    length(numpath+2*(i-1)+2)=index(path(numpath+2*(i-1)+2),' ')-1
  end do


  ! Determine number of available nested domains
  !*********************************************

30   numbnests=i-1

  close(unitpath)
  return

998   write(*,*) ' #### TRAJECTORY MODEL ERROR! ERROR WHILE     #### '
  write(*,*) ' #### READING FILE PATHNAMES.                 #### '
  stop

999   write(*,*) ' #### TRAJECTORY MODEL ERROR! FILE "pathnames"#### '
  write(*,*) ' #### CANNOT BE OPENED IN THE CURRENT WORKING #### '
  write(*,*) ' #### DIRECTORY.                              #### '
  stop

end subroutine readpaths
