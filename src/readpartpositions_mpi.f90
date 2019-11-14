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

subroutine readpartpositions

  !*****************************************************************************
  !                                                                            *
  !   This routine opens the particle dump file and reads all the particle     *
  !   positions from a previous run to initialize the current run.             *
  !                                                                            *
  !                                                                            *
  !     Author: A. Stohl                                                       *
  !                                                                            *
  !     24 March 2000                                                          *
  !                                                                            *
  !   CHANGES                                                                  *
  !     12/2014 eso: MPI version                                               *
  !                  Root process reads positions and distributes the data     *
  !                  :TODO: The above..                                        *
  !*****************************************************************************
  !                                                                            *
  ! Variables:                                                                 *
  !                                                                            *
  !*****************************************************************************

  use par_mod
  use com_mod
  use random_mod, only: ran1
  use mpi_mod !, only: mp_seed

  implicit none

  integer :: ibdatein,ibtimein,nspecin,itimein,numpointin,i,j,ix,ip
  integer :: id1,id2,it1,it2
  integer :: addone,numparticlecount_all,numpart_all,lbnd,ubnd
  real :: xlonin,ylatin,topo,hmixi,pvi,qvi,rhoi,tri,tti
  character :: specin*7
  real(kind=dp) :: julin,julpartin,juldate

  integer :: idummy = -8

  ! These variables are allocated at the root process for all particles in file.
  real,dimension(maxpart) :: xtra1_all,ytra1_all,ztra1_all
  real,dimension(maxpart,maxspec) :: xmass1_all
  integer,dimension(maxpart) :: npoint_all,itramem_all

  ! Different seed for each process
  idummy=idummy+mp_seed

  numparticlecount=0
  numparticlecount_all=0
  numpart_all=0

  ! Open header file of dumped particle data
  ! Each MPI process sequentially access file (just in case)
  !*********************************************************

  do ip=0, mp_partgroup_np-1
    call mpif_mpi_barrier
    open(unitpartin,file=path(2)(1:length(2))//'header', &
         form='unformatted',err=998)

    read(unitpartin) ibdatein,ibtimein
    read(unitpartin)
    read(unitpartin)

    read(unitpartin)
    read(unitpartin)
    read(unitpartin) nspecin
    nspecin=nspecin/3
    if ((ldirect.eq.1).and.(nspec.ne.nspecin)) goto 997

    do i=1,nspecin
      read(unitpartin)
      read(unitpartin)
      read(unitpartin) j,specin
      if ((ldirect.eq.1).and.(species(i)(1:7).ne.specin)) goto 996
    end do
    
    read(unitpartin) numpointin
    if (numpointin.ne.numpoint) then ! goto 995

! eso 2016: moved this warning here to avoid out-of-block goto
!995   write(*,*) ' #### FLEXPART MODEL WARNING IN READPARTPOSITIONS#### '
      write(*,*) ' #### FLEXPART MODEL WARNING IN READPARTPOSITIONS#### '
      write(*,*) ' #### NUMBER OF RELEASE LOCATIONS DOES NOT     #### '
      write(*,*) ' #### AGREE WITH CURRENT SETTINGS!             #### '
!  stop
      goto 999 
    end if

999 continue 
    do i=1,numpointin
      read(unitpartin)
      read(unitpartin)
      read(unitpartin)
      read(unitpartin)
      do j=1,nspec
        read(unitpartin)
        read(unitpartin)
        read(unitpartin)
      end do
    end do
    read(unitpartin)
    read(unitpartin)

    do ix=0,numxgrid-1
      read(unitpartin)
    end do
    close(unitpartin)


  ! Open data file of dumped particle data
  ! All processes read the whole file
  !***************************************

    open(unitpartin,file=path(2)(1:length(2))//'partposit_end', &
         form='unformatted',err=998)
    
100 read(unitpartin,end=99) itimein
    i=0
200 i=i+1
    read(unitpartin) npoint_all(i),xlonin,ylatin,ztra1_all(i),itramem_all(i), &
         topo,pvi,qvi,rhoi,hmixi,tri,tti,(xmass1_all(i,j),j=1,nspec)
    
    if (xlonin.eq.-9999.9) goto 100
    xtra1_all(i)=(xlonin-xlon0)/dx
    ytra1_all(i)=(ylatin-ylat0)/dy
    numparticlecount_all=max(numparticlecount_all,npoint(i))
    goto 200

99  numpart_all=i-1
    
    close(unitpartin)

  end do


  ! Distribute particles among processes
  !**************************************
  lbnd=1
  ubnd=0

  do ip=0, mp_partgroup_np-1
    
! Extra particle distributed in case remainder in the division
    if (ip.lt.mod(numpart_all,mp_partgroup_np)) then
      addone=1
    else
      addone=0
    end if

    if (ip==0) then 
      ubnd=ubnd + numpart_all/mp_partgroup_np + addone
    else 
      ubnd=lbnd + numpart_all/mp_partgroup_np + addone - 1
    end if
    
    if (ip==mp_pid) then
      
      numparticlecount=numparticlecount_all/mp_partgroup_np+addone
      numpart=numpart_all/mp_partgroup_np+addone

      xtra1(1:numpart) = xtra1_all(lbnd:ubnd)
      ytra1(1:numpart) = ytra1_all(lbnd:ubnd)
      ztra1(1:numpart) = ztra1_all(lbnd:ubnd)
      xmass1(1:numpart,:) = xmass1_all(lbnd:ubnd,:)
      npoint(1:numpart) = npoint_all(lbnd:ubnd)
      itramem(1:numpart) = itramem_all(lbnd:ubnd)

    end if

    lbnd=ubnd+1

  end do

  julin=juldate(ibdatein,ibtimein)+real(itimein,kind=dp)/86400._dp
  if (abs(julin-bdate).gt.1.e-5) goto 994
  do i=1,numpart
    julpartin=juldate(ibdatein,ibtimein)+ &
         real(itramem(i),kind=dp)/86400._dp
    nclass(i)=min(int(ran1(idummy)*real(nclassunc))+1, &
         nclassunc)
    idt(i)=mintime
    itra1(i)=0
    itramem(i)=nint((julpartin-bdate)*86400.)
    itrasplit(i)=ldirect*itsplit
  end do

  return


994   write(*,*) ' #### FLEXPART MODEL ERROR IN READPARTPOSITIONS#### '
  write(*,*) ' #### ENDING TIME OF PREVIOUS MODEL RUN DOES   #### '
  write(*,*) ' #### NOT AGREE WITH STARTING TIME OF THIS RUN.#### '
  call caldate(julin,id1,it1)
  call caldate(bdate,id2,it2)
  write(*,*) 'julin: ',julin,id1,it1
  write(*,*) 'bdate: ',bdate,id2,it2
  stop

996   write(*,*) ' #### FLEXPART MODEL ERROR IN READPARTPOSITIONS#### '
  write(*,*) ' #### SPECIES NAMES TO BE READ IN DO NOT       #### '
  write(*,*) ' #### AGREE WITH CURRENT SETTINGS!             #### '
  stop

997   write(*,*) ' #### FLEXPART MODEL ERROR IN READPARTPOSITIONS#### '
  write(*,*) ' #### THE NUMBER OF SPECIES TO BE READ IN DOES #### '
  write(*,*) ' #### NOT AGREE WITH CURRENT SETTINGS!         #### '
  stop

998   write(*,*) ' #### FLEXPART MODEL ERROR!   THE FILE         #### '
  write(*,*) ' #### '//path(2)(1:length(2))//'grid'//' #### '
  write(*,*) ' #### CANNOT BE OPENED. IF A FILE WITH THIS    #### '
  write(*,*) ' #### NAME ALREADY EXISTS, DELETE IT AND START #### '
  write(*,*) ' #### THE PROGRAM AGAIN.                       #### '
  stop

end subroutine readpartpositions
