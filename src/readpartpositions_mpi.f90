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
  use mpi_mod, only: mp_seed

  implicit none

  integer :: ibdatein,ibtimein,nspecin,itimein,numpointin,i,j,ix
  integer :: id1,id2,it1,it2
  real :: xlonin,ylatin,topo,hmixi,pvi,qvi,rhoi,tri,tti
  character :: specin*7
  real(kind=dp) :: julin,julpartin,juldate

  integer :: idummy = -8

  ! Different seed for each process
  idummy=idummy+mp_seed

  numparticlecount=0

  ! Open header file of dumped particle data
  !*****************************************

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
  if (numpointin.ne.numpoint) goto 995
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


  ! Open data file of dumped particle data
  !***************************************

  close(unitpartin)
  open(unitpartin,file=path(2)(1:length(2))//'partposit_end', &
       form='unformatted',err=998)


100   read(unitpartin,end=99) itimein
    i=0
200   i=i+1
    read(unitpartin) npoint(i),xlonin,ylatin,ztra1(i),itramem(i), &
         topo,pvi,qvi,rhoi,hmixi,tri,tti,(xmass1(i,j),j=1,nspec)

    if (xlonin.eq.-9999.9) goto 100
    xtra1(i)=(xlonin-xlon0)/dx
    ytra1(i)=(ylatin-ylat0)/dy
    numparticlecount=max(numparticlecount,npoint(i))
    goto 200

99   numpart=i-1

  close(unitpartin)

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

!995   write(*,*) ' #### FLEXPART MODEL ERROR IN READPARTPOSITIONS#### '
995   write(*,*) ' #### FLEXPART MODEL WARNING IN READPARTPOSITIONS#### '
  write(*,*) ' #### NUMBER OF RELEASE LOCATIONS DOES NOT     #### '
  write(*,*) ' #### AGREE WITH CURRENT SETTINGS!             #### '
!  stop
  goto 999 

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
