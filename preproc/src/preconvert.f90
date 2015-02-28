PROGRAM PRECONVERT
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !                                                                 !
  ! PROGRAM PRECONVERT - PREPARES INPUT DATA FOR POP MODEL METEOR-  !
  !                      OLOGICAL PREPROCESSOR                      !
  !                                                                 !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !                                                                 !
  ! CALCULATION OF ETAPOINT ON A REGULAR LAMDA/PHI GRID AND WRITING !
  ! U,V,ETAPOINT,T,PS,Q,SD,MSL,TCC,10U, 10V, 2T,2D,LSP,CP,SSHF,SSR, !
  ! EWSS,NSSS TO AN OUTPUT FILE (GRIB 1 or 2 FORMAT).               ! 
  !                                                                 !
  ! AUTHORS: L. HAIMBERGER, G. WOTAWA, 1994-04                      !
  !                     adapted: A. BECK                            !
  !                     2003-05-11                                  !
  !          L. Haimberger 2006-12    V2.0                          !
  !                    modified to handle arbitrary regular grids   !
  !                    and T799 resolution data                     !
  !          L. Haimberger 2010-03    V4.0                          !
  !                    modified to grib edition 2 fields            !
  !                    and T1279 resolution data                    !
  !                                                                 !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !                                                                 !
  ! DESCRIPTION OF NEEDED INPUT:                                    !
  !                                                                 !
  !       FILE      PARAMETER(S)    DATA REPRESENTATION             !
  !                                                                 !
  !       fort.10   U,V             spherical harmonics             !
  !       fort.11   T               regular lamda phi grid          !
  !       fort.12   LNSP            spherical harmonics             !
  !       fort.13   D               spherical harmonics             !
  !       fort.14   SD,MSL,TCC,10U,                                 !
  !                 10V,2T,2D       regular lamda phi grid          !
  !       fort.18   Q               regular lamda phi grid          !
  !                                                                 !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !                                                                 !
  ! DESCRIPTION OF OUTPUT:                                          !
  !                                                                 !
  ! UNIT  FILE      PARAMETER(S)    DATA REPRESENTATION             !
  !                                                                 !
  ! 15    fort.15   U,V,ETA,T,PS,                                   !
  !                 Q,SD,MSL,TCC,                                   !
  !                 10U,10V,2T,2D,  regular lamda phi grid          !
  !                 LSP,CP,SSHF,                                    !
  !                 SSR,EWSS,NSSS                                   !
  !                                                                 !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !

  USE PHTOGR
  USE GRTOPH
  USE FTRAFO
  USE RWGRIB2
  USE GRIB_API

  IMPLICIT NONE

  REAL, ALLOCATABLE, DIMENSION (:,:)   :: Z
  REAL, ALLOCATABLE, DIMENSION (:,:,:) :: T, UV 
  REAL, ALLOCATABLE, DIMENSION (:,:,:) :: DIV, ETA
  REAL, ALLOCATABLE, DIMENSION (:,:)   :: DPSDL, DPSDM
  REAL, ALLOCATABLE, DIMENSION (:,:,:) :: PS,DPSDT
  REAL, ALLOCATABLE, DIMENSION (:)     :: WSAVE,H
  REAL, ALLOCATABLE, DIMENSION (:)     :: BREITE, GBREITE
  REAL, ALLOCATABLE, DIMENSION (:)     :: AK, BK, PV
  INTEGER                              :: NPV

  ! Arrays for Gaussian grid calculations

  REAL  :: X1,X2,RMS,MW,SIG,LAM
  REAL,ALLOCATABLE :: CUA(:,:,:),CVA(:,:,:)

  REAL, ALLOCATABLE, DIMENSION (:,:) :: P,PP !,P2
  REAL, ALLOCATABLE, DIMENSION (:,:) :: XMN,HILFUV
  REAL, ALLOCATABLE, DIMENSION (:)   :: LNPMN,LNPMN2,LNPMN3
  REAL, ALLOCATABLE, DIMENSION (:)   :: WEIGHT
  REAL, ALLOCATABLE, DIMENSION (:,:) :: UGVG
  REAL, ALLOCATABLE, DIMENSION (:,:) :: DG, ETAG
  REAL, ALLOCATABLE, DIMENSION (:,:) :: GWSAVE
  REAL, ALLOCATABLE, DIMENSION (:)   :: PSG,HILF

  ! end arrays for Gaussian grid calculations

  INTEGER, ALLOCATABLE, DIMENSION (:) :: MLAT
  INTEGER, ALLOCATABLE :: GIFAX(:,:)

  REAL, PARAMETER   :: PI=ACOS(-1.D0)

  REAL COSB,DAK,DBK,P00
  REAL URLAR8,JMIN1,LLLAR8,MAXBMIN1,PIR8,DCOSB

  INTEGER I,J,K,L,IERR,M,MK,NGI,NGJ
  INTEGER LUNIT,LUNIT_OUT

  INTEGER MAXL, MAXB, MLEVEL, LEVOUT,LEVMIN,LEVMAX
  INTEGER MGAUSS,MSMOOTH, MNAUF,META
  INTEGER MDPDETA,METAPAR
  REAL RLO0, RLO1, RLA0, RLA1
  CHARACTER*300 MLEVELIST

  INTEGER MAUF, MANF,IFAX(10)

  INTEGER igrib,iret

  CHARACTER*80 FILENAME

  NAMELIST /NAMGEN/ &
       MAXL, MAXB,  &
       MLEVEL,MLEVELIST,MNAUF,METAPAR, &
       RLO0, RLO1, RLA0, RLA1, &
       MGAUSS,MSMOOTH,META,&
       MDPDETA

  read (4,NAMGEN)

  PRINT*, 'MAXL= ', MAXL, ' RLO0 = ', RLO0, ' RLO1 = ', RLO1
  MAUF=INT(360.*(REAL(MAXL)-1.)/(RLO1-RLO0)+0.0001)
  PRINT*, 'MAUF= ', MAUF

  MANF=INT(REAL(MAUF)/360.*(360.+RLO0)+1.0001)
  PRINT*, 'MANF= ', MANF

  IF(MANF .gt. MAUF) MANF=MANF-MAUF


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !     GAUSS STUFF !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! computation of etadot on gaussian grid

! a new grib message is loaded from an existing sample.
! we want a regular gaussian grid
! set the environment variable GRIB_SAMPLES_PATH
  write (filename, "(A13,I0,A3,I0)") "EI_regular_gg",MNAUF+1,"_ml", MLEVEL
  call grib_new_from_samples(igrib, FILENAME)

  call grib_get(igrib,'numberOfPointsAlongAMeridian', NGJ)
  ALLOCATE (MLAT(NGJ))

  !     get as a integer
  call grib_get(igrib,'pl', MLAT)
  print*,'Number of points along a meridian  = ', NGJ
  NGI=SUM(MLAT)

  call grib_get(igrib,'numberOfVerticalCoordinateValues',mk)
  print*, 'numberOfVerticalCoordinateValues = ', mk

  IF(mk/2-1 .ne. MLEVEL) THEN 
     WRITE(*,*) 'FATAL: Number of model levels',mk, &
          ' does not agree with', MLEVEL,' in namelist'
     STOP
  ENDIF

  call grib_get_size(igrib,'pv',npv)
  allocate(pv(npv))
  call grib_get(igrib,'pv',pv)
  ALLOCATE(AK(NPV/2))
  ALLOCATE(BK(NPV/2))
  AK=pv(:NPV/2)
  BK=pv(NPV/2+1:)
  deallocate(pv)
! END GAUSS INFO



  ! Initialization of Legendre transform on LAT/LON grid.

  ALLOCATE (BREITE(MAXB))
  ALLOCATE (Z(0:((MNAUF+3)*(MNAUF+4))/2,MAXB))
  !$OMP PARALLEL DO
  DO  J=1,MAXB
     BREITE(J)=SIN((RLA1-(J-1.D0)*(RLA1-RLA0)/(MAXB-1))* PI/180.D0)
     CALL PLGNFA(MNAUF,BREITE(J),Z(0,J))
  ENDDO
  !$OMP END PARALLEL DO


  ! Initialisation of fields for  FFT and Legendre transformation
  !	to  Gaussian grid and back to phase space
  ALLOCATE (GBREITE(NGJ),WEIGHT(NGJ))
  ALLOCATE (P(0:((MNAUF+3)*(MNAUF+4))/2,NGJ/2))
  ALLOCATE (PP(NGJ/2,0:((MNAUF+3)*(MNAUF+4))/2))
  X1=-1.D0
  X2=1.D0
  CALL GAULEG(X1,X2,GBREITE,WEIGHT,NGJ)

  !$OMP PARALLEL DO PRIVATE(M)
  DO J=1,NGJ/2
     CALL PLGNFA(MNAUF,GBREITE(J),P(:,J))
     DO M=0,(MNAUF+3)*(MNAUF+4)/2
        PP(J,M)=P(M,J)
     ENDDO
  ENDDO
  !$OMP END PARALLEL DO


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   READ LNSP
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  FILENAME='fort.12'
  ALLOCATE (LNPMN(0:(MNAUF+1)*(MNAUF+2)-1))
  CALL READSPECTRAL(FILENAME,LNPMN,MNAUF,1,MLEVEL,(/152/))
  ! Call set99 an initialization routine that must be called once
  ! before a sequence of calls to the fft routines
  ALLOCATE (WSAVE(4*MAUF+15))
  ALLOCATE (PS(MAXL, MAXB,1))
  CALL SET99(WSAVE,IFAX,mauf)
  CALL PHGCUT(LNPMN,PS,WSAVE,IFAX,Z, &
       MNAUF,MNAUF,MAUF,MANF,MAXL,MAXB,1)
  CALL STATIS(MAXL,MAXB,1,EXP(PS),RMS,MW,SIG)
  WRITE(*,'(A12,3F12.4)') 'STATISTICS : ',RMS,MW,SIG

  ALLOCATE (GWSAVE(8*NGJ+15,NGJ/2))
  ALLOCATE (GIFAX(10,NGJ))
  DO J=1,NGJ/2
     CALL SET99(GWSAVE(1,J),GIFAX(1,J),MLAT(J))
  ENDDO
  ALLOCATE (PSG(NGI),HILF(NGI))
  ALLOCATE (LNPMN2(0:(MNAUF+1)*(MNAUF+2)-1))
  CALL PHGR213(LNPMN,HILF,GWSAVE,GIFAX,P,MLAT,MNAUF,NGI,NGJ,1)
  PSG=HILF
  CALL GRPH213(LNPMN2,PSG,GWSAVE,GIFAX,PP,WEIGHT,MLAT, &
       MNAUF,NGI,NGJ,1)
  CALL PHGR213(LNPMN2,HILF,GWSAVE,GIFAX,P,MLAT,MNAUF,NGI,NGJ,1)


  HILF=exp(PSG)-exp(HILF)

  CALL STATIS(NGI,1,1,HILF,RMS,MW,SIG)
  WRITE(*,'(A12,3F11.4)') 'STATISTICS: ',RMS,MW,SIG

  PSG=EXP(PSG)
  HILF=PSG
  CALL STATIS(NGI,1,1,HILF,RMS,MW,SIG)
  WRITE(*,'(A12,3F11.4)') 'STATISTICS: ',RMS,MW,SIG

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   READ U/V
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  FILENAME='fort.10'
  ALLOCATE (XMN(0:(MNAUF+1)*(MNAUF+2)-1, 2*MLEVEL))
  CALL READSPECTRAL(FILENAME, &
       XMN,MNAUF,2*MLEVEL,MLEVEL,(/131,132/)) 
  ! Transform winds on gaussian grid
  ALLOCATE (UGVG(NGI, 2*MLEVEL))
  CALL PHGR213(XMN,UGVG,GWSAVE,GIFAX,P,MLAT,MNAUF,NGI,NGJ,2*MLEVEL)
  ALLOCATE (CUA(2,4,MLEVEL))
  ALLOCATE (CVA(2,4,MLEVEL))
  DO K=1,MLEVEL
     ! North Pole
     CALL JSPPOLE(XMN(:,K),1,MNAUF,.TRUE.,CUA(:,:,K))
     CALL JSPPOLE(XMN(:,MLEVEL+K),1,MNAUF,.TRUE.,CVA(:,:,K))
     ! South Pole
     CALL JSPPOLE(XMN(:,K),-1,MNAUF,.TRUE.,CUA(:,3:4,K))
     CALL JSPPOLE(XMN(:,MLEVEL+K),-1,MNAUF,.TRUE.,CVA(:,3:4,K))
  ENDDO

  IF(MSMOOTH .ne. 0) THEN
    DO K=1,2*MLEVEL
       CALL SPFILTER(XMN(:,K),MNAUF,MSMOOTH)
    ENDDO
  ENDIF
  ALLOCATE (UV(MAXL, MAXB, 2*MLEVEL))
  CALL PHGCUT(XMN,UV,WSAVE,IFAX,Z, &
       MNAUF,MNAUF,MAUF,MANF,MAXL,MAXB,2*MLEVEL)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   READ Divergence
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  FILENAME='fort.13'
  CALL READSPECTRAL(FILENAME,XMN,MNAUF,MLEVEL,MLEVEL,(/155/))
  ! Tranform horizontal divergence on gaussian grid
  ALLOCATE (DG(NGI,MLEVEL))
  CALL PHGR213(XMN,DG,GWSAVE,GIFAX,P,MLAT,MNAUF,NGI,NGJ,MLEVEL)


  ! compute gradient of LNSP (log of surface pressure) on gaussian grid 
  ALLOCATE (H(0:(MNAUF+2)*(MNAUF+3)/2))
  ALLOCATE (DPSDL(NGI,1),DPSDM(NGI,1))
  CALL PHGRAD(LNPMN,DPSDL,DPSDM,GWSAVE,GIFAX,P,H,MLAT,MNAUF,NGI,NGJ,1)

  ! Compute vertical weed spind on gaussian grid
  ALLOCATE (ETAG(NGI,MLEVEL))
  CALL CONTGL(HILF,DPSDL,DPSDM,DG,UGVG(:,1),UGVG(:,MLEVEL+1), &
       GBREITE,ETAG,MLAT,AK,BK,NGI,NGJ,MLEVEL)


  CALL GRPH213(XMN,ETAG,GWSAVE,GIFAX,PP,WEIGHT,MLAT, &
       MNAUF,NGI,NGJ,MLEVEL)
  IF(MSMOOTH .ne. 0) THEN
    DO K=1,MLEVEL
       CALL SPFILTER(XMN(:,K),MNAUF,MSMOOTH)
    ENDDO
  ENDIF
  ALLOCATE (ETA(MAXL,MAXB,MLEVEL))
  CALL PHGCUT(XMN,ETA,WSAVE,IFAX,Z,MNAUF,MNAUF,MAUF,MANF,MAXL,MAXB,MLEVEL)

  CALL GRPH213(XMN,HILF,GWSAVE,GIFAX,PP,WEIGHT,MLAT, MNAUF,NGI,NGJ,1)

  IF(MSMOOTH .ne. 0) CALL SPFILTER(XMN(:,1),MNAUF,MSMOOTH)
  ALLOCATE (DPSDT(MAXL, MAXB,1))
  CALL PHGCUT(XMN,DPSDT,WSAVE,IFAX,Z,MNAUF,MNAUF,MAUF,MANF,MAXL,MAXB,1)

  CALL STATIS(MAXL,MAXB,1,DPSDT,RMS,MW,SIG)
  WRITE(*,'(A12,3F11.4)') 'STATISTICS DPSDT: ',RMS,MW,SIG

  CALL GRPH213(XMN,PSG,GWSAVE,GIFAX,PP,WEIGHT,MLAT,MNAUF,NGI,NGJ,1)
  CALL PHGCUT(XMN,PS,WSAVE,IFAX,Z,MNAUF,MNAUF,MAUF,MANF,MAXL,MAXB,1)

  CALL STATIS(MAXL,MAXB,1,PS,RMS,MW,SIG)
  WRITE(*,'(A12,3F11.4)') 'STATISTICS: ',RMS,MW,SIG

114 DEALLOCATE(HILF,PSG,DPSDL,DPSDM,ETAG,DG,LNPMN)

  DEALLOCATE(PP,P,UGVG,MLAT,GBREITE,WEIGHT,GWSAVE,XMN)

  ! CREATE FILE VERTICAL.EC NEEDED BY POP MODEL

  open(21,file='VERTICAL.EC')
  write(21,'(a)')
  write(21,'(a)') 'VERTICAL DISCRETIZATION OF POP MODEL'
  write(21,'(a)')
  write(21,'(i3,a)') MLEVEL,'   number of layers'
  write(21,'(a)')
  write(21,'(a)') '* A(NLEV+1)'
  write(21,'(a)')
  do  i=1,MLEVEL+1
     write(21,'(f18.12)') AK(I)
  enddo
  write(21,'(a)')
  write(21,'(a)') '* B(NLEV+1)'
  write(21,'(a)')
  do  i=1,MLEVEL+1
     write(21,'(f18.12)') BK(I)
  enddo
  close(21)


  !     Calculation of etadot in CONTGL needed scaled winds (ucosphi,vcosphi)
  !     Now we are transforming back to the usual winds. 
  ALLOCATE (HILFUV(2*MAXL,2))
  DO K=1,MLEVEL
     DO J=2,MAXB-1
        COSB=SQRT(1.0-(BREITE(J))*(BREITE(J)))
        UV(:,J,K)=UV(:,J,K)/COSB
        UV(:,J,MLEVEL+K)=UV(:,J,MLEVEL+K)/COSB
     ENDDO
     ! special treatment for poles, if necessary. 
     DO J=1,MAXB,MAXB-1
        COSB=SQRT(1.0-(BREITE(J))*(BREITE(J)))
        if(1.0-BREITE(J)*BREITE(J) .gt. 0 .OR. MGAUSS .NE. 1) then
           IF(RLA0 .EQ. -90.0 .AND. J .EQ. MAXB .OR. &
                RLA1 .EQ. 90.0 .AND. J .EQ. 1) then
              UV(:,J,K)=UV(:,J,K)*1.D6
              UV(:,J,MLEVEL+K)=UV(:,J,MLEVEL+K)*1.D6
           else
              UV(:,J,K)=UV(:,J,K)/COSB
              UV(:,J,MLEVEL+K)=UV(:,J,MLEVEL+K)/COSB
           endif
        else
           HILFUV(5:MAXL,:)=0.
           HILFUV(1:2,:)=0.
           IF(J.EQ.MAXB) THEN
              ! Suedpol
              HILFUV(3:4,1)=CUA(:,4,K)
              HILFUV(3:4,2)=CVA(:,4,K)
           ELSE
              ! Nordpol
              HILFUV(3:4,1)=CUA(:,2,K)
              HILFUV(3:4,2)=CVA(:,2,K)
           ENDIF
           CALL RFOURTR(HILFUV(:,1),WSAVE,IFAX,MAXL/2-1,MAXL)
           DO I=0,MAXL-1
              IF(MANF+I.LE.MAXL) THEN
                 UV(I+1,J,K)=HILFUV(MANF+I,1)
              ELSE
                 UV(I+1,J,K)=HILFUV(MANF-MAXL+I,1)
              ENDIF
           ENDDO
           CALL RFOURTR(HILFUV(:,2),WSAVE,IFAX,MAXL/2-1,MAXL)
           DO I=0,MAXL-1
              IF(MANF+I.LE.MAXL) THEN
                 UV(I+1,J,MLEVEL+K)=HILFUV(MANF+I,2)
              ELSE
                 UV(I+1,J,MLEVEL+K)=HILFUV(MANF-MAXL+I,2)
              ENDIF
           ENDDO
        endif
     ENDDO
  ENDDO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !                      READING OF T                      !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ALLOCATE (T(MAXL, MAXB, MLEVEL))
  FILENAME='fort.11'
  CALL READLATLON(FILENAME,T,MAXL,MAXB,MLEVEL,(/130/))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !                    WRITE MODEL LEVEL DATA TO fort.15            !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! open output file
  call grib_open_file(LUNIT,'fort.15','w')

  ! we use temperature on lat/lon on model levels as template for model level data
  LUNIT_OUT=0
  call grib_open_file(LUNIT_OUT,'fort.11','r')
  call grib_new_from_file(LUNIT_OUT,igrib, iret)
  call grib_close_file(LUNIT_OUT)


  CALL WRITELATLON(LUNIT,igrib,UV(:,:,1),MAXL,MAXB,MLEVEL,MLEVELIST,131)

  CALL WRITELATLON(LUNIT,igrib,UV(:,:,MLEVEL+1),MAXL,MAXB,MLEVEL,MLEVELIST,132)

  CALL WRITELATLON(LUNIT,igrib,ETA,MAXL,MAXB,MLEVEL,MLEVELIST,METAPAR)

  CALL WRITELATLON(LUNIT,igrib,T,MAXL,MAXB,MLEVEL,MLEVELIST,130)

  CALL WRITELATLON(LUNIT,igrib,PS,MAXL,MAXB,1,'1',134)

  call grib_close_file(LUNIT)


2000 STOP 'SUCCESSFULLY FINISHED CONVERT_PRE: CONGRATULATIONS'
3000 STOP 'ROUTINE CONVERT_PRE: ERROR'
9999 stop 'ROUTINE CONVERT_PRE: ERROR'

CONTAINS

  SUBROUTINE STATIS (NI,NJ,NK,PHI,RMS,MW,SIG)
    IMPLICIT NONE

    integer, intent(in) :: NI, NJ, NK

    REAL PHI(NI,NJ,NK),SIG,MW,RMS,P
    INTEGER N

    N=NI*NJ*NK

    RMS=0.
    MW=0.

    DO  I=1,NI
       DO  J=1,NJ
          DO  K=1,NK
             P=PHI(I,J,K)
             RMS=RMS+P*P
             MW=MW+P
          ENDDO
       ENDDO
    ENDDO

    RMS=SQRT(RMS/N)
    MW=MW/N

    IF(RMS*RMS-MW*MW.LT.0.) THEN	  
       SIG=0.0
    ELSE
       SIG=SQRT(RMS*RMS-MW*MW)
    ENDIF

    RETURN
  END SUBROUTINE STATIS

END PROGRAM PRECONVERT
