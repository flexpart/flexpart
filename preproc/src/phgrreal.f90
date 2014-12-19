MODULE PHTOGR

  INTEGER, PARAMETER :: MAXAUF=36000

CONTAINS

  SUBROUTINE PHGR213(CXMN,FELD,WSAVE,IFAX,Z,MLAT,MNAUF, &
    &MAXL,MAXB,MLEVEL)

    !     DIE ROUTINE F]HRT EINE TRANSFORMATION EINER
    !     FELDVARIABLEN VOM PHASENRAUM IN  DEN PHYSIKALISCHEN
    !     RAUM AUF DAS REDUZIERTE GAUSS'SCHE GITTER DURCH
    !
    !     CXMN   = SPEKTRALKOEFFIZIENTEN IN DER REIHENFOLGE
    !              CX00,CX01,CX11,CX02,....CXMNAUFMNAUF
    !     FELD   = FELD DER METEOROLOGISCHEN VARIABLEN
    !	WSAVE  = Working Array fuer Fouriertransformation
    !     Z 	 = LEGENDREFUNKTIONSWERTE
    !
    !     MNAUF    ANZAHL DER FOURIERKOEFFIZIENTEN
    !     MAXL     ANZAHL DER FUER DAS GITTER BENUTZTEN LAENGEN
    !     MAXB     ANZAHL DER FUER DAS GITTER BENOETIGTEN BREITEN
    !     MLEVEL   ANZAHL DER LEVELS, DIE TRANSFORMIERT WERDEN
    !
    IMPLICIT NONE

    !			Anzahl der Gitterpunkte auf jedem Breitenkreis
    INTEGER MLAT(MAXB/2)
    INTEGER K,MAXL,MAXB,MLEVEL,MNAUF
    INTEGER IND(MAXB)


    !    FELD DER LEGENDREPOLYNOME FUER EINE BREITE
    REAL Z(0:((MNAUF+3)*(MNAUF+4))/2,MAXB/2)

    REAL, INTENT(IN) ::  CXMN(0:(MNAUF+1)*(MNAUF+2)-1,MLEVEL)
    REAL FELD(MAXL,MLEVEL)
    REAL WSAVE(8*MAXB+15,MAXB/2)
    INTEGER :: IFAX(10,MAXB)

    IND(1)=0
    DO  K=2,MAXB/2
       IND(K)=IND(K-1)+MLAT(K-1)
    ENDDO

    !$OMP PARALLEL DO SCHEDULE(DYNAMIC)
    DO  K=1,MAXB/2
       CALL PHSYM(K,IND,CXMN,FELD,Z,WSAVE,IFAX,MLAT, &
       &MNAUF,MAXL,MAXB,MLEVEL)

    ENDDO
    !$OMP END PARALLEL DO

  END SUBROUTINE PHGR213
  !
  !
  SUBROUTINE PHSYM(K,IND,CXMN,FELD,Z,WSAVE,IFAX,MLAT, &
       &MNAUF,MAXL,MAXB,MLEVEL)

    IMPLICIT NONE

    INTEGER MLAT(MAXB/2)
    INTEGER K,L,I,J,LLS,LLPS,LL,LLP,MAXL,MAXB,MLEVEL,MNAUF
    INTEGER IND(MAXB)
    INTEGER :: IFAX(10,MAXB)


    !    FELD DER FOURIERKOEFFIZIENTEN
    REAL :: CXMS(0:MAXAUF-1),CXMA(0:MAXAUF-1)

    !    FELD DER LEGENDREPOLYNOME FUER EINE BREITE
    REAL Z(0:((MNAUF+3)*(MNAUF+4))/2,MAXB/2)
    REAL ACR,ACI,SCR,SCI

    REAL, INTENT(IN) ::  CXMN(0:(MNAUF+1)*(MNAUF+2)-1,MLEVEL)
    REAL FELD(MAXL,MLEVEL)
    REAL WSAVE(8*MAXB+15,MAXB/2)

    DO  L=1,MLEVEL
       LL=0
       LLP=0
       DO  I=0,MNAUF
          SCR=0.D0
          SCI=0.D0
          ACR=0.D0
          ACI=0.D0
          LLS=LL
          LLPS=LLP
          IF(2*I+1.LT.MLAT(K)) THEN
             !	Innerste Schleife aufgespalten um if-Abfrage zu sparen
             DO J=I,MNAUF,2
                SCR=SCR+Z(LLP,K)*CXMN(2*LL,L)
                SCI=SCI+Z(LLP,K)*CXMN(2*LL+1,L)
                LL=LL+2            
                LLP=LLP+2            
             ENDDO
             LL=LLS+1
             LLP=LLPS+1
             DO  J=I+1,MNAUF,2
                ACR=ACR+Z(LLP,K)*CXMN(2*LL,L)
                ACI=ACI+Z(LLP,K)*CXMN(2*LL+1,L)
                LL=LL+2            
                LLP=LLP+2            
             ENDDO
          ENDIF
          LL=LLS+(MNAUF-I+1)
          LLP=LLPS+(MNAUF-I+3)
          CXMS(2*I)=SCR+ACR
          CXMS(2*I+1)=SCI+ACI
          CXMA(2*I)=SCR-ACR
          CXMA(2*I+1)=SCI-ACI
       ENDDO
       CALL RFOURTR(CXMS,WSAVE(:,K),IFAX(:,K),MNAUF,&
            &MLAT(K))
       FELD(IND(k)+1:IND(K)+MLAT(K),L)=CXMS(0:MLAT(K)-1)
       CALL RFOURTR(CXMA,&
            &WSAVE(:,K),IFAX(:,K),MNAUF,MLAT(K))
       FELD(MAXL-IND(k)-MLAT(K)+1:MAXL-IND(k),L)=CXMA(0:MLAT(K)-1)

    ENDDO

  END SUBROUTINE PHSYM

  SUBROUTINE PHGCUT(CXMN,FELD,WSAVE,IFAX,Z,&
       &                  MNAUF,MMAX,MAUF,MANF,MAXL,MAXB,MLEVEL)

    !     DIE ROUTINE FUEHRT EINE TRANSFORMATION EINER
    !     FELDVARIABLEN VOM PHASENRAUM IN  DEN PHYSIKALISCHEN
    !     RAUM AUF KUGELKOORDINATEN DURCH. Es kann ein Teilausschnitt
    !			Der Erde angegeben werden. Diese Routine ist langsamer als
    !			phgrph
    !
    !     CXMN   = SPEKTRALKOEFFIZIENTEN IN DER REIHENFOLGE
    !              CX00,CX01,CX11,CX02,....CXMNAUFMNAUF
    !     FELD   = FELD DER METEOROLOGISCHEN VARIABLEN
    !     Z      = SINUS DER GEOGRAFISCHEN BREITEN
    !
    !     MNAUF    ANZAHL DER FOURIERKOEFFIZIENTEN
    !     MAUF     ANZAHL DER LAENGEN UND DER FOURIERKOEFFIZIENTEN
    !     MANF     ANFANG DES LAENGENBEREICHS FUER DAS GITTER,
    !              AUF DAS INTERPOLIERT WERDEN SOLL
    !     MAXL     ANZAHL DER FUER DAS GITTER BENUTZTEN LAENGEN
    !     MAXB     ANZAHL DER FUER DAS GITTER BENOETIGTEN BREITEN
    !     MLEVEL   ANZAHL DER LEVELS, DIE TRANSFORMIERT WERDEN
    !
    IMPLICIT NONE
    INTEGER, INTENT(IN)  :: MNAUF, MMAX, MAUF, MANF, MAXL, MAXB, MLEVEL
    REAL, INTENT(OUT)    :: FELD(MAXL,MAXB,MLEVEL)
    INTEGER, INTENT(INOUT) :: IFAX(10)

    !    FELD DER LEGENDREPOLYNOME FUER EINE BREITE
    REAL, INTENT(IN)     :: Z(0:((MMAX+3)*(MMAX+4))/2,MAXB)
    REAL, INTENT(IN)     :: CXMN(0:(MMAX+1)*(MMAX+2)-1,MLEVEL)
    REAL, INTENT(INOUT)  :: WSAVE(4*MAUF+15)  ! work array previously
                                              ! initialized with set99 
                                              ! (will be used in fft99)
    INTEGER              :: J
    LOGICAL              :: SYM


    IF(MAUF.LE.MNAUF) WRITE(*,*) 'TOO COARSE LONGITUDE RESOLUTION'
    IF((MANF.LT.1).OR.(MAXL.LT.1).OR. &
         &       (MANF.GT.MAUF).OR.(MAXL.GT.MAUF)) THEN
       WRITE(*,*) 'WRONG LONGITUDE RANGE',MANF,MAXL
       STOP
    ENDIF

    ! Pruefe, ob Ausgabegitter symmetrisch zum Aequator ist
    ! Wenn ja soll Symmetrie der Legendrepolynome ausgenutzt werden
    ! Check if the output grid is symmetrical to the equator
    ! If yes we exploit it

    IF(MAXB .GT. 4) THEN
       SYM=.TRUE.
       DO  J=5,5
          IF(ABS(ABS(Z(100,J))-ABS(Z(100,MAXB+1-J))).GT.1E-11) THEN
             SYM=.FALSE.
          ENDIF
       ENDDO
       WRITE(*,*) 'Symmetrisch: ',SYM
    ELSE
       SYM=.FALSE.
    ENDIF


    IF(SYM) THEN
       !$OMP PARALLEL DO 
       DO J=1,(MAXB+1)/2
          CALL PHSYMCUT(J,CXMN,FELD,Z,WSAVE,IFAX, &
               &MAUF,MNAUF,MAXL,MAXB,MLEVEL,MANF)

       ENDDO
       !$OMP END PARALLEL DO
    ELSE
       !$OMP PARALLEL DO 
       DO J=1,MAXB
          CALL PHGPNS(CXMN,FELD,Z,WSAVE,IFAX, &
               &J,MNAUF,MAUF,MANF,MAXL,MAXB,MLEVEL)
       ENDDO
       !$OMP END PARALLEL DO

    ENDIF


    RETURN
  END SUBROUTINE PHGCUT

  SUBROUTINE PHSYMCUT(J,CXMN,FELD,Z,WSAVE,IFAX,&
       &MAUF,MNAUF,MAXL,MAXB,MLEVEL,MANF)

    IMPLICIT NONE

    INTEGER, INTENT(IN)     :: J
    REAL,    INTENT(IN)     :: CXMN(0:(MNAUF+1)*(MNAUF+2)-1,MLEVEL)
    REAL,    INTENT(OUT)    :: FELD(MAXL,MAXB,MLEVEL)
    REAL,    INTENT(IN)     :: Z(0:((MNAUF+3)*(MNAUF+4))/2,MAXB)
    REAL,    INTENT(INOUT)  :: WSAVE(4*MAUF+15)
    INTEGER, INTENT(IN)     :: IFAX(10)
    INTEGER, INTENT(IN)     :: MAUF,MNAUF,MAXL,MAXB,MLEVEL,MANF
    
    !    FELD DER FOURIERKOEFFIZIENTEN

    REAL :: CXM(0:MAXAUF-1),CXMA(0:MAXAUF-1)

    !    FELD DER LEGENDREPOLYNOME FUER EINE BREITE
    REAL SCR,SCI,ACR,ACI


    INTEGER            :: L, LL, LLP, LLS, LLPS, I, K

    DO L=1,MLEVEL
       LL=0
       LLP=0
       DO  I=0,MNAUF
          SCR=0.D0
          SCI=0.D0
          ACR=0.D0
          ACI=0.D0
          LLS=LL
          LLPS=LLP
          !	Innerste Schleife aufgespalten um if-Abfrage zu sparen
          DO  K=I,MNAUF,2
             SCR=SCR+Z(LLP,J)*CXMN(2*LL,L)
             SCI=SCI+Z(LLP,J)*CXMN(2*LL+1,L)
             LL=LL+2            
             LLP=LLP+2            
          ENDDO
          LL=LLS+1
          LLP=LLPS+1
          DO  K=I+1,MNAUF,2
             ACR=ACR+Z(LLP,J)*CXMN(2*LL,L)
             ACI=ACI +Z(LLP,J)*CXMN(2*LL+1,L)
             LL=LL+2            
             LLP=LLP+2            
          ENDDO
          LL=LLS+MNAUF-I+1
          LLP=LLPS+MNAUF-I+3
          CXM(2*I)=SCR+ACR
          CXM(2*I+1)=SCI+ACI
          CXMA(2*I)=SCR-ACR
          CXMA(2*I+1)=SCI-ACI
       ENDDO

       CALL RFOURTR(CXM,WSAVE,IFAX,MNAUF,MAUF)
       DO  I=0,MAXL-1
          IF(MANF+I.LE.MAUF) THEN
             FELD(I+1,J,L)=CXM(MANF+I-1)
          ELSE
             FELD(I+1,J,L)=CXM(MANF-MAUF+I-1)
          ENDIF
       ENDDO
       CALL RFOURTR(CXMA,WSAVE,IFAX,MNAUF,MAUF)
       DO  I=0,MAXL-1
          IF(MANF+I.LE.MAUF) THEN
             FELD(I+1,MAXB+1-J,L)=CXMA(MANF+I-1)
          ELSE
             FELD(I+1,MAXB+1-J,L)=CXMA(MANF-MAUF+I-1)
          ENDIF
       ENDDO
    ENDDO

  END SUBROUTINE PHSYMCUT

  SUBROUTINE PHGPNS(CXMN,FELD,Z,WSAVE,IFAX, &
       J,MNAUF,MAUF,MANF,MAXL,MAXB,MLEVEL)

    IMPLICIT NONE
    REAL,    intent(in)    :: CXMN(0:(MNAUF+1)*(MNAUF+2)-1,MLEVEL)
    REAL,    intent(out)   :: FELD(MAXL,MAXB,MLEVEL)
    REAL,    intent(in)    :: Z(0:((MNAUF+3)*(MNAUF+4))/2,MAXB)
    REAL,    intent(inout) :: WSAVE(4*MAUF+15)
    INTEGER, intent(in)    :: IFAX(10)
    INTEGER, intent(in)    :: MNAUF,MAUF,MANF,J,MAXL,MAXB,MLEVEL

    REAL :: CXM(0:MAXAUF-1)

    INTEGER I,L

    DO L=1,MLEVEL
       CALL LEGTR(CXMN(:,L),CXM,Z(:,J),MNAUF,MAUF)
       CALL RFOURTR(CXM,WSAVE,IFAX,MNAUF,MAUF)

       DO I=0,MAXL-1
          IF(MANF+I.LE.MAUF) THEN
             FELD(I+1,J,L)=CXM(MANF+I-1)
          ELSE
             FELD(I+1,J,L)=CXM(MANF-MAUF+I-1)
          ENDIF
       ENDDO
    ENDDO
  END SUBROUTINE PHGPNS
  !
  SUBROUTINE LEGTR(CXMN,CXM,Z,MNAUF,MAUF)
    IMPLICIT NONE
    INTEGER MNAUF,MAUF,LL,LLP,I,J
    REAL CXM(0:MAXAUF-1)
    REAL CXMN(0:(MNAUF+1)*(MNAUF+2)-1)
    REAL Z(0:((MNAUF+3)*(MNAUF+4))/2)
    REAL CI,CR
    !
    !     DIESE ROUTINE BERECHNET DIE FOURIERKOEFFIZIENTEN CXM
    !
    LL=0
    LLP=0
    DO  I=0,MNAUF
       CR=0.D0
       CI=0.D0
       DO  J=I,MNAUF
          CR=CR+Z(LLP)*CXMN(2*LL)
          CI=CI+Z(LLP)*CXMN(2*LL+1)
          LL=LL+1
          LLP=LLP+1
       ENDDO
       LLP=LLP+2
       CXM(2*I)=CR
       CXM(2*I+1)=CI
    ENDDO

  END SUBROUTINE LEGTR
  !
  !     
  SUBROUTINE RFOURTR(CXM,TRIGS,IFAX,MNAUF,MAXL)
    !     BERECHNET DIE FOURIERSUMME MIT EINEM FFT-ALGORITHMUS
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: MNAUF, MAXL
    INTEGER,intent(IN)  :: IFAX(10)
    REAL, intent(inout) :: CXM(0:MAXAUF-1)
    REAL                :: WSAVE(2*MAXL),TRIGS(2*MAXL)

    INTEGER :: I

    DO I=MNAUF+1,MAXL-1
       CXM(2*I)=0.0
       CXM(2*I+1)=0.0
    ENDDO
    CALL FFT99(CXM,WSAVE,TRIGS,IFAX,1,1,MAXL,1,1)
    DO I=0,MAXL-1
       CXM(I)=CXM(I+1)
    ENDDO

    RETURN
  END SUBROUTINE RFOURTR
  !     
  !    
  SUBROUTINE GAULEG(X1,X2,X,W,N)
  ! From numerical recipes
  ! Given the lower and upper limits of integration X1 and X2, 
  ! this routine returns arrayx X and W of length N containing the
  ! abscissa and weights of the Gauss-Legendre N-point quadrature
  ! formula. The parameter EPS is the relative precision. Note that
  ! internal computations are done in double precision.
    IMPLICIT NONE

    INTEGER, INTENT(IN)             :: N
    REAL, INTENT(IN)                :: X1, X2
    REAL, DIMENSION(N), INTENT(OUT) :: X, W

    REAL                            :: PI
    REAL, PARAMETER                 :: EPS=3.D-14
    INTEGER                         :: I,J,M 
    REAL                            :: P1, P2, P3
    REAL                            :: PP
    REAL                            :: Z1, Z
    REAL                            :: XM, XL

    PI=ACOS(-1.D0)
   
    M=(N+1)/2
    XM=0.5D0*(X2+X1)
    XL=0.5D0*(X2-X1)
    DO  I=1,M
       Z=DCOS(3.141592654D0*(I-.25D0)/(N+.5D0))
       DO 
          P1=1.D0
          P2=0.D0
          DO  J=1,N
             P3=P2
             P2=P1
             P1=((2.D0*J-1.D0)*Z*P2-(J-1.D0)*P3)/J
          ENDDO
          PP=N*(Z*P1-P2)/(Z*Z-1.D0)
          Z1=Z
          Z=Z1-P1/PP
          IF (ABS(Z-Z1).LE.EPS) EXIT
       ENDDO

       X(I)=XM-XL*Z
       X(N+1-I)=XM+XL*Z
       W(I)=2.D0*XL/((1.D0-Z*Z)*PP*PP)
       W(N+1-I)=W(I)
    ENDDO
  END SUBROUTINE GAULEG
  !
  !
  SUBROUTINE PLGNFA(LL,X,Z)
    !
    ! PLGNFA is associated with normalization of 
    ! legendre functions OF P00(X) to PLL(X)
    ! and write result in Z
    ! The polynomials are indexed as the ECMWF, i.e.
    ! P00, P10, P11, P20, P21, P22, ... 
    ! This routine is analogous to PLGNDN
    ! X is the cosinus of the zenith angle or
    ! sinus of latitude
    !
    IMPLICIT NONE

    INTEGER, INTENT(IN)                                :: LL
    REAL, INTENT(IN)                                   :: X
    REAL, DIMENSION (0:((LL+3)*(LL+4))/2), INTENT(OUT) :: Z

    REAL                 :: FACT, POT, SOMX2
    REAL                 :: DJ, DK, DDK
    INTEGER              :: I, J, K, L 

    !
    L=LL+2
    I=1
    Z(0)=1.D0
    FACT=1.D0
    POT=1.D0
    SOMX2=DSQRT(1.D0-X*X)
    DO  J=0,L
       DJ=DBLE(J)
       IF(J.GT.0) THEN
          FACT=FACT*(2.D0*DJ-1.D0)/(2.D0*DJ)
          POT=POT*SOMX2
          Z(I)=DSQRT((2.D0*DJ+1.D0)*FACT)*POT
          I=I+1
       ENDIF
       IF(J.LT.L) THEN
          Z(I)=X* &
               &DSQRT((4.D0*DJ*DJ+8.D0*DJ+3.D0)/(2.D0*DJ+1.D0))*Z(I-1)
          I=I+1
       ENDIF
       DK=DJ+2.D0
       DO  K=J+2,L
          DDK=(DK*DK-DJ*DJ)
          Z(I)=X*DSQRT((4.D0*DK*DK-1.D0)/DDK)*Z(I-1)-       &
               &    DSQRT(((2.D0*DK+1.D0)*(DK-DJ-1.D0)*(DK+DJ-1.D0))/ &
               &    ((2.D0*DK-3.D0)*DDK))*Z(I-2)
          DK=DK+1.D0
          I=I+1
       ENDDO
    ENDDO
    RETURN
  END SUBROUTINE PLGNFA


  SUBROUTINE DPLGND(MNAUF,Z,DZ)
    !
    ! DPLGND BERECHNET DIE ABLEITUNG DER NORMIERTEN ASSOZIIERTEN
    ! LEGENDREFUNKTIONEN VON P00(X) BIS PLL(X)
    ! UND SCHREIBT SIE IN DAS FELD DZ
    ! DIE REIHENFOLGE IST
    ! P00(X),P01(X),P11(X),P02(X),P12(X),P22(X),..PLL(X)
    !
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: MNAUF
    REAL Z(0:((MNAUF+3)*(MNAUF+4))/2)
    REAL DZ(0:((MNAUF+2)*(MNAUF+3))/2)

    INTEGER   :: I, J, LLP, LLH
    REAL      :: WURZELA, WURZELB
    !
    IF(Z(0).NE.1.D0) THEN
       WRITE(*,*) 'DPLGND: Z(0) must be 1.0'
       STOP
    ENDIF

    LLP=0
    LLH=0
    DO  I=0,MNAUF+1
       DO  J=I,MNAUF+1
          IF(I.EQ.J) THEN
             WURZELA=DSQRT(DBLE((J+1)*(J+1)-I*I)/DBLE(4*(J+1)*(J+1)-1))
             DZ(LLH)=DBLE(J)*WURZELA*Z(LLP+1)
          ELSE
             WURZELB=DSQRT(DBLE((J+1)*(J+1)-I*I)/DBLE(4*(J+1)*(J+1)-1))
             DZ(LLH)=DBLE(J)*WURZELB*Z(LLP+1)-DBLE(J+1)*WURZELA*Z(LLP-1)
             WURZELA=WURZELB
          ENDIF
          LLH=LLH+1              
          LLP=LLP+1
       ENDDO
       LLP=LLP+1
    ENDDO
  END SUBROUTINE DPLGND


  ! Spectral Filter of Sardeshmukh and Hoskins (1984, MWR)
  ! MM=Spectral truncation of field
  ! MMAX= Spectral truncation of filter
  !
  SUBROUTINE SPFILTER(FELDMN,MM,MMAX)

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: MM,MMAX
    REAL, INTENT(INOUT) :: FELDMN(0:(MM+1)*(MM+2)-1)
! local variables
    INTEGER             :: I,J,K,L
    REAL                :: KMAX,SMAX,FAK

    SMAX=0.1
    KMAX=-ALOG(SMAX)
    KMAX=KMAX/(float(MMAX)*float(MMAX+1))**2
    l=0
    do i=0,MM
       do j=i,MM
          if(j .le. MMAX) then
             fak=1.0
             feldmn(2*l)=feldmn(2*l)*fak
             feldmn(2*l+1)=feldmn(2*l+1)*fak
          else
             feldmn(2*l)=0.
             feldmn(2*l+1)=0.
          endif
          l=l+1
       enddo
    enddo
  END SUBROUTINE SPFILTER

END MODULE PHTOGR
    


































