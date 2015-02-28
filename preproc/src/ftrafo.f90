MODULE FTRAFO

CONTAINS

  ! Berechnung des Gradienten eines Skalars aus dem Feld des
  !	Skalars XMN im Phasenraum. Zurueckgegeben werden die Felder der
  !	Komponenten des horizontalen Gradienten XLAM,XPHI auf dem Gauss'schen Gitter.
  !	GWSAVE ist ein Hilfsfeld fuer die FFT
  !	P enthaelt die assoziierten Legendrepolynome, H deren Ableitung
  !	MLAT enthaelt die Anzahl der Gitterpunkte pro Breitenkreis
  !	MNAUF gibt die spektrale Aufloesung an, 
  !	NI = Anzahl der Gauss'schen Gitterpunkte,
  !	NJ = Anzahl der Gauss'schen Breiten,
  !	NK = Anzahl der Niveaus

  SUBROUTINE PHGRAD(XMN,XLAM,XPHI,GWSAVE,IFAX,P,H,MLAT, &
       MNAUF,NI,NJ,NK)

    USE PHTOGR
    IMPLICIT NONE
    INTEGER   J,K,M,N,NI,NJ,NK,MNAUF,GGIND,LL,LLP,LLH,LLS,LLPS,LLHS
    INTEGER		MLAT(NJ),IFAX(10,NJ)
    REAL     	UFOUC(0:MAXAUF),MUFOUC(0:MAXAUF)
    REAL 			VFOUC(0:MAXAUF),MVFOUC(0:MAXAUF)
    REAL			XMN(0:(MNAUF+1)*(MNAUF+2)-1,NK)
    REAL		P(0:(MNAUF+3)*(MNAUF+4)/2,NJ/2)
    REAL		H(0:(MNAUF+2)*(MNAUF+3)/2)
    REAL			XLAM(NI,NK),XPHI(NI,NK)
    REAL			GWSAVE(8*NJ+15,NJ/2)
    REAL      ERAD
    REAL 		SCR,SCI,ACR,ACI,MUSCR,MUSCI,MUACR,MUACI,RT,IT

    ERAD = 6367470.0

    GGIND=0
    DO  J = 1,NJ/2
       CALL DPLGND(MNAUF,P(0,J),H)
       DO  K = 1,NK
	  LL=0
          LLP=0
	  LLH=0
	  DO  M = 0,MNAUF
             SCR=0.D0
             SCI=0.D0
             ACR=0.D0
             ACI=0.D0
             MUSCR=0.D0
             MUSCI=0.D0
             MUACR=0.D0
             MUACI=0.D0
             LLS=LL
             LLPS=LLP
             LLHS=LLH
             IF(2*M+1.LT.MLAT(J)) THEN
                DO  N = M,MNAUF,2
                   RT=XMN(2*LL,K)
                   IT=XMN(2*LL+1,K)
                   SCR =SCR+ RT*P(LLP,J)
                   SCI =SCI+ IT*P(LLP,J)
                   MUACR =MUACR+RT*H(LLH)
                   MUACI =MUACI+ IT*H(LLH)
                   LL=LL+2
                   LLP=LLP+2
                   LLH=LLH+2
                ENDDO
                LL=LLS+1
                LLP=LLPS+1
                LLH=LLHS+1
                DO N = M+1,MNAUF,2
                   RT=XMN(2*LL,K)
                   IT=XMN(2*LL+1,K)
                   ACR =ACR+ RT*P(LLP,J)
                   ACI =ACI+ IT*P(LLP,J)
                   MUSCR =MUSCR+ RT*H(LLH)
                   MUSCI =MUSCI+ IT*H(LLH)
                   LL=LL+2
                   LLP=LLP+2
                   LLH=LLH+2
                ENDDO
             ENDIF
             LL=LLS+(MNAUF-M+1)
             LLP=LLPS+(MNAUF-M+3)
             LLH=LLHS+(MNAUF-M+2)

             UFOUC(2*M)=-M*(SCI-ACI)/ERAD
             UFOUC(2*M+1)=M*(SCR-ACR)/ERAD
             VFOUC(2*M)=-M*(SCI+ACI)/ERAD
             VFOUC(2*M+1)=M*(SCR+ACR)/ERAD

             MUFOUC(2*M)=-(MUSCR-MUACR)/ERAD
             MUFOUC(2*M+1)=-(MUSCI-MUACI)/ERAD
             MVFOUC(2*M)=-(MUSCR+MUACR)/ERAD
             MVFOUC(2*M+1)=-(MUSCI+MUACI)/ERAD
          ENDDO

          CALL RFOURTR(VFOUC,GWSAVE(:,J),IFAX(:,J),MNAUF,MLAT(J))
          XLAM(GGIND+1:GGIND+MLAT(J),K)=VFOUC(0:MLAT(J)-1)
          CALL RFOURTR(UFOUC,GWSAVE(:,J),IFAX(:,J),MNAUF,MLAT(J))
          XLAM(NI-GGIND-MLAT(J)+1:NI-GGIND,K)=UFOUC(0:MLAT(J)-1)

          CALL RFOURTR(MVFOUC, GWSAVE(:,J),IFAX(:,J),MNAUF,MLAT(J))
          XPHI(GGIND+1:GGIND+MLAT(J),K)=MVFOUC(0:MLAT(J)-1)
          CALL RFOURTR(MUFOUC,GWSAVE(:,J),IFAX(:,J),MNAUF,MLAT(J))
          XPHI(NI-GGIND-MLAT(J)+1:NI-GGIND,K)=MUFOUC(0:MLAT(J)-1)

       ENDDO
       GGIND=GGIND+MLAT(J)
    ENDDO

  END SUBROUTINE PHGRAD

  ! Berechnung der Divergenz aus dem Windfeld (U,V)
  !	im Phasenraum. Zurueckgegeben werden die Felder der
  !	Komponenten des horizontalen Gradienten XLAM,XPHI auf dem Gauss'schen Gitter.
  !	GWSAVE ist ein Hilfsfeld fuer die FFT
  !	P enthaelt die assoziierten Legendrepolynome, H deren Ableitung
  !	MLAT enthaelt die Anzahl der Gitterpunkte pro Breitenkreis
  !	MNAUF gibt die spektrale Aufloesung an, 
  !	NI = Anzahl der Gauss'schen Gitterpunkte,
  !	NJ = Anzahl der Gauss'schen Breiten,
  !	NK = Anzahl der Niveaus
  ! Beachte, dass das Windfeld eine um 1 erhoehte Aufloesung in mu-Richtung hat.

  SUBROUTINE CONTGL(PS,DPSDL,DPSDM,DIV,U,V,BREITE,ETA,MLAT,A,B,NI,NJ,NK)

    IMPLICIT NONE

    INTEGER NI,NJ,NK,I,J,K,MLAT(NJ),L

    REAL A(NK+1),B(NK+1)
    REAL PS(NI),DPSDL(NI),DPSDM(NI)
    REAL DIV(NI,NK),U(NI,NK),V(NI,NK),ETA(NI,NK)
    REAL BREITE(NJ)

    REAL DIVT1,DIVT2,POB,PUN,DPSDT,COSB

    L=0
    DO  J=1,NJ
       COSB=(1.0-BREITE(J)*BREITE(J))
       DO  I=1,MLAT(J)
          L=L+1
          DIVT1=0.0
          DIVT2=0.0
          DO  K=1,NK
             POB=A(K)+B(K)*PS(L)
             PUN=A(K+1)+B(K+1)*PS(L)

             DIVT1=DIVT1+DIV(L,K)*(PUN-POB)
             if(cosb .gt. 0.) then
                DIVT2=DIVT2+(B(K+1)-B(K))*PS(L)* &
                (U(L,K)*DPSDL(L)+V(L,K)*DPSDM(L))/COSB
             endif

             ETA(L,K)=-DIVT1-DIVT2
          ENDDO

          DPSDT=(-DIVT1-DIVT2)/PS(L)

          DO  K=1,NK
             ETA(L,K)=ETA(L,K)-DPSDT*B(K+1)*PS(L)
          ENDDO
          PS(L)=DPSDT*PS(L)
       ENDDO
    ENDDO
  END SUBROUTINE CONTGL

END MODULE FTRAFO
