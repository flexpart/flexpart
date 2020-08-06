! SPDX-FileCopyrightText: FLEXPART 1998-2019, see flexpart_license.txt
! SPDX-License-Identifier: GPL-3.0-or-later

!**************************************************************************
!****                       SUBROUTINE CONVECT                        *****
!****                          VERSION 4.3c                           *****
!****                          20 May, 2002                           *****
!****                          Kerry Emanuel                          *****
!**************************************************************************
!
  SUBROUTINE CONVECT &
         (ND,  NL,   DELT, IFLAG, &
         PRECIP, WD,   TPRIME, QPRIME, CBMF    )
  !
  !-cv *************************************************************************
  !-cv C. Forster, November 2003 - May 2004:
  !-cv
  !-cv The subroutine has been downloaded from Kerry Emanuel's homepage,
  !-cv where further infos on the convection scheme can be found
  !-cv http://www-paoc.mit.edu/~emanuel/home.html
  !-cv
  !-cv The following changes have been made to integrate this subroutine
  !-cv into FLEXPART
  !-cv
  !-cv Putting most of the variables in a new common block
  !-cv renaming eps to eps0 because there is some eps already in includepar
  !-cv
  !-cv removing the arrays U,V,TRA and related arrays
  !-cv
  !-cv renaming the original arrays T,Q,QS,P,PH to
  !-cv TCONV,QCONV,QSCONV,PCONV_HPA,PHCONV_HPA
  !-cv
  !-cv Initialization of variables has been put into parameter statements
  !-cv instead of assignment of values at each call, in order to save 
  !-cv computation time.
  !***************************************************************************
  !
  !-----------------------------------------------------------------------------
  !    *** On input:      ***
  !
  !T:   Array of absolute temperature (K) of dimension ND, with first
  !      index corresponding to lowest model level. Note that this array
  !      will be altered by the subroutine if dry convective adjustment
  !      occurs and if IPBL is not equal to 0.
  !
  !Q:   Array of specific humidity (gm/gm) of dimension ND, with first
  !       index corresponding to lowest model level. Must be defined
  !       at same grid levels as T. Note that this array will be altered
  !       if dry convective adjustment occurs and if IPBL is not equal to 0.
  !
  !QS:  Array of saturation specific humidity of dimension ND, with first
  !       index corresponding to lowest model level. Must be defined
  !       at same grid levels as T. Note that this array will be altered
  !       if dry convective adjustment occurs and if IPBL is not equal to 0.
  !
  !U:   Array of zonal wind velocity (m/s) of dimension ND, witth first
  !       index corresponding with the lowest model level. Defined at
  !       same levels as T. Note that this array will be altered if
  !       dry convective adjustment occurs and if IPBL is not equal to 0.
  !
  !V:   Same as U but for meridional velocity.
  !
  !TRA: Array of passive tracer mixing ratio, of dimensions (ND,NTRA),
  !       where NTRA is the number of different tracers. If no
  !       convective tracer transport is needed, define a dummy
  !       input array of dimension (ND,1). Tracers are defined at
  !       same vertical levels as T. Note that this array will be altered
  !       if dry convective adjustment occurs and if IPBL is not equal to 0.
  !
  !P:   Array of pressure (mb) of dimension ND, with first
  !       index corresponding to lowest model level. Must be defined
  !       at same grid levels as T.
  !
  !PH:  Array of pressure (mb) of dimension ND+1, with first index
  !       corresponding to lowest level. These pressures are defined at
  !       levels intermediate between those of P, T, Q and QS. The first
  !       value of PH should be greater than (i.e. at a lower level than)
  !       the first value of the array P.
  !
  !ND:  The dimension of the arrays T,Q,QS,P,PH,FT and FQ
  !
  !NL:  The maximum number of levels to which convection can
  !       penetrate, plus 1.
  !       NL MUST be less than or equal to ND-1.
  !
  !NTRA:The number of different tracers. If no tracer transport
  !       is needed, set this equal to 1. (On most compilers, setting
  !       NTRA to 0 will bypass tracer calculation, saving some CPU.)
  !
  !DELT: The model time step (sec) between calls to CONVECT
  !
  !----------------------------------------------------------------------------
  !    ***   On Output:         ***
  !
  !IFLAG: An output integer whose value denotes the following:
  !
  !           VALUE                        INTERPRETATION
  !           -----                        --------------
  !             0               No moist convection; atmosphere is not
  !                             unstable, or surface temperature is less
  !                             than 250 K or surface specific humidity
  !                             is non-positive.
  !
  !             1               Moist convection occurs.
  !
  !             2               No moist convection: lifted condensation
  !                             level is above the 200 mb level.
  !
  !             3               No moist convection: cloud base is higher
  !                             then the level NL-1.
  !
  !             4               Moist convection occurs, but a CFL condition
  !                             on the subsidence warming is violated. This
  !                             does not cause the scheme to terminate.
  !
  !FT:   Array of temperature tendency (K/s) of dimension ND, defined at same
  !        grid levels as T, Q, QS and P.
  !
  !FQ:   Array of specific humidity tendencies ((gm/gm)/s) of dimension ND,
  !        defined at same grid levels as T, Q, QS and P.
  !
  !FU:   Array of forcing of zonal velocity (m/s^2) of dimension ND,
  !        defined at same grid levels as T.
  !
  !FV:   Same as FU, but for forcing of meridional velocity.
  !
  !FTRA: Array of forcing of tracer content, in tracer mixing ratio per
  !        second, defined at same levels as T. Dimensioned (ND,NTRA).
  !
  !PRECIP: Scalar convective precipitation rate (mm/day).
  !
  !WD:    A convective downdraft velocity scale. For use in surface
  !        flux parameterizations. See convect.ps file for details.
  !
  !TPRIME: A convective downdraft temperature perturbation scale (K).
  !         For use in surface flux parameterizations. See convect.ps
  !         file for details.
  !
  !QPRIME: A convective downdraft specific humidity
  !         perturbation scale (gm/gm).
  !         For use in surface flux parameterizations. See convect.ps
  !         file for details.
  !
  !CBMF:   The cloud base mass flux ((kg/m**2)/s). THIS SCALAR VALUE MUST
  !         BE STORED BY THE CALLING PROGRAM AND RETURNED TO CONVECT AT
  !         ITS NEXT CALL. That is, the value of CBMF must be "remembered"
  !         by the calling program between calls to CONVECT.
  !
  !-----------------------------------------------------------------------------
  !
  !    ***  THE PARAMETER NA SHOULD IN GENERAL BE GREATER THAN   ***
  !    ***                OR EQUAL TO  ND + 1                    ***
  !
  !
  use par_mod
  use conv_mod

  implicit none
  !
  !-cv====>Begin Module CONVECT    File convect.f      Undeclared variables
  !
  !Argument variables
  !
  integer :: iflag, nd, nl
  !
  real :: cbmf, delt, precip, qprime, tprime, wd
  !
  !Local variables
  !
  integer :: i, icb, ihmin, inb, inb1, j, jtt, k
  integer :: nk
  !
  real :: ad, afac, ahmax, ahmin, alt, altem
  real :: am, amp1, anum, asij, awat, b6, bf2, bsum, by
  real :: byp, c6, cape, capem, cbmfold, chi, coeff
  real :: cpinv, cwat, damps, dbo, dbosum
  real :: defrac, dei, delm, delp, delt0, delti, denom, dhdp
  real :: dpinv, dtma, dtmin, dtpbl, elacrit, ents
  real :: epmax, fac, fqold, frac, ftold
  real :: plcl, qp1, qsm, qstm, qti, rat
  real :: rdcp, revap, rh, scrit, sigt, sjmax
  real :: sjmin, smid, smin, stemp, tca
  real :: tvaplcl, tvpplcl, tvx, tvy, wdtrain

  !integer jc,jn
  !real alvnew,a2,ahm,alv,rm,sum,qnew,dphinv,tc,thbar,tnew,x

  real :: FUP(NA),FDOWN(NA)
  !
  !-cv====>End Module   CONVECT    File convect.f

  INTEGER :: NENT(NA)
  REAL :: M(NA),MP(NA),MENT(NA,NA),QENT(NA,NA),ELIJ(NA,NA)
  REAL :: SIJ(NA,NA),TVP(NA),TV(NA),WATER(NA)
  REAL :: QP(NA),EP(NA),TH(NA),WT(NA),EVAP(NA),CLW(NA)
  REAL :: SIGP(NA),TP(NA),CPN(NA)
  REAL :: LV(NA),LVCP(NA),H(NA),HP(NA),GZ(NA),HM(NA)
  !REAL TOLD(NA)
  !
  ! -----------------------------------------------------------------------
  !
  !   ***                     Specify Switches                         ***
  !
  !   ***   IPBL: Set to zero to bypass dry adiabatic adjustment       ***
  !   ***    Any other value results in dry adiabatic adjustment       ***
  !   ***     (Zero value recommended for use in models with           ***
  !   ***                   boundary layer schemes)                    ***
  !
  !   ***   MINORIG: Lowest level from which convection may originate  ***
  !   ***     (Should be first model level at which T is defined       ***
  !   ***      for models using bulk PBL schemes; otherwise, it should ***
  !   ***      be the first model level at which T is defined above    ***
  !   ***                      the surface layer)                      ***
  !
    INTEGER,PARAMETER :: IPBL=0
    INTEGER,PARAMETER :: MINORIG=1
  !
  !------------------------------------------------------------------------------
  !
  !   ***                    SPECIFY PARAMETERS                        ***
  !
  !   *** ELCRIT IS THE AUTOCONVERSION THERSHOLD WATER CONTENT (gm/gm) ***
  !   ***  TLCRIT IS CRITICAL TEMPERATURE BELOW WHICH THE AUTO-        ***
  !   ***       CONVERSION THRESHOLD IS ASSUMED TO BE ZERO             ***
  !   ***     (THE AUTOCONVERSION THRESHOLD VARIES LINEARLY            ***
  !   ***               BETWEEN 0 C AND TLCRIT)                        ***
  !   ***   ENTP IS THE COEFFICIENT OF MIXING IN THE ENTRAINMENT       ***
  !   ***                       FORMULATION                            ***
  !   ***  SIGD IS THE FRACTIONAL AREA COVERED BY UNSATURATED DNDRAFT  ***
  !   ***  SIGS IS THE FRACTION OF PRECIPITATION FALLING OUTSIDE       ***
  !   ***                        OF CLOUD                              ***
  !   ***        OMTRAIN IS THE ASSUMED FALL SPEED (P/s) OF RAIN       ***
  !   ***     OMTSNOW IS THE ASSUMED FALL SPEED (P/s) OF SNOW          ***
  !   ***  COEFFR IS A COEFFICIENT GOVERNING THE RATE OF EVAPORATION   ***
  !   ***                          OF RAIN                             ***
  !   ***  COEFFS IS A COEFFICIENT GOVERNING THE RATE OF EVAPORATION   ***
  !   ***                          OF SNOW                             ***
  !   ***     CU IS THE COEFFICIENT GOVERNING CONVECTIVE MOMENTUM      ***
  !   ***                         TRANSPORT                            ***
  !   ***    DTMAX IS THE MAXIMUM NEGATIVE TEMPERATURE PERTURBATION    ***
  !   ***        A LIFTED PARCEL IS ALLOWED TO HAVE BELOW ITS LFC      ***
  !   ***    ALPHA AND DAMP ARE PARAMETERS THAT CONTROL THE RATE OF    ***
  !   ***                 APPROACH TO QUASI-EQUILIBRIUM                ***
  !   ***   (THEIR STANDARD VALUES ARE  0.20 AND 0.1, RESPECTIVELY)    ***
  !   ***                   (DAMP MUST BE LESS THAN 1)                 ***
  !
    REAL,PARAMETER :: ELCRIT=.0011
    REAL,PARAMETER :: TLCRIT=-55.0
    REAL,PARAMETER :: ENTP=1.5
    REAL,PARAMETER :: SIGD=0.05
    REAL,PARAMETER :: SIGS=0.12
    REAL,PARAMETER :: OMTRAIN=50.0
    REAL,PARAMETER :: OMTSNOW=5.5
    REAL,PARAMETER :: COEFFR=1.0
    REAL,PARAMETER :: COEFFS=0.8
    REAL,PARAMETER :: CU=0.7
    REAL,PARAMETER :: BETA=10.0
    REAL,PARAMETER :: DTMAX=0.9
    REAL,PARAMETER :: ALPHA=0.025  !original 0.2
    REAL,PARAMETER :: DAMP=0.1
  !
  !   ***        ASSIGN VALUES OF THERMODYNAMIC CONSTANTS,        ***
  !   ***            GRAVITY, AND LIQUID WATER DENSITY.           ***
  !   ***             THESE SHOULD BE CONSISTENT WITH             ***
  !   ***              THOSE USED IN CALLING PROGRAM              ***
  !   ***     NOTE: THESE ARE ALSO SPECIFIED IN SUBROUTINE TLIFT  ***
  !
  REAL,PARAMETER :: CPD=1005.7
  REAL,PARAMETER :: CPV=1870.0
  REAL,PARAMETER :: CL=2500.0
  REAL,PARAMETER :: RV=461.5
  REAL,PARAMETER :: RD=287.04
  REAL,PARAMETER :: LV0=2.501E6
  REAL,PARAMETER :: G=9.81
  REAL,PARAMETER :: ROWL=1000.0
  !
  REAL,PARAMETER :: CPVMCL=CL-CPV
  REAL,PARAMETER :: EPS0=RD/RV
  REAL,PARAMETER :: EPSI=1./EPS0
  REAL,PARAMETER :: GINV=1.0/G
  REAL,PARAMETER :: EPSILON=1.e-20

  ! EPSILON IS A SMALL NUMBER USED TO EXCLUDE MASS FLUXES OF ZERO
  !
  DELTI=1.0/DELT
  !
  !      ***  INITIALIZE OUTPUT ARRAYS AND PARAMETERS  ***
  !

    DO I=1,NL+1
     FT(I)=0.0
     FQ(I)=0.0
     FDOWN(I)=0.0
     SUB(I)=0.0
     FUP(I)=0.0
     M(I)=0.0
     MP(I)=0.0
    DO J=1,NL+1
     FMASS(I,J)=0.0
     MENT(I,J)=0.0
    END DO
    END DO
    DO I=1,NL+1
     RDCP=(RD*(1.-QCONV(I))+QCONV(I)*RV)/ &
          (CPD*(1.-QCONV(I))+QCONV(I)*CPV)
     TH(I)=TCONV(I)*(1000.0/PCONV_HPA(I))**RDCP
    END DO
    PRECIP=0.0
    WD=0.0
    TPRIME=0.0
    QPRIME=0.0
    IFLAG=0
  !
  !  IF(IPBL.NE.0)THEN
  !
  !***            PERFORM DRY ADIABATIC ADJUSTMENT            ***
  !
  !  JC=0
  !  DO 30 I=NL-1,1,-1
  !   JN=0
  !    SUM=TH(I)*(1.+QCONV(I)*EPSI-QCONV(I))
  !   DO 10 J=I+1,NL
  !    SUM=SUM+TH(J)*(1.+QCONV(J)*EPSI-QCONV(J))
  !    THBAR=SUM/REAL(J+1-I)
  !    IF((TH(J)*(1.+QCONV(J)*EPSI-QCONV(J))).LT.THBAR)JN=J
  !  10    CONTINUE
  !   IF(I.EQ.1)JN=MAX(JN,2)
  !   IF(JN.EQ.0)GOTO 30
  !  12    CONTINUE
  !   AHM=0.0
  !   RM=0.0
  !   DO 15 J=I,JN
  !    AHM=AHM+(CPD*(1.-QCONV(J))+QCONV(J)*CPV)*TCONV(J)*
  !    +   (PHCONV_HPA(J)-PHCONV_HPA(J+1))
  !    RM=RM+QCONV(J)*(PHCONV_HPA(J)-PHCONV_HPA(J+1))
  !  15    CONTINUE
  !   DPHINV=1./(PHCONV_HPA(I)-PHCONV_HPA(JN+1))
  !   RM=RM*DPHINV
  !   A2=0.0
  !   DO 20 J=I,JN
  !    QCONV(J)=RM
  !    RDCP=(RD*(1.-QCONV(J))+QCONV(J)*RV)/
  !    1     (CPD*(1.-QCONV(J))+QCONV(J)*CPV)
  !    X=(0.001*PCONV_HPA(J))**RDCP
  !    TOLD(J)=TCONV(J)
  !    TCONV(J)=X
  !    A2=A2+(CPD*(1.-QCONV(J))+QCONV(J)*CPV)*X*
  !    1    (PHCONV_HPA(J)-PHCONV_HPA(J+1))
  !  20    CONTINUE
  !   DO 25 J=I,JN
  !    TH(J)=AHM/A2
  !    TCONV(J)=TCONV(J)*TH(J)
  !    TC=TOLD(J)-273.15
  !    ALV=LV0-CPVMCL*TC
  !    QSCONV(J)=QSCONV(J)+QSCONV(J)*(1.+QSCONV(J)*(EPSI-1.))*ALV*
  !    1    (TCONV(J)- TOLD(J))/(RV*TOLD(J)*TOLD(J))
  ! if (qslev(j) .lt. 0.) then
  !   write(*,*) 'qslev.lt.0 ',j,qslev
  ! endif
  !  25    CONTINUE
  !   IF((TH(JN+1)*(1.+QCONV(JN+1)*EPSI-QCONV(JN+1))).LT.
  !    1    (TH(JN)*(1.+QCONV(JN)*EPSI-QCONV(JN))))THEN
  !    JN=JN+1
  !    GOTO 12
  !   END IF
  !   IF(I.EQ.1)JC=JN
  !  30   CONTINUE
  !
  !   ***   Remove any supersaturation that results from adjustment ***
  !
  !IF(JC.GT.1)THEN
  ! DO 38 J=1,JC
  !    IF(QSCONV(J).LT.QCONV(J))THEN
  !     ALV=LV0-CPVMCL*(TCONV(J)-273.15)
  !     TNEW=TCONV(J)+ALV*(QCONV(J)-QSCONV(J))/(CPD*(1.-QCONV(J))+
  !    1      CL*QCONV(J)+QSCONV(J)*(CPV-CL+ALV*ALV/(RV*TCONV(J)*TCONV(J))))
  !     ALVNEW=LV0-CPVMCL*(TNEW-273.15)
  !     QNEW=(ALV*QCONV(J)-(TNEW-TCONV(J))*(CPD*(1.-QCONV(J))
  !    1     +CL*QCONV(J)))/ALVNEW
  !     PRECIP=PRECIP+24.*3600.*1.0E5*(PHCONV_HPA(J)-PHCONV_HPA(J+1))*
  !    1      (QCONV(J)-QNEW)/(G*DELT*ROWL)
  !     TCONV(J)=TNEW
  !     QCONV(J)=QNEW
  !     QSCONV(J)=QNEW
  !    END IF
  !  38  CONTINUE
  !END IF
  !
  !END IF
  !
  !  *** CALCULATE ARRAYS OF GEOPOTENTIAL, HEAT CAPACITY AND STATIC ENERGY
  !
    GZ(1)=0.0
    CPN(1)=CPD*(1.-QCONV(1))+QCONV(1)*CPV
    H(1)=TCONV(1)*CPN(1)
    LV(1)=LV0-CPVMCL*(TCONV(1)-273.15)
    HM(1)=LV(1)*QCONV(1)
    TV(1)=TCONV(1)*(1.+QCONV(1)*EPSI-QCONV(1))
    AHMIN=1.0E12
    IHMIN=NL
    DO I=2,NL+1
      TVX=TCONV(I)*(1.+QCONV(I)*EPSI-QCONV(I))
      TVY=TCONV(I-1)*(1.+QCONV(I-1)*EPSI-QCONV(I-1))
      GZ(I)=GZ(I-1)+0.5*RD*(TVX+TVY)*(PCONV_HPA(I-1)-PCONV_HPA(I))/ &
           PHCONV_HPA(I)
      CPN(I)=CPD*(1.-QCONV(I))+CPV*QCONV(I)
      H(I)=TCONV(I)*CPN(I)+GZ(I)
      LV(I)=LV0-CPVMCL*(TCONV(I)-273.15)
      HM(I)=(CPD*(1.-QCONV(I))+CL*QCONV(I))*(TCONV(I)-TCONV(1))+ &
           LV(I)*QCONV(I)+GZ(I)
      TV(I)=TCONV(I)*(1.+QCONV(I)*EPSI-QCONV(I))
  !
  !  ***  Find level of minimum moist static energy    ***
  !
      IF(I.GE.MINORIG.AND.HM(I).LT.AHMIN.AND.HM(I).LT.HM(I-1))THEN
       AHMIN=HM(I)
       IHMIN=I
      END IF
    END DO
    IHMIN=MIN(IHMIN, NL-1)
  !
  !  ***     Find that model level below the level of minimum moist       ***
  !  ***  static energy that has the maximum value of moist static energy ***
  !
    AHMAX=0.0
  !  ***  bug fixed: need to assign an initial value to NK
  !  HSO, 05.08.2009
    NK=MINORIG
    DO I=MINORIG,IHMIN
     IF(HM(I).GT.AHMAX)THEN
      NK=I
      AHMAX=HM(I)
     END IF
    END DO
  !
  !  ***  CHECK WHETHER PARCEL LEVEL TEMPERATURE AND SPECIFIC HUMIDITY   ***
  !  ***                          ARE REASONABLE                         ***
  !  ***      Skip convection if HM increases monotonically upward       ***
  !
    IF(TCONV(NK).LT.250.0.OR.QCONV(NK).LE.0.0.OR.IHMIN.EQ.(NL-1)) &
         THEN
     IFLAG=0
     CBMF=0.0
     RETURN
    END IF
  !
  !   ***  CALCULATE LIFTED CONDENSATION LEVEL OF AIR AT PARCEL ORIGIN LEVEL ***
  !   ***       (WITHIN 0.2% OF FORMULA OF BOLTON, MON. WEA. REV.,1980)      ***
  !
    RH=QCONV(NK)/QSCONV(NK)
    CHI=TCONV(NK)/(1669.0-122.0*RH-TCONV(NK))
    PLCL=PCONV_HPA(NK)*(RH**CHI)
    IF(PLCL.LT.200.0.OR.PLCL.GE.2000.0)THEN
     IFLAG=2
     CBMF=0.0
     RETURN
    END IF
  !
  !   ***  CALCULATE FIRST LEVEL ABOVE LCL (=ICB)  ***
  !
    ICB=NL-1
    DO I=NK+1,NL
     IF(PCONV_HPA(I).LT.PLCL)THEN
      ICB=MIN(ICB,I)
     END IF
    END DO
    IF(ICB.GE.(NL-1))THEN
     IFLAG=3
     CBMF=0.0
     RETURN
    END IF
  !
  !   *** FIND TEMPERATURE UP THROUGH ICB AND TEST FOR INSTABILITY           ***
  !
  !   *** SUBROUTINE TLIFT CALCULATES PART OF THE LIFTED PARCEL VIRTUAL      ***
  !   ***  TEMPERATURE, THE ACTUAL TEMPERATURE AND THE ADIABATIC             ***
  !   ***                   LIQUID WATER CONTENT                             ***
  !
    CALL TLIFT(GZ,ICB,NK,TVP,TP,CLW,ND,NL,1)
    DO I=NK,ICB
     TVP(I)=TVP(I)-TP(I)*QCONV(NK)
    END DO
  !
  !   ***  If there was no convection at last time step and parcel    ***
  !   ***       is stable at ICB then skip rest of calculation        ***
  !
    IF(CBMF.EQ.0.0.AND.TVP(ICB).LE.(TV(ICB)-DTMAX))THEN
     IFLAG=0
     RETURN
    END IF
  !
  !   ***  IF THIS POINT IS REACHED, MOIST CONVECTIVE ADJUSTMENT IS NECESSARY ***
  !
    IF(IFLAG.NE.4)IFLAG=1
  !
  !   ***  FIND THE REST OF THE LIFTED PARCEL TEMPERATURES          ***
  !
    CALL TLIFT(GZ,ICB,NK,TVP,TP,CLW,ND,NL,2)
  !
  !   ***  SET THE PRECIPITATION EFFICIENCIES AND THE FRACTION OF   ***
  !   ***          PRECIPITATION FALLING OUTSIDE OF CLOUD           ***
  !   ***      THESE MAY BE FUNCTIONS OF TP(I), PCONV_HPA(I) AND CLW(I)     ***
  !
    DO I=1,NK
     EP(I)=0.0
     SIGP(I)=SIGS
    END DO
    DO I=NK+1,NL
     TCA=TP(I)-273.15
     IF(TCA.GE.0.0)THEN
      ELACRIT=ELCRIT
     ELSE
      ELACRIT=ELCRIT*(1.0-TCA/TLCRIT)
     END IF
     ELACRIT=MAX(ELACRIT,0.0)
     EPMAX=0.999
     EP(I)=EPMAX*(1.0-ELACRIT/MAX(CLW(I),1.0E-8))
     EP(I)=MAX(EP(I),0.0)
     EP(I)=MIN(EP(I),EPMAX)
     SIGP(I)=SIGS
    END DO
  !
  !   ***       CALCULATE VIRTUAL TEMPERATURE AND LIFTED PARCEL     ***
  !   ***                    VIRTUAL TEMPERATURE                    ***
  !
    DO I=ICB+1,NL
     TVP(I)=TVP(I)-TP(I)*QCONV(NK)
    END DO
    TVP(NL+1)=TVP(NL)-(GZ(NL+1)-GZ(NL))/CPD
  !
  !   ***        NOW INITIALIZE VARIOUS ARRAYS USED IN THE COMPUTATIONS       ***
  !
    DO I=1,NL+1
     HP(I)=H(I)
     NENT(I)=0
     WATER(I)=0.0
     EVAP(I)=0.0
     WT(I)=OMTSNOW
     LVCP(I)=LV(I)/CPN(I)
     DO J=1,NL+1
      QENT(I,J)=QCONV(J)
      ELIJ(I,J)=0.0
      SIJ(I,J)=0.0
     END DO
    END DO
    QP(1)=QCONV(1)
    DO I=2,NL+1
     QP(I)=QCONV(I-1)
    END DO
  !
  !  ***  FIND THE FIRST MODEL LEVEL (INB1) ABOVE THE PARCEL'S      ***
  !  ***          HIGHEST LEVEL OF NEUTRAL BUOYANCY                 ***
  !  ***     AND THE HIGHEST LEVEL OF POSITIVE CAPE (INB)           ***
  !
    CAPE=0.0
    CAPEM=0.0
    INB=ICB+1
    INB1=INB
    BYP=0.0
    DO I=ICB+1,NL-1
     BY=(TVP(I)-TV(I))*(PHCONV_HPA(I)-PHCONV_HPA(I+1))/PCONV_HPA(I)
     CAPE=CAPE+BY
     IF(BY.GE.0.0)INB1=I+1
     IF(CAPE.GT.0.0)THEN
      INB=I+1
      BYP=(TVP(I+1)-TV(I+1))*(PHCONV_HPA(I+1)-PHCONV_HPA(I+2))/ &
           PCONV_HPA(I+1)
      CAPEM=CAPE
     END IF
    END DO
    INB=MAX(INB,INB1)
    CAPE=CAPEM+BYP
    DEFRAC=CAPEM-CAPE
    DEFRAC=MAX(DEFRAC,0.001)
    FRAC=-CAPE/DEFRAC
    FRAC=MIN(FRAC,1.0)
    FRAC=MAX(FRAC,0.0)
  !
  !   ***   CALCULATE LIQUID WATER STATIC ENERGY OF LIFTED PARCEL   ***
  !
    DO I=ICB,INB
     HP(I)=H(NK)+(LV(I)+(CPD-CPV)*TCONV(I))*EP(I)*CLW(I)
    END DO
  !
  !   ***  CALCULATE CLOUD BASE MASS FLUX AND RATES OF MIXING, M(I),  ***
  !   ***                   AT EACH MODEL LEVEL                       ***
  !
    DBOSUM=0.0
  !
  !   ***     INTERPOLATE DIFFERENCE BETWEEN LIFTED PARCEL AND      ***
  !   ***  ENVIRONMENTAL TEMPERATURES TO LIFTED CONDENSATION LEVEL  ***
  !
    TVPPLCL=TVP(ICB-1)-RD*TVP(ICB-1)*(PCONV_HPA(ICB-1)-PLCL)/ &
         (CPN(ICB-1)*PCONV_HPA(ICB-1))
    TVAPLCL=TV(ICB)+(TVP(ICB)-TVP(ICB+1))*(PLCL-PCONV_HPA(ICB))/ &
         (PCONV_HPA(ICB)-PCONV_HPA(ICB+1))
    DTPBL=0.0
    DO I=NK,ICB-1
     DTPBL=DTPBL+(TVP(I)-TV(I))*(PHCONV_HPA(I)-PHCONV_HPA(I+1))
    END DO
    DTPBL=DTPBL/(PHCONV_HPA(NK)-PHCONV_HPA(ICB))
    DTMIN=TVPPLCL-TVAPLCL+DTMAX+DTPBL
    DTMA=DTMIN
  !
  !   ***  ADJUST CLOUD BASE MASS FLUX   ***
  !
  CBMFOLD=CBMF
  ! *** C. Forster: adjustment of CBMF is not allowed to depend on FLEXPART timestep
  DELT0=DELT/3.
  DAMPS=DAMP*DELT/DELT0
  CBMF=(1.-DAMPS)*CBMF+0.1*ALPHA*DTMA
  CBMF=MAX(CBMF,0.0)
  !
  !   *** If cloud base mass flux is zero, skip rest of calculation  ***
  !
  IF(CBMF.EQ.0.0.AND.CBMFOLD.EQ.0.0)THEN
   RETURN
  END IF

  !
  !   ***   CALCULATE RATES OF MIXING,  M(I)   ***
  !
  M(ICB)=0.0
  DO I=ICB+1,INB
   K=MIN(I,INB1)
   DBO=ABS(TV(K)-TVP(K))+ &
        ENTP*0.02*(PHCONV_HPA(K)-PHCONV_HPA(K+1))
   DBOSUM=DBOSUM+DBO
   M(I)=CBMF*DBO
  END DO
  DO I=ICB+1,INB
   M(I)=M(I)/DBOSUM
  END DO
  !
  !   ***  CALCULATE ENTRAINED AIR MASS FLUX (MENT), TOTAL WATER MIXING  ***
  !   ***     RATIO (QENT), TOTAL CONDENSED WATER (ELIJ), AND MIXING     ***
  !   ***                        FRACTION (SIJ)                          ***
  !
    DO I=ICB+1,INB
     QTI=QCONV(NK)-EP(I)*CLW(I)
     DO J=ICB,INB
      BF2=1.+LV(J)*LV(J)*QSCONV(J)/(RV*TCONV(J)*TCONV(J)*CPD)
      ANUM=H(J)-HP(I)+(CPV-CPD)*TCONV(J)*(QTI-QCONV(J))
      DENOM=H(I)-HP(I)+(CPD-CPV)*(QCONV(I)-QTI)*TCONV(J)
      DEI=DENOM
      IF(ABS(DEI).LT.0.01)DEI=0.01
      SIJ(I,J)=ANUM/DEI
      SIJ(I,I)=1.0
      ALTEM=SIJ(I,J)*QCONV(I)+(1.-SIJ(I,J))*QTI-QSCONV(J)
      ALTEM=ALTEM/BF2
      CWAT=CLW(J)*(1.-EP(J))
      STEMP=SIJ(I,J)
      IF((STEMP.LT.0.0.OR.STEMP.GT.1.0.OR. &
           ALTEM.GT.CWAT).AND.J.GT.I)THEN
       ANUM=ANUM-LV(J)*(QTI-QSCONV(J)-CWAT*BF2)
       DENOM=DENOM+LV(J)*(QCONV(I)-QTI)
       IF(ABS(DENOM).LT.0.01)DENOM=0.01
       SIJ(I,J)=ANUM/DENOM
       ALTEM=SIJ(I,J)*QCONV(I)+(1.-SIJ(I,J))*QTI-QSCONV(J)
       ALTEM=ALTEM-(BF2-1.)*CWAT
      END IF
      IF(SIJ(I,J).GT.0.0.AND.SIJ(I,J).LT.0.9)THEN
       QENT(I,J)=SIJ(I,J)*QCONV(I)+(1.-SIJ(I,J))*QTI
       ELIJ(I,J)=ALTEM
       ELIJ(I,J)=MAX(0.0,ELIJ(I,J))
       MENT(I,J)=M(I)/(1.-SIJ(I,J))
       NENT(I)=NENT(I)+1
      END IF
      SIJ(I,J)=MAX(0.0,SIJ(I,J))
      SIJ(I,J)=MIN(1.0,SIJ(I,J))
     END DO
  !
  !   ***   IF NO AIR CAN ENTRAIN AT LEVEL I ASSUME THAT UPDRAFT DETRAINS  ***
  !   ***   AT THAT LEVEL AND CALCULATE DETRAINED AIR FLUX AND PROPERTIES  ***
  !
     IF(NENT(I).EQ.0)THEN
      MENT(I,I)=M(I)
      QENT(I,I)=QCONV(NK)-EP(I)*CLW(I)
      ELIJ(I,I)=CLW(I)
      SIJ(I,I)=1.0
     END IF
    END DO
    SIJ(INB,INB)=1.0
  !
  !   ***  NORMALIZE ENTRAINED AIR MASS FLUXES TO REPRESENT EQUAL  ***
  !   ***              PROBABILITIES OF MIXING                     ***
  !
    DO I=ICB+1,INB
    IF(NENT(I).NE.0)THEN
     QP1=QCONV(NK)-EP(I)*CLW(I)
     ANUM=H(I)-HP(I)-LV(I)*(QP1-QSCONV(I))
     DENOM=H(I)-HP(I)+LV(I)*(QCONV(I)-QP1)
     IF(ABS(DENOM).LT.0.01)DENOM=0.01
     SCRIT=ANUM/DENOM
     ALT=QP1-QSCONV(I)+SCRIT*(QCONV(I)-QP1)
     IF(ALT.LT.0.0)SCRIT=1.0
     SCRIT=MAX(SCRIT,0.0)
     ASIJ=0.0
     SMIN=1.0
     DO J=ICB,INB
      IF(SIJ(I,J).GT.0.0.AND.SIJ(I,J).LT.0.9)THEN
       IF(J.GT.I)THEN
        SMID=MIN(SIJ(I,J),SCRIT)
        SJMAX=SMID
        SJMIN=SMID
        IF(SMID.LT.SMIN.AND.SIJ(I,J+1).LT.SMID)THEN
         SMIN=SMID
         SJMAX=MIN(SIJ(I,J+1),SIJ(I,J),SCRIT)
         SJMIN=MAX(SIJ(I,J-1),SIJ(I,J))
         SJMIN=MIN(SJMIN,SCRIT)
        END IF
       ELSE
        SJMAX=MAX(SIJ(I,J+1),SCRIT)
        SMID=MAX(SIJ(I,J),SCRIT)
        SJMIN=0.0
        IF(J.GT.1)SJMIN=SIJ(I,J-1)
        SJMIN=MAX(SJMIN,SCRIT)
       END IF
       DELP=ABS(SJMAX-SMID)
       DELM=ABS(SJMIN-SMID)
       ASIJ=ASIJ+(DELP+DELM)*(PHCONV_HPA(J)-PHCONV_HPA(J+1))
       MENT(I,J)=MENT(I,J)*(DELP+DELM)* &
            (PHCONV_HPA(J)-PHCONV_HPA(J+1))
      END IF
     END DO
     ASIJ=MAX(1.0E-21,ASIJ)
     ASIJ=1.0/ASIJ
     DO J=ICB,INB
      MENT(I,J)=MENT(I,J)*ASIJ
     END DO
     BSUM=0.0
     DO J=ICB,INB
      BSUM=BSUM+MENT(I,J)
     END DO
     IF(BSUM.LT.1.0E-18)THEN
      NENT(I)=0
      MENT(I,I)=M(I)
      QENT(I,I)=QCONV(NK)-EP(I)*CLW(I)
      ELIJ(I,I)=CLW(I)
      SIJ(I,I)=1.0
     END IF
    END IF
    END DO
  !
  !   ***  CHECK WHETHER EP(INB)=0, IF SO, SKIP PRECIPITATING    ***
  !   ***             DOWNDRAFT CALCULATION                      ***
  !
    IF(EP(INB).LT.0.0001)GOTO 405
  !
  !   ***  INTEGRATE LIQUID WATER EQUATION TO FIND CONDENSED WATER   ***
  !   ***                AND CONDENSED WATER FLUX                    ***
  !
    JTT=2
  !
  !    ***                    BEGIN DOWNDRAFT LOOP                    ***
  !
    DO I=INB,1,-1
  !
  !    ***              CALCULATE DETRAINED PRECIPITATION             ***
  !
    WDTRAIN=G*EP(I)*M(I)*CLW(I)
    IF(I.GT.1)THEN
     DO J=1,I-1
     AWAT=ELIJ(J,I)-(1.-EP(I))*CLW(I)
     AWAT=MAX(0.0,AWAT)
       WDTRAIN=WDTRAIN+G*AWAT*MENT(J,I)
     END DO
    END IF
  !
  !    ***    FIND RAIN WATER AND EVAPORATION USING PROVISIONAL   ***
  !    ***              ESTIMATES OF QP(I)AND QP(I-1)             ***
  !
  !
  !  ***  Value of terminal velocity and coefficient of evaporation for snow   ***
  !
    COEFF=COEFFS
    WT(I)=OMTSNOW
  !
  !  ***  Value of terminal velocity and coefficient of evaporation for rain   ***
  !
    IF(TCONV(I).GT.273.0)THEN
     COEFF=COEFFR
     WT(I)=OMTRAIN
    END IF
    QSM=0.5*(QCONV(I)+QP(I+1))
    AFAC=COEFF*PHCONV_HPA(I)*(QSCONV(I)-QSM)/ &
         (1.0E4+2.0E3*PHCONV_HPA(I)*QSCONV(I))
    AFAC=MAX(AFAC,0.0)
    SIGT=SIGP(I)
    SIGT=MAX(0.0,SIGT)
    SIGT=MIN(1.0,SIGT)
    B6=100.*(PHCONV_HPA(I)-PHCONV_HPA(I+1))*SIGT*AFAC/WT(I)
    C6=(WATER(I+1)*WT(I+1)+WDTRAIN/SIGD)/WT(I)
    REVAP=0.5*(-B6+SQRT(B6*B6+4.*C6))
    EVAP(I)=SIGT*AFAC*REVAP
    WATER(I)=REVAP*REVAP
  !
  !    ***  CALCULATE PRECIPITATING DOWNDRAFT MASS FLUX UNDER     ***
  !    ***              HYDROSTATIC APPROXIMATION                 ***
  !
    IF(I.EQ.1)GOTO 360
    DHDP=(H(I)-H(I-1))/(PCONV_HPA(I-1)-PCONV_HPA(I))
    DHDP=MAX(DHDP,10.0)
    MP(I)=100.*GINV*LV(I)*SIGD*EVAP(I)/DHDP
    MP(I)=MAX(MP(I),0.0)
  !
  !   ***   ADD SMALL AMOUNT OF INERTIA TO DOWNDRAFT              ***
  !
    FAC=20.0/(PHCONV_HPA(I-1)-PHCONV_HPA(I))
    MP(I)=(FAC*MP(I+1)+MP(I))/(1.+FAC)
  !
  !    ***      FORCE MP TO DECREASE LINEARLY TO ZERO                 ***
  !    ***      BETWEEN ABOUT 950 MB AND THE SURFACE                  ***
  !
      IF(PCONV_HPA(I).GT.(0.949*PCONV_HPA(1)))THEN
       JTT=MAX(JTT,I)
       MP(I)=MP(JTT)*(PCONV_HPA(1)-PCONV_HPA(I))/(PCONV_HPA(1)- &
            PCONV_HPA(JTT))
      END IF
  360   CONTINUE
  !
  !    ***       FIND MIXING RATIO OF PRECIPITATING DOWNDRAFT     ***
  !
    IF(I.EQ.INB)GOTO 400
    IF(I.EQ.1)THEN
     QSTM=QSCONV(1)
    ELSE
     QSTM=QSCONV(I-1)
    END IF
    IF(MP(I).GT.MP(I+1))THEN
      RAT=MP(I+1)/MP(I)
      QP(I)=QP(I+1)*RAT+QCONV(I)*(1.0-RAT)+100.*GINV* &
           SIGD*(PHCONV_HPA(I)-PHCONV_HPA(I+1))*(EVAP(I)/MP(I))
     ELSE
      IF(MP(I+1).GT.0.0)THEN
        QP(I)=(GZ(I+1)-GZ(I)+QP(I+1)*(LV(I+1)+TCONV(I+1)*(CL-CPD))+ &
             CPD*(TCONV(I+1)-TCONV(I)))/(LV(I)+TCONV(I)*(CL-CPD))
      END IF
    END IF
    QP(I)=MIN(QP(I),QSTM)
    QP(I)=MAX(QP(I),0.0)
400 CONTINUE
    END DO
  !
  !   ***  CALCULATE SURFACE PRECIPITATION IN MM/DAY     ***
  !
    PRECIP=PRECIP+WT(1)*SIGD*WATER(1)*3600.*24000./(ROWL*G)
  !
  405   CONTINUE
  !
  !   ***  CALCULATE DOWNDRAFT VELOCITY SCALE AND SURFACE TEMPERATURE AND  ***
  !   ***                    WATER VAPOR FLUCTUATIONS                      ***
  !
  WD=BETA*ABS(MP(ICB))*0.01*RD*TCONV(ICB)/(SIGD*PCONV_HPA(ICB))
  QPRIME=0.5*(QP(1)-QCONV(1))
  TPRIME=LV0*QPRIME/CPD
  !
  !   ***  CALCULATE TENDENCIES OF LOWEST LEVEL POTENTIAL TEMPERATURE  ***
  !   ***                      AND MIXING RATIO                        ***
  !
    DPINV=0.01/(PHCONV_HPA(1)-PHCONV_HPA(2))
    AM=0.0
    IF(NK.EQ.1)THEN
     DO K=2,INB
       AM=AM+M(K)
     END DO
    END IF
  ! save saturated upward mass flux for first level
    FUP(1)=AM
    IF((2.*G*DPINV*AM).GE.DELTI)IFLAG=4
    FT(1)=FT(1)+G*DPINV*AM*(TCONV(2)-TCONV(1)+(GZ(2)-GZ(1))/CPN(1))
    FT(1)=FT(1)-LVCP(1)*SIGD*EVAP(1)
    FT(1)=FT(1)+SIGD*WT(2)*(CL-CPD)*WATER(2)*(TCONV(2)- &
         TCONV(1))*DPINV/CPN(1)
    FQ(1)=FQ(1)+G*MP(2)*(QP(2)-QCONV(1))* &
         DPINV+SIGD*EVAP(1)
    FQ(1)=FQ(1)+G*AM*(QCONV(2)-QCONV(1))*DPINV
    DO J=2,INB
     FQ(1)=FQ(1)+G*DPINV*MENT(J,1)*(QENT(J,1)-QCONV(1))
    END DO
  !
  !   ***  CALCULATE TENDENCIES OF POTENTIAL TEMPERATURE AND MIXING RATIO  ***
  !   ***               AT LEVELS ABOVE THE LOWEST LEVEL                   ***
  !
  !   ***  FIRST FIND THE NET SATURATED UPDRAFT AND DOWNDRAFT MASS FLUXES  ***
  !   ***                      THROUGH EACH LEVEL                          ***
  !
    DO I=2,INB
    DPINV=0.01/(PHCONV_HPA(I)-PHCONV_HPA(I+1))
    CPINV=1.0/CPN(I)
    AMP1=0.0
    AD=0.0
    IF(I.GE.NK)THEN
     DO K=I+1,INB+1
       AMP1=AMP1+M(K)
     END DO
    END IF
    DO K=1,I
    DO J=I+1,INB+1
     AMP1=AMP1+MENT(K,J)
    END DO
    END DO
  ! save saturated upward mass flux
    FUP(I)=AMP1
    IF((2.*G*DPINV*AMP1).GE.DELTI)IFLAG=4
    DO K=1,I-1
    DO J=I,INB
     AD=AD+MENT(J,K)
    END DO
    END DO
  ! save saturated downward mass flux
    FDOWN(I)=AD
    FT(I)=FT(I)+G*DPINV*(AMP1*(TCONV(I+1)-TCONV(I)+(GZ(I+1)-GZ(I))* &
         CPINV)-AD*(TCONV(I)-TCONV(I-1)+(GZ(I)-GZ(I-1))*CPINV)) &
         -SIGD*LVCP(I)*EVAP(I)
    FT(I)=FT(I)+G*DPINV*MENT(I,I)*(HP(I)-H(I)+ &
         TCONV(I)*(CPV-CPD)*(QCONV(I)-QENT(I,I)))*CPINV
    FT(I)=FT(I)+SIGD*WT(I+1)*(CL-CPD)*WATER(I+1)* &
         (TCONV(I+1)-TCONV(I))*DPINV*CPINV
    FQ(I)=FQ(I)+G*DPINV*(AMP1*(QCONV(I+1)-QCONV(I))- &
         AD*(QCONV(I)-QCONV(I-1)))
    DO K=1,I-1
     AWAT=ELIJ(K,I)-(1.-EP(I))*CLW(I)
     AWAT=MAX(AWAT,0.0)
     FQ(I)=FQ(I)+G*DPINV*MENT(K,I)*(QENT(K,I)-AWAT-QCONV(I))
    END DO
    DO K=I,INB
     FQ(I)=FQ(I)+G*DPINV*MENT(K,I)*(QENT(K,I)-QCONV(I))
    END DO
    FQ(I)=FQ(I)+SIGD*EVAP(I)+G*(MP(I+1)* &
         (QP(I+1)-QCONV(I))-MP(I)*(QP(I)-QCONV(I-1)))*DPINV
    END DO
  !
  !   *** Adjust tendencies at top of convection layer to reflect  ***
  !   ***       actual position of the level zero CAPE             ***
  !
    FQOLD=FQ(INB)
    FQ(INB)=FQ(INB)*(1.-FRAC)
    FQ(INB-1)=FQ(INB-1)+FRAC*FQOLD*((PHCONV_HPA(INB)- &
         PHCONV_HPA(INB+1))/ &
         (PHCONV_HPA(INB-1)-PHCONV_HPA(INB)))*LV(INB)/LV(INB-1)
    FTOLD=FT(INB)
    FT(INB)=FT(INB)*(1.-FRAC)
    FT(INB-1)=FT(INB-1)+FRAC*FTOLD*((PHCONV_HPA(INB)- &
         PHCONV_HPA(INB+1))/ &
         (PHCONV_HPA(INB-1)-PHCONV_HPA(INB)))*CPN(INB)/CPN(INB-1)
  !
  !   ***   Very slightly adjust tendencies to force exact   ***
  !   ***     enthalpy, momentum and tracer conservation     ***
  !
    ENTS=0.0
    DO I=1,INB
     ENTS=ENTS+(CPN(I)*FT(I)+LV(I)*FQ(I))* &
          (PHCONV_HPA(I)-PHCONV_HPA(I+1))
    END DO
    ENTS=ENTS/(PHCONV_HPA(1)-PHCONV_HPA(INB+1))
    DO I=1,INB
     FT(I)=FT(I)-ENTS/CPN(I)
    END DO

  ! ************************************************
  ! **** DETERMINE MASS DISPLACEMENT MATRIX
  ! ***** AND COMPENSATING SUBSIDENCE
  ! ************************************************

  ! mass displacement matrix due to saturated up-and downdrafts
  ! inside the cloud and determine compensating subsidence
  ! FUP(I) (saturated updrafts), FDOWN(I) (saturated downdrafts) are assumed to be
  ! balanced by  compensating subsidence (SUB(I))
  ! FDOWN(I) and SUB(I) defined positive downwards

  ! NCONVTOP IS THE TOP LEVEL AT WHICH CONVECTIVE MASS FLUXES ARE DIAGNOSED
  ! EPSILON IS A SMALL NUMBER

   SUB(1)=0.
   NCONVTOP=1
   do i=1,INB+1
   do j=1,INB+1
    if (j.eq.NK) then
     FMASS(j,i)=FMASS(j,i)+M(i)
    endif
     FMASS(j,i)=FMASS(j,i)+MENT(j,i)
     IF (FMASS(J,I).GT.EPSILON) NCONVTOP=MAX(NCONVTOP,I,J)
   end do
   if (i.gt.1) then
    SUB(i)=FUP(i-1)-FDOWN(i)
   endif
   end do
   NCONVTOP=NCONVTOP+1

    RETURN
  !
END SUBROUTINE CONVECT
!
! ---------------------------------------------------------------------------
!
SUBROUTINE TLIFT(GZ,ICB,NK,TVP,TPK,CLW,ND,NL,KK)
  !
  !-cv
  use par_mod
  use conv_mod

  implicit none
  !-cv
  !====>Begin Module TLIFT      File convect.f      Undeclared variables
  !
  !Argument variables
  !
  integer :: icb, kk, nd, nk, nl
  !
  !Local variables
  !
  integer :: i, j, nsb, nst
  !
  real :: ah0, ahg, alv, cpinv, cpp, denom
  real :: es, qg, rg, s, tc, tg
  !
  !====>End Module   TLIFT      File convect.f

    REAL :: GZ(ND),TPK(ND),CLW(ND)
    REAL :: TVP(ND)
  !
  !   ***   ASSIGN VALUES OF THERMODYNAMIC CONSTANTS     ***
  !
    REAL,PARAMETER :: CPD=1005.7
    REAL,PARAMETER :: CPV=1870.0
    REAL,PARAMETER :: CL=2500.0
    REAL,PARAMETER :: RV=461.5
    REAL,PARAMETER :: RD=287.04
    REAL,PARAMETER :: LV0=2.501E6
  !
    REAL,PARAMETER :: CPVMCL=CL-CPV
    REAL,PARAMETER :: EPS0=RD/RV
    REAL,PARAMETER :: EPSI=1./EPS0
  !
  !   ***  CALCULATE CERTAIN PARCEL QUANTITIES, INCLUDING STATIC ENERGY   ***
  !
    AH0=(CPD*(1.-QCONV(NK))+CL*QCONV(NK))*TCONV(NK)+QCONV(NK)* &
         (LV0-CPVMCL*( &
         TCONV(NK)-273.15))+GZ(NK)
    CPP=CPD*(1.-QCONV(NK))+QCONV(NK)*CPV
    CPINV=1./CPP
  !
    IF(KK.EQ.1)THEN
  !
  !   ***   CALCULATE LIFTED PARCEL QUANTITIES BELOW CLOUD BASE   ***
  !
    DO I=1,ICB-1
     CLW(I)=0.0
    END DO
    DO I=NK,ICB-1
     TPK(I)=TCONV(NK)-(GZ(I)-GZ(NK))*CPINV
     TVP(I)=TPK(I)*(1.+QCONV(NK)*EPSI)
    END DO
    END IF
  !
  !    ***  FIND LIFTED PARCEL QUANTITIES ABOVE CLOUD BASE    ***
  !
    NST=ICB
    NSB=ICB
    IF(KK.EQ.2)THEN
     NST=NL
     NSB=ICB+1
    END IF
    DO I=NSB,NST
     TG=TCONV(I)
     QG=QSCONV(I)
     ALV=LV0-CPVMCL*(TCONV(I)-273.15)
     DO J=1,2
      S=CPD+ALV*ALV*QG/(RV*TCONV(I)*TCONV(I))
      S=1./S
      AHG=CPD*TG+(CL-CPD)*QCONV(NK)*TCONV(I)+ALV*QG+GZ(I)
      TG=TG+S*(AH0-AHG)
      TG=MAX(TG,35.0)
      TC=TG-273.15
      DENOM=243.5+TC
      IF(TC.GE.0.0)THEN
       ES=6.112*EXP(17.67*TC/DENOM)
      ELSE
       ES=EXP(23.33086-6111.72784/TG+0.15215*LOG(TG))
      END IF
      QG=EPS0*ES/(PCONV_HPA(I)-ES*(1.-EPS0))
     END DO
     ALV=LV0-CPVMCL*(TCONV(I)-273.15)
     TPK(I)=(AH0-(CL-CPD)*QCONV(NK)*TCONV(I)-GZ(I)-ALV*QG)/CPD
     CLW(I)=QCONV(NK)-QG
     CLW(I)=MAX(0.0,CLW(I))
     RG=QG/(1.-QCONV(NK))
     TVP(I)=TPK(I)*(1.+RG*EPSI)
    END DO
    RETURN
END SUBROUTINE TLIFT
