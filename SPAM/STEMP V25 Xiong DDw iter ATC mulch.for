C=======================================================================
C  COPYRIGHT 1998-2010 The University of Georgia, Griffin, Georgia
C                      University of Florida, Gainesville, Florida
C                      Iowa State University, Ames, Iowa
C                      International Center for Soil Fertility and
C                       Agricultural Development, Muscle Shoals, Alabama
C                      University of Guelph, Guelph, Ontario
C  ALL RIGHTS RESERVED
C=======================================================================
C=======================================================================
C  STEMP, Subroutine
C
C  Determines soil temperature by layer
C-----------------------------------------------------------------------
C  Revision history
C  12/01/1980     Originally based on EPIC soil temperature routines
C  12/01/1999 CHP Combined SOILT and INSOILT into STEMP.for for modular
C                 format.
C  01/01/2000 AJG Added surface temperature for the CENTURY-based
C                 SOM/soil-N module.
C  07/01/2000 GH  Incorporated in CROPGRO
C  06/07/2002 GH  Modified for crop rotations
C  06/07/2002 GH  Moved TAMP and TAV to IPWTH
C  09/17/2002 CHP Added computation for ISWWAT = 'N' (necessary for potato)
C  07/15/2003 CHP No re-initialization for sequenced runs.
C  01/14/2005 CHP Added METMP = 3: Corrected water content in temp. eqn.
!  07/24/2006 CHP Use MSALB instead of SALB (includes mulch and soil
!                 water effects on albedo)
!  12/09/2008 CHP Remove METMP and code for old (incorrect) soil water effect
C-----------------------------------------------------------------------
C  Called : Main
C  Calls  : SOILT
C=======================================================================

      SUBROUTINE STEMP(CONTROL, ISWITCH,
     &    SOILPROP, SRAD, SW, TAVG, TMAX, XLAT, TAV, TAMP,!Input
     &    EOP, TRWUP, XHLAI, VPD,TDEW, ES,EP,WINDSP,CANHT,!Input
     &    MULCH,                                          !Input
     &    SRFTEMP, ST)                                    !Output

C-----------------------------------------------------------------------
      USE ModuleDefs     !Definitions of constructed variable types,
                         ! which contain control information, soil
                         ! parameters, hourly weather data.
      IMPLICIT  NONE
      EXTERNAL YR_DOY,ERROR,FIND,SOILTBK,OPSTEMPBK,MULCHWATER
      SAVE

      CHARACTER*1  RNMODE, ISWWAT,MEEVP !, IDETL (CSVC ADD MEEVP)
      CHARACTER*6  SECTION
      CHARACTER*6, PARAMETER :: ERRKEY = "STEMP "
      CHARACTER*30 FILEIO

      INTEGER DOY, DYNAMIC, I, L, NLAYR,J,M
      INTEGER RUN, YRDOY, YEAR
      INTEGER ERRNUM, FOUND, LNUM, LUNIO

      REAL ABD, ALBEDO, ATOT, B, CUMDPT
      REAL DP, FX, HDAY, ICWD, PESW,MSALB,SRAD,SRFTEMP
      REAL TAMP, TAV, TAVG, TBD, TMAX, XLAT, WW
      REAL TDL, TLL, TSW, TA, DT
      REAL TMA(5),DelA(5),DelTOT
      REAL EOP,TRWUP,XHLAI,VPD, TDEW, SWFAC,EP1,AVP,SVP,ES,EP,WINDSP
      REAL, DIMENSION(NL) :: BD, DLAYR, DS, DUL, LL, ST, SW, SWI, DSMID,
     &                     CLAY,SILT,SAND,OC
      REAL CLAYV(NL),SILTV(NL),SANDV(NL),OMV(NL),TotSolid(NL),POR(NL)
      REAL CLAYC(NL),SILTC(NL),SANDC(NL),OMC(NL),OM(NL),Totmineral(NL)
      REAL ClayFrac(NL), SiltFrac(NL), SandFrac(NL), OMFrac(NL)
      REAL TcondDry(NL), TcondS(NL), TcondSat(NL), SWREL(NL)
      REAL  PX, QX, RX, SX
      REAL HeatCap(NL),STBot(NL),AMP(NL),DSMIDV2(NL)
      REAL STCOND(NL)
      REAL Omega, Del,STa(NL),STboti(NL),AMPi(NL)
      REAL ASTCOND,AHeatCap, DampDa,DampDw
      REAL SolAvg,WINDmps,ESWatt,EPWatt,X1S,X1P
      REAL X2S,X2P,XEPS,Xsky,Plht,AEROra,RiNo,MEK,PHI,XG
      REAL AEROraSM1,AEROraS0,AEROraPM1,AEROraP0
      REAL XcubeS,XcubeP,DelS,DelP
      REAL Ga,GwSM1,GwS0,GwPM1,GwP0,TSM1,TS0,TPM1,TP0
      REAL RADSM1,RADS0,RADPM1,RADP0,HSM1,HS0,HPM1,HP0
      REAL FSM1,FS0,FPM1,FP0,FLAGS,FLAGP,CANHT,STEP
      REAL MULCHMASS,MULCHCOVER,MULCHTHICK,MULCHWAT,MULCHEVAP,MULCHALB
      REAL Sfrac,Pfrac,Mfrac,TM0,TMBot

!-----------------------------------------------------------------------
      TYPE (ControlType) CONTROL
      TYPE (SoilType)    SOILPROP
      TYPE (SwitchType)  ISWITCH
      TYPE (MulchType)   MULCH
!      TYPE (WeatherType) WEATHER

!     Check for output verbosity
!     IDETL  = ISWITCH % IDETL
!     IDETL = 'N', '0' (zero), suppress output
!     IF (INDEX('N0',IDETL) > 0) RETURN

!     Transfer values from constructed data types into local variables.
      DYNAMIC = CONTROL % DYNAMIC
      YRDOY   = CONTROL % YRDOY

      ISWWAT = ISWITCH % ISWWAT
      MEEVP  = ISWITCH % MEEVP  !SVC

      BD     = SOILPROP % BD     ! bulk den g/cm3
      DLAYR  = SOILPROP % DLAYR  ! thickness of soil layer cm
      DS     = SOILPROP % DS     ! cumulative depth of soil layer cm
      DUL    = SOILPROP % DUL    ! drained upper limit (cm3/cm3)
      LL     = SOILPROP % LL     ! lower limit (cm3/cm3)
      NLAYR  = SOILPROP % NLAYR  ! number of soil layers
      MSALB  = SOILPROP % MSALB  ! soil soil albedo with mulch
                                 !  and soil water effects
      CLAY   = SOILPROP % CLAY   ! clay (% by weight)
      SILT   = SOILPROP % SILT   ! silt (% by weight)
      SAND   = SOILPROP % SAND   ! sand (% by weight)
      OC     = SOILPROP % OC     ! organic carbon (g C/g soil)
      
      MULCHMASS  = MULCH % MULCHMASS
      MULCHCOVER = MULCH % MULCHCOVER
      MULCHTHICK = MULCH % MULCHTHICK
      MULCHWAT   = MULCH % MULCHWAT
      MULCHEVAP  = MULCH % MULCHEVAP
      MULCHALB   = MULCH % MULCHALB
      MULCHTHICK = MULCH % MULCHTHICK
      
!      SRAD   = WEATHER % SRAD    ! Solar radiation (MJ/(m2 day)
!      WINDSP = WEATHER % WINDSP  ! wind speed (km/day)
!-----------------------------------------------------------------------
      CALL YR_DOY(YRDOY, YEAR, DOY)
      
!-------------------------------------------------------------
!      Compute Water Stress Factor       
! ------------------------------------------------------------
          SWFAC  = 1.0
          IF(ISWWAT.NE.'N') THEN
             IF (EOP .GT. 0.0) THEN
                EP1 = EOP * 0.1 ! mm/day to cm/day
                IF (EP1 .GE. TRWUP) THEN
                  SWFAC = TRWUP / EP1
                ENDIF
             ENDIF
          ENDIF
!      Calculate vapor pressure deficit of air
!         Following Buck (1981, J. Applied Met. 20:1527-1532. (kPa)
      IF(TAVG .GE. 0.) THEN
            AVP = 0.61121*EXP(17.368*TDEW/(238.88 + TDEW))  ! for water
          ELSE
            AVP = 0.61115*EXP(22.452*TDEW/(272.55 + TDEW))  ! for ice
      ENDIF
!     Saturation vapor pressure (kPa)
           IF(TAVG .GE. 0.) THEN
            SVP = 0.61121*EXP(17.368*TAVG/(238.88 + TAVG))  ! for water
          ELSE
            SVP = 0.6115*EXP(22.452*TAVG/(272.55 + TAVG))  ! for ice
      ENDIF
          VPD = SVP - AVP   ! kPa
!
!***********************************************************************
!***********************************************************************
!     Run initialization - run once per simulation
!***********************************************************************
!      IF (DYNAMIC .EQ. RUNINIT) THEN
!-----------------------------------------------------------------------
!***********************************************************************
!***********************************************************************
!     Seasonal initialization - run once per season
!***********************************************************************
!      ELSEIF (DYNAMIC .EQ. SEASINIT) THEN
      IF (DYNAMIC .EQ. SEASINIT) THEN
!-----------------------------------------------------------------------
      FILEIO  = CONTROL % FILEIO
      LUNIO   = CONTROL % LUNIO
      RUN     = CONTROL % RUN
      RNMODE  = CONTROL % RNMODE

      IF (RUN .EQ. 1 .OR. INDEX('QF',RNMODE) .LE. 0) THEN

        IF (ISWWAT .NE. 'N') THEN
!         Read inital soil water values from FILEIO
!         (not yet done in WATBAL, so need to do here)
          OPEN (LUNIO, FILE = FILEIO, STATUS = 'OLD', IOSTAT=ERRNUM)
          IF (ERRNUM .NE. 0) CALL ERROR(ERRKEY,ERRNUM,FILEIO,0)
          SECTION = '*INITI'
          CALL FIND(LUNIO, SECTION, LNUM, FOUND)
          IF (FOUND .EQ. 0) CALL ERROR(SECTION, 42, FILEIO, LNUM)

!         Initial depth to water table (not currently used)
          READ(LUNIO,'(40X,F6.0)',IOSTAT=ERRNUM) ICWD ; LNUM = LNUM + 1
          IF (ERRNUM .NE. 0) CALL ERROR(ERRKEY,ERRNUM,FILEIO,LNUM)

          DO L = 1, NLAYR
            READ(LUNIO,'(9X,F5.3)',IOSTAT=ERRNUM) SWI(L)
            LNUM = LNUM + 1
            IF (ERRNUM .NE. 0) CALL ERROR(ERRKEY,ERRNUM,FILEIO,LNUM)
            IF (SWI(L) .LT. LL(L)) SWI(L) = LL(L)
          ENDDO

          CLOSE (LUNIO)
        ELSE
          SWI = DUL
        ENDIF

        IF (XLAT .LT. 0.0) THEN
          HDAY =  20.0           !DOY (hottest) for southern hemisphere
        ELSE
          HDAY = 200.0           !DOY (hottest) for northern hemisphere
        ENDIF
       
        TBD = 0.0
        TLL = 0.0
        TSW = 0.0
        TDL = 0.0
        CUMDPT = 0.0
        DO L = 1, NLAYR
!           BAK correction to have DSMID and CUMDPT in cm
!          DSMID(L) = CUMDPT + DLAYR(L)* 5.0
          DSMID(L) = CUMDPT + DLAYR(L)/5.0
!          CUMDPT   = CUMDPT + DLAYR(L)*10.0
          CUMDPT   = CUMDPT + DLAYR(L)
          TBD = TBD + BD(L)  * DLAYR(L)       !CHP
          TLL = TLL + LL(L)  * DLAYR(L)
          TSW = TSW + SWI(L) * DLAYR(L)
          TDL = TDL + DUL(L) * DLAYR(L)
        END DO

        IF (ISWWAT .EQ. 'Y') THEN
          PESW = AMAX1(0.0, TSW - TLL)      !cm
        ELSE
          !If water not being simulated, use DUL as water content
          PESW = AMAX1(0.0, TDL - TLL)
        ENDIF

        ABD    = TBD / DS(NLAYR)                   !CHP
        FX     = ABD/(ABD+686.0*EXP(-5.63*ABD))
        DP     = 1000.0 + 2500.0*FX
        WW     = 0.356  - 0.144*ABD
        B      = ALOG(500.0/DP)
        ALBEDO = MSALB

! CVF: difference in soil temperatures occur between different optimization
!     levels in compiled versions.
! Keep only 4 decimals. chp 06/03/03
!     Prevents differences between release & debug modes:
        DO I = 1, 5
!          TMA(I) = NINT(TAVG*10000.)/10000.   !chp
          TMA(I) = TAVG
          DelA(I) = 0.0
        END DO
        ATOT = TMA(1) * 5.0
        DelTOT = 0.0
        
        DO L = 1, NLAYR
          ST(L) = TAVG
        END DO
        ! also initialize SRFTEMP   BAK 2024 07 16
        SRFTEMP = TAVG

        DO I = 1, 8  ! spin 8 times
          CALL SOILTBK (
     &        ALBEDO,B,CUMDPT,DOY,DP,HDAY,NLAYR,           !Input
     &        PESW, SRAD, TAMP, TAV, TAVG, TMAX, WW, DSMID,!Input
!         added by BAK on 8 July 2024          
     &    BD,DLAYR,DS,DUL,LL,MSALB,CLAY,SILT,SAND,        !Input
     &    OC,SW,AVP,XHLAI,ES,EP,WINDSP,MULCHMASS,         !Input
     &    MULCHCOVER,MULCHTHICK,MULCHWAT,MULCHEVAP,       !Input
     &    MULCHALB,                            !Input
     &        ATOT, TMA, SRFTEMP, ST, DelA,DelTOT,        !Output
!            added by BAK 2023 11 29 for testing
     &    TA,DT,POR,                                      !Output
     &    SWREL,TcondDry, TcondSat, STCOND,HeatCap,       !Output
     &    DampDa,DampDw,CLAYFrac,SILTFrac,SANDFrac,OMFrac,
     &      ASTCOND,AHeatCap,                         !Output
     &    Del,STa,SolAvg,WINDmps,ESWatt,EPWatt,X1S,X1P, !Output
     &    X2S,X2P,XEPS,Xsky,Plht,
     &    AEROraSM1,AEROraS0,AEROraPM1,AEROraP0,
     &    RiNo,MEK,PHI,XG,
     &    XcubeS,XcubeP,DelS,DelP,J,M,
     &    Ga,GwSM1,GwS0,GwPM1,GwP0,TSM1,TS0,TPM1,TP0,
     &    RADSM1,RADS0,RADPM1,RADP0,HSM1,HS0,HPM1,HP0,
     &    FSM1,FS0,FPM1,FP0,FLAGS,FLAGP,CANHT,STEP,
     &    Sfrac,Pfrac,Mfrac,TM0,TMBot)                  !Output
          END DO
      ENDIF

!     Print soil temperature data in STEMP.OUT
          IF (MEEVP .NE. 'Z')    !CSVC
     & CALL OPSTEMPBK(CONTROL, ISWITCH, DOY, SRFTEMP, ST, TAV, TAMP,
!        added following outputs BAK 2023 11 29
     &   TMA,ATOT,TA,DT,
     &   DS,CLAY,SILT,SAND,OC,BD,SW,SWREL,POR,
     &     CLAYFrac,SILTFrac,SANDFrac,OMFrac,STCOND,HeatCap,
     &   TcondDry, TcondSat, ASTCOND,AHeatCap,DampDa,DampDw,
     &   Del,STa,SolAvg,WINDmps,ESWatt,EPWatt,X1S,X1P, !Output
     &    X2S,X2P,XEPS,Xsky,Plht,
     &    AEROraSM1,AEROraS0,AEROraPM1,AEROraP0,
     &    RiNo,MEK,PHI,XG,
     &    XcubeS,XcubeP,DelS,DelP,SRAD,WINDSP,ES,EP,TAVG,XHLAI,AVP,
     &    J,M,Ga,GwSM1,GwS0,GwPM1,GwP0,TSM1,TS0,TPM1,TP0,
     &    RADSM1,RADS0,RADPM1,RADP0,HSM1,HS0,HPM1,HP0,
     &    FSM1,FS0,FPM1,FP0,FLAGS,FLAGP,CANHT,STEP,
     &    Sfrac,Pfrac,Mfrac,TM0,TMBot)
!***********************************************************************
!***********************************************************************
!     Daily rate calculations
!***********************************************************************
      ELSEIF (DYNAMIC .EQ. RATE) THEN
!-----------------------------------------------------------------------
      TBD = 0.0
      TLL = 0.0
      TSW = 0.0
      DO L = 1, NLAYR
        TBD = TBD + BD(L) * DLAYR(L)
        TDL = TDL + DUL(L)* DLAYR(L)
        TLL = TLL + LL(L) * DLAYR(L)
        TSW = TSW + SW(L) * DLAYR(L)
      ENDDO

      ABD    = TBD / DS(NLAYR)                    !CHP
      FX     = ABD/(ABD+686.0*EXP(-5.63*ABD))
      DP     = 1000.0 + 2500.0*FX   !DP in mm
      WW     = 0.356  - 0.144*ABD   !vol. fraction
      B      = ALOG(500.0/DP)
      ALBEDO = MSALB

      IF (ISWWAT .EQ. 'Y') THEN
        PESW = MAX(0.0, TSW - TLL)      !cm
      ELSE
        !If water not being simulated, use DUL as water content
        PESW = AMAX1(0.0, TDL - TLL)    !cm
      ENDIF

      CALL SOILTBK (
     &    ALBEDO,B,CUMDPT,DOY,DP,HDAY,NLAYR,              !Input
     &    PESW, SRAD, TAMP, TAV, TAVG, TMAX, WW, DSMID,   !Input
!         added by BAK on 8 July 2024          
     &    BD,DLAYR,DS,DUL,LL,MSALB,CLAY,SILT,SAND,        !Input
     &    OC,SW,AVP,XHLAI,ES,EP,WINDSP,MULCHMASS,         !Input
     &    MULCHCOVER,MULCHTHICK,MULCHWAT,MULCHEVAP,       !Input
     &    MULCHALB,                           !Input
     &    ATOT, TMA, SRFTEMP, ST,DelA,DelTOT,             !Output
!            added by BAK 2023 11 29 for testing
     &    TA,DT,POR,                                      !Output
     &    SWREL,TcondDry, TcondSat, STCOND,HeatCap,       !Output
     &    DampDa,DampDw,CLAYFrac,SILTFrac,SANDFrac,OMFrac,
     &     ASTCOND,AHeatCap,                              !Output
     &    Del,STa,SolAvg,WINDmps,ESWatt,EPWatt,X1S,X1P,   !Output
     &    X2S,X2P,XEPS,Xsky,Plht,
     &    AEROraSM1,AEROraS0,AEROraPM1,AEROraP0,
     &    RiNo,MEK,PHI,XG,
     &    XcubeS,XcubeP,DelS,DelP,J,M,
     &    Ga,GwSM1,GwS0,GwPM1,GwP0,TSM1,TS0,TPM1,TP0,
     &    RADSM1,RADS0,RADPM1,RADP0,HSM1,HS0,HPM1,HP0,
     &    FSM1,FS0,FPM1,FP0,FLAGS,FLAGP,CANHT,STEP,      !Output
     &    Sfrac,Pfrac,Mfrac,TM0,TMBot)
!***********************************************************************
!***********************************************************************
!     Output & Seasonal summary
!***********************************************************************
      ELSEIF (DYNAMIC .EQ. OUTPUT .OR. DYNAMIC .EQ. SEASEND) THEN
!-----------------------------------------------------------------------
          IF (MEEVP .NE. 'Z')    !CSVC
     & CALL OPSTEMPBK(CONTROL, ISWITCH, DOY, SRFTEMP, ST, TAV, TAMP,
!        added following outputs BAK 2023 11 29
     &   TMA,ATOT,TA,DT,
     &   DS,CLAY,SILT,SAND,OC,BD,SW,SWREL,POR,
     &     CLAYFrac,SILTFrac,SANDFrac,OMFrac,STCOND,HeatCap,
     &   TcondDry, TcondSat, ASTCOND,AHeatCap,DampDa,DampDw,
     &   Del,STa,SolAvg,WINDmps,ESWatt,EPWatt,X1S,X1P, !Output
     &    X2S,X2P,XEPS,Xsky,Plht,
     &    AEROraSM1,AEROraS0,AEROraPM1,AEROraP0,
     &    RiNo,MEK,PHI,XG,
     &    XcubeS,XcubeP,DelS,DelP,SRAD,WINDSP,ES,EP,TAVG,XHLAI,AVP,
     &    J,M,Ga,GwSM1,GwS0,GwPM1,GwP0,TSM1,TS0,TPM1,TP0,
     &    RADSM1,RADS0,RADPM1,RADP0,HSM1,HS0,HPM1,HP0,
     &    FSM1,FS0,FPM1,FP0,FLAGS,FLAGP,CANHT,STEP,
     &    Sfrac,Pfrac,Mfrac,TM0,TMBot)
!***********************************************************************
!***********************************************************************
!     END OF DYNAMIC IF CONSTRUCT
!***********************************************************************
      ENDIF
!***********************************************************************
      RETURN
      END !SUBROUTINE STEMP
!=======================================================================


C=======================================================================
C  SOILTBK, Subroutine
C  Determines soil temperature by layer
C-----------------------------------------------------------------------
C  Revision history
C  02/09/1933 PWW Header revision and minor changes.
C  12/09/1999 CHP Revisions for modular format.
C  01/01/2000 AJG Added surface temperature for the CENTURY-based
C                SOM/soil-N module.
C  01/14/2005 CHP Added METMP = 3: Corrected water content in temp. eqn.
!  12/07/2008 CHP Removed METMP -- use only corrected water content
C-----------------------------------------------------------------------
C  Called : STEMP
C  Calls  : None
C=======================================================================

      SUBROUTINE SOILTBK (
     &    ALBEDO,B,CUMDPT,DOY,DP,HDAY,NLAYR,            !Input
     &    PESW, SRAD, TAMP, TAV, TAVG, TMAX, WW, DSMID, !Input
!         added by BAK on 8 July 2024          
     &    BD,DLAYR,DS,DUL,LL,MSALB,CLAY,SILT,SAND,      !Input
     &    OC, SW,AVP,XHLAI,ES,EP,WINDSP,MULCHMASS,      !Input
     &    MULCHCOVER,MULCHTHICK,MULCHWAT,MULCHEVAP,     !Input
     &    MULCHALB,                          !Input
     &    ATOT, TMA, SRFTEMP, ST,DelA,DelTOT,           !Output
!          added by BAK 2023 11 29 for testing
     &    TA,DT,POR,                                    !Output
     &    SWREL,TcondDry, TcondSat, STCOND,HeatCap,     !Output
     &    DampDa,DampDw,CLAYFrac,SILTFrac,SANDFrac,OMFrac,
     &       ASTCOND,AHeatCap,                          !Output
     &    Del,STa,SolAvg,WINDmps,ESWatt,EPWatt,X1S,X1P, !Output
     &    X2S,X2P,XEPS,Xsky,Plht,
     &    AEROraSM1,AEROraS0,AEROraPM1,AEROraP0,
     &    RiNo,MEK,PHI,XG,
     &    XcubeS,XcubeP,DelS,DelP,J,M,
     &    GaS,GwSM1,GwS0,GwPM1,GwP0,TSM1,TS0,TPM1,TP0,
     &    RADSM1,RADS0,RADPM1,RADP0,HSM1,HS0,HPM1,HP0,
     &    FSM1,FS0,FPM1,FP0,FLAGS,FLAGP,CANHT,STEP,
     &    Sfrac,Pfrac,Mfrac,TM0,TMBot)

!     ------------------------------------------------------------------
      USE ModuleDefs     !Definitions of constructed variable types,
                         ! which contain control information, soil
                         ! parameters, hourly weather data.
!     NL defined in ModuleDefs.for

      IMPLICIT  NONE
      SAVE

      INTEGER  K, L, DOY, NLAYR, J, M

      REAL ALBEDO, ALX, ATOT, B, CUMDPT, DD, DP, DT, FX
      REAL HDAY, PESW, SRAD, SRFTEMP, TA, TAMP, TAV, TAVG, TMAX
      REAL WC, WW, ZD,AVP,XHLAI
      REAL TMA(5),DelA(5),DelTOT
      REAL DSMID(NL),DlAYR(NL)
      REAL ST(NL)
      REAL SW(NL),DS(NL)
      REAL CLAYV(NL),SILTV(NL),SANDV(NL),OMV(NL),TotSolid(NL),POR(NL)
      REAL CLAYC(NL),SILTC(NL),SANDC(NL),OMC(NL),OM(NL),Totmineral(NL)
      REAL ClayFrac(NL), SiltFrac(NL), SandFrac(NL), OMFrac(NL)
      REAL TcondDry(NL), TcondSol(NL), TcondSat(NL), SWREL(NL)
      REAL  PX, QX, RX, SX
      REAL HeatCap(NL),STBot(NL),AMP(NL),DSMIDV2(NL)
      REAL BD(NL),CLAY(NL),SILT(NL),SAND(NL),OC(NL),STCOND(NL)
      REAL DUL(NL),LL(NL),MSALB
      REAL TSC,THC,DampDa,DampDw,Omega, Del,STa(NL),STboti(NL),AMPi(NL)
      REAL ASTCOND,AHeatCap, DDa(NL),DDw(NL)
      REAL STBZ,ALBS,ALBP,EPSS,EPSP,Rho,CP,Lamda,ES,EP
      REAL SolAvg,ESWatt,EPWatt,X1S,X1P,X2S,X2p,XEPS,Xsky,XG
      REAL XcubeS,XcubeP,DelS,DelP
      REAL PlHt,Z0,Disp,Zratio,RiNo,MEK,PHI,WINDSP,WINDmps,AEROra
      REAl GaS,GwSM1,GwS0,GwPM1,GwP0,TSM1,TS0,TPM1,TP0
      REAL RADSM1,RADS0,RADPM1,RADP0,HSM1,HS0,HPM1,HP0
      REAL FSM1,FS0,FPM1,FP0,FLAGS,FLAGP,CANHT,STEP
      REAL AEROraSM1,AEROraS0,AEROraPM1,AEROraP0
      REAL STbota(NL), AMPa(NL), AMPw(NL)
      REAL MULCHMASS,MULCHCOVER,MULCHTHICK,MULCHWAT,MULCHEVAP
      REAL MULCHALB,Pfrac,Sfrac,Mfrac,EPSM,EMWatt,X1M,GaM
      REAL MULCHgpcm3,MULCHcm3pcm3,MULCHpor,MULCHdryThCon,MPX,MQX,MSWREL
      REAL MULCHsatThCon,MSREL,MSX,MRX,MTCond,MHeatCap,DDMa,DDMw
      REAL TMM1,RADMM1,AEROraMM1,HMM1,GwMM1,FMM1
      REAL TM0,RADM0,AEROraM0,HM0,GwM0,FM0,DelM,TMBota,TMBot,SRFT

!-----------------------------------------------------------------------
!
      STBZ = 5.6697E-8 ! Stefan-Boltzmann constant [W/(m2 K)]
      ALBS = 0.2       ! albedo of soil; e.g.Campell,G.S. and Diac,G.R.
                       !  2005. p. 70 in G.L. Harfield and J.M. Baker,
                       ! Micrometeorology in Agricultural Sysatems,
                       ! Am. Soc.Agron., Crop Sci.Soc. Am., & Soil Sci.
                       ! Soc. Am, Madison, Wisconsin, USA.
      ALBP = 0.22      ! albedo and emissivity of plants; e.g. p.88-89
      EPSP = 0.98      ! in Monteith, J.
                       ! and Unsworth, M., 2008. Priciples of
                       ! Environmental Physics, Elsevier, Amsterdam.
      EPSS = 0.95      ! emissivity of soil; e.g. Idso et al. 1969.
                       ! Ecology 50(5):899-902.
      EPSM = 0.98      ! assume mulch emissivity same as plants
      Rho  = 1.204     ! Density of dry air (kg/m3) at sea level and 20C
      CP   = 1004.67   ! and heat capacity [J/(kg K)] in
                       ! Ham, J.H. 2005. p. 535 in G.L. Harfield
                       ! and J.M. Baker,
                       ! Micrometeorology in Agricultural Sysatems,
                       ! Am. Soc.Agron., Crop Sci.Soc. Am., & Soil Sci.
                       ! Soc. Am, Madison, Wisconsin, USA.
      Lamda = 2.501E6 - 2361*TAVG ! Latent heat of vaporization (J/kg),
                       ! also from Ham, p. 541
                       ! convert DOY to radians wi annual temp cycle        
      ALX = (FLOAT(DOY) - HDAY) * 0.0174
      
      !
!  Calculate fractonal area coverages of bare soil, mulch, and crop plants
      IF(XHLAI .GT. 3.0) THEN
          Pfrac = 1.0
      ELSE
          Pfrac = XHLAI/3.0
      END IF
          Mfrac = MULCHCOVER*(1.0 - Pfrac)
          Sfrac = 1.0 - Pfrac - Mfrac
      
      
!    11/28/2023 BAK Inserting Xiong (2023) alternative method for
!    simulating soil thermal conductivity and ultimatly damping depth.
!    [Xiong, K. et al. 2023. Scientific Reports.
!        https://doi.org/10.1038/s41598-023-37413-5]
      DO L = 1, NLAYR
!          correct % for presence of organic matter
          OM(L) = OC(L)/0.4 ! assume organic matter is 40% C
          Totmineral(L)=CLAY(L)+SILT(L)+SAND(L)
          IF(Totmineral(L)+OM(L) .GT. 100.) THEN
              CLAYC(L)=CLAY(L)*Totmineral(L)/(Totmineral(L)+OM(L))
              SILTC(L)=SILT(L)*Totmineral(L)/(Totmineral(L)+OM(L))
              SANDC(L)=SAND(L)*Totmineral(L)/(Totmineral(L)+OM(L))
              OMC(L)  =  OM(L)*Totmineral(L)/(Totmineral(L)+OM(L))
          ENDIF
!          convert from % by weight to cm3/cm3 volumes
!     2.65 is the particle density of sand, silt, and clay
!     1.3 is the partcle density of organic matter (DeVries 1963,1975)     
          CLAYV(L) = (CLAYC(L)/100.)*BD(L)/2.65     
          SILTV(L) = (SILTC(L)/100.)*BD(L)/2.65
          SANDV(L) = (SANDC(L)/100.)*BD(L)/2.65
          OMV(L)   = (OMC(L)/100)*BD(L)/1.3
          TotSolid(L) = CLAYV(L) + SILTV(L) + SANDV(L) + OMV(L)
          IF(TotSolid(L) .GT. 1.) THEN
              TotSolid(L) = 1.0
              END IF
          POR(L) = 1. - TotSolid(L)  ! Porosity
          END DO       
!
      DO L = 1, NLAYR
          ClayFrac(L) = CLAYV(L)/TotSolid(L)
          IF(ClayFrac(L) < 0.0) THEn
              ClayFrac(L) = 0.0
              END IF
          SiltFrac(L) = SILTV(L)/TotSolid(L)
          IF(SiltFrac(L) < 0.0) THEN
              SiltFrac(L) = 0.0
              END IF         
          SandFrac(L) = SANDV(L)/TotSolid(L)
          IF(SandFrac(L) < 0.0) THEN
              SandFrac(L) = 0.0
              END IF
          OMFrac(L)   = OMV(L)/TotSolid(L)
          IF(OMFrac(L) < 0.0) THEn
              OMFrac(L) = 0.0
              END IF
      END DO
   
!
      DO L = 1, NLAYR
!
! Calculate dry thermal conductivity (W m-1 C-1)
!     Equation was fitted to mineral soils
       TCondDry(L) = -0.6*POR(L) + 0.51
        PX = TCondDry(L)
!       
! Geometric mean thermal conductivity of solid materials for mineral soils
! from Xiong et al (2023)
        If(SandFrac(L) .GT. 0.2) THEN            
           TcondSol(L) = (7.7**SandFrac(L))*
     &         (2.0**(SiltFrac(L) + ClayFrac(L)))*(0.25**OMFrac(L))
        ELSE
         TcondSol(L) = (7.7**SandFrac(L))*
     &         (3.0**(SiltFrac(L) + ClayFrac(L)))*(0.25**OMFrac(L))
        END IF
! where 7.7 (W m-1 C-1) = thermal conductivity of quartz (sand) and
!       2.0 (W m-1 C-1) = thermal conductivity of other soil minerals
!               if sand > 0.2; 
!       3.0 otherwise
!       0.25 thermal conductivity of organic matter (DeVries 1963, 1975)
!
! Geometric mean thermal conductivity of soil solids and water at saturation
        TcondSat(L) = (TcondSol(L)**(1. - POR(L)))*(0.594**POR(L))
!   where 0.594 W m-1 C-1 is the thermal conductivity of wqter at 20C (DeVries 1963, 1975)
        QX = TcondSat(L) - TcondDry(L)
      
! Relative soil water content (cm3/cm3) compared to saturation when pores are full of water
        SWREL(L) = SW(L)/POR(L)
        IF(SWREL(L) .LT. 0.00001) SWREL(L) = 0.00001
        IF(SWREL(L) .GT. 1.0) SWREL(L) = 1.0
         SX = 1.5*(SWREL(L) - SWREL(L)**2)
!
!     R vs Sand fitted by BAK to Table 1 of Xiong et al. (2023)         
         RX = -1.2125*SANDV(L) + 1.8935
!
! Compute thermal conductivity of the soil (W m-1 C-1)
        STCond(L) = PX + QX*(SWREL(L)**RX)
     &       + SX*EXP(SWREL(L)*(1. - SWREL(L)))
!
!   Calculate Heat Capacity (J m-3 C-1) following DeVries (1963, 1975)
        HeatCap(L) = (SANDV(L)+SILTV(L)+CLAYV(L))*2.0E6 + OMV(L)*2.5E6
     &                  + SW(L)*4.2E6
        END DO
!
!    Calcualte average thermal conductivity and heat capacity for
!      whole soil profile
	TSC = 0.0
	THC = 0.0
	DO L = 1,NLAYR
		TSC = TSC + STCond(L)*DLAYR(L)
          THC = THC + HeatCap(L)*DLAYR(L)
      END DO
      ASTCond = TSC/DS(NLAYR)
      AHeatCap = THC/DS(NLAYR)
        
!     Calculate avg annual and weather front Damping depths (cm)
        Omega = 2.0*3.14159/(365.0*24.0*3600.0)    ! radians/s
        DampDa = 100.*SQRT(2.*ASTCOND/(AHeatCap*Omega)) ! annual
        Omega = 2.0*3.14159/(5.0*24.0*3600.0)    ! radians/s
        DampDw = 100.*SQRT(2.*ASTCOND/(AHeatCap*Omega)) ! 5 day
        
!  Calculte annual and weather damping depths for each layer (cm)
        DO L = 1,NLAYR
         Omega = 2.0*3.14159/(365.0*24.0*3600.0)    ! rad/s annual
         DDa(L) = 100.*SQRT(2.*STCond(L)/(HeatCap(L)*Omega)) ! annual
         Omega = 2.0*3.14159/(5.0*24.0*3600.0)    ! rad/s 5 day
         DDw(L) = 100.*SQRT(2.*STCond(L)/(HeatCap(L)*Omega)) ! annual
      END DO
!
!
!   Get ready for energy balances
!         
!     Calculate average solar rad for day From MJ/(m2 day) to W/m2
      SolAvg = SRAD*1.0E6/(24.0*3600.0)
!     Convert wind in km/day to average m/s
      WINDmps= WINDSP*0.01157
!     Convert ES and EP from mm/day to W/m2
      ESWatt = ES*28.36
      EPWatt = EP*28.36
      EMWatt = MULCHEVAP*28.36
!
!      Calculate sky radiation term following Prata, A.I.,
!      1996. Q.J.R. Meteorological Soc. 122:1127-1151,
!      doi:10.1002/qj.49712253306
      XEPS = 465.*AVP/(TAVG+273.15)
      Xsky = (1. - (1. + XEPS)*EXP(-SQRT(1.2 + 3.0*XEPS)))*
     &         STBZ*(TAVG+273.15)**4
!      
!     Caculate net solar for mulch, soil, and plants
      X1M = SolAVG*(1. - MULChALB)
      X1S = SolAvg*(1. - ALBS)
      X1P = SolAvg*(1. - ALBP)
     
!     Compute surface soil flux for annual wave per Novak (2005, Eq. 6;
!     pp. 105-129 in J.L. Hatfield and J.M. Baker. 2005. Micrometeorology
!     in Agricultural Systems, #47 in Agronomy Series, Am. Soc.Agron.,
!     Crop Sci. Soc. Am., and Soil Sci. Sci. Am., Madison, WI)
!     or Kimball and Jackson (1979, Eq. 3.4-18; pp. 211-229 in B.J.
!     Barfield and J.F. Gerber, 1979, Modification of the Environment
!     of Plants, ASAE Monograph, AM. Soc. Ag. Eng. St. Joseph, MI)
      GaS = (TAMP/2.)*SQRT((2.0*3.14159/(365.0*24.0*3600.0))*
     &     STCond(1)*HeatCap(1))*COS(ALX + 3.1416/4.)
!
!     Similarly compute surface mulch heat flux for annual wave
      GaM = (TAMP/2.)*SQRT((2.0*3.14159/(365.0*24.0*3600.))*
     &  MTCond*MheatCap)*COS(ALX + 3.1416/4.)
            
!     Check plant height
      If(CANHT < 0.01) THEN
          PlHt = 0.01 ! assume bare soil and mulch have roughness 
!                      elements 1 cm high
      ELSE
          PlHt = CANHT
      END IF
      Z0 = 0.13*PlHt  ! roughness length from Monteith
      Disp = 0.63*PlHt ! displacement height from Monteith
      Zratio = (PlHt + 1. - Disp + Z0)/Z0  ! assume wind measured at
!        1.0 m above whatever the plant height is
      
!  *********** START MULCH SECTION ******************
! ** Compute thermal conductivity of mulch **
!  following the soil thermal conductivity model of Xiong et al. (2023)
      
      IF(Mfrac .GT. 0.001) THEN 

!    Compute mulch density (kg/ha to g/cm3)
      IF(MULCHTHICK .LT. 0.0001) THEN
          MULCHgpcm3 = 0.0
          MULCHcm3pcm3 = 0.0
          MULCHpor = 1.0
      ELSE
!     MULCHMASS in kg/ha and MULCHTHICK in mm
      MULCHgpcm3 = MULCHMASS*1.E-3/(MULCHTHICK*0.1*1.E8)
! 1.3 is the partcle density of organic matter (DeVries 1963, 1975)      
      MULCHcm3pcm3 = MULCHgpcm3/1.3
      MULCHpor     = 1. - MULCHcm3pcm3 ! porosity
          END IF           
!
!     Dry thermal conductivity
!  0.25 (W/(m C) thermal conductivity of organic matter (DeVries 1963,1975)
! Assume it decreases linearly to 0.0 with porosity going from 1.0 to 0.0
      MULCHdryThCon = 0.25 - 0.25*MULCHpor
      MPX = MULCHdryThCon
      
!     Saturated thermal conductivity
! Geometric mean thermal conductivity of soil solids and water at saturation
        MULCHsatThCon = (0.25**(1. - MULCHpor))*(0.594**MULCHpor)
!   where 0.594 W m-1 C-1 is the thermal conductivity of wqter at 20C (DeVries 1963,1975)
        MQX = MULCHsatThCon - MULCHdryThCon
                
! Relative soil water content (cm3/cm3) compared to saturation
        MSWREL = MULCHWAT/MULCHpor
        IF(MSWREL .LT. 0.00001) MSWREL = 0.00001
        IF(MSWREL .GT. 1.0) MSWREL = 1.0
        MSX = 1.5*(MSWREL - MSWREL**2)
        
!     R vs Sand fitted by BAK to Table 1 of Xiong et al. (2023)
!     However, all the soils in Table 1 are mineral soils so may
!     not be correct for organic matter.
!        Taking the value for zero sand
         MRX = 1.8935
         
! Compute thermal conductivity of the mulch (W m-1 C-1)
        MTCond = MPX + MQX*(MSWREL**MRX)
     &       + MSX*EXP(MSWREL*(1. - MSWREL))
!
!   Calculate Mulch Heat Capacity (J m-3 C-1) following DeVries (1963,1975)
        MHeatCap = MULCHcm3pcm3*2.5E6 + MULCHWAT*4.2E6
        
!  Calculte annual and weather damping depths for mulch layer (cm)
         Omega = 2.0*3.14159/(365.0*24.0*3600.0)    ! rad/s annual
         DDMa = 100.*SQRT(2.*MTCond/(MHeatCap*Omega)) ! annual
         Omega = 2.0*3.14159/(5.0*24.0*3600.0)    ! rad/s 5 day
         DDMw = 100.*SQRT(2.*MTCond/(MHeatCap*Omega)) ! annual

      
!      **************** Start Mulch Surface temperature loop ****************
		
				TMM1 = TAVG ! Initial guess for mulch surfaCE temperature


			RADMM1=EPSS*STBZ*(TMM1+273.15)**4	! Upwelling canopy radiation
!     Calculate aerodynamic resistance following Kimball et al.,
!     2015. Agronomy J. 107(1):129-141.
!     doi:10.2134/agronj14.0109
!     using Mahrt and Ek, 1984. J. Clim. Appl. Meteorology
!     23:222-234.
!      doi:10.1175/1520-0450(1984)023<0222:TIOASO>2.0.CO;2
      

      IF(WINDmps .LT. 0.1 .AND. ABS(TMM1 - TAVG) .LT. 0.1) THEN
          AEROraMM1 = RHO*CP/2.32
          ELSE IF(WINDmps .LT. 0.1 .AND. ABS(TSM1 - TAVG) .GE. 0.1) THEN
              AEROraMM1 = RHO*CP/(5.*(ABS(TMM1 - TAVG))**0.33)
      ELSE
          
      PlHt = 0.01 ! assume bare soil and mulch have roughness 
!                      elements 1 cm high
      Z0 = 0.13*PlHt  ! roughness length from Monteith
      Disp = 0.63*PlHt ! displacement height from Monteith
      Zratio = (PlHt + 1. - Disp + Z0)/Z0  ! assume wind measured at
!        2.0 m 
            RiNo = 9.8*(TAVG-TMM1)*(2.0-Disp)/
     &            ((TAVG+273.15)*WINDmps**2.)
!         where 9.8 is acceleration of gravity
!         and assume reference 1m above PlHt
!     Mahrt & Ek K  where 0.4 is von Karmen's constant
           MEK = 75.*0.4*0.4*SQRT(Zratio)/(LOG(Zratio))**2
          IF(TMM1 < TAVG) THEN
           PHI = (1. + 15.*RiNo)*SQRT(1. + 5.*RiNo)  ! stable conditions
              ELSE
!             unstable conditions          
          PHI = 1./(1. - (15.*RiNo)/(1. + MEK*SQRT(-RiNo)))
              END IF
      AEROraMM1 = ((1./WINDmps)*((1./0.4)*LOG(Zratio))**2.)*PHI
      END IF

!         Calcualte sensible heat
      HMM1 = (Rho*CP/AEROraMM1)*(TMM1-TAVG)
      
!         Calculate "weather" soil heat flux wave ! 5 day
      GwMM1 = (TMM1 - TAVG)*SQRT((2.0*3.14159/(5.0*24.0*3600.0))*
     &     MTCond*MHeatCap)
      
			FMM1 = -X1S +RADMM1 -XSky +HMM1 +GaM +GwMM1 +EMWatt



!		*** Start MULCH surface temperature iteration loop ***
		STEP=SIGN(1.,TMM1)
		TM0=TMM1 + STEP

		FLAGS=0

		DO j = 1, 10000
		
              RADM0=EPSM*STBZ*(TM0+273.15)**4
!      Upwelling canopy radiation
!     Calculate aerodynamic resistance following Kimball et al.,
!     2015. Agronomy J. 107(1):129-141.
!     doi:10.2134/agronj14.0109
!     using Mahrt and Ek, 1984. J. Clim. Appl. Meteorology
!     23:222-234.
!      doi:10.1175/1520-0450(1984)023<0222:TIOASO>2.0.CO;2     

      IF(WINDmps .LT. 0.1 .AND. ABS(TM0 - TAVG) .LT. 0.1) THEN
          AEROraM0 = RHO*CP/2.32
          ELSE IF(WINDmps .LT. 0.1 .AND. ABS(TM0 - TAVG) .GE. 0.1) THEN
              AEROraM0 = RHO*CP/(5.*(ABS(TM0 - TAVG))**0.33)
          ELSE
!          !      Richardson No.
            RiNo = 9.8*(TAVG-TM0)*(2.0-Disp)/
     &            ((TAVG+273.15)*WINDmps**2.)
!         where 9.8 is acceleration of gravity
!         and assume reference 1m above PlHt
!     Mahrt & Ek K  where 0.4 is von Karmen's constant
           MEK = 75.*0.4*0.4*SQRT(Zratio)/(LOG(Zratio))**2
          IF(TM0 < TAVG) THEN
           PHI = (1. + 15.*RiNo)*SQRT(1. + 5.*RiNo)  ! stable conditions
              ELSE
!             unstable conditions          
          PHI = 1./(1. - (15.*RiNo)/(1. + MEK*SQRT(-RiNo)))
              END IF
      AEROraM0 = ((1./WINDmps)*((1./0.4)*LOG(Zratio))**2.)*PHI
      END IF

      !         Calcualte sensible heat
      HM0 = (Rho*CP/AEROraM0)*(TM0-TAVG)
      
!         Calculate "weather" mulch heat flux wave
      GwM0 = (TM0 - TAVG)*SQRT((2.0*3.14159/(5.0*24.0*3600.0))* ! 5 day
     &       MTCond*MheatCap)
      
			FM0 = -X1S +RADM0 -XSky +HM0 +GaM +GwM0 +EMWatt



          IF(ABS(FM0)<0.0001 .OR. ABS(FM0-FMM1)<0.00001 
     &            .OR. ABS(TS0-TMM1)<0.0001
     &	        .AND. ABS(HM0-HMM1)<0.0001 .AND. ABS(GwM0-GwMM1)<0.0001
     &			.AND. ABS(RADM0-RADMM1)<0.0001
     &			.AND. j>2) THEN			! Have convergence											
				EXIT									

				ELSE IF(ABS(FM0-FMM1)< 1.D-12) THEN
				! need to avoid divide by zero, so assume convergence		
				EXIT
		
				ELSE IF(SIGN(1.,FM0)/=SIGN(1.,FMM1)) THEN	
						STEP=-0.5*STEP
                          
                  ELSE IF(SIGN(1.,FM0)==SIGN(1.,FMM1) .AND.
     &              ABS(FM0)>ABS(FMM1)) THEN !going wrong way,
                          STEP=-STEP         !  so need to go back
                          IF(ABS(STEP)<1.E-8) THEN
                              STEP=0.5*(TMM1+TM0)
                              END IF
                  END IF
                  
                  TMM1 = TM0
				TM0 = TM0 + STEP
                  FMM1 = FM0
                  AEROraMM1=AEROraM0
                  HMM1 = HM0
                  GwMM1 = GwM0
                  RADMM1 = RADM0
                  

!				IF(IOUT==1) THEN
!					WRITE(7,105) j, TA(i), TSM1, TS0, TC1, SVPTC, FM1, FK0, &
!					 RADCM1, RADC0, RADNETM1, RADNET0, &
!					 RichNo(i), PHI(i), KCON, HM1, EM1, HvM1, &
!					 RAM1, RAK(i), H0, E0, Hv0
!					 105 FORMAT(' ', I5, 5F9.3, 17E14.5)
!					END IF

	
			IF(j==9999) THEN	! convergence not achieved. Set flag and
							! use air temperature as soil surface temperature
				TM0 = TAVG
				END IF

              ENDDO	! *********** End of Mulch Iteration Loop **********
              DelM = TM0 - TAVG
              
!     Compute temperature at bottom of mulch layer
!       For annual temperature wave
        ZD = -(MULCHTHICK/10)/DDMa
        TMBota = TAV + ((TAMP/2.0)*COS(ALX+ZD))*EXP(ZD)
! Add attenuated weather + surface deviaiton from ann air temp        
        ZD = -(MULCHTHICK/10)/DDMw
        TMBot = TMBota + DelM*EXP(ZD)
        
      END IF   ! End of Mulch Section

!      **************** Start Soil Surface temperature loop ****************
	IF (Sfrac .GT. 0.001) THEN	
				TSM1 = TAVG ! Initial guess for soil surfaCE temperature


			RADSM1=EPSS*STBZ*(TSM1+273.15)**4	! Upwelling canopy radiation
!     Calculate aerodynamic resistance following Kimball et al.,
!     2015. Agronomy J. 107(1):129-141.
!     doi:10.2134/agronj14.0109
!     using Mahrt and Ek, 1984. J. Clim. Appl. Meteorology
!     23:222-234.
!      doi:10.1175/1520-0450(1984)023<0222:TIOASO>2.0.CO;2
      

      IF(WINDmps .LT. 0.1 .AND. ABS(TSM1 - TAVG) .LT. 0.1) THEN
          AEROraSM1 = RHO*CP/2.32
          ELSE IF(WINDmps .LT. 0.1 .AND. ABS(TSM1 - TAVG) .GE. 0.1) THEN
              AEROraSM1 = RHO*CP/(5.*(ABS(TSM1 - TAVG))**0.33)
          ELSE
!          !      Richardson No.
            RiNo = 9.8*(TAVG-TSM1)*(PlHt+1.0-Disp)/
     &            ((TAVG+273.15)*WINDmps**2.)
!         where 9.8 is acceleration of gravity
!         and assume reference 1m above PlHt
!     Mahrt & Ek K  where 0.4 is von Karmen's constant
           MEK = 75.*0.4*0.4*SQRT(Zratio)/(LOG(Zratio))**2
          IF(TSM1 < TAVG) THEN
           PHI = (1. + 15.*RiNo)*SQRT(1. + 5.*RiNo)  ! stable conditions
              ELSE
!             unstable conditions          
          PHI = 1./(1. - (15.*RiNo)/(1. + MEK*SQRT(-RiNo)))
              END IF
      AEROraSM1 = ((1./WINDmps)*((1./0.4)*LOG(Zratio))**2.)*PHI
      END IF

!         Calcualte sensible heat
      HSM1 = (Rho*CP/AEROraSM1)*(TSM1-TAVG)
      
!         Calculate "weather" soil heat flux wave ! 5 day
      GwSM1 = (TSM1 - TAVG)*SQRT((2.0*3.14159/(5.0*24.0*3600.0))*
     &     STCond(1)*HeatCap(1))
      
			FSM1 = -X1S +RADSM1 -XSky +HSM1 +GaS +GwSM1 +ESWatt



!		*** Start Soil surface temperature iteration loop ***
		STEP=SIGN(1.,TSM1)
		TS0=TSM1 + STEP

		FLAGS=0

		DO j = 1, 10000
		
              RADS0=EPSS*STBZ*(TS0+273.15)**4
!      Upwelling canopy radiation
!     Calculate aerodynamic resistance following Kimball et al.,
!     2015. Agronomy J. 107(1):129-141.
!     doi:10.2134/agronj14.0109
!     using Mahrt and Ek, 1984. J. Clim. Appl. Meteorology
!     23:222-234.
!      doi:10.1175/1520-0450(1984)023<0222:TIOASO>2.0.CO;2     

      IF(WINDmps .LT. 0.1 .AND. ABS(TS0 - TAVG) .LT. 0.1) THEN
          AEROraS0 = RHO*CP/2.32
          ELSE IF(WINDmps .LT. 0.1 .AND. ABS(TS0 - TAVG) .GE. 0.1) THEN
              AEROraS0 = RHO*CP/(5.*(ABS(TS0 - TAVG))**0.33)
          ELSE
!          !      Richardson No.
            RiNo = 9.8*(TAVG-TS0)*(PlHt+1.0-Disp)/
     &            ((TAVG+273.15)*WINDmps**2.)
!         where 9.8 is acceleration of gravity
!         and assume reference 1m above PlHt
!     Mahrt & Ek K  where 0.4 is von Karmen's constant
           MEK = 75.*0.4*0.4*SQRT(Zratio)/(LOG(Zratio))**2
          IF(TS0 < TAVG) THEN
           PHI = (1. + 15.*RiNo)*SQRT(1. + 5.*RiNo)  ! stable conditions
              ELSE
!             unstable conditions          
          PHI = 1./(1. - (15.*RiNo)/(1. + MEK*SQRT(-RiNo)))
              END IF
      AEROraS0 = ((1./WINDmps)*((1./0.4)*LOG(Zratio))**2.)*PHI
      END IF

      !         Calcualte sensible heat
      HS0 = (Rho*CP/AEROraS0)*(TS0-TAVG)
      
!         Calculate "weather" soil heat flux wave
      GwS0 = (TS0 - TAVG)*SQRT((2.0*3.14159/(5.0*24.0*3600.0))* ! 5 day
     &     STCond(1)*HeatCap(1))
      
			FS0 = -X1S +RADS0 -XSky +HS0 +GaS +GwS0 +ESWatt



          IF(ABS(FS0)<0.0001 .OR. ABS(FS0-FSM1)<0.00001 
     &            .OR. ABS(TS0-TSM1)<0.0001
     &	        .AND. ABS(HS0-HSM1)<0.0001 .AND. ABS(GwS0-GwSM1)<0.0001
     &			.AND. ABS(RADS0-RADSM1)<0.0001
     &			.AND. j>2) THEN			! Have convergence											
				EXIT									

				ELSE IF(ABS(FS0-FSM1)< 1.D-12) THEN
				! need to avoid divide by zero, so assume convergence		
				EXIT
		
				ELSE IF(SIGN(1.,FS0)/=SIGN(1.,FSM1)) THEN	
						STEP=-0.5*STEP
                          
                  ELSE IF(SIGN(1.,FS0)==SIGN(1.,FSM1) .AND.
     &              ABS(FS0)>ABS(FSM1)) THEN !going wrong way,
                          STEP=-STEP         !  so need to go back
                          IF(ABS(STEP)<1.E-8) THEN
                              STEP=0.5*(TSM1+TS0)
                              END IF
                  END IF
                  
                  TSM1 = TS0
				TS0 = TS0 + STEP
                  FSM1 = FS0
                  AEROraSM1=AEROraS0
                  HSM1 = HS0
                  GwSM1 = GwS0
                  RADSM1 = RADS0
                  

!				IF(IOUT==1) THEN
!					WRITE(7,105) j, TA(i), TSM1, TS0, TC1, SVPTC, FM1, FK0, &
!					 RADCM1, RADC0, RADNETM1, RADNET0, &
!					 RichNo(i), PHI(i), KCON, HM1, EM1, HvM1, &
!					 RAM1, RAK(i), H0, E0, Hv0
!					 105 FORMAT(' ', I5, 5F9.3, 17E14.5)
!					END IF

	
			IF(j==9999) THEN	! convergence not achieved. Set flag and
							! use air temperature as soil surface temperature
				TS0 = TAVG
				END IF

              END DO	! *********** End of Soil Surface Loop **********!
      END IF ! End of soil section

              
!      ******* Start Plant Surface temperature loop ****************
	IF(Pfrac .GT. 0.001) THEN
          
				TPM1 = TAVG		! Initial guess for Plant temperature


			RADPM1=EPSP*STBZ*(TPM1+273.15)**4	! Upwelling canopy radiation
!     Calculate aerodynamic resistance following Kimball et al.,
!     2015. Agronomy J. 107(1):129-141.
!     doi:10.2134/agronj14.0109
!     using Mahrt and Ek, 1984. J. Clim. Appl. Meteorology
!     23:222-234.
!      doi:10.1175/1520-0450(1984)023<0222:TIOASO>2.0.CO;2
      

      IF(WINDmps .LT. 0.1 .AND. ABS(TPM1 - TAVG) .LT. 0.1) THEN
          AEROraPM1 = RHO*CP/2.32
          ELSE IF(WINDmps .LT. 0.1 .AND. ABS(TPM1 - TAVG) .GE. 0.1) THEN
              AEROraPM1 = RHO*CP/(5.*(ABS(TPM1 - TAVG))**0.33)
          ELSE
!          !      Richardson No.
            RiNo = 9.8*(TAVG-TPM1)*(PlHt+1.0-Disp)/
     &            ((TAVG+273.15)*WINDmps**2.)
!         where 9.8 is acceleration of gravity
!         and assume reference 1m above PlHt
!     Mahrt & Ek K  where 0.4 is von Karmen's constant
           MEK = 75.*0.4*0.4*SQRT(Zratio)/(LOG(Zratio))**2
          IF(TPM1 < TAVG) THEN
           PHI = (1. + 15.*RiNo)*SQRT(1. + 5.*RiNo)  ! stable conditions
              ELSE
!             unstable conditions          
          PHI = 1./(1. - (15.*RiNo)/(1. + MEK*SQRT(-RiNo)))
              END IF
      AEROraPM1 = ((1./WINDmps)*((1./0.4)*LOG(Zratio))**2.)*PHI
      END IF

!         Calcualte sensible heat
      HPM1 = (Rho*CP/AEROraPM1)*(TPM1-TAVG)
      
!         Calculate "weather" plant heat flux wave ! 5 day
      GwPM1 = (TPM1 - TAVG)*SQRT((2.0*3.14159/(5.0*24.0*3600.0))*
     &     STCond(1)*HeatCap(1))
      
			FPM1 = -X1P +RADPM1 -XSky +HPM1 +GaS +GwPM1 +EPWatt



!		*** Start Plant surface temperature iteration loop ***
		STEP=SIGN(1.,TPM1)
		TP0=TPM1 + STEP

		FLAGP=0

		DO M = 1, 10000
		
              RADP0=EPSP*STBZ*(TP0+273.15)**4
              ! Upwelling canopy radiation
!     Calculate aerodynamic resistance following Kimball et al.,
!     2015. Agronomy J. 107(1):129-141.
!     doi:10.2134/agronj14.0109
!     using Mahrt and Ek, 1984. J. Clim. Appl. Meteorology
!     23:222-234.
!      doi:10.1175/1520-0450(1984)023<0222:TIOASO>2.0.CO;2     

      IF(WINDmps .LT. 0.1 .AND. ABS(TP0 - TAVG) .LT. 0.1) THEN
          AEROraP0 = RHO*CP/2.32
          ELSE IF(WINDmps .LT. 0.1 .AND. ABS(TP0 - TAVG) .GE. 0.1) THEN
              AEROraP0 = RHO*CP/(5.*(ABS(TP0 - TAVG))**0.33)
          ELSE
!          !      Richardson No.
            RiNo = 9.8*(TAVG-TP0)*(PlHt+1.0-Disp)/
     &            ((TAVG+273.15)*WINDmps**2.)
!         where 9.8 is acceleration of gravity
!         and assume reference 1m above PlHt
!     Mahrt & Ek K  where 0.4 is von Karmen's constant
           MEK = 75.*0.4*0.4*SQRT(Zratio)/(LOG(Zratio))**2
          IF(TP0 < TAVG) THEN
           PHI = (1. + 15.*RiNo)*SQRT(1. + 5.*RiNo)  ! stable conditions
              ELSE
!             unstable conditions          
          PHI = 1./(1. - (15.*RiNo)/(1. + MEK*SQRT(-RiNo)))
              END IF
      AEROraP0 = ((1./WINDmps)*((1./0.4)*LOG(Zratio))**2.)*PHI
      END IF

      !         Calcualte sensible heat
      HP0 = (Rho*CP/AEROraP0)*(TP0-TAVG)
      
!         Calculate "weather" plant heat flux wave
      GwP0 = (TP0 - TAVG)*SQRT((2.0*3.14159/(5.0*24.0*3600.0))* ! 5 day
     &     STCond(1)*HeatCap(1))
      
			FP0 = -X1P +RADP0 -XSky +HP0 +GaS +GwP0 +EPWatt



			IF(ABS(FP0)<0.0001 .OR. ABS(FP0-FPM1)<0.00001 
     &            .OR. ABS(TP0-TPM1)<0.0001
     &	        .AND. ABS(HP0-HPM1)<0.0001 .AND. ABS(GwP0-GwPM1)<0.0001
     &			.AND. ABS(RADP0-RADPM1)<0.0001
     &			.AND. M>2) THEN			! Have convergence											
				EXIT									

				ELSE IF(ABS(FP0-FPM1)< 1.D-12) THEN
				! need to avoid divide by zero, so assume convergence		
				EXIT
		
                  ELSE IF(SIGN(1.,FP0)/=SIGN(1.,FPM1)) THEN	
						STEP=-0.5*STEP
                          
                  ELSE IF(SIGN(1.,FP0)==SIGN(1.,FPM1) .AND.
     &              ABS(FP0)>ABS(FPM1)) THEN !going wrong way so,
                          STEP=-STEP         ! need to go back
                          IF(ABS(STEP)<1.E-8) THEN
                              STEP=0.5*(TPM1+TP0)
                              END IF
                  END IF
                  
                  TPM1 = TP0
				TP0 = TP0 + STEP
                  FPM1 = FP0
                  AEROraPM1=AEROraP0
                  HPM1 = HP0
                  GwPM1 = GwP0
                  RADPM1 = RADP0
                  

!				IF(IOUT==1) THEN
!					WRITE(7,105) j, TA(i), TSM1, TS0, TC1, SVPTC, FM1, FK0, &
!					 RADCM1, RADC0, RADNETM1, RADNET0, &
!					 RichNo(i), PHI(i), KCON, HM1, EM1, HvM1, &
!					 RAM1, RAK(i), H0, E0, Hv0
!					 105 FORMAT(' ', I5, 5F9.3, 17E14.5)
!					END IF

	
			IF(M==9999) THEN	! convergence not achieved. Set flag and
							! use air temperature as plant surface temperature
				TP0 = TAVG
				END IF

              END DO	! ******** End of Plant Surface Loop **********
              
          END IF    ! End of plant section      
              
!      Calculate weighted average differece between mulch, soil surface, 
!        and plant surface  and the average air temperature for today
          SRFT = Sfrac*TS0 + Pfrac*TP0 + Mfrac*TMBot
!      Calculate standard annual temperature for the day
          TA = TAV + (TAMP/2.0) * COS(ALX) ! air temp from annual curve   
          Del = SRFT - TA
          
      
!     Compute average air temp for last 5 days    
      ATOT   = ATOT - TMA(5)
      DelTOT = DelTOT - DelA(5)

      DO K = 5, 2, -1
        TMA(K) = TMA(K-1)
        DelA(K)= DelA(K-1)
      END DO
!
!      Get rid of solar radiation stuff and just use TAVG      
!      TMA(1) = (1.0 - ALBEDO) * (TAVG + (TMAX - TAVG) *
!     &      SQRT(SRAD * 0.03)) + ALBEDO * TMA(1)
       TMA(1) = TAVG
       DelA(1)= Del
!     Instead of using air temperature alone or modified
!      by solar radiation, use the soil surface temperature
!      averaged over the last 5 days
       
!     Prevents differences between release & debug modes:
!       Keep only 4 decimals. chp 06/03/03
      TMA(1) = NINT(TMA(1)*10000.)/10000.  !chp       
        ATOT = ATOT + TMA(1)
        DelTOT = DelTOT + DelA(1)
        DT = DelTOT/5.
        SRFTEMP = (ATOT/5.) + DT


!-----------------------------------------------------------------------
!      !Water content function - compare old and new
!      SELECT CASE (METMP)
!      CASE ('O')  !Old, uncorrected equation
!        !OLD EQUATION (used in DSSAT v3.5 CROPGRO, SUBSTOR, CERES-Maize
!         WC = AMAX1(0.01, PESW) / (WW * CUMDPT * 10.0)
!
!      CASE ('E')  !Corrected method (EPIC)
!        !NEW (CORRECTED) EQUATION
!        !chp 11/24/2003 per GH and LAH
!        WC = AMAX1(0.01, PESW) / (WW * CUMDPT) * 10.0
!     frac =              cm   / (    mm     ) * mm/cm
        !WC (ratio)
        !PESW (cm)
        !WW (dimensionless)
        !CUMDPT (mm)
!      END SELECT
!-----------------------------------------------------------------------

!      FX = EXP(B * ((1.0 - WC) / (1.0 + WC))**2)

!      DD = FX * DP                                  !DD in mm
!     JWJ, GH 12/9/2008
!     Checked damping depths against values from literature and
!       values are reasonable (after fix to WC equation).
!     Hillel, D. 2004. Introduction to Environmental Soil Physics.
!       Academic Press, San Diego, CA, USA.
!
!   *** Calcualte soil temperatures for "ideal" annual cosine curve
        ! using average annual damping depth for whole profile
!     DO L = 1, NLAYR
!        ZD    = -DSMID(L) / DampDa
!        STa(L) = TAV + ((TAMP/2.0) * COS(ALX + ZD)) * EXP(ZD)
!     END DO
           
!   *** calculate soil temperature accounting for weather fronts
      ! using 5-day damping depth averaged for whole profile
!      DO L = 1, NLAYR
!          ZD = -DSMID(L)/DampDw
!          ST(L) = STa(L) + DT*EXP(ZD)
!         ST(L) = NINT(ST(L) * 1000.) / 1000. !debug vs release fix
!     END DO
      
      !     Added: soil T for surface litter layer.
!     NB: this should be done by adding array element 0 to ST(L). Now
!     temporarily done differently.
!      SRFTEMP = TAV + (TAMP / 2. * COS(ALX) + DT)
!     Note: ETPHOT calculates TSRF(3), which is surface temperature by
!     canopy zone.  1=sunlit leaves.  2=shaded leaves.  3= soil.  Should
!     we combine these variables?  At this time, only SRFTEMP is used
!     elsewhere. - chp 11/27/01
!
! ** Calculate temperature at top, middle and bottom of each layer **
!
! Calculate annual temperature at middle of top soil layer
        ZD = -DSMID(1)/DDa(1)
       STa(1) = TAV + ((TAMP/2.0)*COS(ALX+ZD))*EXP(ZD)
! Add attenuated weather + surface deviaiton from ann air temp
        ZD = -DSMID(1)/DDw(1)
        ST(1) = STa(1) + DT*EXP(ZD)
! Calculate annual temp and amp at bottom of top soil layer
        ZD = -DLAYR(1)/DDa(1)
        STBota(1) = TAV + ((TAMP/2.0)*COS(ALX+ZD))*EXP(ZD)
       AMPa(1) = (TAMP/2.0)*EXP(ZD)
! Add attenuated weather + surface deviaiton from ann air temp        
        ZD = -DLAYR(1)/DDw(1)
        STBoti(1) = STBota(1) + DT*EXP(ZD)
        AMPw(1) = DT*exp(ZD)
!
!       compute soil temp at middele and bottom of each layer
!        and amplitude
      
      IF (NLAYR .GT. 1) THEN
          DO L = 2,NLAYR
! Calculate annual temperature at middle of Lth soil layer
             ZD = -DSMID(L)/DDa(L)
              STa(L) = TAV + (AMPa(L-1)*COS(ALX+ZD))*EXP(ZD)
! Add attenuated weather + surface deviaiton from ann air temp
              ZD = -DSMID(L)/DDw(L)
              ST(L) = STa(L) + AMPw(L-1)*EXP(ZD)
! Calculate annual temp and amp at bottom of Lth soil layer
              ZD = -DLAYR(L)/DDa(L) 
              STBota(L) = TAV + (AMPa(L-1)*COS(ALX+ZD))*EXP(ZD)
              AMPa(L) = AMPa(L-1)*EXP(ZD)
! Add attenuated weather + surface deviaiton from ann air temp              
              ZD = -DLAYR(L)/DDw(L)
              STBoti(L) = STBota(L) + AMPw(L-1)*EXP(ZD)
              AMPw(L) = AMPw(L-1)*EXP(ZD)
          END DO
      END IF        
        
!-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE SOILTBK
C=======================================================================


!=======================================================================
! STEMP and SOILT Variable definitions - updated 2/15/2004
!=======================================================================
! ABD      Average bulk density for soil profile (g [soil] / cm3 [soil])
! AHeatCap Average soil heat capactity for whole profile (J/(m3 C))
! ALBEDO   Reflectance of soil-crop surface (fraction)
! ALBP     Albedo of plant canopy
! ALBS     Albedo of bare soil
! ALX      = (Day of year - hotest day)* PI/180 to convert deg to rad
! ASTCond  Average soil thermal conductivity for whole profile (W/(m C))
! ATOT     Sum of TMA array (last 5 days soil temperature) ( C)
! B        Exponential decay factor (Parton and Logan) (in subroutine
!            HTEMP)
! BD(L)    Bulk density, soil layer L (g [soil] / cm3 [soil])
! CLAY(L)  Clay (% by weight)
! CONTROL  Composite variable containing variables related to control
!            and/or timing of simulation.    See Appendix A.
! CP       Heat capaciy of dry air at sea level and 20 C (J/(kg C))
! CUMDPT   Cumulative depth of soil profile (mm)
! DampDa   Annual damping depth using average thermal conductivity
!            and heat capacity for whole profile (cm)
! DampDw   5-day damping depth for weather using average thermal
!            conductivity and heat capacity for whole profile (cm)
! DDa(L)   Annual damping depth for Lth layer (cm)
! DDw(L)   5-day damping depth for Lth layer (cm)
! Del      Deviation of daily surface temp (SRFT) from
!              daily avg air temp (TAVG) (C)
! DLAYR(L) Thickness of soil layer L (cm)
! DOY      Current day of simulation (d)
! DP
! DS(L)    Cumulative depth in soil layer L (cm)
! DSMID    Depth to midpoint of soil layer L (cm)
! DT       Deviation of surface temp (SRFTEMP) from
!             annual air temp for the day (TA) (C)
! DUL(L)   Volumetric soil water content at Drained Upper Limit in soil
!            layer L (cm3[water]/cm3[soil])
! EPSM     Emissivity of mulch
! EPSP     Emissivity of plant surface
! EPSS     Emissivity of bare soil surface
! ERRNUM   Error number for input
! EMWatt   Daily Average latent heat flux from mulch surface (W/m2)
! EPWatt   Daily average latent heat flux from plant canopy (W/m2)
! ESWatt   Daily average latent heat flux from soil surface (W/m2)      
! FILEIO   Filename for input file (e.g., IBSNAT35.INP)
! FOUND    Indicator that good data was read from file by subroutine FIND
!            (0 - End-of-file encountered, 1 - NAME was found)
! FX
! GaM      Mulch surface heat flux for annual wave (W/m2)
! GaS      Soil surface heat flux for annual wave per Novak (W/m2)
! HDAY
! ICWD     Initial water table depth (cm)
! ISWITCH  Composite variable containing switches which control flow of
!            execution for model.  The structure of the variable
!            (SwitchType) is defined in ModuleDefs.for.
! ISWWAT   Water simulation control switch (Y or N)
! Lamda    Latent heat of vaporization of water (J/kg)
! LINC     Line number of input file
! LL(L)    Volumetric soil water content in soil layer L at lower limit
!           (cm3 [water] / cm3 [soil])
! LNUM     Current line number of input file
! LUNIO    Logical unit number for FILEIO
! Mfrac    Fraction of surface that is mulch
! MSG      Text array containing information to be written to WARNING.OUT
!            file.
! MSGCOUNT Number of lines of message text to be sent to WARNING.OUT
 !
!   *** MULCH stuff from Cheryl ***
 !     Data construct for mulch layer
! Organic mulch variables are stored in a composite variable MULCH,
!      which contains: 
!     Data construct for mulch layer
!      TYPE MulchType
!        REAL MULCHMASS    !Mass of surface mulch layer (kg[dry mat.]/ha)
!        REAL MULCHALB     !Albedo of mulch layer
!        REAL MULCHCOVER   !Coverage of mulch layer (frac. of surface)
!        REAL MULCHTHICK   !Thickness of mulch layer (mm)
!        REAL MULCHEVAP    !Evaporation from mulch layer (mm/d)
!        REAL MULCHSAT     !Saturation water content of mulch (mm3/mm3)
!        REAL MULCHN       !N content of mulch layer (kg[N]/ha)
!        REAL MULCHP       !P content of mulch layer (kg[P]/ha)
!        REAL NEWMULCH     !Mass of new surface mulch (kg[dry mat.]/ha)
!        REAL NEWMULCHWAT  !Water content of new mulch ((mm3/mm3)
!        REAL MULCH_AM     !Area covered / dry weight of residue (ha/kg)
!        REAL MUL_EXTFAC   !Light extinction coef. for mulch layer
!        REAL MUL_WATFAC   !Saturation water content (mm[water] ha kg-1)
!      END TYPE MulchType

!Mulch mass, N, and P are initialized and modified by the two organic matter routines. Surface organic
!matter can be applied, tilled, or decomposed, so values change daily with the Century model.
!The CERES model does not decompose surface mulch.

!Subroutine MULCHLAYER determines the thickness and areal fractional coverage of the mulch layer.
!This is called daily by both organic matter routines.

!Subroutine MULCHWATER determines the water content of the mulch layer. Called by WatBal. Rainfall
!and surface irrigation wet the mulch layer preferentially, then the soil layers over the portion
!of soil covered by mulch.

!MULCHEVAP calculates the mulch evaporation, which is prioritized over soil evaporation.
Called by SPAM.

!I think that the effects of a surface mulch layer on soil temperature could be done completely
!in the soil temperature routine, using the daily mulch state variables, MulchMass, MulchAlb, 
!MulchCover, MulchThick, and MulchWat. Can it be handled as another soil layer? There is a potential
!for partial soil coverage which may be tricky to handle.

           
! NLAYR    Actual number of soil layers
! Omega    Temporary varialbe holding angular velocity (radians/s)
! PESW     Potential extractable soil water (= SW - LL) summed over root
!            depth (cm)
! Pfrac    Fraction of surface that is plant canopy
! RHO      Density of dry air at sea level and 20 C (kg/m3)
! RNMODE    Simulation run mode (I=Interactive, A=All treatments,
!             B=Batch mode, E=Sensitivity, D=Debug, N=Seasonal, Q=Sequence)
! RUN      Change in date between two observations for linear interpolation
! MSALB    Soil albedo with mulch and soil water effects (fraction)
! SAND(L)  Sand (% by weight)
! SECTION  Section name in input file
! Sfrac    Fraction of surface that is bare soil
! SILT(L)  Silt 
! SOILPROP Composite variable containing soil properties including bulk
!            density, drained upper limit, lower limit, pH, saturation
!            water content.  Structure defined in ModuleDefs.
! SolAvg   Average solar radiation for day (W/m2)
! SRAD     Solar radiation (MJ/m2-d)
! SRFT     Daily average of soil surface, plant surface, 
!             and mulch bottom temperatures
! SRFTEMP  Surface temperature = Average air temperature for last 5 days
!            plus average 5-day difference between soil-plant-mulch-bottom
!            temperature and daily average air temperature ( C).
! ST(L)    Soil temperature at middle of soil layer L ( C)
! STBZ     Stefan-Boltzmann constant 5.6697E-8 (W/(m2 K))
! SW(L)    Volumetric soil water content in layer L
!           (cm3 [water] / cm3 [soil])
! SWI(L)   Initial soil water content (cm3[water]/cm3[soil])
! TA       Daily normal temperature ( C)
! TAMP     Amplitude of temperature function used to calculate soil
!            temperatures ( C)
! TAV      Average annual soil temperature, used with TAMP to calculate
!            soil temperature. ( C)
! TAVG     Average daily temperature ( C)
! TBD      Sum of bulk density over soil profile
! TDL      Total water content of soil at drained upper limit (cm)
! TLL      Total soil water in the profile at the lower limit of
!            plant-extractable water (cm)
! TM0      Mulch surface temperature (C)
! TMA(I)   Array of previous 5 days of average soil temperatures. ( C)
! TMAX     Maximum daily temperature ( C)
! TMBot    Mulch bottom temperature (C)
! TP0      Plane surface temperature (C)
! TS0      Bare soil surface temperature (C)
! TSW      Total soil water in profile (cm)
! WC
! WW
! X1M      Net solar radiation for mulch (W/m2)
! X1P      Net solar radiation for crop canopy (W/m2)
! X1S      Net solar radiaiton for bare soil surface (W/m2)
! XHLAI    Leaf area index (cm2/cm2)
! XLAT     Latitude (deg.)
! Xsky     Long-wave average daily sky radiation following Prata (W/m2)
! YEAR     Year of current date of simulation
! YRDOY    Current day of simulation (YYYYDDD)
! ZD       Temporary variable for soil depth/damping depth
!=======================================================================
