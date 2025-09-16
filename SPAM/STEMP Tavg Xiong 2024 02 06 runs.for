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
     &      EOP, TRWUP, XHLAI,                        !Input
     &    SRFTEMP, ST)                                    !Output

C-----------------------------------------------------------------------
      USE ModuleDefs     !Definitions of constructed variable types,
                         ! which contain control information, soil
                         ! parameters, hourly weather data.
      IMPLICIT  NONE
      EXTERNAL YR_DOY,ERROR,FIND,SOILTBK,OPSTEMP
      SAVE

      CHARACTER*1  RNMODE, ISWWAT,MEEVP !, IDETL (CSVC ADD MEEVP)
      CHARACTER*6  SECTION
      CHARACTER*6, PARAMETER :: ERRKEY = "STEMP "
      CHARACTER*30 FILEIO

      INTEGER DOY, DYNAMIC, I, L, NLAYR
      INTEGER RUN, YRDOY, YEAR
      INTEGER ERRNUM, FOUND, LNUM, LUNIO

      REAL ABD, ALBEDO, ATOT, B, CUMDPT
      REAL DP, FX, HDAY, ICWD, PESW, MSALB, SRAD, SRFTEMP
      REAL TAMP, TAV, TAVG, TBD, TMAX, XLAT, WW
      REAL TDL, TLL, TSW, TA, DT
      REAL TMA(5)
      REAL EOP, TRWUP,XHLAI, VPD, TDEW, SWFAC, EP1,AVP,SVP
      REAL, DIMENSION(NL) :: BD, DLAYR, DS, DUL, LL, ST, SW, SWI, DSMID
      REAL DLI(NL),TEXSUM(NL)
      REAL CLAYV(NL),SILTV(NL),SANDV(NL),OMV(NL),TotSolid(NL),POR(NL)
      REAL CLAYC(NL),SILTC(NL),SANDC(NL),OMC(NL),OM(NL),Totmineral(NL)
      REAL ClayFrac(NL), SiltFrac(NL), SandFrac(NL), OMFrac(NL)
      REAL TcondDry(NL), TcondS(NL), TcondSat(NL), SWREL(NL),DSI(NL)
      REAL  PX, QX, RX, SX
      REAL HeatCap(NL),STBot(NL),AMP(NL),DSMIDV2(NL)
      REAL CLAY(NL),SILT(NL),SAND(NL),OC(NL),STCOND(NL)
      REAL DampD(NL),DampDd(NL),Omega, Del,STa(NL),STboti(NL),AMPi(NL)

!-----------------------------------------------------------------------
      TYPE (ControlType) CONTROL
      TYPE (SoilType) SOILPROP
      TYPE (SwitchType) ISWITCH

!     Check for output verbosity
!     IDETL  = ISWITCH % IDETL
!     IDETL = 'N', '0' (zero), suppress output
!     IF (INDEX('N0',IDETL) > 0) RETURN

!     Transfer values from constructed data types into local variables.
      DYNAMIC = CONTROL % DYNAMIC
      YRDOY   = CONTROL % YRDOY

      ISWWAT = ISWITCH % ISWWAT
      MEEVP  = ISWITCH % MEEVP  !SVC

      BD     = SOILPROP % BD
      DLAYR  = SOILPROP % DLAYR
      DS     = SOILPROP % DS
      DUL    = SOILPROP % DUL
      LL     = SOILPROP % LL
      NLAYR  = SOILPROP % NLAYR
      MSALB  = SOILPROP % MSALB
      CLAY   = SOILPROP % CLAY   ! clay (% by weight)
      SILT   = SOILPROP % SILT   ! silt (% by weight)
      SAND   = SOILPROP % SAND   ! sand (% by weight)
      OC     = SOILPROP % OC     ! organic carbon (g C/g soil)
      DSI    = SOILPROP % DS
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
          DSMID(L) = CUMDPT + DLAYR(L)* 5.0
          CUMDPT   = CUMDPT + DLAYR(L)*10.0
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
          TMA(I) = NINT(TAVG*10000.)/10000.   !chp
        END DO
        ATOT = TMA(1) * 5.0

        DO L = 1, NLAYR
          ST(L) = TAVG
        END DO

        DO I = 1, 8
          CALL SOILTBK (
     &        ALBEDO, B, CUMDPT, DOY, DP, HDAY, NLAYR,    !Input
     &        PESW, SRAD, TAMP, TAV, TAVG, TMAX, WW, DSMID,!Input
     &    SOILPROP,                                     !Input
     &        ATOT, TMA, SRFTEMP, ST,                     !Output
!            added by BAK 2023 11 29 for testing
     &    TA,DT,                                          !Output
     &    SWREL,TcondDry, TcondSat, STCOND,HeatCap,       !Output
     &    DampD,DampDd,STBot,AMP,                         !Output
     &    CLAYV,SILTV,SANDV,OMV,Del,STa,STboti,AMPi)      !Output
          END DO
      ENDIF

!     Print soil temperature data in STEMP.OUT
          IF (MEEVP .NE. 'Z')    !CSVC
     & CALL OPSTEMP(CONTROL, ISWITCH, DOY, SRFTEMP, ST, TAV, TAMP)

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
     &    ALBEDO, B, CUMDPT, DOY, DP, HDAY, NLAYR,        !Input
     &    PESW, SRAD, TAMP, TAV, TAVG, TMAX, WW, DSMID,   !Input
     &    SOILPROP,                                       !Input
     &    ATOT, TMA, SRFTEMP, ST,                         !Output
!            added by BAK 2023 11 29 for testing
     &    TA,DT,                                          !Output
     &    SWREL,TcondDry, TcondSat, STCOND,HeatCap,       !Output
     &    DampD,DampDd,STBot,AMP,                         !Output
     &    CLAYV,SILTV,SANDV,OMV,Del,STa,STboti,AMPi)      !Output
!***********************************************************************
!***********************************************************************
!     Output & Seasonal summary
!***********************************************************************
      ELSEIF (DYNAMIC .EQ. OUTPUT .OR. DYNAMIC .EQ. SEASEND) THEN
!-----------------------------------------------------------------------
          IF (MEEVP .NE. 'Z')    !CSVC
     & CALL OPSTEMP(CONTROL, ISWITCH, DOY, SRFTEMP, ST, TAV, TAMP)

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
     &    ALBEDO, B, CUMDPT, DOY, DP, HDAY, NLAYR,      !Input
     &    PESW, SRAD, TAMP, TAV, TAVG, TMAX, WW, DSMID, !Input
     &    SOILPROP,                                     !Input
     &    ATOT, TMA, SRFTEMP, ST,                       !Output
!          added by BAK 2023 11 29 for testing
     &    TA,DT,                          !Output
     &    SWREL,TcondDry, TcondSat, STCOND,HeatCap,     !Output
     &    DampD,DampDd,STBot,AMP,                       !Output
     &    CLAYV,SILTV,SANDV,OMV,Del,STa,STboti,AMPi)    !Output

!     ------------------------------------------------------------------
      USE ModuleDefs     !Definitions of constructed variable types,
                         ! which contain control information, soil
                         ! parameters, hourly weather data.
!     NL defined in ModuleDefs.for

      IMPLICIT  NONE
      SAVE

      INTEGER  K, L, DOY, NLAYR

      REAL ALBEDO, ALX, ATOT, B, CUMDPT, DD, DP, DT, FX
      REAL HDAY, PESW, SRAD, SRFTEMP, TA, TAMP, TAV, TAVG, TMAX
      REAL WC, WW, ZD
      REAL TMA(5)
      REAL DSMID(NL)
      REAL ST(NL)
      REAL SW(NL),DS(NL),DLI(NL)
      REAL CLAYV(NL),SILTV(NL),SANDV(NL),OMV(NL),TotSolid(NL),POR(NL)
      REAL CLAYC(NL),SILTC(NL),SANDC(NL),OMC(NL),OM(NL),Totmineral(NL)
      REAL ClayFrac(NL), SiltFrac(NL), SandFrac(NL), OMFrac(NL)
      REAL TcondDry(NL), TcondS(NL), TcondSat(NL), SWREL(NL)
      REAL  PX, QX, RX, SX
      REAL HeatCap(NL),STBot(NL),AMP(NL),DSMIDV2(NL)
      REAL BD(NL),CLAY(NL),SILT(NL),SAND(NL),OC(NL),STCOND(NL)
      REAL DampD(NL),DampDd(NL),Omega, Del,STa(NL),STboti(NL),AMPi(NL)

      TYPE (SoilType) SOILPROP

      BD     = SOILPROP % BD
      DS     = SOILPROP % DS
      NLAYR  = SOILPROP % NLAYR
      CLAY   = SOILPROP % CLAY   ! clay (% by weight)
      SILT   = SOILPROP % SILT   ! silt (% by weight)
      SAND   = SOILPROP % SAND   ! sand (% by weight)
      OC     = SOILPROP % OC     ! organic carbon (g C/g soil)
!-----------------------------------------------------------------------
      ALX    = (FLOAT(DOY) - HDAY) * 0.0174
      ATOT   = ATOT - TMA(5)

      DO K = 5, 2, -1
        TMA(K) = TMA(K-1)
      END DO
      
!    11/28/2023 BAK Inserting Xiong (2023) alternative method for
!    simulating soil thermal conductivity and ultimatly damping depth.
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
!         2.65 is the particle density of sand, silt, and clay
!         1.3 is the partcle density of organic matter (DeVries 1963, 1975)     
          CLAYV(L) = (CLAYC(L)/100.)*BD(L)/2.65     
          SILTV(L) = (SILTC(L)/100.)*BD(L)/2.65
          SANDV(L) = (SANDC(L)/100.)*BD(L)/2.65
          OMV(L)   = (OMC(L)/100)*BD(L)/1.3
          TotSolid(L) = CLAYV(L) + SILTV(L) + SANDV(L) + OMV(L)
          POR(L) = 1. - TotSolid(L)  ! Porosity
!
          ClayFrac(L) = CLAYV(L)/TotSolid(L)
          SiltFrac(L) = SILTV(L)/TotSolid(L)
          SandFrac(L) = SANDV(L)/TotSolid(L)
          OMFrac(L)   = OMV(L)/TotSolid(L)
!
! Calculate dry thermal conductivity (W m-1 C-1)
!     Equation was fitted to mineral soils
        TCondDry(L) = -0.6*POR(L) + 0.51
        PX = TCondDry(L)
!       
! Geometric mean thermal conductivity of solid materials for mineral soils
! from Xiong et al (2023)
        If(SandFrac(L) .GT. 0.2) THEN            
           TcondS(L) = (7.7**SandFrac(L))*
     &         (2.0**(SiltFrac(L) + ClayFrac(L)))*(0.25**OMFrac(L))
        ELSE
          TcondS(L) = (7.7**SandFrac(L))*
     &         (3.0**(SiltFrac(L) + ClayFrac(L)))*(0.25**OMFrac(L))
        END IF
! where 7.7 (W m-1 C-1) = thermal conductivity of quartz (sand) and
!       2.0 (W m-1 C-1) = thermal conductivity of other soil minerals
!               if sand > 0.2; 
!       3.0 otherwise
!       0.25 thermal conductivity of organic matter (DeVries 1963, 1975)
!
! Geometric mean thermal conductivity of soil solids and water at saturation
        TcondSat(L) = (TcondS(L)**(1. - POR(L)))*(0.594**POR(L))
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
     &                  + SW(L)*1.0E6
        
!     Calculate annual and daily Damping depths (cm)
        Omega = 2.0*3.14159/(365.0*24.0*3600.0)    ! radians/s
        DampD(L) = 100.*SQRT(2.*STCOND(L)/(HeatCap(L)*Omega)) ! annual
        Omega = 2.0*3.14159/(24.0*3600.0)    ! radians/s
        DampDd(L) = 100.*SQRT(2.*STCOND(L)/(HeatCap(L)*Omega)) ! daily 
!
       END DO       
!
!      Get rid of solar radiation stuff and just use TAVG      
!      TMA(1) = (1.0 - ALBEDO) * (TAVG + (TMAX - TAVG) *
!     &      SQRT(SRAD * 0.03)) + ALBEDO * TMA(1)
       TMA(1) = TAVG
       
!     Prevents differences between release & debug modes:
!       Keep only 4 decimals. chp 06/03/03
      TMA(1) = NINT(TMA(1)*10000.)/10000.  !chp
      ATOT = ATOT + TMA(1)

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
        WC = AMAX1(0.01, PESW) / (WW * CUMDPT) * 10.0
!     frac =              cm   / (    mm     ) * mm/cm
        !WC (ratio)
        !PESW (cm)
        !WW (dimensionless)
        !CUMDPT (mm)
!      END SELECT
!-----------------------------------------------------------------------

      FX = EXP(B * ((1.0 - WC) / (1.0 + WC))**2)

      DD = FX * DP                                  !DD in mm
!     JWJ, GH 12/9/2008
!     Checked damping depths against values from literature and
!       values are reasonable (after fix to WC equation).
!     Hillel, D. 2004. Introduction to Environmental Soil Physics.
!       Academic Press, San Diego, CA, USA.

      TA = TAV + TAMP * COS(ALX) / 2.0
      DT = ATOT / 5.0 - TA

      DO L = 1, NLAYR
        ZD    = -DSMID(L) / DD
        ST(L) = TAV + (TAMP / 2.0 * COS(ALX + ZD) + DT) * EXP(ZD)
        ST(L) = NINT(ST(L) * 1000.) / 1000.   !debug vs release fix
      END DO

!     Added: soil T for surface litter layer.
!     NB: this should be done by adding array element 0 to ST(L). Now
!     temporarily done differently.
      SRFTEMP = TAV + (TAMP / 2. * COS(ALX) + DT)
!     Note: ETPHOT calculates TSRF(3), which is surface temperature by
!     canopy zone.  1=sunlit leaves.  2=shaded leaves.  3= soil.  Should
!     we combine these variables?  At this time, only SRFTEMP is used
!     elsewhere. - chp 11/27/01

!-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE SOILTBK
C=======================================================================


!=======================================================================
! STEMP and SOILT Variable definitions - updated 2/15/2004
!=======================================================================
! ABD      Average bulk density for soil profile (g [soil] / cm3 [soil])
! ALBEDO   Reflectance of soil-crop surface (fraction)
! ALX
! ATOT     Sum of TMA array (last 5 days soil temperature) (°C)
! B        Exponential decay factor (Parton and Logan) (in subroutine
!            HTEMP)
! BD(L)    Bulk density, soil layer L (g [soil] / cm3 [soil])
! CONTROL  Composite variable containing variables related to control
!            and/or timing of simulation.    See Appendix A.
! CUMDPT   Cumulative depth of soil profile (mm)
! DD
! DLAYR(L) Thickness of soil layer L (cm)
! DOY      Current day of simulation (d)
! DP
! DS(L)    Cumulative depth in soil layer L (cm)
! DSMID    Depth to midpoint of soil layer L (cm)
! DT
! DUL(L)   Volumetric soil water content at Drained Upper Limit in soil
!            layer L (cm3[water]/cm3[soil])
! ERRNUM   Error number for input
! FILEIO   Filename for input file (e.g., IBSNAT35.INP)
! FOUND    Indicator that good data was read from file by subroutine FIND
!            (0 - End-of-file encountered, 1 - NAME was found)
! FX
! HDAY
! ICWD     Initial water table depth (cm)
! ISWITCH  Composite variable containing switches which control flow of
!            execution for model.  The structure of the variable
!            (SwitchType) is defined in ModuleDefs.for.
! ISWWAT   Water simulation control switch (Y or N)
! LINC     Line number of input file
! LL(L)    Volumetric soil water content in soil layer L at lower limit
!           (cm3 [water] / cm3 [soil])
! LNUM     Current line number of input file
! LUNIO    Logical unit number for FILEIO
! MSG      Text array containing information to be written to WARNING.OUT
!            file.
! MSGCOUNT Number of lines of message text to be sent to WARNING.OUT
! NLAYR    Actual number of soil layers
! PESW     Potential extractable soil water (= SW - LL) summed over root
!            depth (cm)
! RNMODE    Simulation run mode (I=Interactive, A=All treatments,
!             B=Batch mode, E=Sensitivity, D=Debug, N=Seasonal, Q=Sequence)
! RUN      Change in date between two observations for linear interpolation
! MSALB    Soil albedo with mulch and soil water effects (fraction)
! SECTION  Section name in input file
! SOILPROP Composite variable containing soil properties including bulk
!            density, drained upper limit, lower limit, pH, saturation
!            water content.  Structure defined in ModuleDefs.
! SRAD     Solar radiation (MJ/m2-d)
! SRFTEMP  Temperature of soil surface litter (°C)
! ST(L)    Soil temperature in soil layer L (°C)
! SW(L)    Volumetric soil water content in layer L
!           (cm3 [water] / cm3 [soil])
! SWI(L)   Initial soil water content (cm3[water]/cm3[soil])
! TA       Daily normal temperature (°C)
! TAMP     Amplitude of temperature function used to calculate soil
!            temperatures (°C)
! TAV      Average annual soil temperature, used with TAMP to calculate
!            soil temperature. (°C)
! TAVG     Average daily temperature (°C)
! TBD      Sum of bulk density over soil profile
! TDL      Total water content of soil at drained upper limit (cm)
! TLL      Total soil water in the profile at the lower limit of
!            plant-extractable water (cm)
! TMA(I)   Array of previous 5 days of average soil temperatures. (°C)
! TMAX     Maximum daily temperature (°C)
! TSW      Total soil water in profile (cm)
! WC
! WW
! XLAT     Latitude (deg.)
! YEAR     Year of current date of simulation
! YRDOY    Current day of simulation (YYYYDDD)
! ZD
!=======================================================================
