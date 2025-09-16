C=======================================================================
C  OPSTEMP, Subroutine, C.H.Porter
C  Generates output for daily soil temperature data
C-----------------------------------------------------------------------
C  REVISION HISTORY
C  11/01/2001 CHP Written
C  06/07/2002 GH  Modified for crop rotations
C-----------------------------------------------------------------------
C  Called from:   STEMP
C  Calls:         None
C=======================================================================
      SUBROUTINE OPSTEMPBK(CONTROL, ISWITCH, DOY,
     &   SRFTEMP, ST,TAV,TAMP,
!        added following outputs BAK 2023 11 29
     &   TMA,ATOT,TA,DT,
     &   DS,CLAY,SILT,SAND,OC,BD,SW,SWREL,POR,
     &    CLAYFrac,SILTFrac,SANDFrac,OMFrac,
     &   TcondDry, TcondSat, ASTCOND,AHeatCap,DampDa,DampDw,
     &   Del,STa)
!-----------------------------------------------------------------------
      USE ModuleDefs
      USE ModuleData
!     VSH
      USE CsvOutput
      USE Linklist
      IMPLICIT NONE
      EXTERNAL GETLUN, HEADER, YR_DOY
      SAVE
!-----------------------------------------------------------------------
      CHARACTER*1  RNMODE
      CHARACTER*12 OUTT

      INTEGER DAS, DOY, DYNAMIC, ERRNUM, FROP, L, N_LYR
      INTEGER MaxN_LYR
      INTEGER NOUTDT, RUN, YEAR, YRDOY, REPNO
      REAL ST(NL), SRFTEMP, TAV, TAMP
!        added following outputs BAK 2023 11 29
      REAL  DS(NL),CLAY(NL),SILT(NL),SAND(NL),OC(NL),BD(NL),SW(NL)
      REAL  SWREL(NL),POR(NL)
      REAL  TcondDry(NL),TcondSat(NL),STCOND(NL),HeatCap(NL)
      REAL  DampDa,DampDw,STBot(NL),AMP(NL)
      REAL  CLAYV(NL),SILTV(NL),SANDV(NL),OMV(NL)
      REAL  ClayFrac(NL),SiltFrac(NL),SandFrac(NL),OMFrac(NL)
      REAL  TMA(5),ATOT,TA,DT,Del,STa(NL),STboti(NL),AMPi(NL)
      REAL  ASTCond,AHeatCap
      LOGICAL FEXIST, DOPRINT

!-----------------------------------------------------------------------
!     The variable "CONTROL" is of constructed type "ControlType" as
!     defined in ModuleDefs.for, and contains the following variables.
!     The components are copied into local variables for use here.
!-----------------------------------------------------------------------
      TYPE (ControlType) CONTROL
      TYPE (SwitchType)  ISWITCH
      TYPE (SoilType)    SOILPROP

      IF (INDEX('N0',ISWITCH % IDETL) > 0) RETURN

      DAS     = CONTROL % DAS
      DYNAMIC = CONTROL % DYNAMIC
      FROP    = CONTROL % FROP
      YRDOY   = CONTROL % YRDOY

      FMOPT   = ISWITCH % FMOPT   ! VSH
!***********************************************************************
!***********************************************************************
!     Seasonal initialization - run once per season
!***********************************************************************
!      ELSEIF (DYNAMIC .EQ. SEASINIT) THEN
      IF (DYNAMIC .EQ. SEASINIT) THEN
!-----------------------------------------------------------------------
      RNMODE  = CONTROL % RNMODE
      REPNO   = CONTROL % REPNO
      RUN     = CONTROL % RUN

      IF (FMOPT == 'A' .OR. FMOPT == ' ') THEN   ! VSH
        CALL GETLUN('OUTT',NOUTDT)
!       Open the output files
        OUTT = 'SoilTemp.OUT'
        INQUIRE (FILE = OUTT, EXIST = FEXIST)
        IF (FEXIST) THEN
          OPEN (UNIT=NOUTDT, FILE=OUTT, STATUS='OLD',
     &      IOSTAT = ERRNUM, POSITION='APPEND')
          !IF (RNMODE .NE. 'Q') THEN
          !ENDIF
        ELSE
          OPEN (UNIT=NOUTDT, FILE=OUTT, STATUS='NEW',
     &      IOSTAT = ERRNUM)
 !        Write headers info to daily output file
          WRITE(NOUTDT,'("*SOIL TEMPERATURE OUTPUT FILE (DAILY)")')
        ENDIF
      END IF   ! VSH
C-----------------------------------------------------------------------
C     Variable heading for SoilTemp.OUT
C-----------------------------------------------------------------------
      IF (RNMODE .NE. 'Q' .OR. RUN .EQ. 1) THEN

        IF (FMOPT == 'A' .OR. FMOPT == ' ') THEN   ! VSH
          !For first run of a sequenced run, use replicate
          ! number instead of run number in header.
          IF (RNMODE .EQ. 'Q') THEN
            CALL HEADER(SEASINIT, NOUTDT, REPNO)
          ELSE
            CALL HEADER(SEASINIT, NOUTDT, RUN)
          ENDIF
        END IF   ! VSH

        CALL GET(SOILPROP)
        N_LYR = MIN(10, MAX(4,SOILPROP%NLAYR))

        IF (FMOPT == 'A' .OR. FMOPT == ' ') THEN   ! VSH
          WRITE (NOUTDT, '("! TAV  =",F8.1,/,"! TAMP =",F8.1)') TAV,TAMP
          WRITE (NOUTDT,
     &      '("!",T17,"Temperature (oC) by soil depth (cm):",
     &      /,"!",T17,"Surface",10A8)')(SoilProp%LayerText(L),L=1,N_LYR)
          
!       To output more data to right of soil temps,
!          Set maximum number of layers to 15
          MaxN_LYR = 15
          IF (N_LYR <= MaxN_LYR) THEN
              If(N_LYR < MaxN_LYR) THEN
              DO L = N_LYR+1, MaxN_LYR
                  ST(L) = 99.9
              END DO
              END IF
            WRITE (NOUTDT,120) ("TS",L,"D",L=1,MaxN_LYR),
     &                         (" DS",L,L=1,2),
     &                         ("CLA",L,L=1,2),
     &                         ("SIL",L,L=1,2),
     &                         ("SAN",L,L=1,2),
     &                         (" OC",L,L=1,2),
     &                         (" BD",L,L=1,2),
     &                         (" SW",L,L=1,2),
     &                         ("SWR",L,L=1,2),
     &                         ("POR",L,L=1,2),
     &                         ("CLF",L,L=1,2),
     &                         ("SIF",L,L=1,2),
     &                         ("SAF",L,L=1,2),
     &                         ("OMF",L,L=1,2),        
     &                        ("TCD",L,L=1,2),
     &                        ("TCS",L,L=1,2),    
     &                        ("STa",L,L=1,2)
   
  120    FORMAT('@YEAR DOY   DAS    TS0D',9("    ",A2,I1,A1),
     &                                    6("   ",A2,I2,A1),
     &         "     TMA     ATO      TA      DT",
     &                                32("    ",A3,I1),
     &         "     DDa     DDw     AST     AHC")
!     &    '    TS1D    TS2D    TS3D    TS4D    TS5D',
!     &    '    TS6D    TS7D    TS8D    TS9D    TS10')
          ELSE
            WRITE (NOUTDT,122) ("ST",L,"D",L=1,9), "    TS10"
  122       FORMAT('@YEAR DOY   DAS    TS0D',9("  ",A2,I1,A1),A8)
          ENDIF
        END IF   ! VSH
      ENDIF

      ENDIF !DYNAMIC

!***********************************************************************
!***********************************************************************
!     Daily Output
!***********************************************************************
      DOPRINT = .FALSE.
      SELECT CASE (DYNAMIC)
!      CASE (SEASINIT)
!        DOPRINT = .TRUE.
      CASE (OUTPUT)
        IF (MOD(DAS, FROP) == 0) THEN
          DOPRINT = .TRUE.
        ENDIF
      CASE (SEASEND)
        IF (MOD(DAS, FROP) /= 0) THEN
          DOPRINT = .TRUE.
        ENDIF
      END SELECT
      IF (DAS == 1) DOPRINT = .TRUE.

      IF (DOPRINT) THEN
        CALL YR_DOY(YRDOY, YEAR, DOY)
        IF (FMOPT == 'A' .OR. FMOPT == ' ') THEN   ! VSH
!         Generate output for file SoilTemp.OUT
          WRITE (NOUTDT,300) YEAR,DOY,DAS,SRFTEMP,(ST(L),L=1,MaxN_LYR),
     &        TMA(1),ATOT,TA,DT,
     &        DS(1),DS(2),CLAY(1),CLAY(2),
     &        SILT(1),SILT(2), 
     &        SAND(1),SAND(2),OC(1),OC(2),
     &        BD(1),BD(2),
     &        SW(1),SW(2),
     &        SWREL(1),SWREL(2),POR(1),POR(2),
     &        CLAYFrac(1),ClayFrac(2),SILTFrac(1),SILTFrac(2),
     &        SANDFrac(1),SandFrac(2),OMFrac(1),OMFrac(2),     
     &        TCondDry(1),TCondDry(2),
     &        TcondSat(1),TCondSat(2),     
     &        STa(1),STa(2),
     &        DampDa,DampDw,ASTCOND,AHeatCap
  300     FORMAT(1X,I4,1X,I3.3,1X,I5,16F8.1,
     &           39F8.2, E15.4)
             ! 10F8.2,21F8.3,4F8.2,2F12.0,E15.4,7F8.2)
        END IF   ! VSH

!       VSH CSV output corresponding to SoilTEMP.OUT
        IF (FMOPT == 'C') THEN ! VSH
          CALL CsvOutTemp_crgro(EXPNAME,CONTROL%RUN, CONTROL%TRTNUM,
     &CONTROL%ROTNUM,CONTROL%REPNO, YEAR, DOY, DAS, SRFTEMP,
     &N_LYR, ST, vCsvlineTemp, vpCsvlineTemp, vlngthTemp)

          CALL LinklstTemp(vCsvlineTemp)
        ENDIF

      ENDIF

!***********************************************************************
!***********************************************************************
!     SEASEND
!***********************************************************************
!      IF (DYNAMIC .EQ. SEASEND) THEN
      IF ((DYNAMIC == SEASEND)
     & .AND. (FMOPT == 'A'.OR.FMOPT == ' ')) THEN ! VSH
!-----------------------------------------------------------------------
        CLOSE (NOUTDT)

!***********************************************************************
!***********************************************************************
!     END OF DYNAMIC IF CONSTRUCT
!***********************************************************************
      ENDIF
!***********************************************************************
      RETURN
      END SUBROUTINE OPSTEMPBK
!***********************************************************************
