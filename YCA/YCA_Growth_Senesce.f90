!***************************************************************************************************************************
! This is the code from the section (DYNAMIC == RATE) lines 4538 - 4642 of the original CSCAS code. The names of the 
! dummy arguments are the same as in the original CSCAS code and the call statement and are declared here. The variables 
! that are not arguments are declared in module YCA_First_Trans_m. Unless identified as by MF, all comments are those of 
! the original CSCAS.FOR code.
!
! Subroutine YCA_Growth_Senesce calculates senescence and remobilization.
!***************************************************************************************************************************
    SUBROUTINE YCA_Growth_Senesce ( &
        ISWNIT      , ISWWAT,     BRSTAGE      & 
        )
    
        USE YCA_First_Trans_m
        USE YCA_LeafControl
    
        IMPLICIT NONE
        
        CHARACTER(LEN=1) ISWNIT      , ISWWAT
        REAL BRSTAGE
    
        !-----------------------------------------------------------------------
        !           Calculate senescence of leaves,stems,etc..
        !-----------------------------------------------------------------------

        ! LAH Notes from original cassava model. May need to take
        ! into account procedure to calculate leaf senescence. 
        ! Leaves are assumed to have a variety-specific maximum 
        ! life, which can be influenced by temperature and shading
        ! by leaves above. Water stress is assumed not to have any
        ! effect on leaf life (Cock, pers. comm.). However, on 
        ! release of stress leaves drop off and are replaced by a 
        ! flush of new leaves. This is not yet built into the 
        ! model.

        PLASP = 0.0
        PLASI = 0.0
        PLASL = 0.0
        PLASS = 0.0

        ! Leaf senescence - phyllochron or real time driven
        LAPSTMP = 0.0

        DO BR = 0, BRSTAGE                                                                                        !LPM 21MAR15
            DO LF = 1, LNUMSIMSTG(BR)         
                IF (node(BR,LF)%LAGETT+TTLFLIFE*EMRGFR <= LLIFGTT+LLIFATT) EXIT                                                     !EQN 371 LPM28MAR15 Deleted LLIFGTT
                IF (node(BR,LF)%LATL3T-node(BR,LF)%LAPS > 0.0) THEN
                    LAPSTMP = AMIN1((node(BR,LF)%LATL3T - node(BR,LF)%LAPS),(node(BR,LF)%LATL3T*(AMIN1((node(BR,LF)%LAGETT+(TTLFLIFE*EMRGFR)-(LLIFGTT+LLIFATT)), &         !EQN 372
                        (TTLFLIFE*EMRGFR))/LLIFSTT)))
                    node(BR,LF)%LAPS = node(BR,LF)%LAPS + LAPSTMP
                    PLASP = PLASP + LAPSTMP                                                                                !EQN 370
                ENDIF
            ENDDO
        ENDDO

        ! Leaf senescence - injury        ! LAH  To add later?
        !PLASI = PLA*(LSENI/100.0)*DU/STDAY  ! May need injury loss

        ! Leaf senescence - water or N stress
        ! LAH Need to accelerated senescence rather than lose leaf
        PLASW = 0.0
        PLASN = 0.0
        IF (ISWWAT /= 'N') THEN
            IF (PLA-SENLA > 0.0.AND.WUPR < WFSU) PLASW = AMAX1(0.0,AMIN1((PLA-SENLA)-PLAS,(PLA-SENLA)*LLOSA))        !EQN 373
        ENDIF
        IF (ISWNIT /= 'N') THEN
            LNCSEN = LNCM + NFSU * (LNCX-LNCM)                                                                         !EQN 374
            IF (PLA-SENLA > 0.0.AND.LANC < LNCSEN) PLASN = AMAX1(0.0,AMIN1((PLA-SENLA)-PLAS,(PLA-SENLA)*LLOSA))
        ENDIF
        ! LAH TMP
        PLASW = 0.0
        PLASN = 0.0
        PLASS = PLASW + PLASN    ! Loss because of stress
              
        ! LPM As a NOTE for further code: if PLASW and PLASN want to be added to the code, should be calculated by cohort
        
        !-----------------------------------------------------------------------
        !        LAI by Cohort
        !-----------------------------------------------------------------------
        node(BR,LF)%LAIByCohort=0.0                               ! DA re-initializing LAIByCohort
        LAI=0.0                                              ! DA re-initializing LAI
        
        DO Bcount=0,BRSTAGE
            BR= BRSTAGE - Bcount                                                        ! DA 28OCT2016 to run the loop to the higher branch to the lowest
            DO Lcount=0,LNUMSIMSTG(BR)-1
                LF=LNUMSIMSTG(BR)-Lcount                                                ! DA to run the loop to the higher leaf to the lowest
                IF (isLeafAlive(node(BR,LF))) THEN                      ! DA if leave is alive
                    IF(BR == INT(BRSTAGE) .AND. LF == INT(LNUMSIMSTG(INT(BRSTAGE)))) THEN                   ! DA if the very first leaf of the top of the highest branch
                        node(BR,LF)%LAIByCohort = (node(BR,LF)%LATL3T-node(BR,LF)%LAPS)*PLTPOP*0.0001                      ! DA calculate LAI of the leaf
                    ELSE                                                                                    ! DA if further leaf
                        node(BR,LF)%LAIByCohort= LAI + (node(BR,LF)%LATL3T-node(BR,LF)%LAPS)*PLTPOP*0.0001                ! DA the LAI calculation is accumulative from the previous cohort LAI
                    ENDIF
                    LAI = node(BR,LF)%LAIByCohort                                                                  ! DA updating LAI
                    
                ENDIF
            ENDDO
        ENDDO
        
        ! Leaf senescence - low light at base of canopy
        ! NB. Just senesces any leaf below critical light fr 
        PLASL = 0.0
        !IF (LAI > LAIXX) THEN
        !    PLASL = (LAI-LAIXX) / (PLTPOP*0.0001)
        !    ! LAH Eliminated! Replaced by accelerated senescence
        !    PLASL = 0.0
        !ENDIF
            
        ! Leaf senescence - overall
        PLAS =  PLASP + PLASI + PLASS + PLASL                                                                          !EQN 369
        ! Overall check to restrict senescence to what available
        PLAS = AMAX1(0.0,AMIN1(PLAS,PLA-SENLA))

        !-----------------------------------------------------------------------
        !           Calculate C and N made available through senescence
        !-----------------------------------------------------------------------

        SENLFG = 0.0
        SENLFGRS = 0.0
        SENNLFG = 0.0
        SENNLFGRS = 0.0
        IF (PLA-SENLA > 0.0) THEN
        ! LAH New algorithms 03/04/13
        SENLFG = AMIN1(LFWT*LWLOS,(AMAX1(0.0,(LFWT*(PLAS/(PLA-SENLA))*LWLOS))))                                        !EQN 375
        SENLFGRS = AMIN1(LFWT*(1.0-LWLOS),(AMAX1(0.0,(LFWT*(PLAS/(PLA-SENLA))*(1.0-LWLOS)))))                          !EQN 376
        ENDIF
  
        IF (ISWNIT /= 'N') THEN
            ! NB. N loss has a big effect if low N
            ! Assumes that all reserve N in leaves
            IF (LFWT > 0.0) LANCRS = (LEAFN+RSN) / LFWT                                                               !EQN 377
            SENNLFG = AMIN1(LEAFN,(SENLFG+SENLFGRS)*LNCM)                                                              !EQN 378
            SENNLFGRS = AMIN1(LEAFN-SENNLFG,(SENLFG+SENLFGRS)*(LANC-LNCM))                                             !EQN 379
        ELSE
            SENNLFG = 0.0
            SENNLFGRS = 0.0
        ENDIF

        !-----------------------------------------------------------------------
        !           Calculate overall senescence loss from tops
        !-----------------------------------------------------------------------

        SENFR = 1.0
        SENTOPLITTERG = 0.0
        SENTOPLITTERG = SENLFG*SENFR                                                                                   !EQN 380
        
    END SUBROUTINE YCA_Growth_Senesce