!**********************************************************************************************************************
! This is the code from the section (DYNAMIC.EQ.RUNINIT) ! Initialization, lines 1827 - 2336 of the original CSCAS code.
! The names of the dummy arguments are the same as in the original CSCAS code and the call statement and are declared 
! here. The variables that are not arguments are declared in module CS_First_Trans_m. Unless identified as by MF, all 
! comments are those of the original CSCAS.FOR code.
!
! Subroutine CS_SeasInit_VarInit initializes state and rate variables.
!**********************************************************************************************************************
    
    SUBROUTINE CS_SeasInit_VarInit( &
        BRSTAGE     , CAID        , CANHT       , DEWDUR      , LAIL        , LAILA       , NFP         , PARIP       , &
        PARIPA      , RESCALG     , RESLGALG    , RESNALG     , RLV         , SENCALG     , SENLALG     , SENNALG     , &
        STGYEARDOY  , TRWUP       , UH2O        , UNH4        , UNO3         &
        )  
        
        USE ModuleDefs
        USE CS_First_Trans_m

        IMPLICIT     NONE
        
        INTEGER STGYEARDOY(0:19)            
        
        REAL    BRSTAGE     , CAID        , CANHT       , DEWDUR      , LAIL(30)    , LAILA(30)   , NFP         , PARIP       
        REAL    PARIPA      , RESCALG(0:NL)             , RESLGALG(0:NL)            , RESNALG(0:NL)             , RLV(NL)     
        REAL    SENCALG(0:NL)             , SENLALG(0:NL)             , SENNALG(0:NL)             , TRWUP       , UH2O(NL)    
        REAL    UNH4(NL)    , UNO3(NL)    
        
        
        !-----------------------------------------------------------------------
        !       Initialize both state and rate variables                   
        !-----------------------------------------------------------------------

        aflf = 0.0
        amtnit = 0.0
        andem = 0.0
        brnumst = 1.0
        caid = 0.0
        canht = 0.0
        canhtg = 0.0
        carboadj = 0.0
        carbobeg = 0.0
        carbobegi = 0.0
        carbobegm = 0.0
        carbobegr = 0.0
        carboc = 0.0
        carboend = 0.0
        carbor = 0.0
        carbot = 0.0
        cdays  = 0
        cflfail = 'N'
        cfllflife = '-'
        cflharmsg = 'N'
        cflsdrsmsg = 'N'
        !crrswt = 0.0
        crwad = 0.0
        crwt = 0.0
        crwtp = 0.0 !LPM 23MAY2015 Added to keep the potential planting stick weight
        cnad = 0.0
        cnadprev = 0.0
        cnadstg = 0.0
        cnam = -99.0
        cnamm = -99.0
        co2cc = 0.0
        co2fp = 1.0
        co2intppm = 0.0
        co2intppmp = 0.0
        co2max = -99.0
        co2pav = -99.0
        co2pc = 0.0
        cumdu = 0.0
        cwad = 0.0
        cwadprev = 0.0
        cwadstg = 0.0
        cwahc = 0.0
        cwahcm = -99.0
        cwam = -99.0
        cwamm = -99.0
        dae = -99
        dap = -99
        dawwp = 0.0 !LPM 06MAR2016 DAWWP added to save Development Age (with stress)
        daylcc = 0.0
        daylpav = -99.0
        daylpc = 0.0
        daylst = 0.0
        daysum = 0.0
        sentoplitter = 0.0
        dewdur = -99.0
        df = 1.0
        dfout = 1.0
        dalf = 0
        dglf = 0
        dslf = 0
        drainc = 0.0
        !dstage = 0.0 !LPM 05JUN2015 DSTAGE is not used
        du = 0.0
        duneed = 0.0
        dynamicprev = 99999
        edap = -99
        edapfr = 0.0
        edapm = -99
        edayfr = 0.0
        emrgfr = 0.0
        emrgfrprev = 0.0
        emflag = 'N'
        eoc = 0.0
        eoebud = 0.0
        eoebudc = 0.0
        eoebudcrp = 0.0
        eoebudcrpc = 0.0
        eoebudcrpco2 = 0.0
        eoebudcrpco2c = 0.0
        eoebudcrpco2h2o = 0.0
        eoebudcrpco2h2oc = 0.0
        eompen = 0.0
        eompenc = 0.0
        eompcrp = 0.0
        eompcrpc = 0.0
        eompcrpco2 = 0.0
        eompcrpco2c = 0.0
        eompcrpco2h2o = 0.0
        eompcrpco2h2oc = 0.0
        eopen = 0.0
        eopenc = 0.0
        eopt = 0.0
        eoptc = 0.0
        epcc   = 0.0
        epsratio = 0.0
        established = 'n'
        etcc   = 0.0
        eyeardoy = -99
        fappline = ' '
        fappnum = 0
        fernitprev = 0.0
        fldap = 0
        fln = 0.0
        gdap = -99
        gdap = -99
        gdapm = -99
        gdayfr = 0.0
        gedayse = 0.0
        gedaysg = 0.0
        germfr = -99.0
        gestage = 0.0
        gestageprev = 0.0
        geucum = 0.0
        grocr = 0.0
        grocradj = 0.0
        grolf = 0.0
        grolfadj = 0.0
        grors = 0.0
        !grosr = 0.0 !LPM 05JUN2105 GROSR or basic growth of storage roots will not be used
        grost = 0.0
        grostadj = 0.0
        grostcr = 0.0
        gyeardoy = -99
        hamt = 0.0
        hbpc = -99.0
        hbpcf = -99.0
        hiad = 0.0
        hiam = -99.0
        hiamm = -99.0
        hind = 0.0
        hinm = -99.0
        hinmm = -99.0
        hnad = 0.0
        hnam = -99.0
        hnamm = -99.0
        hnpcm = -99.0
        hnpcmm = -99.0
        hnumam = -99.0
        hnumamm = -99.0
        hnumber = -99
        hnumgm = -99.0
        hnumgmm = -99.0
        hop = ' '
        hpc = -99.0
        hpcf = -99.0
        hstage = 0.0
        hwam = -99.0 
        hwamm = -99.0
        hwum = -99.0
        hwummchar = ' -99.0'
        hyrdoy = -99
        hyeardoy = -99
        idetgnum = 0
        irramtc = 0.0
        lagett = 0.0
        !lagep = 0.0
        lagl = 0.0
        lagl3 = 0.0 !LPM 15NOV15 added to save leaf area growing by cohort (considering assimilates restriction)
        lagl3t = 0.0 !LPM 15NOV15 added to save leaf area growing by cohort (considering assimilates restriction)
        laglt = 0.0 !LPM 15NOV15 added to save leaf area growing by cohort
        lai = 0.0
        laiprev = 0.0
        lail = 0.0
        laila = 0.0
        laistg = 0.0
        laix = 0.0
        laixm = -99.0
        lanc = 0.0
        lap = 0.0
        laphc = 0.0
        lapp = 0.0
        laps = 0.0
        latl = 0.0
        latl2 = 0.0
        latl2t = 0.0 !LPM 15NOV15 added to save leaf area by cohort
        latl3 = 0.0 
        latl3t = 0.0 !LPM 15NOV15 added to save leaf area by cohort (considering assimilates restriction)
        latl4 = 0.0
        !lcnum = 0 !LPM 28MAR15 Non necessary variables
        !lcoa = 0.0 !LPM 28MAR15 These variables are not necessary
        !lcoas = 0.0
        leafn = 0.0
        lfwt = 0.0
        lfwtm = 0.0
        llnad = 0.0
        !llrswad = 0.0
        !llrswt = 0.0 !LPM 21MAY2015 The reserves distribution will not be included, it needs to be reviewed
        llwad = 0.0
        lncr = 0.0
        lncx = 0.0
        lndem = 0.0
        lnphc = 0.0
        lnum = 0.0
        lnumsimstg = 0.0   !LPM 09AGO2015 To initialize the variable lnumsimstg
        lnumsimtostg = 0.0 !LPM 09AGO2015 To initialize the variable lnumsimtostg
        !lnumend = 0.0
        lnumg = 0.0
        lnumprev = 0.0
        lnumsg = 0
        lnumsm = -99.0
        lnumsmm = -99.0
        lnumstg = 0.0
        lnuse = 0.0
        lseed = -99
        lpeai = 0.0
        !lperswad = 0.0 !LPM 21MAY2015 The reserves distribution will not be included, it needs to be reviewed
        lperswt = 0.0
        lpewad = 0.0
        lstage = 0.0
        lwphc = 0.0
        mdap = -99
        mdat = -99
        mdayfr = -99
        mdoy = -99
        nfg = 1.0
        nfgcc = 0.0
        nfgcc = 0.0
        nfgpav = 1.0
        nfgpc = 0.0
        nflf = 1.0
        nflf2 = 0.0
        nflfp = 1.0
        nfp = 1.0
        nfpcav = 1.0
        nfpcc = 0.0
        nfppav = 1.0
        nfppc = 0.0
        NODEWTGB = 0.0      !LPM 11APR15 New variables of node growth
        NODEWT = 0.0       !LPM 11APR15 New variables of node growth
        nsdays = 0
        nuf = 1.0
        nupac = 0.0
        nupad = 0.0
        nupap = 0.0
        nupapcsm = 0.0
        nupapcsm1 = 0.0
        nupapcrp = 0.0
        nupc = 0.0
        nupd = 0.0
        nupratio = 0.0
        pari = 0.0
        pari1 = 0.0
        parip = -99.0
        paripa = -99.0
        pariue = 0.0
        parmjc = 0.0
        parmjic = 0.0
        paru = 0.0
        pdadj = -99.0
        pdays  = 0
        photqr = 0.0
        pla = 0.0
        !plags2 = 0.0     !LPM 23MAR15 non necessary PLAGSB2 considers all the branches and shoots
        plagsb2 = 0.0
        plagsb3 = 0.0
        plagsb4 = 0.0
        plas = 0.0
        plasi = 0.0
        plasl = 0.0
        plasp = 0.0
        plass = 0.0
        plax = 0.0
        pltpop = 0.0
        plyear = -99
        plyeardoy = 9999999
        psdap = -99
        psdapm = -99
        psdat = -99
        psdayfr = 0.0
        !ptf = 0.0 !LPM 19MAY2015 PTF is not considered in the model 
        rainc = 0.0
        raincc = 0.0
        rainpav = -99.0
        rainpc = 0.0
        ranc = 0.0
        rescal = 0.0
        rescalg = 0.0
        reslgal = 0.0
        reslgalg = 0.0
        resnal = 0.0
        resnalg = 0.0
        respc = 0.0
        resprc = 0.0
        resptc = 0.0
        reswal = 0.0
        reswalg = 0.0
        rlf = 0.0
        rlfc = 0.0
        rlv = 0.0
        rnad = 0.0
        rnam = -99.0
        rnamm = -99.0
        rncr = 0.0
        rndem = 0.0
        rnuse = 0.0
        rootn = 0.0
        rootns = 0.0
        rscd = 0.0
        rscm = 0.0
        rscx = 0.0
        rsfp = 1.0
        rsn = 0.0
        rsnph = 0.0
        rsnphc = 0.0
        rsnad = 0.0
        rsnused = 0.0
        brstage = 0.0
        brstageprev = 0.0
        rswad = 0.0
        rswam = -99.0
        rswamm = -99.0
        rswph = 0.0
        rswphc = 0.0
        rswt = 0.0
        rswtm = 0.0
        rswtx = 0.0
        rtdep = 0.0
        rtdepg = 0.0
        rtnsl = 0.0
        rtresp = 0.0
        rtrespadj = 0.0
        rtslxdate = -99 
        rtwt = 0.0
        rtwtal = 0.0
        rtwtg = 0.0
        rtwtgadj = 0.0
        rtwtgl = 0.0
        rtwtl = 0.0
        rtwtm = 0.0
        rtwtsl = 0.0
        runoffc = 0.0
        rwad = 0.0
        rwam = -99.0
        rwamm = -99.0
        said = 0.0
        sanc = 0.0
        sancout = 0.0
        sdnad = 0.0
        sdnc = 0.0
        sdwad = 0.0
        sdwam = -99.0
        seeduse = 0.0
        seeduser = 0.0
        seeduset = 0.0
        sencags = 0.0
        sencalg = 0.0
        sencas = 0.0
        sencl = 0.0
        sencs = 0.0
        senla = 0.0
        senlalitter = 0.0
        senlags = 0.0
        senlalg = 0.0
        senlas = 0.0
        senlfg = 0.0
        senlfgrs = 0.0
        senll = 0.0
        senls = 0.0
        sennags = 0.0
        sennal = 0.0
        sennalg = 0.0
        sennas = 0.0
        sennatc = -99.0
        sennatcm = -99.0
        sennl = 0.0
        sennlfg = 0.0
        sennlfgrs = 0.0
        senns = 0.0
        senrtg = 0.0
        SENTOPLITTERG = 0.0
        senwacm = -99.0
        senwacmm = -99.0
        senwags = 0.0
        senwal = 0.0
        senwalg = 0.0
        senwl = 0.0
        senroot = 0.0
        senroota = 0.0
        shdat = 0
        shla = 0.0
        shlag2 = 0.0
        shlag2b = 0.0 !LPM 23MAR15 add new variable
        shlas = 0.0
        shnum = 0.0
        shnumad = 0.0
        shnuml = 1.0
        shrtd = 0.0
        shrtm = 0.0
        sla = -99.0
        snad = 0.0
        sncr = 0.0
        sno3profile = 0.0
        sno3profile = 0.0
        sno3rootzone = 0.0
        snh4rootzone = 0.0
        snph = 0.0
        snphc = 0.0
        snuse = 0.0
        sradc = 0.0
        sradcav  = -99.0
        sradcc = 0.0
        sradd = 0.0
        sradpav  = -99.0
        sradpc = 0.0
        sradprev = 0.0
        srnam = -99.0
        sranc = 0.0
        srndem = 0.0
        srnoad = 0.0
        srnoam = 0.0
        srnogm = 0.0
        srnopd = 0.0
        srnuse = 0.0
        srootn = 0.0
        srwt = 0.0
        srwtgrs = 0.0
        srwud = 0.0
        srwum = 0.0
        srwum = 0.0
        stai = 0.0
        staig = 0.0
        stais = 0.0
        stemn = 0.0
        stemnn = 0.0!LPM 23MAY2015 Added to keep nitrogen content by node
        stgedat = 0
        stgyeardoy = 9999999
        !strswt = 0.0
        stwad = 0.0
        stwt = 0.0
        stwtp = 0.0 !LPM 23MAY2015 Added to keep the potential stem weight
        stwtm = 0.0
        swphc = 0.0
        tcan = 0.0
        tdifav = -99.0
        tdifnum = 0
        tdifsum = 0.0
        tfd = 0.0
        tfg = 1.0
        tfglf = 0.0
        tfdlf = 0.0
        tfp = 1.0
        tlchc = 0.0
        tmaxcav  = -99.0
        tmaxcc = 0.0
        tmaxm = -99.0
        tmaxpav  = -99.0
        tmaxpc = 0.0
        tmaxsum = 0.0
        tmaxx = -99.0
        tmean = -99.0
        tmeanav = -99.0
        tmeancc = 0.0
        tmeand = 0.0
        tmeane  = 0.0
        tmeanec = 0.0
        tmeang  = 0.0
        tmeangc = 0.0
        tmeannum = 0.0
        tmeanpc = 0.0
        tmeansum = 0.0
        tmincav  = -99.0
        tmincc = 0.0
        tminm = 999.0
        tminn = 99.0
        tminpav  = -99.0
        tminpc = 0.0
        tminsum = 0.0
        tnad = 0.0
        tnoxc = 0.0
        tofixc = 0.0
        tominc = 0.0
        tominfomc = 0.0
        tominsom1c = 0.0
        tominsom2c = 0.0
        tominsom3c = 0.0
        tominsomc = 0.0
        tratio = 0.0
        trwup = 0.0
        tt = 0.0
        tt20 = -99.0
        ttgem = 0.0
        ttcum = 0.0
        ttcumws = 0.0 !LPM 31JUL2015 Added to have a new clock with water stress
        ttcumls = 0.0 !LPM 12JUL2015 Added to have a new clock using a different optimum temperature for leaf size
        ttd = 0.0
        ttnext = 0.0
        twad = 0.0
        uh2o = 0.0
        unh4 = 0.0
        uno3 = 0.0
        vanc = 0.0
        vcnc = 0.0
        vmnc = 0.0
        vnad = 0.0
        vnam = -99.0
        vnamm = -99.0
        vnpcm = -99.0
        vnpcmm = -99.0
        vpdfp = 1.0
        vwad = 0.0
        vwam = -99.0
        vwamm = -99.0
        wfg = 1.0
        wfgcc = 0.0
        wfgpav = 1.0
        wfgpc = 0.0
        wflf = 0.0
        wfp = 1.0
        wfpcav = 1.0
        wfppav = 1.0
        wfpcc = 0.0
        wfppc = 0.0
        wsdays = 0
        wupr = 1.0
        
        h2ocf = -99.0
        no3cf = -99.0
        nh4cf = -99.0
        rtno3 = -99.0
        rtnh4 = -99.0
        no3mn = -99.0
        nh4mn = -99.0

       END SUBROUTINE CS_SeasInit_VarInit
