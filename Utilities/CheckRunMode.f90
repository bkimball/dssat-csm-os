Subroutine CheckRunMode(RNMODE)

Implicit None
INTEGER COUNT
CHARACTER*1 RNMODE
CHARACTER*120 MSG(66)

DATA MSG / &
"-----------------------------------------------------------------------------", &
"The command line arguments you provided are invalid.                         ", & 
"Please use the following syntax:                                             ", &
"-------------------------------------------------------                      ", &
"                                                                             ", &
"  Model_binary <Model> Runmode <argA> <argB> <FileCTR>                       ", &
"                                                                             ", &
"-----------------------------------------------------------------------------", &
"Details:                                                                     ", &
"  <Model>   - optional                                                       ", &
"            - 8-character name of crop model (e.g., MZIXM047 or WHAPS047).   ", &
"            - If model name is blank or invalid, the default will be used.   ", &
"                                                                             ", &
"  Runmode   - required                                                       ", &
"            - 1-character run mode code                                      ", &
"            - see table below for valid values and required arguments        ", &
"Run                                                                          ", &
"mode argA       argB  Description                                            ", &
"---- ---------  ----- ------------------------------------------------------ ", &
" A   FileX      NA    All: Run all treatments in the specified FileX.        ", &
" B   BatchFile  NA    Batch: Batchfile lists experiments and treatments.     ", &
" C   FileX      TrtNo Command line: Run single FileX and treatment #.        ", &
" D   TempFile   NA    Debug: Skip input module and use existing TempFile.    ", &
" E   BatchFile  NA    Sensitivity: Batchfile lists FileX and TrtNo.          ", &
" F   BatchFile  NA    Farm model: Batchfile lists experiments and treatments.", &
" G   FileX      TrtNo Gencalc: Run single FileX and treatment #.             ", &
" I   NA         NA    Interactive: Interactively select FileX and TrtNo.     ", &
" L   BatchFile  NA    Gene-based model (Locus): Batchfile for FileX and TrtNo", &
" N   BatchFile  NA    Seasonal analysis: Batchfile lists FileX and TrtNo.    ", &
" Q   BatchFile  NA    Sequence analysis: Batchfile lists FileX & rotation #. ", &
" S   BatchFile  NA    Spatial: Batchfile lists experiments and treatments.   ", &
" T   BatchFile  NA    Gencalc: Batchfile lists experiments and treatments.   ", &
"                                                                             ", &
"  BatchFile - Name of DSSAT batch file with list of exeriments and treatments", &
"                (e.g., DSSBATCH.v47)                                         ", &
"            - Current directory, 30 characters maximum                       ", &
"                                                                             ", &
"  FileX     - Name of Experimental file (e.g., UFGA7801.SBX)                 ", &
"            - Current directory, 12-character DSSAT naming convention        ", &
"                                                                             ", &
"  TempFile  - Name of temporary I/O file, normally generated by the input    ", &
"                module (e.g., DSSAT47.INP)                                   ", &
"            - Current directory, 30 characters maximum                       ", &
"                                                                             ", &
"  TrtNo     - Treatment # (integer) in specified FileX to be simulated       ", &
"                                                                             ", &
"  <FileCTR> - optional                                                       ", &
"            - path + filename of external file which contains overrides for  ", &
"                simulation controls.                                         ", &
"            - This option is available with all run modes except D and I.    ", &
"            - Default file (DSCSM047.CTR) is found in the root directory.    ", &
"            - 120 characters maximum.                                        ", &
"-----------------------------------------------------------------------------", &
" Example #1:                                                                 ", &
" DSCSM047 B DSSBATCH.V47                                                     ", &
" Effect: Run in batch mode. Name of the batch file is DSSBATCH.V47.          ", &
"                                                                             ", &
" Example #2:                                                                 ", &
" DSCSM047 MZIXM047 A UFGA8201.MZX                                            ", &
" Effect: Run all treatments in experiment UFGA8201.MZX using IXIM model.     ", &
"                                                                             ", &
" Example #3:                                                                 ", &
" DSCSM047 Q DSSBATCH.V47 DSCSM047.CTR                                        ", &
" Effect: Run sequence simulation listed in DSSBATCH.V47 using the            ", &
"           simulation control options specified by DSCSM047.CTR              ", &
"-----------------------------------------------------------------------------"/

  COUNT = SIZE(MSG)

  IF (INDEX('ABCDEFGILNQSTabcdefginlqst',RNMODE) .GT. 0) RETURN

  WRITE(*,'(100(/,A))') MSG
  CALL WARNING(COUNT, "CSM", MSG)
  CALL ERROR ("CSM",90,"",0)
  RETURN
End Subroutine CheckRunMode