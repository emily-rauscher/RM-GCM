C**********************************************************               
C             SUBROUTINE INIVARPARAM                                          
C**********************************************************               
      subroutine INIVARPARAM                                                 
      use tracerclds, only : ADYE, RHODYE, PDYEFIX, PDYEUPPER,
     &                       TRELAXORB, KDYERAD, KDYESPEC, LDYERAD,
     &                       KDYEMODE, TCONDDYE, TEVAPDYE, LDYEFIX,
     &                       LBENCH, PDYEBUD, NDYEBUD,
     &                       L1COLRT, ZEN1COL, NPAIRMX, NDYERAD,
     &                       KSLOTRAD, ISPRAD, LDYELOGN
      use datapaths, only : DATA_DIR
C-----------------------------------------------------------------------  
C     Subroutine to initialise various model parameters            
C-----------------------------------------------------------------------  
C                                                                         
C     Determines model resolution                                         
C                                                                         
      include 'params.i'                                                                                    
C                                                                         
C                                                                         
C     Sets basic constants, especially those needed for array dimensions  
C                                                                         
      PARAMETER(MH=2,PI=3.14159265359,PI2=2.0*PI                          
     +,NNP=NN+1,MGPP=MG+2,JGP=JG+1,JGG=JG*NHEM,JGGP=JGG+1,MJP=NWJ2+NWJ2   
     +,NLM=NL-1,NLP=NL+1,NLPP=NL+2,NLA=NL+3,NLB=NL+4,NL2=NL*NL            
     +,IDA=(MG+MG+MG)/2+1,IDB=NWJ2*NL,IDC=IDB+IDB,IDD=MGPP*NL             
     +,IDE=NL2*NN,IDF=NCRAY*(MG+1),IDG=JG*NL,IDH=JG*MG                    
     +,IDI=NNP/2,IDJ=IDI*IDI,IDK=NL*IDI,IDL=MGPP/2,IDM=NNP/2,IDN=IDM*NL   
     +,NWW=1+(MM-1)/MOCT)                                                 
      PARAMETER(IGA=NWJ2*NHEM,IGB=IDB*NHEM,IGC=MGPP*NHEM,IGD=IDD*NHEM     
     +,IGG=IDG*NHEM,IGL=IDL*NHEM,IGM=IDM*NHEM,IGN=IDN*NHEM                
     +,IGO=IGA+IGA,IGP=IGB+IGB,NFTWG=(5+NTRAC)*NL+3                       
     +,NFTGW=(6+3*NTRAC)*NL+2,NFTGD=(3+NTRAC)*NL,NLTR=NL*NTRAC)           
C     Number of 2D (surface) output fields. This value is                 
C     Doubled due to averaged and instantaneous fields.                   
      PARAMETER (N2DFLD=21,NGRPAD=N2DFLD*2*IGC)                           
C
C      
C     Basic planetary parameters for run plus information about           
C     vertical grid structure                                                                                                                      
C     Note that RD and GASCON are identical and CPD is set from RD,AKAP. 

      COMMON        SQ(NNP),RSQ(NNP),SIGMAH(NLM),SIGMA(NL)                
     +              ,T01S2(NLM),T0(NL),ALPHA(NL),DSIGMA(NL),RDSIG(NL)     
     +              ,TKP(NL),C(NL2),SQH(NNP)                              
     +              ,MF,MFP,JZF,NF                                    
     +              ,AKAP,GA,GASCON,RADEA,WW,PFAC,EZ,AIOCT             
     +              ,RD,RV,CPD,CLATNT                                     
     +              ,P0,LRSTRT,LSHORT,LTVEC,LSTRETCH                         
     +              ,LFLUX                                                
     +              ,LBALAN,LRESTIJ                                       
     +              ,LCLIM, LPERPET, L22L,LOROG ,LCSFCT                   
     +              ,LNOISE,NFP                                               
      COMPLEX EZ,AIOCT                                                    
      LOGICAL LRSTRT,LSHORT,LTVEC,LSTRETCH,LBALAN,LRESTIJ                 
     +       ,LFLUX,LNOISE                                                
     +       ,LCLIM, LPERPET, L22L,LOROG,LCSFCT                           
                                                           
C     Constant arrays and variables associated with time and vertical     
C     differencing. Also counters.                                        

      COMMON/BATS/  BEGDAY,CTRA(NTRAC),BM1(IDE),AK(NNP),AQ(NL2),G(NL2)              
     +              ,TAU(NL2),KOUNT,KITS,KSTART,KTOTAL,KRUN,ITSPD           
     +              ,DELT,DELT2,CV,CG,CT,CQ,PNU,PNU2,PNU21                
     +              ,NTRACO,KOLOUR(NTRAC),RGG(NL2)            
     +              ,BEGDOY,DOY                                           
C

       COMMON/VARPARAM/OOM_IN, LPLOTMAP,NLPLOTMAP_IN,RFCOEFF_IN, 
     & NTSTEP_IN, NSKIP_IN, BOTRELAXTIME, FBASEFLUX, FORCE1DDAYS, 
     & OPACIR_POWERLAW, OPACIR_REFPRES, SOLC_IN, TOAALB, 
     & PORB, OBLIQ, ECCEN, TAULIMIT
      
       LOGICAL LPLOTMAP

C      Radiatively active tracer cloud. LDYERAD=F leaves the RT on the prescribed
C      ropprmulti field, and is off by default. KDYEMODE=1 has one tracer per species:
C      KDYERAD picks the active one (2 = first dye) and the rest stay passive
C      (i.e., they don't contribute to the radiative transfer). With KDYEMODE=2
C      the dyes are vapour/solid pairs and KDYERAD must be a solid slot.

       COMMON/MAG/ BFIELD,TDRAG_MIN,RAMPUP,LBDRAG
       LOGICAL LBDRAG

       COMMON/BINVAL/PORBST,ECCPL,ECCST,SMAPL,SMAST,STMASS1,
     & STMASS2,STRAD1,STRAD2,STTEMP1,STTEMP2,LBIN

       LOGICAL LBIN

       NAMELIST/INVARPARAM/OOM_IN, LPLOTMAP,NLPLOTMAP_IN,RFCOEFF_IN,
     & NTSTEP_IN, NSKIP_IN, BOTRELAXTIME, FBASEFLUX, FORCE1DDAYS,
     & OPACIR_POWERLAW, OPACIR_REFPRES, SOLC_IN, TOAALB,
     & PORB, OBLIQ, ECCEN,
     & ADYE, PDYEFIX, PDYEUPPER, TRELAXORB,
     & LDYERAD, KDYERAD, KDYESPEC,
     & KDYEMODE, TCONDDYE, TEVAPDYE, LDYEFIX, LBENCH,
     & PDYEBUD, NDYEBUD, L1COLRT, ZEN1COL, LDYELOGN, DATA_DIR

       NAMELIST/INMAG/ LBDRAG,BFIELD,TDRAG_MIN,RAMPUP

       NAMELIST/INBINVAL/LBIN,PORBST,ECCPL,ECCST,SMAPL,SMAST,STMASS1,
     & STMASS2, STRAD1,STRAD2,STTEMP1,STTEMP2

       LPLOTMAP=.TRUE.

C      Defaults before the namelist read, so a fort.7 predating these
C      keys still gives a defined configuration.
       DO 5 KK=1,NDYE+1
          ADYE(KK)=1.0E-6
    5  CONTINUE
       LDYEFIX=.TRUE.
       LDYELOGN=.TRUE.
       LBENCH=.FALSE.
       L1COLRT=.FALSE.
       ZEN1COL=0.0
       PDYEBUD=1.0E3
       NDYEBUD=0
C      Guard value only; cloud_properties_set_up.f replaces
       RHODYE=3.0E3
       PDYEFIX=1.0E5
       PDYEUPPER=1.0E6
       TRELAXORB=0.1
       LDYERAD=.FALSE.
       KDYERAD=2
C      One species per vapour/solid pair (KDYEMODE=2); KDYEMODE=1 and
C      a one-pair run read only the first, so "KDYESPEC = 7" still works.
       KDYESPEC=0
       KDYESPEC(1)=7
       KDYEMODE=1
       TCONDDYE=1.0E3
C      Negative means "same as TCONDDYE", resolved after the read.
       TEVAPDYE=-1.0

       DATA_DIR='../'

       READ (7,INVARPARAM)
       IF (TEVAPDYE.LT.0.0) TEVAPDYE=TCONDDYE
C      DATA_DIR gets a trailing '/', so '..' and '../' both work; blank
C      means the run directory itself.
       IF (LEN_TRIM(DATA_DIR) .EQ. 0) DATA_DIR='./'
       IF (DATA_DIR(LEN_TRIM(DATA_DIR):LEN_TRIM(DATA_DIR)) .NE. '/')
     &    DATA_DIR=TRIM(DATA_DIR)//'/'
       WRITE(*,'(2A)') ' Input data directory (DATA_DIR): ', TRIM(DATA_DIR)
       WRITE(2,INVARPARAM)

C      Single-column RT: its clouds come from the condensation curve, as
C      the column has no dye state.
       IF (L1COLRT .AND. LDYERAD) THEN
          WRITE(2,*) 'INIVARPARAM: L1COLRT needs LDYERAD=F'
          STOP 'L1COLRT needs LDYERAD=F'
       ENDIF

C      Cloud scheme. Two-phase dyes come in vapour/solid pairs filling
C      slots 2..NDYE+1, each pair with its own species.
       IF (KDYEMODE.NE.1 .AND. KDYEMODE.NE.2) THEN
          WRITE(2,*) 'INIVARPARAM: KDYEMODE must be 1 or 2:',
     &               KDYEMODE
          STOP 'KDYEMODE out of range'
       ENDIF
       IF (KDYEMODE.EQ.2 .AND. MOD(NDYE,2).NE.0) THEN
          WRITE(2,*) 'INIVARPARAM: KDYEMODE=2 needs NDYE even ',
     &               '(2 per vapour/solid pair):', NDYE
          STOP 'KDYEMODE=2 needs even NDYE'
       ENDIF

C      Species in use: one per pair in KDYEMODE=2, KDYESPEC(1) for every
C      dye in KDYEMODE=1. Each must be a cloud species 1-NCLOUDS, and no
C      species may be used by two pairs.
       IF (KDYEMODE.EQ.2) THEN
          NSPDY=NDYE/2
       ELSE
          NSPDY=1
       ENDIF
       DO 30 KP=1,NSPDY
          IF (KDYESPEC(KP).LT.1 .OR. KDYESPEC(KP).GT.NCLOUDS) THEN
             WRITE(2,*) 'INIVARPARAM: KDYESPEC(',KP,') must be a ',
     &                  'cloud species 1 -',NCLOUDS,' (one per pair):',
     &                  KDYESPEC(KP)
             WRITE(*,*) 'INIVARPARAM: KDYESPEC(',KP,') must be a ',
     &                  'cloud species 1 -',NCLOUDS,' (one per pair):',
     &                  KDYESPEC(KP)
             STOP 'KDYESPEC out of range'
          ENDIF
          DO 25 KQ=1,KP-1
             IF (KDYESPEC(KQ).EQ.KDYESPEC(KP)) THEN
                WRITE(2,*) 'INIVARPARAM: species ',KDYESPEC(KP),
     &                     ' is in KDYESPEC twice'
                WRITE(*,*) 'INIVARPARAM: species ',KDYESPEC(KP),
     &                     ' is in KDYESPEC twice'
                STOP 'KDYESPEC repeats a species'
             ENDIF
   25     CONTINUE
   30  CONTINUE

C      The radiatively active dyes: every pair's solid in KDYEMODE=2,
C      the single tracer KDYERAD in KDYEMODE=1.
       NDYERAD=0
       ISPRAD=0
       IF (LDYERAD) THEN
          IF (KDYEMODE.EQ.2) THEN
             NDYERAD=NSPDY
             DO 35 KP=1,NSPDY
                KSLOTRAD(KP)=2*KP+1
   35        CONTINUE
          ELSE
             IF (KDYERAD.LT.2 .OR. KDYERAD.GT.NDYE+1) THEN
                WRITE(2,*) 'INIVARPARAM: KDYERAD must be 2..NDYE+1:',
     &                     KDYERAD
                STOP 'KDYERAD out of range'
             ENDIF
             NDYERAD=1
             KSLOTRAD(1)=KDYERAD
          ENDIF
          DO 40 KP=1,NDYERAD
             ISPRAD(KDYESPEC(KP))=KP
   40     CONTINUE
       ENDIF

       READ (7,INMAG)
       WRITE(2,INMAG)
       READ(7,INBINVAL)
       WRITE(2,INBINVAL)
       
       TAULIMIT=2.*LOG(1.e-2)/1.66/(10.**(-1.*OOM_IN/NL)-10.**(OOM_IN/NL))

       ! Isaac is guessing
       !TAULIMIT=100

      END                                                                 
