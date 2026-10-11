C**********************************************************               
C             SUBROUTINE INIVARPARAM                                          
C**********************************************************               
      subroutine INIVARPARAM                                                 
      use tracerclds, only : ADYE, RHODYE, PDYEFIX, PDYEUPPER,
     &                       TRELAXORB, KDYERAD, KDYESPEC, LDYERAD,
     &                       KDYEMODE, TCONDDYE, TEVAPDYE
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
     & KDYEMODE, TCONDDYE, TEVAPDYE

       NAMELIST/INMAG/ LBDRAG,BFIELD,TDRAG_MIN,RAMPUP

       NAMELIST/INBINVAL/LBIN,PORBST,ECCPL,ECCST,SMAPL,SMAST,STMASS1,
     & STMASS2, STRAD1,STRAD2,STTEMP1,STTEMP2

       LPLOTMAP=.TRUE.

C      Defaults before the namelist read, so a fort.7 predating these
C      keys still gives a defined configuration.
       DO 5 KK=1,NTRAC
          ADYE(KK)=1.0E-6
    5  CONTINUE
C      Guard value only; cloud_properties_set_up.f replaces
       RHODYE=3.0E3
       PDYEFIX=1.0E5
       PDYEUPPER=1.0E6
       TRELAXORB=0.1
       LDYERAD=.FALSE.
       KDYERAD=2
       KDYESPEC=7
       KDYEMODE=1
       TCONDDYE=1.0E3
C      Negative means "same as TCONDDYE", resolved after the read.
       TEVAPDYE=-1.0

       READ (7,INVARPARAM)
       IF (TEVAPDYE.LT.0.0) TEVAPDYE=TCONDDYE
       WRITE(2,INVARPARAM)


C      Check that the dye species is in range, and that the radiatively active
C      tracer is a valid species. The dye species is used for both settling
C      and radiative feedback
       IF (KDYESPEC.LT.1 .OR. KDYESPEC.GT.13) THEN
          WRITE(2,*) 'INIVARPARAM: KDYESPEC out of range 1-13:',
     &               KDYESPEC
          STOP 'KDYESPEC out of range'
       ENDIF
       IF (LDYERAD) THEN
          IF (KDYERAD.LT.2 .OR. KDYERAD.GT.NTRAC) THEN
             WRITE(2,*) 'INIVARPARAM: KDYERAD must be 2..NTRAC:',
     &                  KDYERAD
             STOP 'KDYERAD out of range'
          ENDIF
       ENDIF

C      Cloud scheme. Two-phase dyes come in vapour/solid pairs filling
C      slots 2..NTRAC, and only a solid (odd) slot may be opaque.
       IF (KDYEMODE.NE.1 .AND. KDYEMODE.NE.2) THEN
          WRITE(2,*) 'INIVARPARAM: KDYEMODE must be 1 or 2:',
     &               KDYEMODE
          STOP 'KDYEMODE out of range'
       ENDIF
       IF (KDYEMODE.EQ.2) THEN
          IF (MOD(NTRAC-1,2).NE.0) THEN
             WRITE(2,*) 'INIVARPARAM: KDYEMODE=2 needs NTRAC odd ',
     &                  '(1 + 2 per vapour/solid pair):', NTRAC
             STOP 'KDYEMODE=2 needs odd NTRAC'
          ENDIF
          IF (LDYERAD .AND. MOD(KDYERAD,2).NE.1) THEN
             WRITE(2,*) 'INIVARPARAM: KDYEMODE=2 needs KDYERAD on a ',
     &                  'solid slot (3,5,...):', KDYERAD
             STOP 'KDYERAD is not a solid slot'
          ENDIF
       ENDIF

       READ (7,INMAG)
       WRITE(2,INMAG)
       READ(7,INBINVAL)
       WRITE(2,INBINVAL)
       
       TAULIMIT=2.*LOG(1.e-2)/1.66/(10.**(-1.*OOM_IN/NL)-10.**(OOM_IN/NL))

       ! Isaac is guessing
       !TAULIMIT=100

      END                                                                 
