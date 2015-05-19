C**********************************************************               
C             SUBROUTINE INIVARPARAM                                          
C**********************************************************               
      subroutine INIVARPARAM                                                 
C-----------------------------------------------------------------------  
C     Subroutine to initialise various model parameters            
C-----------------------------------------------------------------------  
C                                                                         
C     Determines model resolution                                         
C                                                                         
      include 'params.i'
C      PARAMETER(NN=21,MM=21,NHEM=2,NL=5,MOCT=1,MG=64,JG=16,NWJ2=121       
C     P         ,NCRAY=64,JGL=JG,NTRAC=1,NLEVRF=1)                         
                                                                          
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
C                                                                         
C     Note that RD and GASCON are identical and CPD is set from RD,AKAP.  
      COMMON        SQ(NNP),RSQ(NNP),SIGMAH(NLM),SIGMA(NL)                
     +              ,T01S2(NLM),T0(NL),ALPHA(NL),DSIGMA(NL),RDSIG(NL)     
     +              ,TKP(NL),C(NL2),SQH(NNP)                              
     +              ,MF,MFP,JZF,NF,NFP                                    
     +              ,AKAP,GA,GASCON,RADEA,WW,P0,PFAC,EZ,AIOCT             
     +              ,RD,RV,CPD,CLATNT                                     
     +              ,LRSTRT,LSHORT,LTVEC,LSTRETCH                         
     +              ,LFLUX                                                
     +              ,LBALAN,LRESTIJ                                       
     +              ,LCLIM, LPERPET, L22L,LOROG ,LCSFCT                   
     +              ,LNOISE                                               
      COMPLEX EZ,AIOCT                                                    
      LOGICAL LRSTRT,LSHORT,LTVEC,LSTRETCH,LBALAN,LRESTIJ                 
     +       ,LFLUX,LNOISE                                                
     +       ,LCLIM, LPERPET, L22L,LOROG,LCSFCT                           
C                                                                         
C                                                                         
C     Constant arrays and variables associated with time and vertical     
C     differencing. Also counters.                                        
C                                                                         
      COMMON/BATS/  BM1(IDE),AK(NNP),AQ(NL2),G(NL2),TAU(NL2)              
     +              ,KOUNT,KITS,KSTART,KTOTAL,KRUN,BEGDAY,ITSPD           
     +              ,DELT,DELT2,CV,CG,CT,CQ,PNU,PNU2,PNU21                
     +              ,NTRACO,CTRA(NTRAC),KOLOUR(NTRAC),RGG(NL2)            
     +              ,BEGDOY,DOY                                           
C                                                                         
C                                                                         
C                                                                         
       COMMON/VARPARAM/OOM_IN, LPLOTMAP,NLPLOTMAP_IN,RFCOEFF_IN, 
     & NTSTEP_IN, NSKIP_IN, BOTRELAXTIME, FBASEFLUX, FORCE1DDAYS, 
     & OPACIR_POWERLAW, OPACIR_REFPRES, SOLC_IN, TOAALB, 
     & PORB, OBLIQ, ECCEN 
      
       LOGICAL LPLOTMAP

       COMMON/MAG/ LBDRAG,BFIELD,TDRAG_MIN,RAMPUP
       LOGICAL LBDRAG

       COMMON/BINVAL/LBIN,PORBST,ECCPL,ECCST,SMAPL,SMAST,STMASS1,
     & STMASS2,STRAD1,STRAD2,STTTEMP1,STTEMP2

       LOGICAL LBIN
C       COMMON/VARPARAM2/LDIUR_IN, OCFLUX_IN, ROUGHLENGTH,NSWALBEDO_FUNC,
C     & SURFDEPTH_TOP, SURFDEPTH_BOT, EFFSOILCONDUCT,
C     & TSURF_INIT_TOP, TSURF_INIT_BOT, SNET_INIT

C       LOGICAL LDIUR_IN

       NAMELIST/INVARPARAM/OOM_IN, LPLOTMAP,NLPLOTMAP_IN,RFCOEFF_IN, 
     & NTSTEP_IN, NSKIP_IN, BOTRELAXTIME, FBASEFLUX, FORCE1DDAYS, 
     & OPACIR_POWERLAW, OPACIR_REFPRES, SOLC_IN, TOAALB, 
     & PORB, OBLIQ, ECCEN 

       NAMELIST/INMAG/ LBDRAG,BFIELD,TDRAG_MIN,RAMPUP

       NAMELIST/INBINVAL/LBIN,PORBST,ECCPL,ECCST,SMAPL,SMAST,STMASS1,
     & STMASS2, STRAD1,STRAD2,STTEMP1,STTEMP2
C       NAMELIST/INVARPARAM2/LDIUR_IN, OCFLUX_IN, ROUGHLENGTH,
C      & NSWALBEDO_FUNC,
C     & SURFDEPTH_TOP, SURFDEPTH_BOT, EFFSOILCONDUCT,
C     & TSURF_INIT_TOP, TSURF_INIT_BOT, SNET_INIT
       
C Switch to 
       LPLOTMAP=.TRUE.

C  240 FORMAT(' SIMPLIFIED RADIATIVE SCHEME EMPLOYED:')      
C  243 FORMAT(' SW MODEL 1: CLEAR SKY + SURFACE REFLECTION'     
C     :,' ALBEDO ALBSW1=',f6.1)  
C  244  FORMAT(' SW MODEL 2: SEMI-INFINITE SCATTERING ATMOSPHERE WITH'     
C     :,' TOTAL OPTICAL THICKNESS ABSSW2=',f6.1
C     :,', SCATTERING COEFF SCATSW2= ',f4.1 
C     :,', ASYMMETRY FACTOR ASYMSW2= ',f4.1)  
C  245 FORMAT(' LW MODEL 1: SINGLE ABSORPT. COEFF. (*PRESS) ABSLW1= ',e8.2)     
C  246  FORMAT(' LW MODEL 2: PLANCK MEAN OPACITY TABLE') 
C 247   FORMAT(' NSWMODEL INDEX INAPPROPRIATE:',I4)
C 248   FORMAT(' NLWMODEL INDEX INAPPROPRIATE:',I4)

C       write(*,*) 'reading fort.7 in inivarparam'
       READ (7,INVARPARAM)                                                     
       WRITE(2,INVARPARAM)             

C       write(*,*) 'second read in inivarparam'
       READ (7,INMAG)
       WRITE(2,INMAG)

       READ(7,INBINVAL)
       WRITE(2,INBINVAL)

C      WRITE(2,240)                                                        

C      IF(NSWMODEL.EQ.1) THEN
C         WRITE(2,243) ALBSW1
C      ELSE IF(NSWMODEL.EQ.2) THEN
C         WRITE(2,244) ABSSW2, SCATSW2, ASYMSW2
C      ELSE
C      WRITE(2,247) NSWMODEL
C      CALL ABORT
C      ENDIF
 
C      IF(NLWMODEL.EQ.1) THEN
C         WRITE(2,245) ABSLW1
C      ELSE IF(NLWMODEL.EQ.2) THEN
C         WRITE(2,246)
C      ELSE
C      WRITE(2,248) NLWMODEL
C      CALL ABORT
C      ENDIF
    
      END                                                                 
