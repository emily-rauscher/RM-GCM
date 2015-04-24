C************************************************************             
C                 SUBROUTINE INIBAL                                       
C************************************************************             
      SUBROUTINE INIBAL                                                   
C                                                                         
C     This subroutine reads data used for balancing and                   
C     calculates arrays needed for balancing.                             
C                                                                         
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
C     Array ordering in SPECTR must correspond to that in GRIDP.          
C                                                                         
      COMMON/SPECTR/Z(IGB),D(IGB),T(IGB),TRA(IGB,NTRAC),SP(IGA),GS(IGA)   
     :              ,SPA(IGA),VP(IGA),DTE(IGB),TT(IGB)                    
     :              ,TRAT(IGB,NTRAC),DT(IGB),ZT(IGB)                      
     :              ,ZMI(IGB),DMI(IGB),TMI(IGB)                           
     :              ,TRAMI(IGB,NTRAC),SPMI(IGA)                           
      COMPLEX Z,D,T,TRA,SP,GS,SPA,VP,DTE,TT,TRAT,DT,ZT                    
     :       ,ZMI,DMI,TMI,TRAMI,SPMI                                      
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
C     Switches counters and constants controlling type and frequency of   
C     model output                                                        
C                                                                         
      COMMON/OUTCON/RNTAPE,NCOEFF,NLAT,INLAT,INSPC                        
     +              ,RNTAPO                                               
     +              ,KOUNTP,KOUNTE,KOUNTH,KOUNTR                          
     +              ,KOUTP,KOUTE,KOUTH,KOUTR,DAY                          
     +              ,SQR2,RSQR2,EAM1,EAM2,TOUT1,TOUT2,RMG                 
     +              ,LSPO(NL),LGPO(NL)                                    
     $              ,LSHIST,LMINIH                                        
      LOGICAL LSHIST,LMINIH                                               
      LOGICAL LSPO,LGPO                                                   
C                                                                         
C                                                                         
C     Constants and arrays needed for balancing                           
C                                                                         
      COMMON/BALAN/BFILT(NL),RGT0(NL),RG(NL2),TMEAN(NL)                   
     +            ,EP1(IGA),EP2(IGA),KBAL,MFTBAL,SRGT0,LTBAL              
      LOGICAL LTBAL                                                       
C                                                                         
      PARAMETER(IDI2=IDI+IDI,NLT=NL+NL)                                   
C                                                                         
      DIMENSION TBM1(NL2),WA1(IDI),WA2(NL)                                
      INTEGER IWA1(IDI),IWA2(NL)                                          
      DIMENSION TBME(NL,NL)                                               
      EQUIVALENCE (TBM1(1),TBME(1,1))                                     
      DIMENSION AU(IDJ),SM(IDJ),PMNRE(IDJ)                                
      DIMENSION AUE(IDI,IDI)                                              
      EQUIVALENCE (AU(1),AUE(1,1))                                        
C                                                                         
      DIMENSION ZDATN(IGN),ZDAT1(IGM)                                     
C                                                                         
      NAMELIST/INPBL/ KBAL,LTBAL,TMEAN                                    
      NAMELIST/TMPSP/ ZDATN,ZDAT1                                         
      NAMELIST/WVORT/ ZDATN                                               
C                                                                         
  219 FORMAT(/' BALANCING FROM TEMPERATURE AND SURFACE PRESSURE TO',      
     +       ' OBTAIN VORTICITY')                                         
  220 FORMAT(/' BALANCING FROM VORTICITY TO OBTAIN TEMPERATURE AND',      
     +       ' SURFACE PRESSURE')                                         
  234 FORMAT(/' ***ABORT*** BALANCING ATTEMPTED WITH OROGRAPHY'/          
     +       ' TEMPERATURE FIELD WOULD CONTAIN 2-GRID WAVE')              
C                                                                         
C     Set default values and override as desired through NAMELIST input   
C                                                                         
      KBAL=0                                                              
      LTBAL=.FALSE.                                                       
      DO 21 L=1,NL                                                        
         TMEAN(L)=250.0                                                   
   21 CONTINUE                                                            
C                                                                         
C      write(*,*) 'now reading fort.7 in inibal'
      READ(7,INPBL)                                                       
      WRITE(2,INPBL)                                                      
C                                                                         
      IF (     LTBAL) WRITE(2,219)                                        
      IF (.NOT.LTBAL) WRITE(2,220)                                        
C                                                                         
C     Make TMEAN dimensionless                                            
C                                                                         
      DO 62 L=1,NL                                                        
        TMEAN(L)=TMEAN(L)/CT                                              
   62 CONTINUE                                                            
C                                                                         
C     Read zonally averaged state from NAMELIST WVORT or TMPSP.           
C     Assumes data is non-dimensionalised spectral coefficients           
C     and that temperature (if used) includes layer mean.                 
C                                                                         
      IF(.NOT.LTBAL) THEN                                                 
C         write(*,*) '2nd reading in inibal'
         READ (7,WVORT)                                                   
         WRITE(2,WVORT)                                                   
         M=0                                                              
         DO 800 IHEM=1,NHEM                                               
            DO 801 L=1,NL                                                 
               K=NWJ2*(IHEM-1)+(L-1)*IGA                                  
               DO 802 I=1,IDM                                             
                  K=K+1                                                   
                  M=M+1                                                   
                  Z(K)=ZDATN(M)                                           
  802          CONTINUE                                                   
  801       CONTINUE                                                      
  800    CONTINUE                                                         
         I=1                                                              
         DO 170 L=1,NL                                                    
            T(I)=T(I)+SQR2*(TMEAN(L)-T0(L))                               
            I=I+IGA                                                       
  170    CONTINUE                                                         
      ELSE                                                                
C         write(*,*) 'third reading in inibal'
         READ (7,TMPSP)                                                   
         WRITE(2,TMPSP)                                                   
         M=0                                                              
         DO 803 IHEM=1,NHEM                                               
            DO 804 L=1,NL                                                 
               K=NWJ2*(IHEM-1)+(L-1)*IGA                                  
               DO 805 I=1,IDM                                             
                  K=K+1                                                   
                  M=M+1                                                   
                  T(K)=ZDATN(M)                                           
  805          CONTINUE                                                   
  804       CONTINUE                                                      
  803    CONTINUE                                                         
         M=0                                                              
         DO 806 IHEM=1,NHEM                                               
            K=NWJ2*(IHEM-1)                                               
            DO 807 I=1,IDM                                                
               K=K+1                                                      
               M=M+1                                                      
               SP(K)=ZDAT1(M)                                             
  807       CONTINUE                                                      
  806    CONTINUE                                                         
         I=1                                                              
         DO 172 L=1,NL                                                    
            T(I)=T(I)-SQR2*T0(L)                                          
            I=I+IGA                                                       
  172    CONTINUE                                                         
      ENDIF                                                               
      IL=1                                                                
      DO 174 L=1,NL                                                       
         Z(IL)=Z(IL)+EZ                                                   
         IL=IL+IGA                                                        
  174 CONTINUE                                                            
C                                                                         
      IF (KBAL.EQ.0) RETURN                                               
C                                                                         
      IF(.NOT.LTBAL) THEN                                                 
C                                                                         
C        Set values required in BALANC.                                   
C        With orography the balanced temperature field contains a         
C        2-grid wave in the vertical. ABORT if this is attempted.         
C                                                                         
         MAXIND=ICAMAX(IGA,GS,1)                                          
         GSMAX=ABS(GS(MAXIND))                                            
         IF(GSMAX.GT.1.0E-10) THEN                                        
           WRITE(2,234)                                                   
           CALL ABORT                                                     
         ENDIF                                                            
         DO 90 L=1,NL2                                                    
            TBM1(L)=G(L)                                                  
   90    CONTINUE                                                         
         CALL MATINV(TBME,NL,NL,IWA2,WA2)                                 
         DO 92 L=1,NL2                                                    
            RG(L)=TBM1(L)                                                 
   92    CONTINUE                                                         
         DO 91 L=1,NL                                                     
            BFILT(L)=0.                                                   
   91    CONTINUE                                                         
         BFILT(1)=1.                                                      
         BFILT(2)=1.0                                                     
         IF (NLM.GE.2) THEN                                               
            TEMPP=1.0                                                     
            DO 93 I=2,NLM                                                 
               DO 94 J=2,I                                                
                  TEMP=BFILT(J)                                           
                  BFILT(J)=TEMP+TEMPP                                     
                  TEMPP=TEMP                                              
   94          CONTINUE                                                   
               BFILT(I+1)=1.0                                             
   93       CONTINUE                                                      
         ENDIF                                                            
         FACT=-1.0                                                        
         DO 95 I=2,NL                                                     
            BFILT(I)=BFILT(I)*FACT                                        
            FACT=-FACT                                                    
   95    CONTINUE                                                         
         SRGT0=0.                                                         
         IG=0                                                             
         DO 98 L=1,NL                                                     
            TRGT0=0.                                                      
            DO 97 M=1,NL                                                  
               IG=IG+1                                                    
               TRGT0=TRGT0+RG(IG)*T0(M)                                   
   97       CONTINUE                                                      
            RGT0(L)=TRGT0                                                 
            SRGT0=SRGT0+TRGT0*BFILT(L)                                    
   98    CONTINUE                                                         
      ELSE                                                                
C                                                                         
C        Set values required in TBAL                                      
C                                                                         
         MWV1=1+MOCT                                                      
         MFTBAL=9                                                         
C         REWIND 9                                                        
C                                                                         
         IE=0                                                             
         DO 410 MP=1,MFTBAL,MOCT                                          
            AM=MP-1                                                       
            AMSQ=AM*AM                                                    
            AN=AM+1.0                                                     
            DO 409 NP=MP,NFP,MH                                           
               IE=IE+1                                                    
               ANSQ=AN*AN                                                 
               EP2(IE)=2.0*SQRT((ANSQ-AMSQ)/                              
     +                 (4.0*ANSQ-1.0))*(1.0-1.0/AN)                       
               AN=AN+1.0                                                  
               ANSQ=AN*AN                                                 
               EP1(IE)=2.0*SQRT((ANSQ-AMSQ)/                              
     +                 (4.0*ANSQ-1.0))*(1.0+1.0/AN)                       
               AN=AN+1.0                                                  
 409        CONTINUE                                                      
 410     CONTINUE                                                         
         IE=(NFP-1)/2+1                                                   
         DO 1220 MP=MWV1,MFTBAL,MOCT                                      
            J2=(NFP-MP)/2+1                                               
            J22=J2*J2                                                     
            DO 1202 I=1,J22                                               
               AU(I)=0.0                                                  
               SM(I)=0.0                                                  
 1202       CONTINUE                                                      
            IJ=1                                                          
            DO 1205 J=1,J2                                                
               IE=IE+1                                                    
               SM(IJ)=EP2(IE)                                             
               AU(IJ)=EP2(IE)*EP2(IE)+EP1(IE)*EP1(IE)                     
               IF(J.GT.1)AU(IJ-1)=EP1(IE-1)*EP2(IE)                       
               IF (J.NE.J2) THEN                                          
                  SM(IJ+1)=EP1(IE)                                        
                  AU(IJ+1)=EP1(IE)*EP2(IE+1)                              
               ENDIF                                                      
               IJ=IJ+J2+1                                                 
 1205       CONTINUE                                                      
            CALL MATINV(AUE,J2,J2,IWA1,WA1)                               
            IJ=0                                                          
            IUS=1                                                         
            DO 1210 I=1,J2                                                
               DO 1208 J=1,J2                                             
                  TAL=0.0                                                 
                  IJ=IJ+1                                                 
                  IU=IUS                                                  
                  IM=J                                                    
                  DO 1207 K=1,J2                                          
                     TAL=TAL+AU(IU)*SM(IM)                                
                     IU=IU+1                                              
                     IM=IM+J2                                             
 1207             CONTINUE                                                
                  PMNRE(IJ)=TAL                                           
 1208          CONTINUE                                                   
               IUS=IUS+J2                                                 
 1210       CONTINUE                                                      
            WRITE(99)(PMNRE(I),I=1,J22)                                   
 1220    CONTINUE                                                         
C                                                                         
         IF(NHEM.EQ.2) THEN                                               
            IE=1+NWJ2                                                     
            AN=2.0                                                        
            DO 420 NP=3,NFP,MH                                            
               IE=IE+1                                                    
               ANSQ=AN*AN                                                 
               EP2(IE)=2.0*SQRT(ANSQ/(4.0*ANSQ-1.0))*(1.0-1.0/AN)         
               AN=AN+1.0                                                  
               ANSQ=AN*AN                                                 
               EP1(IE)=2.0*SQRT(ANSQ/(4.0*ANSQ-1.0))*(1.0+1.0/AN)         
               AN=AN+1.0                                                  
  420       CONTINUE                                                      
            DO 430 MP=MWV1,MFTBAL,MOCT                                    
               AM=MP-1                                                    
               AMSQ=AM*AM                                                 
               AN=AM+1.0                                                  
               DO 429 NP=MP,NFP,MH                                        
                  IE=IE+1                                                 
                  ANSQ=AN*AN                                              
                  EP1(IE)=2.0*SQRT((ANSQ-AMSQ)/                           
     +                    (4.0*ANSQ-1.0))*(1.0+1.0/AN)                    
                  AN=AN+1.0                                               
                  ANSQ=AN*AN                                              
                  EP2(IE)=2.0*SQRT((ANSQ-AMSQ)/                           
     +                    (4.0*ANSQ-1.0))*(1.0-1.0/AN)                    
                  AN=AN+1.0                                               
  429          CONTINUE                                                   
  430       CONTINUE                                                      
            IE=1+NWJ2                                                     
            J2=IDM                                                        
            J22=J2*J2                                                     
            J2L=J2-1                                                      
            J22L=J2*J2L                                                   
            DO 1222 I=1,J22                                               
               AU(I)=0.0                                                  
               SM(I)=0.0                                                  
 1222       CONTINUE                                                      
            IJ=1                                                          
            DO 1225 J=1,J2L                                               
               IE=IE+1                                                    
               IJB=IJ+J-1                                                 
               SM(IJB)=EP2(IE)                                            
               AU(IJ)=EP1(IE)*EP1(IE) + EP2(IE)*EP2(IE)                   
               IF(J.GT.1) AU(IJ-1)=EP1(IE-1)*EP2(IE)                      
               SM(IJB+1)=EP1(IE)                                          
               IF(J.LT.J2L) AU(IJ+1)=EP1(IE)*EP2(IE+1)                    
               IJ=IJ+J2                                                   
 1225       CONTINUE                                                      
            CALL MATINV(AUE,J2L,J2L,IWA1,WA1)                             
            IJ=0                                                          
            IUS=1                                                         
            DO 1230 I=1,J2L                                               
               DO 1228 J=1,J2                                             
                  TAL=0.0                                                 
                  IJ=IJ+1                                                 
                  IU=IUS                                                  
                  IM=J                                                    
                  DO 1227 K=1,J2L                                         
                     TAL=TAL + AU(IU)*SM(IM)                              
                     IU=IU+1                                              
                     IM=IM+J2                                             
 1227             CONTINUE                                                
                  PMNRE(IJ)=TAL                                           
 1228          CONTINUE                                                   
               IUS=IUS+J2L                                                
 1230       CONTINUE                                                      
            WRITE(99) (PMNRE(I),I=1,J22L)                                 
         ENDIF                                                            
C                                                                         
      ENDIF                                                               
C                                                                         
      END                                                                 
