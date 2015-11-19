C**********************************************************               
C             SUBROUTINE ENERGY                                           
C**********************************************************               
      SUBROUTINE ENERGY                                                   
C                                                                         
C     Calculates various global diagnostic quantities                     
C     every itstp timesteps.                                              
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
C     Constant arrays and variables associated with time and vertical     
C     differencing. Also counters.                                        
C                                                                         
      COMMON/BATS/  BEGDAY,CTRA(NTRAC),BM1(IDE),AK(NNP),AQ(NL2),G(NL2)              
     +              ,TAU(NL2),KOUNT,KITS,KSTART,KTOTAL,KRUN,ITSPD           
     +              ,DELT,DELT2,CV,CG,CT,CQ,PNU,PNU2,PNU21                
     +              ,NTRACO,KOLOUR(NTRAC),RGG(NL2)            
     +              ,BEGDOY,DOY                                           
C                                                                         
      COMPLEX Q(IGB)                                                      
      EQUIVALENCE (Q(1),TRA(1,1))                                         
      COMPLEX TIG,SPIP,CTIG,CRGS                                          
      COMMON/STATS/GMSP0,GMSPMI,LMASCOR,LMASOLD,LMASPRT                   
      LOGICAL LMASCOR,LMASOLD,LMASPRT                                     
C                                                                         
C-----------------------------------------------------------------------  
C                                                                         
C     First remove planetary vorticity so Z contains relative vorticity   
C                                                                         
      I=1                                                                 
      DO 101 L=1,NL                                                       
         Z(I)=Z(I)-EZ                                                     
         I=I+IGA                                                          
  101 CONTINUE                                                            
C                                                                         
C     Calculate means - PSITOT RMS vorticity                              
C                       CHITOT RMS divergence                             
C                       TMPTOT RMS temperature                            
C                       TOTP  IE+PE potential energy                      
C                       AMSP mean surface pressure                        
C                                                                         
      PSITOT=0.0                                                          
      CHITOT=0.0                                                          
      TMPTOT=0.0                                                          
      TOTP=0.0                                                            
      TOTI=0.0                                                            
      TOTQ=0.0 ! Total moisture diagnostic (as in BGCM5).                 
      AMQ=0.0  ! Accounts for the p0 term excluded from SPA(1) in MGRMLT. 
      IL=1                                                                
      ST2B=0.                                                             
      ST=0.                                                               
      IOFS=0                                                              
      DO 800 IHEM=1,NHEM                                                  
         IG=IOFS                                                          
         DO 30 L=1,NL                                                     
            TPSITT=0.                                                     
            TCHITT=0.                                                     
            TTMPTT=0.                                                     
            TTOTI=0.                                                      
            TTOTQ=0.                                                      
            DSIG=DSIGMA(L)                                                
            DSIGH=0.5*DSIG                                                
            IP=IOFS                                                       
            DO 10 JP=1,NFP,MH                                             
               IG=IG+1                                                    
               IP=IP+1                                                    
               RTIG=REAL(T(IG))                                           
               RZIG=REAL(Z(IG))                                           
               RDIG=REAL(D(IG))                                           
               TPSITT=TPSITT+RZIG*RZIG                                    
               TCHITT=TCHITT+RDIG*RDIG                                    
               TTMPTT=TTMPTT+RTIG*RTIG                                    
               RSPIP=REAL(SPA(IP))                                        
               IF (L.EQ.1) THEN                                           
                  TOTP=TOTP+0.5*RSPIP*REAL(GS(IP))                        
               ENDIF                                                      
               TTOTI=TTOTI+RSPIP*RTIG                                     
               TTOTQ=TTOTQ+RSPIP*REAL(Q(IG))                              
   10       CONTINUE                                                      
            PSITOT=PSITOT+DSIGH*TPSITT                                    
            CHITOT=CHITOT+DSIGH*TCHITT                                    
            TMPTOT=TMPTOT+DSIGH*TTMPTT                                    
            TOTI=TOTI+DSIGH*TTOTI                                         
            TOTQ=TOTQ+DSIGH*TTOTQ                                         
            TPSITT=0.                                                     
            TCHITT=0.                                                     
            TTMPTT=0.                                                     
            TTOTI=0.                                                      
            TTOTQ=0.                                                      
            DO 25 M=MOCT,MF,MOCT                                          
               DO 20 JP=M,NF,MH                                           
                  IG=IG+1                                                 
                  IP=IP+1                                                 
                  TIG=T(IG)                                               
                  CTIG=CONJG(TIG)                                         
                  SPIP=SPA(IP)                                            
                  TPSITT=TPSITT+REAL(Z(IG)*CONJG(Z(IG)))                  
                  TCHITT=TCHITT+REAL(D(IG)*CONJG(D(IG)))                  
                  TTMPTT=TTMPTT+REAL(TIG*CTIG)                            
                  IF (L.EQ.1) THEN                                        
                     CRGS=CONJG(GS(IP))                                   
                     TOTP=TOTP+REAL(SPIP*CRGS)                            
                  ENDIF                                                   
                  TTOTI=TTOTI+REAL(SPIP*CTIG)                             
                  TTOTQ=TTOTQ+REAL(SPIP*CONJG(Q(IG)))                     
   20          CONTINUE                                                   
   25       CONTINUE                                                      
            PSITOT=PSITOT+DSIG*TPSITT                                     
            CHITOT=CHITOT+DSIG*TCHITT                                     
            TMPTOT=TMPTOT+DSIG*TTMPTT                                     
            TOTI=TOTI+DSIG*TTOTI                                          
            TOTQ=TOTQ+DSIG*TTOTQ                                          
            IF (IHEM.EQ.1) THEN                                           
               RTL=REAL(T(IL))                                            
               ST2B=ST2B+T0(L)*RTL*DSIG                                   
               ST=ST+RTL*DSIG                                             
               AMQ=AMQ+DSIG*REAL(Q(IL))                                   
               IL=IL+IGA                                                  
            ENDIF                                                         
            IG=IG+IGA-NWJ2                                                
   30    CONTINUE                                                         
         IOFS=NWJ2                                                        
  800 CONTINUE                                                            
      AMSP=1.0+REAL(SPA(1))*RSQR2                                         
      AMQ=CQ*RSQR2*AMQ                                                    
      TOTQ=CQ*TOTQ+AMQ                                                    
      PSITOT=SQRT(PSITOT)                                                 
      CHITOT=SQRT(CHITOT)                                                 
      TMPTOT=SQRT(TMPTOT+TOUT2+ST2B*SQR2)                                 
      TOTP=TOTP+RSQR2*REAL(GS(1))+(AMSP*TOUT1+TOTI+RSQR2*ST)/AKAP         
C                                                                         
      IF (KOUNT .EQ. 0) WRITE(2,40)                                       
        IF (KOUTP.EQ.0) THEN                                              
           WRITE (2,50) KOUNT,PSITOT,CHITOT,TMPTOT,TOTP,TOTQ,AMSP         
        ENDIF                                                             
   40 FORMAT(/3X,'KOUNT',3X,'RMSVORT',8X,'RMSDIV',8X,'RMSTEMP'            
     +       ,8X,'PE+IE',8X,'TOTQ',9X,'MSP')                              
   50 FORMAT(I6,1X,4E15.6,2F13.10)                                        
C                                                                         
C     Restore Z to absolute vorticity                                     
C                                                                         
      I=1                                                                 
      DO 102 L=1,NL                                                       
         Z(I)=Z(I)+EZ                                                     
         I=I+IGA                                                          
  102 CONTINUE                                                            
C                                                                         
      RETURN                                                              
      END                                                                 
