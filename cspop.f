C**********************************************************               
C             SUBROUTINE SPOP                                             
C**********************************************************               
      SUBROUTINE SPOP                                                     
C                                                                         
C     Controls diagnostic output from model run.                          
C     Outputs spectral coefficients.                                      
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
C     Array ordering in GRIDP must correspond to that in SPECTR.          
C     Real arrays: multi-level arrays are 1-dimensional.                  
C                                                                         
      COMMON/GRIDP/ CHIG(IGD),SFG(IGD),UG(IGD),VG(IGD)                    
     *              ,ZG(IGD),DG(IGD),TG(IGD)                              
     +              ,TRAG(IGD,NTRAC)                                      
     *              ,PLG(IGC),PJG(IGC),PMG(IGC)                           
     *              ,SPG(IGC),VPG(IGC),EG(IGD)                            
     +              ,TNLG(IGD),TRANLG(IGD,NTRAC),FUG(IGD),FVG(IGD)        
     +              ,UTG(IGD),UTRAG(IGD,NTRAC)                            
     +              ,VTG(IGD),VTRAG(IGD,NTRAC),FVGT(IGD),FUGT(IGD)        
     $              ,GRPAD(NGRPAD)                                        
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
C     Legendre polynomials and information about gaussian latitudes       
C                                                                         
      COMMON/LEGAU/ ALPJ(MJP),DALPJ(MJP)                                  
     +              ,ALP(NWJ2,2,JGL),DALP(NWJ2,2,JGL)                     
     +              ,RLP(NWJ2,2,JGL),RDLP(NWJ2,2,JGL)                     
     +              ,SI(JGG),CS(JGG),SISQ(JGG),CSSQ(JGG),SECSQ(JGG)       
     +              ,ALAT(JGG),GWT(JGG),AW(JGG),JH,JL,JINC                
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
 200  FORMAT(/' NUMBER OF TIME STEPS COMPLETED =',I5)                     
 202  FORMAT(' SPECTRAL COEFFICIENTS (COEFF ; AMPLITUDE ; PHASE)')        
 204  FORMAT(' VORTICITY AT LEVEL',I2)                                    
 206  FORMAT(' DIVERGENCE AT LEVEL',I2)                                   
 208  FORMAT(' PERTURBATION TEMPERATURE AT LEVEL',I2)                     
 211  FORMAT(' LOG(SURFACE PRESSURE)')                                    
 215  FORMAT(' TRACER FIELD ',I4,' AT LEVEL',I2)                          
C                                                                         
      IF (NCOEFF.EQ.0) RETURN                                             
C                                                                         
C     Spectral coeficients are wanted                                     
C                                                                         
      WRITE (2,200) KOUNT                                                 
      WRITE (2,202)                                                       
C                                                                         
C     Absolute vorticity                                                  
C                                                                         
      DO 18 L=1,NL                                                        
         IF (LSPO(L)) THEN                                                
            WRITE (2,204) L                                               
            CALL WRSPA(Z(1+(L-1)*IGA),1)                                  
         ENDIF                                                            
 18   CONTINUE                                                            
C                                                                         
C     Divergence                                                          
C                                                                         
      DO 28 L=1,NL                                                        
         IF (LSPO(L)) THEN                                                
            WRITE (2,206) L                                               
            CALL WRSPA(D(1+(L-1)*IGA),2)                                  
         ENDIF                                                            
 28   CONTINUE                                                            
C                                                                         
C     Temperature                                                         
C                                                                         
      DO 38 L=1,NL                                                        
         IF (LSPO(L)) THEN                                                
            WRITE (2,208) L                                               
            CALL WRSPA(T(1+(L-1)*IGA),2)                                  
         ENDIF                                                            
 38   CONTINUE                                                            
C                                                                         
C     Tracers                                                             
C                                                                         
      DO 48 KK=1,NTRAC                                                    
         DO 49 L=1,NL                                                     
            IF(LSPO(L)) THEN                                              
               WRITE(2,215) KK,L                                          
               CALL WRSPA(TRA(1+(L-1)*IGA,KK),2)                          
            ENDIF                                                         
   49    CONTINUE                                                         
   48 CONTINUE                                                            
C                                                                         
C     Log (Surface Pressure)                                              
C                                                                         
      WRITE(2,211)                                                        
      CALL WRSPA(SP(1),2)                                                 
C                                                                         
      RETURN                                                              
      END                                                                 
