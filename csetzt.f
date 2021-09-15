C**********************************************************               
C             SUBROUTINE SETZT                                            
C**********************************************************               
      SUBROUTINE SETZT                                                    
C                                                                         
C     This subroutine sets up restoration temperature field.              
C     The temperature at SIGMA = 1 is TGR, entered in Kelvin.             
C     a lapse rate of ALR k/m is assumed under the tropopause and         
C     zero above. The actual profile tends to this away from the          
C     tropopause, with smooth interpolation depending on DTTRP            
C     at the model tropopause. The height of                              
C     the tropopause is given as ZTROP m.                                 
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
C     Restoration temperature field and constants which determine it,     
C     also contains timescales                                            
C                                                                         
      COMMON/RESTIJ/TTRES(IGB)                                            
     + ,DTNS,DTEP,DTTRP,FAC(NL),DDAMP(NL),TFRC(NL),YRLEN,TRS(NL)          
     +  ,RESTTT(NL),REDTEP(NL)
     +  ,ALR,ZTROP,TGR                                                    
      COMPLEX TTRES                                                       
C                    
      COMMON/NEWFORC/TTRESN(IGB)
      COMPLEX TTRESN
CC                                                                         
C     Constant arrays and variables associated with time and vertical     
C     differencing. Also counters.                                        
C                                                                         
      COMMON/BATS/  BEGDAY,CTRA(NTRAC),BM1(IDE),AK(NNP),AQ(NL2),G(NL2)              
     +              ,TAU(NL2),KOUNT,KITS,KSTART,KTOTAL,KRUN,ITSPD           
     +              ,DELT,DELT2,CV,CG,CT,CQ,PNU,PNU2,PNU21                
     +              ,NTRACO,KOLOUR(NTRAC),RGG(NL2)            
     +              ,BEGDOY,DOY                                           
C                                                                         
      DO 100 I=1,IGB                                                      
         TTRES(I)=(0.0,0.0)                                               
100   CONTINUE                                                            
      DTTRP=DTTRP*CT                                                      
      SIGPREV=1.                                                          
      TPREV=TGR                                                           
      ZPREV=0.                                                            
C ER Modif: can use RESTTT to define TRS instead of using default calculation below
      IF (RESTTT(1) .EQ. 0.) THEN
         DO 150 L=NL,1,-1                                                    
            ZP=ZPREV+(GASCON*TPREV/GA)*LOG(SIGPREV/SIGMA(L))                 
            TP=TGR-ZTROP*ALR                                                 
            TP=TP+SQRT((.5*ALR*(ZP-ZTROP))**2+DTTRP**2)                      
            TP=TP-.5*ALR*(ZP-ZTROP)                                          
            TPM=.5*(TPREV+TP)                                                
            ZPP=ZPREV+(GASCON*TPM/GA)*LOG(SIGPREV/SIGMA(L))                  
            TPP=TGR-ZTROP*ALR                                                
            TPP=TPP+SQRT((.5*ALR*(ZPP-ZTROP))**2+DTTRP**2)                   
            TPP=TPP-.5*ALR*(ZPP-ZTROP)                                       
            TRS(L)=TPP                                                       
            ZPREV=ZPREV+(.5*(TPP+TPREV)*GASCON/GA)*LOG(SIGPREV/SIGMA(L))     
            TPREV=TPP                                                        
            SIGPREV=SIGMA(L)                                                 
 150     CONTINUE                                                            
      ELSE
         DO 160 L=1,NL
            TRS(L)=RESTTT(L)
 160     CONTINUE
      ENDIF

C                                                                         
      WRITE(2,2000)                                                       
      WRITE(2,2010) TRS                                                   
 2000 FORMAT(/' RESTORATION TEMPERATURE STRATIFICATION IN K ')            
 2010 FORMAT(10F7.2)                                                      
C                                                                         
      DO 170 L=1,NL                                                       
         TRS(L)=TRS(L)/CT                                                 
 170  CONTINUE                                                            
C                                                                         
C     Now the latitudinal variation in TTRES is set up                    
C     (this being in terms of a deviation from T0 which                   
C     is usually constant with height)                                    
C                                                                         
C      DO 200 L=1,NL                                                       
C         I=(L-1)*IGA                                                      
C         TTRES(I+1)=SQRT(2.)*(TRS(L)-T0(L))                               
C         TTRES(I+2)=-2./3.*SQRT(0.4)*DTEP*FAC(L)                          
C         TTRES(I+NWJ2+1)=(1./SQRT(6.))*DTNS*FAC(L)                        
C200   CONTINUE                                                            
      DTTRP=DTTRP/CT                                                      
C                                                                         
C KM Modif to use New Forcing
      TTRES=TTRESN
C
      END                                                                 
