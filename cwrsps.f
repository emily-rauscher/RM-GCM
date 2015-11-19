C**********************************************************               
C             SUBROUTINE WRSPS                                            
C**********************************************************               
      SUBROUTINE WRSPS(A,IA)                                              
C                                                                         
C     Prints spectral coefficients                                        
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
      PARAMETER(RAD=180./PI)                                              
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
      COMPLEX A(NWJ2)                                                     
      CHARACTER COEFF(NWJ2,2)*8                                           
      REAL FP(IGA)                                                        
      INTEGER IP(IGA)                                                     
      SAVE COEFF,FP,IP                                                    
      COMPLEX POLAR,Z                                                     
      POLAR(Z)=CMPLX(ABS(Z),ATAN2(AIMAG(Z),REAL(Z)+1.0E-20)*RAD)          
      DATA COEFF/MJP*' (  ,  )'/                                          
      IG=0                                                                
      I=0                                                                 
      DO 200 MP=1,MFP,MOCT                                                
C                                                                         
C        Only print if MP<NCOEFF                                          
C                                                                         
         IF (MP.GT.NCOEFF) RETURN                                         
C                                                                         
         IBEG=IG                                                          
         DO 100 JP=MP,NCOEFF,MH                                           
            IG=IG+1                                                       
            I=I+1                                                         
            IP(I)=IG                                                      
            IF (MP.EQ.1) THEN                                             
               FP(I)=1.0                                                  
            ELSE                                                          
               FP(I)=2.0                                                  
            ENDIF                                                         
            WRITE(COEFF(I,1)(3:4),'(I2)')MP-1                             
            WRITE(COEFF(I,2)(3:4),'(I2)')MP-1                             
            WRITE(COEFF(I,1)(6:7),'(I2)')JP                               
            WRITE(COEFF(I,2)(6:7),'(I2)')JP-1                             
 100     CONTINUE                                                         
         IG=IBEG+(NFP-MP+2)/MH                                            
 200  CONTINUE                                                            
      RETURN                                                              
C******************************                                           
      ENTRY WRSPA(A,IA)                                                   
      IF (NHEM.EQ.1) THEN                                                 
         WRITE(2,1000)(COEFF(I,IA),POLAR(A(IP(I))*FP(I)),I=1,INSPC)       
      ELSE                                                                
         WRITE(2,1000)(COEFF(I,2),POLAR(A(IP(I)+(2-IA)*NWJ2)*FP(I)),      
     1   COEFF(I,1),POLAR(A(IP(I)+(IA-1)*NWJ2)*FP(I)),I=1,INSPC)          
      ENDIF                                                               
 1000 FORMAT(3(A8,1X,1PE9.2,0PF7.1))                                      
      RETURN                                                              
      END                                                                 
