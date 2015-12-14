C**************************************************************           
C                    SUBROUTINE INIGAU                                    
C**************************************************************           
      SUBROUTINE INIGAU                                                   
C                                                                         
C     This subroutine calculates gaussian weights and latitudes           
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
C     Constant arrays and variables associated with time and vertical     
C     differencing. Also counters.                                        
C                                                                         
      COMMON/BATS/  BEGDAY,CTRA(NTRAC),BM1(IDE),AK(NNP),AQ(NL2),G(NL2)              
     +              ,TAU(NL2),KOUNT,KITS,KSTART,KTOTAL,KRUN,ITSPD           
     +              ,DELT,DELT2,CV,CG,CT,CQ,PNU,PNU2,PNU21                
     +              ,NTRACO,KOLOUR(NTRAC),RGG(NL2)            
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
  230 FORMAT(/' *** ABORT *** PARAMETER JGL IS DIFFERENT FROM 1 OR JG')   
  202 FORMAT(' GAUSSIAN LATITUDES')                                       
  201 FORMAT(10F7.2)                                                      
C                                                                         
      IF(JGL.NE.1.AND.JGL.NE.JG) THEN                                     
         WRITE(2,230)                                                     
         CALL ABORT                                                       
      ENDIF                                                               
      IF(JGL.EQ.1) JINC=0                                                 
      IF(JGL.EQ.JG)JINC=1                                                 
      JL=1                                                                
      DO 8 J=1,JG                                                         
         JH=J                                                             
         CALL GWTLT(SI(J),WEIGHT,J,JG)                                    
         SISQ(J)=SI(J)*SI(J)                                              
         CSSQ(J)=1.-SISQ(J)                                               
         SECSQ(J)=1./CSSQ(J)                                              
         CS(J)=SQRT(CSSQ(J))                                              
         ALAT(J)=ATAN(SI(J)/CS(J))*180.0/PI                               
         GWT(J)=WEIGHT/REAL(NHEM)                                         
         AW(J)=WEIGHT*2.0*SECSQ(J)                                        
C                                                                         
C        Compute Legendre functions at the current latitude.              
C                                                                         
         CALL LGNDRE(NN,MM,MOCT,ALPJ,DALPJ,MJP,1,SI(JH),CS(JH))           
C                                                                         
C        Reorder Legendre functions, separating even/odd functions.       
C                                                                         
         DO 58 K=1,2                                                      
            I=0                                                           
            II=K-2                                                        
            DO 57 M=0,MM-1,MOCT                                           
               DO 56 N=M,NN-1,2                                           
                  I=I+1                                                   
                  II=II+2                                                 
                  ALP(I,K,JL)=ALPJ(II)                                    
                  DALP(I,K,JL)=DALPJ(II)                                  
                  RLP(I,K,JL)=-RSQ(N+K)*ALP(I,K,JL)                       
                  RDLP(I,K,JL)=-RSQ(N+K)*DALP(I,K,JL)                     
   56          CONTINUE                                                   
   57       CONTINUE                                                      
   58    CONTINUE                                                         
C                                                                         
         IF(JGL.EQ.1) WRITE(25) ALP,DALP,RLP,RDLP                         
C                                                                         
         JL=JL+JINC                                                       
    8 CONTINUE                                                            
C                                                                         
      IF (NHEM.EQ.2) THEN                                                 
CDIR$    IVDEP                                                            
         DO 59 J=1,JG                                                     
            SI(JGGP-J)=-SI(J)                                             
            CS(JGGP-J)=CS(J)                                              
            SISQ(JGGP-J)=SISQ(J)                                          
            CSSQ(JGGP-J)=CSSQ(J)                                          
            SECSQ(JGGP-J)=SECSQ(J)                                        
            ALAT(JGGP-J)=-ALAT(J)                                         
            GWT(JGGP-J)=GWT(J)                                            
            AW(JGGP-J)=AW(J)                                              
   59    CONTINUE                                                         
      ENDIF                                                               
C                                                                         
C     Output the Gaussian latitudes                                       
C                                                                         
      WRITE(2,202)                                                        
      WRITE(2,201)(ALAT(J),J=1,JGG)                                       
C                                                                         
      END                                                                 
