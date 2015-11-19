C**********************************************************               
C             SUBROUTINE BALANC                                           
C**********************************************************               
      SUBROUTINE BALANC                                                   
C                                                                         
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
C     Constants and arrays needed for balancing                           
C                                                                         
      COMMON/BALAN/BFILT(NL),RGT0(NL),RG(NL2),TMEAN(NL)                   
     +            ,EP1(IGA),EP2(IGA),KBAL,MFTBAL,SRGT0,LTBAL              
      LOGICAL LTBAL                                                       
C                                                                         
      DIMENSION GP(NL)                                                    
      COMPLEX TA,GP,GSI1,VPS,TK,SRGT                                      
C                                                                         
C**********************************************************************   
C     2-grid vertical temperature wave is not removed if orography        
C     is included. Program aborts in INITAL if attempted.                 
C**********************************************************************   
C                                                                         
      I1=0                                                                
      DO 800 IHEM=1,NHEM                                                  
         DO 3 MP=1,MFP,MOCT                                               
            DO 4 IN=MP,NFP,MH                                             
               I1=I1+1                                                    
               INR=IN+IHEM-1                                              
               IF (INR.GT.1) THEN                                         
                  VPS=VP(I1)                                              
                  K=I1                                                    
                  GSI1=GS(I1)                                             
                  IL=0                                                    
                  DO 10 L=1,NL                                            
                     TA=(0.,0.)                                           
                     KK=I1                                                
                     DO 9 M=1,NL                                          
                        IL=IL+1                                           
                        TA=TA+G(IL)*TT(KK)                                
                        KK=KK+IGA                                         
    9                CONTINUE                                             
                     TA=(T0(L)*VPS-TA)*DELT - RSQ(INR)*DT(K)              
                     GP(L)=TA-GSI1                                        
                     K=K+IGA                                              
   10             CONTINUE                                                
                  IL=0                                                    
                  K=I1                                                    
                  SRGT=(0.,0.)                                            
                  DO 12 L=1,NL                                            
                     TK=(0.,0.)                                           
                     DO 11 M=1,NL                                         
                        IL=IL+1                                           
                        TK=TK+RG(IL)*GP(M)                                
   11                CONTINUE                                             
                     T(K)=TK                                              
                     SRGT=SRGT+BFILT(L)*TK                                
                     K=K+IGA                                              
   12             CONTINUE                                                
                  SRGT=SRGT/SRGT0                                         
                  SP(I1)=SRGT                                             
                  K=I1                                                    
                  DO 13 L=1,NL                                            
                     T(K)=T(K)-RGT0(L)*SRGT                               
                     K=K+IGA                                              
   13             CONTINUE                                                
               ENDIF                                                      
    4       CONTINUE                                                      
    3    CONTINUE                                                         
         I1=NWJ2                                                          
  800 CONTINUE                                                            
C                                                                         
      IF (KOUNT.EQ.0) THEN                                                
         DO 2 I=1,IGA                                                     
            SPMI(I)=SP(I)                                                 
    2    CONTINUE                                                         
         DO 5 J=1,IGB                                                     
            ZMI(J)=Z(J)                                                   
            DMI(J)=D(J)                                                   
            TMI(J)=T(J)                                                   
    5    CONTINUE                                                         
      ENDIF                                                               
C                                                                         
      RETURN                                                              
      END                                                                 
