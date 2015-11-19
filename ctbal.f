C**********************************************************               
C             SUBROUTINE TBAL                                             
C**********************************************************               
      SUBROUTINE TBAL                                                     
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
C     Constants and arrays needed for balancing                           
C                                                                         
      COMMON/BALAN/BFILT(NL),RGT0(NL),RG(NL2),TMEAN(NL)                   
     +            ,EP1(IGA),EP2(IGA),KBAL,MFTBAL,SRGT0,LTBAL              
      LOGICAL LTBAL                                                       
C                                                                         
      DIMENSION PMNRE(IDJ)                                                
      COMPLEX ERR(IDK),VPS,GSI1,TA,SRGT                                   
C                                                                         
C***********************************************************************  
C     Note that this scheme does not converge if wavenumbers M are        
C     present for which there are only a small number of modes in the     
C     truncation. For a small number of iterations the problem is not     
C     serious but it may be removed altogether by limiting the range of   
C     the 20 and 80 loops. For a given M, 7 modes are sufficient and 3    
C     are insufficient. The loop terminator MFTBAL is set in INITAL.      
C***********************************************************************  
C                                                                         
      REWIND 99                                                           
      I1=0                                                                
      DO 20 MP=1,MFTBAL,MOCT                                              
         J2=(NFP-MP)/2+1                                                  
         J22=J2*J2                                                        
         IJ=0                                                             
         DO 3 IN=MP,NFP,MH                                                
            I1=I1+1                                                       
            IF (IN.GT.1) THEN                                             
               VPS=VP(I1)                                                 
               K=I1                                                       
               GSI1=GS(I1)                                                
               IL=0                                                       
               DO 10 L=1,NL                                               
                  TA=(0.,0.)                                              
                  KK=I1                                                   
                  SRGT=0.                                                 
                  DO 9 M=1,NL                                             
                     IL=IL+1                                              
                     TA=TA+G(IL)*TT(KK)                                   
                     SRGT=SRGT+G(IL)*T(KK)                                
                     KK=KK+IGA                                            
    9             CONTINUE                                                
                  IJ=IJ+1                                                 
                  ERR(IJ)=-DT(K)-SQ(IN)*                                  
     +                    (SRGT+GSI1+T0(L)*SP(I1)+DELT*(TA-T0(L)*VPS))    
                  K=K+IGA                                                 
   10          CONTINUE                                                   
            ENDIF                                                         
    3    CONTINUE                                                         
         IF (MP.GT.1) THEN                                                
            READ(99) (PMNRE(I),I=1,J22)                                   
            IJ=0                                                          
            DO 24 J=1,J2                                                  
               IZJ=IZJ+1                                                  
               IZ=IZJ                                                     
               DO 23 L=1,NL                                               
                  IK=L                                                    
                  DO 22 K=1,J2                                            
                     IJ=IJ+1                                              
                     Z(IZ)=Z(IZ)+PMNRE(IJ)*ERR(IK)                        
                     IK=IK+NL                                             
   22             CONTINUE                                                
                  IZ=IZ+IGA                                               
                  IJ=IJ-J2                                                
   23          CONTINUE                                                   
               IJ=IJ+J2                                                   
   24       CONTINUE                                                      
         ELSE                                                             
            IJS=(J2-2)*NL                                                 
            IL=J2                                                         
            DO 28 L=1,NL                                                  
               IE=J2                                                      
               IZ=IL                                                      
               SRGT=0.                                                    
               IJ=IJS+L                                                   
               DO 26 J=2,J2                                               
                  SRGT=(ERR(IJ)-EP2(IE)*SRGT)/EP1(IE-1)                   
                  IJ=IJ-NL                                                
                  IZ=IZ-1                                                 
                  Z(IZ)=Z(IZ)+SRGT                                        
                  IE=IE-1                                                 
   26          CONTINUE                                                   
               IL=IL+IGA                                                  
   28       CONTINUE                                                      
            IZJ=J2                                                        
         ENDIF                                                            
   20 CONTINUE                                                            
C                                                                         
      IF (NHEM.EQ.2) THEN                                                 
         I1=NWJ2                                                          
         DO 80 MP=1,MFTBAL,MOCT                                           
            J2=(NFP-MP)/2+1                                               
            J2L=J2-1                                                      
            J22L=J2*J2L                                                   
            IJ=0                                                          
            DO 60 IN=MP,NFP,MH                                            
               I1=I1+1                                                    
               VPS=VP(I1)                                                 
               GSI1=GS(I1)                                                
               K=I1                                                       
               IL=0                                                       
               DO 58 L=1,NL                                               
                  TA=(0.0,0.0)                                            
                  SRGT=(0.0,0.0)                                          
                  KK=I1                                                   
                  DO 56 M=1,NL                                            
                     IL=IL+1                                              
                     TA=TA + G(IL)*TT(KK)                                 
                     SRGT=SRGT + G(IL)*T(KK)                              
                     KK=KK+IGA                                            
   56             CONTINUE                                                
                  IJ=IJ+1                                                 
                  ERR(IJ)=-DT(K)-SQ(IN+1)*(SRGT+GSI1+T0(L)*SP(I1)+        
     +                    DELT*(TA-T0(L)*VPS))                            
                  K=K+IGA                                                 
   58          CONTINUE                                                   
   60       CONTINUE                                                      
            IF (MP.EQ.1) THEN                                             
               READ(99) (PMNRE(I),I=1,J22L)                               
               IZJ=1+NWJ2                                                 
               IJ=0                                                       
               DO 64 J=1,J2L                                              
                  IZJ=IZJ+1                                               
                  IZ=IZJ                                                  
                  DO 63 L=1,NL                                            
                     IK=L                                                 
                     DO 62 K=1,J2                                         
                        IJ=IJ+1                                           
                        Z(IZ)=Z(IZ) + PMNRE(IJ)*ERR(IK)                   
                        IK=IK+NL                                          
   62                CONTINUE                                             
                     IZ=IZ+IGA                                            
                     IJ=IJ-J2                                             
   63             CONTINUE                                                
                  IJ=IJ+J2                                                
   64          CONTINUE                                                   
            ELSE                                                          
               IJS=(J2-1)*NL                                              
               IL=I1                                                      
               DO 78 L=1,NL                                               
                  IE=I1                                                   
                  IZ=IL                                                   
                  SRGT=(0.0,0.0)                                          
                  IJ=IJS+L                                                
                  DO 76 J=1,J2                                            
                     SRGT=(ERR(IJ)-EP2(IE)*SRGT)/EP1(IE)                  
                     Z(IZ)=Z(IZ)+SRGT                                     
                     IJ=IJ-NL                                             
                     IZ=IZ-1                                              
                     IE=IE-1                                              
   76             CONTINUE                                                
                  IL=IL+IGA                                               
   78          CONTINUE                                                   
            ENDIF                                                         
   80    CONTINUE                                                         
      ENDIF                                                               
C                                                                         
      IF(KOUNT.EQ.0)THEN                                                  
         DO 4 I=1,IGA                                                     
            SPMI(I)=SP(I)                                                 
    4    CONTINUE                                                         
         DO 8 J=1,IGB                                                     
            ZMI(J)=Z(J)                                                   
            DMI(J)=D(J)                                                   
            TMI(J)=T(J)                                                   
    8    CONTINUE                                                         
      ENDIF                                                               
C                                                                         
      RETURN                                                              
      END                                                                 
