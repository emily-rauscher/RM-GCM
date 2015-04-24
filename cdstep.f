C**********************************************************               
C             SUBROUTINE DSTEP                                            
C**********************************************************               
      SUBROUTINE DSTEP                                                    
C                                                                         
C     Diabatic part of timestep. Completion of time-filter.               
C     Note that only Z,D,T have diabatic tendencies at present.           
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
      IF(KOUNT.GT.KITS) THEN                                              
C                                                                         
C       Ordinary centred timestep                                         
C                                                                         
        DO 10 I=1,IGB                                                     
           Z(I)=Z(I)+DELT2*ZT(I)                                          
           T(I)=T(I)+DELT2*TT(I)                                          
           D(I)=D(I)+DELT2*DT(I)                                          
           ZMI(I)=ZMI(I)+PNU*Z(I)                                         
           TMI(I)=TMI(I)+PNU*T(I)                                         
           DMI(I)=DMI(I)+PNU*D(I)                                         
   10   CONTINUE                                                          
        DO 11 KK=1,NTRAC                                                  
           DO 11 I=1,IGB                                                  
              TRA(I,KK)=TRA(I,KK)+DELT2*TRAT(I,KK)                        
              TRAMI(I,KK)=TRAMI(I,KK)+PNU*TRA(I,KK)                       
   11   CONTINUE                                                          
        DO 20 I=1,IGA                                                     
           SPMI(I)=SPMI(I)+PNU*SP(I)                                      
20      CONTINUE                                                          
        RETURN                                                            
      ELSE                                                                
C                                                                         
C       Short initial timestep                                            
C                                                                         
        DO 60 I=1,IGB                                                     
           Z(I)=Z(I)+DELT2*ZT(I)                                          
           T(I)=T(I)+DELT2*TT(I)                                          
           D(I)=D(I)+DELT2*DT(I)                                          
   60   CONTINUE                                                          
        DO 61 KK=1,NTRAC                                                  
           DO 61 I=1,IGB                                                  
             TRA(I,KK)=TRA(I,KK)+DELT2*TRAT(I,KK)                         
   61   CONTINUE                                                          
        DELT=DELT2                                                        
        DELT2=DELT2+DELT2                                                 
        RETURN                                                            
      ENDIF                                                               
C                                                                         
      END                                                                 
