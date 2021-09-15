C**********************************************************               
C             SUBROUTINE SETRES                                           
C**********************************************************               
      SUBROUTINE SETRES                                                   
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
C     Restoration fields and timescale                                    
C                                                                         
      COMMON/RESTOR/ZRES(IGN),DRES(IGN),TRES(IGN),SPRES(IGM),DAMP         
C                                                                         
C                                                                         
C     Set up restoration state from the KOUNT=0 zonally averaged          
C     state and write this to FT13 for future use.                        
C     This is only done when KOUNT=0 and DAMP.GT.0.0.                     
C                                                                         
 2200 FORMAT(/' RESTORATION RECORD WRITTEN TO CHANNEL ',I3)               
C                                                                         
      IF (DAMP.LE.0) RETURN                                               
C                                                                         
      DO 840 IHEM=1,NHEM                                                  
         I=NWJ2*(IHEM-1)                                                  
         IR=IDM*(IHEM-1)                                                  
         DO 100 J=1,IDM                                                   
            I=I+1                                                         
            IR=IR+1                                                       
            SPRES(IR)=SP(I)                                               
  100    CONTINUE                                                         
         DO 850 L=1,NL                                                    
            I=NWJ2*(IHEM-1)+(L-1)*IGA                                     
            IR=IDM*(IHEM-1)+(L-1)*IGM                                     
            DO 860 J=1,IDM                                                
               I=I+1                                                      
               IR=IR+1                                                    
               ZRES(IR)=Z(I)                                              
               DRES(IR)=D(I)                                              
               TRES(IR)=T(I)                                              
  860       CONTINUE                                                      
  850    CONTINUE                                                         
  840 CONTINUE                                                            
C                                                                         
      REWIND 13                                                           
      WRITE(13)ZRES,DRES,TRES,SPRES                                       
      WRITE(2,2200)13                                                     
C                                                                         
      END                                                                 
