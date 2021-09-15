C**********************************************************               
C             SUBROUTINE LTIDT                                            
C**********************************************************               
      SUBROUTINE LTIDT                                                    
C                                                                         
C     Inverse Legendre transform, from spectral to Fourier space          
C     for the derivatives of temperature required to calculate PV.        
C     HEXP transforms fields having the same symmetry and type of         
C     Legendre function.                                                  
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
C     Array ordering in GRIDP must correspond to that in SPECTR.          
C     Complex arrays: multi-level arrays are 2-dimensional.               
C                                                                         
      COMMON/GRIDP/ CHIG(IGL,NL),SFG(IGL,NL),UG(IGL,NL),VG(IGL,NL)        
     :              ,ZG(IGL,NL),DG(IGL,NL),TG(IGL,NL)                     
     :              ,TRAG(IGL,NL,NTRAC)                                   
     :              ,PLG(IGL),PJG(IGL),PMG(IGL)                           
     :              ,SPG(IGL),VPG(IGL),EG(IGL,NL)                         
     :              ,TNLG(IGL,NL),TRANLG(IGL,NL,NTRAC),FUG(IGL,NL)        
     :              ,FVG(IGL,NL),UTG(IGL,NL),UTRAG(IGL,NL,NTRAC)          
     :              ,VTG(IGL,NL),VTRAG(IGL,NL,NTRAC)                      
     :              ,FVGT(IGL,NL),FUGT(IGL,NL)                            
     :              ,GRPAD(NGRPAD)                                        
      COMPLEX CHIG,SFG,UG,VG,ZG,DG,TG,PLG,PJG,PMG                         
     :       ,SPG,VPG,EG,TNLG,FUG,FVG,UTG,VTG,FVGT,FUGT                   
     :       ,TRAG,TRANLG,UTRAG,VTRAG                                     
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
C     Polynomial used to aid vectorization of Legendre transforms         
C                                                                         
      COMMON/POLYNO/POLY(NWJ2,2,4),CMPA(IGL)                              
      COMPLEX CMPA                                                        
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
C     Calls to HEXPFL give following Fourier fields :                     
C        CHIG :   x-derivative of temperature field                       
C        SFG  :   y-derivative of temperature field                       
C     These are only the temporary contents for PV calculation.           
C                                                                         
      DO 10 L=1,NL                                                        
         DO 20 I=1,IGL                                                    
            CHIG(I,L)=0.                                                  
            SFG(I,L)=0.                                                   
 20      CONTINUE                                                         
 10   CONTINUE                                                            
      CALL HEXP(T,CHIG,NL,2)                                              
      CALL HEXP(T,SFG,NL,4)                                               
      DO 30 L=1,NL                                                        
         DO 40 I=1,IGL                                                    
            CHIG(I,L)=CMPA(I)*CHIG(I,L)                                   
 40      CONTINUE                                                         
 30   CONTINUE                                                            
C                                                                         
      RETURN                                                              
      END                                                                 
