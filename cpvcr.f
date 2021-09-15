C**********************************************************               
C             SUBROUTINE PVCR                                             
C**********************************************************               
      SUBROUTINE PVCR(TFIELD,TTYPE)                                       
C                                                                         
C     Calculate potential temperature and Ertel potential vorticity       
C     on model levels.                                                    
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
      COMMON/COMGRM/DUDLSG(IGC,NL),DVDLSG(IGC,NL),DTDLSG(IGC,NL)          
C                                                                         
C                                                                         
C     Array ordering in GRIDP must correspond to that in SPECTR.          
C     Real arrays: multi-level arrays are 2-dimensional.                  
C                                                                         
      COMMON/GRIDP/ CHIG(IGC,NL),SFG(IGC,NL),UG(IGC,NL),VG(IGC,NL)        
     :              ,ZG(IGC,NL),DG(IGC,NL),TG(IGC,NL)                     
     :              ,TRAG(IGC,NL,NTRAC)                                   
     :              ,PLG(IGC),PJG(IGC),PMG(IGC)                           
     :              ,SPG(IGC),VPG(IGC),EG(IGC,NL)                         
     :              ,TNLG(IGC,NL),TRANLG(IGC,NL,NTRAC),FUG(IGC,NL)        
     :              ,FVG(IGC,NL),UTG(IGC,NL),UTRAG(IGC,NL,NTRAC)          
     :              ,VTG(IGC,NL),VTRAG(IGC,NL,NTRAC)                      
     :              ,FVGT(IGC,NL),FUGT(IGC,NL)                            
     :              ,GRPAD(NGRPAD)                                        
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
      REAL PSG(IGC,NL),TFIELD(IGC,NL)                                     
      INTEGER TTYPE                                                       
C                                                                         
C     Calculate pressure on model levels.                                 
C                                                                         
      DO 10 L=1,NL                                                        
         DO 11 II=1,IGC                                                   
            PSG(II,L)=SIGMA(L)*PLG(II)                                    
 11      CONTINUE                                                         
 10   CONTINUE                                                            
C                                                                         
      DO 100 IHEM=1,NHEM                                                  
         IOF=(IHEM-1)*MGPP                                                
C                                                                         
C     Calculate potential temperature and PV on model levels.             
C     Generalised formula using array of model level pressures is used    
C     to calculate potential temperature.  But note that the reference    
C     pressure for adiabatic processes is assumed to be P0, the non-      
C     dimensionalising pressure.                                          
C                                                                         
C     TFIELD is the tracer field to be initialised.                       
C     TTYPE=1   : Initialise tracer as potential temperature.             
C     TTYPE=2   : Initialise tracer as Ertel PV.                          
C                                                                         
C     Remember that at this point:-                                       
C     CHIG contains the x-derivative of T.                                
C     SFG  contains the y-derivative of T.                                
C                                                                         
         DO 20 L=1,NL                                                     
            DO 30 I=1,MG                                                  
               II=I+IOF                                                   
               PMKG=PSG(II,L)**(-AKAP)                                    
               PM1KG=PMKG/PSG(II,L)                                       
               AKT=AKAP*TG(II,L)                                          
               TXP=CHIG(II,L)-AKT*PMG(II)                                 
               TYP=SFG(II,L)-AKT*PJG(II)                                  
               IF(TTYPE.EQ.1) THEN                                        
                  TFIELD(II,L)=TG(II,L)*PMKG                              
               ELSE                                                       
                  TFIELD(II,L)=PM1KG*(ZG(II,L)*(-DTDLSG(II,L)+AKT)        
     :                 +SECSQ(JH)*(-DUDLSG(II,L)*TYP+DVDLSG(II,L)*TXP))   
               ENDIF                                                      
 30         CONTINUE                                                      
 20      CONTINUE                                                         
 100  CONTINUE                                                            
C                                                                         
      RETURN                                                              
      END                                                                 
