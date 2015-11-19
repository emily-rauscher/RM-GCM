C**********************************************************               
C             SUBROUTINE HEXPV                                            
C**********************************************************               
      SUBROUTINE HEXPV(Z,D,T,Q,SP,CHIG,SFG,UG,VG,ZG,DG,TG,QG,PLG,PJG)     
C                                                                         
C     Perform all the indirect Legendre transforms for the adiabatic      
C     part of the timestep at the current latitude (pair), in place       
C     of separate calls to HEXP for individual transform types.           
C                                                                         
C     The following types of Legendre function and thence types of        
C     transform may be used:                                              
C        ITYPE=1,2  :  ALP   :  ALPN(,,,1)   :  normal transform.         
C        ITYPE=3,4  :  DALP  :  ALPN(,,,2)   :  y-derivative.             
C        ITYPE=5,6  :  RLP   :  ALPN(,,,3)   :  del(-2).                  
C        ITYPE=7,8  :  RDLP  :  ALPN(,,,4)   :  y-derivative of del(-2).  
C     An even/odd value of ITYPE denotes a spectral field of even/odd     
C     symmetry.                                                           
C                                                                         
C     Maximum vector efficiency is achieved by chaining all multi-level   
C     transforms in one loop and by chaining all single-level transforms  
C     in a second loop.                                                   
C                                                                         
C     All dummy argument arrays are declared complex.                     
C     All array dimensions are parameters.                                
C     Multi-level arrays are 3-dimensional.                               
C                                                                         
C     Tracers (QNLG,QT) are included conditionally inside inner loops.    
C     This does *not* affect vectorisation, since the number of tracers,  
C     NTRAC, is a parameter and the relevant IF constructs are removed    
C     at compilation.  The same applies to code conditional on NHEM.      
C     (Only tested on Cray J90 with cf77).                                
C     Version for RSGUP3.                     Mike Blackburn,  12.01.95.  
C     Include tracers for IGCM2.              Piers Forster,   ??.04.97.  
C     Conditional separate loops for tracers. Mike Blackburn,  04.09.98.  
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
C     Legendre polynomials and information about gaussian latitudes       
C                                                                         
      COMMON/LEGAU/ ALPJ(MJP),DALPJ(MJP)                                  
     +              ,ALP(NWJ2,2,JGL),DALP(NWJ2,2,JGL)                     
     +              ,RLP(NWJ2,2,JGL),RDLP(NWJ2,2,JGL)                     
     +              ,SI(JGG),CS(JGG),SISQ(JGG),CSSQ(JGG),SECSQ(JGG)       
     +              ,ALAT(JGG),GWT(JGG),AW(JGG),JH,JL,JINC                
C                                                                         
C                                                                         
      COMPLEX     Z(NWJ2,NHEM,NL),D(NWJ2,NHEM,NL),T(NWJ2,NHEM,NL)         
     *           ,Q(NWJ2,NHEM,NL*NTRAC)                                   
     *           ,SP(NWJ2,NHEM)                                           
      COMPLEX     CHIG(IDL,NHEM,NL),SFG(IDL,NHEM,NL)                      
     *           ,UG(IDL,NHEM,NL),VG(IDL,NHEM,NL)                         
     *           ,ZG(IDL,NHEM,NL),DG(IDL,NHEM,NL),TG(IDL,NHEM,NL)         
     *           ,QG(IDL,NHEM,NL*NTRAC)                                   
     *           ,PLG(IDL,NHEM),PJG(IDL,NHEM)                             
      COMPLEX     TEMP                                                    
      REAL        ALPN(NWJ2,2,JGL,4)                                      
      EQUIVALENCE (ALPN(1,1,1,1),ALP(1,1,1))                              
C                                                                         
C     Transform multi-level fields to create even and odd contributions   
C     to the Fourier coefficients at the northern hemisphere latitude.    
C     Second subscript of ALPN denotes odd or even subset of Legendre     
C     functions, and depends of symmetry of spectral field.               
C     Fourth subscript of ALPN is Legendre function type, (ITYPE+1)/2.    
C                                                                         
      IM=0                                                                
      IP=0                                                                
      DO 10 M=0,MM-1,MOCT                                                 
         IM=IM+1                                                          
         DO 10 N=M,NN-1,2                                                 
            IP=IP+1                                                       
            DO 10 L=1,NL                                                  
C              Velocity potential    : type 6.                            
               CHIG(IM,1,L)=CHIG(IM,1,L) + ALPN(IP,1,JL,3)*D(IP,1,L)      
C              Streamfunction        : type 5.                            
               SFG (IM,1,L)=SFG (IM,1,L) + ALPN(IP,2,JL,3)*Z(IP,1,L)      
C              Zonal (rot) wind      : type 7.                            
               UG  (IM,1,L)=UG  (IM,1,L) + ALPN(IP,2,JL,4)*Z(IP,1,L)      
C              Merid (div) wind      : type 8.                            
               VG  (IM,1,L)=VG  (IM,1,L) + ALPN(IP,1,JL,4)*D(IP,1,L)      
C              (Relative) vorticity  : type 1.                            
               ZG  (IM,1,L)=ZG  (IM,1,L) + ALPN(IP,2,JL,1)*Z(IP,1,L)      
C              Divergence            : type 2.                            
               DG  (IM,1,L)=DG  (IM,1,L) + ALPN(IP,1,JL,1)*D(IP,1,L)      
C              Temperature           : type 2.                            
               TG  (IM,1,L)=TG  (IM,1,L) + ALPN(IP,1,JL,1)*T(IP,1,L)      
               IF (NTRAC.EQ.1) THEN                                       
C                 Tracers            : type 2.                            
                  QG (IM,1,L)=QG(IM,1,L) + ALPN(IP,1,JL,1)*Q(IP,1,L)      
               ENDIF                                                      
               IF (NHEM.EQ.2) THEN                                        
C                 Velocity potential    : type 6.                         
                  CHIG(IM,2,L)=CHIG(IM,2,L) + ALPN(IP,2,JL,3)*D(IP,2,L)   
C                 Streamfunction        : type 5.                         
                  SFG (IM,2,L)=SFG (IM,2,L) + ALPN(IP,1,JL,3)*Z(IP,2,L)   
C                 Zonal (rot) wind      : type 7.                         
                  UG  (IM,2,L)=UG  (IM,2,L) + ALPN(IP,1,JL,4)*Z(IP,2,L)   
C                 Merid (div) wind      : type 8.                         
                  VG  (IM,2,L)=VG  (IM,2,L) + ALPN(IP,2,JL,4)*D(IP,2,L)   
C                 (Relative) vorticity  : type 1.                         
                  ZG  (IM,2,L)=ZG  (IM,2,L) + ALPN(IP,1,JL,1)*Z(IP,2,L)   
C                 Divergence            : type 2.                         
                  DG  (IM,2,L)=DG  (IM,2,L) + ALPN(IP,2,JL,1)*D(IP,2,L)   
C                 Temperature           : type 2.                         
                  TG  (IM,2,L)=TG  (IM,2,L) + ALPN(IP,2,JL,1)*T(IP,2,L)   
                  IF (NTRAC.EQ.1) THEN                                    
C                    Tracers            : type 2.                         
                     QG (IM,2,L)=QG(IM,2,L) + ALPN(IP,2,JL,1)*Q(IP,2,L)   
                  ENDIF                                                   
               ENDIF                                                      
   10 CONTINUE                                                            
C                                                                         
C     Transform tracers, treated as a single (NL*NTRAC)-level field.      
C                                                                         
      IF (NTRAC.GT.1) THEN                                                
         IM=0                                                             
         IP=0                                                             
         DO 12 M=0,MM-1,MOCT                                              
            IM=IM+1                                                       
            DO 12 N=M,NN-1,2                                              
               IP=IP+1                                                    
               DO 12 L=1,NL*NTRAC                                         
C                 Tracers            : type 2.                            
                  QG(IM,1,L)=QG(IM,1,L)+ALPN(IP,1,JL,1)*Q(IP,1,L)         
                  IF (NHEM.EQ.2) THEN                                     
                     QG(IM,2,L)=QG(IM,2,L)+ALPN(IP,2,JL,1)*Q(IP,2,L)      
                  ENDIF                                                   
   12    CONTINUE                                                         
      ENDIF                                                               
C                                                                         
C     Transform single-level fields to create even and odd contributions  
C     to the Fourier coefficients at the northern hemisphere latitude.    
C     Vectorisation is over total wavenumber for each zonal wavenumber.   
C                                                                         
      IM=0                                                                
      IP=0                                                                
      DO 20 M=0,MM-1,MOCT                                                 
         IM=IM+1                                                          
         DO 20 N=M,NN-1,2                                                 
            IP=IP+1                                                       
C           Log (surface pressure)    : type 2.                           
            PLG(IM,1)=PLG(IM,1) + ALPN(IP,1,JL,1)*SP(IP,1)                
C           Merid gradient of ln (ps) : type 4.                           
            PJG(IM,1)=PJG(IM,1) + ALPN(IP,1,JL,2)*SP(IP,1)                
            IF (NHEM.EQ.2) THEN                                           
C              Log (surface pressure)    : type 2.                        
               PLG(IM,2)=PLG(IM,2) + ALPN(IP,2,JL,1)*SP(IP,2)             
C              Merid gradient of ln (ps) : type 4.                        
               PJG(IM,2)=PJG(IM,2) + ALPN(IP,2,JL,2)*SP(IP,2)             
            ENDIF                                                         
   20 CONTINUE                                                            
C                                                                         
C     For a global run, sum and difference even and odd contributions     
C     to give the complete Fourier transforms at the northern and         
C     southern latitude rows: N=E+O, S=E-O.                               
C     For symmetric Fourier fields, even (IM,1) precedes odd (IM,2).      
C     For asymmetric Fourier fields, odd (IM,1) precedes even (IM,2).     
C                                                                         
      IF (NHEM.EQ.2) THEN                                                 
C                                                                         
C        Multi-level fields.                                              
C                                                                         
         DO 30 IM=1,NWW                                                   
            DO 30 L=1,NL                                                  
C              Velocity potential : symmetric.                            
               TEMP=CHIG(IM,1,L)                                          
               CHIG(IM,1,L)=TEMP+CHIG(IM,2,L)                             
               CHIG(IM,2,L)=TEMP-CHIG(IM,2,L)                             
C              Streamfunction : anti-symmetric.                           
               TEMP=SFG(IM,1,L)                                           
               SFG(IM,1,L)=SFG(IM,2,L)+TEMP                               
               SFG(IM,2,L)=SFG(IM,2,L)-TEMP                               
C              Zonal (rotational) wind : symmetric.                       
               TEMP=UG(IM,1,L)                                            
               UG(IM,1,L)=TEMP+UG(IM,2,L)                                 
               UG(IM,2,L)=TEMP-UG(IM,2,L)                                 
C              Meridional (divergent) wind : anti-symmetric.              
               TEMP=VG(IM,1,L)                                            
               VG(IM,1,L)=VG(IM,2,L)+TEMP                                 
               VG(IM,2,L)=VG(IM,2,L)-TEMP                                 
C              Vorticity : anti-symmetric.                                
               TEMP=ZG(IM,1,L)                                            
               ZG(IM,1,L)=ZG(IM,2,L)+TEMP                                 
               ZG(IM,2,L)=ZG(IM,2,L)-TEMP                                 
C              Divergence : symmetric.                                    
               TEMP=DG(IM,1,L)                                            
               DG(IM,1,L)=TEMP+DG(IM,2,L)                                 
               DG(IM,2,L)=TEMP-DG(IM,2,L)                                 
C              Temperature : symmetric.                                   
               TEMP=TG(IM,1,L)                                            
               TG(IM,1,L)=TEMP+TG(IM,2,L)                                 
               TG(IM,2,L)=TEMP-TG(IM,2,L)                                 
               IF (NTRAC.EQ.1) THEN                                       
C                 Tracers : symmetric.                                    
                  TEMP=QG(IM,1,L)                                         
                  QG(IM,1,L)=TEMP+QG(IM,2,L)                              
                  QG(IM,2,L)=TEMP-QG(IM,2,L)                              
               ENDIF                                                      
   30    CONTINUE                                                         
C                                                                         
C        Tracers, treated as a single (NL*NTRAC)-level field.             
C                                                                         
         IF (NTRAC.GT.1) THEN                                             
            DO 32 IM=1,NWW                                                
               DO 32 L=1,NL*NTRAC                                         
C                 Tracers : symmetric.                                    
                  TEMP=QG(IM,1,L)                                         
                  QG(IM,1,L)=TEMP+QG(IM,2,L)                              
                  QG(IM,2,L)=TEMP-QG(IM,2,L)                              
   32       CONTINUE                                                      
         ENDIF                                                            
C                                                                         
C        Single-level fields.                                             
C                                                                         
         DO 40 IM=1,NWW                                                   
C           Log (surface pressure) : symmetric.                           
            TEMP=PLG(IM,1)                                                
            PLG(IM,1)=TEMP+PLG(IM,2)                                      
            PLG(IM,2)=TEMP-PLG(IM,2)                                      
C           Meridional gradient of ln(ps) : anti-symmetric.               
            TEMP=PJG(IM,1)                                                
            PJG(IM,1)=PJG(IM,2)+TEMP                                      
            PJG(IM,2)=PJG(IM,2)-TEMP                                      
   40    CONTINUE                                                         
C                                                                         
      ENDIF                                                               
C                                                                         
      RETURN                                                              
      END                                                                 
