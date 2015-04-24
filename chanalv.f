C**********************************************************               
C             SUBROUTINE HANALV                                           
C**********************************************************               
      SUBROUTINE HANALV(SPA,VP,DTE,TT,QT,DT,ZT,SPG,VPG,EG                 
     *                 ,TNLG,QNLG,FUG,FVG,VTG,VQG,FVGT,FUGT)              
C                                                                         
C     Perform all the direct Legendre transforms for the adiabatic        
C     part of the timestep at the current latitude (pair), in place       
C     of separate calls to HANAL for individual transform types.          
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
C                                                                         
C     NOTE: The y-derivative transforms use integration by parts and      
C           are valid only if the input field has zero zonal mean at      
C           both poles.                                                   
C                                                                         
C     NOTE: *** THE INPUT FOURIER FIELDS ARE MODIFIED IF GLOBAL ***       
C                                                                         
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
C     Polynomial used to aid vectorization of Legendre transforms         
C                                                                         
      COMMON/POLYNO/POLY(NWJ2,2,4),CMPA(IGL)                              
      COMPLEX CMPA                                                        
C                                                                         
C                                                                         
      COMPLEX     SPA(NWJ2,NHEM),VP(NWJ2,NHEM),DTE(NWJ2,NHEM,NL)          
     *           ,TT(NWJ2,NHEM,NL),DT(NWJ2,NHEM,NL),ZT(NWJ2,NHEM,NL)      
     *           ,QT(NWJ2,NHEM,NL*NTRAC)                                  
      COMPLEX     SPG(IDL,NHEM),VPG(IDL,NHEM),EG(IDL,NHEM,NL)             
     *           ,TNLG(IDL,NHEM,NL),FUG(IDL,NHEM,NL),FVG(IDL,NHEM,NL)     
     *           ,VTG(IDL,NHEM,NL),FVGT(IDL,NHEM,NL),FUGT(IDL,NHEM,NL)    
     *           ,QNLG(IDL,NHEM,NL*NTRAC),VQG(IDL,NHEM,NL*NTRAC)          
      COMPLEX     TEMP                                                    
      REAL        ALPN(NWJ2,2,JGL,4)                                      
      EQUIVALENCE (ALPN(1,1,1,1),ALP(1,1,1))                              
C                                                                         
C     For a global run, sum and difference the complete Fourier           
C     transforms at the northern and southern latitude rows to give       
C     the even and odd contributions : E=(N+S)/2, O=(N-S)/2.              
C     For Fourier fields symmetric about equator  : even precedes odd.    
C     For Fourier fields asymmetric about equator : odd precedes even.    
C                                                                         
      IF (NHEM.EQ.2) THEN                                                 
C                                                                         
C        Single-level fields.                                             
C                                                                         
         DO 10 IM=1,NWW                                                   
C           Surface pressure : symmetric.                                 
            TEMP=SPG(IM,1)                                                
            SPG(IM,1)=0.5*(TEMP+SPG(IM,2))                                
            SPG(IM,2)=0.5*(TEMP-SPG(IM,2))                                
C           Surface pressure tendency : symmetric.                        
            TEMP=VPG(IM,1)                                                
            VPG(IM,1)=0.5*(TEMP+VPG(IM,2))                                
            VPG(IM,2)=0.5*(TEMP-VPG(IM,2))                                
   10    CONTINUE                                                         
C                                                                         
C        Multi-level fields.                                              
C                                                                         
         DO 20 IM=1,NWW                                                   
            DO 20 L=1,NL                                                  
C              Divergence tendency : energy term : symmetric.             
               TEMP=EG(IM,1,L)                                            
               EG(IM,1,L)=0.5*(TEMP+EG(IM,2,L))                           
               EG(IM,2,L)=0.5*(TEMP-EG(IM,2,L))                           
C              Temperature tendency : main + d/dx part : symmetric.       
               TEMP=TNLG(IM,1,L)                                          
               TNLG(IM,1,L)=0.5*(TEMP+TNLG(IM,2,L))                       
               TNLG(IM,2,L)=0.5*(TEMP-TNLG(IM,2,L))                       
               IF (NTRAC.EQ.1) THEN                                       
C                 Tracer tendency : main + d/dx part : symmetric.         
                  TEMP=QNLG(IM,1,L)                                       
                  QNLG(IM,1,L)=0.5*(TEMP+QNLG(IM,2,L))                    
                  QNLG(IM,2,L)=0.5*(TEMP-QNLG(IM,2,L))                    
               ENDIF                                                      
C              Divergence tendency : d/dx part : symmetric.               
               TEMP=FUG(IM,1,L)                                           
               FUG(IM,1,L)=0.5*(TEMP+FUG(IM,2,L))                         
               FUG(IM,2,L)=0.5*(TEMP-FUG(IM,2,L))                         
C              Vorticity tendency : d/dx part : anti-symmetric.           
               TEMP=FVG(IM,1,L)                                           
               FVG(IM,1,L)=0.5*(TEMP-FVG(IM,2,L))                         
               FVG(IM,2,L)=0.5*(TEMP+FVG(IM,2,L))                         
C              Temperature tendency : d/dy part : anti-symmetric.         
               TEMP=VTG(IM,1,L)                                           
               VTG(IM,1,L)=0.5*(TEMP-VTG(IM,2,L))                         
               VTG(IM,2,L)=0.5*(TEMP+VTG(IM,2,L))                         
               IF (NTRAC.EQ.1) THEN                                       
C                 Tracer tendency : d/dy part : anti-symmetric.           
                  TEMP=VQG(IM,1,L)                                        
                  VQG(IM,1,L)=0.5*(TEMP-VQG(IM,2,L))                      
                  VQG(IM,2,L)=0.5*(TEMP+VQG(IM,2,L))                      
               ENDIF                                                      
C              Divergence tendency : d/dy part : anti-symmetric.          
               TEMP=FVGT(IM,1,L)                                          
               FVGT(IM,1,L)=0.5*(TEMP-FVGT(IM,2,L))                       
               FVGT(IM,2,L)=0.5*(TEMP+FVGT(IM,2,L))                       
C              Vorticity tendency : d/dy part : symmetric.                
               TEMP=FUGT(IM,1,L)                                          
               FUGT(IM,1,L)=0.5*(TEMP+FUGT(IM,2,L))                       
               FUGT(IM,2,L)=0.5*(TEMP-FUGT(IM,2,L))                       
   20    CONTINUE                                                         
C                                                                         
C        Tracers, treated as a single (NL*NTRAC)-level field.             
C                                                                         
         IF (NTRAC.GT.1) THEN                                             
            DO 22 IM=1,NWW                                                
               DO 22 L=1,NL*NTRAC                                         
C                 Tracer tendency : main + d/dx part : symmetric.         
                  TEMP=QNLG(IM,1,L)                                       
                  QNLG(IM,1,L)=0.5*(TEMP+QNLG(IM,2,L))                    
                  QNLG(IM,2,L)=0.5*(TEMP-QNLG(IM,2,L))                    
C                 Tracer tendency : d/dy part : anti-symmetric.           
                  TEMP=VQG(IM,1,L)                                        
                  VQG(IM,1,L)=0.5*(TEMP-VQG(IM,2,L))                      
                  VQG(IM,2,L)=0.5*(TEMP+VQG(IM,2,L))                      
   22       CONTINUE                                                      
         ENDIF                                                            
C                                                                         
      ENDIF                                                               
C                                                                         
C     Set up the appropriate Gaussian weight for the current latitude,    
C     dependent on transform type.                                        
C     Assumes JH in /LEGAU/ contains latitude counter from calling loop.  
C                                                                         
      AW1256=AW(JH)*CSSQ(JH)                                              
      AW3478=-AW(JH)                                                      
C                                                                         
C     Calculate POLY array in a vector loop before the main transforms,   
C     for the required Legendre Function types.                           
C     Both even and odd functions are required, irrespective of NHEM.     
C     Second subscript of ALPN denotes odd or even subset of Legendre     
C     functions, and depends of symmetry of spectral field.               
C     Fourth subscript of ALPN is Legendre function type, (ITYPE+1)/2.    
C                                                                         
      DO 30 IHEM=1,2                                                      
         DO 30 IP=1,NWJ2                                                  
            POLY(IP,IHEM,1)=AW1256*ALPN(IP,IHEM,JL,1)                     
            POLY(IP,IHEM,2)=AW3478*ALPN(IP,IHEM,JL,2)                     
   30 CONTINUE                                                            
C                                                                         
C     Transform single-level fields.                                      
C     Vectorisation is over total wavenumber for each zonal wavenumber.   
C                                                                         
      IM=0                                                                
      IP=0                                                                
      DO 40 M=0,MM-1,MOCT                                                 
         IM=IM+1                                                          
         DO 40 N=M,NN-1,2                                                 
            IP=IP+1                                                       
C           Surface pressure          : type 2.                           
            SPA(IP,1)=SPA(IP,1) + POLY(IP,1,1)*SPG(IM,1)                  
C           Surface pressure tendency : type 2.                           
            VP (IP,1)=VP (IP,1) + POLY(IP,1,1)*VPG(IM,1)                  
            IF (NHEM.EQ.2) THEN                                           
C              Surface pressure          : type 2.                        
               SPA(IP,2)=SPA(IP,2) + POLY(IP,2,1)*SPG(IM,2)               
C              Surface pressure tendency : type 2.                        
               VP (IP,2)=VP (IP,2) + POLY(IP,2,1)*VPG(IM,2)               
            ENDIF                                                         
   40 CONTINUE                                                            
C                                                                         
C     Transform multi-level fields.                                       
C     Inner loop vectorisation is over total wavenumber, to access        
C     spectral memory sequentially, avoiding skip distances being a       
C     multiple of 8 (which causes memory bank conflicts on Cray vector    
C     machines).                                                          
C                                                                         
      IM=0                                                                
      IP=0                                                                
      DO 50 M=0,MM-1,MOCT                                                 
         IM=IM+1                                                          
         IPM=IP                                                           
         DO 50 L=1,NL                                                     
            IP=IPM                                                        
            DO 50 N=M,NN-1,2                                              
               IP=IP+1                                                    
C              Divergence tendency  : energy term      : type 2.          
               DTE(IP,1,L)=DTE(IP,1,L) + POLY(IP,1,1)*EG  (IM,1,L)        
C              Temperature tendency : main + d/dx part : type 2.          
C              Temperature tendency : d/dy part        : type 4.          
               TT (IP,1,L)=TT (IP,1,L) + POLY(IP,1,1)*TNLG(IM,1,L)        
     *                                 + POLY(IP,1,2)*VTG (IM,1,L)        
               IF (NTRAC.EQ.1) THEN                                       
C                 Tracer tendency   : main + d/dx part : type 2.          
C                 Tracer tendency   : d/dy part        : type 4.          
                  QT(IP,1,L)=QT(IP,1,L) +POLY(IP,1,1)*QNLG(IM,1,L)        
     *                                  +POLY(IP,1,2)*VQG (IM,1,L)        
               ENDIF                                                      
C              Divergence tendency  : d/dx part        : type 2.          
C              Divergence tendency  : d/dy part        : type 4.          
               DT (IP,1,L)=DT (IP,1,L) + POLY(IP,1,1)*FUG (IM,1,L)        
     *                                 + POLY(IP,1,2)*FVGT(IM,1,L)        
C              Vorticity tendency   : d/dx part        : type 1.          
C              Vorticity tendency   : d/dy part        : type 3.          
               ZT (IP,1,L)=ZT (IP,1,L) + POLY(IP,2,1)*FVG (IM,1,L)        
     *                                 + POLY(IP,2,2)*FUGT(IM,1,L)        
               IF (NHEM.EQ.2) THEN                                        
C                 Divergence tendency  : energy term      : type 2.       
                  DTE(IP,2,L)=DTE(IP,2,L) + POLY(IP,2,1)*EG  (IM,2,L)     
C                 Temperature tendency : main + d/dx part : type 2.       
C                 Temperature tendency : d/dy part        : type 4.       
                  TT (IP,2,L)=TT (IP,2,L) + POLY(IP,2,1)*TNLG(IM,2,L)     
     *                                    + POLY(IP,2,2)*VTG (IM,2,L)     
                  IF (NTRAC.EQ.1) THEN                                    
C                    Tracer tendency   : d/dy part        : type 4.       
                     QT(IP,2,L)=QT(IP,2,L) +POLY(IP,2,1)*QNLG(IM,2,L)     
     *                                     +POLY(IP,2,2)*VQG (IM,2,L)     
                  ENDIF                                                   
C                 Divergence tendency  : d/dx part        : type 2.       
C                 Divergence tendency  : d/dy part        : type 4.       
                  DT (IP,2,L)=DT (IP,2,L) + POLY(IP,2,1)*FUG (IM,2,L)     
     *                                    + POLY(IP,2,2)*FVGT(IM,2,L)     
C                 Vorticity tendency   : d/dx part        : type 1.       
C                 Vorticity tendency   : d/dy part        : type 3.       
                  ZT (IP,2,L)=ZT (IP,2,L) + POLY(IP,1,1)*FVG (IM,2,L)     
     *                                    + POLY(IP,1,2)*FUGT(IM,2,L)     
               ENDIF                                                      
   50 CONTINUE                                                            
C                                                                         
C     Transform tracers, treated as a single (NL*NTRAC)-level field.      
C                                                                         
      IF (NTRAC.GT.1) THEN                                                
         IM=0                                                             
         IP=0                                                             
         DO 60 M=0,MM-1,MOCT                                              
            IM=IM+1                                                       
            IPM=IP                                                        
            DO 60 L=1,NL*NTRAC                                            
               IP=IPM                                                     
               DO 60 N=M,NN-1,2                                           
                  IP=IP+1                                                 
C                 Tracer tendency   : main + d/dx part : type 2.          
C                 Tracer tendency   : d/dy part        : type 4.          
                  QT(IP,1,L)=QT(IP,1,L)+POLY(IP,1,1)*QNLG(IM,1,L)         
     *                                 +POLY(IP,1,2)*VQG (IM,1,L)         
                  IF (NHEM.EQ.2) THEN                                     
                     QT(IP,2,L)=QT(IP,2,L)+POLY(IP,2,1)*QNLG(IM,2,L)      
     *                                    +POLY(IP,2,2)*VQG (IM,2,L)      
                  ENDIF                                                   
   60    CONTINUE                                                         
      ENDIF                                                               
C                                                                         
      RETURN                                                              
      END                                                                 
