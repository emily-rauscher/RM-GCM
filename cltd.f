C**********************************************************               
C             SUBROUTINE LTD                                              
C**********************************************************               
      SUBROUTINE LTD                                                      
C                                                                         
C     Direct Legendre transform for the adiabatic part of the timestep.   
C     Transforms from Fourier to spectral space at the current latitude   
C     (pair).  In a global run the input arrays are complete (even+odd)   
C     Fourier coefficients at the northern & southern hemisphere rows.    
C                                                                         
C     Includes the option either to call the modular routine HANAL for    
C     each field to be transformed, or to call the fast vectorising       
C     routine HANALV to perform all transforms together.  The choice      
C     is controlled by logical LTVEC.                                     
C                                                                         
C     Each call to HANAL transforms fields having the same symmetry       
C     and type of Legendre Function.  HANAL1 is a separate routine        
C     with improved efficiency for single-level transforms.               
C                                                                         
C     The Fourier work array passed to HANAL must be dimensioned with     
C     (at least) the maximum number of levels used in the HANAL calls.    
C                                                                         
C     Version for RSGUP3.                     Mike Blackburn,  12.01.95.  
C     ANSI work arrays for HANAL,SPDEL2.      Mike Blackburn,  04.09.96.  
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
      REAL    FILT(0:NN)                                                  
      COMPLEX GWORK(IGL,(3+NTRAC)*NL+2)                                   
C                                                                         
C     Prepare Fourier arrays:                                             
C     - change sign of terms which contribute negatively to tendency,     
C     - apply (1-mu**2) weighting,                                        
C     - take zonal derivatives,                                           
C     - make copies of effective momentum tendencies.                     
C                                                                         
      DO 10 L=1,NL                                                        
         DO 10 I=1,IGL                                                    
            EG(I,L)=-0.5*EG(I,L)*SECSQ(JH)                                
            VTG(I,L)=-VTG(I,L)                                            
            TNLG(I,L)=TNLG(I,L)-CMPA(I)*UTG(I,L)*SECSQ(JH)                
            FVGT(I,L)=FVG(I,L)                                            
            FUGT(I,L)=-FUG(I,L)                                           
            FUG(I,L)=CMPA(I)*FUG(I,L)*SECSQ(JH)                           
            FVG(I,L)=CMPA(I)*FVG(I,L)*SECSQ(JH)                           
   10 CONTINUE                                                            
C                                                                         
C     Tracer transforms depend on the form of advection used.             
C                                                                         
      IF (LFLUX) THEN                                                     
         AFAC=-1.0                                                        
      ELSE                                                                
         AFAC=1.0                                                         
      ENDIF                                                               
      DO 11 KK=1,NTRAC                                                    
         DO 11 L=1,NL                                                     
            DO 11 I=1,IGL                                                 
               TRANLG(I,L,KK)=AFAC*TRANLG(I,L,KK)                         
     *                        -CMPA(I)*UTRAG(I,L,KK)*SECSQ(JH)            
               VTRAG(I,L,KK)=-VTRAG(I,L,KK)                               
   11 CONTINUE                                                            
C                                                                         
      IF (LTVEC) THEN                                                     
C                                                                         
C        Call single routine to perform all transforms with maximum       
C        vector efficiency.                                               
C        *** NOTE : THE INPUT FOURIER FIELDS ARE MODIFIED IF GLOBAL ***   
C                                                                         
         CALL HANALV(SPA,VP,DTE,TT,TRAT,DT,ZT,SPG,VPG,EG                  
     +              ,TNLG,TRANLG,FUG,FVG,VTG,VTRAG,FVGT,FUGT)             
C                                                                         
      ELSE                                                                
C                                                                         
C        Main transform of even fields:                                   
C        SPG to FUG (and SPA to DT) must be contiguous in common.         
C                                                                         
         CALL HANAL(SPG,GWORK,SPA,(3+NTRAC)*NL+2,2)                       
C                                                                         
C        Remaining transforms: VTG,FVGT (and TT,DT) must be contiguous.   
C                                                                         
         CALL HANAL(FVG,GWORK,ZT,NL,1)                                    
         CALL HANAL(VTG,GWORK,TT,(2+NTRAC)*NL,4)                          
         CALL HANAL(FUGT,GWORK,ZT,NL,3)                                   
C                                                                         
      ENDIF                                                               
C                                                                         
C     At the last latitude, take del**2 of the energy term and add to     
C     the divergence tendency.                                            
C                                                                         
      IF (JH.EQ.JG) THEN                                                  
         CALL SPDEL2(DTE,FILT,NWJ2,NN,MM,MOCT,NHEM,NL,0,2)                
         DO 20 I=1,IGB                                                    
            DT(I)=DT(I)+DTE(I)                                            
   20    CONTINUE                                                         
      ENDIF                                                               
C                                                                         
      RETURN                                                              
      END                                                                 
