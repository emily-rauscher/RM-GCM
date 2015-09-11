C**********************************************************               
C             SUBROUTINE LTI                                              
C**********************************************************               
      SUBROUTINE LTI                                                      
C                                                                         
C     Inverse Legendre transform for the adiabatic part of the timestep.  
C     Transforms from spectral to Fourier space at the current latitude   
C     (pair).  In a global run the resulting arrays are complete          
C     (i.e. even+odd) Fourier coefficients at the northern & southern     
C     hemisphere rows.                                                    
C                                                                         
C     Includes the option either to call the modular routine HEXP for     
C     each field to be transformed, or to call the fast vectorising       
C     routine HEXPV to perform all transforms together.  The choice       
C     is controlled by logical LTVEC.                                     
C                                                                         
C     Each call to HEXP transforms fields having the same symmetry        
C     and type of Legendre Function.  HEXP1 is a separate routine         
C     with improved efficiency for single-level transforms.               
C                                                                         
C     Version for RSGUP3.                     Mike Blackburn,  12.01.95.  
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
      COMPLEX CEZN                                                        
C                                                                         
C     Preset Fourier arrays.                                              
C                                                                         
      DO 10 L=1,NL                                                        
         DO 10 I=1,IGL                                                    
            CHIG(I,L)=0.                                                  
            SFG(I,L)=0.                                                   
            UG(I,L)=0.                                                    
            VG(I,L)=0.                                                    
            ZG(I,L)=0.                                                    
            DG(I,L)=0.                                                    
            TG(I,L)=0.                                                    
   10 CONTINUE                                                            
C                                                                         
      DO 11 KK=1,NTRAC                                                    
         DO 11 L=1,NL                                                     
            DO 11 I=1,IGL                                                 
               TRAG(I,L,KK)=0.                                            
   11 CONTINUE                                                            
C                                                                         
      DO 20 I=1,IGL                                                       
         PLG(I)=0.                                                        
         PJG(I)=0.                                                        
   20 CONTINUE                                                            
C                                                                         
C     Remove planetary vorticity in spectral space, so all transforms     
C     use relative vorticity.                                             
C                                                                         
      DO 30 I=1,IGB,IGA                                                   
         Z(I)=Z(I)-EZ                                                     
   30 CONTINUE                                                            
C                                                                         
      IF (LTVEC) THEN                                                     
C                                                                         
C        Call single routine to perform all transforms with maximum       
C        vector efficiency.                                               
C                                                                         
         CALL HEXPV(Z,D,T,TRA,SP,CHIG,SFG,UG,VG,ZG,DG,TG,TRAG,PLG,PJG)    
C                                                                         
      ELSE                                                                
C                                                                         
C        Transform prognostic fields & meridional derivative of ln(ps).   
C        D,T,TRA,SP (and DG,TG,TRAG,PLG) must be contiguous in common.    
C                                                                         
         CALL HEXP(Z , ZG,  NL  ,1)                                       
         CALL HEXP(D,DG,(2+NTRAC)*NL+1,2)                                 
         CALL HEXP1(SP,PJG,    1,4)                                       
C                                                                         
C        Wind components: calls to HEXP give following Fourier fields:    
C           SFG  :   streamfunction.                                      
C           CHIG :   velocity potential.                                  
C           UG   :   -U(rotational).                                      
C           VG   :   V(divergent).                                        
C                                                                         
         CALL HEXP(Z,SFG ,NL,5)                                           
         CALL HEXP(D,CHIG,NL,6)                                           
         CALL HEXP(Z,UG  ,NL,7)                                           
         CALL HEXP(D,VG  ,NL,8)                                           
C                                                                         
      ENDIF                                                               
C                                                                         
C     Restore planetary vorticity in spectral space.                      
C                                                                         
      DO 40 I=1,IGB,IGA                                                   
         Z(I)=Z(I)+EZ                                                     
   40 CONTINUE                                                            
C                                                                         
C     Convert from relative to absolute vorticity in Fourier space:       
C     (real part of) m=0 coefficient only.                                
C                                                                         
      CEZN=CMPLX(2.0*SI(JH),0.0)                                          
      DO 50 L=1,NL                                                        
         ZG(1,L)=ZG(1,L)+CEZN                                             
         IF (NHEM.EQ.2) ZG(1+IDL,L)=ZG(1+IDL,L)-CEZN                      
   50 CONTINUE                                                            
C                                                                         
C     Sum to give total winds.  CMPA takes x-derivative.                  
C                                                                         
      DO 60 L=1,NL                                                        
         DO 60 I=1,IGL                                                    
            UG(I,L)=CMPA(I)*CHIG(I,L)-UG(I,L)                             
            VG(I,L)=CMPA(I)* SFG(I,L)+VG(I,L)                             
   60 CONTINUE                                                            
C                                                                         
C     Zonal gradient of ln(ps).                                           
C                                                                         
      DO 70 I=1,IGL                                                       
         PMG(I)=CMPA(I)*PLG(I)                                            
   70 CONTINUE                                                            
C                                                                         
      RETURN
      write(*,*) 'ctli line 207, plg:',PLG                                                              
      END                                                                 
