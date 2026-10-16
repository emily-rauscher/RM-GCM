C**********************************************************               
C             SUBROUTINE HEXP                                             
C**********************************************************               
      SUBROUTINE HEXP(SV,GV,NLS,ITYPE)                                    
C                                                                         
C     Performs an indirect Legendre transform for a (set of) field(s)     
C     having a total of NLS levels, from spectral to Fourier space.       
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
C     Version for RSGUP3.                     Mike Blackburn,  10.01.95.  
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
C     Legendre polynomials and information about gaussian latitudes       
C                                                                         
      COMMON/LEGAU/ ALPJ(MJP),DALPJ(MJP)                                  
     +              ,ALP(NWJ2,2,JGL),DALP(NWJ2,2,JGL)                     
     +              ,RLP(NWJ2,2,JGL),RDLP(NWJ2,2,JGL)                     
     +              ,SI(JGG),CS(JGG),SISQ(JGG),CSSQ(JGG),SECSQ(JGG)       
     +              ,ALAT(JGG),GWT(JGG),AW(JGG),JH,JL,JINC                
C                                                                         
      COMPLEX     SV(IGA,NLS),GV(IGL,NLS),TEMP                            
      REAL        ALPN(NWJ2,2,JGL,4)                                      
      EQUIVALENCE (ALPN(1,1,1,1),ALP(1,1,1))                              
C                                                                         
 6900 FORMAT(/' ***ABORT : HEXP CALLED WITH INVALID TYPE =',I5)           
C                                                                         
C     Use ITYPE to define transform type and symmetry labels.             
C     ISPAR is symmetry of spectral field    = 0 for D,T,SP etc.          
C                                            = 1 for Z.                   
C     IGPAR is symmetry of Fourier field: same as ISPAR unless transform  
C                                         involves a d/dy.                
C                                                                         
      IF (ITYPE.LE.0.OR.ITYPE.GE.9) THEN                                  
         WRITE(2,6900) ITYPE                                              
         CALL ABORT                                                       
      ENDIF                                                               
      IALP=(ITYPE+1)/2                                                    
      ISPAR=MOD(ITYPE,2)                                                  
      IGPAR=ISPAR                                                         
      IF (IALP.EQ.2.OR.IALP.EQ.4) IGPAR=1-ISPAR                           
C                                                                         
C     Perform inverse Legendre transform from spectral space to form
C     the even and odd contributions to the Fourier transforms, then
C     (global run) sum and difference them to give the complete
C     Fourier transforms at the northern and southern latitude rows:
C        IGPAR=0 : even (first half) precedes odd (second half).
C        IGPAR=1 : odd precedes even.
C     Fields L are independent, so they are shared among OpenMP
C     threads. Each GV(IM,L) still sums its N terms in the original
C     order, so the result does not depend on the thread count. Small
C     calls stay serial (IF clause) to avoid thread start-up costs.
C
      IA=ISPAR+1
      IB=2-ISPAR
!$OMP PARALLEL DO IF(NLS.GE.16) DEFAULT(SHARED)
!$OMP& PRIVATE(L,M,N,IM,IG,IP,K,KB,TEMP)
      DO 50 L=1,NLS
         IM=0
         IP=0
         IF (NHEM.EQ.1) THEN
            DO 10 M=0,MM-1,MOCT
               IM=IM+1
               DO 10 N=M,NN-1,2
                  IP=IP+1
                  GV(IM,L)=GV(IM,L)+ALPN(IP,IA,JL,IALP)*SV(IP,L)
   10       CONTINUE
         ELSE
            DO 20 M=0,MM-1,MOCT
               IM=IM+1
               IG=IM+IDL
               DO 20 N=M,NN-1,2
                  IP=IP+1
                  GV(IM,L)=GV(IM,L)+ALPN(IP,IA,JL,IALP)*SV(IP     ,L)
                  GV(IG,L)=GV(IG,L)+ALPN(IP,IB,JL,IALP)*SV(IP+NWJ2,L)
   20       CONTINUE
            IF (IGPAR.EQ.0) THEN
               DO 30 K=1,NWW
                  KB=K+IDL
                  TEMP=GV(K,L)
                  GV(K,L)=TEMP+GV(KB,L)
                  GV(KB,L)=TEMP-GV(KB,L)
   30          CONTINUE
            ELSE
               DO 40 K=1,NWW
                  KB=K+IDL
                  TEMP=GV(K,L)
                  GV(K,L)=GV(KB,L)+TEMP
                  GV(KB,L)=GV(KB,L)-TEMP
   40          CONTINUE
            ENDIF
         ENDIF
   50 CONTINUE
!$OMP END PARALLEL DO
C
      RETURN                                                              
      END                                                                 
