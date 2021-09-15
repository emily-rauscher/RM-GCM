C**********************************************************               
C             SUBROUTINE MGRMLT                                           
C**********************************************************               
      SUBROUTINE MGRMLT                                                   
C                                                                         
C     Computes nonlinear tendencies in grid point space                   
C     for the present latitude                                            
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
      DIMENSION SDOTP(MG,NLM),SUMD(MG),TPTA(MG),TPTB(MG)                  
      DIMENSION CC(NL,NL)                                                 
      EQUIVALENCE (CC(1,1),C(1))                                          
C                                                                         
      IOFM=0                                                              
C                                                                         
C     Loop over hemispheres                                               
C                                                                         
      DO 800 IHEM=1,NHEM                                                  
         DO 100 I=1,MG                                                    
            J=I+IOFM                                                      
            SUMD(I)=0.0                                                   
            VPG(J)=0.0                                                    
C                                                                         
C           Change from Ln(PSTAR) to PSTAR                                
C                                                                         
            SPG(J)=EXP(PLG(J))-1.0                                        
C                                                                         
 100     CONTINUE                                                         
         DO 110 L=1,NLM                                                   
            DO 120 I=1,MG                                                 
              J=I+IOFM                                                    
              SUMD(I)=SUMD(I)+DSIGMA(L)*DG(J,L)                           
              VPG(J)=VPG(J)+DSIGMA(L)*SECSQ(JH)*(UG(J,L)*PMG(J)+          
     1               VG(J,L)*PJG(J))                                      
              SDOTP(I,L)=SUMD(I)+VPG(J)                                   
 120        CONTINUE                                                      
 110     CONTINUE                                                         
         DO 121 I=1,MG                                                    
            J=I+IOFM                                                      
            SUMD(I)=SUMD(I)+DSIGMA(NL)*DG(J,NL)                           
            VPG(J)=VPG(J)+DSIGMA(NL)*SECSQ(JH)*(UG(J,NL)*PMG(J)+          
     1             VG(J,NL)*PJG(J))                                       
 121     CONTINUE                                                         
         DO 130 L=1,NLM                                                   
            DO 140 I=1,MG                                                 
               J=I+IOFM                                                   
               SDOTP(I,L)=SIGMAH(L)*(SUMD(I)+VPG(J))-SDOTP(I,L)           
 140        CONTINUE                                                      
 130     CONTINUE                                                         
         DO 150 I=1,MG                                                    
            SUMD(I)=0.0                                                   
 150     CONTINUE                                                         
         DO 160 L=1,NL                                                    
            DO 170 I=1,MG                                                 
               TPTA(I)=0.0                                                
               TPTB(I)=0.0                                                
 170        CONTINUE                                                      
            DO 180 LL=1,L                                                 
               DO 190 I=1,MG                                              
                  J=I+IOFM                                                
                  VGPG=SECSQ(JH)*(UG(J,LL)*PMG(J)+VG(J,LL)*PJG(J))        
                  TPTA(I)=TPTA(I)+CC(LL,L)*VGPG                           
                  TPTB(I)=TPTB(I)+CC(LL,L)*(VGPG+DG(J,LL))                
 190           CONTINUE                                                   
 180        CONTINUE                                                      
            DO 200 I=1,MG                                                 
               J=I+IOFM                                                   
               UTG(J,L)=UG(J,L)*TG(J,L)                                   
               VTG(J,L)=VG(J,L)*TG(J,L)                                   
               EG(J,L)=UG(J,L)*UG(J,L)+VG(J,L)*VG(J,L)                    
 200        CONTINUE                                                      
            DO 202 KK=1,NTRAC                                             
               DO 201 I=1,MG                                              
                  J=I+IOFM                                                
                  UTRAG(J,L,KK)=UG(J,L)*TRAG(J,L,KK)                      
                  VTRAG(J,L,KK)=VG(J,L)*TRAG(J,L,KK)                      
  201          CONTINUE                                                   
  202       CONTINUE                                                      
            IF (L.GT.1.AND.L.LT.NL) THEN                                  
               DO 210 I=1,MG                                              
                  J=I+IOFM                                                
                  TSUM=SECSQ(JH)*(UG(J,L)*PMG(J)+VG(J,L)*PJG(J))          
                  SUMD(I)=SUMD(I)+TSUM*DSIGMA(L)                          
                  TNLG(J,L)=TG(J,L)*DG(J,L)+AKAP*TG(J,L)*(TSUM-TPTB(I))+  
     1                      TKP(L)*(TSUM-TPTA(I))-                        
     2                      RDSIG(L)*(SDOTP(I,L)*(TG(J,L+1)-TG(J,L))+     
     3                      SDOTP(I,L-1)*(TG(J,L)-                        
     4                      TG(J,L-1))+VPG(J)*(T01S2(L)*SIGMAH(L)+        
     5                      T01S2(L-1)*SIGMAH(L-1))-                      
     6                      SUMD(I)*(T01S2(L-1)+T01S2(L))+                
     7                      TSUM*DSIGMA(L)*T01S2(L-1))                    
                  FVG(J,L)=-UG(J,L)*ZG(J,L)-PJG(J)*TG(J,L)-               
     1                     RDSIG(L)*(SDOTP(I,L)*(VG(J,L+1)-VG(J,L))+      
     2                     SDOTP(I,L-1)*(VG(J,L)-VG(J,L-1)))              
                  FUG(J,L)=VG(J,L)*ZG(J,L)-PMG(J)*TG(J,L)-                
     1                     RDSIG(L)*(SDOTP(I,L)*(UG(J,L+1)-UG(J,L))+      
     2                     SDOTP(I,L-1)*(UG(J,L)-UG(J,L-1)))              
 210           CONTINUE                                                   
               DO 213 KK=1,NTRAC                                          
                  IF (.NOT.LFLUX) THEN                                    
                     DO 211 I=1,MG                                        
                        J=I+IOFM                                          
                        TRANLG(J,L,KK)=TRAG(J,L,KK)*DG(J,L)-RDSIG(L)*     
     1                    (SDOTP(I,L  )*(TRAG(J,L+1,KK)-TRAG(J,L  ,KK))   
     2                    +SDOTP(I,L-1)*(TRAG(J,L  ,KK)-TRAG(J,L-1,KK)))  
  211                CONTINUE                                             
                  ELSE                                                    
                     DO 212 I=1,MG                                        
                        J=I+IOFM                                          
                        TRANLG(J,L,KK)=RDSIG(L)*                          
     1                    (SDOTP(I,L  )*(TRAG(J,L+1,KK)+TRAG(J,L  ,KK))   
     2                    -SDOTP(I,L-1)*(TRAG(J,L  ,KK)+TRAG(J,L-1,KK)))  
  212                CONTINUE                                             
                  ENDIF                                                   
  213          CONTINUE                                                   
            ELSE                                                          
               FAC=1.0                                                    
               FAC2=1.                                                    
               K=L                                                        
               IF (L.EQ.NL) THEN                                          
                  K=L-1                                                   
                  FAC=0.0                                                 
                  FAC2=-1.                                                
               ENDIF                                                      
               DO 220 I=1,MG                                              
                  J=I+IOFM                                                
                  TSUM=SECSQ(JH)*(UG(J,L)*PMG(J)+VG(J,L)*PJG(J))          
                  SUMD(I)=SUMD(I)+TSUM*DSIGMA(L)*FAC                      
                  TNLG(J,L)=TG(J,L)*DG(J,L)+AKAP*TG(J,L)*(TSUM-TPTB(I))+  
     1                      TKP(L)*(TSUM-TPTA(I))-                        
     2                      RDSIG(L)*(SDOTP(I,K)*(TG(J,K+1)-TG(J,K))+     
     3                      T01S2(K)*(SIGMAH(K)*VPG(J)-SUMD(I)))          
                  FVG(J,L)=-UG(J,L)*ZG(J,L)-PJG(J)*TG(J,L)-               
     1                     RDSIG(L)*(SDOTP(I,K)*(VG(J,K+1)-VG(J,K)))      
                  FUG(J,L)=VG(J,L)*ZG(J,L)-PMG(J)*TG(J,L)-                
     1                     RDSIG(L)*(SDOTP(I,K)*(UG(J,K+1)-UG(J,K)))      
 220           CONTINUE                                                   
               DO 223 KK=1,NTRAC                                          
                  IF (.NOT.LFLUX) THEN                                    
                     DO 221 I=1,MG                                        
                        J=I+IOFM                                          
                        TRANLG(J,L,KK)=TRAG(J,L,KK)*DG(J,L)-RDSIG(L)*     
     1                         SDOTP(I,K)*(TRAG(J,K+1,KK)-TRAG(J,K,KK))   
  221                CONTINUE                                             
                  ELSE                                                    
                     DO 222 I=1,MG                                        
                        J=I+IOFM                                          
                        TRANLG(J,L,KK)=RDSIG(L)*FAC2*                     
     1                         SDOTP(I,K)*(TRAG(J,K+1,KK)+TRAG(J,K,KK))   
  222                CONTINUE                                             
                  ENDIF                                                   
  223          CONTINUE                                                   
            ENDIF                                                         
 160     CONTINUE                                                         
         IOFM=MGPP                                                        
 800  CONTINUE                                                            
      RETURN                                                              
      END                                                                 
