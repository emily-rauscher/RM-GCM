C**********************************************************               
C             SUBROUTINE CBADJ                                            
C**********************************************************               
      SUBROUTINE CBADJ(NCRB,NCRT,J,IHEM)                                  
C                                                                         
C     PRECIPITATING DEEP CONVECTION SCHEME.                               
C     ADJUSTMENT TOWARDS A SUBSATURATED MOIST ADIABAT.                    
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
C     Array ordering in GRIDP must correspond to that in SPECTR.          
C     Real arrays: multi-level arrays are 2-dimensional.                  
C     the variables have been renamed to coincide with                    
C     variable names in bgcm5 DGRMLT                                      
C                                                                         
C swapped round VLNG and UNLG to check                                    
      COMMON/GRIDP/ CHIG(IGC,NL),SFG(IGC,NL),UG(IGC,NL),VG(IGC,NL)        
     :              ,TTVD(IGC,NL),QTVD(IGC,NL),TG(IGC,NL)                 
     :              ,TRAG(IGC,NL,NTRAC)                                   
     :              ,PLG(IGC),TYBL(IGC),TXBL(IGC)                         
     :              ,SPG(IGC),VPG(IGC),TTRD(IGC,NL)                       
     :              ,TNLG(IGC,NL),TRANLG(IGC,NL,NTRAC),UNLG(IGC,NL)       
     :              ,VNLG(IGC,NL),TTLR(IGC,NL),UTRAG(IGC,NL,NTRAC)        
     :              ,TTCR(IGC,NL),VTRAG(IGC,NL,NTRAC)                     
     :              ,UTVD(IGC,NL),VTVD(IGC,NL)                            
     :         ,ASSBL(IGC),ASHBL(IGC),ASLBL(IGC),ARRCR(IGC),ARRLR(IGC)    
     :         ,arflux(igc,6),asfld(igc,6),acld(igc,4)                    
     :         ,SSBL(IGC),SHBL(IGC),SLBL(IGC),RRCR(IGC),RRLR(IGC)         
     :         ,rflux(igc,6),sfld(igc,6),cld(igc,4)                       
C                                                                         
      COMMON/PHYS/  LBL,LVD,LCR,LLR,LRD,LCUBM,LCBADJ                      
     :              ,TSTAR(IGC,JG),QSTAR(IGC,JG),FRAD(JG,NHEM)            
     :              ,TSTARO(IGC,JG),TDEEPO(IGC,JG),smstar(igc,jg)         
     :              ,tdeep(igc,jg),hsnow(igc,jg),sqstar(igc,jg)           
     :              ,SALB(IGC,JG),SBAL(IGC,JG),BLCD(IGC)                  
     :              ,SVEGE(IGC,JG),CD,DRAG,BLVAD,BLA,BLRH,BLVB(IGC)       
     :              ,AKVV,AKTV,AKQV,ESCONA,ESCONB,EPSIQ,CTQ,CCC           
     : ,ctqi,sdsn,shcs,shcsp,shcsn,skse,sksn,slhf,sd1,sd2,sdw             
     :        ,ssmc,sdsnd,LSL,sasnow,saice,shsstar,shsmax                 
     :     ,LOC,SHCO,SHCI,LNOICE,LOLDBL,LCOND,LNNSK,ITSLL,ITSLO           
     :              ,CCR,RCON,DTBUOY,TSLA,TSLB,TSLC,TSLD,CUT1,CUT2        
     :              ,NLCR,NCUTOP,CURHM,AKTC,AKQC,CUBMT,CBADJT,CBADJP      
     :              ,SKAP(NL),SK(NLM),FWS(NL),CLR(NL),FB(NLM)             
     :              ,TTDC(NL),QTDC(NL),TTMC(NL),QTMC(NL),TC(NL),QC(NL)    
     :              ,CTCR(NL,NHEM),CTLR(NL,NHEM),NAVRD,NAVWT,DELT2C       
      LOGICAL LBL,LVD,LCR,LLR,LRD,LCUBM,LCBADJ,LSL,LOC                    
     :       ,LNOICE,LOLDBL,LCOND,LNNSK                                   
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
       COMMON/CPIERS/ICFLAG(IGC,5,2),CFRAC(IGC,5),PNET(IGC,JG)            
     :               ,SNET(IGC,JG),RRFLUX(IGC,JG,6)                       
     :               ,TTSW(IGC,NL),TTLW(IGC,NL)                           
      REAL PL(NL),PSL(NL),PSK(NL),TSL(NL),TR(NL),QR(NL),CBC(NL),GH1(NL)   
C                                                                         
C     Setup moisture variables by equivilencing them to                   
C     Tracer No. 1                                                        
C                                                                         
      REAL QG(IGC,NL),QNLG(IGC,NL),QTLR(IGC,NL),QTCR(IGC,NL)              
      EQUIVALENCE (QG(1,1),TRAG(1,1,1)) , (QNLG(1,1),TRANLG(1,1,1))       
     *            ,(QTLR(1,1),UTRAG(1,1,1)),(QTCR(1,1),VTRAG(1,1,1))      
C                                                                         
C*      2001 FORMAT(' CBADJ DEEP SWAP: I,J,NCRB,NCRT,QRINT,QGINT = '      
C*     F      ,4I5,2F10.4)                                                
C                                                                         
C     Choose appropriate value of latent heat based on ground T           
c     ie, we assume it snows if the ground is colder than 273.15K         
      if (tstar(j,jh).gt.0.363224029) then                                
         ctquse=ctq                                                       
      else                                                                
         ctquse=ctqi                                                      
      endif                                                               
      NLEV=1+NCRB-NCRT                                                    
      EPSIT=CTQUSE*EPSIQ                                                  
C                                                                         
      DO 10 L=NCRT,NCRB                                                   
        PL(L)=SIGMA(L)*PLG(J)                                             
        PSL(L)=PL(L)+CBADJP                                               
        PSK(L)=(PSL(L)/PL(L))**AKAP                                       
        CBC(L)=CTQUSE*PSK(L)/PSL(L)                                       
        TSL(L)=TC(L)*PSK(L)                                               
        QR(L)=PQSAT(TSL(L))/PSL(L)                                        
   10 CONTINUE                                                            
C                                                                         
      SINT=SSUM(NLEV,DSIGMA(NCRT),1)                                      
      TGINT=SDOT(NLEV,TG(J,NCRT),IGC,DSIGMA(NCRT),1)                      
      QGINT=SDOT(NLEV,QG(J,NCRT),IGC,DSIGMA(NCRT),1)                      
      TRINT=SDOT(NLEV,TC(NCRT),1,DSIGMA(NCRT),1)                          
      QRINT=SDOT(NLEV,QR(NCRT),1,DSIGMA(NCRT),1)                          
      HGINT=TGINT+CTQUSE*QGINT                                            
      DHL=(HGINT-TRINT-CTQUSE*QRINT)/SINT                                 
      DO 15 L=NCRT,NCRB                                                   
 15   GH1(L)=(TC(L)+CTQUSE*QR(L)+DHL)*PSK(L)                              
C                                                                         
      NIT=0                                                               
   20 NIT=NIT+1                                                           
      DTMAX=0.0                                                           
      DO 30 L=NCRT,NCRB                                                   
        QN=CBC(L)*PQSAT(TSL(L))                                           
        DTSL=(TSL(L)+QN-GH1(L))/(1.0+ESCONB*QN/(TSL(L)*TSL(L)))           
        DTMAX=MAX(ABS(DTSL),DTMAX)                                        
        TSL(L)=TSL(L)-DTSL                                                
   30 CONTINUE                                                            
      IF(DTMAX.LT.EPSIT) GOTO 40                                          
      IF(NIT.GT.10) GOTO 40                                               
      GOTO 20                                                             
   40 CONTINUE                                                            
C                                                                         
      DO 50 L=NCRT,NCRB                                                   
        TR(L)=TSL(L)/PSK(L)                                               
        QR(L)=PQSAT(TSL(L))/PSL(L)                                        
   50 CONTINUE                                                            
C                                                                         
      QRINT=SDOT(NLEV,QR(NCRT),1,DSIGMA(NCRT),1)                          
      IF(QRINT.GT.QGINT) THEN                                             
C*      WRITE(2,2001) J,JH,NCRB,NCRT,CQ*QRINT/SINT,CQ*QGINT/SINT          
        IF(LCUBM)      CALL CUBM (NCRB,NCRT,J,IHEM)                       
        IF(.NOT.LCUBM) CALL CUDIF(NCRB,NCRT,J,IHEM)                       
      ELSE                                                                
         write(*,*) 'MOIST CONVECTION using CBADJT!!!'
        TRINT=SDOT(NLEV,TR(NCRT),1,DSIGMA(NCRT),1)                        
        DHL=(HGINT-TRINT-CTQUSE*QRINT)/SINT                               
        DO 60 L=NCRT,NCRB                                                 
          TR(L)=TR(L)+DHL                                                 
          TTEND=(TR(L)-TG(J,L))/CBADJT                                    
          QTEND=(QR(L)-QG(J,L))/CBADJT                                    
          TTMC(L)=TTMC(L)+TTEND                                           
          QTMC(L)=QTMC(L)+QTEND                                           
          TG(J,L)=TG(J,L)+DELT2C*TTEND                                    
          QG(J,L)=QG(J,L)+DELT2C*QTEND                                    
          CTCR(L,IHEM)=CTCR(L,IHEM)+1.0                                   
   60   CONTINUE                                                          
        CTCR(NCRB,IHEM)=CTCR(NCRB,IHEM)-1.0                               
        RRCR(J)=RRCR(J)-RCON*PLG(J)*(QRINT-QGINT)/CBADJT                  
C deep convection parameterised by precip.:  P is mm/day                  
        P=-REAL(ITSPD)*RCON*PLG(J)*(QRINT-QGINT)/CBADJT                   
        IF (P.GT.0.15) THEN                                               
          CFRAC(J,5)=0.245+0.125*LOG(P)                                   
          CFRAC(J,5)=MIN(1.0,CFRAC(J,5))                                  
          ICFLAG(J,5,1)=NCRB                                              
          ICFLAG(J,5,2)=NCRT                                              
        ENDIF                                                             
      ENDIF                                                               
C                                                                         
      RETURN                                                              
      END                                                                 
