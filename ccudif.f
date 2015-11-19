C**********************************************************               
C             SUBROUTINE CUDIF                                            
C**********************************************************               
      SUBROUTINE CUDIF(NCRB,NCT,J,IHEM)                                   
C                                                                         
C     NON-PRECIPITATING CONVECTION. DIFFUSION TOWARDS MIXING-LINE THETA   
C     AND CONSTANT Q. CLOUD-TOP FLUXES AND VARIABLE K'S WHEN SHALLOW.     
C     OPTIONAL RAINOUT OVER A SPECIFIED REL-HUM.                          
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
C                                                                         
      REAL PL(NL),PFK(NL),TR(NL)                                          
C                                                                         
C     Setup moisture variables by equivilencing them to                   
C     Tracer No. 1                                                        
C                                                                         
      REAL QG(IGC,NL),QNLG(IGC,NL),QTLR(IGC,NL),QTCR(IGC,NL)              
      EQUIVALENCE (QG(1,1),TRAG(1,1,1)) , (QNLG(1,1),TRANLG(1,1,1))       
     *            ,(QTLR(1,1),UTRAG(1,1,1)),(QTCR(1,1),VTRAG(1,1,1))      
C                                                                         
C                                                                         
 2001 FORMAT(' CUDIF QMIN APPLIED AT BAS: I,J,NCRB,NCRT,QG =',5I5,F10.4)  
 2002 FORMAT(' CUDIF QMIN APPLIED FOR QM: I,J,NCRB,NCRT,QM =',4I5,F10.4)  
 2004 FORMAT(' CUDIF MID LEVEL SWAP IGNORED: I,J,NCRB,NCRT =',4I5)        
 2103 FORMAT(' CUBM LAPSE RATE LIMIT IMPOSED: I,J,NCRB,NCRT = ',4I5)      
 2104 FORMAT(' CUBM LAPSE RATE LIMIT IMPOSED: I,J,NCRB,NCRT,TLIM          
     +                                                     = ',4I5,F7.4)  
C                                                                         
                                                                          
C                                                                         
      FTSL(T,Q,P)=TSLA+TSLB/(TSLC*LOG(T)-LOG(Q*P)+TSLD)                   
C                                                                         
C     Minimum moisture for mixing line calculations.                      
      QGMIN=1.E-6                                                         
C Stable limit to mixing line slope d(theta)/dp=-10^-3 K/Pa               
      TLIM=-1.0E-3*P0/CT                                                  
C                                                                         
      ESCON=1./PLG(J)                                                     
      IF(NCT.GE.NCUTOP) THEN                                              
        IF(BLVB(J).EQ.0.0)RETURN                                          
        NCRT=NCT-1                                                        
        AKTCU=BLVB(J)*BLCD(J)*CUT2/CUT1                                   
      ELSE IF (NCRB.GE.NCUTOP) THEN                                       
C        Limit the convective top for swap points to the maximum          
C        height of shallow convection.                                    
C        Note that NCRT includes the level above cloud top.               
         NCRT=NCUTOP-1                                                    
         AKTCU=AKTC                                                       
      ELSE                                                                
C        No adjustment for mid-level swap points.                         
C         WRITE(2,2004) J,JH,NCRB,NCT                                     
         RETURN                                                           
      ENDIF                                                               
      AKQCU=AKTCU*AKQC/AKTC                                               
      NLEV=1+NCRB-NCRT                                                    
      NCRTP=NCRT+1                                                        
      NCRBM=NCRB-1                                                        
C                                                                         
      RKAP=1.0/AKAP                                                       
      DO 50 L=NCRT,NCRB                                                   
        PL(L)=SIGMA(L)*PLG(J)                                             
        PFK(L)=PL(L)**AKAP                                                
   50 CONTINUE                                                            
      THB=TG(J,NCRB)/PFK(NCRB)                                            
      IF (QG(J,NCRB).LT.QGMIN) THEN                                       
         WRITE(2,2001) J,JH,L,NCRB,NCRT,CQ*QG(J,NCRB)                     
      ENDIF                                                               
C  Apply a consistent minimum moisture criterion in the mixing line       
C  calculations.                                                          
      QQG=MAX(QG(J,NCRB),QGMIN)                                           
      TSB=FTSL(TG(J,NCRB),QQG,PL(NCRB))                                   
      PSLB=PL(NCRB)*(TSB/TG(J,NCRB))**RKAP                                
      THM=0.5*(THB+TG(J,NCRT)/PFK(NCRT))                                  
      QM =0.5*(QG(J,NCRB)+QG(J,NCRT))                                     
      IF (QM.LT.QGMIN) THEN                                               
         WRITE(2,2002) J,JH,NCRB,NCRT,CQ*QM                               
      ENDIF                                                               
      QM=MAX(QM,QGMIN)                                                    
      TMT=THM*PFK(NCRT)                                                   
      TSLM=FTSL(TMT,QM,PL(NCRT))                                          
      PSLM=PL(NCRT)*(TSLM/TMT)**RKAP                                      
      TLAPSM=(THM-THB)/(PSLM-PSLB)                                        
      TLAPSM=MAX(TLAPSM,TLIM) ! limit stability                           
      TLAPSM=MIN(TLAPSM,0.0) ! prevents instability                       
      IF (TLAPSM.GT.0.0) THEN                                             
         WRITE(2,2103) J,JH,NCRB,NCRT                                     
      ELSEIF (TLAPSM.LT.TLIM)THEN                                         
         WRITE(2,2104) J,JH,NCRB,NCRT,TLIM                                
      ENDIF                                                               
                                                                          
      DO 60 L=NCRT,NCRB                                                   
        THRL=THB+(PL(L)-PL(NCRB))*TLAPSM                                  
        TR(L)=THRL*PFK(L)                                                 
   60 CONTINUE                                                            
C                                                                         
      TLPH=TG(J,NCRT)+TG(J,NCRTP)                                         
      FPH=FB(NCRT)/(TLPH*TLPH)                                            
      DQPH=QG(J,NCRTP)-QG(J,NCRT)                                         
      DTPH=SK(NCRT)*(TG(J,NCRTP)-TR(NCRTP))-(TG(J,NCRT)-TR(NCRT))         
      TTMC(NCRT)=AKTCU*FPH*DTPH/DSIGMA(NCRT)                              
      QTMC(NCRT)=AKQCU*FPH*DQPH/DSIGMA(NCRT)                              
C                                                                         
      IF(NLEV.GT.2) THEN                                                  
        DO 10 L=NCRTP,NCRBM                                               
          LP=L+1                                                          
          FMH=FPH                                                         
          DQMH=DQPH                                                       
          DTMH=DTPH/SK(L-1)                                               
          TLPH=TG(J,LP)+TG(J,L)                                           
          FPH=FB(L)/(TLPH*TLPH)                                           
          DQPH=QG(J,LP)-QG(J,L)                                           
          DTPH=SK(L)*(TG(J,LP)-TR(LP))-(TG(J,L)-TR(L))                    
          TTMC(L)=AKTCU*(FPH*DTPH-FMH*DTMH)/DSIGMA(L)                     
          QTMC(L)=AKQCU*(FPH*DQPH-FMH*DQMH)/DSIGMA(L)                     
          CTCR(L,IHEM)=CTCR(L,IHEM)+1.0                                   
   10   CONTINUE                                                          
      ENDIF                                                               
C                                                                         
      TTMC(NCRB)=-AKTCU*FPH*DTPH/(SK(NCRBM)*DSIGMA(NCRB))                 
      QTMC(NCRB)=-AKQCU*FPH*DQPH/DSIGMA(NCRB)                             
C                                                                         
      SINT=SSUM(NLEV,DSIGMA(NCRT),1)                                      
      TTCOR=SDOT(NLEV,TTMC(NCRT),1,DSIGMA(NCRT),1)/SINT                   
      DO 20 L=NCRT,NCRB                                                   
   20 TTMC(L)=TTMC(L)-TTCOR                                               
C                                                                         
      IF(CURHM.LT.1.00001) THEN                                           
        DO 30 L=NCRT,NCRB                                                 
          TGP=TG(J,L)+DELT2C*TTMC(L)                                      
          QGP=QG(J,L)+DELT2C*QTMC(L)                                      
          QS=CURHM*ESCON*PQSAT(TGP)/SIGMA(L)                              
          QEXS=QGP-QS                                                     
          QEX=QEXS/(1.0+CTQ*QS*ESCONB/(TGP*TGP))                          
          IF (QEXS.LT.0.0) THEN                                           
            QEX=0.0                                                       
          END IF                                                          
          QTCOR=-QEX/DELT2C                                               
          TTMC(L)=TTMC(L)-CTQ*QTCOR                                       
          QTMC(L)=QTMC(L)+QTCOR                                           
          RRCR(J)=RRCR(J)-QTCOR*CLR(L)*PLG(J)                             
   30   CONTINUE                                                          
      ENDIF                                                               
C                                                                         
      DO 70 L=NCRT,NCRB                                                   
        TG(J,L)=TG(J,L)+DELT2C*TTMC(L)                                    
        QG(J,L)=QG(J,L)+DELT2C*QTMC(L)                                    
   70 CONTINUE                                                            
C                                                                         
      RETURN                                                              
      END                                                                 
