C**********************************************************               
C             SUBROUTINE DGRMLT                                           
C**********************************************************               
      SUBROUTINE DGRMLT(IH)                                                   
C                                                                         
C     COMPUTE DIABATIC TENDENCIES IN GRID POINT SPACE FOR PRESENT LAT.    
C     ACCUMULATE TIME AVERAGES FOR PRINTED OUTPUT AND HISTORY             
C     ASSUMES OUTPUT AND HISTORY TIMESTEPS ARE SYNCHRONISED               
C                                                                         
C     Full Physics                                                        
C     Piers 6/2/97 Version 2.0                                            
C     Water vapour (things beginning with Q have changed                  
C     into tracer no 1)                                                   
C     and converted from Volume mixing  ratio to Mass mixing ratio        
C     As water vapour advection is now done in Flux form                  
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
C     Switches counters and constants controlling type and frequency of   
C     model output                                                        
C                                                                         
      COMMON/OUTCON/RNTAPE,NCOEFF,NLAT,INLAT,INSPC                        
     +              ,RNTAPO                                               
     +              ,KOUNTP,KOUNTE,KOUNTH,KOUNTR                          
     +              ,KOUTP,KOUTE,KOUTH,KOUTR,DAY                          
     +              ,SQR2,RSQR2,EAM1,EAM2,TOUT1,TOUT2,RMG                 
     +              ,LSPO(NL),LGPO(NL)                                    
     $              ,LSHIST,LMINIH                                        
      LOGICAL LSHIST,LMINIH                                               
      LOGICAL LSPO,LGPO                                                   
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
C                                                                         
C     Constant arrays and variables associated with time and vertical     
C     differencing. Also counters.                                        
C                                                                         
      COMMON/BATS/  BEGDAY,CTRA(NTRAC),BM1(IDE),AK(NNP),AQ(NL2),G(NL2)              
     +              ,TAU(NL2),KOUNT,KITS,KSTART,KTOTAL,KRUN,ITSPD           
     +              ,DELT,DELT2,CV,CG,CT,CQ,PNU,PNU2,PNU21                
     +              ,NTRACO,KOLOUR(NTRAC),RGG(NL2)            
     +              ,BEGDOY,DOY                                           
C                                                                         
            COMMON/PHYS/  CCR,RCON,DTBUOY,TSLA,TSLB,TSLC,TSLD,CUT1,CUT2
     :              ,TSTAR(IGC,JG),QSTAR(IGC,JG),FRAD(JG,NHEM)            
     :              ,TSTARO(IGC,JG),TDEEPO(IGC,JG),smstar(igc,jg)         
     :              ,tdeep(igc,jg),hsnow(igc,jg),sqstar(igc,jg)           
     :              ,SALB(IGC,JG),SBAL(IGC,JG),BLCD(IGC)                  
     :              ,SVEGE(IGC,JG),CD,DRAG,BLVAD,BLA,BLRH,BLVB(IGC)       
     :              ,AKVV,AKTV,AKQV,ESCONA,ESCONB,EPSIQ,CTQ,CCC           
     : ,ctqi,sdsn,shcs,shcsp,shcsn,skse,sksn,slhf,sd1,sd2,sdw             
     :        ,ssmc,sdsnd,sasnow,saice,shsstar,shsmax
     :     ,LOC,LNOICE,LOLDBL,LCOND,LNNSK
     :              ,NLCR,CURHM,AKTC,AKQC,CUBMT,CBADJT,CBADJP
     :              ,SKAP(NL),SK(NLM),FWS(NL),CLR(NL),FB(NLM)             
     :              ,TTDC(NL),QTDC(NL),TTMC(NL),QTMC(NL),TC(NL),QC(NL)    
     :              ,CTCR(NL,NHEM),CTLR(NL,NHEM)
     :              ,LBL,LVD,LCR,LLR,LRD,LCUBM,LCBADJ
     :              ,LSL,NAVRD,NAVWT,DELT2C,SHCO,SHCI,ITSLL,ITSLO,NCUTOP
                LOGICAL LBL,LVD,LCR,LLR,LRD,LCUBM,LCBADJ,LSL,LOC 
     :                 ,LNOICE,LOLDBL,LCOND,LNNSK                       
C                                                                         
      COMMON/PTENDZ/UTVDZ(IGG),VTVDZ(IGG),TTVDZ(IGG),QTVDZ(IGG)           
     :              ,TTCRZ(IGG),QTCRZ(IGG),TTLRZ(IGG),QTLRZ(IGG)          
     :              ,TTRDZ(IGG),CTCRZ(IGG),CTLRZ(IGG)                     
     :              ,UTBLZ(JGG),VTBLZ(JGG),TTBLZ(JGG),QTBLZ(JGG)          
     :              ,AUTVDZ(IGG),AVTVDZ(IGG),ATTVDZ(IGG),AQTVDZ(IGG)      
     :              ,ATTCRZ(IGG),AQTCRZ(IGG),ATTLRZ(IGG),AQTLRZ(IGG)      
     :              ,ATTRDZ(IGG),ACTCRZ(IGG),ACTLRZ(IGG)                  
     :              ,AUTBLZ(JGG),AVTBLZ(JGG),ATTBLZ(JGG),AQTBLZ(JGG)      
C                                                                         
      LOGICAL LSUM,LWRITE                                                 
       COMMON/CPIERS/ICFLAG(IGC,5,2),CFRAC(IGC,5),PNET(IGC,JG)            
     :               ,SNET(IGC,JG),RRFLUX(IGC,JG,6)                       
     :               ,TTSW(IGC,NL),TTLW(IGC,NL)                           
      REAL UTBL(IGC),VTBL(IGC),TTBL(IGC),QTBL(IGC)                        
       COMMON/GSG/GSG(IGC,JG)                                             
      INTEGER IFIRST                                                      
      REAL TROPHT(MG,NHEM,JG)                                             
c                                                                         
      common/gridsss/assbl1(igc,jg),ashbl1(igc,jg),aslbl1(igc,jg),        
     *               arrcr1(igc,jg),arrlr1(igc,jg),                       
     *               arflux1(igc,6,jg),asfld1(igc,6,jg),                  
     *               acld1(igc,4,jg)                                      



C     Setup moisture variables by equivilencing them to                   
C     Tracer No. 1                                                        
C                                                                         
      REAL QG(IGC,NL),QNLG(IGC,NL),QTLR(IGC,NL),QTCR(IGC,NL)              
      EQUIVALENCE (QG(1,1),TRAG(1,1,1)) , (QNLG(1,1),TRANLG(1,1,1))       
     *            ,(QTLR(1,1),UTRAG(1,1,1)),(QTCR(1,1),VTRAG(1,1,1))      
C                                                                         
      EQUIVALENCE (UTBL(1),UNLG(1,NL)),(VTBL(1),VNLG(1,NL))               
     :           ,(TTBL(1),TNLG(1,NL)),(QTBL(1),TRANLG(1,NL,1))           
C                                                                         
      SAVE IFIRST                                                         
      DATA IFIRST/1/                                                      
      LSUM=.FALSE.                                                        
      LWRITE=.FALSE.                                                      
      IF(KOUTH.GE.1.AND.KOUTH.LE.KOUNTH) LSUM=.TRUE.                      
      IF(KOUTH.EQ.KOUNTH) LWRITE=.TRUE.                                   
      RKP=1.0/KOUNTH                                                      
      DELT2C=DELT2                                                        
C                                                                         
      IF(LSUM) THEN                                                       
         IF(KOUTH.EQ.1) THEN                                               
            do k=1,6                                                         
               do j=1,igc                                                      
                  arflux(j,k)=0.                                                
                  asfld(j,k)=0.                                                 
               enddo                                                           
            enddo                                                            
            do k=1,4                                                         
               do j=1,igc                                                      
                  acld(j,k)=0.                                                  
               enddo                                                           
            enddo                                                            
            DO 10 J=1,IGC                                                    
               ASSBL(J)=0.0                                                   
               ASHBL(J)=0.0                                                   
               ASLBL(J)=0.0                                                   
               ARRCR(J)=0.0                                                   
               ARRLR(J)=0.0                                                   
 10         CONTINUE                                                         
         ELSE                                                              
            do i=1,igc                                                       
               assbl(i)=assbl1(i,jh)                                          
               ashbl(i)=ashbl1(i,jh)                                          
               aslbl(i)=aslbl1(i,jh)                                          
               arrcr(i)=arrcr1(i,jh)                                          
               arrlr(i)=arrlr1(i,jh)                                          
            end do                                                           
            do k=1,6                                                         
               do i=1,igc                                                     
                  arflux(i,k)=arflux1(i,k,jh)                                 
                  asfld(i,k)=asfld1(i,k,jh)                                   
               end do                                                         
            end do                                                           
            do k=1,4                                                         
               do i=1,igc                                                     
                  acld(i,k)=acld1(i,k,jh)                                     
               end do                                                         
            end do                                                           
c                                                                         
c        READ(NAVRD)ASSBL,ASHBL,ASLBL,ARRCR,ARRLR                         
c     $        ,arflux,asfld,acld                                         
c     
         ENDIF                                                             
      ENDIF                                                               
C                                                                         
      DO 20 L=1,NL                                                        
         DO 20 J=1,IGC                                                      
            UNLG(J,L)=0.0                                                     
            VNLG(J,L)=0.0                                                     
            TNLG(J,L)=0.0                                                     
            QNLG(J,L)=0.0                                                     
            QTVD(J,L)=0.0                                                     
            UTVD(J,L)=0.0                                                     
            VTVD(J,L)=0.0                                                     
            TTVD(J,L)=0.0                                                     
            TTCR(J,L)=0.0                                                     
            QTCR(J,L)=0.0                                                     
            TTLR(J,L)=0.0                                                     
            QTLR(J,L)=0.0                                                     
            TTRD(J,L)=0.0                                                     
            TG(J,L)=TG(J,L)+T0(L)                                             
 20   CONTINUE                                                            
      DO 30 J=1,IGC                                                       
C sets ICFLAG and CFRAC to ZERO at each timestep                          
         DO L=1,5                                                          
C bottom level is complex for ICFLAG                                      
C need to make sure that max and min commands in MORC.632                 
C change things if cloud is there                                         
C PMF bug fix 11-5-98                                                     
            ICFLAG(J,L,1)=2                                                
            ICFLAG(J,L,2)=NLP                                              
            CFRAC(J,L)=0.0                                                 
         ENDDO                                                             
         TXBL(J)=0.0                                                       
         TYBL(J)=0.0                                                       
         SSBL(J)=0.0                                                       
         SHBL(J)=0.0                                                       
         SLBL(J)=0.0                                                       
         RRCR(J)=0.0                                                       
         RRLR(J)=0.0                                                       
 30   CONTINUE                                                            
      DO 60 IHEM=1,NHEM                                                   
         DO 40 L=1,NL                                                      
            CTCR(L,IHEM)=0.0                                                
 40      CTLR(L,IHEM)=0.0                                                  
         IOFM=(IHEM-1)*MGPP                                                
         DO 50 I=1,MG                                                      
            J=I+IOFM                                                        
            PLG(J)=EXP(PLG(J))                                              
 50      CONTINUE                                                          
         IF (LFLUX) THEN                                                   
C  Convert from volume mixing ratio to mass mixing ratio.                 
            DO 52 L=1,NL                                                    
               DO 51 I=1,MG                                                  
                  J=I+IOFM                                                   
                  QG(J,L)=QG(J,L)/PLG(J)                                     
 51            CONTINUE                                                      
 52         CONTINUE                                                        
         ENDIF                                                             
 60   CONTINUE                                                            
      IF (LCSFCT) THEN                                                    
         IF (.NOT.LPERPET.OR.IFIRST.EQ.1)                                   
     +        CALL SFCT(PLG,JH,IFIRST,TROPHT)                              
         IF (JH.EQ.JG) IFIRST=0                                             
      ENDIF                                                               
C                                                                         
      if (LOLDBL) then                                                    
         if (LBL) call BLAYER                                             
      else                                                                
         if (LBL) call BLSURF                                             
      endif                                                               
      IF(LVD) CALL VDIFF
      IF(LCR) CALL CONVEC            
      IF(LLR) CALL LSCRN

      IF(LRD) CALL RADIATION(TROPHT,IH)
      if (LBL.AND.(.NOT.LOLDBL)) CALL SURFM                               
C                                                                         
      do j=1,igc                                                          
         sfld(j,1)=salb(j,jh)                                             
         sfld(j,2)=tstar(j,jh)                                            
         sfld(j,3)=tdeep(j,jh)                                            
         sfld(j,4)=qstar(j,jh)                                            
         sfld(j,5)=smstar(j,jh)                                           
         sfld(j,6)=hsnow(j,jh)                                            
         cld(j,1)=cfrac(j,1)                                              
         cld(j,2)=cfrac(j,2)                                              
         cld(j,3)=cfrac(j,3)                                              
         cld(j,4)=cfrac(j,4)+cfrac(j,5)                                   
         do k=1,6                                                         
            rflux(j,k)=rrflux(j,jh,k)                                     
         enddo                                                            
C rflux(j,2) is +downwards surface flux Wm-2                              
        rflux(j,2)=rflux(j,1)-rflux(j,2)+rflux(j,3)-rflux(j,4)-           
     :    CV*1.0E5*(shbl(j)+slbl(j))                                      
      enddo                                                               
      IF(LSUM) THEN                                                       
      do k=1,6                                                            
         do j=1,igc                                                       
            arflux(j,k)=arflux(j,k)+rflux(j,k)*RKP                        
            asfld(j,k)=asfld(j,k)+sfld(j,k)*RKP                           
         enddo                                                            
      enddo                                                               
      do k=1,4                                                            
         do j=1,igc                                                       
            acld(j,k)=acld(j,k)+cld(j,k)*RKP                              
         enddo                                                            
      enddo                                                               
      DO 100 J=1,IGC                                                      
       SSBL(J)=SQRT(TXBL(J)*TXBL(J)+TYBL(J)*TYBL(J))                      
       ASSBL(J)=ASSBL(J)+SSBL(J)*RKP                                      
       ASHBL(J)=ASHBL(J)+SHBL(J)*RKP                                      
       ASLBL(J)=ASLBL(J)+SLBL(J)*RKP                                      
       ARRCR(J)=ARRCR(J)+RRCR(J)*RKP                                      
       ARRLR(J)=ARRLR(J)+RRLR(J)*RKP                                      
  100 CONTINUE                                                            
      DO 140 IHEM=1,NHEM                                                  
       IOF=(IHEM-1)*MGPP+1                                                
       K1=JH*(2-IHEM)+(JGG+1-JH)*(IHEM-1)                                 
       K=K1                                                               
       DO 120 L=1,NL                                                      
        UTVDZ(K)=SSUM(MG,UTVD(IOF,L),1)*RMG                               
        VTVDZ(K)=SSUM(MG,VTVD(IOF,L),1)*RMG                               
        TTVDZ(K)=SSUM(MG,TTVD(IOF,L),1)*RMG                               
        QTVDZ(K)=SSUM(MG,QTVD(IOF,L),1)*RMG                               
        TTCRZ(K)=SSUM(MG,TTCR(IOF,L),1)*RMG                               
        QTCRZ(K)=SSUM(MG,QTCR(IOF,L),1)*RMG                               
        TTLRZ(K)=SSUM(MG,TTLR(IOF,L),1)*RMG                               
        QTLRZ(K)=SSUM(MG,QTLR(IOF,L),1)*RMG                               
        TTRDZ(K)=SSUM(MG,TTRD(IOF,L),1)*RMG                               
        CTCRZ(K)=CTCR(L,IHEM)                                             
        CTLRZ(K)=CTLR(L,IHEM)                                             
  120  K=K+JGG                                                            
       UTBLZ(K1)=SSUM(MG,UTBL(IOF),1)*RMG                                 
       VTBLZ(K1)=SSUM(MG,VTBL(IOF),1)*RMG                                 
       TTBLZ(K1)=SSUM(MG,TTBL(IOF),1)*RMG                                 
       QTBLZ(K1)=SSUM(MG,QTBL(IOF),1)*RMG                                 
  140 CONTINUE                                                            
c      IF(LWRITE) THEN                                                    
c        WRITE(NAVWT)ASSBL,ASHBL,ASLBL,ARRCR,ARRLR                        
c     $        ,arflux,asfld,acld                                         
c     :  ,SSBL,SHBL,SLBL,RRCR,RRLR                                        
c     $        ,rflux,sfld,cld                                            
c      ELSE                                                               
c        WRITE(NAVWT)ASSBL,ASHBL,ASLBL,ARRCR,ARRLR                        
c     $        ,arflux,asfld,acld                                         
c      ENDIF                                                              
c                                                                         
      do i=1,igc                                                          
         assbl1(i,jh)=assbl(i)                                            
         ashbl1(i,jh)=ashbl(i)                                            
         aslbl1(i,jh)=aslbl(i)                                            
         arrcr1(i,jh)=arrcr(i)                                            
         arrlr1(i,jh)=arrlr(i)                                            
      end do                                                              
      do k=1,6                                                            
         do i=1,igc                                                       
            arflux1(i,k,jh)=arflux(i,k)                                   
            asfld1(i,k,jh)=asfld(i,k)                                     
         end do                                                           
      end do                                                              
      do k=1,4                                                            
         do i=1,igc                                                       
            acld1(i,k,jh)=acld(i,k)                                       
         end do                                                           
      end do                                                              
      ENDIF ! IF(LSUM)                                                    
C                                                                         
      DO 200 L=1,NL                                                       
       DO 200 J=1,IGC                                                     
        UNLG(J,L)=UNLG(J,L)+UTVD(J,L)                                     
        VNLG(J,L)=VNLG(J,L)+VTVD(J,L)                                     
        TNLG(J,L)=TNLG(J,L)+TTVD(J,L)+TTCR(J,L)+TTLR(J,L)+TTRD(J,L)       
        QNLG(J,L)=QNLG(J,L)+QTVD(J,L)+QTCR(J,L)+QTLR(J,L)
  200 CONTINUE  
C                                                                        
      IF (LFLUX) THEN                                                     
C        Convert to volume mixing ratio from mass mixing ratio.           
         DO 211 L=1,NL                                                    
            DO 210 J=1,IGC                                                
               QNLG(J,L)=PLG(J)*QNLG(J,L)                                 
  210       CONTINUE                                                      
  211    CONTINUE                                                         
      ENDIF                                                               
C                                                                         
      RETURN                                                              
      END                                                                 
