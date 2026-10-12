C**********************************************************               
C             SUBROUTINE DGRMLT                                           
C**********************************************************               
      SUBROUTINE DGRMLT(IH,IPASS)
      use tracerclds, only : ADYE, PDYEFIX, TRELAXORB,
     &                       KDYEMODE, TCONDDYE, TEVAPDYE, QSATDYE,
     &                       SETTLECOL, FCONDYE
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
      INTEGER IPASS
C
C     Needed for PORB, OBLIQ (nightside settling source term below).
       COMMON/VARPARAM/OOM_IN, LPLOTMAP,NLPLOTMAP_IN,RFCOEFF_IN,
     & NTSTEP_IN, NSKIP_IN, BOTRELAXTIME, FBASEFLUX, FORCE1DDAYS,
     & OPACIR_POWERLAW, OPACIR_REFPRES, SOLC_IN, TOAALB,
     & PORB, OBLIQ, ECCEN
       LOGICAL LPLOTMAP
       DIMENSION QSED(NL)
C     One column for SETTLECOL, and the two-phase post-exchange state.
       DIMENSION PCSED(NL),TCSED(NL),QCSED(NL),QVSED(NL)
C     Per-level, grain-independent settling inputs, shared by all dyes.
       DIMENSION GLAMC(NL),GRHOC(NL),GETAC(NL),FSCLC(NL)
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

      if (LOLDBL) then
         if (LBL) call BLAYER                                             
      else                                                                
         if (LBL) call BLSURF                                             
      endif


      IF(LVD) CALL VDIFF
      IF(LCR) CALL CONVEC
      IF(LLR) CALL LSCRN

C     IPASS=1 is the pre-pass for RADIATION_ALLLATS: it runs the T0
C     addition (l.242), the PLG EXP (l.271), VDIFF and CONVEC, then
C     returns so the radiation pre-pass sees exactly the TG/PLG that
C     RADIATION saw when it lived inside this routine.  Returning here
C     also stays clear of the accumulator store-back further below.
      IF (IPASS.EQ.1) RETURN
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
      ENDIF
C                                                                         
      DO 200 L=1,NL                                                       
       DO 200 J=1,IGC                                                     
        UNLG(J,L)=UNLG(J,L)+UTVD(J,L)                                     
        VNLG(J,L)=VNLG(J,L)+VTVD(J,L)
        TNLG(J,L)=TNLG(J,L)+TTVD(J,L)+TTCR(J,L)+TTLR(J,L)+TTRD(J,L)

        QNLG(J,L)=QNLG(J,L)+QTVD(J,L)+QTCR(J,L)+QTLR(J,L)
  200 CONTINUE
C
C     Deep reservoir: relax toward 1 below PDYEFIX at all local times.
C     Must be a tendency, not a reset of TRAG, which is rebuilt from
C     spectral TRA every call. RLXDYE follows the DAMP=1/(PI2*RESTIM)
C     convention of cinires.f. TRELAXORB is in orbits and PORB the
C     period in rotation days, hence TRELAXORB*PORB; PORB=0 means no
C     orbital motion, so fall back to TRELAXORB rotation days.
      IF (PORB.NE.0.0) THEN
         TRELAXDYE=TRELAXORB*PORB
      ELSE
         TRELAXDYE=TRELAXORB
      ENDIF
      RLXDYE=1.0/(PI2*TRELAXDYE)
C
C     Two-phase exchange: the fraction of the gap to equilibrium closed
C     in DELT2, from exact integration of dq/dt=(qeq-q)/tau. TCONDDYE
C     and TEVAPDYE are seconds, so tau*WW is model time; tau<=0 means
C     instantaneous equilibrium.
      IF (KDYEMODE.EQ.2) THEN
         FCDYE=1.0
         FEDYE=1.0
         IF (TCONDDYE.GT.0.0) FCDYE=1.0-EXP(-DELT2/(TCONDDYE*WW))
         IF (TEVAPDYE.GT.0.0) FEDYE=1.0-EXP(-DELT2/(TEVAPDYE*WW))
      ENDIF
C
C     Gravitational settling of the dyes (tracers 2..NTRAC):
C     dq/dt = (T/P) d/dz[ (P/T) q VFALL ]. For a mass mixing ratio the
C     density factor is P/T INSIDE the derivative and T/P outside;
C     reversed, a uniform q gives dq/dt>0, not the physical decay. No
C     height coordinate here, so dz=-(GASCON*T/GA) dlnP per sigma
C     layer. ADYE(KK) and RHODYE come from the tracerclds module;
C     CALC_V_FALL is
C     SI, as are GA, SIGMA*PLG*P0 and TG*CT, so VFALL is m/s for /WW.
C     SEDSLN/SEDSLT use the PORB/OBLIQ formula of NIKOSRAD; night is
C     where the local zenith-angle cosine is negative.
      IF (PORB.NE.0.0) THEN
         SEDSLN=(1./PORB-1.)*KOUNT*360./ITSPD
         SEDSLN=MOD(SEDSLN,360.)
      ELSE
         SEDSLN=0.0
      ENDIF
      IF (OBLIQ.EQ.0.0) THEN
         SEDSLT=0.0
      ELSE
         SEDSLT=ASIN(SIN(OBLIQ*PI/180.)*SIN(PI2*KOUNT/ITSPD/PORB))
     +          *180./PI
      ENDIF
C
      IOFM=0
      DO 320 IHEM=1,NHEM
         IF (IHEM.EQ.1) THEN
            ALATLC=ALAT(JH)
         ELSE
            ALATLC=-ALAT(JH)
         ENDIF
         DO 310 I=1,MG
            J=I+IOFM
            ALON=360.0*REAL(I-1)/REAL(MG)
            COSZLC=SIN(ALATLC/360.*PI2)*SIN(SEDSLT/360.*PI2)
     +            +COS(ALATLC/360.*PI2)*COS(SEDSLT/360.*PI2)
     +            *COS((ALON-SEDSLN)/360.*PI2)
C           Calculate gas properties for each level
            DO 290 L=1,NL
               PCSED(L)=SIGMA(L)*PLG(J)*P0
               TCSED(L)=TG(J,L)*CT
               CALL SED_GAS(PCSED(L),TCSED(L),GLAMC(L),GRHOC(L),
     +                      GETAC(L))
               IF (KDYEMODE.EQ.1) THEN
                  FSCLC(L)=FCONDYE(PCSED(L),TCSED(L))
               ELSE
                  FSCLC(L)=1.0
               ENDIF
  290       CONTINUE
C           With LFLUX the dyes are advected in flux form, so TRAG holds
C           p*q (PLG=p*/P0 by now, see the EXP above), not q. Everything
C           below works in q: divide by PSDY on the way in, multiply the
C           target by PSDY on the way out, as QG/QNLG are for tracer 1.
            IF (LFLUX) THEN
               PSDY=PLG(J)
            ELSE
               PSDY=1.0
            ENDIF
            IF (KDYEMODE.EQ.1) THEN
C           Tracers 2..NTRAC are the dyes: identical but for ADYE(KK),
C           so one run compares grain sizes. Tracer 1 is water vapour.
C           Each dye is total species material, and only its condensed
C           fraction (FCONDYE) falls.
            DO 308 KK=2,NTRAC
C              TRAG strides levels by IGC, so copy out the column.
               DO 301 L=1,NL
                  QCSED(L)=TRAG(J,L,KK)/PSDY
  301          CONTINUE
               CALL SETTLECOL(QCSED,PCSED,TCSED,SIGMA,DSIGMA,ADYE(KK),
     +                        GA,GASCON,WW,DELT2,FSCLC,GLAMC,GRHOC,
     +                        GETAC,QSED)
               DO 302 L=1,NL
                  TRANLG(J,L,KK)=TRANLG(J,L,KK)
     +                 +(PSDY*QSED(L)-TRAG(J,L,KK))/DELT2
  302          CONTINUE
C              Deep reservoir, relaxing q toward 1.
               DO 305 L=1,NL
                  IF (PCSED(L).GT.PDYEFIX) THEN
                     TRANLG(J,L,KK)=TRANLG(J,L,KK)
     +                    +(PSDY-TRAG(J,L,KK))*RLXDYE
                  ENDIF
  305          CONTINUE
  308       CONTINUE
            ELSE IF (KDYEMODE.EQ.2) THEN
C           Vapour/solid pairs: vapour KV=2P, solid KC=2P+1. First the
C           pair relaxes toward equilibrium, solid = max(total-QSATDYE,0),
C           exactly integrated over DELT2 so any timescale is stable and
C           the total is untouched; then the solid alone settles.
            DO 318 KV=2,NTRAC-1,2
               KC=KV+1
               DO 312 L=1,NL
                  QVDY=TRAG(J,L,KV)/PSDY
                  QCDY=TRAG(J,L,KC)/PSDY
                  QTDY=QVDY+QCDY
                  QEDY=MAX(QTDY-QSATDYE(PCSED(L),TCSED(L)),0.0)
                  IF (QEDY.GT.QCDY) THEN
                     QCDY=QCDY+(QEDY-QCDY)*FCDYE
                  ELSE
                     QCDY=QCDY+(QEDY-QCDY)*FEDYE
                  ENDIF
                  QCSED(L)=QCDY
                  QVSED(L)=QTDY-QCDY
  312          CONTINUE
               CALL SETTLECOL(QCSED,PCSED,TCSED,SIGMA,DSIGMA,ADYE(KC),
     +                        GA,GASCON,WW,DELT2,FSCLC,GLAMC,GRHOC,
     +                        GETAC,QSED)
               DO 314 L=1,NL
                  TRANLG(J,L,KV)=TRANLG(J,L,KV)
     +                 +(PSDY*QVSED(L)-TRAG(J,L,KV))/DELT2
                  TRANLG(J,L,KC)=TRANLG(J,L,KC)
     +                 +(PSDY*QSED(L)-TRAG(J,L,KC))/DELT2
  314          CONTINUE
C              Deep reservoir: relax the pair's total q toward 1,
C              supplied as vapour; the exchange above decides its phase.
               DO 316 L=1,NL
                  IF (PCSED(L).GT.PDYEFIX) THEN
                     TRANLG(J,L,KV)=TRANLG(J,L,KV)
     +                    +(PSDY-TRAG(J,L,KV)-TRAG(J,L,KC))*RLXDYE
                  ENDIF
  316          CONTINUE
  318       CONTINUE
            ELSE
               WRITE(*,*)'ERROR: KDYEMODE must be 1 or 2 (cdgrmlt.f)'
               STOP
            ENDIF
  310    CONTINUE
         IOFM=MGPP
  320 CONTINUE
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
