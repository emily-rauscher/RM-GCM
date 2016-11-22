C**********************************************************               
C             SUBROUTINE CONVEC                                           
C**********************************************************               
      SUBROUTINE CONVEC                                                   
C                                                                         
C     CONVECTION SCHEME. THIS ROUTINE FINDS UNSTABLE LAYERS AND           
C     CALLS RELEVANT ROUTINES TO CALCULATE TENDENCIES ETC.                
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
     :       ,LNOICE,LOLDBL,LCOND,LNNSK                                   
C                                                                         
C                                                                         
      LOGICAL LCDRY,LCWET                                                 
      REAL ESCON(NL),SIGMIN,TMIN                                          
C                                                                         
C     Setup moisture variables by equivilencing them to                   
C     Tracer No. 1                                                        
C                                                                         
      REAL QG(IGC,NL),QNLG(IGC,NL),QTLR(IGC,NL),QTCR(IGC,NL)              
      EQUIVALENCE (QG(1,1),TRAG(1,1,1)) , (QNLG(1,1),TRANLG(1,1,1))       
     *            ,(QTLR(1,1),UTRAG(1,1,1)),(QTCR(1,1),VTRAG(1,1,1))      
C                                                                         
C                                                                         
c 28-05-97 PMF Water vapour fix to stop too low values                    
c takes water vapour from level below                                     
c this produces a slight vertical advection of Q                          
c assumes no -ve water vapour in bottom layer                             
c                                                                         
      REAL QGO(NL)                                                        
c                                                                         
      PRESSMIN=0.20000E+29  !  Minimum parcel pressure for                    
                        !  which moist convection is attempted.           
      TMIN=0.31000418E+29   !  Minimum parcel temperature (dedim.) for        
                        !  which moist convection is attempted.           
      LTOP=NL-NLCR      !  NLCR gives highest level for which             
                        !  convection is attempted.                       



      DO 800 IHEM=1,NHEM                                                  
        IOFM=(IHEM-1)*MGPP                                                
        DO 700 I=1,MG                                                     
          J=I+IOFM                                                        
C         IPH=J                                                           
          DO 10 L=1,NL                                                    
            TTDC(L)=0.0                                                   
            QTDC(L)=0.0                                                   
            TTMC(L)=0.0                                                   
            QTMC(L)=0.0                                                   
            ESCON(L)=1./(PLG(J)*SIGMA(L))                                 
            QGO(L)=QG(J,L)                                                
   10     CONTINUE                                                        
c Water vapour fix                                                        
          IFL=0                                                           
          DO L=1,NL-1                                                     
            IF (QG(J,L).LT.1.0E-6) THEN                                   
              IFL=1                                                       
              QEXS=1.0E-6-QG(J,L)                                         
              IF (L.EQ.1) THEN                                            
                QEXS=QEXS*SIGMAH(1)                                       
              ELSE                                                        
                QEXS=QEXS*(SIGMAH(L)-SIGMAH(L-1))                         
              ENDIF                                                       
              IF (L.NE.NL-1) THEN                                         
                QEX=QEXS/(SIGMAH(L+1)-SIGMAH(L))                          
              ELSE                                                        
                QEX=QEXS/(1.0-SIGMAH(L))                                  
              ENDIF                                                       
c assumes no large scale rain here in layer                               
              QTLR(J,L)=(1.0E-6-QGO(L))/DELT2C                            
              QTLR(J,L+1)=QTLR(J,L+1)-QEX/DELT2C                          
              QG(J,L)=1.0E-6                                              
              QG(J,L+1)=QG(J,L+1)-QEX                                     
            ENDIF                                                         
          ENDDO                                                           
C                                                                         
C     DRY ADJUSTMENT                                                      
C                                                                         
          L=NL                                                            
C     BEGIN NEW PARCEL CURVE                                              
   20     NCRB=L                                                          
          LCDRY=.FALSE.                                                   
          TCLP=TG(J,L)                                                    
          QCL=QG(J,L)                                                     
C     DRY ADIABAT TO NEXT LEVEL                                           
   30     L=L-1                                                           
          IF(L.LT.LTOP) GOTO 40                                           
          TCL=TCLP*SK(L)                                                  
          IF(TCL.LE.TG(J,L)) GOTO 40  !If stable layer                    
C ER Modif: commenting out saturation check
C     BUOYANT - IF SATURATED - MOIST CONVECTION - IGNORE                  
c          QSL=ESCON(L)*PQSAT(TCL)                                         
c          IF ( QCL.GE.QSL ) THEN   ! Is parcel supersaturated?            
c!KM Bypass supersaturation issue and Tmin criterion avoided
c            IF ( ( (SIGMA(NCRB)*PLG(J)).GE.PRESSMIN )                     
c     :        .AND.( TG(J,NCRB).GE.TMIN ) ) GOTO 40                       
C If parcel is in stratosphere, do dry rather than moist convection       
C as the dry and moist adiabats are approximately parallel and            
C dry convection is a lot cheaper.                                        
c          ENDIF                                                           
C     DRY CONVECTION - CONTINUE PARCEL CURVE UP                           
   35     LCDRY=.TRUE.                                                    
          TCLP=TCL                                                        
          GOTO 30                                                         
   40     CONTINUE                                                        
C     STABLE LAYER OR MODEL TOP - ADJUST ANY UNSTABLE LAYER BELOW         
          IF(LCDRY) THEN                                                  
            NCRT=L+1            
!!       write(*,*) 'Calling Dry Convection between levels:',NCRB,NCRT     
            CALL DRYADJ(NCRB,NCRT,J,IHEM)                                 
          ENDIF                                                           
          IF(L.GT.LTOP) GOTO 20                                           

! KM Bypass Moist convection
!          GOTO 900
C                                                                         
C     MOIST CONVECTION                                                    
C         
          L=NL                                                            
C         BEGIN NEW PARCEL CURVE                                          
  100     NCRB=L                                                          
          LCWET=.FALSE. 
          IF ( (SIGMA(NCRB)*PLG(J)).LT.PRESSMIN ) THEN                    
            L=LTOP
            GOTO 150                                                      
          ELSEIF ( TG(J,NCRB).LT.TMIN ) THEN                              
            L=L-1                                                         
            GOTO 100                                                      
          ENDIF 
          TC(L)=TG(J,L)                                                   
          QC(L)=QG(J,L)                                                   
C     UP ONE LEVEL - DRY ADIABAT AS FIRST GUESS                           
  110     LP=L                                                            
          L=L-1                                                           
          IF(L.LT.LTOP) GOTO 150                                          
          TC(L)=TC(LP)*SK(L)                                              
          QC(L)=QC(LP)
          QSL=ESCON(L)*PQSAT(TC(L))                                       
C     IF SATURATED MAY BE MOIST CONVECTION                                
          IF(QC(L).GE.QSL) GOTO 120                                       
C     TEST FOR DRY INSTABILITY (SHOULD BE STABLE)                         
          GOTO 150  !No moist convection                                  
C     POSSIBLE MOIST CONVECTION - ITERATE FOR THETAE (T,Q)                
  120     NIT=0                                                           
  130     NIT=NIT+1                                                       
          DQN=QC(L)-QSL                                                   
          DQ=DQN/(1.0+CCC*QSL/(TC(L)*TC(L)))                              
          QC(L)=QC(L)-DQ                                                  
          TC(L)=TC(L)+CTQ*DQ                                              
          IF(ABS(DQN).LT.EPSIQ) GOTO 140                                  
          IF(NIT.GT.10) GOTO 140                                          
          QSL=ESCON(L)*PQSAT(TC(L))                                       
          GOTO 130                                                        
C     BUOYANCY TEST - IF STABLE GO TO 150                                 
  140     IF(TC(L).LT.TG(J,L)) GOTO 150                                   
C     MOIST CONVECTION - CONTINUE PARCEL CURVE UP                         
          LCWET=.TRUE.                                                    
          GOTO 110                                                        
  150     CONTINUE                                                        
C     STABLE LAYER OR MODEL TOP - ADJUST ANY UNSTABLE LAYER BELOW         
          IF(LCWET) THEN                                                  
            NCRT=LP    

            IF(NCRT.LT.NCUTOP) THEN                                       
              IF(.NOT.LCBADJ) CALL CBCON(NCRB,NCRT,J,IHEM)                
              IF(LCBADJ)      CALL CBADJ(NCRB,NCRT,J,IHEM)                
            ELSE                                                          
              NCT=NCRT                                                    
              IF(.NOT.LCUBM)  CALL CUDIF(NCRB,NCT,J,IHEM)                 
              IF(LCUBM)       CALL CUBM (NCRB,NCT,J,IHEM)                 
            ENDIF                                                         
          ENDIF                                                           
          IF(L.GT.LTOP) GOTO 100                                          
C                                                                         
C     COLLECT TENDENCIES AND INCREMENT T,Q FOR MOIST CONVECTION           
C                                                                         
          DO 200 L=1,NL                                                   
            TTCR(J,L)=TTDC(L)+TTMC(L)                                     
            QTCR(J,L)=QTDC(L)+QTMC(L)+QTLR(J,L)                           
            QTLR(J,L)=0.0                                                 
  200     CONTINUE                                                        
  700   CONTINUE                                                          
  800 CONTINUE                                                            
C                                                                         
 900  CONTINUE
      RETURN                                                              
      END                                                                 
