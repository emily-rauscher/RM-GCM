C**********************************************************               
C             SUBROUTINE BLAYER                                           
C**********************************************************               
      SUBROUTINE BLAYER                                                   
C                                                                         
C     SURFACE LAYER FLUXES FOR PRESENT LATITUDE.                          
C     BULK AERODYNAMIC FORMULATION FOR LOWEST LEVEL TENDENCIES.           
C                                                                         
C                                                                         
C     Determines model resolution                                         
C                                                                         
      include 'params.i'
                                                                          
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
C     Setup moisture variables by equivilencing them to                   
C     Tracer No. 1                                                        
C                                                                         
      REAL QG(IGC,NL),QNLG(IGC,NL),QTLR(IGC,NL),QTCR(IGC,NL)              
      EQUIVALENCE (QG(1,1),TRAG(1,1,1)) , (QNLG(1,1),TRANLG(1,1,1))       
     *            ,(QTLR(1,1),UTRAG(1,1,1)),(QTCR(1,1),VTRAG(1,1,1))      
C                                                                         
C                                                                         
      RCSJ=SECSQ(JH)                                                      
      SQRC=SQRT(RCSJ)                                                     
      DO 800 IHEM=1,NHEM                                                  
        IOFM=(IHEM-1)*MGPP                                                
        DO 10 I=1,MG                                                      
          J=I+IOFM                                                        
          BLCD(J)=CD                                                      
          CTAU=CD/TSTAR(J,JH)                                             
          CSH =CTAU/AKAP                                                  
          DRAGJ=CD*DRAG/TSTAR(J,JH)                                       
                                                                          
          VM=SQRT(RCSJ*(UG(J,NL)*UG(J,NL)+VG(J,NL)*VG(J,NL)))             
          VMP=VM+BLVAD                                                    
          CVM=DRAGJ*VMP                                                   
          UNLG(J,NL)=-CVM*UG(J,NL)                                        
          VNLG(J,NL)=-CVM*VG(J,NL)                                        
          CVM=CTAU*PLG(J)*VMP*SQRC                                        
          TXBL(J)=CVM*UG(J,NL)                                            
          TYBL(J)=CVM*VG(J,NL)                                            
          THNL=TG(J,NL)/SKAP(NL)                                          
C     The next line prevents condensation onto a warm surface,            
C     which doesn't happen in reality due to the surface                  
C     effectively becoming saturated after a very small amount            
C     of condensation has occurred                                        
          QSURF=QSTAR(J,JH)                                               
          IF (.NOT.LCOND) THEN                                            
            QSURF=MIN(SQSTAR(J,JH),MAX(QSTAR(J,JH),QG(J,NL)))             
          ENDIF                                                           
          DTH=TSTAR(J,JH)-THNL                                            
          DQ=BLRH*(QSURF-QG(J,NL))                                        
          BLVB(J)=BLA*SQRT(ABS(DTH)/(0.5*(THNL+TSTAR(J,JH))))             
          IF (DTH.LT.0.0) THEN                                            
            BLVB(J)=0.0                                                   
          END IF                                                          
          VMP=VM+BLVB(J)                                                  
          CVM=0.2*DRAGJ*VMP                                               
          TNLG(J,NL)=CVM*DTH                                              
          QNLG(J,NL)=CVM*DQ                                               
          CVM=0.2*CSH*PLG(J)*VMP                                          
          SHBL(J)=CVM*DTH                                                 
C     Choose appropriate value of latent heat based on ground T           
c     ie, we assume it snows if the ground is colder than 273.15K         
          IF (tstar(j,jh).gt.0.363224029) THEN                            
           slbl(j)=cvm*ctq*dq                                             
          ELSE                                                            
           slbl(j)=cvm*ctqi*dq                                            
          ENDIF                                                           
          UG(J,NL)=UG(J,NL)+DELT2C*UNLG(J,NL)                             
          VG(J,NL)=VG(J,NL)+DELT2C*VNLG(J,NL)                             
          TG(J,NL)=TG(J,NL)+DELT2C*TNLG(J,NL)                             
          QG(J,NL)=QG(J,NL)+DELT2C*QNLG(J,NL)                             
   10   CONTINUE                                                          
  800 CONTINUE                                                            
      RETURN                                                              
      END                                                                 
