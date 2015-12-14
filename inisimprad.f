C**********************************************************               
C             SUBROUTINE INISIMPRAD                                          
C**********************************************************               
      subroutine INISIMPRAD                                                 
C-----------------------------------------------------------------------  
C     Subroutine to initialise the simplified radiation scheme             
C-----------------------------------------------------------------------  
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
C     Constant arrays and variables associated with time and vertical     
C     differencing. Also counters.                                        
C                                                                         
      COMMON/BATS/  BEGDAY,CTRA(NTRAC),BM1(IDE),AK(NNP),AQ(NL2),G(NL2)              
     +              ,TAU(NL2),KOUNT,KITS,KSTART,KTOTAL,KRUN,ITSPD           
     +              ,DELT,DELT2,CV,CG,CT,CQ,PNU,PNU2,PNU21                
     +              ,NTRACO,KOLOUR(NTRAC),RGG(NL2)            
     +              ,BEGDOY,DOY                                           
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
       COMMON/CPIERS/ICFLAG(IGC,5,2),CFRAC(IGC,5),PNET(IGC,JG)            
     :               ,SNET(IGC,JG),RRFLUX(IGC,JG,6)                       
     :               ,TTSW(IGC,NL),TTLW(IGC,NL)                           
C                                                                         
C     Legendre polynomials and information about gaussian latitudes       
C                                                                         
      COMMON/LEGAU/ ALPJ(MJP),DALPJ(MJP)                                  
     +              ,ALP(NWJ2,2,JGL),DALP(NWJ2,2,JGL)                     
     +              ,RLP(NWJ2,2,JGL),RDLP(NWJ2,2,JGL)                     
     +              ,SI(JGG),CS(JGG),SISQ(JGG),CSSQ(JGG),SECSQ(JGG)       
     +              ,ALAT(JGG),GWT(JGG),AW(JGG),JH,JL,JINC                
C                                                                         
       COMMON/GSG/GSG(IGC,JG)                                             
C                
       COMMON/SIMPIRRAD/LLOGPLEV,LFLUXDIAG,L1DZENITH,LDIUR,ABSMMRLW,
     & JSKIPLON,JSKIPLAT,NSWMODEL,NLWMODEL,ABSSW1,ABSSTRAT,PRMIN,ALBSW1,
     & ABSSW2,SCATSW2,ASYMSW2,ABSLW1,NEWTB,NEWTE
       LOGICAL LLOGPLEV,LFLUXDIAG,L1DZENITH,LDIUR

       NAMELIST/INSIMPRAD/LLOGPLEV,LFLUXDIAG,L1DZENITH,LDIUR,ABSMMRLW,
     &      JSKIPLON,JSKIPLAT,NSWMODEL,NLWMODEL,ABSSW1,ABSSTRAT,PRMIN,
     &      ALBSW1,ABSSW2,SCATSW2,ASYMSW2,ABSLW1,NEWTB,NEWTE

C Switch to enable proper sub-layer calculation in LW scheme given log(P) 
C distribution of vertical levels
       LPLOTLEV=.FALSE.
C Switch to output in fort.63 the various flux and heating rates diagnostics
       LFLUXDIAG=.TRUE.
C Switch to force irradiation at zenith around sphere (1D equivalent run)
       L1DZENITH=.FALSE.
C ER modif: switch to use diurnally averaged forcing instead of hemispheric
       LDIUR=.FALSE.

C Value of Mass Mixing Ratio (MMR) of absorber in LW scheme
C Should be 1.0 but allows overall scaling up/down of asborber amount  
       ABSMMRLW=1.0

C Indices to be used to skip RT calculation every JSKIP index in LON and LAT 
       JSKIPLON=1
       JSKIPLAT=1

C Indices to select the SW and LW models used in RT calculation 
C NSWMODEL=1 is SW clear sky with surface reflection specified by albedo
C NSWMODEL=2 is SW semi-infinite scattering atmosphere (e.g Schneider & Liu 08)
C NLWMODEL=1 is LW unique absorption coeff (*Press) like Schneider & Liu 08
C NLWMODEL=2 is LW with Planck opacities from Freedman, Lodders & Marley 08  
       NSWMODEL=1   
       NLWMODEL=1

C Values of parameters entering the above models
       ABSSW1=0.3  ! optical thickness at surface (* Press)
       ALBSW1=0.0   ! 0.0 for full absorption at bottom
       ABSSW2=0.3  ! optical thickness at surface (* Press)
       SCATSW2= 0.8 ! 0.8 in Schneider & Liu 08 for Jupiter
       ASYMSW2= 0.204 ! 0.204 in Schneider & Liu 08 for Jupiter
       ABSLW1= 3e-5 !in cm^2/g (* Press, hence  *Press^2 for optical depth)
       NEWTB=0
       NEWTE=0

  240 FORMAT(' SIMPLIFIED RADIATIVE SCHEME EMPLOYED:')      
        
  241 FORMAT(' LW SUB-LAYERS TREATED ASSUMING -LOG- PRESSURE LEVELS')         
  242 FORMAT(' LW SUB-LAYERS TREATED ASSUMING -LINEAR- PRESSURE LEVELS')
  243 FORMAT(' SW MODEL 1: CLEAR SKY + SURFACE REFLECTION'     
     :,' ALBEDO ALBSW1=',f6.1)  
  244  FORMAT(' SW MODEL 2: SEMI-INFINITE SCATTERING ATMOSPHERE WITH'     
     :,' TOTAL OPTICAL THICKNESS ABSSW2=',f6.1
     :,', SCATTERING COEFF SCATSW2= ',f4.1 
     :,', ASYMMETRY FACTOR ASYMSW2= ',f4.1)  
  245 FORMAT(' LW MODEL 1: SINGLE ABSORPT. COEFF. ABSLW1= ',e9.2)     
  246 FORMAT(' LW MODEL 2: PLANCK MEAN OPACITY TABLE') 
 247  FORMAT(' NSWMODEL INDEX INAPPROPRIATE:',I4)
 248  FORMAT(' NLWMODEL INDEX INAPPROPRIATE:',I4)
 249  FORMAT('Newtonian heating only on levels: ',I4,' to ',I4)

C      write(*,*) 'reading fort.7 in inisimprad'
      READ (7,INSIMPRAD)                                                     
      WRITE(2,INSIMPRAD)                                                     

      WRITE(2,240)                                                        
      IF(LLOGPLEV) THEN 
         WRITE(2,241)                                        
      ELSE 
         WRITE(2,242) 
      ENDIF

      IF(NSWMODEL.EQ.1) THEN

C Scale bottom of atmosphere optical depth for RADSW, using P0 and GA
C Factor of 10 to scale ABSSW1 from CGS to code units (like ABSLW1)
         ABSSW1=P0/GA*ABSSW1/10.0
         ABSSTRAT=P0/GA*ABSSTRAT/10.0

         WRITE(2,243) ALBSW1
      ELSE IF(NSWMODEL.EQ.2) THEN
         WRITE(2,244) ABSSW2, SCATSW2, ASYMSW2
      ELSE
      WRITE(2,247) NSWMODEL
      CALL ABORT
      ENDIF
 
      IF(NLWMODEL.EQ.1) THEN
         WRITE(2,245) ABSLW1
      ELSE IF(NLWMODEL.EQ.2) THEN
         WRITE(2,246)
      ELSE
      WRITE(2,248) NLWMODEL
      CALL ABORT
      ENDIF
    
      IF(.NOT.(NEWTB.EQ.NEWTE)) WRITE(2,249),NEWTB,NEWTE
                                                                          
      END                                                                 
