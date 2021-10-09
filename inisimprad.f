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
      REAL AVG,LO,MWTOT,RUNIV,TAU_KMAMGALL,KMAMGperBAR,WVO,TAU_KMAMGH2
     +,mwh2,mwhe,mwh2o,rfhe,rfh2o,fh2,RAYPERBARCONS                                                                   
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
       COMMON/SIMPIRRAD/LLOGPLEV,LFLUXDIAG,L1DZENITH,LDIUR,
     & JSKIPLON,JSKIPLAT, DOSWRAD, DOLWRAD, LWSCAT,
     & FLXLIMDIF,SURFEMIS, RAYSCAT, RAYSCATLAM, AEROSOLS,ABSSW, ABSLW,
     & ALBSW, NEWTB, NEWTE,RAYPERBARCONS, TEMPERATURE_INTERNAL, TEMPERATURE_IRRAD

       REAL SURFEMIS,RAYSCATLAM,ABSSW,ABSLW,ALBSW
       LOGICAL LLOGPLEV,LFLUXDIAG,L1DZENITH,LDIUR,DOSWRAD,DOLWRAD
     + ,LWSCAT, FLXLIMDIF, RAYSCAT,AEROSOLS
       
       NAMELIST/INSIMPRAD/LLOGPLEV,LFLUXDIAG,L1DZENITH,LDIUR,
     & JSKIPLON,JSKIPLAT, DOSWRAD, DOLWRAD, LWSCAT,FLXLIMDIF,SURFEMIS,
     & RAYSCAT, RAYSCATLAM,AEROSOLS,ABSSW,ABSLW, ALBSW, NEWTB, NEWTE, TEMPERATURE_INTERNAL,
     & TEMPERATURE_IRRAD

c
C Switch to enable proper sub-layer calculation in LW scheme given log(P) 
C distribution of vertical levels
       LPLOTLEV=.FALSE.
C Switch to output in fort.63 the various flux and heating rates diagnostics
       LFLUXDIAG=.TRUE.
C Switch to force irradiation at zenith around sphere (1D equivalent run)
       L1DZENITH=.FALSE.
C ER modif: switch to use diurnally averaged forcing instead of hemispheric
       LDIUR=.FALSE.

C Indices to be used to skip RT calculation every JSKIP index in LON and LAT 
       JSKIPLON=1
       JSKIPLAT=1

C The following are related to the flags in the radiative transfer suite
!     DOSWRAD= ISL        - do solar calculations when = 1
!     DOLWRAD= IR        - do infrared calculations when = 1
!     LWSCAT= IRS        - do infrared scattering when = 1
!     FLXLIMDIF  - do flux limited diffusion correction = 1 ~MTR
!     SURFEMIS=  EMISIR     - SURFACE IR EMISSIVITY
      DOSWRAD      = .TRUE.
      DOLWRAD      = .TRUE.
      LWSCAT       = .TRUE.
      FLXLIMDIF    = .TRUE.
      SURFEMIS     = 1.

!     RAYSCAT    - include rayleigh scattering
!     RAYSCATWL  - wavelength at which to compute scattering optical
!     depth in double gray regime
      RAYSCAT      =.FALSE.
      RAYSCATLAM    = 0.      
      
!     INCLUDE AEROSOLS 
      AEROSOLS  =.FALSE.
      AEROSOLMODEL ='GENERATE'
      PI0AERSW = 1.
      ASYMSW = 1.
      EXTFACTLW= .1
      PI0AERLW= .8
      ASYMLW= 0.01 
!     GASEOUS ABSORPTION
      ABSSW =1.0
      ABSLW =1.0      
      ALBSW =0.0   ! 0.0 for full absorption at bottom

        LLOGPLEV        = .FALSE.
        LFLUXDIAG       = .TRUE.
        L1DZENITH       = .FALSE.
        LDIUR           = .FALSE.
        JSKIPLON        = 1
        JSKIPLAT        = 1
        DOSWRAD        = .TRUE.
        DOLWRAD         = .TRUE.
        LWSCAT          = .TRUE.
        FLXLIMDIF       = .TRUE.
        EMISIR          = 1.
        RAYSCAT         = .FALSE.
        RAYSCATLAM      = 0.
        AEROSOLS        = .FALSE.
        ABSSW           = 8.14e-4
        ABSLW           = 1E-2
        ALBSW           = 0.0
        NEWTB           = 0
        NEWTE           = 0


c       ABSSW2=0.3  ! optical thickness at surface (* Press)
c       SCATSW2= 0.8 ! 0.8 in Schneider & Liu 08 for Jupiter
c       ASYMSW2= 0.204 ! 0.204 in Schneider & Liu 08 for Jupiter
!MTR       ABSLW1= 3e-5 !in cm^2/g (* Press, hence  *Press^2 for optical depth)
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
! 249  FORMAT('Newtonian heating only on levels: ',I4,' to ',I4)
 250  FORMAT('MATRIX SW FROM TOON ET AL. 1989: '
     &,'SW ABSCOEFF AT REFERENCE PRESSURE=',f10.5)
 251  FORMAT('MATRIX LW & SOURCE FUNCTION FROM TOON ET AL. 1989: '
     &,'LW ABSCOEFF AT REFERENCE PRESSURE=',f10.5) 

C      write(*,*) 'reading fort.7 in inisimprad'
      READ (7,INSIMPRAD)                                                     
C      WRITE(2,INSIMPRAD)                                                     
      WRITE(2,240)                                                        
      IF(LLOGPLEV) THEN 
         WRITE(2,241)                                        
      ELSE 
         WRITE(2,242) 
      ENDIF

!      IF(NSWMODEL.EQ.1) THEN

C Scale bottom of atmosphere optical depth for RADSW, using P0 and GA
C Factor of 10 to scale ABSSW1 from CGS to code units (like ABSLW1)
!MTR         ABSSW1=P0/GA*ABSSW1/10.0
!MTR         ABSSTRAT=P0/GA*ABSSTRAT/10.0

!         WRITE(2,243) ALBSW1
!      ELSE IF(NSWMODEL.EQ.2) THEN
!         WRITE(2,244) ABSSW2, SCATSW2, ASYMSW2
!      ELSE IF(NSWMODEL.EQ.3) THEN
         WRITE(2,250) ABSSW
!      ELSE
!      WRITE(2,247) NSWMODEL
!      CALL ABORT
!      ENDIF
 
!      IF(NLWMODEL.EQ.1) THEN
!         WRITE(2,245) ABSLW1
!      ELSE IF(NLWMODEL.EQ.2) THEN
!         WRITE(2,246)
!      ELSE IF(NLWMODEL.EQ.3) THEN
         WRITE(2,251) ABSLW
         
!      ELSE
!      WRITE(2,248) NLWMODEL
!      CALL ABORT
!      ENDIF
!      COMPUTER KM-AMG PER BAR OF PRESSURE
       AVG          = 6.0221409E23 !#
       LO           = 2.686763E25  !#/m^3
       RUNIV        = 8.3143  !J/mol/K 
       MWTOT        = RUNIV/GASCON * 1000. !Kg/mol
       KMAMGperBAR  = AVG*1.E5/(LO*GA*MWTOT)
!      Now compute the rayleigh scattering optical depth per km-amg for
!      H2 gas using the expression from Dalgarno and Williams (ApJ ,
!      1962).
       WVO          = RAYSCATLAM*1E4 ! ANGSTROMS
       TAU_KMAMGH2  = 2.687*(8.14E11/WVO/WVO/WVO/WVO 
     &   +1.28E18/WVO/WVO/WVO/WVO/WVO/WVO 
     &   +1.61E24/WVO/WVO/WVO/WVO/WVO/WVO/WVO/WVO)
!      However, this is just for H2 gas.  The molecular weight chosen by
!      choice of gas constant implies a gas with heavier elements, so it
!      is not technically self consistent.  To make it self consistent,
!      a fraction of heavier elements can be computed from the molecular
!      weight if we choose some values. 
!      Assumption: Let's arbitrarily assume the remainder (not H2)
!      is 90% HE and 10% H2O (though one could just as arbitrarily use
!      CH4 to increase the rayleigh scattering given CH4's high index of
!      refraction); then the resulting H2 fraction is:
       mwh2=2.0
       mwhe=4.0
       mwh2o=18.0 
       rfhe=0.8
       rfh2o=0.2
       fh2=(MWTOT-rfhe*mwhe-rfh2o*mwh2o)/(mwh2-rfhe*mwhe-rfh2o*mwh2o)
!      Each molecular component contributes its mole fraction x a
!      product that depends on its index of refraction and that of
!      hydrogen:  indexfact= (n-1)**2 / (nh2 -1)**2
!      for HE, indexfact = 0.0641; for H2O = 3.3690; for H2 = 1.0,
!      CH4=10.1509,etc...
       TAU_KMAMGALL =TAU_KMAMGH2* 
     &   (fh2*1.0 + (1.-fh2)*rfhe*0.0641 + (1.-fh2)*rfh2o*3.3690)
!      Finally, the rayleight scattering optical depth per bar is...
       RAYPERBARCONS = TAU_KMAMGH2 * KMAMGperBAR
!      NOW WRITE SOME GENERAL RADIATIVE TRANSFER SCHEME DETAILS TO FILE 
       WRITE(60,*) 'RADIATIVE TRANSFER SUMMARY'
       WRITE(60,*) ''
       WRITE(60,*) 'Two-stream, plane-parallel columns' 
       write(60,*) 'Quadrature closure in SW,'
       write(60,*) 'Hemispheric Mean followed by source function in LW'
       write(60,*) '(Toon et al.,1989). Potential scattering and' 
       write(60,*) 'absorption in SW and LW by gas and aerosols.'
       write(60,*) '(by default, the gas does not scatter in longwave)'
     
       write(60,*) '' 
       WRITE(60,*) 'GAS TRANSMISSIONS'
       write(60,*) ' kappa SW: ',ABSSW,'(cm^2/g)'
       write(60,*) ' kappa LW: ',ABSLW,'(cm^2/g)'
       write(60,*) ' tau per bar (absorption):' 
       write(60,*) '       SW: ',ABSSW*1e6/GA/100.
       write(60,*) '       LW: ',ABSLW*1e6/GA/100.
       WRITE(60,*)'Pressure (bars) of SW clear gas tau = 2/3 @ mu0=1',
     &               (2./3.)/(ABSSW*1e6/GA/100.)

       WRITE(60,*)'Pressure (bars) of LW clear gas tau = 2/3 @ mu0=1',
     &               (2./3.)/(ABSLW*1e6/GA/100.)

       WRITE(60,*) 'RAYLEIGH SCATTERING'
       WRITE(60,*) 'RAYLEIGH OPTICAL DEPTH PER BAR: ',RAYPERBARCONS 
       WRITE(60,*) 'KM-AMAGATS OF GAS PER BAR: ',KMAMGperBAR
       WRITE(60,*) ' Pressure of SW Two-way gaseous (inc.ray) tau= 1: ',
     &               .5/(ABSSW*1e6/GA/100. + RAYPERBARCONS) 
       WRITE(60,*) ' N.B. Assumes molecular composition with:'
       write(60,*) '      H2 fraction = ',fh2
       write(60,*) '      He fraction = ',rfhe
       write(60,*) '      H2O fraction= ',rfh2o
       write(60,*) ' ...To match the molecular weight. (see inisimprad)'
       WRITE(60,*) ''
       write(60,*) 'Surface Albedo in SW = ',ALBSW
       write(60,*) 'Surface Emissivity in LW = ',EMISIR
       
       write(60,*) 'Flux limited diffusion?, ',FLXLIMDIF 
       
       write(60,*)'INSIMPRAD/'
       write(60,*)'LLOGPLEV',LLOGPLEV
       write(60,*)'LFLUXDIAG',LFLUXDIAG
       write(60,*)'L1DZENITH',L1DZENITH
       write(60,*)'LDIUR',LDIUR
       write(60,*)'JSKIPLON',JSKIPLON
       write(60,*)'JSKIPLAT',JSKIPLAT
       write(60,*) 'DOSWRAD',DOSWRAD 
       write(60,*)'DOLWRAD',DOLWRAD
       write(60,*)'LWSCAT',LWSCAT
       write(60,*)'FLXLIMDIF',FLXLIMDIF
       write(60,*)'FLD TAULIMIT',TAULIMIT
       write(60,*)'SURFEMIS',SURFEMIS
       write(60,*)'RAYSCAT',RAYSCAT
       write(60,*)'RAYSCATLAM',RAYSCATLAM
       write(60,*)'AEROSOLS',AEROSOLS
       write(60,*)'ABSSW',ABSSW
       write(60,*)'ABSLW',ABSLW
       write(60,*)'ALBSW',ALBSW

       ! MALSKY SHOULD THESE LINES BE HERE????
       write(60,*)'TEMPERATURE_INTERNAL', TEMPERATURE_INTERNAL
       write(60,*)'TEMPERATURE_IRRAD',TEMPERATURE_IRRAD
      END                                                                 
