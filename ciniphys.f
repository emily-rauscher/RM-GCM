C**********************************************************               
C             SUBROUTINE INIPHYS                                          
C**********************************************************               
      SUBROUTINE INIPHYS                                                  
C                                                                         
C     Sets up PHYSICS variables and arrays. Sets NAMELIST variables       
C     to their default settings, then reads NAMELIST                      
C                 (Piers Forster, 23/10/96)                               
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
C     Legendre polynomials and information about gaussian latitudes       
C                                                                         
      COMMON/LEGAU/ ALPJ(MJP),DALPJ(MJP)                                  
     +              ,ALP(NWJ2,2,JGL),DALP(NWJ2,2,JGL)                     
     +              ,RLP(NWJ2,2,JGL),RDLP(NWJ2,2,JGL)                     
     +              ,SI(JGG),CS(JGG),SISQ(JGG),CSSQ(JGG),SECSQ(JGG)       
     +              ,ALAT(JGG),GWT(JGG),AW(JGG),JH,JL,JINC                
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
C physics namelist                                                        
C                                                                         
      NAMELIST/INPHYS/LBL,LVD,LCR,LLR,LRD,LCUBM,LCBADJ                    
     :,CD,BLRH,AKVV,AKTV,AKQV,AKTC,AKQC,NLCR,NCUTOP,FW,FRADC              
     :,CURHM,CUBMT,CBADJT,CBADJP                                          
     $     ,LSL,LOC,LNOICE,LOLDBL,LCOND,LNNSK,ITSLL,ITSLO                 
C                                                                         
      LSL=.TRUE.                                                          
      LOC=.TRUE.                                                          
      LNOICE=.FALSE.                                                      
      LOLDBL=.FALSE.                                                      
      LCOND=.TRUE.                                                        
      LNNSK=.TRUE.                                                        
      ITSLL=6.                                                            
      ITSLO=6.                                                            
C                                                                         
  240 FORMAT(' PHYSICAL PROCESSES INCLUDED ARE AS FOLLOWS:')              
  241 FORMAT(' SURFACE FLUXES : CD, RH-FACTOR =',F10.4,F10.1,'%')         
  242 FORMAT(' CONVECTION THROUGH',I4,' LAYERS FROM LOWEST LEVEL'/        
     :' DRY CONVECTION   : ADJUSTMENT TO NEUTRALITY IN SINGLE TIMESTEP')  
  246 FORMAT(' PRECIPITATING CONVECTION : NON-ENTRAINING CLOUD MODEL'     
     :,' WITH ENVIRONMENTAL SDOT =',E8.1)                                 
  247 FORMAT(' PRECIPITATING CONVECTION : ADJUSTMENT TO SUBSATURATED'     
     :,' MOIST ADIABAT WITH SCRIPT-P(MB) TAU(HOURS) =',2F6.1)             
  248 FORMAT(' NON-PRECIPITATING DIFFUSIVE MOIST CONVECTION WITH'         
     :,' BACKGROUND KT KQ =',2F6.1,' M2/S')                               
  249 FORMAT(' NON-PRECIPITATING BETTS-MILLER MOIST CONVECTION WITH'      
     :,' BACKGROUND TIMESCALE ',F6.1,' HOURS')                            
  250 FORMAT(' FORMALLY NON-PRECIPITATING CONVECTION INCLUDES RAINOUT'    
     :,' EXCESS OVER',F6.1,'% RH')                                        
  243 FORMAT(' LARGE-SCALE CONDENSATION TO SATURATION IN TIMESTEP')       
  244 FORMAT(' RADIATION: MORCRETTE CODE(SMR+PMF), DIAGNOSTIC CLOUD       
     : VERSION 2.0')                                                      
  245 FORMAT(' VERTICAL DIFFUSION OF MOMENTUM HEAT AND MOISTURE WITH',    
     :       ' CONSTANT COEFFS AKVV AKTV AKQV =',3F8.2)                   
C                                                                         
C     Preset namelist values.                                             
C                                                                         
      LBL=.FALSE.                                                         
      LVD=.FALSE.                                                         
      LCR=.FALSE.                                                         
      LLR=.FALSE.                                                         
      LRD=.FALSE.                                                         
      LCUBM=.FALSE.                                                       
      LCBADJ=.FALSE.                                                      
      CD=0.001                                                            
      BLRH=100.0                                                          
      AKVV=1.0                                                            
      AKTV=1.0                                                            
      AKQV=1.0                                                            
      AKTC=10.0                                                           
      AKQC=10.0                                                           
      FW=1.0E-6                                                           
      FRADC=1.25                                                          
      CURHM=200.0                                                         
      CUBMT=3.0                                                           
      CBADJT=3.0                                                          
      CBADJP=-30.0                                                        
                                                                          
C     INITIALISE PARAMETERS FOR PHYSICAL PROCESSES                        
C                                                                         
      NLCR=NLM
      NCNT=0                                                              
      DO 188 L=1,NL                                                       
       IF(SIGMA(L).LT.0.75) GOTO 188                                      
       NCNT=NCNT+1                                                        
  188 CONTINUE                                                            
      NCUTOP=NLP-NCNT                                                     
C      write(*,*) 'reading fort.7 from iniphys'
      READ (7,INPHYS)                                                     
      WRITE(2,INPHYS)                                                     
      WRITE(2,240)                                                        
      IF(LBL) WRITE(2,241)CD,BLRH                                         
      IF(LCR) WRITE(2,242)NLCR                                            
      IF(LCR.AND..NOT.LCBADJ)WRITE(2,246)FW                               
      IF(LCR.AND.LCBADJ)     WRITE(2,247)CBADJP,CBADJT                    
      IF(LCR.AND..NOT.LCUBM) WRITE(2,248)AKTC,AKQC                        
      IF(LCR.AND.LCUBM)      WRITE(2,249)CUBMT                            
      IF(LCR.AND.(CURHM.LT.100.01))WRITE(2,250)CURHM                      
      IF(LLR) WRITE(2,243)                                                
      IF(LRD) WRITE(2,244)                                                
      IF(LVD) WRITE(2,245)AKVV,AKTV,AKQV                                  
      CLATNTI=2.834E6                                                     
      SFAC=0.5**KITS                                                      
      IF(LRSTRT.AND..NOT.LSHORT) SFAC=1.0                                 
      CTQ=CLATNT/(CPD*CT)                                                 
      CTQI=CLATNTI/(CPD*CT)                                               
      if (LNOICE) CTQI=CTQ                                                
      CCC=CLATNT*CLATNT/(RV*CPD*CT*CT)                                    
      ESCONB=CLATNT/RV                                                    
      ESCONA=RD*EXP(ESCONB/273.15)*610.7/(RV*P0)                          
      ESCONB=ESCONB/CT                                                    
      EPSIQ=0.01/CQ                                                       
C     Used to be multiplied by old (constant) CD. Now done in BLAYER      
      DRAG=GA*RADEA/(RD*CT*DSIGMA(NL))                                    
      BLVAD=3.0/CV                                                        
      BLA=500.0/CV ! A=500m/s                                             
      BLRH=BLRH/100.0                                                     
      FW=FW/WW                                                            
      DTBUOY=1.0E-5/CT                                                    
      TSLA=55.0/CT                                                        
      TSLB=2840.0/CT                                                      
      TSLC=3.5                                                            
      TSLD=3.5*LOG(CT)-LOG(P0)-0.67485                                    
      CBADJT=CBADJT*3600.0*WW                                             
      CBADJP=CBADJP*100.0/P0                                              
      CURHM=CURHM/100.0                                                   
      CUBMT=CUBMT*3600.0*WW                                               
C     Used to be divided by old CD, now done in CUDIF and CUBM            
      CUT1=1.0E4*DSIGMA(NL)/(RADEA*0.2)                                   
      CUT2=1.0E4*DSIGMA(NLM)/RADEA                                        
      CUT2=CUT2*CUT2*0.5                                                  
      AKTC=AKTC/(RADEA*CV)                                                
      AKQC=AKQC/(RADEA*CV)                                                
      CDIF=GA/(CV*WW)                                                     
      CDIF=CDIF*CDIF                                                      
      RCON=P0*DELT/(SFAC*GA)                                              
      CCR=RCON*FW                                                         
C     Initialise values for the soil model                                
      sdsnd=300.             ! density of snow (kg/m3) dimensional        
      sdsn=sdsnd*(RD*CT)/P0  ! density of snow dedimensionalised          
      sdw=1000.*(RD*CT)/P0   ! density of water (kg/m3)                   
      shcs=2085.0E3*CT/P0    ! rho*c for 50% saturated soil               
      shcsp=15.312E6*CT/P0   ! equ. heat cap. for freezing of soil water  
      shcsn=627.0E3*CT/P0    ! rho*c for snow                             
      skse=1.072818*CT/(RADEA*CV*P0)  ! effective soil conductivity       
      sksn=0.24*CT/(RADEA*CV*P0)  ! thermal conductivity for snow         
      slhf=3.5E5/(RD*CT)     ! latent heat of fusion of ice               
      sd1=0.06/RADEA         ! depth of upper soil level                  
      sd2=2.2999/RADEA       ! depth of lower soil level                  
      ssmc=0.5/RADEA         ! depth equivalent to capacity of model      
C-------------------------------                                          
C     Albedos for snow and sea ice, depth for snow albedo function        
C-------------------------------                                          
      sasnow=0.8 ! deep snow albedo                                       
      saice=0.6  ! ice albedo                                             
      shsstar=0.30/RADEA  ! depth of snow for albedo transfer function    
C     max depth of snow for heat capacity of top soil layer               
C     Note - equal to e-folding depth/sqrt(2)                             
      shsmax=1.4/RADEA                                                    
      shco=25.*1000.*4190.*CT/P0/RADEA ! heat capacity of ocean mixed     
                                       ! layer. Mixed layer depth 25m.    
      write(2,*)'Mixed Layer Depth is 25m'                                
      shci=2.*900.*2100.*CT/P0/RADEA                                      
      DO 190 L=1,NLM                                                      
        CLR(L)=RCON*DSIGMA(L)                                             
        FWS(L) =-FW*RDSIG(L)                                              
        FB(L)=CDIF*4.0*SIGMAH(L)*SIGMAH(L)/(SIGMA(L+1)-SIGMA(L))          
        SKAP(L)=SIGMA(L)**AKAP                                            
  190 SK(L)  =(SIGMA(L)/SIGMA(L+1))**AKAP                                 
      CLR(NL)=RCON*DSIGMA(NL)                                             
      FWS(NL)=-FW*RDSIG(NL)                                               
      SKAP(NL)=SIGMA(NL)**AKAP                                            
      FR=180.0/PI                                                         
      FC=FRADC/(PI2*CT)                                                   
      DO 193 IHEM=1,NHEM                                                  
      DO 193 J=1,JG  ! Cooling rate now identical at all latitudes.       
  193 FRAD(J,IHEM)=-FC                                                    
      AKVV=AKVV/(RADEA*CV)                                                
      AKTV=AKTV/(RADEA*CV)                                                
      AKQV=AKQV/(RADEA*CV)                                                
      NAVRD=20                                                            
      NAVWT=21                                                            
C                                                                         
      RETURN                                                              
      END                                                                 
