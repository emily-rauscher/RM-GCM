C*************************************************************            
C                   SUBROUTINE INISET                                     
C*************************************************************            
      SUBROUTINE INISET                                                   
C                                                                         
C     Sets up various variables and arrays. Sets NAMELIST variables       
C     to their default settings, then reads NAMELIST                      
C                                                                         
C                                                                         
C     Determines model resolution                                         
C                                                                         
      include 'params.i'
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
C     Array ordering in SPECTR must correspond to that in GRIDP.          
C                                                                         
      COMMON/SPECTR/Z(IGB),D(IGB),T(IGB),TRA(IGB,NTRAC),SP(IGA),GS(IGA)   
     :              ,SPA(IGA),VP(IGA),DTE(IGB),TT(IGB)                    
     :              ,TRAT(IGB,NTRAC),DT(IGB),ZT(IGB)                      
     :              ,ZMI(IGB),DMI(IGB),TMI(IGB)                           
     :              ,TRAMI(IGB,NTRAC),SPMI(IGA)                           
      COMPLEX Z,D,T,TRA,SP,GS,SPA,VP,DTE,TT,TRAT,DT,ZT                    
     :       ,ZMI,DMI,TMI,TRAMI,SPMI                                      
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
C     Constants and arrays needed for the fast Fourier transforms         
C                                                                         
      COMMON/COMFFT/NTGW,NRSTGW,NTWG,NRSTWG,NTGD,NRSTGD                   
     +              ,TRIG(IDA),WORK(IDF),IFAX(10)                         
C                                                                         
C                                                                         
C     Polynomial used to aid vectorization of Legendre transforms         
C                                                                         
      COMMON/POLYNO/POLY(NWJ2,2,4),CMPA(IGL)                              
      COMPLEX CMPA                                                        
C                                                                         
C                                                                         
C     Restoration fields and timescale                                    
C                                                                         
      COMMON/RESTOR/ZRES(IGN),DRES(IGN),TRES(IGN),SPRES(IGM),DAMP         
C                                                                         
      COMMON/STATS/GMSP0,GMSPMI,LMASCOR,LMASOLD,LMASPRT

      COMMON/CLOUDY/AEROSOLMODEL,AERTOTTAU,CLOUDBASE,
     &   CLOUDTOP,CLDFRCT,AERHFRAC,PI0AERSW,ASYMSW,EXTFACTLW,PI0AERLW,
     &   ASYMLW,DELTASCALE,SIG_AREA,PHI_LON,TAUAEROSOL,AEROPROF,
     &   MAXTAU,MAXTAULOC,TCON,AEROSOLCOMP,MTLX,METALLICITY,HAZES,PICKET_FENCE_CLOUDS,MOLEF,AERLAYERS
      NAMELIST/INCLOUDY/AEROSOLMODEL,AEROSOLCOMP,MTLX,METALLICITY,HAZES,PICKET_FENCE_CLOUDS,MOLEF,AERLAYERS,
     &  AERTOTTAU,CLOUDBASE,CLOUDTOP,AERHFRAC,PI0AERSW,
     &  ASYMSW,EXTFACTLW,PI0AERLW,ASYMLW,DELTASCALE,SIG_AREA,PHI_LON
      CHARACTER(30) :: AEROSOLMODEL
      CHARACTER(30) :: AEROSOLCOMP
      REAL TAUAEROSOL(nl+1,mg,2,jg),AEROPROF(NL+1),MAXTAU,TCON(NL+1)
      REAL MOLEF(13)
      REAL MTLX, METALLICITY
      INTEGER AERLAYERS
      LOGICAL DELTASCALE, HAZES, PICKET_FENCE_CLOUDS

      LOGICAL LMASCOR,LMASOLD,LMASPRT                                     
C                                                                         
C-----------------------------------------------------------------------  
C                                                                         
      REAL*8 AKK1                                                         
      NAMELIST/INPPL/ GA,GASCON,RADEA,AKAP,WW,P0,RV,CPD,CLATNT            
      NAMELIST/INPRN/ KRUN,BEGDAY,TSPD,KITS,PNU,TDISS                     
     +         ,LFLUX, BEGDOY                                             
     + ,NDEL,T0,LRSTRT,LSTRETCH,LSHORT,LTVEC,LBALAN,LRESTIJ               
     +   ,LCLIM, LPERPET, L22L,LOROG,LCSFCT                               
     + ,KOLOUR,LNOISE,LMASCOR,LMASOLD,LMASPRT                             
      NAMELIST/INPOP/RNTAPE,KOUNTH,KOUNTR,KOUNTP,KOUNTE                   
     + ,NCOEFF,NLAT,LGPO,LSPO                                             
     + ,RNTAPO,NTRACO                                                     
     $     ,LSHIST,LMINIH                                                 
C
C
C                                                                         
  205 FORMAT(/' *****RNTAPE*****',F12.3)                                  
  207 FORMAT(' PRINTED OUTPUT EVERY ',I3,' TIMESTEPS'/                    
     +       ' RMS QUANTITIES OUTPUT EVERY ',I3,' TIMESTEPS'/             
     +       ' HISTORY RECORD WRITTEN EVERY ',I3,' TIMESTEPS'/            
     +       ' RESTART RECORD WRITTEN EVERY ',I3,' TIMESTEPS')            
  209 FORMAT(' INTEGRATION WITH',I3,                                      
     +' LEVELS IN THE VERTICAL (NL=',I3,')'/                              
     +' JAGGED TRIANGULAR/TRAPEZOIDAL TRUNCATION AT TOTAL WAVENO. ',I3/   
     +' AND ZONAL WAVENO. ',I3,' (NN=',I3,' MM=',I3,')')                  
  221 FORMAT(' NO LATERAL DISSIPATION')                                   
  222 FORMAT(' DEL',I2,' LATERAL DISSIPATION ON VORTICITY, DIVERGENCE'/   
     + ' AND TEMPERATURE WITH DIFFUSION COEFFICIENT  ',E11.4,             
     + ' m**',I1,'/s'/' THE E-FOLDING TIME FOR SMALLEST RESOLVED',        
     + ' SCALE IS ',F6.3,' DAYS')                                         
  223 FORMAT(' NO TIME FILTER')                                           
  224 FORMAT(' ROBERT TIME FILTER WITH PARAMETER PNU ',F5.2)              
  280 FORMAT(' GLOBAL DOMAIN: BOTH EVEN AND ODD COEFFICIENTS INCLUDED',   
     +' (NHEM=',I1,')')                                                   
  281 FORMAT(' HEMISPHERIC DOMAIN: ONLY EVEN DIVERGENCE TEMPERATURE'/     
     +' SURFACE PRESSURE AND ODD VORTICITY COEFFICIENTS INCLUDED',        
     +' (NHEM=',I1,')')                                                   
  210 FORMAT(' ',I2,'-FOLD SYMMETRY IN LONGITUDE IMPOSED AND ONLY'/       
     +' 1 /',I2,' OF THE DOMAIN USED (MOCT=',I2,')')                      
  211 FORMAT(' NON LINEAR TERMS EVALUATED ON GRID OF ',I3,                
     +' GAUSSIAN LATITUDES '/' AND ',I3,                                  
     +' EVENLY SPACED LONGITUDES (JG=',I3,' MG=',I3,')')                  
  212 FORMAT(' ECMWF ANGULAR MOMENTUM CONSERVING VERTICAL SCHEME')        
  213 FORMAT(' GLOBAL MASS CORRECTION SWITCHED ON')                       
  214 FORMAT(' GLOBAL MASS CORRECTION SWITCHED OFF')                      
  225 FORMAT(/' ***ABORT*** CORRECT VALUE OF NWJ2',I5,' VALUE GIVEN',I5)  
  232 FORMAT(/' ***ABORT***  NLAT IS GREATER THAN JG*NHEM')               
  233 FORMAT(/' ***ABORT***NCOEFF IS GREATER THAN NN')                    
  240 FORMAT(/' ***ABORT*** NTRACO IS GREATER THAN NTRAC')                
C                                                                         
C     Set default values and override as desired through NAMELIST input   
C      Mass correction namelist defaults                                  
                                                                          
      LMASCOR=.TRUE.  ! mass correction switched on                       
      LMASOLD=.TRUE.  ! masses read from restart file                     
      LMASPRT=.FALSE. ! mass correction printed every KOUNTP steps        
C                                                                         
C     THE FOLLOWING PARAMETERS ARE ALL DEFINED AS INPUTS IN FORT.7, 
C     WHICH IS READ IN BELOW. ALL THESE VALUES ARE THEREFORE
C     PLACEHOLDERS
                                                                         
      GA=0.0                                                          
      GASCON=0.0                                                        
      RADEA=0.0                                                     
      AKAP=0.0                                                          
      WW=0.0                                                         
      P0=0.0                                                         
      RV=0.0                                                           
      CLATNT=0.0                                                        
      CPD=0.0 !GASCON/AKAP                                                     
C     
      KRUN=0                                                              
      BEGDAY=0.0                                                          
      BEGDOY=0.0                                                          
      TSPD=0.0 !24.0                                                           
      KITS=0.0 !3                                                              
      PNU=0.02                                                            
      TDISS=0.0 !0.25                                                          
      NDEL=0.0  !6                                                              
      DO 17 L=1,NL                                                        
         T0(L)=00.0                                                      
   17 CONTINUE                                                            
      LRSTRT=.FALSE.                                                      
      LSTRETCH =.FALSE.                                                   
      LPERPET=.TRUE.                                                      
      LCLIM=.FALSE.                                                       
      L22L=.TRUE.                                                         
      LOROG=.FALSE.                                                       
      LCSFCT=.TRUE.                                                       
      LSHORT=.FALSE.                                                      
      LTVEC =.TRUE.                                                       
      LFLUX=.TRUE.                                                        
      LBALAN=.FALSE.                                                      
      LRESTIJ=.FALSE.                                                     
      LNOISE=.FALSE.                                                      
C                                                                         
      RNTAPE=0.0                                                          
      RNTAPO=-999.0                                                       
      NTRACO=0                                                            
      KOUNTH=0                                                            
      KOUNTR=0                                                            
      KOUNTP=0                                                            
      KOUNTE=0                                                            
      NCOEFF=0                                                            
      NLAT=MIN(16,JGG)                                                    
      LSHIST=.TRUE.                                                       
      LMINIH=.TRUE.                                                       
      DO 18 I=1,NL                                                        
         LSPO(I)=.FALSE.                                                  
         LGPO(I)=.FALSE.                                                  
   18 CONTINUE                                                            
C                                                                         
C     Read NAMELISTs, overwrite defaults and write them out               
C
C      write(*,*) 'reading fort.7 in iniset'
      READ(7,INPPL)                                                       
      WRITE(2,INPPL)                                                      
C      write(*,*) 'second read in iniset'
      READ(7,INPRN)                                                       
      WRITE(2,INPRN)                                                      
C      write(*,*) 'third read in iniset'
      READ(7,INPOP)                                                       
      IF (RNTAPO.EQ.-999.0) RNTAPO=RNTAPE                                 
      WRITE(2,INPOP)

C                                                                         
C     Set remaining physical constants.                                   
C                                                                         
      RD=GASCON                                                           
      CPD=GASCON/AKAP  !No longer defined explicitly in fort.7                                                                     
C     Write out details of model run                                      
C                                                                         
      WRITE(2,205)RNTAPE                                                  
      WRITE(2,209)NL,NL,NN,MM,NN,MM                                       
      IF(NHEM.EQ.2) WRITE(2,280) NHEM                                     
      IF(NHEM.EQ.1) WRITE(2,281) NHEM                                     
      WRITE(2,210) MOCT,MOCT,MOCT                                         
      WRITE(2,211)JG,MG,JG,MG                                             
      WRITE(2,212)                                                        
      IF (LMASCOR) THEN                                                   
         WRITE(2,213)                                                     
      ELSE                                                                
         WRITE(2,214)                                                     
      ENDIF


C                                                                         
C     Set resolution dependent quantities                                 
C                                                                         
      AMH=MH                                                              
      MF=MM-1                                                             
      MFP=MF+1                                                            
      MFPP=MFP+1                                                          
      NF=NN-1                                                             
      NFP=NF+1                                                            
      NFPP=NFP+1                                                          
      AIOCT=(0.,1.)*MOCT                                                  
C                                                                         
C     Set various spectral limits and coefficients which                  
C     depend on wavenumber                                                
C                                                                         
      NW=1+MF/MOCT                                                        
      NWP=NW+1                                                            
      MJPP=MJP+NW                                                         
      MGP=MG+1                                                            
      MG2=MG/2                                                            
      RMG=1./REAL(MG)                                                     
      JZF=MGPP-NW-NW                                                      
      DO 4 NP=1,NFPP                                                      
         SQ(NP)=NP*(NP-1)                                                 
         SQH(NP)=0.5*SQ(NP)                                               
         IF (NP.GT.1) THEN                                                
            RSQ(NP)=1./SQ(NP)                                             
         ELSE                                                             
            RSQ(1)=0.                                                     
         ENDIF                                                            
    4 CONTINUE                                                            
C                                                                         
C     Compute internal diffusion parameter                                
C                                                                         
      IF(TDISS.EQ.0.0) THEN                                               
         AKK=0.0                                                          
      ELSE                                                                
c                                                                         
c        Slightly more complicated version to cope with r4 precision:     
c        AKK=WW*(RADEA**NDEL)/(2.0*PI*TDISS*((NN*(NN+1))**REAL(NDEL/2)))  
c        AKK=AKK/(WW*(RADEA**NDEL))                                       
c                                                                         
         AKK=NN*(NN+1)                                                    
         AKK=AKK**REAL(NDEL/2)                                            
         AKK=1.0/(2.0*PI*TDISS*AKK)                                       
         AKK1=RADEA                                                       
         AKK1=WW*AKK1**NDEL                                               
         AKK1=AKK*AKK1                                                    
      END IF                                                              
      IF(AKK.EQ.0.0) WRITE(2,221)                                         
      IF(AKK.NE.0.0) WRITE(2,222) NDEL,AKK1,NDEL,TDISS                     
      NDELH=NDEL/2                                                        
      DO 3 NP=1,NNP                                                       
         AK(NP)=AKK*(SQ(NP)**NDELH)                                       
    3 CONTINUE                                                            
C                                                                         
C     Set time variables and counters                                     
C                                                                         
      IF(PNU.EQ.0.0)WRITE(2,223)                                          
      IF(PNU.NE.0.0)WRITE(2,224)PNU                                       
      WRITE(2,207)KOUNTP,KOUNTE,KOUNTH,KOUNTR                             
C                                                                         
      IF(KOUNTP.EQ.0) KOUNTP=-999                                         
      IF(KOUNTE.EQ.0) KOUNTE=-999                                         
      IF(KOUNTH.EQ.0) KOUNTH=-999                                         
      IF(KOUNTR.EQ.0) KOUNTR=-999                                         
      DELT=PI2/TSPD                                                       
      PNU2=PNU+PNU                                                        
      PNU21=1.0-PNU2                                                      
      ITSPD=NINT(TSPD)                                                    
C                                                                         
C     Check variables make sense                                          
C                                                                         
      IF (NLAT.GT.JGG) THEN                                               
         WRITE(2,232)                                                     
         CALL ABORT                                                       
      ENDIF                                                               
      IF (NCOEFF.GT.NN) THEN                                              
         WRITE(2,233)                                                     
         CALL ABORT                                                       
      ENDIF                                                               
      NWJCH=0                                                             
      DO 310 MP=1,MFP,MOCT                                                
         DO 320 JP=MP,NFP,MH                                              
            NWJCH=NWJCH+1                                                 
 320     CONTINUE                                                         
 310  CONTINUE                                                            
      IF (NWJ2.NE.NWJCH) THEN                                             
         WRITE(2,225) NWJCH,NWJ2                                          
         CALL ABORT                                                       
      ENDIF                                                               
      IF (NTRACO.GT.NTRAC) THEN                                           
         WRITE(2,240)                                                     
         CALL ABORT                                                       
      ENDIF                                                               
C                                                                         
C     Set dimensionalising factors                                        
C                                                                         
      EZ=1.0/SQRT(.375)                                                   
      CV=RADEA*WW                                                         
      CG=CV*CV                                                            
      CT=CG/GASCON                                                        
      CQ=1000.0                                                           
      PFAC=0.5*CV*CV*1.0E5/GA                                             
      SQR2=SQRT(2.0)                                                      
      RSQR2=1.0/SQR2                                                      
      EAM1=SQR2/3.                                                        
      EAM2=SQRT(2./45.)                                                   
C                                                                         
C     Make T0 dimensionless                                               
C                                                                         
      DO 61 L=1,NL                                                        
        T0(L)=T0(L)/CT                                                    
   61 CONTINUE                                                            
      DO 30 KK=1,NTRAC                                                    
c Hard wires water to Tracer Number 1                                     
C water has "colour" of 3.                                                
         IF (KK.EQ.1) THEN                                                
            KOLOUR(KK)=3                                                  
c convert from g/Kg                                                       
            CTRA(KK)=CQ                                                   
C         ELSE                                                            
          IF(KOLOUR(KK).EQ.1)  CTRA(KK)=CT                                
          IF(KOLOUR(KK).EQ.2)  CTRA(KK)=1.E6*CT*WW*GA/P0                  
         ENDIF                                                            
   30 CONTINUE                                                            
C                                                                         
C     Set up arrays and variables for use in FFT routines                 
C                                                                         
      NTRWG=NFTWG*NHEM                                                    
      NTRGW=NFTGW*NHEM                                                    
      NTRGD=NFTGD*NHEM                                                    
      NTWG=(NTRWG-1)/NCRAY                                                
      NTGW=(NTRGW-1)/NCRAY                                                
      NTGD=(NTRGD-1)/NCRAY                                                
      NRSTWG=NTRWG-NCRAY*NTWG                                             
      NRSTGW=NTRGW-NCRAY*NTGW                                             
      NRSTGD=NTRGD-NCRAY*NTGD                                             
C                                                                         
C     Calculate auxiliary values required by FFT991                       
C                                                                         
      CALL SET99(TRIG,IFAX,MG)                                            
C                                                                         
C     Set output control variables and initialise WRSPS                   
C                                                                         
      INSPC=0                                                             
      DO 28 MP=1,NCOEFF,MOCT                                              
         DO 27 JP=MP,NCOEFF,MH                                            
            INSPC=INSPC+1                                                 
   27    CONTINUE                                                         
   28 CONTINUE                                                            
      CALL WRSPS(Z(1),1)                                                  
      IF (NLAT.NE.0) THEN                                                 
         INLAT=JGG/NLAT                                                   
      ELSE                                                                
         INLAT=0                                                          
      ENDIF                                                               
C                                                                         
C     Set up CMPA array to calculate x-derivative of half transforms      
C                                                                         
      DO 40 I=1,IGL                                                       
         CMPA(I)=0.0                                                      
   40 CONTINUE                                                            
      NROW=0                                                              
      DO 41 MP=1,MFP,MOCT                                                 
         NROW=NROW+1                                                      
         CMPA(NROW)=CMPLX(0.,REAL(MP-1))                                  
   41 CONTINUE                                                            
      IF(NHEM.EQ.2)THEN                                                   
         NROW=0                                                           
CDIR$    IVDEP                                                            
         DO 42 MP=1,MFP,MOCT                                              
            NROW=NROW+1                                                   
            CMPA(NROW+IDL)=CMPA(NROW)                                     
   42    CONTINUE                                                         
      END IF
      ! Kennedy changing where cloud stuff gets read in from fort.7
      ! These get overwritten immediately, the following are just placeholders
      AEROSOLMODEL = 'Global'
      AEROSOLCOMP  = 'standard'
      MTLX         = .1
      METALLICITY  = 1.0
      HAZES = .False.
      PICKET_FENCE_CLOUDS = .False.
      MOLEF        = (/ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 /)
      AERLAYERS    = 20.0
      AERTOTTAU    = 100.
      CLOUDBASE    = 10.
      CLOUDTOP     = 0.0001
      AERHFRAC     = 1.
      PI0AERSW     = 1.0
      ASYMSW       = 0.376
      EXTFACTLW    = 0.01
      PI0AERLW     = 0.0452
      ASYMLW       = 0.005896
      DELTASCALE   = .True.
      SIG_AREA     = 25.
      PHI_LON      = 65.
      GRAYCLDV     = .False.
      READ(7,INCLOUDY)
      CALL get_cloud_scattering_properties_wrapper


      END                                                                 
