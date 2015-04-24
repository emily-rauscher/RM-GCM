C**********************************************************               
C             SUBROUTINE SW                                               
C**********************************************************               
      SUBROUTINE SW                                                       
     :  (IAER,PAER,ZEPSCQ,ZEPSCT,ZEELOG,                                  
     :   PSCT,PRMU0,PCARDI,PPSOL,PWV,                                     
     :   PCLDSW,POZ,PPMB,PTAVE,PFDOWN,PFUP,                               
     :   PCG,POMEGA,PTAU,PALBS)                                           
                                                                          
C                                                                         
C**** *SW* - COMPUTES THE SHORTWAVE RADIATION FLUXES.                     
C                                                                         
C     PURPOSE.                                                            
C     --------                                                            
C           COMPUTES THE SHORTWAVE RADIATION FLUXES IN TWO SPECTRAL       
C     INTERVALS FOLLOWING FOUQUART AND BONNEL (1980).                     
C                                                                         
                                                                          
C                                                                         
C          1. COMPUTES ABSORBER AMOUNTS WITH TEMPERATURE AND PRESSURE     
C     SCALING.                                                            
C          2. COMPUTES UPWARD AND DOWNWARD FLUXES IN THE 0.25-0.68        
C     MICRON SPECTRAL INTERVAL.                                           
C          3. COMPUTES UPWARD AND DOWNWARD FLUXES IN THE 0.68-4.0         
C     MICRON SPECTRAL INTERVAL.                                           
C                                                                         
C PMF 23-2-98                                                             
C     ----------------------------------------------------------------    
C                                                                         
C                                                                         
C                                                                         
C     Determines model resolution                                         
C                                                                         
      include 'params.i'
C      PARAMETER(NN=21,MM=21,NHEM=2,NL=5,MOCT=1,MG=64,JG=16,NWJ2=121       
C     P         ,NCRAY=64,JGL=JG,NTRAC=1,NLEVRF=1)                         
                                                                          
C                                                                         
      INTEGER NLON,NLEV,NLM,KDLON,JLON,JK,IAER,KFLEV                      
     : ,JKL,JAE,NLP,JKP1,JKLP1                                            
      PARAMETER(NLON=1,NLEV=NL+1)                                         
      PARAMETER (NLM=NLEV-1,KDLON=NLON                                    
     :          ,KFLEV=NLEV,NLP=KFLEV+1)                                  
c passed stuff                                                            
      REAL ZEPSCQ,ZEPSCT,ZEELOG                                           
      REAL PALBS(NLON,2),PFUP(NLON,NLP),PFDOWN(NLON,NLP)                  
     :,     PPMB(NLON,NLP),PRMU0(NLON)                                    
      REAL PTAVE(NLON,NLEV),POZ(NLON,NLEV),PCLDSW(NLON,NLEV)              
     :,    PTAU(NLON,2,NLEV),PCG(NLON,2,NLEV),PWV(NLON,NLEV)              
     :,    POMEGA(NLON,2,NLEV),PAER(NLON,NLEV,5)                          
      REAL PPSOL(NLON)                                                    
      REAL PSCT,PCARDI                                                    
C*********************                                                    
      REAL   ZAKI  (NLON , 2 )        ,                                   
     S       ZCGAZ (NLON , NLEV)     ,                                    
     S       ZC1I  (NLON , NLP)    ,                                      
     S       ZDSIG (NLON , NLEV)     ,                                    
     S       ZFACT (NLON)             ,                                   
     S       ZFD   (NLON , NLP)    ,                                      
     S       ZFU   (NLON , NLP)    ,                                      
     S       ZG    (NLON)             ,                                   
     S       ZGG   (NLON)             ,                                   
     S       ZPIZAZ(NLON , NLEV)     ,                                    
     S       ZRAYL (NLON),  ZRAY1 (NLON , NLP)                            
      REAL                                                                
     S       ZRAY2 (NLON , NLP)    ,                                      
     S       ZREF  (NLON)             ,                                   
     S       ZREFZ (NLON , 2 , NLP),                                      
     S       ZRE1  (NLON)             ,                                   
     S       ZRE2  (NLON)             ,                                   
     S       ZRJ   (NLON , 6 , NLP),                                      
     S       ZRK   (NLON , 6 , NLP),                                      
     S       ZRL   (NLON , 8)         ,                                   
     S       ZRMU  (NLON)             ,                                   
     S       ZRMUE (NLON , NLP)    ,                                      
     S       ZRMUZ (NLON)                                                 
      REAL                                                                
     S       ZRNEB (NLON)             ,                                   
     S       ZRUEF (NLON , 8)         ,                                   
     S       ZR1   (NLON)             ,                                   
     S       ZR21  (NLON)             ,                                   
     S       ZR22  (NLON)             ,                                   
     S       ZR23  (NLON)             ,                                   
     S       ZS    (NLON)             ,                                   
     S       ZSEC  (NLON)             ,                                   
     S       ZSS1  (NLON)             ,                                   
     S       ZTAUAZ(NLON , NLEV)     ,                                    
     S       ZTO1  (NLON)             ,                                   
     S       ZTR   (NLON , 2 , NLP),                                      
     S       ZTRA1 (NLON , NLP)    ,                                      
     S       ZTRA2 (NLON , NLP)                                           
      REAL                                                                
     S       ZTR1  (NLON)             ,                                   
     S       ZTR2  (NLON)             ,                                   
     S       ZUD   (NLON , 3 , NLP),                                      
     S       ZUM   (NLON , NLP)    ,                                      
     S       ZU1D  (NLON)             ,                                   
     S       ZU2D  (NLON)             ,                                   
     S       ZW    (NLON)             ,                                   
     S       ZN175 (NLON)             ,                                   
     S       ZN190 (NLON)             ,                                   
     S       ZO175 (NLON)             ,                                   
     S       ZO190 (NLON)             ,                                   
     S       ZP75  (NLON)             ,                                   
     S       ZP90  (NLON)             ,                                   
     S       ZSIGN (NLON)             ,                                   
     S       ZSIGO (NLON)                                                 
      INTEGER in,ja,k,j,i,inu,jkm1,jaj,jkki,jn2j,jkkp4,kref               
      INTEGER jajp,jabs,jn                                                
      REAL xl2o,zpsig,ztray,zgar,zratio,zff,zfaoa,zfaoc                   
      REAL zcorae,zcorcd,zmue,zgap,zbmu0,zww,zto,zden                     
      REAL zmu1,zbmu1,zden1,zre11,zrmum1,zaa,zrki                         
      REAL XL2, XL3, ALOS,SOLC,ZRT                                        
      REAL ZWH2O,ZDSCO2,ZDSH2O                                            
      REAL FOXDWN(NLON, NLP)                                              
      REAL CH2O,CCO2                                                      
C SUSW and SUAER vairables                                                
      REAL SUN(2),D(2,3),APAD(2,3,7),BPAD(2,3,7),CRAY(2,6)                
      REAL TAUA(2,5),PIZA(2,5),CGA(2,5)                                   
C                                                                         
C                                                                         
      CH2O   = 5.3669274E-03                                              
      CCO2   = 5.8269497E-03                                              
C SW aerosols init PMF 23-2-98                                            
C                                                                         
C                                                                         
C      -------------------------------------------------------------      
C                                                                         
C*       1.    SHORTWAVE COEFFICIENTS                                     
C              ----------------------                                     
C                                                                         
      DATA ((TAUA(IN,JA),JA=1,5),IN=1,2) /                                
     S .730719, .912819, .725059, .745405, .682188 ,                      
     S .730719, .912819, .725059, .745405, .682188 /                      
      DATA ((PIZA(IN,JA),JA=1,5),IN=1,2) /                                
     S .872212, .982545, .623143, .944887, .997975 ,                      
     S .872212, .982545, .623143, .944887, .997975 /                      
      DATA ((CGA (IN,JA),JA=1,5),IN=1,2) /                                
     S .647596, .739002, .580845, .662657, .624246 ,                      
     S .647596, .739002, .580845, .662657, .624246 /                      
C      -------------------------------------------------------------      
C Init SW coeffs PMF 23-2-98                                              
C                                                                         
C      ----------------------------------------------------------------   
C                                                                         
C*       1.    SET VALUES.                                                
C              -----------                                                
C    COEFFICIENTS FOR THE SHORTWAVE RADIATION SUBROUTINES                 
C                                                                         
C APAD  :  PADE APPROXIMANTS NUMERATOR                                    
C BPAD  :  PADE APPROXIMANTS DENOMINATOR                                  
C D     :  TRANSMISSION LIMIT FOR INFINITE ABSORBER AMOUNT                
C CRAY  :  RAYLEIGH SCATTERING COEFFICIENTS                               
C SUN   :  SOLAR FRACTION IN SPECTRAL INTERVALS                           
C                                                                         
      DATA SOLC/1376.0/                                                   
      DATA SUN(1) / 0.441130 /                                            
      DATA (D(1,K),K = 1,3) / 0.00, 0.00, 0.00 /                          
c old without O2                                                          
c      DATA ((APAD(1,I,J),I=1,3),J=1,7) /                                 
c     S 0.000000000E-00, 0.000000000E-00, 0.925887084E-04,                
c     S 0.000000000E-00, 0.000000000E-00, 0.129353723E-01,                
c     S 0.000000000E-00, 0.000000000E-00, 0.800821928E+00,                
c     S 0.000000000E-00, 0.000000000E-00, 0.242715973E+02,                
c     S 0.000000000E-00, 0.000000000E-00, 0.878331486E+02,                
c     S 0.000000000E-00, 0.000000000E-00, 0.191559725E+02,                
c     S 0.000000000E-00, 0.000000000E-00, 0.000000000E+00 /               
C                                                                         
c      DATA ((BPAD(1,I,J),I=1,3),J=1,7) /                                 
c     S 0.000000000E-00, 0.000000000E-00, 0.925887084E-04,                
c     S 0.000000000E-00, 0.000000000E-00, 0.131812683E-01,                
c     S 0.000000000E-00, 0.000000000E-00, 0.812706117E+00,                
c     S 0.000000000E-00, 0.000000000E-00, 0.249863591E+02,                
c     S 0.000000000E-00, 0.000000000E-00, 0.931071925E+02,                
c     S 0.000000000E-00, 0.000000000E-00, 0.252233437E+02,                
c     S 0.000000000E-00, 0.000000000E-00, 0.100000000E+01 /               
C new with O2                                                             
C                                                                         
      DATA ((APAD(1,I,J),I=1,3),J=1,7) /                                  
     S 0.000000000E-00, 0.000000000E-00,  0.200763774E+00,                
     S 0.000000000E-00, 0.000000000E-00,  0.273755555E+02,                
     S 0.000000000E-00, 0.000000000E-00,  0.759011354E+02,                
     S 0.000000000E-00, 0.000000000E-00, -0.542785355E+03,                
     S 0.000000000E-00, 0.000000000E-00,  0.771898951E+03,                
     S 0.000000000E-00, 0.000000000E-00,  0.240490730E+02,                
     S 0.000000000E-00, 0.000000000E-00,  0.000000000E+00 /               
C                                                                         
      DATA ((BPAD(1,I,J),I=1,3),J=1,7) /                                  
     S 0.000000000E-00, 0.000000000E-00,  0.200763774E+00,                
     S 0.000000000E-00, 0.000000000E-00,  0.278643596E+02,                
     S 0.000000000E-00, 0.000000000E-00,  0.827455632E+02,                
     S 0.000000000E-00, 0.000000000E-00, -0.561842163E+03,                
     S 0.000000000E-00, 0.000000000E-00,  0.768518240E+03,                
     S 0.000000000E-00, 0.000000000E-00,  0.727341645E+02,                
     S 0.000000000E-00, 0.000000000E-00,  0.100000000E+01 /               
C                                                                         
      DATA (CRAY(1,K),K=1,6) /                                            
     S .428937E-01, .890743E+00,-.288555E+01,                             
     S .522744E+01,-.469173E+01, .161645E+01/                             
C                                                                         
c      DATA SUN(2) / 0.558324 /                                           
      DATA SUN(2) / 0.557633 /                                            
      DATA (D(2,K),K=1,3) / 0.317735694, 0.775570952, 0.800000000 /       
c      DATA ((APAD(2,I,J),I=1,3),J=1,7) /                                 
c     S 0.822745535E-02, 0.145321703E-04, 0.410177786E+03,                
c     S 0.705825794E+01, 0.175741897E-01, 0.672595424E+02,                
c     S 0.348747605E+03, 0.259696276E+01, 0.000000000E-00,                
c     S 0.174921268E+04, 0.599852834E+02, 0.000000000E-00,                
c     S 0.100138721E+04, 0.203510317E+03, 0.000000000E-00,                
c     S 0.518496206E+02, 0.757222990E+02, 0.000000000E-00,                
c     S 0.000000000E+00, 0.000000000E+00, 0.000000000E+00 /               
C                                                                         
c      DATA ((BPAD(2,I,J),I=1,3),J=1,7) /                                 
c     S 0.822745535E-02, 0.145321703E-04, 0.410177786E+03,                
c     S 0.719686994E+01, 0.176779370E-01, 0.731185438E+02,                
c     S 0.381961022E+03, 0.265802733E+01, 0.100000000E+01,                
c     S 0.219460901E+04, 0.634292876E+02, 0.000000000E+00,                
c     S 0.159321079E+04, 0.228829763E+03, 0.000000000E+00,                
c     S 0.138748279E+03, 0.992586506E+02, 0.000000000E+00,                
c     S 0.100000000E+01, 0.100000000E+01, 0.000000000E+00 /               
C                                                                         
      DATA ((APAD(2,I,J),I=1,3),J=1,7) /                                  
     S 0.822745535E-02, 0.145321703E-04, 0.262335915E+03,                 
     S 0.705825794E+01, 0.175741897E-01, 0.368136730E+03,                 
     S 0.348747605E+03, 0.259696276E+01, 0.000000000E-00,                 
     S 0.174921268E+04, 0.599852834E+02, 0.000000000E-00,                 
     S 0.100138721E+04, 0.203510317E+03, 0.000000000E-00,                 
     S 0.518496206E+02, 0.757222990E+02, 0.000000000E-00,                 
     S 0.000000000E+00, 0.000000000E+00, 0.000000000E+00 /                
C                                                                         
      DATA ((BPAD(2,I,J),I=1,3),J=1,7) /                                  
     S 0.822745535E-02, 0.145321703E-04, 0.262335915E+03,                 
     S 0.719686994E+01, 0.176779370E-01, 0.368956332E+03,                 
     S 0.381961022E+03, 0.265802733E+01, 0.100000000E+01,                 
     S 0.219460901E+04, 0.634292876E+02, 0.000000000E+00,                 
     S 0.159321079E+04, 0.228829763E+03, 0.000000000E+00,                 
     S 0.138748279E+03, 0.992586506E+02, 0.000000000E+00,                 
     S 0.100000000E+01, 0.100000000E+01, 0.000000000E+00 /                
                                                                          
      DATA (CRAY(2,K),K=1,6) /                                            
     S .697200E-02, .173297E-01,-.850903E-01,                             
     S .248261E+00,-.302031E+00, .129662E+00/                             
C                                                                         
C*         1.     COMPUTES AMOUNTS OF ABSORBERS                           
C                 -----------------------------                           
C                                                                         
C*         1.1    INITIALIZES QUANTITIES                                  
C                 ----------------------                                  
C                                                                         
      DO 111 JLON = 1 , KDLON                                             
        ZC1I(JLON,KFLEV+1)=0.                                             
        ZUD(JLON,1,KFLEV+1)=0.                                            
        ZUD(JLON,2,KFLEV+1)=0.                                            
        ZUD(JLON,3,KFLEV+1)=0.                                            
C        ZFACT(JLON)= PRMU0(JLON) * PSCT * RDAYL(JLON)                    
        ZFACT(JLON)= PRMU0(JLON) * PSCT                                   
        ZRMU(JLON)=SQRT(1224.* PRMU0(JLON) * PRMU0(JLON) + 1.) / 35.      
        ZSEC(JLON)=1./ZRMU(JLON)                                          
 111  CONTINUE                                                            
C                                                                         
C                                                                         
C*         1.2    OZONE FOR DOWNWARD LOOKING PATH                         
C                 -------------------------------                         
C                                                                         
 120  CONTINUE                                                            
C                                                                         
      DO 122 JK = 1 , KFLEV                                               
       JKL = KFLEV+1 - JK                                                 
       JKLP1 = JKL + 1                                                    
       DO 121 JLON = 1 , KDLON                                            
        ZUD(JLON,3,JKL)=ZUD(JLON,3,JKLP1)+POZ(JLON,JKL)*ZSEC(JLON)        
 121   CONTINUE                                                           
 122  CONTINUE                                                            
C                                                                         
c          SCHUMANN-RUNGE and HERZBERG CONTINUUM (0.12 - 0.25 MICRON)     
c          ----------------------------------------------------------     
c                                                                         
      ALOS=2.687E19                                                       
      DO 50 JLON = 1 , KDLON                                              
        FOXDWN(JLON,1)=0.0                                                
        IF (PRMU0(JLON).GT.0.0001) THEN                                   
          XL3=0.0                                                         
          DO 55 JK=1,KFLEV                                                
            XL3=XL3+POZ(JLON,JK)*ZSEC(JLON)*ALOS                          
   55     CONTINUE                                                        
C                                                                         
C set XL2O as mesopause value if too small                                
          XL2O=4.442E19*1.E-4*ZSEC(JLON)                                  
C                                                                         
          DO 60 JK=1,KFLEV                                                
            XL3=XL3-POZ(JLON,JK)*ZSEC(JLON)*ALOS                          
C            XL2=4.442E19*PPMB(JLON,JK+1)*ZSEC(JLON)*100.0                
            XL2=MAX(4.442E19*PPMB(JLON,JK+1)*ZSEC(JLON)*100.0, XL2O)      
C                                                                         
            FOXDWN(JLON,JK+1) =                                           
     *      (  0.913 *EXP(-5.500E-24*XL2-6.215E-18*XL3)                   
     *    + 0.6308 *EXP(-1.342E-26*XL2-1.656E-18*XL3)                     
     *    + 0.4*2.30E-3*EXP(-3.159E-20*XL2-1.857E-06*XL2**0.33837)        
     *    + 0.4*3.00E-3*EXP(-2.261E-24*XL2-1.917E-14*XL2**0.78903)        
     *    + 0.3*3.50E-3*EXP(-5.399E-24*XL2-2.466E-17*XL2**0.92105)        
     *    + 0.3*2.80E-3*EXP(-3.406E-24*XL2-7.787E-18*XL2**0.9269 )        
     *    + 0.1*0.1466 *EXP(-2.915E-25*XL2**1.06135                       
     *                   - 1.866E-07*XL2**0.2895) )/SOLC                  
C                                                                         
   60     CONTINUE                                                        
        ELSE                                                              
          DO 65 JK=KFLEV,1,-1                                             
            FOXDWN(JLON,JK)=0.0                                           
   65     CONTINUE                                                        
        ENDIF                                                             
   50 CONTINUE                                                            
C                                                                         
C*         1.3    OZONE FOR UPWARD LOOKING PATH AND OTHER ABSORBERS       
C                 -------------------------------------------------       
C                                                                         
 130  CONTINUE                                                            
C                                                                         
      DO 131 JLON = 1 , KDLON                                             
         ZUM(JLON,1) = ZUD(JLON,3,1)                                      
         ZU1D(JLON) = 0.                                                  
         ZU2D(JLON) = 0.                                                  
         ZPSIG = PPSOL(JLON) / 101325.                                    
         ZP75(JLON) = PPSOL(JLON) * ZPSIG ** 0.75                         
         ZP90(JLON) = PPSOL(JLON) * ZPSIG ** 0.90                         
         ZO175(JLON) = 1.0                                                
         ZO190(JLON) = 1.0                                                
         ZSIGO(JLON) = 1.0                                                
 131  CONTINUE                                                            
C                                                                         
      DO 133 JK = 1 , KFLEV                                               
         JKP1 = JK + 1                                                    
         JKL = KFLEV+1 - JK                                               
         DO 132 JLON = 1 , KDLON                                          
            ZUM(JLON,JKP1) = ZUM(JLON,JK) + POZ(JLON,JK) * 1.66           
            ZRT = 273.15 / PTAVE(JLON,JK)                                 
            ZWH2O = MAX(PWV(JLON,JK) , ZEPSCQ )                           
            ZSIGN(JLON) = 100. * PPMB(JLON,JKP1) / PPSOL(JLON)            
            ZDSIG(JLON,JK) = ZSIGO(JLON) - ZSIGN(JLON)                    
            ZN175(JLON) = ZSIGN(JLON) ** 1.75                             
            ZN190(JLON) = ZSIGN(JLON) ** 1.90                             
            ZDSCO2 = ZO175(JLON) - ZN175(JLON)                            
            ZDSH2O = ZO190(JLON) - ZN190(JLON)                            
            ZUD(JLON,1,JK) = ZP90(JLON) * ZDSH2O*CH2O*ZWH2O  * ZRT **0.45 
            ZUD(JLON,2,JK) = ZP75(JLON) * ZDSCO2*CCO2*PCARDI * ZRT **0.375
            ZU1D(JLON) = ZU1D(JLON) + ZUD(JLON,1,JK)                      
            ZU2D(JLON) = ZU2D(JLON) + ZUD(JLON,2,JK)                      
            ZSIGO(JLON) = ZSIGN(JLON)                                     
            ZO175(JLON) = ZN175(JLON)                                     
            ZO190(JLON) = ZN190(JLON)                                     
 132     CONTINUE                                                         
 133  CONTINUE                                                            
C                                                                         
C*         1.4    COMPUTES CLEAR-SKY GREY ABSORPTION COEFFICIENTS         
C                 -----------------------------------------------         
C                                                                         
 140  CONTINUE                                                            
C                                                                         
      DO 141 JLON = 1 , KDLON                                             
         ZU1D(JLON) = ZU1D(JLON) * ZSEC(JLON)                             
         ZU2D(JLON) = ZU2D(JLON) * ZSEC(JLON)                             
         ZW(JLON) = ZU1D(JLON)                                            
 141  CONTINUE                                                            
C                                                                         
      CALL SWTT ( 2, 1, APAD,BPAD,D,ZW,ZR1)                               
C                                                                         
      DO 142 JLON = 1 , KDLON                                             
         ZAKI(JLON,1) = -LOG( ZR1 (JLON)) / ZU1D(JLON)                    
         ZW(JLON) = ZU2D(JLON)                                            
 142  CONTINUE                                                            
C                                                                         
      CALL SWTT ( 2, 2, APAD,BPAD,D,ZW,ZR1)                               
C                                                                         
      DO 143 JLON = 1 , KDLON                                             
         ZAKI(JLON,2) = -LOG( ZR1 (JLON)) / ZU2D(JLON)                    
 143  CONTINUE                                                            
C                                                                         
C     ---------------------------------------------------------------     
C                                                                         
C*         2.     FIRST SPECTRAL INTERVAL (0.25-0.68 MICRON)              
C                 ----------------------- ------------------              
C                                                                         
 200  CONTINUE                                                            
C                                                                         
      INU = 1                                                             
C                                                                         
C*         2.1    OPTICAL THICKNESS FOR RAYLEIGH SCATTERING               
C                 -----------------------------------------               
C                                                                         
 210  CONTINUE                                                            
C                                                                         
      DO 211 JLON = 1 , KDLON                                             
          ZRAYL(JLON) = CRAY(INU,1) + ZRMU(JLON)                          
     :         * (CRAY(INU,2) + ZRMU(JLON)                                
     S         * (CRAY(INU,3) + ZRMU(JLON) * (CRAY(INU,4)+ZRMU(JLON)      
     S         * (CRAY(INU,5) + ZRMU(JLON) *   CRAY(INU,6)))))            
 211  CONTINUE                                                            
C                                                                         
C*         2.2    OPTICAL PARAMETERS FOR AEROSOLS AND RAYLEIGH            
C                 --------------------------------------------            
C                                                                         
 220  CONTINUE                                                            
C                                                                         
      DO 225 JK = 1 , KFLEV                                               
         DO 221 JLON = 1 , KDLON                                          
            ZCGAZ(JLON,JK) = 0.                                           
            ZPIZAZ(JLON,JK) =  0.                                         
            ZTAUAZ(JLON,JK) = 0.                                          
 221     CONTINUE                                                         
         DO 223 JAE=1,5                                                   
            DO 222 JLON = 1 , KDLON                                       
               ZTAUAZ(JLON,JK)=ZTAUAZ(JLON,JK)                            
     S                      +PAER(JLON,JK,JAE)*TAUA(INU,JAE)              
               ZPIZAZ(JLON,JK)=ZPIZAZ(JLON,JK)+PAER(JLON,JK,JAE)          
     S                   * TAUA(INU,JAE)*PIZA(INU,JAE)                    
               ZCGAZ(JLON,JK) =  ZCGAZ(JLON,JK) +PAER(JLON,JK,JAE)        
     S                   * TAUA(INU,JAE)*PIZA(INU,JAE)*CGA(INU,JAE)       
 222        CONTINUE                                                      
 223     CONTINUE                                                         
C                                                                         
         DO 224 JLON = 1 , KDLON                                          
            IF(IAER.EQ.0) THEN                                            
              ZCGAZ(JLON,JK) = 0.0                                        
              ZPIZAZ(JLON,JK) = 1.0                                       
            ELSE                                                          
              ZCGAZ(JLON,JK) = ZCGAZ(JLON,JK) / ZPIZAZ(JLON,JK)           
              ZPIZAZ(JLON,JK) = ZPIZAZ(JLON,JK) / ZTAUAZ(JLON,JK)         
            ENDIF                                                         
            ZTRAY = ZRAYL(JLON) * ZDSIG(JLON,JK)                          
            ZRATIO = ZTRAY / (ZTRAY + ZTAUAZ(JLON,JK))                    
            ZGAR = ZCGAZ(JLON,JK)                                         
            ZFF = ZGAR * ZGAR                                             
            ZTAUAZ(JLON,JK)=ZTRAY+ZTAUAZ(JLON,JK)                         
     :                      *(1.-ZPIZAZ(JLON,JK)*ZFF)                     
            ZCGAZ(JLON,JK) = ZGAR * (1. - ZRATIO) / (1. + ZGAR)           
            ZPIZAZ(JLON,JK) =ZRATIO+(1.-ZRATIO)*ZPIZAZ(JLON,JK)           
     S                   *(1.-ZFF) / (1. - ZPIZAZ(JLON,JK) * ZFF)         
 224     CONTINUE                                                         
 225  CONTINUE                                                            
                                                                          
C                                                                         
C*         2.3    TOTAL EFFECTIVE CLOUDINESS ABOVE A GIVEN LEVEL          
C                 ----------------------------------------------          
C                                                                         
 230  CONTINUE                                                            
C                                                                         
      DO 231 JLON = 1 , KDLON                                             
         ZR23(JLON) = 0.                                                  
         ZC1I(JLON,KFLEV+1) = 0.                                          
 231  CONTINUE                                                            
C                                                                         
      DO 233 JK = 1 , KFLEV                                               
         JKL = KFLEV+1 - JK                                               
         JKLP1 = JKL + 1                                                  
         DO 232 JLON = 1 , KDLON                                          
           ZFAOA = 1.-ZPIZAZ(JLON,JKL)*ZCGAZ(JLON,JKL)*ZCGAZ(JLON,JKL)    
           ZFAOC = 1. - POMEGA(JLON,INU,JKL) * PCG(JLON,INU,JKL)          
     S                                       * PCG(JLON,INU,JKL)          
           ZCORAE = ZFAOA * ZTAUAZ(JLON,JKL) * ZSEC(JLON)                 
           ZCORCD = ZFAOC * PTAU(JLON,INU,JKL) * ZSEC(JLON)               
           ZR21(JLON) = EXP(-ZCORAE   )                                   
           ZR22(JLON) = EXP(-ZCORCD   )                                   
           ZSS1(JLON) = PCLDSW(JLON,JKL)*(1.0-ZR21(JLON)*ZR22(JLON))      
     S               + (1.0-PCLDSW(JLON,JKL))*(1.0-ZR21(JLON))            
           ZC1I(JLON,JKL) = 1.0-(1.0-ZSS1(JLON))*(1.0-ZC1I(JLON,JKLP1))   
 232     CONTINUE                                                         
 233  CONTINUE                                                            
                                                                          
C                                                                         
C*         2.4    REFLECTIVITY/TRANSMISSIVITY FOR PURE SCATTERING         
C                 -----------------------------------------------         
C                                                                         
 240  CONTINUE                                                            
C                                                                         
      DO 241 JLON = 1 , KDLON                                             
         ZREFZ(JLON,2,1) = PALBS(JLON,INU)                                
         ZREFZ(JLON,1,1) = PALBS(JLON,INU)                                
 241  CONTINUE                                                            
C                                                                         
      DO 246 JK = 2 , KFLEV+1                                             
         JKM1 = JK-1                                                      
         DO 242 JLON = 1 , KDLON                                          
            ZRNEB(JLON)= PCLDSW(JLON,JKM1)                                
             ZRE1(JLON)=0.                                                
            ZTR1(JLON)=0.                                                 
            ZRE2(JLON)=0.                                                 
            ZTR2(JLON)=0.                                                 
C                                                                         
C*         2.4.1  EQUIVALENT ZENITH ANGLE                                 
C                 -----------------------                                 
C                                                                         
                                                                          
 2410 CONTINUE                                                            
C                                                                         
            ZMUE = (1.-ZC1I(JLON,JK)) * ZSEC(JLON)                        
     S            + ZC1I(JLON,JK) * 1.66                                  
            ZRMUE(JLON,JK) = 1./ZMUE                                      
C                                                                         
C*         2.4.2  REFLECT./TRANSMISSIVITY DUE TO RAYLEIGH AND AEROSOLS    
C                 ----------------------------------------------------    
C                                                                         
                                                                          
 2420 CONTINUE                                                            
C                                                                         
            ZGAP = ZCGAZ(JLON,JKM1)                                       
            ZBMU0 = 0.5 - 0.75 * ZGAP / ZMUE                              
            ZWW = ZPIZAZ(JLON,JKM1)                                       
            ZTO = ZTAUAZ(JLON,JKM1)                                       
            ZDEN = 1. + (1. - ZWW + ZBMU0 * ZWW) * ZTO * ZMUE             
     S           + (1-ZWW) * (1. - ZWW +2.*ZBMU0*ZWW)*ZTO*ZTO*ZMUE*ZMUE   
            ZRAY1(JLON,JKM1) = ZBMU0 * ZWW * ZTO * ZMUE / ZDEN            
            ZTRA1(JLON,JKM1) = 1. / ZDEN                                  
C                                                                         
            ZMU1 = 0.5                                                    
            ZBMU1 = 0.5 - 0.75 * ZGAP * ZMU1                              
            ZDEN1= 1. + (1. - ZWW + ZBMU1 * ZWW) * ZTO / ZMU1             
     S           + (1-ZWW) * (1. - ZWW +2.*ZBMU1*ZWW)*ZTO*ZTO/ZMU1/ZMU1   
            ZRAY2(JLON,JKM1) = ZBMU1 * ZWW * ZTO / ZMU1 / ZDEN1           
            ZTRA2(JLON,JKM1) = 1. / ZDEN1                                 
C                                                                         
C*         2.4.3  EFFECT OF CLOUD LAYER                                   
C                 ---------------------                                   
C                                                                         
 2430 CONTINUE                                                            
C                                                                         
            PTAU(JLON,INU,JKM1) = MAX( PTAU(JLON,INU,JKM1) , ZEPSCT )     
            ZW(JLON) = POMEGA(JLON,INU,JKM1)                              
            ZTO1(JLON) = PTAU(JLON,INU,JKM1)/ZW(JLON)                     
     S            + ZTAUAZ(JLON,JKM1)/ZPIZAZ(JLON,JKM1)                   
            ZR21(JLON) = PTAU(JLON,INU,JKM1) + ZTAUAZ(JLON,JKM1)          
            ZR22(JLON) = PTAU(JLON,INU,JKM1) / ZR21(JLON)                 
            ZGG(JLON) = ZR22(JLON) * PCG(JLON,INU,JKM1)                   
     S              + (1. - ZR22(JLON)) * ZCGAZ(JLON,JKM1)                
            ZW(JLON) = ZR21(JLON) / ZTO1(JLON)                            
            ZREF(JLON) = ZREFZ(JLON,1,JKM1)                               
            ZRMUZ(JLON) = ZRMUE(JLON,JK)                                  
 242     CONTINUE                                                         
C                                                                         
         CALL DEDD(ZGG,ZREF,ZRMUZ,ZTO1,ZW,ZRE1,ZTR1,ZRE2,ZTR2)            
C                                                                         
                                                                          
         DO 245 JLON = 1 , KDLON                                          
C                                                                         
             ZREFZ(JLON,1,JK) = (1.-ZRNEB(JLON)) * (ZRAY1(JLON,JKM1)      
     S                     + ZREFZ(JLON,1,JKM1) * ZTRA1(JLON,JKM1)        
     S                     * ZTRA2(JLON,JKM1)                             
     S                     /(1.-ZRAY2(JLON,JKM1)*ZREFZ(JLON,1,JKM1)))     
     S                     + ZRNEB(JLON) * ZRE2(JLON)                     
C                                                                         
            ZTR(JLON,1,JKM1) = ZRNEB(JLON) * ZTR2(JLON) +                 
     S                     (ZTRA1(JLON,JKM1)/(1.-ZRAY2(JLON,JKM1)         
     S                     *ZREFZ(JLON,1,JKM1)))* (1.-ZRNEB(JLON))        
C                                                                         
 245     CONTINUE                                                         
 246  CONTINUE                                                            
C                                                                         
C*         2.5    REFLECT./TRANSMISSIVITY BETWEEN SURFACE AND LEVEL       
C                 -------------------------------------------------       
C                                                                         
                                                                          
 250  CONTINUE                                                            
C                                                                         
      JAJ = 2                                                             
      DO 251 JLON = 1 , KDLON                                             
         ZRJ(JLON,JAJ,KFLEV+1) = 1.                                       
         ZRK(JLON,JAJ,KFLEV+1) = ZREFZ(JLON, 1,KFLEV+1)                   
 251  CONTINUE                                                            
C                                                                         
      DO 253 JK = 1 , KFLEV                                               
         JKL = KFLEV+1 - JK                                               
         JKLP1 = JKL + 1                                                  
         DO 252 JLON = 1 , KDLON                                          
            ZRE11= ZRJ(JLON,JAJ,JKLP1) * ZTR(JLON, 1,JKL)                 
            ZRJ(JLON,JAJ,JKL) = ZRE11                                     
            ZRK(JLON,JAJ,JKL) = ZRE11 * ZREFZ(JLON, 1,JKL)                
 252     CONTINUE                                                         
 253  CONTINUE                                                            
C                                                                         
C*         2.6    OZONE ABSORPTION AND FLUXES                             
C                 ---------------------------                             
C                                                                         
 260  CONTINUE                                                            
C                                                                         
      DO 264 JK = 1 , KFLEV+1                                             
         JKL = KFLEV+1 - JK + 1                                           
         DO 262 JLON = 1 , KDLON                                          
            ZW(JLON) = ZUD(JLON,3,JKL)                                    
 262     CONTINUE                                                         
C                                                                         
         CALL SWTT ( INU, 3, APAD,BPAD,D,ZW,ZR1)                          
C                                                                         
         DO 263 JLON = 1 , KDLON                                          
           ZFD(JLON,JKL) = ZR1(JLON)*ZRJ(JLON,JAJ,JKL)*SUN(INU)           
 263     CONTINUE                                                         
 264  CONTINUE                                                            
C                                                                         
      DO 265 JLON = 1 , KDLON                                             
         ZFU(JLON,1) = PALBS(JLON,INU) * ZFD(JLON,1)                      
 265  CONTINUE                                                            
C                                                                         
      DO 268 JK = 1 , KFLEV+1                                             
         DO 266 JLON = 1 , KDLON                                          
            ZW(JLON) = ZUM(JLON,JK)                                       
 266     CONTINUE                                                         
C                                                                         
         CALL SWTT ( INU, 3, APAD,BPAD,D,ZW,ZR1)                          
C                                                                         
         DO 267 JLON = 1 , KDLON                                          
           ZFU(JLON,JK) = ZR1(JLON) * ZRK(JLON,JAJ,JK) * SUN(INU)         
 267     CONTINUE                                                         
 268  CONTINUE                                                            
C                                                                         
C                                                                         
 878  FORMAT(1X,' LOOP ',I4,' SUCCESSFUL IN SW    INDEX=',I3)             
C     ------------------------------------------------------------------  
C                                                                         
C*         3.     SECOND SPECTRAL INTERVAL (0.68-4.00 MICRON)             
C                 ----------------------- -------------------             
C                                                                         
 300  CONTINUE                                                            
C                                                                         
      INU = 2                                                             
C                                                                         
C*         3.1    OPTICAL THICKNESS FOR RAYLEIGH SCATTERING               
C                 -----------------------------------------               
C                                                                         
 310  CONTINUE                                                            
C                                                                         
      DO 311 JLON = 1 , KDLON                                             
         ZRMUM1 = 1. - ZRMU(JLON)                                         
         ZRAYL(JLON) = CRAY(INU,1) + ZRMUM1   * (CRAY(INU,2) + ZRMUM1     
     S            * (CRAY(INU,3) + ZRMUM1   * (CRAY(INU,4) + ZRMUM1       
     S            * (CRAY(INU,5) + ZRMUM1   *  CRAY(INU,6)       ))))     
 311  CONTINUE                                                            
C                                                                         
C*         3.2    OPTICAL PARAMETERS FOR AEROSOLS AND RAYLEIGH            
C                 --------------------------------------------            
C                                                                         
 320  CONTINUE                                                            
C                                                                         
      DO 325 JK = 1 , KFLEV                                               
         DO 321 JLON = 1 , KDLON                                          
            ZCGAZ(JLON,JK) = 0.                                           
            ZPIZAZ(JLON,JK) = 0.                                          
            ZTAUAZ(JLON,JK) = 0.                                          
 321     CONTINUE                                                         
         DO 323 JAE=1,5                                                   
            DO 322 JLON = 1 , KDLON                                       
               ZTAUAZ(JLON,JK) = ZTAUAZ(JLON,JK) + PAER(JLON,JK,JAE)      
     S                       * TAUA(INU,JAE)                              
               ZPIZAZ(JLON,JK) = ZPIZAZ(JLON,JK) + PAER(JLON,JK,JAE)      
     S                       * TAUA(INU,JAE) * PIZA(INU,JAE)              
               ZCGAZ(JLON,JK) =  ZCGAZ(JLON,JK) + PAER(JLON,JK,JAE)       
     S                   * TAUA(INU,JAE)*PIZA(INU,JAE)*CGA(INU,JAE)       
 322        CONTINUE                                                      
 323     CONTINUE                                                         
         DO 324 JLON = 1 , KDLON                                          
            IF(IAER.EQ.0) THEN                                            
              ZCGAZ(JLON,JK) = 0.0                                        
              ZPIZAZ(JLON,JK) = 1.0                                       
            ELSE                                                          
              ZCGAZ(JLON,JK) = ZCGAZ(JLON,JK) / ZPIZAZ(JLON,JK)           
              ZPIZAZ(JLON,JK) = ZPIZAZ(JLON,JK) / ZTAUAZ(JLON,JK)         
            ENDIF                                                         
            ZTRAY = ZRAYL(JLON) * ZDSIG(JLON,JK)                          
            ZRATIO = ZTRAY / (ZTRAY + ZTAUAZ(JLON,JK))                    
            ZGAR = ZCGAZ(JLON,JK)                                         
            ZFF = ZGAR * ZGAR                                             
            ZTAUAZ(JLON,JK)=ZTRAY+ZTAUAZ(JLON,JK)                         
     :                           *(1.-ZPIZAZ(JLON,JK)*ZFF)                
            ZCGAZ(JLON,JK) = ZGAR * (1. - ZRATIO) / (1. + ZGAR)           
            ZPIZAZ(JLON,JK) = ZRATIO+(1. - ZRATIO)*ZPIZAZ(JLON,JK)        
     S                   *(1.-ZFF) / (1. - ZPIZAZ(JLON,JK) * ZFF)         
 324     CONTINUE                                                         
 325  CONTINUE                                                            
C                                                                         
C*         3.3    TOTAL EFFECTIVE CLOUDINESS ABOVE A GIVEN LEVEL          
C                 ----------------------------------------------          
C                                                                         
 330  CONTINUE                                                            
C                                                                         
      DO 331 JLON = 1 , KDLON                                             
         ZR23(JLON) = 0.                                                  
         ZC1I(JLON,KFLEV+1) = 0.                                          
 331  CONTINUE                                                            
      DO 333 JK = 1 , KFLEV                                               
         JKL = KFLEV+1 - JK                                               
         JKLP1 = JKL + 1                                                  
         DO 332 JLON = 1 , KDLON                                          
           ZFAOA =1.-ZPIZAZ(JLON,JKL)*ZCGAZ(JLON,JKL)*ZCGAZ(JLON,JKL)     
           ZFAOC =1. - POMEGA(JLON,INU,JKL) * PCG(JLON,INU,JKL)           
     S                                       * PCG(JLON,INU,JKL)          
           ZCORAE = ZFAOA * ZTAUAZ(JLON,JKL) * ZSEC(JLON)                 
           ZCORCD = ZFAOC * PTAU(JLON,INU,JKL) * ZSEC(JLON)               
           ZR21(JLON) = EXP(-ZCORAE   )                                   
           ZR22(JLON) = EXP(-ZCORCD   )                                   
           ZSS1(JLON) = PCLDSW(JLON,JKL)*(1.0-ZR21(JLON)*ZR22(JLON))      
     S               + (1.0-PCLDSW(JLON,JKL))*(1.0-ZR21(JLON))            
           ZC1I(JLON,JKL) = 1.0-(1.0-ZSS1(JLON))*(1.0-ZC1I(JLON,JKLP1))   
 332     CONTINUE                                                         
 333  CONTINUE                                                            
C                                                                         
C*         3.4    REFLECTIVITY/TRANSMISSIVITY FOR PURE SCATTERING         
C                 -----------------------------------------------         
C                                                                         
 340  CONTINUE                                                            
C                                                                         
      DO 341 JLON = 1 , KDLON                                             
         ZREFZ(JLON,2,1) = PALBS(JLON,INU)                                
         ZREFZ(JLON,1,1) = PALBS(JLON,INU)                                
 341  CONTINUE                                                            
C                                                                         
      DO 346 JK = 2 , KFLEV+1                                             
         JKM1 = JK - 1                                                    
         DO 342 JLON = 1 , KDLON                                          
            ZRNEB(JLON) = PCLDSW(JLON,JKM1)                               
            ZRE1(JLON)=0.                                                 
            ZTR1(JLON)=0.                                                 
            ZRE2(JLON)=0.                                                 
            ZTR2(JLON)=0.                                                 
C                                                                         
C*         3.4.1  EQUIVALENT ZENITH ANGLE                                 
C                 -----------------------                                 
C                                                                         
 3410 CONTINUE                                                            
C                                                                         
            ZMUE = (1.-ZC1I(JLON,JK)) * ZSEC(JLON)                        
     S           + ZC1I(JLON,JK) * 1.66                                   
            ZRMUE(JLON,JK) = 1./ZMUE                                      
C                                                                         
C*         3.4.2  REFLECT./TRANSMISSIVITY DUE TO RAYLEIGH AND AEROSOLS    
C                 ----------------------------------------------------    
C                                                                         
 3420 CONTINUE                                                            
C                                                                         
            ZGAP = ZCGAZ(JLON,JKM1)                                       
            ZBMU0 = 0.5 - 0.75 * ZGAP / ZMUE                              
            ZWW = ZPIZAZ(JLON,JKM1)                                       
            ZTO = ZTAUAZ(JLON,JKM1)                                       
            ZDEN = 1. + (1. - ZWW + ZBMU0 * ZWW) * ZTO * ZMUE             
     S           + (1-ZWW)*(1.-ZWW+2.*ZBMU0*ZWW)*ZTO*ZTO*ZMUE*ZMUE        
            ZRAY1(JLON,JKM1) = ZBMU0 * ZWW * ZTO * ZMUE / ZDEN            
            ZTRA1(JLON,JKM1) = 1. / ZDEN                                  
C                                                                         
            ZMU1 = 0.5                                                    
            ZBMU1 = 0.5 - 0.75 * ZGAP * ZMU1                              
            ZDEN1= 1. + (1. - ZWW + ZBMU1 * ZWW) * ZTO / ZMU1             
     S           + (1.-ZWW)*(1.-ZWW+2.*ZBMU1*ZWW)*ZTO*ZTO/ZMU1/ZMU1       
            ZRAY2(JLON,JKM1) = ZBMU1 * ZWW * ZTO / ZMU1 / ZDEN1           
            ZTRA2(JLON,JKM1) = 1. / ZDEN1                                 
C                                                                         
C*         3.4.3  EFFECT OF CLOUD LAYER                                   
C                 ---------------------                                   
C                                                                         
 3430 CONTINUE                                                            
C                                                                         
            PTAU(JLON,INU,JKM1) = MAX( PTAU(JLON,INU,JKM1) , ZEPSCT )     
            ZW(JLON) = POMEGA(JLON,INU,JKM1)                              
            ZTO1(JLON) = PTAU(JLON,INU,JKM1)/ZW(JLON)                     
     S               + ZTAUAZ(JLON,JKM1)/ZPIZAZ(JLON,JKM1)                
            ZR21(JLON) = PTAU(JLON,INU,JKM1) + ZTAUAZ(JLON,JKM1)          
            ZR22(JLON) = PTAU(JLON,INU,JKM1) / ZR21(JLON)                 
            ZGG(JLON) = ZR22(JLON) * PCG(JLON,INU,JKM1)                   
     S              + (1. - ZR22(JLON)) * ZCGAZ(JLON,JKM1)                
             ZW(JLON) = ZR21(JLON) / ZTO1(JLON)                           
            ZREF(JLON)=ZREFZ(JLON,1,JKM1)                                 
            ZRMUZ(JLON)=ZRMUE(JLON,JK)                                    
 342     CONTINUE                                                         
C                                                                         
      CALL DEDD(ZGG,ZREF,ZRMUZ,ZTO1,ZW,ZRE1,ZTR1,ZRE2,ZTR2)               
C                                                                         
                                                                          
         DO 345 JLON = 1 , KDLON                                          
C                                                                         
            ZREFZ(JLON,2,JK) = (1.-ZRNEB(JLON)) * (ZRAY1(JLON,JKM1)       
     S                     + ZREFZ(JLON,2,JKM1) * ZTRA1(JLON,JKM1)        
     S                     * ZTRA2(JLON,JKM1) )                           
     S                     + ZRNEB(JLON) * ZRE1(JLON)                     
C                                                                         
            ZTR(JLON,2,JKM1) = ZRNEB(JLON) * ZTR1(JLON)                   
     S                     + ZTRA1(JLON,JKM1) * (1.-ZRNEB(JLON))          
C                                                                         
            ZREFZ(JLON,1,JK) = (1.-ZRNEB(JLON)) * (ZRAY1(JLON,JKM1)       
     S                     + ZREFZ(JLON,1,JKM1) * ZTRA1(JLON,JKM1)        
     S                     * ZTRA2(JLON,JKM1)                             
     S                     / (1.-ZRAY2(JLON,JKM1)*ZREFZ(JLON,1,JKM1)))    
     S                     + ZRNEB(JLON)*ZRE2(JLON)                       
C                                                                         
            ZTR(JLON,1,JKM1) = ZRNEB(JLON) * ZTR2(JLON)                   
     S                     + (ZTRA1(JLON,JKM1)/ (1.-ZRAY2(JLON,JKM1)      
     S                     * ZREFZ(JLON,1,JKM1)))* (1.-ZRNEB(JLON))       
C                                                                         
 345     CONTINUE                                                         
 346  CONTINUE                                                            
C                                                                         
C*         3.5    REFLECT./TRANSMISSIVITY BETWEEN SURFACE AND LEVEL       
C                 -------------------------------------------------       
C                                                                         
 350  CONTINUE                                                            
C                                                                         
      DO 354 JABS = 1 , 2                                                 
         DO 351 JLON = 1 , KDLON                                          
            ZRJ(JLON,JABS,KFLEV+1) = 1.                                   
            ZRK(JLON,JABS,KFLEV+1) = ZREFZ(JLON,JABS,KFLEV+1)             
 351     CONTINUE                                                         
C                                                                         
         DO 353 JK = 1 , KFLEV                                            
            JKL = KFLEV+1 - JK                                            
            JKLP1 = JKL + 1                                               
            DO 352 JLON = 1 , KDLON                                       
               ZRE11 = ZRJ(JLON,JABS,JKLP1) * ZTR(JLON,JABS,JKL)          
               ZRJ(JLON,JABS,JKL) = ZRE11                                 
               ZRK(JLON,JABS,JKL) = ZRE11 * ZREFZ(JLON,JABS,JKL)          
 352        CONTINUE                                                      
 353     CONTINUE                                                         
 354  CONTINUE                                                            
C                                                                         
C*         3.6    REFLECT./TRANSMISSIVITY WITH GREY ABSORPTION            
C                 --------------------------------------------            
C                                                                         
 360  CONTINUE                                                            
C                                                                         
      JN = 2                                                              
C                                                                         
      DO 369 JABS=1,2                                                     
C                                                                         
         DO 361 JLON = 1 , KDLON                                          
            ZREFZ(JLON,2,1) = PALBS(JLON,INU)                             
            ZREFZ(JLON,1,1) = PALBS(JLON,INU)                             
 361     CONTINUE                                                         
C                                                                         
      DO 364 JK = 2 , KFLEV+1                                             
         JKM1 = JK - 1                                                    
         DO 362 JLON = 1 , KDLON                                          
            ZRNEB(JLON) = PCLDSW(JLON,JKM1)                               
            ZAA = ZUD(JLON,JABS,JKM1)                                     
            ZRKI = ZAKI(JLON,JABS)                                        
            ZS(JLON) = EXP(-ZRKI * ZAA * 1.66)                            
            ZG(JLON) = EXP(-ZRKI * ZAA / ZRMUE(JLON,JK))                  
            ZTR1(JLON) = 0.                                               
            ZRE1(JLON) = 0.                                               
            ZTR2(JLON) = 0.                                               
            ZRE2(JLON) = 0.                                               
C                                                                         
C*         3.6.1  INTRODUCING CLOUD EFFECTS                               
C                 -------------------------                               
C                                                                         
 3610 CONTINUE                                                            
C                                                                         
            PTAU(JLON,INU,JKM1) = MAX( PTAU(JLON,INU,JKM1) , ZEPSCT )     
            ZW(JLON)= POMEGA(JLON,INU,JKM1)                               
            ZTO1(JLON) = PTAU(JLON,INU,JKM1) / ZW(JLON)                   
     S               + ZTAUAZ(JLON,JKM1) / ZPIZAZ(JLON,JKM1)              
     S               + ZAA * ZRKI                                         
            ZR21(JLON) = PTAU(JLON,INU,JKM1) + ZTAUAZ(JLON,JKM1)          
            ZR22(JLON) = PTAU(JLON,INU,JKM1) / ZR21(JLON)                 
            ZGG(JLON) = ZR22(JLON) * PCG(JLON,INU,JKM1)                   
     S              + (1. - ZR22(JLON)) * ZCGAZ(JLON,JKM1)                
            ZW(JLON) = ZR21(JLON) / ZTO1(JLON)                            
            ZREF(JLON) = ZREFZ(JLON,1,JKM1)                               
            ZRMUZ(JLON) = ZRMUE(JLON,JK)                                  
 362     CONTINUE                                                         
C                                                                         
         CALL DEDD(ZGG,ZREF,ZRMUZ,ZTO1,ZW,ZRE1,ZTR1,ZRE2,ZTR2)            
C                                                                         
         DO 363 JLON = 1 , KDLON                                          
C                                                                         
            ZREFZ(JLON,2,JK) = (1.-ZRNEB(JLON)) * (ZRAY1(JLON,JKM1)       
     S                     + ZREFZ(JLON,2,JKM1) * ZTRA1(JLON,JKM1)        
     S                     * ZTRA2(JLON,JKM1) ) * ZG(JLON) * ZS(JLON)     
     S                     + ZRNEB(JLON) * ZRE1(JLON)                     
C                                                                         
            ZTR(JLON,2,JKM1)=ZRNEB(JLON)*ZTR1(JLON)                       
     S              + (ZTRA1(JLON,JKM1)) * ZG(JLON) * (1.-ZRNEB(JLON))    
C                                                                         
            ZREFZ(JLON,1,JK)=(1.-ZRNEB(JLON))*(ZRAY1(JLON,JKM1)           
     S             +ZREFZ(JLON,1,JKM1)*ZTRA1(JLON,JKM1)*ZTRA2(JLON,JKM1)  
     S             /(1.-ZRAY2(JLON,JKM1)*ZREFZ(JLON,1,JKM1)))             
     S             *ZG(JLON)*ZS(JLON)+ ZRNEB(JLON) * ZRE2(JLON)           
C                                                                         
            ZTR(JLON,1,JKM1)= ZRNEB(JLON) * ZTR2(JLON)                    
     S                    + (ZTRA1(JLON,JKM1)/(1.-ZRAY2(JLON,JKM1)        
     S                    * ZREFZ(JLON,1,JKM1)))                          
     S                    * ZG(JLON) * (1. -ZRNEB(JLON))                  
C                                                                         
C                                                                         
 363        CONTINUE                                                      
 364     CONTINUE                                                         
C                                                                         
C*         3.6.2  REFLECT./TRANSMISSIVITY BETWEEN SURFACE AND LEVEL       
C                 -------------------------------------------------       
C                                                                         
 3620 CONTINUE                                                            
C                                                                         
         DO 368 KREF=1,2                                                  
C                                                                         
            JN = JN + 1                                                   
C                                                                         
            DO 365 JLON = 1 , KDLON                                       
               ZRJ(JLON,JN,KFLEV+1) = 1.                                  
               ZRK(JLON,JN,KFLEV+1) = ZREFZ(JLON,KREF,KFLEV+1)            
 365        CONTINUE                                                      
C                                                                         
            DO 367 JK = 1 , KFLEV                                         
               JKL = KFLEV+1 - JK                                         
               JKLP1 = JKL + 1                                            
               DO 366 JLON = 1 , KDLON                                    
                  ZRE11 = ZRJ(JLON,JN,JKLP1) * ZTR(JLON,KREF,JKL)         
                  ZRJ(JLON,JN,JKL) = ZRE11                                
                  ZRK(JLON,JN,JKL) = ZRE11 * ZREFZ(JLON,KREF,JKL)         
 366           CONTINUE                                                   
 367        CONTINUE                                                      
 368     CONTINUE                                                         
 369  CONTINUE                                                            
C                                                                         
C*         3.7    UPWARD (ZRK) AND DOWNWARD (ZRJ) PSEUDO-FLUXES           
C                 ---------------------------------------------           
C                                                                         
 370  CONTINUE                                                            
C                                                                         
      DO 374 JK = 1 , KFLEV+1                                             
         DO 373 JAJ = 1 , 5 , 2                                           
            JAJP = JAJ + 1                                                
            DO 372 JLON = 1 , KDLON                                       
               ZRJ(JLON,JAJ,JK)= ZRJ(JLON,JAJ,JK) - ZRJ(JLON,JAJP,JK)     
               ZRK(JLON,JAJ,JK)= ZRK(JLON,JAJ,JK) - ZRK(JLON,JAJP,JK)     
               ZRJ(JLON,JAJ,JK)= MAX( ZRJ(JLON,JAJ,JK) , ZEELOG )         
               ZRK(JLON,JAJ,JK)= MAX( ZRK(JLON,JAJ,JK) , ZEELOG )         
 372        CONTINUE                                                      
 373     CONTINUE                                                         
 374  CONTINUE                                                            
C                                                                         
      DO 377 JK = 1 , KFLEV+1                                             
         DO 376 JAJ = 2 , 6 , 2                                           
            DO 375 JLON = 1 , KDLON                                       
               ZRJ(JLON,JAJ,JK)= MAX( ZRJ(JLON,JAJ,JK) , ZEELOG )         
               ZRK(JLON,JAJ,JK)= MAX( ZRK(JLON,JAJ,JK) , ZEELOG )         
 375        CONTINUE                                                      
 376     CONTINUE                                                         
 377  CONTINUE                                                            
C                                                                         
C*         3.8    EFFECTIVE ABSORBER AMOUNTS BY INVERSE LAPLACE           
C                 ---------------------------------------------           
C                                                                         
 380  CONTINUE                                                            
C                                                                         
      DO 387 JK = 1 , KFLEV+1                                             
         JKKI = 1                                                         
         DO 385 JAJ = 1 , 2                                               
            DO 384 JN = 1 , 2                                             
               JN2J = JN + 2 * JAJ                                        
               JKKP4 = JKKI + 4                                           
C                                                                         
C*         3.8.1  EFFECTIVE ABSORBER AMOUNTS                              
C                 ---------------------------------------------           
C                                                                         
 3810 CONTINUE                                                            
C                                                                         
C                                                                         
               DO 3811 JLON = 1 , KDLON                                   
                 ZW(JLON) = LOG( ZRJ(JLON,JN,JK) / ZRJ(JLON,JN2J,JK))     
     S                   / ZAKI(JLON,JAJ)                                 
 3811          CONTINUE                                                   
C                                                                         
C*         3.8.2  TRANSMISSION FUNCTION                                   
C                 ---------------------                                   
C                                                                         
 3820 CONTINUE                                                            
C                                                                         
                CALL SWTT ( INU, JAJ,APAD,BPAD,D,ZW,ZR1)                  
C                                                                         
                DO 3821 JLON = 1 , KDLON                                  
                   ZRL(JLON,JKKI) = ZR1(JLON)                             
                   ZRUEF(JLON,JKKI) = ZW(JLON)                            
                   ZW(JLON) = LOG( ZRK(JLON,JN,JK) / ZRK(JLON,JN2J,JK))   
     S                    / ZAKI(JLON,JAJ)                                
 3821           CONTINUE                                                  
C                                                                         
                CALL SWTT ( INU, JAJ, APAD,BPAD,D,ZW,ZR1)                 
C                                                                         
                DO 383 JLON = 1 , KDLON                                   
                   ZRL(JLON,JKKP4) = ZR1(JLON)                            
                   ZRUEF(JLON,JKKP4) = ZW(JLON)                           
 383            CONTINUE                                                  
C                                                                         
                JKKI=JKKI+1                                               
 384         CONTINUE                                                     
 385      CONTINUE                                                        
C                                                                         
C*         3.8.3  UPWARD AND DOWNWARD FLUXES WITH H2O AND UMG ABSORPTION  
C                 ------------------------------------------------------  
C                                                                         
 3830 CONTINUE                                                            
C                                                                         
          DO 386 JLON = 1 , KDLON                                         
            PFDOWN(JLON,JK) = ZRJ(JLON,1,JK)*ZRL(JLON,1)*ZRL(JLON,3)      
     S                      + ZRJ(JLON,2,JK)*ZRL(JLON,2)*ZRL(JLON,4)      
            PFUP(JLON,JK)   = ZRK(JLON,1,JK)*ZRL(JLON,5)*ZRL(JLON,7)      
     S                      + ZRK(JLON,2,JK)*ZRL(JLON,6)*ZRL(JLON,8)      
 386      CONTINUE                                                        
 387  CONTINUE                                                            
C                                                                         
C*         3.9    INTRODUCTION OF OZONE ABSORPTION                        
C                 --------------------------------                        
C                                                                         
 390  CONTINUE                                                            
C                                                                         
      JABS=3                                                              
      DO 395 JK = 1 , KFLEV+1                                             
         DO 392 JLON = 1 , KDLON                                          
            ZW(JLON) = ZUD(JLON,JABS,JK)                                  
 392     CONTINUE                                                         
C                                                                         
         CALL SWTT ( INU, JABS, APAD,BPAD,D,ZW,ZR1)                       
C                                                                         
         DO 393 JLON = 1 , KDLON                                          
            PFDOWN(JLON,JK) = ZR1(JLON) * PFDOWN(JLON,JK) * SUN(INU)      
            ZW(JLON) = ZUM(JLON,JK)                                       
 393     CONTINUE                                                         
C                                                                         
         CALL SWTT ( INU, JABS,APAD,BPAD,D,ZW,ZR1)                        
C                                                                         
         DO 394 JLON = 1 , KDLON                                          
            PFUP(JLON,JK) = ZR1(JLON) * PFUP(JLON,JK) * SUN(INU)          
 394     CONTINUE                                                         
 395  CONTINUE                                                            
C                                                                         
C     ------------------------------------------------------------        
C                                                                         
C*         4.     NET TOTAL SHORTWAVE FLUXES                              
C                 --------------------------                              
C                                                                         
C 400  CONTINUE                                                           
CZ      DO 402 JK = 1 , KFLEV+1                                           
C         DO 401 JLON = 1 , KDLON                                         
C           PFUP(JLON,JK) = (PFUP(JLON,JK) +ZFU(JLON,JK)) * ZFACT(JLON)   
C           PFDOWN(JLON,JK)= (PFDOWN(JLON,JK)+ZFD(JLON,JK))*ZFACT(JLON)   
C 401     CONTINUE                                                        
C 402  CONTINUE                                                           
C                                                                         
C*         4.     NET TOTAL SHORTWAVE FLUXES                              
C                 --------------------------                              
C                                                                         
 400  CONTINUE                                                            
      DO 402 JK = 1 , KFLEV+1                                             
         DO 401 JLON = 1 , KDLON                                          
c            IF(ITASK.EQ.1)THEN                                           
c              PFUP(JLON,JK)   = ZFU(JLON,JK) * ZFACT(JLON)               
c              PFDOWN(JLON,JK) = ZFD(JLON,JK) * ZFACT(JLON)               
c            ELSE IF(ITASK.EQ.2)THEN                                      
c              PFUP(JLON,JK)   = PFUP(JLON,JK) * ZFACT(JLON)              
c              PFDOWN(JLON,JK) = PFDOWN(JLON,JK) * ZFACT(JLON)            
c            ELSE IF(ITASK.EQ.3)THEN                                      
              PFUP(JLON,JK)=(PFUP(JLON,JK)+ZFU(JLON,JK))*ZFACT(JLON)      
              PFDOWN(JLON,JK) = (PFDOWN(JLON,JK) + ZFD(JLON,JK)           
     *                      + FOXDWN(JLON,JK)) * ZFACT(JLON)              
c            ENDIF                                                        
 401     CONTINUE                                                         
 402  CONTINUE                                                            
C                                                                         
      RETURN                                                              
      END                                                                 
