      SUBROUTINE MAKECLOUDS(poslats,mg,jg,nl,p0,sigma,tauaerosol)

!      *****************************************************************
!      * This routine generates a an array of aerosol optical depths   *
!      * It returns this lat x lon x layer array and writes it to file *
!      * The file also includes pi0 and asymmetry parameter at both    *
!      * short wave and long wave channels.  ~2016 MTR                 *
!      *****************************************************************

       REAL PRESSURE(NL),PABSDIF1(NL),PABSDIF2(NL),
     &      VERTPROF(NL),LONGYS(MG),LATYS(JG),THELAT,TheTAU,
     &      TOTAL,TAUPA,TAUPB,TAUAEROSOL(nl,mg,2,jg),TAUC,SIGC,
     &      deltalonc,deltalonc360,SIGMA(NL),themin(1),POSLATS(JG)
       INTEGER PINDEXBOTC,PINDEXTOPC,NLAT,NLON,NLEV,NHEM
     &         ,IHEM, ILAT,ILEV,IL,ILON,TheCounter

       COMMON/CLOUDY/AEROSOLMODEL,AERTOTTAU,CLOUDBASE,
     &   CLOUDTOP,AERHFRAC,PI0AERSW,ASYMSW,EXTFACTLW,PI0AERLW,
     &   ASYMLW,SIG_AREA,PHI_LON,AERO4LAT,AEROPROF
       CHARACTER(15) :: AEROSOLMODEL

       NAMELIST/INCLOUDY/AEROSOLMODEL,AERTOTTAU,CLOUDBASE,
     &  CLOUDTOP,AERHFRAC,PI0AERSW,ASYMSW,EXTFACTLW,PI0AERLW,ASYMLW,
     &  SIG_AREA,PHI_LON
      

       READ (7,INCLOUDY)
       
!      DIMENSIONS--------------------
       nlat     =  jg
       nlon     =  mg
       nlev     =  nl
       nhem     =   2
       pressure = sigma*P0
!            Define longitude array:
           DO    ilon  = 1,NLON
           LONGYS(ILON) = REAL(ilon-1)/REAL(NLON)*360.0
           ENDDO
!            Define latitude array:
           LATYS=POSLATS 
!      CLOUD DISTRIBUTION------------

!         1. VERTICAL
!          cloudbase and cloudtop pressures read in from fort.7 cloudy namelist
!          Find the closest sigma-pressure levels to these pressures:
          DO IL=1,NL 
            PABSDIF1(IL)=ABS(PRESSURE(IL)*1e-5 - cloudbase)
            PABSDIF2(IL)=ABS(PRESSURE(IL)*1e-5 - cloudtop) 
            VERTPROF(IL) = 0
          ENDDO
           themin=MINLOC(PABSDIF1)
           PINDEXBOTC=themin(1)
           themin=MINLOC(PABSDIF2)
           PINDEXTOPC=themin(1)
!         Now create a vertical array that equals one at the cloud base pressure,
!         scales with pressure (~uniform mixing ratio) to the power of
!         aerHfrac (e.g. aerHfrac = 1 is scales 1:1 with pressure) until the cloud top
!         pressure, and is zero everywhere else.
          TOTAL=0.0
          DO IL=PINDEXTOPC,PINDEXBOTC
           VERTPROF(IL) = (PRESSURE(IL)/PRESSURE(PINDEXBOTC))**AERHFRAC
             TOTAL=TOTAL+VERTPROF(IL)
          ENDDO
!         Normalize it to preserve the total integrated optical.
          DO IL=PINDEXTOPC,PINDEXBOTC
           VERTPROF(IL)=VERTPROF(IL)/TOTAL
          ENDDO
!         Now we have a dimensionless array that integrates to one and defines
!         the shape of the cloud in the vertical dimension.
 

!         2. HORIZONTAL
!      Specify scattering optical depth as a function of latitude, longitude---and using the array above--height.

!      For Kepler 7-b, a 2-d parameterization of Munoz & Isaak (2015) is used:
!             tau(lon,lat;tauc,sigc,deltalonc)=tc*exp[-({lon-deltalonc}^2+lat^2)/(2*sigidlc^2)]
!             deltalonc is the eastward offset relative to the substeller point, 
!             tauc is cloud thickness, cloud width is sigc.
       tauc=AERTOTTAU
       sigc=SIG_AREA
       deltalonc=360+(-PHI_LON)
       dellonc360=deltalonc+360.
       
 

!      PREPARE TO WRITE THIS CLOUD INFORMATION TO FILE FOR THE RECORD, THOUGH
!      NOT FOR THE SUBSEQUENT COMPUTATIONS. FILE WILL BE FORT.77 
       WRITE(77,*) 'CLOUD MODEL'
       WRITE(77,*) '' 
       WRITE(77,*) 'NAME: ',AEROSOLMODEL
       WRITE(77,*) 'Parameters:'
       WRITE(77,*) 'cloudbase,    cloudtop (bars),     powerlaw:'
       WRITE(77,7010) cloudbase,cloudtop,aerHfrac
 7010  FORMAT(1x,F11.4,2x,F11.4,3x,F7.3)
       WRITE(77,*) 'Short Wave Scattering--' 
       WRITE(77,*) '    Total optical depth:', TAUC
       WRITE(77,*) '    Single scattering albedo (PI0):',PI0AERSW
       WRITE(77,*) '    Asymmetry parameter:',ASYMSW
       WRITE(77,*) 'Long Wave Scattering--'
       WRITE(77,*) '    Extinction Ratio (LW/SW):',EXTFACTLW
       WRITE(77,*) '    Single scattering albedo (PI0):',PI0AERLW
       WRITE(77,*) '    Asymmetry parameter:',ASYMLW
       WRITE(77,*) 'Other parameters:'
       WRITE(77,*) 'sig_area =',sigc
       write(77,*) 'phi_lon =',phi_lon
       write(77,*) ''
       write(77,*) 'NORMALIZED VERTICAL PROFILE:'
       WRITE(77,*) '   PRESSURE(MBAR)        NORMALIZED AEROSOL TAU'
              DO IL = 1, NLEV
              WRITE(77,7712) PRESSURE(IL)*1e-5,VERTPROF(IL)
 7712         FORMAT(1x,F11.6,3x,F12.7)
              ENDDO
       WRITE(77,*)''
       write(77,*) 'FULL GRID:'
       WRITE(77,*) '' 
        thecounter = 0 !counter
       DO ILAT   =1,JG
          DO IHEM  =1,2
!            DEFINE THE LAT, GIVEN THE LAT INDEX AND HEMISPHERE
             THELAT=(LATYS(ILAT)*(-1.)*(-1.)**(IHEM))
             DO ILON = 1,MG
              
!            FIRST WRITE THE LAT & LON
              WRITE(77,*) 'LATITUDE, LONGITUDE:',THELAT,LONGYS(ILON)
              WRITE(77,771) '1)PRESS(BARS)','2)AEROSOL_SW_TAU',
     &         '3)SW_PI0','4)SW_ASYM',
     &         '5)AEROSOL_LW_TAU','6)LW_PI0','7)LW_ASYM'
 771       FORMAT(1X,A13,1X,A16,4X,A8,4X,A9,4x,A16,4x,A8,4x,A9)  
                DO ILEV= 1,NL

!    HERE'S WHERE THE CHOICE OF MODEL SPECIFIED IN FORT.7 IS
!    USED TO CALCULATE THE CLOUD SPATIAL COVERAGE
               
!         KEPLER 7B DAYSIDE CLOUDS
                IF(AEROSOLMODEL.EQ.'Kepler7b') THEN
                
                  IF (thecounter.EQ.0) THEN 
                   write(*,*) 'Using Cloudmodel: Kepler7b'
                   thecounter=1.
                  ENDIF
                TAUpa=(TAUC
     & *EXP(-(((LONGYS(ILON)-DELTALONC)*(LONGYS(ILON)-DELTALONC))
     &     +(THELAT*THELAT))/(2.*SIGC*SIGC)))*VERTPROF(ILEV)
! ONCE AGAIN TO FILL OUT THE CIRCLE THAT THE INDEXING WRAPAROUND MISSES
                TAUpb=(TAUC
     & *EXP(-(((LONGYS(ILON)-DELLONC360)*(LONGYS(ILON)-DELLONC360))
     &     +(THELAT*THELAT))/(2.*SIGC*SIGC)))*VERTPROF(ILEV)
! SUM IT UP
                     TheTAU=TAUpa+TAUpb
                     IF (LONGYS(ILON).EQ.360) THEN
                     TheTAU=TAUpb
                     ENDIF
              
!                ENDIF FIRST CONDITIONAL, KEPLER 7B
!          NIGHTSIDE CLOUDS
                ELSE IF (AEROSOLMODEL.EQ.'Nightside') THEN
                  IF (thecounter.EQ.0) THEN   
                    WRITE(*,*) 'Using Cloudmodel: Nightside'
                   thecounter=1.
                  ENDIF

                    IF(LONGYS(ILON).GE.90.) THEN
                       IF(LONGYS(ILON).LE.270.) THEN
                          TheTAU=TAUC*VERTPROF(ILEV)
                       ENDIF 
                    ENDIF
           
!               END ELSE NIGHTSIDE CLOUDS
!          GLOBAL CLOUD COVERAGE
                ELSE IF (AEROSOLMODEL.EQ.'Global') THEN
                      IF (thecounter.EQ.0) THEN
                       WRITE(*,*) 'Using Cloudmodel: Global'
                       thecounter=1.
                      ENDIF

                          TheTAU=TAUC*VERTPROF(ILEV)
!               END ELSE GLOBAL CLOUDS
!          NOTHING SPECIFIED 
                ELSE                 
                    WRITE(*,*)'NO VALID CLOUDMODEL SPECIFIED! STOPPING'
                    STOP
                END IF

 
                     TAUAEROSOL(ILEV,ILON,IHEM,ILAT)=TheTAU
             WRITE(77,772) PRESSURE(ILEV)*1E-5,TheTAU,PI0AERSW,ASYMSW
     &                         TheTAU*EXTFACTLW,PI0AERLW,ASYMLW
 772       FORMAT(2X,F11.6,3X,F12.7,3X,F7.4,3X,F7.4,3X,F11.6
     &             ,3X,F7.4,3X,F7.4)
                ENDDO  
             WRITE(77,*)''
             ENDDO
          ENDDO
       ENDDO

       END
