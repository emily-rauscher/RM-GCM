      SUBROUTINE MAKECLOUDS(poslats,p0,sigma)
!      SUBROUTINE MAKECLOUDS(poslats,mg,jg,nl,p0,sigma)
!      *****************************************************************
!      * This routine generates an array of aerosol optical depths.    *
!      * It returns this (layer x lon x hem x lat) array and writes it *
!      * to file (fort.60). The file also includes pi0 and asymmetry   *
!      * parameter at both short wave and long wave channels.          *
!      *****************************************************************
       include 'params.i'
 
       REAL PRESSURE(NL+1),PABSDIF1(NL+1),PABSDIF2(NL+1),
     &     VERTPROF(NL+1),PSIGMA(NL),LONGYS(MG),LATYS(JG),THELAT,TheTAU,
     &     TOTAL,TAUPA,TAUPB,TAUC,SIGC,XFACTSW,XFACTLW,nightedge,TAUpaN,
     &     deltalonc,deltalonc360,SIGMA(NL),themin(1),POSLATS(JG)
       INTEGER PINDEXBOTC,PINDEXTOPC,NLAT,NLON,NLEV,NHEM
     &         ,IHEM, ILAT,ILEV,IL,ILON,TheCounter,LD

       CHARACTER(30) :: AEROSOLMODEL
       REAL TAUAEROSOL(nl+1,mg,2,jg),AEROPROF(NL+1)
       LOGICAL DELTASCALE
       COMMON/CLOUDY/AEROSOLMODEL,AERTOTTAU,CLOUDBASE,
     &   CLOUDTOP,CLDFRCT,AERHFRAC,PI0AERSW,ASYMSW,EXTFACTLW,PI0AERLW,
     &   ASYMLW,DELTASCALE,SIG_AREA,PHI_LON,TAUAEROSOL,AEROPROF

       NAMELIST/INCLOUDY/AEROSOLMODEL,AERTOTTAU,CLOUDBASE,
     &  CLOUDTOP,CLDFRCT,AERHFRAC,PI0AERSW,ASYMSW,EXTFACTLW,PI0AERLW,
     &  ASYMLW,DELTASCALE,SIG_AREA,PHI_LON
      

       READ (7,INCLOUDY)
       
!      DIMENSIONS--------------------
       nlat     =  jg
       nlon     =  mg
       nlev     =  nl+1 !To define the layer above the top boundary 
!       nhem     =   2
       psigma   = sigma*P0
       ! LAYER EDGES AT WHICH FLUXES ARE COMPUTED, PRFLUX.
             DO LD    = 1,NL-1
             pressure(LD+1)=(psigma(LD)+psigma(LD+1))/2.
             ENDDO
             pressure(NL+1)=p0
             pressure(1)=psigma(1)*0.5
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
          DO IL=1,NLEV
            PABSDIF1(IL)=ABS(PRESSURE(IL)*1e-5 - cloudbase)
            PABSDIF2(IL)=ABS(PRESSURE(IL)*1e-5 - cloudtop) 
            VERTPROF(IL) = 0
          ENDDO
           themin=MINLOC(PABSDIF1)
           PINDEXBOTC=themin(1)
           themin=MINLOC(PABSDIF2)
           PINDEXTOPC=themin(1)
           IF (PINDEXTOPC.GT.PINDEXBOTC) THEN 
       WRITE(*,*) 'CLOUD TOP MUST BE LOWER PRESSURE THAN BASE!'
       WRITE(*,*) 'ABORTING! PLEASE CHOOSE GREATER RANGE IN FORT.7'     
           STOP
           ENDIF  
!         Now create a vertical array that equals one at the cloud base pressure,
!         scales with pressure (~uniform mixing ratio) to the power of
!         aerHfrac (e.g. aerHfrac = 1 is scales 1:1 with pressure) until the cloud top
!         pressure, and is zero everywhere else.
          TOTAL=0.0
          DO IL=PINDEXTOPC,PINDEXBOTC
           VERTPROF(IL) = (PRESSURE(IL)/PRESSURE(PINDEXBOTC))**AERHFRAC
          ENDDO

!         Smooth it by averaging with neighbors
          DO IL=PINDEXTOPC,PINDEXBOTC
           ABOVE        = MAX(IL-1,1) 
           BELOW        = MIN(IL+1,NLEV)
       VERTPROF(IL) = (VERTPROF(ABOVE)+VERTPROF(IL)+VERTPROF(BELOW))/3.
           TOTAL        = TOTAL+VERTPROF(IL)
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
       deltalonc=360.+(-PHI_LON)
       dellonc360=deltalonc-360.
       nightedge =360.-90.

!      PREPARE TO WRITE THIS CLOUD INFORMATION TO FILE FOR THE RECORD, THOUGH
!      NOT FOR THE SUBSEQUENT COMPUTATIONS. FILE WILL BE FORT.61 
       WRITE(61,*) 'CLOUD MODEL'
       WRITE(61,*) '' 
       WRITE(61,*) 'NAME: ',AEROSOLMODEL
       WRITE(61,*) 'Parameters:'
       WRITE(61,*) 'cloudbase, cloudtop(bars), cloudfraction, powerlaw:'
       WRITE(61,7010) cloudbase,cloudtop,cldfrct,aerHfrac
 7010  FORMAT(1x,F11.4,2x,F11.4,3x,F7.3,3x,F7.3)
       WRITE(61,*) 'AEROSOL SCATTERING PARAMETERS:'
       WRITE(61,*) 'Short Wave Scattering--' 
       WRITE(61,*) '    Total optical depth:', TAUC
       WRITE(61,*) '    Single scattering albedo (PI0):',PI0AERSW
       WRITE(61,*) '    Asymmetry parameter:',ASYMSW
       WRITE(61,*) 'Long Wave Scattering--'
       WRITE(61,*) '    Extinction Ratio (LW/SW):',EXTFACTLW
       WRITE(61,*) '    Single scattering albedo (PI0):',PI0AERLW
       WRITE(61,*) '    Asymmetry parameter:',ASYMLW
       WRITE(61,*) 'NOTE: These are aerosol values only,(i.e. no gas)'
       write(61,*) '      not total layer values.'
            IF(DELTASCALE) THEN 
            WRITE(61,*) 'Delta-scaling applied'
            ELSE
            WRITE(61,*) 'No Delta-scaling applied'
            ENDIF
       WRITE(61,*) 'Other parameters:'
       WRITE(61,*) 'sig_area =',sigc
       write(61,*) 'phi_lon =',phi_lon
       write(61,*) ''
       write(61,*) 'NORMALIZED VERTICAL PROFILE:'
       WRITE(61,*) '   PRESSURE(MBAR)        NORMALIZED AEROSOL TAU'
              DO IL = 1, NLEV
              WRITE(61,2112) PRESSURE(IL)*1e-5,VERTPROF(IL)
 2112         FORMAT(1x,F11.6,3x,F12.7)
              ENDDO
       WRITE(61,*)''
       write(61,*) 'FULL GRID:'
       WRITE(61,*) '' 
        thecounter = 0 !counter
       DO ILAT   =1,JG
          DO IHEM  =1,2
!            DEFINE THE LAT, GIVEN THE LAT INDEX AND HEMISPHERE
             THELAT=(LATYS(ILAT)*(-1.)*(-1.)**(IHEM))
             DO ILON = 1,MG
              
!            FIRST WRITE THE LAT & LON
              WRITE(61,*) 'LATITUDE, LONGITUDE:',THELAT,LONGYS(ILON)
              WRITE(61,211) '1)PRESS(BARS)','2)AEROSOL_SW_TAU',
     &         '3)SW_PI0','4)SW_ASYM',
     &         '5)AEROSOL_LW_TAU','6)LW_PI0','7)LW_ASYM'
 211       FORMAT(1X,A13,1X,A16,4X,A8,4X,A9,4x,A16,4x,A8,4x,A9)  
                DO ILEV= 1,NLEV

!    HERE'S WHERE THE CHOICE OF MODEL SPECIFIED IN FORT.7 IS
!    USED TO CALCULATE THE CLOUD SPATIAL COVERAGE
               
!CASE #
!#########
!#########
!CAS1

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
                 TAUpb=0.0
!             write(*,*) (LONGYS(ILON)-DELTALONC),lONGYS(ILON),DELTALONC
                TAUpb=(TAUC
     & *EXP(-(((LONGYS(ILON)-DELLONC360)*(LONGYS(ILON)-DELLONC360))
     &     +(THELAT*THELAT))/(2.*SIGC*SIGC)))*VERTPROF(ILEV)
! SUM IT UP
                     TheTAU=TAUpa+TAUpb
                     IF (LONGYS(ILON).EQ.360) THEN
                     TheTAU=TAUpb
                     ENDIF
              
!                ENDIF FIRST CONDITIONAL, KEPLER 7B


!#########
!CASE 2
!          NOW NIGHTSIDE CLOUDS
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
!               END OF THE NIGHTSIDE



!##########
!CASE 3
!               NOW NIGHTSIDE, BUT WITH A GENTLER TRANSITION

                ELSE IF (AEROSOLMODEL.EQ.'Nightside_soft') THEN
                  IF (thecounter.EQ.0) THEN
                    WRITE(*,*) 'Using Cloudmodel: Nightside_soft'
                   thecounter=1.
                  ENDIF

                    IF(LONGYS(ILON).GT.90.) THEN
                       IF(LONGYS(ILON).LT.270.) THEN

                   DELTALONC = 180.
                   SIGC=45.
                   TAUpa=(TAUC
     & *EXP(-(((LONGYS(ILON)-DELTALONC)*(LONGYS(ILON)-DELTALONC))
     &     +(THELAT*THELAT))/(2.*SIGC*SIGC)))*VERTPROF(ILEV)
                   TheTAU=taupa
                       ENDIF
                    ENDIF
           
!               END ELSE NIGHTSIDE CLOUDS


!##########
!CASE 4
!          GLOBAL CLOUD COVERAGE
                ELSE IF (AEROSOLMODEL.EQ.'Global') THEN
                      IF (thecounter.EQ.0) THEN
                       WRITE(*,*) 'Using Cloudmodel: Global'
                       thecounter=1.
                      ENDIF

                          TheTAU=TAUC*VERTPROF(ILEV)

!               END ELSE GLOBAL CLOUDS


!###########
!CASE 5
!            MUNOZ K7B + NIGHSIDE
                ELSE IF (AEROSOLMODEL.EQ.'Kepler7b_nightsidex') THEN
                     IF (thecounter.EQ.0) THEN
                      WRITE(*,*) 'Using Cloudmodel: Kepler7b_nightsidex'
                       thecounter=1.
                     ENDIF

! For this case, first create the Munoz distribution
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
! OK, now replace the night side completely
                    IF(LONGYS(ILON).GT.90.) THEN
                       IF(LONGYS(ILON).LE.DELTALONC) THEN
!                          TheTAU=TAUC*VERTPROF(ILEV)
!                   TAUpaN=(TAUC
!     & *EXP(-(((LONGYS(ILON)-180.)*(LONGYS(ILON)-180.))
!     &     +(THELAT*THELAT))/(2.*40.*40.)))*VERTPROF(ILEV)
!                   TheTAU=min(TaupaN+TheTau,tauc)
                    TheTau=(TAUC
     & *EXP(-(THELAT*THELAT)/(2.*SIGC*SIGC)))*VERTPROF(ILEV)   
!     &     +(TAUC
!     & *EXP(-(((270.-DELLONC360)*(270-DELLONC360))
!     &     +(THELAT*THELAT))/(2.*SIGC*SIGC)))*VERTPROF(ILEV) 
                       ENDIF
                    ENDIF
! Finally, now remove a chunk for the antisymmetric western limb of
! night side.
                  TAUpa=(TAUC
     & *EXP(-(((LONGYS(ILON)-nightedge+180.)*
     &         (LONGYS(ILON)-nightedge+180.))
     &     +(THELAT*THELAT))/(2.*SIGC*SIGC)))*VERTPROF(ILEV)

!                  TAUpa=(TAUC
!     & *EXP(-(((LONGYS(ILON)-nightedge+180.)*
!     &         (LONGYS(ILON)-nightedge+180.))
!     &     +(THELAT*THELAT))/(2.*SIGC*SIGC)))*VERTPROF(ILEV)
! ONCE AGAIN TO FILL OUT THE CIRCLE THAT THE INDEXING WRAPAROUND MISSES
!                TAUpb=(TAUC
!     & *EXP(-(((LONGYS(ILON)-DELLONC360+180.)*
!     &        (LONGYS(ILON)-DELLONC360+180.))
!     &     +(THELAT*THELAT))/(2.*SIGC*SIGC)))*VERTPROF(ILEV) 

! Subtract this shifted cloud from the night side                   
                    IF(LONGYS(ILON).GE.90.) THEN
                       IF(LONGYS(ILON).LE.270.) THEN
                          TheTAU=TheTau-TAUpa
                          TheTAU=MAX(TheTAU,0.0)
                       ENDIF
                    ENDIF
               




!###############
!CASE 6
!
!            UNIFORM... AEROSOLS AT ALL HEIGHTS AND LOCATIONS
                ELSE IF (AEROSOLMODEL.EQ.'Uniform') THEN
                     IF (thecounter.EQ.0) THEN
                      WRITE(*,*) 'Using Cloudmodel: Uniform'
                       thecounter=1.
                     ENDIF       
                         TheTAU=TAUC   !all heights

!#############
!NO CASE

!          NOTHING SPECIFIED 
                ELSE                 
                    WRITE(*,*)'NO VALID CLOUDMODEL SPECIFIED! STOPPING'
                    STOP
                END IF

!          HERE WE WRITE THE CLOUD ARRAY TO MEMORY 
             TAUAEROSOL(ILEV,ILON,IHEM,ILAT)=TheTAU
!          AND WRITE IT TO FILE FORT.61
             XFACTSW = 1.0
             XFACTLW = 1.0
             IF (TheTAU.LT.1e-10) THEN 
             XFACTSW = 0.0 
             ENDIF
             IF (TheTAU*EXTFACTLW.LT.1e-6) THEN
             XFACTLW = 0.0
             ENDIF             

             WRITE(61,212) PRESSURE(ILEV)*1E-5,TheTAU,PI0AERSW*XFACTSW,
     &                     ASYMSW*XFACTSW,TheTAU*EXTFACTLW,
     &                     PI0AERLW*XFACTLW,ASYMLW*XFACTLW
 212       FORMAT(2X,F11.6,3X,F12.7,3X,F7.4,3X,F7.4,3X,F11.6
     &             ,3X,F7.4,3X,F7.4)
                ENDDO  
             WRITE(61,*)''
             ENDDO
          ENDDO
       ENDDO
!      WRITE TO FORT.60 FOR ALL INCLUSIVE RADIATIVE TRANSFER SUMMARY
!       write(*,*) 'where is this writing?'
!       write(60,*) 'ok'
!       WRITE(60,*) 'NAMELIST/INCLOUDY/'
!       WRITE(60,*) 'AEROSOLMODEL',AEROSOLMODEL
!       WRITE(60,*)'AERTOTTAU',AERTOTTAU
!       WRITE(60,*)'CLOUDBASE',CLOUDBASE
!       WRITE(60,*) 'CLOUDTOP',CLOUDTOP
!       WRITE(60,*) 'AERHFRAC',AERHFRAC
!       WRITE(60,*) 'PI0AERSW',PI0AERSW
!       WRITE(60,*)'ASYMSW',ASYMSW
!       WRITE(60,*)'EXTFACTLW',EXTFACTLW
!      WRITE(60,*)'PI0AERLW',PI0AERLW
!       WRITE(60,*)'ASYMLW',ASYMLW
!       WRITE(60,*)'DELTASCALE',DELTASCALE
!       WRITE(60,*)'SIG_AREA',SIG_AREA
!       WRITE(60,*)'PHI_LON',PHI_LON
       END
