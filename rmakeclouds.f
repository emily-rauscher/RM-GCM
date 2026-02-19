      SUBROUTINE MAKECLOUDS(poslats,p0,sigma)
!      SUBROUTINE MAKECLOUDS(poslats,mg,jg,nl,p0,sigma)
!      *****************************************************************
!      * This routine generates an array of aerosol optical depths.    *
!      * It returns this (layer x lon x hem x lat) array and writes it *
!      * to file (fort.60). The file also includes pi0 and asymmetry   *
!      * parameter at both short wave and long wave channels.          *
!      * Current cloud models:
!      *    Kepler7b
!      *    Nightside
!      *    Kepler7b_nightsidex
!      *    Nightside_soft
!      *    Global
!      *    Uniform
!      *****************************************************************
       include 'params.i'

       REAL PRESSURE(NL+1),PABSDIF1(NL+1),PABSDIF2(NL+1),
     &     VERTPROF(NL+1),PSIGMA(NL),LONGYS(MG),LATYS(JG),THELAT,TheTAU,
     &     TOTAL,TAUPA,TAUPB,TAUC,SIGC,XFACTSW,XFACTLW,nightedge,TAUpaN,
     &     deltalonc,deltalonc360,SIGMA(NL),themin(1),POSLATS(JG)
       INTEGER PINDEXBOTC,PINDEXTOPC,NLAT,NLON,NLEV,NHEM
     &         ,IHEM, ILAT,ILEV,IL,ILON,TheCounter,LD,MAXTAULOC

       CHARACTER(30) :: AEROSOLMODEL
       CHARACTER(30) :: AEROSOLCOMP
       REAL TAUAEROSOL(nl+1,mg,2,jg),AEROPROF(NL+1),MAXTAU,TCON(NL+1)
       REAL TCONMnS(NL+1),TCONSiO2(NL+1),TCONMg2SiO4(NL+1)
       REAL TCONAl2O3(NL+1),TCONFe(NL+1),TCONMgSiO3(NL+1),MOLEF(13)
       REAL MTLX, METALLICITY
       INTEGER AERLAYERS
       LOGICAL DELTASCALE, HAZES, PICKET_FENCE_CLOUDS
       COMMON/CLOUDY/AEROSOLMODEL,AERTOTTAU,CLOUDBASE,
     &   CLOUDTOP,CLDFRCT,AERHFRAC,PI0AERSW,ASYMSW,EXTFACTLW,PI0AERLW,
     &   ASYMLW,DELTASCALE,SIG_AREA,PHI_LON,TAUAEROSOL,AEROPROF,
     &   MAXTAU,MAXTAULOC,TCON,AEROSOLCOMP,MTLX,METALLICITY,HAZES,PICKET_FENCE_CLOUDS,MOLEF,AERLAYERS

       NAMELIST/INCLOUDY/AEROSOLMODEL,AEROSOLCOMP,MTLX,METALLICITY,HAZES,PICKET_FENCE_CLOUDS,MOLEF,AERLAYERS,
     &  AERTOTTAU,CLOUDBASE,CLOUDTOP,CLDFRCT,AERHFRAC,PI0AERSW,
     &  ASYMSW,EXTFACTLW,PI0AERLW,ASYMLW,DELTASCALE,SIG_AREA,PHI_LON

!     MnS
!      DATA TCONMnS /1113.43, 1120.07,1126.64,1133.21,1139.9,1146.46,
!     &        1153.15,  1159.65, 1166.25, 1173.01, 1180.09, 1187.73,
!     &        1194.69, 1202.12, 1209.26, 1216.50, 1223.92, 1231.43,
!     &       1239.19, 1247.17, 1255.39,1263.78, 1272.32,1281.40,
!     &        1289.69, 1297.76, 1306.00, 1315.24, 1324.18, 1333.27,
!     &        1342.44, 1351.39, 1360.54, 1369.99, 1379.72, 1389.42,
!     &        1399.22, 1409.04,  1418.99,1428.77, 1438.60, 1449.11,
!     &        1459.19, 1469.78, 1481.06, 1492.70, 1504.21, 1515.49,
!     &            1527.84,1540.17,1545.90/

!    SiO2
!      DATA TCONSiO2 /1334.63,1342.58,1350.30,1358.48,1366.64,1374.85,
!     &             1383.15,1391.59,1400.10,1408.68, 1417.25,1425.87,
!     &             1434.53, 1443.14, 1451.71,1460.28,1468.90,1477.44,
!     &             1486.12,1494.77,1503.91,1513.72 ,1524.07, 1534.57,
!     &             1544.98,1555.25, 1565.35,1575.40, 1585.44,1595.55,
!     &             1606.04,1617.06, 1628.55,1640.25,1651.98, 1664.07,
!     &             1676.82,1689.97,1703.53, 1717.16, 1731.36, 1746.01,
!     &             1761.10,1776.94, 1793.70,1811.36,1829.60, 1848.37,
!     &             1867.54, 1887.00, 1896.46/

!    MgSiO3
!      DATA   TCONMgSiO3 /1307.11,1314.43,1321.68,1329.11,1336.79,
!     &             1344.72,1352.04,1359.52,1367.46,1375.60,1383.27,
!     &             1391.60,1400.24,1408.51,1416.97,1425.28,1434.14,
!     &             1442.49,1451.31,1460.37,1469.43,1478.40,1487.67,
!     &             1497.14,1507.25,1517.04,1527.28,1537.29,1546.58,
!     &             1556.41,1567.50,1578.63,1588.50,1598.45,1609.41,
!     &             1620.62,1631.41,1642.37,1654.32,1666.68,1678.34,
!     &             1689.61,1702.01,1715.51,1727.89,1739.21,1751.75,
!     &             1765.34, 1778.58,1790.56,1796.07/

!    Mg2SiO4
!      DATA   TCONMg2SiO4 /1370.00,1378.31,1386.62,1395.03,1403.56,
!     &           1412.29,
!     &           1421.17, 1430.18, 1439.17, 1448.16, 1457.24,1466.52,
!     &           1475.94,1485.57, 1495.23, 1505.09, 1515.04, 1525.21,
!     &           1535.33,1545.80, 1556.37, 1567.12,1578.02, 1589.13,
!     &           1600.29, 1611.67,1623.20, 1634.84,1646.61,1658.58,
!!     &           1670.74, 1683.03, 1695.59, 1708.44,1721.29,1734.41,
!     &           1747.86, 1761.64, 1775.79, 1789.95, 1804.36,1819.11,
!     &           1834.34, 1850.41, 1867.40, 1885.24,1903.85,1923.11,
!     &           1943.09, 1963.67, 1974.17/

!    Al2O3
!       DATA   TCONAl2O3 /1685.09, 1693.08, 1700.92,1708.84, 1716.78,
!     &              1724.79,1732.91, 1741.42,1750.37,1759.75,
!     &              1769.31,1779.06,1789.37,1799.94, 1810.65,1820.94,
!!     &              1830.99,1841.08,1851.02, 1861.02,1871.26,1881.74,
 !    &              1892.36,1903.03,1913.69,1924.39,1935.05,1945.95,
 !!!    &              1957.45,1969.45,1980.72,1989.12,1995.54, 2002.05,
!     &              2010.76,2021.69,2033.01,2044.36, 2055.67,2066.99,
!     &              2078.33,2089.65,2100.99,2112.62,2124.88,2137.43,
!     &              2150.32,2163.28, 2176.89, 2191.32,2198.76/
!    Fe
!      DATA   TCONFe /1362.14,1373.43, 1384.60,1395.83, 1407.00,1418.26,
!     &              1430.00,1442.48, 1455.61, 1469.10, 1482.58,1496.08,
!     &              1509.58,1523.08,1536.57,1550.07,1563.57,1577.13,
!     &              1590.57,1604.74,1619.69, 1635.41, 1651.59, 1667.75,
!!     &             1684.03,1700.80, 1718.31, 1736.60, 1755.27, 1773.98,
 !    &             1792.93, 1812.32,1832.10, 1852.28,1872.54, 1892.90,
 !    &             1913.24, 1934.27, 1956.41,1978.37, 2008.05, 2030.80,
 !    &             2051.13, 2081.89, 2103.84, 2132.13,2157.98, 2190.91,
 !    &              2221.92, 2247.77, 2263.48/
      !  READ (7,INCLOUDY)

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
          DO IL=PINDEXTOPC,PINDEXBOTC-1
           ABOVE        = MAX(IL-1,1)
           BELOW        = MIN(IL+1,NLEV)
       VERTPROF(IL) = (VERTPROF(ABOVE)+VERTPROF(IL)+VERTPROF(BELOW))/3.
           TOTAL        = TOTAL+VERTPROF(IL)
          ENDDO
           TOTAL        = TOTAL+VERTPROF(PINDEXBOTC)
!         Normalize it to preserve the total integrated optical.
          DO IL=PINDEXTOPC,PINDEXBOTC
           VERTPROF(IL)=VERTPROF(IL)/TOTAL
          ENDDO

!         Now we have a dimensionless array that integrates to one and defines
!         the shape of the cloud in the vertical dimension.
!         For future use, let's evaluate the maximum optical depth of
!         any given layer and the location of that cloud depth.
          MAXTAU=MAXVAL(VERTPROF,1)
          MAXTAULOC=MAXLOC(VERTPROF,1)


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
       WRITE(61,*) 'AEROSOL: ',AEROSOLCOMP
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
       WRITE(61,*) '   PRESSURE(BAR)        NORMALIZED AEROSOL TAU'
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
             IF (TheTAU*EXTFACTLW.LT.1e-7) THEN
             XFACTLW = 0.0
             ENDIF

             WRITE(61,212) PRESSURE(ILEV)*1E-5,TheTAU,PI0AERSW*XFACTSW,
     &                     ASYMSW*XFACTSW,TheTAU*EXTFACTLW,
     &                     PI0AERLW*XFACTLW,ASYMLW*XFACTLW
 212       FORMAT(2X,F11.6,3X,F14.5,3X,F7.4,3X,F7.4,3X,F11.6
     &             ,3X,F7.4,3X,F7.4)
                ENDDO
             WRITE(61,*)''
             ENDDO
          ENDDO
       ENDDO
!

      write(*,*) 'In rmakeclouds, aerosolcomp =',aerosolcomp

      IF ((AEROSOLCOMP.EQ. 'standard')) THEN
          DO J = 1, NL+1
              TCON(J) = 0.
          END DO
      ELSE
          WRITE(*,*)'NO VALID AEROSOLCOMP SPECIFIED! STOPPING'
          STOP
      END IF

       END
