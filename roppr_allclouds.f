      SUBROUTINE OPPR
!
!     **************************************************************
!     *  Purpose             :  CaLculates optical properties      *
!     *                         such as single scattering albedo,  *
!     *                         asymmetry parameter, etc.          *
!     *                         This routine is case dependent and *
!     *                         wiLL have to be repLaced by the    *
!     *                         user.                              *
!     *  Subroutines Called  :  None                               *
!     *  Input               :  PAH2O, RDH2O, CO2, O3, ETC         *
!     *  Output              :  TAUL, W0, G0, OPD, Y3              *
!     * ************************************************************
!
      include 'rcommons.h'

!     W0(NWAVE,NLAYER) : SINGLE SCATTERING ALBEDO *** delta scaled ***
!     G0(NWAVE,NLAYER) : ASYMMETRY PARAMETER *** delta scaled ***
!     OPD(NWAVE,NLAYER): cumulative OPTICAL DEPTH *** delta scaled ***
!     SFL(NWAVE)       : SOLAR FLUX
!    uW0(NWAVE,NLAYER)  : unscaled SINGLE SCATTERING ALBEDO 
!    uG0(NWAVE,NLAYER)  : unscaled ASYMMETRY PARAMETER 
!    uTAUL(NWAVE,NLAYER): unscaled OPTICAL DEPTH of layer
!
!     ASSUME THAT P IS SAME ON ALL SIGMA LEVELS. IF PSURFACE
!     VARIES A LOT, THEN WE MAY NEED TO CALCULATE TAUO3,
!     TAUCO2, TAUO2 FOR EACH POINT.
!
!     NOTE : THE TOP LAYER IS A DUMMY. IT CONTAINS A DEFINED GAS
!            AMOUNT. DIFFERENT MODELS WILL REQUIRE DIFFERENT
!            TREATMENT OF THIS.
!     CALCULATE TOTAL OPTICAL DEPTH INCLUDING GASES. THEN
!     GIVEN THE AEROSOL OPTICAL DEPTHS AND CLOUD OPTICAL DEPTHS,
!     CALCULATE FINAL OPTICAL PROPERTIES. WE USE A DELTA
!     TWO STREAM APPROACH TO FIND W0, SINGLE SCATTERING ALBEDO,
!     G0, ASYMMMETRY PARAMETER, TAUL, LAYER OPTICAL DEPTH,
!     OPD, CUMULATIVE OPTICAL DEPTH TO BASE OF LAYER.
!     open( 23, file='tau.dat', form='formatted', status='unknown' )   
!     write(23,*) nlayer, ntotal
!
!   NOTE: THIS IS DOUBLE GRAY SPECIFIC; WOULD REQUIRE GENERALIZATION

!   Simple cloud test:
      REAL TCONSiO2(NL+1),CONDFACTSiO2(NL+1),CLOUDLOCSiO2(NL+1)
      REAL TCONMnS(NL+1),CONDFACTMnS(NL+1),CLOUDLOCMnS(NL+1)
      REAL TCONMg2SiO4(NL+1),CONDFACTMg2SiO4(NL+1)
      REAL CLOUDLOCMg2SiO4(NL+1)
      REAL TCONAl2O3(NL+1),CONDFACTAl2O3(NL+1),CLOUDLOCAl2O3(NL+1)
      REAL TCONFe(NL+1),CONDFACTFe(NL+1),CLOUDLOCFe(NL+1)
      REAL TCONMgSiO3(NL+1),CONDFACTMgSiO3(NL+1),CLOUDLOCMgSiO3(NL+1)
      INTEGER BASELEV,TOPLEV

!     MnS
      DATA TCONMnS /1113.43, 1120.07,1126.64,1133.21,1139.9,1146.46,
     &        1153.15,  1159.65, 1166.25, 1173.01, 1180.09, 1187.73,
     &        1194.69, 1202.12, 1209.26, 1216.50, 1223.92, 1231.43,
     &       1239.19, 1247.17, 1255.39,1263.78, 1272.32,1281.40,
     &        1289.69, 1297.76, 1306.00, 1315.24, 1324.18, 1333.27,
     &        1342.44, 1351.39, 1360.54, 1369.99, 1379.72, 1389.42,
     &        1399.22, 1409.04,  1418.99,1428.77, 1438.60, 1449.11,
     &        1459.19, 1469.78, 1481.06, 1492.70, 1504.21, 1515.49,
     &            1527.84,1540.17,1545.90/

!    SiO2
      DATA TCONSiO2 /1334.63,1342.58,1350.30,1358.48,1366.64,1374.85,
     &             1383.15,1391.59,1400.10,1408.68, 1417.25,1425.87,
     &             1434.53, 1443.14, 1451.71,1460.28,1468.90,1477.44,
     &             1486.12,1494.77,1503.91,1513.72 ,1524.07, 1534.57,
     &             1544.98,1555.25, 1565.35,1575.40, 1585.44,1595.55,
     &             1606.04,1617.06, 1628.55,1640.25,1651.98, 1664.07,
     &             1676.82,1689.97,1703.53, 1717.16, 1731.36, 1746.01,
     &             1761.10,1776.94, 1793.70,1811.36,1829.60, 1848.37,
     &             1867.54, 1887.00, 1896.46/

!    MgSiO3
      DATA   TCONMgSiO3 /1307.11,1314.43,1321.68,1329.11,1336.79, 
     &             1344.72,1352.04,1359.52,1367.46,1375.60,1383.27,
     &             1391.60,1400.24,1408.51,1416.97,1425.28,1434.14,
     &             1442.49,1451.31,1460.37,1469.43,1478.40,1487.67, 
     &             1497.14,1507.25,1517.04,1527.28,1537.29,1546.58,
     &             1556.41,1567.50,1578.63,1588.50,1598.45,1609.41,
     &             1620.62,1631.41,1642.37,1654.32,1666.68,1678.34,
     &             1689.61,1702.01,1715.51,1727.89,1739.21,1751.75, 
     &             1765.34, 1778.58,1790.56,1796.07/     

!    Mg2SiO4
      DATA   TCONMg2SiO4 /1370.00,1378.31,1386.62,1395.03,1403.56,
     &           1412.29,
     &           1421.17, 1430.18, 1439.17, 1448.16, 1457.24,1466.52,
     &           1475.94,1485.57, 1495.23, 1505.09, 1515.04, 1525.21,
     &           1535.33,1545.80, 1556.37, 1567.12,1578.02, 1589.13,
     &           1600.29, 1611.67,1623.20, 1634.84,1646.61,1658.58,
     &           1670.74, 1683.03, 1695.59, 1708.44,1721.29,1734.41,
     &           1747.86, 1761.64, 1775.79, 1789.95, 1804.36,1819.11,
     &           1834.34, 1850.41, 1867.40, 1885.24,1903.85,1923.11,
     &           1943.09, 1963.67, 1974.17/
        
!    Al2O3
       DATA   TCONAl2O3 /1685.09, 1693.08, 1700.92,1708.84, 1716.78,
     &              1724.79,1732.91, 1741.42,1750.37,1759.75, 
     &              1769.31,1779.06,1789.37,1799.94, 1810.65,1820.94,
     &              1830.99,1841.08,1851.02, 1861.02,1871.26,1881.74,
     &              1892.36,1903.03,1913.69,1924.39,1935.05,1945.95,
     &              1957.45,1969.45,1980.72,1989.12,1995.54, 2002.05,
     &              2010.76,2021.69,2033.01,2044.36, 2055.67,2066.99,
     &              2078.33,2089.65,2100.99,2112.62,2124.88,2137.43,
     &              2150.32,2163.28, 2176.89, 2191.32,2198.76/
!    Fe
      DATA   TCONFe /1362.14,1373.43, 1384.60,1395.83, 1407.00,1418.26, 
     &              1430.00,1442.48, 1455.61, 1469.10, 1482.58,1496.08,
     &              1509.58,1523.08,1536.57,1550.07,1563.57,1577.13,
     &              1590.57,1604.74,1619.69, 1635.41, 1651.59, 1667.75,
     &             1684.03,1700.80, 1718.31, 1736.60, 1755.27, 1773.98,
     &             1792.93, 1812.32,1832.10, 1852.28,1872.54, 1892.90,
     &             1913.24, 1934.27, 1956.41,1978.37, 2008.05, 2030.80,
     &             2051.13, 2081.89, 2103.84, 2132.13,2157.98, 2190.91,
     &              2221.92, 2247.77, 2263.48/
  
!  110   CONTINUE
      DO 180  J          =   1,NLAYER
       CONDFACTMgSiO3(J)       = min(max(TCONMgSiO3(J)-TT(J),0.0),1.0)
       CLOUDLOCMgSiO3(J)       = J*CONDFACTMgSiO3(J)
       CONDFACTMnS(J)       = min(max(TCONMnS(J)-TT(J),0.0),1.0)
       CLOUDLOCMnS(J)       = J*CONDFACTMnS(J)      
       CONDFACTMg2SiO4(J)     = min(max(TCONMg2SiO4(J)-TT(J),0.0),1.0)
       CLOUDLOCMg2SiO4(J)       = J*CONDFACTMg2SiO4(J) 
       CONDFACTAl2O3(J)       = min(max(TCONAl2O3(J)-TT(J),0.0),1.0)
       CLOUDLOCAl2O3(J)       = J*CONDFACTAl2O3(J)
       CONDFACTFe(J)       = min(max(TCONFe(J)-TT(J),0.0),1.0)
       CLOUDLOCFe(J)       = J*CONDFACTFe(J)        
180   CONTINUE
! Now go through each, cloud by cloud

      IF (SUM(CONDFACTMnS).gt.0) THEN
          BASELEV         = MAXVAL(CLOUDLOCMnS,1)
          TOPLEV          = max(BASELEV-5,1)
      DO 181  J          = BASELEV,TOPLEV,-1
       TAUAER(1,J)       = AEROPROF(J)*CONDFACTMnS(J)*.008
       WOL(1,J)          = .9999
       GOL(1,J)          = .23
!       TAUS(1,J)         = PI0AERSW*AEROPROF(J)*CONDFACTMnS(J)
181   CONTINUE
       END IF

      IF (SUM(CONDFACTAl2O3).gt.0) THEN
          BASELEV         = MAXVAL(CLOUDLOCAl2O3,1)
          TOPLEV          = max(BASELEV-5,1)
      DO 182  J          = BASELEV,TOPLEV,-1
       TAUAER(1,J)     = AEROPROF(J)*CONDFACTAl2O3(J)*.043+TAUAER(1,J)
       WOL(1,J)          = 0.87
       GOL(1,J)          = 0.66
!       TAUS(1,J)         = PI0AERSW*AEROPROF(J)*CONDFACTAl2O3(J)
182   CONTINUE
       END IF

      IF (SUM(CONDFACTFe).gt.0) THEN
          BASELEV         = MAXVAL(CLOUDLOCFe,1)
          TOPLEV          = max(BASELEV-5,1)
      DO 183  J          = BASELEV,TOPLEV,-1
       TAUAER(1,J)       = AEROPROF(J)*CONDFACTFe(J)*.47+TAUAER(1,J)
       WOL(1,J)          = 0.69
       GOL(1,J)          = 0.42
!       TAUS(1,J)         = PI0AERSW*AEROPROF(J)*CONDFACTFe(J)
183   CONTINUE
       END IF

      IF (SUM(CONDFACTMg2SiO4).gt.0) THEN
          BASELEV         = MAXVAL(CLOUDLOCMg2SiO4,1)
          TOPLEV          = max(BASELEV-5,1)
      DO 184  J          = BASELEV,TOPLEV,-1
       TAUAER(1,J)     = AEROPROF(J)*CONDFACTMg2SiO4(J)*.72+TAUAER(1,J)
       WOL(1,J)          = .99999
       GOL(1,J)          = 0.6
!       TAUS(1,J)         = PI0AERSW*AEROPROF(J)*CONDFACTMg2SiO4(J)
184   CONTINUE
       END IF


      IF (SUM(CONDFACTMgSiO3).gt.0) THEN
          BASELEV         = MAXVAL(CLOUDLOCMgSiO3,1)
          TOPLEV          = max(BASELEV-5,1) 
      DO 185  J          = BASELEV,TOPLEV,-1  
       TAUAER(1,J)       = AEROPROF(J)*CONDFACTMgSiO3(J)+TAUAER(1,J)
       WOL(1,J)          = 0.99999 
       GOL(1,J)          = 0.6
!       TAUS(1,J)         = PI0AERSW*AEROPROF(J)*CONDFACTMgSiO3(J)
185   CONTINUE 
       END IF

      DO 201  J          = NLAYER+1,NDBL
!       DO 130 L          = 1,NTOTAL
!    VIS
       TAUAER(1,J)       = 0.
       WOL(1,J)          = 0. 
       GOL(1,J)          = 0. 
       TAUS(1,J)         = 0. 
201   CONTINUE
      
!           IF (TAUS(2,J) .NE. 0) THEN 
!            TTAS         = TAUS(1,J)
!           ELSE
!            TTAS         = 1.
!           ENDIF

!       GOL(1,J)          = G01(1,J)/TTAS

!    IR
               k         =  1
      DO 202  J          =   1,NDBL,2
              JJ         = J     
       TAUAER(2,JJ)       = TAUAER(1,K)*EXTFACTLW
       WOL(2,JJ)          = PI0AERLW
       GOL(2,JJ)          = ASYMLW
       TAUS(2,JJ)         = PI0AERLW*TAUAER(2,K)
              JJ          = J+1
       TAUAER(2,JJ)       = TAUAER(1,K)*EXTFACTLW
       WOL(2,JJ)          = PI0AERLW
       GOL(2,JJ)          = ASYMLW
       TAUS(2,JJ)         = PI0AERLW*TAUAER(2,K)
             k            = k+1
202   continue

!       DO J = 1,NDBL
!       write(*,*),TAUAER(1,J),TAUAER(2,J)
!       enddo
!       STOP 
!         if (TAUS(2,J).ne.0) then
!            TTAS         = TAUS(2,J)
!         else
!            TTAS         = 1.
!         ENDIF

!       GOL(2,J)          = G01(1,J)/TTAS



!            do j  = 1,nlayer
!            write(*,*)'j,Tauray/Taugas',j,TAURAY(1,J)/TAUGAS(1,J)
!            enddo

!     iradgas = 0: no gas in radiative xfer!
      iradgas = 1     

      DO 500 J           = 1,NLAYER
          j1             = max( 1, j-1 )
!
!     First the solar at standard resolution
          DO 400 L       = LLS,NSOLP
!M          DO 400 L       = LLS,LLA
 
!MTR              TAUL(L,J)   = TAUH2O(L,J)+TAUGAS(L,J)+   &
!MTR                               PARAY(L,J)+TAUAER(L,J)+TAUCLD(L,J)
!              THE CODE ORIGINALLY INCLUDED TAUCLD AS WELL AS TAUAER. I
!              DO NOT SEE THE POINT OF HAVING BOTH OPTICAL DEPTHS IN OUR
!              MODEL SO I AM OMMITTING TAUCLD, W0CLD, GCLD, ETC

          TAUL(L,J) = TAUGAS(L,J)+TAURAY(L,J)+TAUAER(L,J)!+TAUCLD(L,J)
             
             if (iradgas.eq.0) then
             tauL(L,j) = tauaer(L,j)
             endif

             if( TAUL(L,J) .lt. EPSILON ) then
             TAUL(L,J) = EPSILON
             endif

             utauL(L,j)  = TAUL(L,J)
          WOT         = (TAURAY(L,J)+TAUAER(L,J)*WOL(L,J))/TAUL(L,J)
             if (iradgas.eq.0) then
              wot = woL(L,j)
             endif
 
             WOT         = min(1.-EPSILON,WOT)
             uw0(L,j)    = WOT
!             write(*,*) 'WOT',WOT
             DENOM       = (TAURAY(L,J)+ TAUAER(L,J)*WOL(L,J))
             if( DENOM .LE. EPSILON ) then 
             DENOM = EPSILON then
             endif
             if( DENOM .GT. EPSILON ) then
               GOT = ( GOL(L,J)* WOL(L,J)*TAUAER(L,J) ) / DENOM
!       print*, j, L, GCLD(l,j), gol(l,j),wol(l,j),got
             else
               GOT = 0.
             endif
             if (iradgas.eq.0) then
             GOT = goL(L,j)
             endif
             ug0(L,j)    = GOT
             uOPD(L,J)   = 0.0
             uOPD(L,J)   = uOPD(L,J1)+uTAUL(L,J)
             IF (deltascale) THEN
             FO          = GOT*GOT
             DEN         = 1.-WOT*FO
             TAUL(L,J)   = TAUL(L,J) * DEN
             W0(L,J)     = (1.-FO)*WOT/DEN
             G0(L,J)     = GOT/(1.+GOT)
             OPD(L,J)    = 0.0
             OPD(L,J)    = OPD(L,J1)+TAUL(L,J)             
             ELSE 
                  W0(L,J)= uw0(L,J)
                  G0(L,J)= ug0(L,J)
                TAUL(L,J)= utaul(L,J)
                 OPD(L,J)= uOPD(L,J)
             ENDIF 
!!!!!!!!!!!!!!!!HERE'S WHERE YOU CAN HARDWIRE VALUES!!!!!!!!!
 
             if( taul(L,j).lt.0. ) then
       write(*,*) 'taul lt 0'
       stop
       endif
                       
 400        CONTINUE
           
!MTR          IF(IR .EQ. 1) THEN
!MTR             DO 450 I        =   1,NGAUSS
!MTR                DO 425 L     =   LLS,LLA
!MTR                   Y3(L,I,J) =   EXP(-TAUL(L,J)/GANGLE(I))
!MTR 425            CONTINUE
!MTR 450         CONTINUE
!          ENDIF
 500  CONTINUE

!      NOW AGAIN FOR THE IR
      DO 501 J           = 1,NDBL
          j1             = max( 1, j-1 )
!
!     First the solar at standard resolution
          DO 401 L       = NSOLP+1,NTOTAL

         TAUL(L,J)  =  TAUGAS(L,J)+TAURAY(L,J)+TAUAER(L,J)!+TAUCLD(L,J)
             if (iradgas.eq.0) then
             tauL(L,j) = tauaer(L,j)
             endif

             if( TAUL(L,J) .lt. EPSILON ) then
             TAUL(L,J) = EPSILON
             endif

             utauL(L,j)  = TAUL(L,J)
          WOT         = (TAURAY(L,J)+TAUAER(L,J)*WOL(L,J))/TAUL(L,J)
             if (iradgas.eq.0) then
              wot = woL(L,j)
             endif

             WOT         = min(1.-EPSILON,WOT)
             uw0(L,j)    = WOT
             DENOM       = (TAURAY(L,J)+ TAUAER(L,J)*WOL(L,J))
             if( DENOM .LE. EPSILON ) then
             DENOM = EPSILON then
             endif
             if( DENOM .GT. EPSILON ) then
               GOT = ( GOL(L,J)* WOL(L,J)*TAUAER(L,J) ) / DENOM        
             else
               GOT = 0.
             endif
             if (iradgas.eq.0) then
             GOT = goL(L,j)
             endif
             ug0(L,j)    = GOT
             uOPD(L,J)   = 0.0
             uOPD(L,J)   = uOPD(L,J1)+uTAUL(L,J)
             IF (deltascale) THEN
             FO          = GOT*GOT
             DEN         = 1.-WOT*FO
             TAUL(L,J)   = TAUL(L,J) * DEN
             W0(L,J)     = (1.-FO)*WOT/DEN
             G0(L,J)     = GOT/(1.+GOT)
             OPD(L,J)    = 0.0
             OPD(L,J)    = OPD(L,J1)+TAUL(L,J)
             ELSE
                  W0(L,J)= uw0(L,J)
                  G0(L,J)= ug0(L,J)
                TAUL(L,J)= utaul(L,J)
                 OPD(L,J)= uOPD(L,J)
             ENDIF

             if( taul(L,j).lt.0. ) then
       write(*,*) 'taul lt 0'
       stop
       endif

 401        CONTINUE

          IF(IR .EQ. 1) THEN
             DO 450 I        =   1,NGAUSS
                DO 425 L     =   LLS,LLA
                   Y3(L,I,J) =   EXP(-TAUL(L,J)/GANGLE(I))
 425            CONTINUE
 450         CONTINUE
          ENDIF
 501  CONTINUE

!      DO J = 1,NDBL 
!      write(*,*) 'W0',J,W0(2,J)
!      ENDDO 
!      DO J = 1,NDBL
!      write(*,*) 'G0',J,G0(2,J)
!      ENDDO
!      DO J = 1,NDBL
!      write(*,*) 'WOL',J,WOL(2,J)
!      ENDDO
!      DO J = 1,NDBL
!      write(*,*) 'GOL',J,GOL(2,J)
!      ENDDO
!      DO J = 1,NDBL
!      write(*,*) 'TAURAY',J,TAURAY(2,J)
!      ENDDO
!      DO J = 1,NDBL
!      write(*,*) 'TAUAER',J,TAUAER(2,J)
!      ENDDO
!      DO J = 1,NDBL
!      write(*,*) 'TAUGAS',J,TAUGAS(2,J)
!      ENDDO
!      DO J = 1,NDBL
!      write(*,*) 'TAUL',J,TAUL(2,J)
!      ENDDO
!      DO J = 1,NDBL
!      write(*,*) 'OPD',J,OPD(2,J) 
!      ENDDO     
!         stop
      RETURN
      END

