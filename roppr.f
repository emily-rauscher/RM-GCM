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
      REAL CONDFACT(NL+1),CLOUDLOC(NL+1)
      INTEGER BASELEV,TOPLEV
!     MnS
!      DATA  TCON /1113.43, 1120.07, 1126.64, 1133.21,1139.9,1146.46,
!     &            1153.15,  1159.65, 1166.25, 1173.01, 1180.09, 1187.73,
!     &            1194.69, 1202.12, 1209.26, 1216.50, 1223.92, 1231.43,
!     &            1239.19, 1247.17, 1255.39,1263.78, 1272.32,1281.40,
!     &            1289.69, 1297.76, 1306.00, 1315.24, 1324.18, 1333.27,
!     &            1342.44, 1351.39, 1360.54, 1369.99, 1379.72, 1389.42,
!     &            1399.22, 1409.04,  1418.99,1428.77, 1438.60, 1449.11,
!     &            1459.19, 1469.78, 1481.06, 1492.70, 1504.21, 1515.49,
!     &            1527.84,1540.17,1545.90/

!    SiO2
!      DATA   TCON /1334.63,1342.58,1350.30,1358.48,1366.64,1374.85,
!     &             1383.15,1391.59,1400.10,1408.68, 1417.25,1425.87,
!     &             1434.53, 1443.14, 1451.71,1460.28,1468.90,1477.44,
!     &             1486.12,1494.77,1503.91,1513.72 ,1524.07, 1534.57,
!     &             1544.98,1555.25, 1565.35,1575.40, 1585.44,1595.55,
!     &             1606.04,1617.06, 1628.55,1640.25,1651.98, 1664.07,
!     &             1676.82,1689.97,1703.53, 1717.16, 1731.36, 1746.01,
!     &             1761.10,1776.94, 1793.70,1811.36,1829.60, 1848.37,
!     &             1867.54, 1887.00, 1896.46/

!    Mg2SiO4
!      DATA   TCON /1370.00, 1378.31, 1386.62, 1395.03, 1403.56,1412.29,
!     &              1421.17, 1430.18, 1439.17, 1448.16, 1457.24,1466.52,
!     &              1475.94,1485.57, 1495.23, 1505.09, 1515.04, 1525.21,
!     &              1535.33,1545.80, 1556.37, 1567.12,1578.02, 1589.13,
!     &              1600.29, 1611.67,1623.20, 1634.84,1646.61,1658.58,
!     &              1670.74, 1683.03, 1695.59, 1708.44,1721.29,1734.41,
!     &              1747.86, 1761.64, 1775.79, 1789.95, 1804.36,1819.11,
!     &              1834.34, 1850.41, 1867.40, 1885.24,1903.85,1923.11,
!     &              1943.09, 1963.67, 1974.17/
      
!    Al2O3
!      DATA   TCON /1677.58,1677.5,1677.58,1677.58,1677.58, 1677.58,
!     &               1677.58,1677.5, 1677.58,1677.58, 1677.58,1677.59,
!     &              1677.59,1677.59, 1677.60,1677.60, 1677.61,1677.63,
!     &              1677.64,1677.66, 1677.69,1677.73,1677.78,1677.84,
!     &              1677.93,1678.04,1678.19,1678.39, 1678.64,1678.98,
!     &              1679.43,1680.02, 1680.80,1681.82, 1683.16,1684.92,
!     &              1687.25 ,1690.31, 1694.34,1699.66, 1706.67,1715.99,
!     &              1728.44,1746.08,1771.48,1808.33, 1855.14,1919.59,
!     &              1999.09, 2115.80,2202.88/

!    Fe
!      DATA   TCON /1358.12,1358.12, 1358.12,1358.12,1358.12,1358.12,
!     &              1358.12, 1358.12,1358.13, 1358.13, 1358.13,1358.13,
!     &              1358.14, 1358.14 ,1358.15, 1358.16, 1358.17,1358.19,
!     &              1358.21, 1358.24, 1358.28, 1358.33, 1358.40,1358.49,
!     &              1358.61, 1358.76, 1358.97, 1359.24, 1359.59,1360.06,
!     &              1360.67, 1361.48, 1362.55,1363.96, 1365.83,1368.27,
!     &              1371.5, 1375.78, 1381.43,1388.85, 1398.61,1411.46,
!     &              1429.01, 1454.56, 1490.12, 1537.02, 1599.25,1695.24,
!     &              1845.34, 2084.92, 2270.26/

   
!  110   CONTINUE
      DO 199  J          =   1,NLAYER
       CONDFACT(J)       = min(max(TCON(J)-TT(J),0.0),1.0)
       CLOUDLOC(J)       = J*CONDFACT(J)
199   CONTINUE
      IF (SUM(CONDFACT).gt.0) THEN
          BASELEV         = MAXVAL(CLOUDLOC,1)
          TOPLEV          =MAX(BASELEV-AERLAYERS,1) 
      DO 200  J          = BASELEV,TOPLEV,-1  
       TAUAER(1,J)       = AEROPROF(J)*CONDFACT(J)
       WOL(1,J)          = PI0AERSW 
       GOL(1,J)          = ASYMSW
!       TAUS(1,J)         = PI0AERSW*AEROPROF(J)
200   CONTINUE 
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
!       TAUS(2,JJ)         = PI0AERLW*TAUAER(2,K)
              JJ          = J+1
       TAUAER(2,JJ)       = TAUAER(1,K)*EXTFACTLW
       WOL(2,JJ)          = PI0AERLW
       GOL(2,JJ)          = ASYMLW
!       TAUS(2,JJ)         = PI0AERLW*TAUAER(2,K)
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

              TAUL(L,J)   = TAUGAS(L,J)+TAURAY(L,J)+TAUAER(L,J)!+TAUCLD(L,J)
             
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

