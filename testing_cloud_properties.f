      SUBROUTINE get_cloud_scattering_properties(NL)
          implicit none
          integer :: J, L, K, NL

          ! Define all the arrays
          ! These are 50 by 50 because that's what the data in CLOUD_DATA is
          ! That can change but use to code from Elsie and Isaac
          real, dimension(50,50) :: Al2O3_rosselandMean_gg
          real, dimension(50,50) :: Al2O3_rosselandMean_qext
          real, dimension(50,50) :: Al2O3_rosselandMean_qscat
          real, dimension(50,50) :: Al2O3_vis_500_gg
          real, dimension(50,50) :: Al2O3_vis_500_qext
          real, dimension(50,50) :: Al2O3_vis_500_qscat
          real, dimension(50,50) :: Al2O3_vis_650_gg
          real, dimension(50,50) :: Al2O3_vis_650_qext
          real, dimension(50,50) :: Al2O3_vis_650_qscat
          real, dimension(50,50) :: Al2O3_vis_800_gg
          real, dimension(50,50) :: Al2O3_vis_800_qext
          real, dimension(50,50) :: Al2O3_vis_800_qscat
          real, dimension(50,50) :: Fe_rosselandMean_gg
          real, dimension(50,50) :: Fe_rosselandMean_qext
          real, dimension(50,50) :: Fe_rosselandMean_qscat
          real, dimension(50,50) :: Fe_vis_500_gg
          real, dimension(50,50) :: Fe_vis_500_qext
          real, dimension(50,50) :: Fe_vis_500_qscat
          real, dimension(50,50) :: Fe_vis_650_gg
          real, dimension(50,50) :: Fe_vis_650_qext
          real, dimension(50,50) :: Fe_vis_650_qscat
          real, dimension(50,50) :: Fe_vis_800_gg
          real, dimension(50,50) :: Fe_vis_800_qext
          real, dimension(50,50) :: Fe_vis_800_qscat
          real, dimension(50,50) :: Mg2SiO4_rosselandMean_gg
          real, dimension(50,50) :: Mg2SiO4_rosselandMean_qext
          real, dimension(50,50) :: Mg2SiO4_rosselandMean_qscat
          real, dimension(50,50) :: Mg2SiO4_vis_500_gg
          real, dimension(50,50) :: Mg2SiO4_vis_500_qext
          real, dimension(50,50) :: Mg2SiO4_vis_500_qscat
          real, dimension(50,50) :: Mg2SiO4_vis_650_gg
          real, dimension(50,50) :: Mg2SiO4_vis_650_qext
          real, dimension(50,50) :: Mg2SiO4_vis_650_qscat
          real, dimension(50,50) :: Mg2SiO4_vis_800_gg
          real, dimension(50,50) :: Mg2SiO4_vis_800_qext
          real, dimension(50,50) :: Mg2SiO4_vis_800_qscat
          real, dimension(50,50) :: MnS_rosselandMean_gg
          real, dimension(50,50) :: MnS_rosselandMean_qext
          real, dimension(50,50) :: MnS_rosselandMean_qscat
          real, dimension(50,50) :: MnS_vis_500_gg
          real, dimension(50,50) :: MnS_vis_500_qext
          real, dimension(50,50) :: MnS_vis_500_qscat
          real, dimension(50,50) :: MnS_vis_650_gg
          real, dimension(50,50) :: MnS_vis_650_qext
          real, dimension(50,50) :: MnS_vis_650_qscat
          real, dimension(50,50) :: MnS_vis_800_gg
          real, dimension(50,50) :: MnS_vis_800_qext
          real, dimension(50,50) :: MnS_vis_800_qscat


          REAL TCONMnS(51)
          REAL TCONMg2SiO4(51)
          REAL TCONFe(51)
          REAL TCONAL2O3(51)

          REAL TCONDS(51,4)

          REAL QE_OPPR(5, 50, 50, 4)
          REAL PI0_OPPR(5, 50, 50, 4)
          REAL G0_OPPR(5, 50, 50, 4)

          REAL DENSITY(4)
          REAL FMOLW(4)
          REAL CORFACT(51)

          real, dimension(50) :: input_particle_size_array_in_meters
          real, dimension(50) :: input_temperature_array
          real, dimension(50) :: particle_size_vs_layer_array_in_meters

          COMMON /CLOUD_PROPERTIES/ TCONDS, QE_OPPR, PI0_OPPR, G0_OPPR,
     &                              DENSITY, FMOLW, CORFACT,
     &                              input_particle_size_array_in_meters,
     &                              input_temperature_array,
     &                              particle_size_vs_layer_array_in_meters

          ! opening the file for reading
          open (1, file='../CLOUD_DATA/Al2O3_rosselandMean_gg.txt')
          read(1,*) Al2O3_rosselandMean_gg
          close(1)
          open (2, file='../CLOUD_DATA/Al2O3_rosselandMean_qext.txt')
          read(2,*) Al2O3_rosselandMean_qext
          close(2)
          open (3, file='../CLOUD_DATA/Al2O3_rosselandMean_qscat.txt')
          read(3,*) Al2O3_rosselandMean_qscat
          close(3)
          open (4, file='../CLOUD_DATA/Al2O3_vis_500_gg.txt')
          read(4,*) Al2O3_vis_500_gg
          close(4)
          open (5, file='../CLOUD_DATA/Al2O3_vis_500_qext.txt')
          read(5,*) Al2O3_vis_500_qext
          close(5)
          open (6, file='../CLOUD_DATA/Al2O3_vis_500_qscat.txt')
          read(6,*) Al2O3_vis_500_qscat
          close(6)
          open (7, file='../CLOUD_DATA/Al2O3_vis_650_gg.txt')
          read(7,*) Al2O3_vis_650_gg
          close(7)
          open (8, file='../CLOUD_DATA/Al2O3_vis_650_qext.txt')
          read(8,*) Al2O3_vis_650_qext
          close(8)
          open (9, file='../CLOUD_DATA/Al2O3_vis_650_qscat.txt')
          read(9,*) Al2O3_vis_650_qscat
          close(9)
          open (10, file='../CLOUD_DATA/Al2O3_vis_800_gg.txt')
          read(10,*) Al2O3_vis_800_gg
          close(10)
          open (11, file='../CLOUD_DATA/Al2O3_vis_800_qext.txt')
          read(11,*) Al2O3_vis_800_qext
          close(11)
          open (12, file='../CLOUD_DATA/Al2O3_vis_800_qscat.txt')
          read(12,*) Al2O3_vis_800_qscat
          close(12)
          open (13, file='../CLOUD_DATA/Fe_rosselandMean_gg.txt')
          read(13,*) Fe_rosselandMean_gg
          close(13)
          open (14, file='../CLOUD_DATA/Fe_rosselandMean_qext.txt')
          read(14,*) Fe_rosselandMean_qext
          close(14)
          open (15, file='../CLOUD_DATA/Fe_rosselandMean_qscat.txt')
          read(15,*) Fe_rosselandMean_qscat
          close(15)
          open (16, file='../CLOUD_DATA/Fe_vis_500_gg.txt')
          read(16,*) Fe_vis_500_gg
          close(16)
          open (17, file='../CLOUD_DATA/Fe_vis_500_qext.txt')
          read(17,*) Fe_vis_500_qext
          close(17)
          open (18, file='../CLOUD_DATA/Fe_vis_500_qscat.txt')
          read(18,*) Fe_vis_500_qscat
          close(18)
          open (19, file='../CLOUD_DATA/Fe_vis_650_gg.txt')
          read(19,*) Fe_vis_650_gg
          close(19)
          open (20, file='../CLOUD_DATA/Fe_vis_650_qext.txt')
          read(20,*) Fe_vis_650_qext
          close(20)
          open (21, file='../CLOUD_DATA/Fe_vis_650_qscat.txt')
          read(21,*) Fe_vis_650_qscat
          close(21)
          open (22, file='../CLOUD_DATA/Fe_vis_800_gg.txt')
          read(22,*) Fe_vis_800_gg
          close(22)
          open (23, file='../CLOUD_DATA/Fe_vis_800_qext.txt')
          read(23,*) Fe_vis_800_qext
          close(23)
          open (24, file='../CLOUD_DATA/Fe_vis_800_qscat.txt')
          read(24,*) Fe_vis_800_qscat
          close(24)
          open (25, file='../CLOUD_DATA/Mg2SiO4_rosselandMean_gg.txt')
          read(25,*) Mg2SiO4_rosselandMean_gg
          close(25)
          open (26, file='../CLOUD_DATA/Mg2SiO4_rosselandMean_qext.txt')
          read(26,*) Mg2SiO4_rosselandMean_qext
          close(26)
          open (27, file='../CLOUD_DATA/Mg2SiO4_rosselandMean_qscat.txt')
          read(27,*) Mg2SiO4_rosselandMean_qscat
          close(27)
          open (28, file='../CLOUD_DATA/Mg2SiO4_vis_500_gg.txt')
          read(28,*) Mg2SiO4_vis_500_gg
          close(28)
          open (29, file='../CLOUD_DATA/Mg2SiO4_vis_500_qext.txt')
          read(29,*) Mg2SiO4_vis_500_qext
          close(29)
          open (30, file='../CLOUD_DATA/Mg2SiO4_vis_500_qscat.txt')
          read(30,*) Mg2SiO4_vis_500_qscat
          close(30)
          open (31, file='../CLOUD_DATA/Mg2SiO4_vis_650_gg.txt')
          read(31,*) Mg2SiO4_vis_650_gg
          close(31)
          open (32, file='../CLOUD_DATA/Mg2SiO4_vis_650_qext.txt')
          read(32,*) Mg2SiO4_vis_650_qext
          close(32)
          open (33, file='../CLOUD_DATA/Mg2SiO4_vis_650_qscat.txt')
          read(33,*) Mg2SiO4_vis_650_qscat
          close(33)
          open (34, file='../CLOUD_DATA/Mg2SiO4_vis_800_gg.txt')
          read(34,*) Mg2SiO4_vis_800_gg
          close(34)
          open (35, file='../CLOUD_DATA/Mg2SiO4_vis_800_qext.txt')
          read(35,*) Mg2SiO4_vis_800_qext
          close(35)
          open (36, file='../CLOUD_DATA/Mg2SiO4_vis_800_qscat.txt')
          read(36,*) Mg2SiO4_vis_800_qscat
          close(36)
          open (37, file='../CLOUD_DATA/MnS_rosselandMean_gg.txt')
          read(37,*) MnS_rosselandMean_gg
          close(37)
          open (38, file='../CLOUD_DATA/MnS_rosselandMean_qext.txt')
          read(38,*) MnS_rosselandMean_qext
          close(38)
          open (39, file='../CLOUD_DATA/MnS_rosselandMean_qscat.txt')
          read(39,*) MnS_rosselandMean_qscat
          close(39)
          open (40, file='../CLOUD_DATA/MnS_vis_500_gg.txt')
          read(40,*) MnS_vis_500_gg
          close(40)
          open (41, file='../CLOUD_DATA/MnS_vis_500_qext.txt')
          read(41,*) MnS_vis_500_qext
          close(41)
          open (42, file='../CLOUD_DATA/MnS_vis_500_qscat.txt')
          read(42,*) MnS_vis_500_qscat
          close(42)
          open (43, file='../CLOUD_DATA/MnS_vis_650_gg.txt')
          read(43,*) MnS_vis_650_gg
          close(43)
          open (44, file='../CLOUD_DATA/MnS_vis_650_qext.txt')
          read(44,*) MnS_vis_650_qext
          close(44)
          open (45, file='../CLOUD_DATA/MnS_vis_650_qscat.txt')
          read(45,*) MnS_vis_650_qscat
          close(45)
          open (46, file='../CLOUD_DATA/MnS_vis_800_gg.txt')
          read(46,*) MnS_vis_800_gg
          close(46)
          open (47, file='../CLOUD_DATA/MnS_vis_800_qext.txt')
          read(47,*) MnS_vis_800_qext
          close(47)
          open (48, file='../CLOUD_DATA/MnS_vis_800_qscat.txt')
          read(48,*) MnS_vis_800_qscat
          close(48)

          input_particle_size_array_in_meters = (/1.00000000E-07, 1.15139540E-07, 1.32571137E-07,
     &    1.52641797E-07, 1.75751062E-07,
     &    2.02358965E-07, 2.32995181E-07, 2.68269580E-07, 3.08884360E-07, 3.55648031E-07, 4.09491506E-07,
     &    4.71486636E-07, 5.42867544E-07, 6.25055193E-07, 7.19685673E-07, 8.28642773E-07, 9.54095476E-07,
     &    1.09854114E-06, 1.26485522E-06, 1.45634848E-06, 1.67683294E-06, 1.93069773E-06, 2.22299648E-06,
     &    2.55954792E-06, 2.94705170E-06, 3.39322177E-06, 3.90693994E-06, 4.49843267E-06, 5.17947468E-06,
     &    5.96362332E-06, 6.86648845E-06, 7.90604321E-06, 9.10298178E-06, 1.04811313E-05, 1.20679264E-05,
     &    1.38949549E-05, 1.59985872E-05, 1.84206997E-05, 2.12095089E-05, 2.44205309E-05, 2.81176870E-05,
     &    3.23745754E-05, 3.72759372E-05, 4.29193426E-05, 4.94171336E-05, 5.68986603E-05, 6.55128557E-05,
     &    7.54312006E-05, 8.68511374E-05, 1.00000000E-04/)

          input_temperature_array = (/500.00000000, 551.02040816, 602.04081633, 653.06122449, 704.08163265,
     &    755.10204082, 806.12244898, 857.14285714, 908.16326531, 959.18367347, 1010.20408163, 1061.2244898,
     &    1112.24489796, 1163.26530612, 1214.28571429, 1265.30612245, 1316.32653061, 1367.34693878, 1418.36734694,
     &    1469.3877551, 1520.40816327, 1571.42857143, 1622.44897959, 1673.46938776, 1724.48979592, 1775.51020408,
     &    1826.53061224, 1877.55102041, 1928.57142857, 1979.59183673, 2030.6122449, 2081.63265306, 2132.65306122,
     &    2183.67346939, 2234.69387755, 2285.71428571, 2336.73469388, 2387.75510204, 2438.7755102, 2489.79591837,
     &    2540.81632653, 2591.83673469, 2642.85714286, 2693.87755102, 2744.89795918, 2795.91836735, 2846.93877551,
     &    2897.95918367, 2948.97959184, 3000.00000000/)

          particle_size_vs_layer_array_in_meters = (/0.1000E-6, 0.1000E-6,0.1000E-6,
     &    0.1000E-6,0.1000E-6,0.1000E-6,0.1000E-6,0.1000E-6,0.1000E-6,
     &    0.1000E-6, 0.1000E-6,0.1000E-6,0.1000E-6,0.1000E-6,0.1000E-6,0.1000E-6,0.1000E-6,0.1000E-6,0.1023E-6,
     &    0.1060E-6, 0.1108E-6,0.1170E-6,0.1250E-6,0.1360E-6,0.1500E-6,0.1950E-6,0.2285E-6,0.2723E-6,0.3300E-6,
     &    0.4060E-6, 0.5060E-6,0.6387E-6,0.8130E-6,1.0430E-6,1.3458E-6,1.7450E-6,2.2710E-6,2.9660E-6,3.8800E-6,
     &    5.0870E-6, 6.6767E-6,8.7720E-6,11.536E-6,15.1780E-6,19.9800E-6,26.3100E-6,34.6500E-6,45.6500E-6,
     &    60.1540E-6, 79.2700E-6/)


          DENSITY = (/3.95e3, 7.9e3, 3.27e3, 4.0e3/)
          FMOLW   = (/43.20,  23.66, 59.61,  36.87/)

          CORFACT = (/0.005,0.018,0.050,0.135,.367,1.000,1.000,1.000,
     &              1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,
     &              1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,
     &              1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,
     &              1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,
     &              1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,
     &              1.000,1.000,1.000/)
!    Fe
          TconFe = (/1362.14,1373.43, 1384.60,1395.83, 1407.00,1418.26,
     &               1430.00,1442.48, 1455.61, 1469.10, 1482.58,1496.08,
     &               1509.58,1523.08,1536.57,1550.07,1563.57,1577.13,
     &               1590.57,1604.74,1619.69, 1635.41, 1651.59, 1667.75,
     &               1684.03,1700.80, 1718.31, 1736.60, 1755.27, 1773.98,
     &               1792.93, 1812.32,1832.10, 1852.28,1872.54, 1892.90,
     &               1913.24, 1934.27, 1956.41,1978.37, 2008.05, 2030.80,
     &               2051.13, 2081.89, 2103.84, 2132.13,2157.98, 2190.91,
     &               2221.92, 2247.77, 2263.48/)

!    Mg2SiO4
          TCONMg2SiO4 = (/1370.00,1378.31,1386.62,1395.03,1403.56, 1412.29,
     &                  1421.17, 1430.18, 1439.17, 1448.16, 1457.24,1466.52,
     &                  1475.94,1485.57, 1495.23, 1505.09, 1515.04, 1525.21,
     &                  1535.33,1545.80, 1556.37, 1567.12,1578.02, 1589.13,
     &                  1600.29, 1611.67,1623.20, 1634.84,1646.61,1658.58,
     &                  1670.74, 1683.03, 1695.59, 1708.44,1721.29,1734.41,
     &                  1747.86, 1761.64, 1775.79, 1789.95, 1804.36,1819.11,
     &                  1834.34, 1850.41, 1867.40, 1885.24,1903.85,1923.11,
     &                  1943.09, 1963.67, 1974.17/)

!    Al2O3
          TCONAl2O3 =  (/1685.09, 1693.08, 1700.92,1708.84, 1716.78,
     &                 1724.79,1732.91, 1741.42,1750.37,1759.75,
     &                 1769.31,1779.06,1789.37,1799.94, 1810.65,1820.94,
     &                 1830.99,1841.08,1851.02, 1861.02,1871.26,1881.74,
     &                 1892.36,1903.03,1913.69,1924.39,1935.05,1945.95,
     &                 1957.45,1969.45,1980.72,1989.12,1995.54, 2002.05,
     &                 2010.76,2021.69,2033.01,2044.36, 2055.67,2066.99,
     &                 2078.33,2089.65,2100.99,2112.62,2124.88,2137.43,
     &                 2150.32,2163.28, 2176.89, 2191.32,2198.76/)

!     MnS
          TCONMnS = (/1113.43, 1120.07,1126.64,1133.21,1139.9,1146.46,
     &              1153.15,  1159.65, 1166.25, 1173.01, 1180.09, 1187.73,
     &              1194.69, 1202.12, 1209.26, 1216.50, 1223.92, 1231.43,
     &              1239.19, 1247.17, 1255.39,1263.78, 1272.32,1281.40,
     &              1289.69, 1297.76, 1306.00, 1315.24, 1324.18, 1333.27,
     &              1342.44, 1351.39, 1360.54, 1369.99, 1379.72, 1389.42,
     &              1399.22, 1409.04,  1418.99,1428.77, 1438.60, 1449.11,
     &              1459.19, 1469.78, 1481.06, 1492.70, 1504.21, 1515.49,
     &              1527.84,1540.17,1545.90/)

      Tconds(1:51,1)=TconAl2O3
      Tconds(1:51,2)=TconFe
      Tconds(1:51,3)=TconMg2SiO4
      Tconds(1:51,4)=TconMnS

      G0_OPPR(1,1:50,1:50, 1)=Al2O3_vis_500_gg
      G0_OPPR(1,1:50,1:50,2)=Fe_vis_500_gg
      G0_OPPR(1,1:50,1:50,3)=Mg2SiO4_vis_500_gg
      G0_OPPR(1,1:50,1:50,4)=MnS_vis_500_gg

      G0_OPPR(2,1:50,1:50,1)=Al2O3_vis_650_gg
      G0_OPPR(2,1:50,1:50,2)=Fe_vis_650_gg
      G0_OPPR(2,1:50,1:50,3)=Mg2SiO4_vis_650_gg
      G0_OPPR(2,1:50,1:50,4)=MnS_vis_650_gg

      G0_OPPR(3,1:50,1:50,1)=Al2O3_vis_800_gg
      G0_OPPR(3,1:50,1:50,2)=Fe_vis_800_gg
      G0_OPPR(3,1:50,1:50,3)=Mg2SiO4_vis_800_gg
      G0_OPPR(3,1:50,1:50,4)=MnS_vis_800_gg

      G0_OPPR(4,1:50,1:50,1)=Al2O3_rosselandMean_gg
      G0_OPPR(4,1:50,1:50,2)=Fe_rosselandMean_gg
      G0_OPPR(4,1:50,1:50,3)=Mg2SiO4_rosselandMean_gg
      G0_OPPR(4,1:50,1:50,4)=MnS_rosselandMean_gg

      G0_OPPR(5,1:50,1:50,1)=Al2O3_rosselandMean_gg
      G0_OPPR(5,1:50,1:50,2)=Fe_rosselandMean_gg
      G0_OPPR(5,1:50,1:50,3)=Mg2SiO4_rosselandMean_gg
      G0_OPPR(5,1:50,1:50,4)=MnS_rosselandMean_gg

      PI0_OPPR(1,1:50,1:50,1)=Al2O3_vis_500_qscat
      PI0_OPPR(1,1:50,1:50,2)=Fe_vis_500_qscat
      PI0_OPPR(1,1:50,1:50,3)=Mg2SiO4_vis_500_qscat
      PI0_OPPR(1,1:50,1:50,4)=MnS_vis_500_qscat

      PI0_OPPR(2,1:50,1:50,1)=Al2O3_vis_650_qscat
      PI0_OPPR(2,1:50,1:50,2)=Fe_vis_650_qscat
      PI0_OPPR(2,1:50,1:50,3)=Mg2SiO4_vis_650_qscat
      PI0_OPPR(2,1:50,1:50,4)=MnS_vis_650_qscat

      PI0_OPPR(3,1:50,1:50,1)=Al2O3_vis_800_qscat
      PI0_OPPR(3,1:50,1:50,2)=Fe_vis_800_qscat
      PI0_OPPR(3,1:50,1:50,3)=Mg2SiO4_vis_800_qscat
      PI0_OPPR(3,1:50,1:50,4)=MnS_vis_800_qscat

      PI0_OPPR(4,1:50,1:50,1)=Al2O3_rosselandMean_qscat
      PI0_OPPR(4,1:50,1:50,2)=Fe_rosselandMean_qscat
      PI0_OPPR(4,1:50,1:50,3)=Mg2SiO4_rosselandMean_qscat
      PI0_OPPR(4,1:50,1:50,4)=MnS_rosselandMean_qscat

      PI0_OPPR(5,1:50,1:50,1)=Al2O3_rosselandMean_qscat
      PI0_OPPR(5,1:50,1:50,2)=Fe_rosselandMean_qscat
      PI0_OPPR(5,1:50,1:50,3)=Mg2SiO4_rosselandMean_qscat
      PI0_OPPR(5,1:50,1:50,4)=MnS_rosselandMean_qscat

      QE_OPPR(1,1:50,1:50,1)=Al2O3_vis_500_qext
      QE_OPPR(1,1:50,1:50,2)=Fe_vis_500_qext
      QE_OPPR(1,1:50,1:50,3)=Mg2SiO4_vis_500_qext
      QE_OPPR(1,1:50,1:50,4)=MnS_vis_500_qext

      QE_OPPR(2,1:50,1:50,1)=Al2O3_vis_650_qext
      QE_OPPR(2,1:50,1:50,2)=Fe_vis_650_qext
      QE_OPPR(2,1:50,1:50,3)=Mg2SiO4_vis_650_qext
      QE_OPPR(2,1:50,1:50,4)=MnS_vis_650_qext

      QE_OPPR(3,1:50,1:50,1)=Al2O3_vis_800_qext
      QE_OPPR(3,1:50,1:50,2)=Fe_vis_800_qext
      QE_OPPR(3,1:50,1:50,3)=Mg2SiO4_vis_800_qext
      QE_OPPR(3,1:50,1:50,4)=MnS_vis_800_qext

      QE_OPPR(4,1:50,1:50,1)=Al2O3_rosselandMean_qext
      QE_OPPR(4,1:50,1:50,2)=Fe_rosselandMean_qext
      QE_OPPR(4,1:50,1:50,3)=Mg2SiO4_rosselandMean_qext
      QE_OPPR(4,1:50,1:50,4)=MnS_rosselandMean_qext

      QE_OPPR(5,1:50,1:50,1)=Al2O3_rosselandMean_qext
      QE_OPPR(5,1:50,1:50,2)=Fe_rosselandMean_qext
      QE_OPPR(5,1:50,1:50,3)=Mg2SiO4_rosselandMean_qext
      QE_OPPR(5,1:50,1:50,4)=MnS_rosselandMean_qext

      END SUBROUTINE get_cloud_scattering_properties