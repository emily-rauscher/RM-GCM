      SUBROUTINE get_cloud_scattering_properties(NL)
          implicit none
          integer :: J, L, K, NL

          ! Define all the arrays
          ! These are 50 by 50 because that's what the data in CLOUD_DATA is
          ! That can change but use to code from Elsie and Isaac
          real, dimension(50,50) :: KCl_rosselandMean_gg
          real, dimension(50,50) :: KCl_rosselandMean_qext
          real, dimension(50,50) :: KCl_rosselandMean_qscat
          real, dimension(50,50) :: KCl_vis_500_gg
          real, dimension(50,50) :: KCl_vis_500_qext
          real, dimension(50,50) :: KCl_vis_500_qscat
          real, dimension(50,50) :: KCl_vis_650_gg
          real, dimension(50,50) :: KCl_vis_650_qext
          real, dimension(50,50) :: KCl_vis_650_qscat
          real, dimension(50,50) :: KCl_vis_800_gg
          real, dimension(50,50) :: KCl_vis_800_qext
          real, dimension(50,50) :: KCl_vis_800_qscat
          
          real, dimension(50,50) :: ZnS_rosselandMean_gg
          real, dimension(50,50) :: ZnS_rosselandMean_qext
          real, dimension(50,50) :: ZnS_rosselandMean_qscat
          real, dimension(50,50) :: ZnS_vis_500_gg
          real, dimension(50,50) :: ZnS_vis_500_qext
          real, dimension(50,50) :: ZnS_vis_500_qscat
          real, dimension(50,50) :: ZnS_vis_650_gg
          real, dimension(50,50) :: ZnS_vis_650_qext
          real, dimension(50,50) :: ZnS_vis_650_qscat
          real, dimension(50,50) :: ZnS_vis_800_gg
          real, dimension(50,50) :: ZnS_vis_800_qext
          real, dimension(50,50) :: ZnS_vis_800_qscat
          
          real, dimension(50,50) :: Na2S_rosselandMean_gg
          real, dimension(50,50) :: Na2S_rosselandMean_qext
          real, dimension(50,50) :: Na2S_rosselandMean_qscat
          real, dimension(50,50) :: Na2S_vis_500_gg
          real, dimension(50,50) :: Na2S_vis_500_qext
          real, dimension(50,50) :: Na2S_vis_500_qscat
          real, dimension(50,50) :: Na2S_vis_650_gg
          real, dimension(50,50) :: Na2S_vis_650_qext
          real, dimension(50,50) :: Na2S_vis_650_qscat
          real, dimension(50,50) :: Na2S_vis_800_gg
          real, dimension(50,50) :: Na2S_vis_800_qext
          real, dimension(50,50) :: Na2S_vis_800_qscat
          
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
          
          real, dimension(50,50) :: Cr_rosselandMean_gg
          real, dimension(50,50) :: Cr_rosselandMean_qext
          real, dimension(50,50) :: Cr_rosselandMean_qscat
          real, dimension(50,50) :: Cr_vis_500_gg
          real, dimension(50,50) :: Cr_vis_500_qext
          real, dimension(50,50) :: Cr_vis_500_qscat
          real, dimension(50,50) :: Cr_vis_650_gg
          real, dimension(50,50) :: Cr_vis_650_qext
          real, dimension(50,50) :: Cr_vis_650_qscat
          real, dimension(50,50) :: Cr_vis_800_gg
          real, dimension(50,50) :: Cr_vis_800_qext
          real, dimension(50,50) :: Cr_vis_800_qscat
          
          real, dimension(50,50) :: SiO2_rosselandMean_gg
          real, dimension(50,50) :: SiO2_rosselandMean_qext
          real, dimension(50,50) :: SiO2_rosselandMean_qscat
          real, dimension(50,50) :: SiO2_vis_500_gg
          real, dimension(50,50) :: SiO2_vis_500_qext
          real, dimension(50,50) :: SiO2_vis_500_qscat
          real, dimension(50,50) :: SiO2_vis_650_gg
          real, dimension(50,50) :: SiO2_vis_650_qext
          real, dimension(50,50) :: SiO2_vis_650_qscat
          real, dimension(50,50) :: SiO2_vis_800_gg
          real, dimension(50,50) :: SiO2_vis_800_qext
          real, dimension(50,50) :: SiO2_vis_800_qscat
          
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
          
          real, dimension(50,50) :: VO_rosselandMean_gg
          real, dimension(50,50) :: VO_rosselandMean_qext
          real, dimension(50,50) :: VO_rosselandMean_qscat
          real, dimension(50,50) :: VO_vis_500_gg
          real, dimension(50,50) :: VO_vis_500_qext
          real, dimension(50,50) :: VO_vis_500_qscat
          real, dimension(50,50) :: VO_vis_650_gg
          real, dimension(50,50) :: VO_vis_650_qext
          real, dimension(50,50) :: VO_vis_650_qscat
          real, dimension(50,50) :: VO_vis_800_gg
          real, dimension(50,50) :: VO_vis_800_qext
          real, dimension(50,50) :: VO_vis_800_qscat
          
          real, dimension(50,50) :: C2H2_rosselandMean_gg
          real, dimension(50,50) :: C2H2_rosselandMean_qext
          real, dimension(50,50) :: C2H2_rosselandMean_qscat
          real, dimension(50,50) :: C2H2_vis_500_gg
          real, dimension(50,50) :: C2H2_vis_500_qext
          real, dimension(50,50) :: C2H2_vis_500_qscat
          real, dimension(50,50) :: C2H2_vis_650_gg
          real, dimension(50,50) :: C2H2_vis_650_qext
          real, dimension(50,50) :: C2H2_vis_650_qscat
          real, dimension(50,50) :: C2H2_vis_800_gg
          real, dimension(50,50) :: C2H2_vis_800_qext
          real, dimension(50,50) :: C2H2_vis_800_qscat
          
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
          
          real, dimension(50,50) :: Fe2SiO4_rosselandMean_gg
          real, dimension(50,50) :: Fe2SiO4_rosselandMean_qext
          real, dimension(50,50) :: Fe2SiO4_rosselandMean_qscat
          real, dimension(50,50) :: Fe2SiO4_vis_500_gg
          real, dimension(50,50) :: Fe2SiO4_vis_500_qext
          real, dimension(50,50) :: Fe2SiO4_vis_500_qscat
          real, dimension(50,50) :: Fe2SiO4_vis_650_gg
          real, dimension(50,50) :: Fe2SiO4_vis_650_qext
          real, dimension(50,50) :: Fe2SiO4_vis_650_qscat
          real, dimension(50,50) :: Fe2SiO4_vis_800_gg
          real, dimension(50,50) :: Fe2SiO4_vis_800_qext
          real, dimension(50,50) :: Fe2SiO4_vis_800_qscat
          
          real, dimension(50,50) :: CaTiO3_rosselandMean_gg
          real, dimension(50,50) :: CaTiO3_rosselandMean_qext
          real, dimension(50,50) :: CaTiO3_rosselandMean_qscat
          real, dimension(50,50) :: CaTiO3_vis_500_gg
          real, dimension(50,50) :: CaTiO3_vis_500_qext
          real, dimension(50,50) :: CaTiO3_vis_500_qscat
          real, dimension(50,50) :: CaTiO3_vis_650_gg
          real, dimension(50,50) :: CaTiO3_vis_650_qext
          real, dimension(50,50) :: CaTiO3_vis_650_qscat
          real, dimension(50,50) :: CaTiO3_vis_800_gg
          real, dimension(50,50) :: CaTiO3_vis_800_qext
          real, dimension(50,50) :: CaTiO3_vis_800_qscat
          
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


          REAL TconKCl(51) 
          REAL TconZnS(51)
          REAL TconNa2S(51)
          REAL TCONMnS(51)
          REAL TconCr2O3(51)
          REAL TCONSiO2(51)
          REAL TCONMg2SiO4(51) 
          REAL TconVO(51)
          REAL TconNi(51)
          REAL TCONFe(51)
          REAL TconCa2SiO4(51)
          REAL TconCaTiO3(51)
          REAL TCONAl2O3(51)

          REAL CORFACT(51)
          REAL TCONDS(51, 13)

          REAL QE_OPPR(5, 50, 50, 13)
          REAL PI0_OPPR(5, 50, 50, 13)
          REAL G0_OPPR(5, 50, 50, 13)

          REAL DENSITY(13)
          REAL FMOLW(13)
          REAL MOLEF(13)

          real, dimension(50) :: input_particle_size_array_in_meters
          real, dimension(50) :: input_temperature_array
          real, dimension(50) :: particle_size_vs_layer_array_in_meters

          COMMON /CLOUD_PROPERTIES/ TCONDS, QE_OPPR, PI0_OPPR, G0_OPPR,
     &                              DENSITY, FMOLW, MOLEF,
     &                              CORFACT,
     &                              input_particle_size_array_in_meters,
     &                              input_temperature_array,
     &                              particle_size_vs_layer_array_in_meters

          ! opening the file for reading
          open (1, file='../CLOUD_DATA/KCl_rosselandMean_gg.txt')
          open (2, file='../CLOUD_DATA/KCl_rosselandMean_qext.txt')
          open (3, file='../CLOUD_DATA/KCl_rosselandMean_qscat.txt')
          open (4, file='../CLOUD_DATA/KCl_vis_500_gg.txt')
          open (5, file='../CLOUD_DATA/KCl_vis_500_qext.txt')
          open (6, file='../CLOUD_DATA/KCl_vis_500_qscat.txt')
          open (7, file='../CLOUD_DATA/KCl_vis_650_gg.txt')
          open (8, file='../CLOUD_DATA/KCl_vis_650_qext.txt')
          open (9, file='../CLOUD_DATA/KCl_vis_650_qscat.txt')
          open (10, file='../CLOUD_DATA/KCl_vis_800_gg.txt')
          open (11, file='../CLOUD_DATA/KCl_vis_800_qext.txt')
          open (12, file='../CLOUD_DATA/KCl_vis_800_qscat.txt')
          
          read(1,*) KCl_rosselandMean_gg
          read(2,*) KCl_rosselandMean_qext
          read(3,*) KCl_rosselandMean_qscat
          read(4,*) KCl_vis_500_gg
          read(5,*) KCl_vis_500_qext
          read(6,*) KCl_vis_500_qscat
          read(7,*) KCl_vis_650_gg
          read(8,*) KCl_vis_650_qext
          read(9,*) KCl_vis_650_qscat
          read(10,*) KCl_vis_800_gg
          read(11,*) KCl_vis_800_qext
          read(12,*) KCl_vis_800_qscat
          
          close(1)
          close(2)
          close(3)
          close(4)
          close(5)
          close(6)
          close(7)
          close(8)
          close(9)
          close(10)
          close(11)
          close(12)
          
          open (1, file='../CLOUD_DATA/ZnS_rosselandMean_gg.txt')
          open (2, file='../CLOUD_DATA/ZnS_rosselandMean_qext.txt')
          open (3, file='../CLOUD_DATA/ZnS_rosselandMean_qscat.txt')
          open (4, file='../CLOUD_DATA/ZnS_vis_500_gg.txt')
          open (5, file='../CLOUD_DATA/ZnS_vis_500_qext.txt')
          open (6, file='../CLOUD_DATA/ZnS_vis_500_qscat.txt')
          open (7, file='../CLOUD_DATA/ZnS_vis_650_gg.txt')
          open (8, file='../CLOUD_DATA/ZnS_vis_650_qext.txt')
          open (9, file='../CLOUD_DATA/ZnS_vis_650_qscat.txt')
          open (10, file='../CLOUD_DATA/ZnS_vis_800_gg.txt')
          open (11, file='../CLOUD_DATA/ZnS_vis_800_qext.txt')
          open (12, file='../CLOUD_DATA/ZnS_vis_800_qscat.txt')
          
          read(1,*) ZnS_rosselandMean_gg
          read(2,*) ZnS_rosselandMean_qext
          read(3,*) ZnS_rosselandMean_qscat
          read(4,*) ZnS_vis_500_gg
          read(5,*) ZnS_vis_500_qext
          read(6,*) ZnS_vis_500_qscat
          read(7,*) ZnS_vis_650_gg
          read(8,*) ZnS_vis_650_qext
          read(9,*) ZnS_vis_650_qscat
          read(10,*) ZnS_vis_800_gg
          read(11,*) ZnS_vis_800_qext
          read(12,*) ZnS_vis_800_qscat
          
          close(1)
          close(2)
          close(3)
          close(4)
          close(5)
          close(6)
          close(7)
          close(8)
          close(9)
          close(10)
          close(11)
          close(12)
          
          open (1, file='../CLOUD_DATA/Na2S_rosselandMean_gg.txt')
          open (2, file='../CLOUD_DATA/Na2S_rosselandMean_qext.txt')
          open (3, file='../CLOUD_DATA/Na2S_rosselandMean_qscat.txt')
          open (4, file='../CLOUD_DATA/Na2S_vis_500_gg.txt')
          open (5, file='../CLOUD_DATA/Na2S_vis_500_qext.txt')
          open (6, file='../CLOUD_DATA/Na2S_vis_500_qscat.txt')
          open (7, file='../CLOUD_DATA/Na2S_vis_650_gg.txt')
          open (8, file='../CLOUD_DATA/Na2S_vis_650_qext.txt')
          open (9, file='../CLOUD_DATA/Na2S_vis_650_qscat.txt')
          open (10, file='../CLOUD_DATA/Na2S_vis_800_gg.txt')
          open (11, file='../CLOUD_DATA/Na2S_vis_800_qext.txt')
          open (12, file='../CLOUD_DATA/Na2S_vis_800_qscat.txt')
          
          read(1,*) Na2S_rosselandMean_gg
          read(2,*) Na2S_rosselandMean_qext
          read(3,*) Na2S_rosselandMean_qscat
          read(4,*) Na2S_vis_500_gg
          read(5,*) Na2S_vis_500_qext
          read(6,*) Na2S_vis_500_qscat
          read(7,*) Na2S_vis_650_gg
          read(8,*) Na2S_vis_650_qext
          read(9,*) Na2S_vis_650_qscat
          read(10,*) Na2S_vis_800_gg
          read(11,*) Na2S_vis_800_qext
          read(12,*) Na2S_vis_800_qscat
          
          close(1)
          close(2)
          close(3)
          close(4)
          close(5)
          close(6)
          close(7)
          close(8)
          close(9)
          close(10)
          close(11)
          close(12)
          
          open (1, file='../CLOUD_DATA/MnS_rosselandMean_gg.txt')
          open (2, file='../CLOUD_DATA/MnS_rosselandMean_qext.txt')
          open (3, file='../CLOUD_DATA/MnS_rosselandMean_qscat.txt')
          open (4, file='../CLOUD_DATA/MnS_vis_500_gg.txt')
          open (5, file='../CLOUD_DATA/MnS_vis_500_qext.txt')
          open (6, file='../CLOUD_DATA/MnS_vis_500_qscat.txt')
          open (7, file='../CLOUD_DATA/MnS_vis_650_gg.txt')
          open (8, file='../CLOUD_DATA/MnS_vis_650_qext.txt')
          open (9, file='../CLOUD_DATA/MnS_vis_650_qscat.txt')
          open (10, file='../CLOUD_DATA/MnS_vis_800_gg.txt')
          open (11, file='../CLOUD_DATA/MnS_vis_800_qext.txt')
          open (12, file='../CLOUD_DATA/MnS_vis_800_qscat.txt')
          
          read(1,*) MnS_rosselandMean_gg
          read(2,*) MnS_rosselandMean_qext
          read(3,*) MnS_rosselandMean_qscat
          read(4,*) MnS_vis_500_gg
          read(5,*) MnS_vis_500_qext
          read(6,*) MnS_vis_500_qscat
          read(7,*) MnS_vis_650_gg
          read(8,*) MnS_vis_650_qext
          read(9,*) MnS_vis_650_qscat
          read(10,*) MnS_vis_800_gg
          read(11,*) MnS_vis_800_qext
          read(12,*) MnS_vis_800_qscat
          
          close(1)
          close(2)
          close(3)
          close(4)
          close(5)
          close(6)
          close(7)
          close(8)
          close(9)
          close(10)
          close(11)
          close(12)
          
          open (1, file='../CLOUD_DATA/Cr_rosselandMean_gg.txt')
          open (2, file='../CLOUD_DATA/Cr_rosselandMean_qext.txt')
          open (3, file='../CLOUD_DATA/Cr_rosselandMean_qscat.txt')
          open (4, file='../CLOUD_DATA/Cr_vis_500_gg.txt')
          open (5, file='../CLOUD_DATA/Cr_vis_500_qext.txt')
          open (6, file='../CLOUD_DATA/Cr_vis_500_qscat.txt')
          open (7, file='../CLOUD_DATA/Cr_vis_650_gg.txt')
          open (8, file='../CLOUD_DATA/Cr_vis_650_qext.txt')
          open (9, file='../CLOUD_DATA/Cr_vis_650_qscat.txt')
          open (10, file='../CLOUD_DATA/Cr_vis_800_gg.txt')
          open (11, file='../CLOUD_DATA/Cr_vis_800_qext.txt')
          open (12, file='../CLOUD_DATA/Cr_vis_800_qscat.txt')
          
          read(1,*) Cr_rosselandMean_gg
          read(2,*) Cr_rosselandMean_qext
          read(3,*) Cr_rosselandMean_qscat
          read(4,*) Cr_vis_500_gg
          read(5,*) Cr_vis_500_qext
          read(6,*) Cr_vis_500_qscat
          read(7,*) Cr_vis_650_gg
          read(8,*) Cr_vis_650_qext
          read(9,*) Cr_vis_650_qscat
          read(10,*) Cr_vis_800_gg
          read(11,*) Cr_vis_800_qext
          read(12,*) Cr_vis_800_qscat
          
          close(1)
          close(2)
          close(3)
          close(4)
          close(5)
          close(6)
          close(7)
          close(8)
          close(9)
          close(10)
          close(11)
          close(12)
          
          open (1, file='../CLOUD_DATA/SiO2_rosselandMean_gg.txt')
          open (2, file='../CLOUD_DATA/SiO2_rosselandMean_qext.txt')
          open (3, file='../CLOUD_DATA/SiO2_rosselandMean_qscat.txt')
          open (4, file='../CLOUD_DATA/SiO2_vis_500_gg.txt')
          open (5, file='../CLOUD_DATA/SiO2_vis_500_qext.txt')
          open (6, file='../CLOUD_DATA/SiO2_vis_500_qscat.txt')
          open (7, file='../CLOUD_DATA/SiO2_vis_650_gg.txt')
          open (8, file='../CLOUD_DATA/SiO2_vis_650_qext.txt')
          open (9, file='../CLOUD_DATA/SiO2_vis_650_qscat.txt')
          open (10, file='../CLOUD_DATA/SiO2_vis_800_gg.txt')
          open (11, file='../CLOUD_DATA/SiO2_vis_800_qext.txt')
          open (12, file='../CLOUD_DATA/SiO2_vis_800_qscat.txt')
          
          read(1,*) SiO2_rosselandMean_gg
          read(2,*) SiO2_rosselandMean_qext
          read(3,*) SiO2_rosselandMean_qscat
          read(4,*) SiO2_vis_500_gg
          read(5,*) SiO2_vis_500_qext
          read(6,*) SiO2_vis_500_qscat
          read(7,*) SiO2_vis_650_gg
          read(8,*) SiO2_vis_650_qext
          read(9,*) SiO2_vis_650_qscat
          read(10,*) SiO2_vis_800_gg
          read(11,*) SiO2_vis_800_qext
          read(12,*) SiO2_vis_800_qscat
          
          close(1)
          close(2)
          close(3)
          close(4)
          close(5)
          close(6)
          close(7)
          close(8)
          close(9)
          close(10)
          close(11)
          close(12)
          
          open (1, file='../CLOUD_DATA/Mg2SiO4_rosselandMean_gg.txt')
          open (2, file='../CLOUD_DATA/Mg2SiO4_rosselandMean_qext.txt')
          open (3, file='../CLOUD_DATA/Mg2SiO4_rosselandMean_qscat.txt')
          open (4, file='../CLOUD_DATA/Mg2SiO4_vis_500_gg.txt')
          open (5, file='../CLOUD_DATA/Mg2SiO4_vis_500_qext.txt')
          open (6, file='../CLOUD_DATA/Mg2SiO4_vis_500_qscat.txt')
          open (7, file='../CLOUD_DATA/Mg2SiO4_vis_650_gg.txt')
          open (8, file='../CLOUD_DATA/Mg2SiO4_vis_650_qext.txt')
          open (9, file='../CLOUD_DATA/Mg2SiO4_vis_650_qscat.txt')
          open (10, file='../CLOUD_DATA/Mg2SiO4_vis_800_gg.txt')
          open (11, file='../CLOUD_DATA/Mg2SiO4_vis_800_qext.txt')
          open (12, file='../CLOUD_DATA/Mg2SiO4_vis_800_qscat.txt')
          
          read(1,*) Mg2SiO4_rosselandMean_gg
          read(2,*) Mg2SiO4_rosselandMean_qext
          read(3,*) Mg2SiO4_rosselandMean_qscat
          read(4,*) Mg2SiO4_vis_500_gg
          read(5,*) Mg2SiO4_vis_500_qext
          read(6,*) Mg2SiO4_vis_500_qscat
          read(7,*) Mg2SiO4_vis_650_gg
          read(8,*) Mg2SiO4_vis_650_qext
          read(9,*) Mg2SiO4_vis_650_qscat
          read(10,*) Mg2SiO4_vis_800_gg
          read(11,*) Mg2SiO4_vis_800_qext
          read(12,*) Mg2SiO4_vis_800_qscat
          
          close(1)
          close(2)
          close(3)
          close(4)
          close(5)
          close(6)
          close(7)
          close(8)
          close(9)
          close(10)
          close(11)
          close(12)
          
          open (1, file='../CLOUD_DATA/VO_rosselandMean_gg.txt')
          open (2, file='../CLOUD_DATA/VO_rosselandMean_qext.txt')
          open (3, file='../CLOUD_DATA/VO_rosselandMean_qscat.txt')
          open (4, file='../CLOUD_DATA/VO_vis_500_gg.txt')
          open (5, file='../CLOUD_DATA/VO_vis_500_qext.txt')
          open (6, file='../CLOUD_DATA/VO_vis_500_qscat.txt')
          open (7, file='../CLOUD_DATA/VO_vis_650_gg.txt')
          open (8, file='../CLOUD_DATA/VO_vis_650_qext.txt')
          open (9, file='../CLOUD_DATA/VO_vis_650_qscat.txt')
          open (10, file='../CLOUD_DATA/VO_vis_800_gg.txt')
          open (11, file='../CLOUD_DATA/VO_vis_800_qext.txt')
          open (12, file='../CLOUD_DATA/VO_vis_800_qscat.txt')
          
          read(1,*) VO_rosselandMean_gg
          read(2,*) VO_rosselandMean_qext
          read(3,*) VO_rosselandMean_qscat
          read(4,*) VO_vis_500_gg
          read(5,*) VO_vis_500_qext
          read(6,*) VO_vis_500_qscat
          read(7,*) VO_vis_650_gg
          read(8,*) VO_vis_650_qext
          read(9,*) VO_vis_650_qscat
          read(10,*) VO_vis_800_gg
          read(11,*) VO_vis_800_qext
          read(12,*) VO_vis_800_qscat
          
          close(1)
          close(2)
          close(3)
          close(4)
          close(5)
          close(6)
          close(7)
          close(8)
          close(9)
          close(10)
          close(11)
          close(12)
          
          open (1, file='../CLOUD_DATA/C2H2_rosselandMean_gg.txt')
          open (2, file='../CLOUD_DATA/C2H2_rosselandMean_qext.txt')
          open (3, file='../CLOUD_DATA/C2H2_rosselandMean_qscat.txt')
          open (4, file='../CLOUD_DATA/C2H2_vis_500_gg.txt')
          open (5, file='../CLOUD_DATA/C2H2_vis_500_qext.txt')
          open (6, file='../CLOUD_DATA/C2H2_vis_500_qscat.txt')
          open (7, file='../CLOUD_DATA/C2H2_vis_650_gg.txt')
          open (8, file='../CLOUD_DATA/C2H2_vis_650_qext.txt')
          open (9, file='../CLOUD_DATA/C2H2_vis_650_qscat.txt')
          open (10, file='../CLOUD_DATA/C2H2_vis_800_gg.txt')
          open (11, file='../CLOUD_DATA/C2H2_vis_800_qext.txt')
          open (12, file='../CLOUD_DATA/C2H2_vis_800_qscat.txt')
          
          read(1,*) C2H2_rosselandMean_gg
          read(2,*) C2H2_rosselandMean_qext
          read(3,*) C2H2_rosselandMean_qscat
          read(4,*) C2H2_vis_500_gg
          read(5,*) C2H2_vis_500_qext
          read(6,*) C2H2_vis_500_qscat
          read(7,*) C2H2_vis_650_gg
          read(8,*) C2H2_vis_650_qext
          read(9,*) C2H2_vis_650_qscat
          read(10,*) C2H2_vis_800_gg
          read(11,*) C2H2_vis_800_qext
          read(12,*) C2H2_vis_800_qscat
          
          close(1)
          close(2)
          close(3)
          close(4)
          close(5)
          close(6)
          close(7)
          close(8)
          close(9)
          close(10)
          close(11)
          close(12)
          
          open (1, file='../CLOUD_DATA/Fe_rosselandMean_gg.txt')
          open (2, file='../CLOUD_DATA/Fe_rosselandMean_qext.txt')
          open (3, file='../CLOUD_DATA/Fe_rosselandMean_qscat.txt')
          open (4, file='../CLOUD_DATA/Fe_vis_500_gg.txt')
          open (5, file='../CLOUD_DATA/Fe_vis_500_qext.txt')
          open (6, file='../CLOUD_DATA/Fe_vis_500_qscat.txt')
          open (7, file='../CLOUD_DATA/Fe_vis_650_gg.txt')
          open (8, file='../CLOUD_DATA/Fe_vis_650_qext.txt')
          open (9, file='../CLOUD_DATA/Fe_vis_650_qscat.txt')
          open (10, file='../CLOUD_DATA/Fe_vis_800_gg.txt')
          open (11, file='../CLOUD_DATA/Fe_vis_800_qext.txt')
          open (12, file='../CLOUD_DATA/Fe_vis_800_qscat.txt')
          
          read(1,*) Fe_rosselandMean_gg
          read(2,*) Fe_rosselandMean_qext
          read(3,*) Fe_rosselandMean_qscat
          read(4,*) Fe_vis_500_gg
          read(5,*) Fe_vis_500_qext
          read(6,*) Fe_vis_500_qscat
          read(7,*) Fe_vis_650_gg
          read(8,*) Fe_vis_650_qext
          read(9,*) Fe_vis_650_qscat
          read(10,*) Fe_vis_800_gg
          read(11,*) Fe_vis_800_qext
          read(12,*) Fe_vis_800_qscat
          
          close(1)
          close(2)
          close(3)
          close(4)
          close(5)
          close(6)
          close(7)
          close(8)
          close(9)
          close(10)
          close(11)
          close(12)
          
          open (1, file='../CLOUD_DATA/Fe2SiO4_rosselandMean_gg.txt')
          open (2, file='../CLOUD_DATA/Fe2SiO4_rosselandMean_qext.txt')
          open (3, file='../CLOUD_DATA/Fe2SiO4_rosselandMean_qscat.txt')
          open (4, file='../CLOUD_DATA/Fe2SiO4_vis_500_gg.txt')
          open (5, file='../CLOUD_DATA/Fe2SiO4_vis_500_qext.txt')
          open (6, file='../CLOUD_DATA/Fe2SiO4_vis_500_qscat.txt')
          open (7, file='../CLOUD_DATA/Fe2SiO4_vis_650_gg.txt')
          open (8, file='../CLOUD_DATA/Fe2SiO4_vis_650_qext.txt')
          open (9, file='../CLOUD_DATA/Fe2SiO4_vis_650_qscat.txt')
          open (10, file='../CLOUD_DATA/Fe2SiO4_vis_800_gg.txt')
          open (11, file='../CLOUD_DATA/Fe2SiO4_vis_800_qext.txt')
          open (12, file='../CLOUD_DATA/Fe2SiO4_vis_800_qscat.txt')
          
          read(1,*) Fe2SiO4_rosselandMean_gg
          read(2,*) Fe2SiO4_rosselandMean_qext
          read(3,*) Fe2SiO4_rosselandMean_qscat
          read(4,*) Fe2SiO4_vis_500_gg
          read(5,*) Fe2SiO4_vis_500_qext
          read(6,*) Fe2SiO4_vis_500_qscat
          read(7,*) Fe2SiO4_vis_650_gg
          read(8,*) Fe2SiO4_vis_650_qext
          read(9,*) Fe2SiO4_vis_650_qscat
          read(10,*) Fe2SiO4_vis_800_gg
          read(11,*) Fe2SiO4_vis_800_qext
          read(12,*) Fe2SiO4_vis_800_qscat
          
          close(1)
          close(2)
          close(3)
          close(4)
          close(5)
          close(6)
          close(7)
          close(8)
          close(9)
          close(10)
          close(11)
          close(12)
          
          open (1, file='../CLOUD_DATA/CaTiO3_rosselandMean_gg.txt')
          open (2, file='../CLOUD_DATA/CaTiO3_rosselandMean_qext.txt')
          open (3, file='../CLOUD_DATA/CaTiO3_rosselandMean_qscat.txt')
          open (4, file='../CLOUD_DATA/CaTiO3_vis_500_gg.txt')
          open (5, file='../CLOUD_DATA/CaTiO3_vis_500_qext.txt')
          open (6, file='../CLOUD_DATA/CaTiO3_vis_500_qscat.txt')
          open (7, file='../CLOUD_DATA/CaTiO3_vis_650_gg.txt')
          open (8, file='../CLOUD_DATA/CaTiO3_vis_650_qext.txt')
          open (9, file='../CLOUD_DATA/CaTiO3_vis_650_qscat.txt')
          open (10, file='../CLOUD_DATA/CaTiO3_vis_800_gg.txt')
          open (11, file='../CLOUD_DATA/CaTiO3_vis_800_qext.txt')
          open (12, file='../CLOUD_DATA/CaTiO3_vis_800_qscat.txt')
          
          read(1,*) CaTiO3_rosselandMean_gg
          read(2,*) CaTiO3_rosselandMean_qext
          read(3,*) CaTiO3_rosselandMean_qscat
          read(4,*) CaTiO3_vis_500_gg
          read(5,*) CaTiO3_vis_500_qext
          read(6,*) CaTiO3_vis_500_qscat
          read(7,*) CaTiO3_vis_650_gg
          read(8,*) CaTiO3_vis_650_qext
          read(9,*) CaTiO3_vis_650_qscat
          read(10,*) CaTiO3_vis_800_gg
          read(11,*) CaTiO3_vis_800_qext
          read(12,*) CaTiO3_vis_800_qscat
          
          close(1)
          close(2)
          close(3)
          close(4)
          close(5)
          close(6)
          close(7)
          close(8)
          close(9)
          close(10)
          close(11)
          close(12)
          
          open (1, file='../CLOUD_DATA/Al2O3_rosselandMean_gg.txt')
          open (2, file='../CLOUD_DATA/Al2O3_rosselandMean_qext.txt')
          open (3, file='../CLOUD_DATA/Al2O3_rosselandMean_qscat.txt')
          open (4, file='../CLOUD_DATA/Al2O3_vis_500_gg.txt')
          open (5, file='../CLOUD_DATA/Al2O3_vis_500_qext.txt')
          open (6, file='../CLOUD_DATA/Al2O3_vis_500_qscat.txt')
          open (7, file='../CLOUD_DATA/Al2O3_vis_650_gg.txt')
          open (8, file='../CLOUD_DATA/Al2O3_vis_650_qext.txt')
          open (9, file='../CLOUD_DATA/Al2O3_vis_650_qscat.txt')
          open (10, file='../CLOUD_DATA/Al2O3_vis_800_gg.txt')
          open (11, file='../CLOUD_DATA/Al2O3_vis_800_qext.txt')
          open (12, file='../CLOUD_DATA/Al2O3_vis_800_qscat.txt')
          
          read(1,*) Al2O3_rosselandMean_gg
          read(2,*) Al2O3_rosselandMean_qext
          read(3,*) Al2O3_rosselandMean_qscat
          read(4,*) Al2O3_vis_500_gg
          read(5,*) Al2O3_vis_500_qext
          read(6,*) Al2O3_vis_500_qscat
          read(7,*) Al2O3_vis_650_gg
          read(8,*) Al2O3_vis_650_qext
          read(9,*) Al2O3_vis_650_qscat
          read(10,*) Al2O3_vis_800_gg
          read(11,*) Al2O3_vis_800_qext
          read(12,*) Al2O3_vis_800_qscat
          
          close(1)
          close(2)
          close(3)
          close(4)
          close(5)
          close(6)
          close(7)
          close(8)
          close(9)
          close(10)
          close(11)
          close(12)

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


          DENSITY = (/1.98e3,4.09e3,1.86e3,4.0e3,5.22e3,2.65e3,3.27e3,5.76e3,8.9e3, 7.9e3,3.34e3,3.98e3,3.95e3/)
          FMOLW   = (/31.59,41.30,33.07,36.87,64.40,25.46,59.61,28.37,24.87, 23.66, 72.99, 50.83, 43.20/)
          MOLEF   = (/1.23e-07,4.06e-08,9.35e-07,3.11e-07,4.4e-07,3.26e-05,
     &                1.745e-05,9.56e-09,1.61e-06,2.94e-05,9.95e-07,7.83e-08,1.385e-06/)

          CORFACT = (/0.005,0.018,0.050,0.135,.367,1.000,1.000,1.000,
     &              1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,
     &              1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,
     &              1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,
     &              1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,
     &              1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,
     &              1.000,1.000,1.000/)



!    KCl
      TconKCl = (/617.032,621.573,626.038,630.552,635.053, 639.555,
     &         644.050,648.556,653.049,657.552,662.043,666.609,671.436,
     &         676.532,681.879,687.233,692.462,697.665,702.916,708.306,
     &         713.767,719.366,725.024,730.775,736.460,742.266,748.065,
     &        753.932,759.779,765.571, 771.346,777.201,783.301,789.715,
     &         796.379,803.117,809.863,816.737,823.798,831.052,838.426,
     &         845.980,853.873,862.074,870.494,878.913,887.351,895.768,
     &           904.198,912.867,917.385/)


!ZnS
      TconZnS =(/708.296,712.010,716.704,719.507,722.309,727.329,
     &                730.718,736.323,739.126,744.731,747.534,753.139,
     &                755.942,761.547,764.350,769.955,772.758,778.363,
     &                783.969,788.079,792.377,797.982,803.587,806.390, 
     &                811.995,815.183,823.206,828.812,834.364,840.022, 
     &                845.628,851.233,856.839,862.444,868.049,873.655,
     &                879.260,884.865,890.471,896.076,901.682,907.287,
     &                915.695,921.897,929.708,936.547,943.722,949.327,
     &                955.732,963.063, 966.143/)


!    Na2S
      TconNa2S= (/776.182,781.957,787.715,793.430,799.135,804.761,
     &              810.359, 815.868, 821.297, 826.731, 832.157,837.761,
     &              843.767,850.179,856.911,863.688,870.458,877.165,
     &              884.006, 890.930,897.948, 905.112,912.216,919.374,
     &              926.464,933.659,940.930,948.253,955.902,963.944,
     &          972.372,980.941,989.400,998.113,1007.19,1016.45,1025.51,
     &          1035.47,1044.52,1054.07,1063.67,1073.49,1083.32,1093.15,
     &              1103.16,1113.41,1123.84,1134.53,1145.92,1158.12,
     &              1164.38/)



!     MnS
      TCONMnS =(/1113.43, 1120.07,1126.64,1133.21,1139.9,1146.46,
     &        1153.15,  1159.65, 1166.25, 1173.01, 1180.09, 1187.73,
     &        1194.69, 1202.12, 1209.26, 1216.50, 1223.92, 1231.43,
     &       1239.19, 1247.17, 1255.39,1263.78, 1272.32,1281.40,
     &        1289.69, 1297.76, 1306.00, 1315.24, 1324.18, 1333.27,
     &        1342.44, 1351.39, 1360.54, 1369.99, 1379.72, 1389.42,
     &        1399.22, 1409.04,  1418.99,1428.77, 1438.60, 1449.11,
     &        1459.19, 1469.78, 1481.06, 1492.70, 1504.21, 1515.49,
     &            1527.84,1540.17,1545.90/)


!    Cr2O3
      TconCr2O3 =(/1213.17,1219.05,1224.98,1231.07,1237.15,1243.21,
     &        1249.26,1255.35,1261.51,1267.83,1274.22,1280.63,1287.04,
     &       1293.56,1300.89,1309.34,1318.81,1329.04,1339.04,1349.79,
     &       1361.13,1373.13,1385.33,1397.48,1409.69,1421.78,1434.01,
     &      1446.16,1458.55,1471.19,1483.84,1496.49,1508.99,1522.14,
     &      1536.54,1552.17,1568.58,1585.09,1601.61,1618.14,1634.62,
     &      1651.02,1667.41,1683.80,1700.20,1716.57,1732.89,1749.26,
     &             1765.73,1783.15,1792.19/)



!    SiO2
      TCONSiO2 =(/1334.63,1342.58,1350.30,1358.48,1366.64,1374.85,
     &             1383.15,1391.59,1400.10,1408.68, 1417.25,1425.87,
     &             1434.53, 1443.14, 1451.71,1460.28,1468.90,1477.44,
     &             1486.12,1494.77,1503.91,1513.72 ,1524.07, 1534.57,
     &             1544.98,1555.25, 1565.35,1575.40, 1585.44,1595.55,
     &             1606.04,1617.06, 1628.55,1640.25,1651.98, 1664.07,
     &             1676.82,1689.97,1703.53, 1717.16, 1731.36, 1746.01,
     &             1761.10,1776.94, 1793.70,1811.36,1829.60, 1848.37,
     &             1867.54, 1887.00, 1896.46/)


 

!    Mg2SiO4
      TCONMg2SiO4 =(/1370.00,1378.31,1386.62,1395.03,1403.56,
     &           1412.29,
     &           1421.17, 1430.18, 1439.17, 1448.16, 1457.24,1466.52,
     &           1475.94,1485.57, 1495.23, 1505.09, 1515.04, 1525.21,
     &           1535.33,1545.80, 1556.37, 1567.12,1578.02, 1589.13,
     &           1600.29, 1611.67,1623.20, 1634.84,1646.61,1658.58,
     &           1670.74, 1683.03, 1695.59, 1708.44,1721.29,1734.41,
     &           1747.86, 1761.64, 1775.79, 1789.95, 1804.36,1819.11,
     &           1834.34, 1850.41, 1867.40, 1885.24,1903.85,1923.11,
     &           1943.09, 1963.67, 1974.17/)

    
!    VO
      TconVO =(/1363.07,1371.53, 1379.98, 1389.36,1399.69,1409.65,
     &            1419.40, 1429.94, 1439.11,1450.21, 1458.82,1470.07, 
     &       1481.33,1492.58,1503.61,1513.72,1526.34,1536.53,1547.03,
     &              1557.31,1569.87,1580.06,1591.07, 1603.14,1616.05,
     &          1627.63, 1639.06, 1652.94, 1664.91, 1677.84, 1690.74, 
     &             1703.56,1716.58,1729.50, 1742.92,1756.97, 1771.03, 
     &            1783.97, 1796.90,1810.39, 1824.44,1838.46, 1852.55, 
     &           1867.09, 1880.66, 1897.51,1914.24, 1929.95, 1945.67, 
     &               1964.18,1972.03/)


!    Ni
      TconNi =(/ 1315.67, 1323.99,1333.24,1343.43,1353.31,1363.35,
     &      1373.31,1383.34,1393.43,1403.31, 1412.86,1421.18,1432.30,
     &      1443.43,1454.55,1465.97,1478.76,1491.77,1504.48, 1515.78,
     &      1529.72,1540.84,1554.77, 1568.23,1580.96,1593.81,1607.69,
     &      1622.14,1635.56,1650.51,1666.23,1681.67,1697.25,1713.65,
     &      1728.42,1746.81,1766.68,1786.23,1800.17,1817.54,1838.73,
     &      1857.11,1878.31, 1896.69,1917.47,1936.26,1957.45,1973.04,
     &      1994.22, 2018.22, 2028.81/)


!    Fe
       TCONFe =(/1362.14,1373.43, 1384.60,1395.83, 1407.00,1418.26,
     &              1430.00,1442.48, 1455.61, 1469.10, 1482.58,1496.08,
     &              1509.58,1523.08,1536.57,1550.07,1563.57,1577.13,
     &              1590.57,1604.74,1619.69, 1635.41, 1651.59, 1667.75,
     &             1684.03,1700.80, 1718.31, 1736.60, 1755.27, 1773.98,
     &             1792.93, 1812.32,1832.10, 1852.28,1872.54, 1892.90,
     &             1913.24, 1934.27, 1956.41,1978.37, 2008.05, 2030.80,
     &             2051.13, 2081.89, 2103.84, 2132.13,2157.98, 2190.91,
     &              2221.92, 2247.77, 2263.48/)
  

!    Ca2SiO4
      TconCa2SiO4 =(/1508.24,1518.14,1527.48,1536.15,1544.81,1556.36,
     &      1565.02, 1573.68, 1585.23, 1593.89, 1604.95,1614.10,1625.65,
     &     1634.31,1645.86, 1656.04, 1666.07, 1677.62,1689.17, 1700.23,
     &      1712.27,1722.50, 1735.37,1746.92,1758.47,1772.69,1784.45,
     &      1796.00, 1810.24, 1821.99,1835.26,1847.97,1862.41,1876.85,
     &       1891.14,1903.65,1919.06,1934.48, 1949.03,1963.46,1980.70,
     &       1995.22,2011.50, 2026.98, 2044.31,2060.60,2078.90,2097.18,
     &         2118.35,2139.54, 2150.11/)


!    CaTiO3
      TconCaTiO3 =(/1600.35,1609.68,1616.74,1626.47,1633.19,
     &     1643.94, 1652.50,1662.08,1672.50,1681.30,1692.50,1703.94,
     &     1712.95,1723.94,1735.02,1744.61,1755.37,1766.81,1778.25,
     &     1788.78,1798.49,1810.72,1821.12,1832.96,1845.32,1855.43,
     &     1867.42,1879.89,1892.37,1904.05,1917.31,1929.79,1941.23,
     &     1952.67,1966.98,1979.67,1992.73,2007.04,2021.35,2035.31,
     &     2047.77,2063.12,2078.59,2096.68,2112.02,2130.14,2144.45,
     &     2161.64,2182.03,2204.63, 2213.67/)


!    Al2O3
      TCONAl2O3 =(/1685.09, 1693.08, 1700.92,1708.84, 1716.78,
     &              1724.79,1732.91, 1741.42,1750.37,1759.75, 
     &              1769.31,1779.06,1789.37,1799.94, 1810.65,1820.94,
     &              1830.99,1841.08,1851.02, 1861.02,1871.26,1881.74,
     &              1892.36,1903.03,1913.69,1924.39,1935.05,1945.95,
     &              1957.45,1969.45,1980.72,1989.12,1995.54, 2002.05,
     &              2010.76,2021.69,2033.01,2044.36, 2055.67,2066.99,
     &              2078.33,2089.65,2100.99,2112.62,2124.88,2137.43,
     &              2150.32,2163.28, 2176.89, 2191.32,2198.76/)




      Tconds(1:51,1)=TconKCl 
      Tconds(1:51,2)=TconZnS
      Tconds(1:51,3)=TconNa2S
      Tconds(1:51,4)=TCONMnS
      Tconds(1:51,5)=TconCr2O3
      Tconds(1:51,6)=TCONSiO2
      Tconds(1:51,7)=TCONMg2SiO4 
      Tconds(1:51,8)=TconVO
      Tconds(1:51,9)=TconNi
      Tconds(1:51,10)=TCONFe
      Tconds(1:51,11)=TconCa2SiO4
      Tconds(1:51,12)=TconCaTiO3
      Tconds(1:51,13)=TCONAl2O3


      G0_OPPR(1,1:50,1:50,1)=KCl_vis_500_gg
      G0_OPPR(2,1:50,1:50,1)=KCl_vis_650_gg
      G0_OPPR(3,1:50,1:50,1)=KCl_vis_800_gg
      G0_OPPR(4,1:50,1:50,1)=KCl_rosselandMean_gg
      G0_OPPR(5,1:50,1:50,1)=KCl_rosselandMean_gg
      
      G0_OPPR(1,1:50,1:50,2)=ZnS_vis_500_gg
      G0_OPPR(2,1:50,1:50,2)=ZnS_vis_650_gg
      G0_OPPR(3,1:50,1:50,2)=ZnS_vis_800_gg
      G0_OPPR(4,1:50,1:50,2)=ZnS_rosselandMean_gg
      G0_OPPR(5,1:50,1:50,2)=ZnS_rosselandMean_gg
      
      G0_OPPR(1,1:50,1:50,3)=Na2S_vis_500_gg
      G0_OPPR(2,1:50,1:50,3)=Na2S_vis_650_gg
      G0_OPPR(3,1:50,1:50,3)=Na2S_vis_800_gg
      G0_OPPR(4,1:50,1:50,3)=Na2S_rosselandMean_gg
      G0_OPPR(5,1:50,1:50,3)=Na2S_rosselandMean_gg
      
      G0_OPPR(1,1:50,1:50,4)=MnS_vis_500_gg
      G0_OPPR(2,1:50,1:50,4)=MnS_vis_650_gg
      G0_OPPR(3,1:50,1:50,4)=MnS_vis_800_gg
      G0_OPPR(4,1:50,1:50,4)=MnS_rosselandMean_gg
      G0_OPPR(5,1:50,1:50,4)=MnS_rosselandMean_gg
      
      G0_OPPR(1,1:50,1:50,5)=Cr_vis_500_gg
      G0_OPPR(2,1:50,1:50,5)=Cr_vis_650_gg
      G0_OPPR(3,1:50,1:50,5)=Cr_vis_800_gg
      G0_OPPR(4,1:50,1:50,5)=Cr_rosselandMean_gg
      G0_OPPR(5,1:50,1:50,5)=Cr_rosselandMean_gg
      
      G0_OPPR(1,1:50,1:50,6)=SiO2_vis_500_gg
      G0_OPPR(2,1:50,1:50,6)=SiO2_vis_650_gg
      G0_OPPR(3,1:50,1:50,6)=SiO2_vis_800_gg
      G0_OPPR(4,1:50,1:50,6)=SiO2_rosselandMean_gg
      G0_OPPR(5,1:50,1:50,6)=SiO2_rosselandMean_gg
      
      G0_OPPR(1,1:50,1:50,7)=Mg2SiO4_vis_500_gg
      G0_OPPR(2,1:50,1:50,7)=Mg2SiO4_vis_650_gg
      G0_OPPR(3,1:50,1:50,7)=Mg2SiO4_vis_800_gg
      G0_OPPR(4,1:50,1:50,7)=Mg2SiO4_rosselandMean_gg
      G0_OPPR(5,1:50,1:50,7)=Mg2SiO4_rosselandMean_gg
      
      G0_OPPR(1,1:50,1:50,8)=VO_vis_500_gg
      G0_OPPR(2,1:50,1:50,8)=VO_vis_650_gg
      G0_OPPR(3,1:50,1:50,8)=VO_vis_800_gg
      G0_OPPR(4,1:50,1:50,8)=VO_rosselandMean_gg
      G0_OPPR(5,1:50,1:50,8)=VO_rosselandMean_gg
      
      G0_OPPR(1,1:50,1:50,9)=C2H2_vis_500_gg
      G0_OPPR(2,1:50,1:50,9)=C2H2_vis_650_gg
      G0_OPPR(3,1:50,1:50,9)=C2H2_vis_800_gg
      G0_OPPR(4,1:50,1:50,9)=C2H2_rosselandMean_gg
      G0_OPPR(5,1:50,1:50,9)=C2H2_rosselandMean_gg
      
      G0_OPPR(1,1:50,1:50,10)=Fe_vis_500_gg
      G0_OPPR(2,1:50,1:50,10)=Fe_vis_650_gg
      G0_OPPR(3,1:50,1:50,10)=Fe_vis_800_gg
      G0_OPPR(4,1:50,1:50,10)=Fe_rosselandMean_gg
      G0_OPPR(5,1:50,1:50,10)=Fe_rosselandMean_gg
      
      G0_OPPR(1,1:50,1:50,11)=Fe2SiO4_vis_500_gg
      G0_OPPR(2,1:50,1:50,11)=Fe2SiO4_vis_650_gg
      G0_OPPR(3,1:50,1:50,11)=Fe2SiO4_vis_800_gg
      G0_OPPR(4,1:50,1:50,11)=Fe2SiO4_rosselandMean_gg
      G0_OPPR(5,1:50,1:50,11)=Fe2SiO4_rosselandMean_gg
      
      G0_OPPR(1,1:50,1:50,12)=CaTiO3_vis_500_gg
      G0_OPPR(2,1:50,1:50,12)=CaTiO3_vis_650_gg
      G0_OPPR(3,1:50,1:50,12)=CaTiO3_vis_800_gg
      G0_OPPR(4,1:50,1:50,12)=CaTiO3_rosselandMean_gg
      G0_OPPR(5,1:50,1:50,12)=CaTiO3_rosselandMean_gg
      
      G0_OPPR(1,1:50,1:50,13)=Al2O3_vis_500_gg
      G0_OPPR(2,1:50,1:50,13)=Al2O3_vis_650_gg
      G0_OPPR(3,1:50,1:50,13)=Al2O3_vis_800_gg
      G0_OPPR(4,1:50,1:50,13)=Al2O3_rosselandMean_gg
      G0_OPPR(5,1:50,1:50,13)=Al2O3_rosselandMean_gg
      
      PI0_OPPR(1,1:50,1:50,1)=KCl_vis_500_qscat
      PI0_OPPR(2,1:50,1:50,1)=KCl_vis_650_qscat
      PI0_OPPR(3,1:50,1:50,1)=KCl_vis_800_qscat
      PI0_OPPR(4,1:50,1:50,1)=KCl_rosselandMean_qscat
      PI0_OPPR(5,1:50,1:50,1)=KCl_rosselandMean_qscat
      
      PI0_OPPR(1,1:50,1:50,2)=ZnS_vis_500_qscat
      PI0_OPPR(2,1:50,1:50,2)=ZnS_vis_650_qscat
      PI0_OPPR(3,1:50,1:50,2)=ZnS_vis_800_qscat
      PI0_OPPR(4,1:50,1:50,2)=ZnS_rosselandMean_qscat
      PI0_OPPR(5,1:50,1:50,2)=ZnS_rosselandMean_qscat
      
      PI0_OPPR(1,1:50,1:50,3)=Na2S_vis_500_qscat
      PI0_OPPR(2,1:50,1:50,3)=Na2S_vis_650_qscat
      PI0_OPPR(3,1:50,1:50,3)=Na2S_vis_800_qscat
      PI0_OPPR(4,1:50,1:50,3)=Na2S_rosselandMean_qscat
      PI0_OPPR(5,1:50,1:50,3)=Na2S_rosselandMean_qscat
      
      PI0_OPPR(1,1:50,1:50,4)=MnS_vis_500_qscat
      PI0_OPPR(2,1:50,1:50,4)=MnS_vis_650_qscat
      PI0_OPPR(3,1:50,1:50,4)=MnS_vis_800_qscat
      PI0_OPPR(4,1:50,1:50,4)=MnS_rosselandMean_qscat
      PI0_OPPR(5,1:50,1:50,4)=MnS_rosselandMean_qscat
      
      PI0_OPPR(1,1:50,1:50,5)=Cr_vis_500_qscat
      PI0_OPPR(2,1:50,1:50,5)=Cr_vis_650_qscat
      PI0_OPPR(3,1:50,1:50,5)=Cr_vis_800_qscat
      PI0_OPPR(4,1:50,1:50,5)=Cr_rosselandMean_qscat
      PI0_OPPR(5,1:50,1:50,5)=Cr_rosselandMean_qscat
      
      PI0_OPPR(1,1:50,1:50,6)=SiO2_vis_500_qscat
      PI0_OPPR(2,1:50,1:50,6)=SiO2_vis_650_qscat
      PI0_OPPR(3,1:50,1:50,6)=SiO2_vis_800_qscat
      PI0_OPPR(4,1:50,1:50,6)=SiO2_rosselandMean_qscat
      PI0_OPPR(5,1:50,1:50,6)=SiO2_rosselandMean_qscat
      
      PI0_OPPR(1,1:50,1:50,7)=Mg2SiO4_vis_500_qscat
      PI0_OPPR(2,1:50,1:50,7)=Mg2SiO4_vis_650_qscat
      PI0_OPPR(3,1:50,1:50,7)=Mg2SiO4_vis_800_qscat
      PI0_OPPR(4,1:50,1:50,7)=Mg2SiO4_rosselandMean_qscat
      PI0_OPPR(5,1:50,1:50,7)=Mg2SiO4_rosselandMean_qscat
      
      PI0_OPPR(1,1:50,1:50,8)=VO_vis_500_qscat
      PI0_OPPR(2,1:50,1:50,8)=VO_vis_650_qscat
      PI0_OPPR(3,1:50,1:50,8)=VO_vis_800_qscat
      PI0_OPPR(4,1:50,1:50,8)=VO_rosselandMean_qscat
      PI0_OPPR(5,1:50,1:50,8)=VO_rosselandMean_qscat
      
      PI0_OPPR(1,1:50,1:50,9)=C2H2_vis_500_qscat
      PI0_OPPR(2,1:50,1:50,9)=C2H2_vis_650_qscat
      PI0_OPPR(3,1:50,1:50,9)=C2H2_vis_800_qscat
      PI0_OPPR(4,1:50,1:50,9)=C2H2_rosselandMean_qscat
      PI0_OPPR(5,1:50,1:50,9)=C2H2_rosselandMean_qscat
      
      PI0_OPPR(1,1:50,1:50,10)=Fe_vis_500_qscat
      PI0_OPPR(2,1:50,1:50,10)=Fe_vis_650_qscat
      PI0_OPPR(3,1:50,1:50,10)=Fe_vis_800_qscat
      PI0_OPPR(4,1:50,1:50,10)=Fe_rosselandMean_qscat
      PI0_OPPR(5,1:50,1:50,10)=Fe_rosselandMean_qscat
      
      PI0_OPPR(1,1:50,1:50,11)=Fe2SiO4_vis_500_qscat
      PI0_OPPR(2,1:50,1:50,11)=Fe2SiO4_vis_650_qscat
      PI0_OPPR(3,1:50,1:50,11)=Fe2SiO4_vis_800_qscat
      PI0_OPPR(4,1:50,1:50,11)=Fe2SiO4_rosselandMean_qscat
      PI0_OPPR(5,1:50,1:50,11)=Fe2SiO4_rosselandMean_qscat
      
      PI0_OPPR(1,1:50,1:50,12)=CaTiO3_vis_500_qscat
      PI0_OPPR(2,1:50,1:50,12)=CaTiO3_vis_650_qscat
      PI0_OPPR(3,1:50,1:50,12)=CaTiO3_vis_800_qscat
      PI0_OPPR(4,1:50,1:50,12)=CaTiO3_rosselandMean_qscat
      PI0_OPPR(5,1:50,1:50,12)=CaTiO3_rosselandMean_qscat
      
      PI0_OPPR(1,1:50,1:50,13)=Al2O3_vis_500_qscat
      PI0_OPPR(2,1:50,1:50,13)=Al2O3_vis_650_qscat
      PI0_OPPR(3,1:50,1:50,13)=Al2O3_vis_800_qscat
      PI0_OPPR(4,1:50,1:50,13)=Al2O3_rosselandMean_qscat
      PI0_OPPR(5,1:50,1:50,13)=Al2O3_rosselandMean_qscat
      
      QE_OPPR(1,1:50,1:50,1)=KCl_vis_500_qext
      QE_OPPR(2,1:50,1:50,1)=KCl_vis_650_qext
      QE_OPPR(3,1:50,1:50,1)=KCl_vis_800_qext
      QE_OPPR(4,1:50,1:50,1)=KCl_rosselandMean_qext
      QE_OPPR(5,1:50,1:50,1)=KCl_rosselandMean_qext
      
      QE_OPPR(1,1:50,1:50,2)=ZnS_vis_500_qext
      QE_OPPR(2,1:50,1:50,2)=ZnS_vis_650_qext
      QE_OPPR(3,1:50,1:50,2)=ZnS_vis_800_qext
      QE_OPPR(4,1:50,1:50,2)=ZnS_rosselandMean_qext
      QE_OPPR(5,1:50,1:50,2)=ZnS_rosselandMean_qext
      
      QE_OPPR(1,1:50,1:50,3)=Na2S_vis_500_qext
      QE_OPPR(2,1:50,1:50,3)=Na2S_vis_650_qext
      QE_OPPR(3,1:50,1:50,3)=Na2S_vis_800_qext
      QE_OPPR(4,1:50,1:50,3)=Na2S_rosselandMean_qext
      QE_OPPR(5,1:50,1:50,3)=Na2S_rosselandMean_qext
      
      QE_OPPR(1,1:50,1:50,4)=MnS_vis_500_qext
      QE_OPPR(2,1:50,1:50,4)=MnS_vis_650_qext
      QE_OPPR(3,1:50,1:50,4)=MnS_vis_800_qext
      QE_OPPR(4,1:50,1:50,4)=MnS_rosselandMean_qext
      QE_OPPR(5,1:50,1:50,4)=MnS_rosselandMean_qext
      
      QE_OPPR(1,1:50,1:50,5)=Cr_vis_500_qext
      QE_OPPR(2,1:50,1:50,5)=Cr_vis_650_qext
      QE_OPPR(3,1:50,1:50,5)=Cr_vis_800_qext
      QE_OPPR(4,1:50,1:50,5)=Cr_rosselandMean_qext
      QE_OPPR(5,1:50,1:50,5)=Cr_rosselandMean_qext
      
      QE_OPPR(1,1:50,1:50,6)=SiO2_vis_500_qext
      QE_OPPR(2,1:50,1:50,6)=SiO2_vis_650_qext
      QE_OPPR(3,1:50,1:50,6)=SiO2_vis_800_qext
      QE_OPPR(4,1:50,1:50,6)=SiO2_rosselandMean_qext
      QE_OPPR(5,1:50,1:50,6)=SiO2_rosselandMean_qext
      
      QE_OPPR(1,1:50,1:50,7)=Mg2SiO4_vis_500_qext
      QE_OPPR(2,1:50,1:50,7)=Mg2SiO4_vis_650_qext
      QE_OPPR(3,1:50,1:50,7)=Mg2SiO4_vis_800_qext
      QE_OPPR(4,1:50,1:50,7)=Mg2SiO4_rosselandMean_qext
      QE_OPPR(5,1:50,1:50,7)=Mg2SiO4_rosselandMean_qext
      
      QE_OPPR(1,1:50,1:50,8)=VO_vis_500_qext
      QE_OPPR(2,1:50,1:50,8)=VO_vis_650_qext
      QE_OPPR(3,1:50,1:50,8)=VO_vis_800_qext
      QE_OPPR(4,1:50,1:50,8)=VO_rosselandMean_qext
      QE_OPPR(5,1:50,1:50,8)=VO_rosselandMean_qext
      
      QE_OPPR(1,1:50,1:50,9)=C2H2_vis_500_qext
      QE_OPPR(2,1:50,1:50,9)=C2H2_vis_650_qext
      QE_OPPR(3,1:50,1:50,9)=C2H2_vis_800_qext
      QE_OPPR(4,1:50,1:50,9)=C2H2_rosselandMean_qext
      QE_OPPR(5,1:50,1:50,9)=C2H2_rosselandMean_qext
      
      QE_OPPR(1,1:50,1:50,10)=Fe_vis_500_qext
      QE_OPPR(2,1:50,1:50,10)=Fe_vis_650_qext
      QE_OPPR(3,1:50,1:50,10)=Fe_vis_800_qext
      QE_OPPR(4,1:50,1:50,10)=Fe_rosselandMean_qext
      QE_OPPR(5,1:50,1:50,10)=Fe_rosselandMean_qext
      
      QE_OPPR(1,1:50,1:50,11)=Fe2SiO4_vis_500_qext
      QE_OPPR(2,1:50,1:50,11)=Fe2SiO4_vis_650_qext
      QE_OPPR(3,1:50,1:50,11)=Fe2SiO4_vis_800_qext
      QE_OPPR(4,1:50,1:50,11)=Fe2SiO4_rosselandMean_qext
      QE_OPPR(5,1:50,1:50,11)=Fe2SiO4_rosselandMean_qext
      
      QE_OPPR(1,1:50,1:50,12)=CaTiO3_vis_500_qext
      QE_OPPR(2,1:50,1:50,12)=CaTiO3_vis_650_qext
      QE_OPPR(3,1:50,1:50,12)=CaTiO3_vis_800_qext
      QE_OPPR(4,1:50,1:50,12)=CaTiO3_rosselandMean_qext
      QE_OPPR(5,1:50,1:50,12)=CaTiO3_rosselandMean_qext
      
      QE_OPPR(1,1:50,1:50,13)=Al2O3_vis_500_qext
      QE_OPPR(2,1:50,1:50,13)=Al2O3_vis_650_qext
      QE_OPPR(3,1:50,1:50,13)=Al2O3_vis_800_qext
      QE_OPPR(4,1:50,1:50,13)=Al2O3_rosselandMean_qext
      QE_OPPR(5,1:50,1:50,13)=Al2O3_rosselandMean_qext


      END SUBROUTINE get_cloud_scattering_properties