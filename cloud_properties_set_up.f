      SUBROUTINE get_cloud_scattering_properties_wrapper
          include 'rcommons.h'
          call get_cloud_scattering_properties(NCLOUDS, NLAYER, NVERT, NIRP, NSOLP, GASCON)
      END SUBROUTINE get_cloud_scattering_properties_wrapper


      SUBROUTINE get_cloud_scattering_properties(NCLOUDS, NLAYER, NVERT, NIRP, NSOLP, GASCON)
          implicit none
          integer :: J, L, K, NL, NCLOUDS, NLAYER, NVERT, NIRP, NSOLP
          real :: GAS_CONSTANT_R, GASCON

          character (len = 40) :: haze_type


          ! Define all the arrays

          ! HAZE ARRAYS ARE DIFFERENT THAN THE OTHER ONES
          real, dimension(50, 50) :: HAZE_RosselandMean_tau_per_bar, HAZE_RosselandMean_pi0, HAZE_RosselandMean_gg
          real, dimension(50, 50) :: HAZE_PlanckMean_tau_per_bar, HAZE_PlanckMean_pi0, HAZE_PlanckMean_gg
          real, dimension(50, 50) :: HAZE_wav_tau_per_bar, HAZE_wav_pi0, HAZE_wav_gg
          real, dimension(100)    :: haze_pressure_array_pascals

          ! These are 50 by 50 because that's what the data in CLOUD_DATA is
          ! That can change but use to code from Elsie and Isaac
          real, dimension(50,50) :: KCl_rosselandMean_gg
          real, dimension(50,50) :: KCl_rosselandMean_qext
          real, dimension(50,50) :: KCl_rosselandMean_pi0

          real, dimension(50,50) :: KCl_PlanckMean_gg
          real, dimension(50,50) :: KCl_PlanckMean_qext
          real, dimension(50,50) :: KCl_PlanckMean_pi0

          real, dimension(50,50) :: KCl_wav_gg
          real, dimension(50,50) :: KCl_wav_qext
          real, dimension(50,50) :: KCl_wav_pi0


          real, dimension(50,50) :: ZnS_rosselandMean_gg
          real, dimension(50,50) :: ZnS_rosselandMean_qext
          real, dimension(50,50) :: ZnS_rosselandMean_pi0

          real, dimension(50,50) :: ZnS_PlanckMean_gg
          real, dimension(50,50) :: ZnS_PlanckMean_qext
          real, dimension(50,50) :: ZnS_PlanckMean_pi0

          real, dimension(50,50) :: ZnS_wav_gg
          real, dimension(50,50) :: ZnS_wav_qext
          real, dimension(50,50) :: ZnS_wav_pi0


          real, dimension(50,50) :: Na2S_rosselandMean_gg
          real, dimension(50,50) :: Na2S_rosselandMean_qext
          real, dimension(50,50) :: Na2S_rosselandMean_pi0

          real, dimension(50,50) :: Na2S_PlanckMean_gg
          real, dimension(50,50) :: Na2S_PlanckMean_qext
          real, dimension(50,50) :: Na2S_PlanckMean_pi0

          real, dimension(50,50) :: Na2S_wav_gg
          real, dimension(50,50) :: Na2S_wav_qext
          real, dimension(50,50) :: Na2S_wav_pi0


          real, dimension(50,50) :: MnS_rosselandMean_gg
          real, dimension(50,50) :: MnS_rosselandMean_qext
          real, dimension(50,50) :: MnS_rosselandMean_pi0

          real, dimension(50,50) :: MnS_PlanckMean_gg
          real, dimension(50,50) :: MnS_PlanckMean_qext
          real, dimension(50,50) :: MnS_PlanckMean_pi0

          real, dimension(50,50) :: MnS_wav_gg
          real, dimension(50,50) :: MnS_wav_qext
          real, dimension(50,50) :: MnS_wav_pi0


          real, dimension(50,50) :: Cr_rosselandMean_gg
          real, dimension(50,50) :: Cr_rosselandMean_qext
          real, dimension(50,50) :: Cr_rosselandMean_pi0

          real, dimension(50,50) :: Cr_PlanckMean_gg
          real, dimension(50,50) :: Cr_PlanckMean_qext
          real, dimension(50,50) :: Cr_PlanckMean_pi0

          real, dimension(50,50) :: Cr_wav_gg
          real, dimension(50,50) :: Cr_wav_qext
          real, dimension(50,50) :: Cr_wav_pi0


          real, dimension(50,50) :: SiO2_rosselandMean_gg
          real, dimension(50,50) :: SiO2_rosselandMean_qext
          real, dimension(50,50) :: SiO2_rosselandMean_pi0

          real, dimension(50,50) :: SiO2_PlanckMean_gg
          real, dimension(50,50) :: SiO2_PlanckMean_qext
          real, dimension(50,50) :: SiO2_PlanckMean_pi0

          real, dimension(50,50) :: SiO2_wav_gg
          real, dimension(50,50) :: SiO2_wav_qext
          real, dimension(50,50) :: SiO2_wav_pi0


          real, dimension(50,50) :: Mg2SiO4_rosselandMean_gg
          real, dimension(50,50) :: Mg2SiO4_rosselandMean_qext
          real, dimension(50,50) :: Mg2SiO4_rosselandMean_pi0

          real, dimension(50,50) :: Mg2SiO4_PlanckMean_gg
          real, dimension(50,50) :: Mg2SiO4_PlanckMean_qext
          real, dimension(50,50) :: Mg2SiO4_PlanckMean_pi0

          real, dimension(50,50) :: Mg2SiO4_wav_gg
          real, dimension(50,50) :: Mg2SiO4_wav_qext
          real, dimension(50,50) :: Mg2SiO4_wav_pi0


          real, dimension(50,50) :: VO_rosselandMean_gg
          real, dimension(50,50) :: VO_rosselandMean_qext
          real, dimension(50,50) :: VO_rosselandMean_pi0

          real, dimension(50,50) :: VO_PlanckMean_gg
          real, dimension(50,50) :: VO_PlanckMean_qext
          real, dimension(50,50) :: VO_PlanckMean_pi0

          real, dimension(50,50) :: VO_wav_gg
          real, dimension(50,50) :: VO_wav_qext
          real, dimension(50,50) :: VO_wav_pi0


          real, dimension(50,50) :: Ni_rosselandMean_gg
          real, dimension(50,50) :: Ni_rosselandMean_qext
          real, dimension(50,50) :: Ni_rosselandMean_pi0

          real, dimension(50,50) :: Ni_PlanckMean_gg
          real, dimension(50,50) :: Ni_PlanckMean_qext
          real, dimension(50,50) :: Ni_PlanckMean_pi0

          real, dimension(50,50) :: Ni_wav_gg
          real, dimension(50,50) :: Ni_wav_qext
          real, dimension(50,50) :: Ni_wav_pi0


          real, dimension(50,50) :: Fe_rosselandMean_gg
          real, dimension(50,50) :: Fe_rosselandMean_qext
          real, dimension(50,50) :: Fe_rosselandMean_pi0

          real, dimension(50,50) :: Fe_PlanckMean_gg
          real, dimension(50,50) :: Fe_PlanckMean_qext
          real, dimension(50,50) :: Fe_PlanckMean_pi0

          real, dimension(50,50) :: Fe_wav_gg
          real, dimension(50,50) :: Fe_wav_qext
          real, dimension(50,50) :: Fe_wav_pi0



          real, dimension(50,50) :: CaSiO4_rosselandMean_gg
          real, dimension(50,50) :: CaSiO4_rosselandMean_qext
          real, dimension(50,50) :: CaSiO4_rosselandMean_pi0

          real, dimension(50,50) :: CaSiO4_PlanckMean_gg
          real, dimension(50,50) :: CaSiO4_PlanckMean_qext
          real, dimension(50,50) :: CaSiO4_PlanckMean_pi0

          real, dimension(50,50) :: CaSiO4_wav_gg
          real, dimension(50,50) :: CaSiO4_wav_qext
          real, dimension(50,50) :: CaSiO4_wav_pi0


          real, dimension(50,50) :: CaTiO3_rosselandMean_gg
          real, dimension(50,50) :: CaTiO3_rosselandMean_qext
          real, dimension(50,50) :: CaTiO3_rosselandMean_pi0

          real, dimension(50,50) :: CaTiO3_PlanckMean_gg
          real, dimension(50,50) :: CaTiO3_PlanckMean_qext
          real, dimension(50,50) :: CaTiO3_PlanckMean_pi0

          real, dimension(50,50) :: CaTiO3_wav_gg
          real, dimension(50,50) :: CaTiO3_wav_qext
          real, dimension(50,50) :: CaTiO3_wav_pi0


          real, dimension(50,50) :: Al2O3_rosselandMean_gg
          real, dimension(50,50) :: Al2O3_rosselandMean_qext
          real, dimension(50,50) :: Al2O3_rosselandMean_pi0

          real, dimension(50,50) :: Al2O3_PlanckMean_gg
          real, dimension(50,50) :: Al2O3_PlanckMean_qext
          real, dimension(50,50) :: Al2O3_PlanckMean_pi0

          real, dimension(50,50) :: Al2O3_wav_gg
          real, dimension(50,50) :: Al2O3_wav_qext
          real, dimension(50,50) :: Al2O3_wav_pi0

          ! SET UP THE CONDENSATION CURVES
          ! SHOULD BE MET DEPENDENT EVENTUALLY

          REAL TconKCl(NLAYER)
          REAL TconZnS(NLAYER)
          REAL TconNa2S(NLAYER)
          REAL TCONMnS(NLAYER)
          REAL TCONCr(NLAYER)
          REAL TCONSiO2(NLAYER)
          REAL TCONMg2SiO4(NLAYER)
          REAL TconVO(NLAYER)
          REAL TconNi(NLAYER)
          REAL TCONFe(NLAYER)
          REAL TconCa2SiO4(NLAYER)
          REAL TconCaTiO3(NLAYER)
          REAL TCONAl2O3(NLAYER)

          REAL CORFACT(51)
          REAL TCONDS(51, 13)

          REAL QE_OPPR(5, 50, 50, 13)
          REAL PI0_OPPR(5, 50, 50, 13)
          REAL G0_OPPR(5, 50, 50, 13)

          REAL DENSITY(13)
          REAL FMOLW(13)
          REAL CLOUD_MOLAR_MASSES(13)

          real, dimension(50) :: input_temperature_array
          real, dimension(50) :: input_pressure_array_cgs

          real, dimension(50) :: input_particle_size_array_in_meters
          real, dimension(50) :: particle_size_vs_layer_array_in_meters

      COMMON /CLOUD_PROPERTIES/ TCONDS, QE_OPPR, PI0_OPPR, G0_OPPR,
     &                              DENSITY, FMOLW,
     &                              CORFACT,
     &                              input_particle_size_array_in_meters,
     &                              input_temperature_array,
     &                              particle_size_vs_layer_array_in_meters,
     &                              input_pressure_array_cgs,
     &                              HAZE_RosselandMean_tau_per_bar, HAZE_RosselandMean_pi0, HAZE_RosselandMean_gg,
     &                              HAZE_PlanckMean_tau_per_bar,HAZE_PlanckMean_pi0, HAZE_PlanckMean_gg,
     &                              HAZE_wav_tau_per_bar,HAZE_wav_pi0, HAZE_wav_gg,
     &                              haze_pressure_array_pascals

          haze_type = 'sulfur'
          if (haze_type .eq. 'soot') THEN
              !write(*,*) "Model being run with soot hazes"
              open (1, file='../CLOUD_DATA/haze_soot_Ross_tauperbar.txt')
              open (2, file='../CLOUD_DATA/haze_soot_wav_tauperbar.txt')
              open (3, file='../CLOUD_DATA/haze_soot_Planck_tauperbar.txt')

              read(1,*) HAZE_RosselandMean_tau_per_bar
              read(2,*) HAZE_wav_tau_per_bar
              read(3,*) HAZE_PlanckMean_tau_per_bar

              close(1)
              close(2)
              close(3)

              open (1, file='../CLOUD_DATA/haze_soot_Ross_pi0.txt')
              open (2, file='../CLOUD_DATA/haze_soot_wav_pi0.txt')
              open (3, file='../CLOUD_DATA/haze_soot_Planck_pi0.txt')

              read(1,*) HAZE_RosselandMean_pi0
              read(2,*) HAZE_wav_pi0
              read(3,*) HAZE_PlanckMean_pi0

              close(1)
              close(2)
              close(3)

              open (1, file='../CLOUD_DATA/haze_soot_Ross_gg.txt')
              open (2, file='../CLOUD_DATA/haze_soot_wav_gg.txt')
              open (3, file='../CLOUD_DATA/haze_soot_Planck_gg.txt')

              read(1,*) HAZE_RosselandMean_gg
              read(2,*) HAZE_RosselandMean_gg
              read(3,*) HAZE_PlanckMean_gg

              close(1)
              close(2)
              close(3)
          else if (haze_type .eq. 'sulfur') THEN
              !write(*,*) "Model being run with sulfur hazes"
              open (1, file='../CLOUD_DATA/haze_sulfur_Ross_tauperbar.txt')
              open (2, file='../CLOUD_DATA/haze_sulfur_wav_tauperbar.txt')
              open (3, file='../CLOUD_DATA/haze_sulfur_Planck_tauperbar.txt')

              read(1,*) HAZE_RosselandMean_tau_per_bar
              read(2,*) HAZE_wav_tau_per_bar
              read(3,*) HAZE_PlanckMean_tau_per_bar

              close(1)
              close(2)
              close(3)

              open (1, file='../CLOUD_DATA/haze_sulfur_Ross_pi0.txt')
              open (2, file='../CLOUD_DATA/haze_sulfur_wav_pi0.txt')
              open (3, file='../CLOUD_DATA/haze_sulfur_Planck_pi0.txt')

              read(1,*) HAZE_RosselandMean_pi0
              read(2,*) HAZE_wav_pi0
              read(3,*) HAZE_PlanckMean_pi0

              close(1)
              close(2)
              close(3)

              open (1, file='../CLOUD_DATA/haze_sulfur_Ross_gg.txt')
              open (2, file='../CLOUD_DATA/haze_sulfur_wav_gg.txt')
              open (3, file='../CLOUD_DATA/haze_sulfur_Planck_gg.txt')

              read(1,*) HAZE_RosselandMean_gg
              read(2,*) HAZE_wav_gg
              read(3,*) HAZE_PlanckMean_gg

              close(1)
              close(2)
              close(3)
          else if (haze_type .eq. 'tholin') THEN
              !write(*,*) "Model being run with sulfur hazes"
              open (1, file='../CLOUD_DATA/haze_tholin_Ross_tauperbar.txt')
              open (2, file='../CLOUD_DATA/haze_tholin_wav_tauperbar.txt')
              open (3, file='../CLOUD_DATA/haze_tholin_Planck_tauperbar.txt')

              read(1,*) HAZE_RosselandMean_tau_per_bar
              read(2,*) HAZE_wav_tau_per_bar
              read(3,*) HAZE_PlanckMean_tau_per_bar

              close(1)
              close(2)
              close(3)

              open (1, file='../CLOUD_DATA/haze_tholin_Ross_pi0.txt')
              open (2, file='../CLOUD_DATA/haze_tholin_wav_pi0.txt')
              open (3, file='../CLOUD_DATA/haze_tholin_Planck_pi0.txt')

              read(1,*) HAZE_RosselandMean_pi0
              read(2,*) HAZE_wav_pi0
              read(3,*) HAZE_PlanckMean_pi0

              close(1)
              close(2)
              close(3)

              open (1, file='../CLOUD_DATA/haze_tholin_Ross_gg.txt')
              open (2, file='../CLOUD_DATA/haze_tholin_wav_gg.txt')
              open (3, file='../CLOUD_DATA/haze_tholin_Planck_gg.txt')

              read(1,*) HAZE_RosselandMean_gg
              read(2,*) HAZE_wav_gg
              read(3,*) HAZE_PlanckMean_gg

              close(1)
              close(2)
              close(3)
          else if (haze_type .eq. 'soot-2xpi0') THEN
              !write(*,*) "Model being run with sulfur hazes"
              open (1, file='../CLOUD_DATA/haze_soot-2xpi0_Ross_tauperbar.txt')
              open (2, file='../CLOUD_DATA/haze_soot-2xpi0_wav_tauperbar.txt')
              open (3, file='../CLOUD_DATA/haze_soot-2xpi0_Planck_tauperbar.txt')

              read(1,*) HAZE_RosselandMean_tau_per_bar
              read(2,*) HAZE_wav_tau_per_bar
              read(3,*) HAZE_PlanckMean_tau_per_bar

              close(1)
              close(2)
              close(3)

              open (1, file='../CLOUD_DATA/haze_soot-2xpi0_Ross_pi0.txt')
              open (2, file='../CLOUD_DATA/haze_soot-2xpi0_wav_pi0.txt')
              open (3, file='../CLOUD_DATA/haze_soot-2xpi0_Planck_pi0.txt')

              read(1,*) HAZE_RosselandMean_pi0
              read(2,*) HAZE_wav_pi0
              read(3,*) HAZE_PlanckMean_pi0

              close(1)
              close(2)
              close(3)

              open (1, file='../CLOUD_DATA/haze_soot-2xpi0_Ross_gg.txt')
              open (2, file='../CLOUD_DATA/haze_soot-2xpi0_wav_gg.txt')
              open (3, file='../CLOUD_DATA/haze_soot-2xpi0_Planck_gg.txt')

              read(1,*) HAZE_RosselandMean_gg
              read(2,*) HAZE_wav_gg
              read(3,*) HAZE_PlanckMean_gg

              close(1)
              close(2)
              close(3)
          else
              write(*,*) "The haze type is being impropertly specified"
          end if




          ! READ IN ALL THE CLOUD FILES
          ! THIS IS A LOT OF FILES

          open (1, file='../CLOUD_DATA/KCl_rosselandMean_gg.txt')
          open (2, file='../CLOUD_DATA/KCl_rosselandMean_qext.txt')
          open (3, file='../CLOUD_DATA/KCl_rosselandMean_pi0.txt')
          open (4, file='../CLOUD_DATA/KCl_PlanckMean_gg.txt')
          open (5, file='../CLOUD_DATA/KCl_PlanckMean_qext.txt')
          open (6, file='../CLOUD_DATA/KCl_PlanckMean_pi0.txt')
          open (7, file='../CLOUD_DATA/KCl_wav_gg.txt')
          open (8, file='../CLOUD_DATA/KCl_wav_qext.txt')
          open (9, file='../CLOUD_DATA/KCl_wav_pi0.txt')
          
          read(1,*) KCl_rosselandMean_gg
          read(2,*) KCl_rosselandMean_qext
          read(3,*) KCl_rosselandMean_pi0
          read(4,*) KCl_PlanckMean_gg
          read(5,*) KCl_PlanckMean_qext
          read(6,*) KCl_PlanckMean_pi0
          read(7,*) KCl_wav_gg
          read(8,*) KCl_wav_qext
          read(9,*) KCl_wav_pi0
          
          close(1)
          close(2)
          close(3)
          close(4)
          close(5)
          close(6)
          close(7)
          close(8)
          close(9)


          open (1, file='../CLOUD_DATA/ZnS_rosselandMean_gg.txt')
          open (2, file='../CLOUD_DATA/ZnS_rosselandMean_qext.txt')
          open (3, file='../CLOUD_DATA/ZnS_rosselandMean_pi0.txt')
          open (4, file='../CLOUD_DATA/ZnS_wav_gg.txt')
          open (5, file='../CLOUD_DATA/ZnS_wav_qext.txt')
          open (6, file='../CLOUD_DATA/ZnS_wav_pi0.txt')
          open (7, file='../CLOUD_DATA/ZnS_PlanckMean_gg.txt')
          open (8, file='../CLOUD_DATA/ZnS_PlanckMean_qext.txt')
          open (9, file='../CLOUD_DATA/ZnS_PlanckMean_pi0.txt')
          
          read(1,*) ZnS_rosselandMean_gg
          read(2,*) ZnS_rosselandMean_qext
          read(3,*) ZnS_rosselandMean_pi0
          read(4,*) ZnS_wav_gg
          read(5,*) ZnS_wav_qext
          read(6,*) ZnS_wav_pi0
          read(7,*) ZnS_PlanckMean_gg
          read(8,*) ZnS_PlanckMean_qext
          read(9,*) ZnS_PlanckMean_pi0
          
          close(1)
          close(2)
          close(3)
          close(4)
          close(5)
          close(6)
          close(7)
          close(8)
          close(9)

          open (1, file='../CLOUD_DATA/Na2S_rosselandMean_gg.txt')
          open (2, file='../CLOUD_DATA/Na2S_rosselandMean_qext.txt')
          open (3, file='../CLOUD_DATA/Na2S_rosselandMean_pi0.txt')
          open (4, file='../CLOUD_DATA/Na2S_wav_gg.txt')
          open (5, file='../CLOUD_DATA/Na2S_wav_qext.txt')
          open (6, file='../CLOUD_DATA/Na2S_wav_pi0.txt')
          open (7, file='../CLOUD_DATA/Na2S_PlanckMean_gg.txt')
          open (8, file='../CLOUD_DATA/Na2S_PlanckMean_qext.txt')
          open (9, file='../CLOUD_DATA/Na2S_PlanckMean_pi0.txt')
          
          read(1,*) Na2S_rosselandMean_gg
          read(2,*) Na2S_rosselandMean_qext
          read(3,*) Na2S_rosselandMean_pi0
          read(4,*) Na2S_wav_gg
          read(5,*) Na2S_wav_qext
          read(6,*) Na2S_wav_pi0
          read(7,*) Na2S_PlanckMean_gg
          read(8,*) Na2S_PlanckMean_qext
          read(9,*) Na2S_PlanckMean_pi0
          
          close(1)
          close(2)
          close(3)
          close(4)
          close(5)
          close(6)
          close(7)
          close(8)
          close(9)

          open (1, file='../CLOUD_DATA/MnS_rosselandMean_gg.txt')
          open (2, file='../CLOUD_DATA/MnS_rosselandMean_qext.txt')
          open (3, file='../CLOUD_DATA/MnS_rosselandMean_pi0.txt')
          open (4, file='../CLOUD_DATA/MnS_wav_gg.txt')
          open (5, file='../CLOUD_DATA/MnS_wav_qext.txt')
          open (6, file='../CLOUD_DATA/MnS_wav_pi0.txt')
          open (7, file='../CLOUD_DATA/MnS_PlanckMean_gg.txt')
          open (8, file='../CLOUD_DATA/MnS_PlanckMean_qext.txt')
          open (9, file='../CLOUD_DATA/MnS_PlanckMean_pi0.txt')
          
          read(1,*) MnS_rosselandMean_gg
          read(2,*) MnS_rosselandMean_qext
          read(3,*) MnS_rosselandMean_pi0
          read(4,*) MnS_wav_gg
          read(5,*) MnS_wav_qext
          read(6,*) MnS_wav_pi0
          read(7,*) MnS_PlanckMean_gg
          read(8,*) MnS_PlanckMean_qext
          read(9,*) MnS_PlanckMean_pi0
          
          close(1)
          close(2)
          close(3)
          close(4)
          close(5)
          close(6)
          close(7)
          close(8)
          close(9)

          open (1, file='../CLOUD_DATA/Cr_rosselandMean_gg.txt')
          open (2, file='../CLOUD_DATA/Cr_rosselandMean_qext.txt')
          open (3, file='../CLOUD_DATA/Cr_rosselandMean_pi0.txt')
          open (4, file='../CLOUD_DATA/Cr_wav_gg.txt')
          open (5, file='../CLOUD_DATA/Cr_wav_qext.txt')
          open (6, file='../CLOUD_DATA/Cr_wav_pi0.txt')
          open (7, file='../CLOUD_DATA/Cr_PlanckMean_gg.txt')
          open (8, file='../CLOUD_DATA/Cr_PlanckMean_qext.txt')
          open (9, file='../CLOUD_DATA/Cr_PlanckMean_pi0.txt')
          
          read(1,*) Cr_rosselandMean_gg
          read(2,*) Cr_rosselandMean_qext
          read(3,*) Cr_rosselandMean_pi0
          read(4,*) Cr_wav_gg
          read(5,*) Cr_wav_qext
          read(6,*) Cr_wav_pi0
          read(7,*) Cr_PlanckMean_gg
          read(8,*) Cr_PlanckMean_qext
          read(9,*) Cr_PlanckMean_pi0
          
          close(1)
          close(2)
          close(3)
          close(4)
          close(5)
          close(6)
          close(7)
          close(8)
          close(9)

          open (1, file='../CLOUD_DATA/SiO2_rosselandMean_gg.txt')
          open (2, file='../CLOUD_DATA/SiO2_rosselandMean_qext.txt')
          open (3, file='../CLOUD_DATA/SiO2_rosselandMean_pi0.txt')
          open (4, file='../CLOUD_DATA/SiO2_wav_gg.txt')
          open (5, file='../CLOUD_DATA/SiO2_wav_qext.txt')
          open (6, file='../CLOUD_DATA/SiO2_wav_pi0.txt')
          open (7, file='../CLOUD_DATA/SiO2_PlanckMean_gg.txt')
          open (8, file='../CLOUD_DATA/SiO2_PlanckMean_qext.txt')
          open (9, file='../CLOUD_DATA/SiO2_PlanckMean_pi0.txt')
          
          read(1,*) SiO2_rosselandMean_gg
          read(2,*) SiO2_rosselandMean_qext
          read(3,*) SiO2_rosselandMean_pi0
          read(4,*) SiO2_wav_gg
          read(5,*) SiO2_wav_qext
          read(6,*) SiO2_wav_pi0
          read(7,*) SiO2_PlanckMean_gg
          read(8,*) SiO2_PlanckMean_qext
          read(9,*) SiO2_PlanckMean_pi0
          
          close(1)
          close(2)
          close(3)
          close(4)
          close(5)
          close(6)
          close(7)
          close(8)
          close(9)

          open (1, file='../CLOUD_DATA/Mg2SiO4_rosselandMean_gg.txt')
          open (2, file='../CLOUD_DATA/Mg2SiO4_rosselandMean_qext.txt')
          open (3, file='../CLOUD_DATA/Mg2SiO4_rosselandMean_pi0.txt')
          open (4, file='../CLOUD_DATA/Mg2SiO4_wav_gg.txt')
          open (5, file='../CLOUD_DATA/Mg2SiO4_wav_qext.txt')
          open (6, file='../CLOUD_DATA/Mg2SiO4_wav_pi0.txt')
          open (7, file='../CLOUD_DATA/Mg2SiO4_PlanckMean_gg.txt')
          open (8, file='../CLOUD_DATA/Mg2SiO4_PlanckMean_qext.txt')
          open (9, file='../CLOUD_DATA/Mg2SiO4_PlanckMean_pi0.txt')
          
          read(1,*) Mg2SiO4_rosselandMean_gg
          read(2,*) Mg2SiO4_rosselandMean_qext
          read(3,*) Mg2SiO4_rosselandMean_pi0
          read(4,*) Mg2SiO4_wav_gg
          read(5,*) Mg2SiO4_wav_qext
          read(6,*) Mg2SiO4_wav_pi0
          read(7,*) Mg2SiO4_PlanckMean_gg
          read(8,*) Mg2SiO4_PlanckMean_qext
          read(9,*) Mg2SiO4_PlanckMean_pi0
          
          close(1)
          close(2)
          close(3)
          close(4)
          close(5)
          close(6)
          close(7)
          close(8)
          close(9)

          open (1, file='../CLOUD_DATA/VO_rosselandMean_gg.txt')
          open (2, file='../CLOUD_DATA/VO_rosselandMean_qext.txt')
          open (3, file='../CLOUD_DATA/VO_rosselandMean_pi0.txt')
          open (4, file='../CLOUD_DATA/VO_wav_gg.txt')
          open (5, file='../CLOUD_DATA/VO_wav_qext.txt')
          open (6, file='../CLOUD_DATA/VO_wav_pi0.txt')
          open (7, file='../CLOUD_DATA/VO_PlanckMean_gg.txt')
          open (8, file='../CLOUD_DATA/VO_PlanckMean_qext.txt')
          open (9, file='../CLOUD_DATA/VO_PlanckMean_pi0.txt')
          
          read(1,*) VO_rosselandMean_gg
          read(2,*) VO_rosselandMean_qext
          read(3,*) VO_rosselandMean_pi0
          read(4,*) VO_wav_gg
          read(5,*) VO_wav_qext
          read(6,*) VO_wav_pi0
          read(7,*) VO_PlanckMean_gg
          read(8,*) VO_PlanckMean_qext
          read(9,*) VO_PlanckMean_pi0
          
          close(1)
          close(2)
          close(3)
          close(4)
          close(5)
          close(6)
          close(7)
          close(8)
          close(9)

          open (1, file='../CLOUD_DATA/Ni_rosselandMean_gg.txt')
          open (2, file='../CLOUD_DATA/Ni_rosselandMean_qext.txt')
          open (3, file='../CLOUD_DATA/Ni_rosselandMean_pi0.txt')
          open (4, file='../CLOUD_DATA/Ni_wav_gg.txt')
          open (5, file='../CLOUD_DATA/Ni_wav_qext.txt')
          open (6, file='../CLOUD_DATA/Ni_wav_pi0.txt')
          open (7, file='../CLOUD_DATA/Ni_PlanckMean_gg.txt')
          open (8, file='../CLOUD_DATA/Ni_PlanckMean_qext.txt')
          open (9, file='../CLOUD_DATA/Ni_PlanckMean_pi0.txt')
          
          read(1,*) Ni_rosselandMean_gg
          read(2,*) Ni_rosselandMean_qext
          read(3,*) Ni_rosselandMean_pi0
          read(4,*) Ni_wav_gg
          read(5,*) Ni_wav_qext
          read(6,*) Ni_wav_pi0
          read(7,*) Ni_PlanckMean_gg
          read(8,*) Ni_PlanckMean_qext
          read(9,*) Ni_PlanckMean_pi0
          
          close(1)
          close(2)
          close(3)
          close(4)
          close(5)
          close(6)
          close(7)
          close(8)
          close(9)

          open (1, file='../CLOUD_DATA/Fe_rosselandMean_gg.txt')
          open (2, file='../CLOUD_DATA/Fe_rosselandMean_qext.txt')
          open (3, file='../CLOUD_DATA/Fe_rosselandMean_pi0.txt')
          open (4, file='../CLOUD_DATA/Fe_wav_gg.txt')
          open (5, file='../CLOUD_DATA/Fe_wav_qext.txt')
          open (6, file='../CLOUD_DATA/Fe_wav_pi0.txt')
          open (7, file='../CLOUD_DATA/Fe_PlanckMean_gg.txt')
          open (8, file='../CLOUD_DATA/Fe_PlanckMean_qext.txt')
          open (9, file='../CLOUD_DATA/Fe_PlanckMean_pi0.txt')
          
          read(1,*) Fe_rosselandMean_gg
          read(2,*) Fe_rosselandMean_qext
          read(3,*) Fe_rosselandMean_pi0
          read(4,*) Fe_wav_gg
          read(5,*) Fe_wav_qext
          read(6,*) Fe_wav_pi0
          read(7,*) Fe_PlanckMean_gg
          read(8,*) Fe_PlanckMean_qext
          read(9,*) Fe_PlanckMean_pi0
          
          close(1)
          close(2)
          close(3)
          close(4)
          close(5)
          close(6)
          close(7)
          close(8)
          close(9)

          open (1, file='../CLOUD_DATA/CaSiO4_rosselandMean_gg.txt')
          open (2, file='../CLOUD_DATA/CaSiO4_rosselandMean_qext.txt')
          open (3, file='../CLOUD_DATA/CaSiO4_rosselandMean_pi0.txt')
          open (4, file='../CLOUD_DATA/CaSiO4_wav_gg.txt')
          open (5, file='../CLOUD_DATA/CaSiO4_wav_qext.txt')
          open (6, file='../CLOUD_DATA/CaSiO4_wav_pi0.txt')
          open (7, file='../CLOUD_DATA/CaSiO4_PlanckMean_gg.txt')
          open (8, file='../CLOUD_DATA/CaSiO4_PlanckMean_qext.txt')
          open (9, file='../CLOUD_DATA/CaSiO4_PlanckMean_pi0.txt')
          
          read(1,*) CaSiO4_rosselandMean_gg
          read(2,*) CaSiO4_rosselandMean_qext
          read(3,*) CaSiO4_rosselandMean_pi0
          read(4,*) CaSiO4_wav_gg
          read(5,*) CaSiO4_wav_qext
          read(6,*) CaSiO4_wav_pi0
          read(7,*) CaSiO4_PlanckMean_gg
          read(8,*) CaSiO4_PlanckMean_qext
          read(9,*) CaSiO4_PlanckMean_pi0
          
          close(1)
          close(2)
          close(3)
          close(4)
          close(5)
          close(6)
          close(7)
          close(8)
          close(9)

          open (1, file='../CLOUD_DATA/CaTiO3_rosselandMean_gg.txt')
          open (2, file='../CLOUD_DATA/CaTiO3_rosselandMean_qext.txt')
          open (3, file='../CLOUD_DATA/CaTiO3_rosselandMean_pi0.txt')
          open (4, file='../CLOUD_DATA/CaTiO3_wav_gg.txt')
          open (5, file='../CLOUD_DATA/CaTiO3_wav_qext.txt')
          open (6, file='../CLOUD_DATA/CaTiO3_wav_pi0.txt')
          open (7, file='../CLOUD_DATA/CaTiO3_PlanckMean_gg.txt')
          open (8, file='../CLOUD_DATA/CaTiO3_PlanckMean_qext.txt')
          open (9, file='../CLOUD_DATA/CaTiO3_PlanckMean_pi0.txt')
          
          read(1,*) CaTiO3_rosselandMean_gg
          read(2,*) CaTiO3_rosselandMean_qext
          read(3,*) CaTiO3_rosselandMean_pi0
          read(4,*) CaTiO3_wav_gg
          read(5,*) CaTiO3_wav_qext
          read(6,*) CaTiO3_wav_pi0
          read(7,*) CaTiO3_PlanckMean_gg
          read(8,*) CaTiO3_PlanckMean_qext
          read(9,*) CaTiO3_PlanckMean_pi0
          
          close(1)
          close(2)
          close(3)
          close(4)
          close(5)
          close(6)
          close(7)
          close(8)
          close(9)

          open (1, file='../CLOUD_DATA/Al2O3_rosselandMean_gg.txt')
          open (2, file='../CLOUD_DATA/Al2O3_rosselandMean_qext.txt')
          open (3, file='../CLOUD_DATA/Al2O3_rosselandMean_pi0.txt')
          open (4, file='../CLOUD_DATA/Al2O3_wav_gg.txt')
          open (5, file='../CLOUD_DATA/Al2O3_wav_qext.txt')
          open (6, file='../CLOUD_DATA/Al2O3_wav_pi0.txt')
          open (7, file='../CLOUD_DATA/Al2O3_PlanckMean_gg.txt')
          open (8, file='../CLOUD_DATA/Al2O3_PlanckMean_qext.txt')
          open (9, file='../CLOUD_DATA/Al2O3_PlanckMean_pi0.txt')
          
          read(1,*) Al2O3_rosselandMean_gg
          read(2,*) Al2O3_rosselandMean_qext
          read(3,*) Al2O3_rosselandMean_pi0
          read(4,*) Al2O3_wav_gg
          read(5,*) Al2O3_wav_qext
          read(6,*) Al2O3_wav_pi0
          read(7,*) Al2O3_PlanckMean_gg
          read(8,*) Al2O3_PlanckMean_qext
          read(9,*) Al2O3_PlanckMean_pi0
          
          close(1)
          close(2)
          close(3)
          close(4)
          close(5)
          close(6)
          close(7)
          close(8)
          close(9)

!          haze_pressure_array_pascals = (/1.000e-01, 1.456e-01, 2.121e-01, 3.089e-01, 4.498e-01, 6.551e-01, 9.541e-01,
!     &                                    1.389e+00, 2.024e+00, 2.947e+00, 4.292e+00, 6.251e+00, 9.103e+00, 1.326e+01,
!     &                                    1.931e+01, 2.812e+01, 4.095e+01, 5.964e+01, 8.685e+01, 1.265e+02, 1.842e+02,
!     &                                    2.683e+02, 3.907e+02, 5.690e+02, 8.286e+02, 1.207e+03, 1.758e+03, 2.560e+03,
!     &                                    3.728e+03, 5.429e+03, 7.906e+03, 1.151e+04, 1.677e+04, 2.442e+04, 3.556e+04,
!     &                                    5.179e+04, 7.543e+04, 1.099e+05, 1.600e+05, 2.330e+05, 3.393e+05, 4.942e+05,
!     &                                    7.197e+05, 1.048e+06, 1.526e+06, 2.223e+06, 3.237e+06, 4.715e+06, 6.866e+06,
!     &                                    1.000e+07/)

          haze_pressure_array_pascals = (/1.000e-01,1.205e-01,1.451e-01,1.748e-01,2.105e-01,2.535e-01,3.054e-01,
     &                                    3.678e-01,4.431e-01,5.337e-01,6.428e-01,7.743e-01,9.326e-01,1.123e+00,
     &                                    1.353e+00,1.630e+00,1.963e+00,2.364e+00,2.848e+00,3.430e+00,4.132e+00,
     &                                    4.977e+00,5.995e+00,7.221e+00,8.697e+00,1.048e+01,1.262e+01,1.520e+01,
     &                                    1.831e+01,2.205e+01,2.656e+01,3.199e+01,3.854e+01,4.642e+01,5.591e+01,
     &                                    6.734e+01,8.111e+01,9.770e+01,1.177e+02,1.417e+02,1.707e+02,2.057e+02,
     &                                    2.477e+02,2.984e+02,3.594e+02,4.329e+02,5.214e+02,6.280e+02,7.565e+02,
     &                                    9.112e+02,1.097e+03,1.322e+03,1.592e+03,1.918e+03,2.310e+03,2.783e+03,
     &                                    3.352e+03,4.037e+03,4.863e+03,5.857e+03,7.055e+03,8.498e+03,1.024e+04,
     &                                    1.233e+04,1.485e+04,1.789e+04,2.154e+04,2.595e+04,3.126e+04,3.765e+04,
     &                                    4.535e+04,5.462e+04,6.579e+04,7.925e+04,9.545e+04,1.150e+05,1.385e+05,
     &                                    1.668e+05,2.009e+05,2.420e+05,2.915e+05,3.511e+05,4.229e+05,5.094e+05,
     &                                    6.136e+05,7.391e+05,8.902e+05,1.072e+06,1.292e+06,1.556e+06,1.874e+06,
     &                                    2.257e+06,2.719e+06,3.275e+06,3.944e+06,4.751e+06,5.722e+06,6.893e+06,
     &                                    8.302e+06,1.000e+07/)


          input_pressure_array_cgs = (/1.0, 1.46, 2.12, 3.09, 4.5, 6.55, 9.54, 13.89, 20.24,
     &    29.47, 42.92, 62.51, 91.03, 132.57, 193.07, 281.18, 409.49, 596.36, 868.51, 1264.86,
     &    1842.07, 2682.7, 3906.94, 5689.87, 8286.43, 12067.93, 17575.11, 25595.48, 37275.94,
     &    54286.75, 79060.43, 115139.54, 167683.29, 244205.31, 355648.03, 517947.47, 754312.01,
     &    1098541.14, 1599858.72, 2329951.81, 3393221.77, 4941713.36, 7196856.73, 10481131.34,
     &    15264179.67, 22229964.83, 32374575.43, 47148663.63, 68664884.5, 100000000.0/)

          input_particle_size_array_in_meters = (/1.000e-07,1.151e-07,1.326e-07,1.526e-07,1.758e-07,
     &                                             2.024e-07,2.330e-07,2.683e-07,3.089e-07,3.556e-07,
     &                                             4.095e-07,4.715e-07,5.429e-07,6.251e-07,7.197e-07,
     &                                             8.286e-07,9.541e-07,1.099e-06,1.265e-06,1.456e-06,
     &                                             1.677e-06,1.931e-06,2.223e-06,2.560e-06,2.947e-06,
     &                                             3.393e-06,3.907e-06,4.498e-06,5.179e-06,5.964e-06,
     &                                             6.866e-06,7.906e-06,9.103e-06,1.048e-05,1.207e-05,
     &                                             1.389e-05,1.600e-05,1.842e-05,2.121e-05,2.442e-05,
     &                                             2.812e-05,3.237e-05,3.728e-05,4.292e-05,4.942e-05,
     &                                             5.690e-05,6.551e-05,7.543e-05,8.685e-05,1.000e-04/)

          input_temperature_array = (/100.0, 179.59183673, 259.18367347, 338.7755102, 418.36734694, 497.95918367,
     &    577.55102041, 657.14285714, 736.73469388, 816.32653061, 895.91836735, 975.51020408, 1055.10204082,
     &    1134.69387755, 1214.28571429, 1293.87755102, 1373.46938776, 1453.06122449, 1532.65306122,
     &    1612.24489796, 1691.83673469, 1771.42857143, 1851.02040816, 1930.6122449, 2010.20408163,
     &    2089.79591837, 2169.3877551, 2248.97959184, 2328.57142857, 2408.16326531, 2487.75510204,
     &    2567.34693878, 2646.93877551, 2726.53061224, 2806.12244898, 2885.71428571, 2965.30612245,
     &    3044.89795918, 3124.48979592, 3204.08163265, 3283.67346939, 3363.26530612, 3442.85714286,
     &    3522.44897959, 3602.04081633, 3681.63265306, 3761.22448980, 3840.81632653, 3920.40816327, 4000.0/)

          particle_size_vs_layer_array_in_meters = (/1e-07, 1e-07, 1e-07, 1e-07, 1e-07, 1e-07, 1e-07, 1e-07,
     &    1e-07, 1e-07, 1e-07, 1e-07, 1e-07, 1e-07, 1e-07, 1e-07, 1e-07, 1e-07, 1e-07, 1e-07, 1e-07, 1e-07,
     &    1e-07, 1e-07, 1e-07, 1.007e-07, 1.048e-07, 1.111e-07, 1.201e-07, 1.333e-07, 1.556e-07, 2.106e-07,
     &    2.655e-07, 3.452e-07, 4.609e-07, 6.302e-07, 8.764e-07, 1.2349e-06, 1.7566e-06, 2.5167e-06,
     &    3.6238e-06, 5.2369e-06, 7.5846e-06, 1.10048e-05, 1.59852e-05, 2.3239e-05, 3.37988e-05,
     &    4.91833e-05, 7.15907e-05, 0.0001042202/)

          DENSITY = (/1.98e3,4.09e3,1.86e3,4.0e3,5.22e3,2.65e3,3.27e3,5.76e3,8.9e3, 7.9e3,3.34e3,3.98e3,3.95e3/)

          ! The molar masses of the different cloud species in grams/mol
          CLOUD_MOLAR_MASSES = (/74.55E-3,    ! KCl
     &                           97.47E-3,   ! ZnS
     &                           78.05E-3,  ! Na2S
     &                           87.00E-3, ! MnS
     &                           52.00E-3,    ! Cr
     &                           60.08E-3,    ! SiO2
     &                           160.95E-3,   ! Mg2Si04
     &                           66.94E-3,  ! VO
     &                           58.69E-3,  ! Ni
     &                           55.85E-3,    ! Fe
     &                           172.23E-3,   ! Ca2Si04
     &                           135.94E-3,   ! CaTiO3
     &                           102.00E-3/)   ! Al2O3


          GAS_CONSTANT_R = 8.314462618 ! This is in SI

          ! Gives the FMOLW in SI.
          DO J = 1, 13
              FMOLW(J) = CLOUD_MOLAR_MASSES(J) / (GAS_CONSTANT_R/GASCON)
          END DO

          ! This is missing that annoying species I can't find
          ! https://arxiv.org/pdf/astro-ph/9807055.pdf

          !  Not Nucleation Limited
          !  1) KCl     || 1.23e-7
          !  2) ZnS     || 4.06e-8
          !  3) Na2S    || 9.35e-7
          !  4) MnS     || 3.11e-7,
          !  5) Cr      || 4.4e-7
          !  6) SiO2    || 3.26e-5
          !  7) Mg2Si04 || 1.745e-5
          !  8) VO      || 9.56e-9
          !  9) Ni      || 1.61e-6
          ! 10) Fe      || 2.94e-5
          ! 11) Ca2Si04 || 1.99e-6
          ! 12) CaTiO3  || 7.83e-8
          ! 13) Al2O3   || 1.385e-6

          ! Nucleation Limited
          ! 1.23e-7,0.0,0.0,0.0,4.4e-7,3.26e-5,1.745e-5,9.56e-9,0.0,0.0,1.99e-6,7.83e-8,1.385e-6,
          !  1) KCl     || 1.23e-7
          !  2) ZnS     || 0.0
          !  3) Na2S    || 0.0
          !  4) MnS     || 0.0
          !  5) Cr203   || 4.4e-7
          !  6) SiO2    || 3.26e-5
          !  7) Mg2Si04 || 1.745e-5
          !  8) VO      || 9.56e-9
          !  9) Ni      || 0.0
          ! 10) Fe      || 0.0
          ! 11) Ca2Si04 || 1.99e-6
          ! 12) CaTiO3  || 7.83e-8
          ! 13) Al2O3   || 1.385e-6


          ! Cloud Free
          ! 0,0,0,0,0,0,0,0,0,0,0,0,0


      CORFACT =   (/0.005,0.018,0.050,0.135,.367,1.000,1.000,1.000,
     &              1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,
     &              1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,
     &              1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,
     &              1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,
     &              1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,
     &              1.000,1.000,1.000/)

!      CORFACT = (/0.0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.40, 0.55, 0.70, 0.85,
!     &            1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,
!     &            1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,
!     &            1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,
!     &            1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000/)


!      CORFACT = (/1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,
!     &            1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,
!     &            1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,
!     &            1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,
!     &            1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000/)

!    KCl
      TconKCl = (/617.032,621.573,626.038,630.552,635.053, 639.555,
     &            644.050,648.556,653.049,657.552,662.043,666.609,671.436,
     &            676.532,681.879,687.233,692.462,697.665,702.916,708.306,
     &            713.767,719.366,725.024,730.775,736.460,742.266,748.065,
     &            753.932,759.779,765.571, 771.346,777.201,783.301,789.715,
     &            796.379,803.117,809.863,816.737,823.798,831.052,838.426,
     &            845.980,853.873,862.074,870.494,878.913,887.351,895.768,
     &            904.198,912.867,917.385/)

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
     &              972.372,980.941,989.400,998.113,1007.19,1016.45,1025.51,
     &              1035.47,1044.52,1054.07,1063.67,1073.49,1083.32,1093.15,
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
!      TconCr2O3 =(/1213.17,1219.05,1224.98,1231.07,1237.15,1243.21,
!     &        1249.26,1255.35,1261.51,1267.83,1274.22,1280.63,1287.04,
!     &       1293.56,1300.89,1309.34,1318.81,1329.04,1339.04,1349.79,
!     &       1361.13,1373.13,1385.33,1397.48,1409.69,1421.78,1434.01,
!     &      1446.16,1458.55,1471.19,1483.84,1496.49,1508.99,1522.14,
!     &      1536.54,1552.17,1568.58,1585.09,1601.61,1618.14,1634.62,
!     &      1651.02,1667.41,1683.80,1700.20,1716.57,1732.89,1749.26,
!     &             1765.73,1783.15,1792.19/)

       TconCr = (/1180, 1189, 1197, 1205, 1213, 1221, 1229, 1237, 1246, 1254,
     &             1262, 1270, 1279, 1289, 1300, 1311, 1321, 1332, 1342, 1353,
     &             1364, 1374, 1385, 1396, 1406, 1417, 1428, 1439, 1453, 1466,
     &             1480, 1493, 1507, 1520, 1534, 1548, 1561, 1575, 1588, 1602,
     &             1617, 1633, 1649, 1665, 1681, 1697, 1713, 1729, 1745, 1761, 1791/)



!    SiO2
      TCONSiO2 = (/1334.63,1342.58,1350.30,1358.48,1366.64,1374.85,
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
!      TconCaTiO3 =(/1600.35,1609.68,1616.74,1626.47,1633.19,
!     &     1643.94, 1652.50,1662.08,1672.50,1681.30,1692.50,1703.94,
!     &     1712.95,1723.94,1735.02,1744.61,1755.37,1766.81,1778.25,
!     &     1788.78,1798.49,1810.72,1821.12,1832.96,1845.32,1855.43,
!     &     1867.42,1879.89,1892.37,1904.05,1917.31,1929.79,1941.23,
!     &     1952.67,1966.98,1979.67,1992.73,2007.04,2021.35,2035.31,
!     &     2047.77,2063.12,2078.59,2096.68,2112.02,2130.14,2144.45,
!     &     2161.64,2182.03,2204.63, 2213.67/)

       TconCaTiO3 = (/1595, 1604, 1612, 1621, 1630, 1638, 1649, 1659, 1669, 1679,
     &                1689, 1699, 1709, 1719, 1730, 1740, 1750, 1760, 1770, 1780,
     &                1790, 1800, 1812, 1825, 1838, 1850, 1863, 1876, 1889, 1902,
     &                1914, 1927, 1940, 1953, 1966, 1979, 1994, 2009, 2023, 2038,
     &                2053, 2068, 2083, 2098, 2113, 2128, 2143, 2158, 2173, 2187, 2200/)


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
      !Tconds(1:51,5)=TconCr2O3
      Tconds(1:51,5)=TconCr
      Tconds(1:51,6)=TCONSiO2
      Tconds(1:51,7)=TCONMg2SiO4 
      Tconds(1:51,8)=TconVO
      Tconds(1:51,9)=TconNi
      Tconds(1:51,10)=TCONFe
      Tconds(1:51,11)=TconCa2SiO4
      Tconds(1:51,12)=TconCaTiO3
      Tconds(1:51,13)=TCONAl2O3


      G0_OPPR(1,1:50,1:50,1)=KCl_wav_gg
      G0_OPPR(2,1:50,1:50,1)=KCl_wav_gg
      G0_OPPR(3,1:50,1:50,1)=KCl_wav_gg
      G0_OPPR(4,1:50,1:50,1)=KCl_PlanckMean_gg
      G0_OPPR(5,1:50,1:50,1)=KCl_rosselandMean_gg
      
      G0_OPPR(1,1:50,1:50,2)=ZnS_wav_gg
      G0_OPPR(2,1:50,1:50,2)=ZnS_wav_gg
      G0_OPPR(3,1:50,1:50,2)=ZnS_wav_gg
      G0_OPPR(4,1:50,1:50,2)=ZnS_PlanckMean_gg
      G0_OPPR(5,1:50,1:50,2)=ZnS_rosselandMean_gg
      
      G0_OPPR(1,1:50,1:50,3)=Na2S_wav_gg
      G0_OPPR(2,1:50,1:50,3)=Na2S_wav_gg
      G0_OPPR(3,1:50,1:50,3)=Na2S_wav_gg
      G0_OPPR(4,1:50,1:50,3)=Na2S_PlanckMean_gg
      G0_OPPR(5,1:50,1:50,3)=Na2S_rosselandMean_gg
      
      G0_OPPR(1,1:50,1:50,4)=MnS_wav_gg
      G0_OPPR(2,1:50,1:50,4)=MnS_wav_gg
      G0_OPPR(3,1:50,1:50,4)=MnS_wav_gg
      G0_OPPR(4,1:50,1:50,4)=MnS_PlanckMean_gg
      G0_OPPR(5,1:50,1:50,4)=MnS_rosselandMean_gg
      
      G0_OPPR(1,1:50,1:50,5)=Cr_wav_gg
      G0_OPPR(2,1:50,1:50,5)=Cr_wav_gg
      G0_OPPR(3,1:50,1:50,5)=Cr_wav_gg
      G0_OPPR(4,1:50,1:50,5)=Cr_PlanckMean_gg
      G0_OPPR(5,1:50,1:50,5)=Cr_rosselandMean_gg

      G0_OPPR(1,1:50,1:50,6)=SiO2_wav_gg
      G0_OPPR(2,1:50,1:50,6)=SiO2_wav_gg
      G0_OPPR(3,1:50,1:50,6)=SiO2_wav_gg
      G0_OPPR(4,1:50,1:50,6)=SiO2_PlanckMean_gg
      G0_OPPR(5,1:50,1:50,6)=SiO2_rosselandMean_gg

      G0_OPPR(1,1:50,1:50,7)=Mg2SiO4_wav_gg
      G0_OPPR(2,1:50,1:50,7)=Mg2SiO4_wav_gg
      G0_OPPR(3,1:50,1:50,7)=Mg2SiO4_wav_gg
      G0_OPPR(4,1:50,1:50,7)=Mg2SiO4_PlanckMean_gg
      G0_OPPR(5,1:50,1:50,7)=Mg2SiO4_rosselandMean_gg
      
      G0_OPPR(1,1:50,1:50,8)=VO_wav_gg
      G0_OPPR(2,1:50,1:50,8)=VO_wav_gg
      G0_OPPR(3,1:50,1:50,8)=VO_wav_gg
      G0_OPPR(4,1:50,1:50,8)=VO_PlanckMean_gg
      G0_OPPR(5,1:50,1:50,8)=VO_rosselandMean_gg
      
      G0_OPPR(1,1:50,1:50,9)=Ni_wav_gg
      G0_OPPR(2,1:50,1:50,9)=Ni_wav_gg
      G0_OPPR(3,1:50,1:50,9)=Ni_wav_gg
      G0_OPPR(4,1:50,1:50,9)=Ni_PlanckMean_gg
      G0_OPPR(5,1:50,1:50,9)=Ni_rosselandMean_gg
      
      G0_OPPR(1,1:50,1:50,10)=Fe_wav_gg
      G0_OPPR(2,1:50,1:50,10)=Fe_wav_gg
      G0_OPPR(3,1:50,1:50,10)=Fe_wav_gg
      G0_OPPR(4,1:50,1:50,10)=Fe_PlanckMean_gg
      G0_OPPR(5,1:50,1:50,10)=Fe_rosselandMean_gg
      
      G0_OPPR(1,1:50,1:50,11)=CaSiO4_wav_gg
      G0_OPPR(2,1:50,1:50,11)=CaSiO4_wav_gg
      G0_OPPR(3,1:50,1:50,11)=CaSiO4_wav_gg
      G0_OPPR(4,1:50,1:50,11)=CaSiO4_PlanckMean_gg
      G0_OPPR(5,1:50,1:50,11)=CaSiO4_rosselandMean_gg
      
      G0_OPPR(1,1:50,1:50,12)=CaTiO3_wav_gg
      G0_OPPR(2,1:50,1:50,12)=CaTiO3_wav_gg
      G0_OPPR(3,1:50,1:50,12)=CaTiO3_wav_gg
      G0_OPPR(4,1:50,1:50,12)=CaTiO3_PlanckMean_gg
      G0_OPPR(5,1:50,1:50,12)=CaTiO3_rosselandMean_gg
      
      G0_OPPR(1,1:50,1:50,13)=Al2O3_wav_gg
      G0_OPPR(2,1:50,1:50,13)=Al2O3_wav_gg
      G0_OPPR(3,1:50,1:50,13)=Al2O3_wav_gg
      G0_OPPR(4,1:50,1:50,13)=Al2O3_PlanckMean_gg
      G0_OPPR(5,1:50,1:50,13)=Al2O3_rosselandMean_gg
      
      PI0_OPPR(1,1:50,1:50,1)=KCl_wav_pi0
      PI0_OPPR(2,1:50,1:50,1)=KCl_wav_pi0
      PI0_OPPR(3,1:50,1:50,1)=KCl_wav_pi0
      PI0_OPPR(4,1:50,1:50,1)=KCl_PlanckMean_pi0
      PI0_OPPR(5,1:50,1:50,1)=KCl_rosselandMean_pi0
      
      PI0_OPPR(1,1:50,1:50,2)=ZnS_wav_pi0
      PI0_OPPR(2,1:50,1:50,2)=ZnS_wav_pi0
      PI0_OPPR(3,1:50,1:50,2)=ZnS_wav_pi0
      PI0_OPPR(4,1:50,1:50,2)=ZnS_PlanckMean_pi0
      PI0_OPPR(5,1:50,1:50,2)=ZnS_rosselandMean_pi0
      
      PI0_OPPR(1,1:50,1:50,3)=Na2S_wav_pi0
      PI0_OPPR(2,1:50,1:50,3)=Na2S_wav_pi0
      PI0_OPPR(3,1:50,1:50,3)=Na2S_wav_pi0
      PI0_OPPR(4,1:50,1:50,3)=Na2S_PlanckMean_pi0
      PI0_OPPR(5,1:50,1:50,3)=Na2S_rosselandMean_pi0
      
      PI0_OPPR(1,1:50,1:50,4)=MnS_wav_pi0
      PI0_OPPR(2,1:50,1:50,4)=MnS_wav_pi0
      PI0_OPPR(3,1:50,1:50,4)=MnS_wav_pi0
      PI0_OPPR(4,1:50,1:50,4)=MnS_PlanckMean_pi0
      PI0_OPPR(5,1:50,1:50,4)=MnS_rosselandMean_pi0
      
      PI0_OPPR(1,1:50,1:50,5)=Cr_wav_pi0
      PI0_OPPR(2,1:50,1:50,5)=Cr_wav_pi0
      PI0_OPPR(3,1:50,1:50,5)=Cr_wav_pi0
      PI0_OPPR(4,1:50,1:50,5)=Cr_PlanckMean_pi0
      PI0_OPPR(5,1:50,1:50,5)=Cr_rosselandMean_pi0
      
      PI0_OPPR(1,1:50,1:50,6)=SiO2_wav_pi0
      PI0_OPPR(2,1:50,1:50,6)=SiO2_wav_pi0
      PI0_OPPR(3,1:50,1:50,6)=SiO2_wav_pi0
      PI0_OPPR(4,1:50,1:50,6)=SiO2_PlanckMean_pi0
      PI0_OPPR(5,1:50,1:50,6)=SiO2_rosselandMean_pi0
      
      PI0_OPPR(1,1:50,1:50,7)=Mg2SiO4_wav_pi0
      PI0_OPPR(2,1:50,1:50,7)=Mg2SiO4_wav_pi0
      PI0_OPPR(3,1:50,1:50,7)=Mg2SiO4_wav_pi0
      PI0_OPPR(4,1:50,1:50,7)=Mg2SiO4_PlanckMean_pi0
      PI0_OPPR(5,1:50,1:50,7)=Mg2SiO4_rosselandMean_pi0
      
      PI0_OPPR(1,1:50,1:50,8)=VO_wav_pi0
      PI0_OPPR(2,1:50,1:50,8)=VO_wav_pi0
      PI0_OPPR(3,1:50,1:50,8)=VO_wav_pi0
      PI0_OPPR(4,1:50,1:50,8)=VO_PlanckMean_pi0
      PI0_OPPR(5,1:50,1:50,8)=VO_rosselandMean_pi0
      
      PI0_OPPR(1,1:50,1:50,9)=Ni_wav_pi0
      PI0_OPPR(2,1:50,1:50,9)=Ni_wav_pi0
      PI0_OPPR(3,1:50,1:50,9)=Ni_wav_pi0
      PI0_OPPR(4,1:50,1:50,9)=Ni_PlanckMean_pi0
      PI0_OPPR(5,1:50,1:50,9)=Ni_rosselandMean_pi0
      
      PI0_OPPR(1,1:50,1:50,10)=Fe_wav_pi0
      PI0_OPPR(2,1:50,1:50,10)=Fe_wav_pi0
      PI0_OPPR(3,1:50,1:50,10)=Fe_wav_pi0
      PI0_OPPR(4,1:50,1:50,10)=Fe_PlanckMean_pi0
      PI0_OPPR(5,1:50,1:50,10)=Fe_rosselandMean_pi0
      
      PI0_OPPR(1,1:50,1:50,11)=CaSiO4_wav_pi0
      PI0_OPPR(2,1:50,1:50,11)=CaSiO4_wav_pi0
      PI0_OPPR(3,1:50,1:50,11)=CaSiO4_wav_pi0
      PI0_OPPR(4,1:50,1:50,11)=CaSiO4_PlanckMean_pi0
      PI0_OPPR(5,1:50,1:50,11)=CaSiO4_rosselandMean_pi0
      
      PI0_OPPR(1,1:50,1:50,12)=CaTiO3_wav_pi0
      PI0_OPPR(2,1:50,1:50,12)=CaTiO3_wav_pi0
      PI0_OPPR(3,1:50,1:50,12)=CaTiO3_wav_pi0
      PI0_OPPR(4,1:50,1:50,12)=CaTiO3_PlanckMean_pi0
      PI0_OPPR(5,1:50,1:50,12)=CaTiO3_rosselandMean_pi0
      
      PI0_OPPR(1,1:50,1:50,13)=Al2O3_wav_pi0
      PI0_OPPR(2,1:50,1:50,13)=Al2O3_wav_pi0
      PI0_OPPR(3,1:50,1:50,13)=Al2O3_wav_pi0
      PI0_OPPR(4,1:50,1:50,13)=Al2O3_PlanckMean_pi0
      PI0_OPPR(5,1:50,1:50,13)=Al2O3_rosselandMean_pi0
      
      QE_OPPR(1,1:50,1:50,1)=KCl_wav_qext
      QE_OPPR(2,1:50,1:50,1)=KCl_wav_qext
      QE_OPPR(3,1:50,1:50,1)=KCl_wav_qext
      QE_OPPR(4,1:50,1:50,1)=KCl_PlanckMean_qext
      QE_OPPR(5,1:50,1:50,1)=KCl_rosselandMean_qext
      
      QE_OPPR(1,1:50,1:50,2)=ZnS_wav_qext
      QE_OPPR(2,1:50,1:50,2)=ZnS_wav_qext
      QE_OPPR(3,1:50,1:50,2)=ZnS_wav_qext
      QE_OPPR(4,1:50,1:50,2)=ZnS_PlanckMean_qext
      QE_OPPR(5,1:50,1:50,2)=ZnS_rosselandMean_qext
      
      QE_OPPR(1,1:50,1:50,3)=Na2S_wav_qext
      QE_OPPR(2,1:50,1:50,3)=Na2S_wav_qext
      QE_OPPR(3,1:50,1:50,3)=Na2S_wav_qext
      QE_OPPR(4,1:50,1:50,3)=Na2S_PlanckMean_qext
      QE_OPPR(5,1:50,1:50,3)=Na2S_rosselandMean_qext
      
      QE_OPPR(1,1:50,1:50,4)=MnS_wav_qext
      QE_OPPR(2,1:50,1:50,4)=MnS_wav_qext
      QE_OPPR(3,1:50,1:50,4)=MnS_wav_qext
      QE_OPPR(4,1:50,1:50,4)=MnS_PlanckMean_qext
      QE_OPPR(5,1:50,1:50,4)=MnS_rosselandMean_qext
      
      QE_OPPR(1,1:50,1:50,5)=Cr_wav_qext
      QE_OPPR(2,1:50,1:50,5)=Cr_wav_qext
      QE_OPPR(3,1:50,1:50,5)=Cr_wav_qext
      QE_OPPR(4,1:50,1:50,5)=Cr_PlanckMean_qext
      QE_OPPR(5,1:50,1:50,5)=Cr_rosselandMean_qext
      
      QE_OPPR(1,1:50,1:50,6)=SiO2_wav_qext
      QE_OPPR(2,1:50,1:50,6)=SiO2_wav_qext
      QE_OPPR(3,1:50,1:50,6)=SiO2_wav_qext
      QE_OPPR(4,1:50,1:50,6)=SiO2_PlanckMean_qext
      QE_OPPR(5,1:50,1:50,6)=SiO2_rosselandMean_qext

      QE_OPPR(1,1:50,1:50,7)=Mg2SiO4_wav_qext
      QE_OPPR(2,1:50,1:50,7)=Mg2SiO4_wav_qext
      QE_OPPR(3,1:50,1:50,7)=Mg2SiO4_wav_qext
      QE_OPPR(4,1:50,1:50,7)=Mg2SiO4_PlanckMean_qext
      QE_OPPR(5,1:50,1:50,7)=Mg2SiO4_rosselandMean_qext
      
      QE_OPPR(1,1:50,1:50,8)=VO_wav_qext
      QE_OPPR(2,1:50,1:50,8)=VO_wav_qext
      QE_OPPR(3,1:50,1:50,8)=VO_wav_qext
      QE_OPPR(4,1:50,1:50,8)=VO_PlanckMean_qext
      QE_OPPR(5,1:50,1:50,8)=VO_rosselandMean_qext
      
      QE_OPPR(1,1:50,1:50,9)=Ni_wav_qext
      QE_OPPR(2,1:50,1:50,9)=Ni_wav_qext
      QE_OPPR(3,1:50,1:50,9)=Ni_wav_qext
      QE_OPPR(4,1:50,1:50,9)=Ni_PlanckMean_qext
      QE_OPPR(5,1:50,1:50,9)=Ni_rosselandMean_qext
      
      QE_OPPR(1,1:50,1:50,10)=Fe_wav_qext
      QE_OPPR(2,1:50,1:50,10)=Fe_wav_qext
      QE_OPPR(3,1:50,1:50,10)=Fe_wav_qext
      QE_OPPR(4,1:50,1:50,10)=Fe_PlanckMean_qext
      QE_OPPR(5,1:50,1:50,10)=Fe_rosselandMean_qext
      
      QE_OPPR(1,1:50,1:50,11)=CaSiO4_wav_qext
      QE_OPPR(2,1:50,1:50,11)=CaSiO4_wav_qext
      QE_OPPR(3,1:50,1:50,11)=CaSiO4_wav_qext
      QE_OPPR(4,1:50,1:50,11)=CaSiO4_PlanckMean_qext
      QE_OPPR(5,1:50,1:50,11)=CaSiO4_rosselandMean_qext
      
      QE_OPPR(1,1:50,1:50,12)=CaTiO3_wav_qext
      QE_OPPR(2,1:50,1:50,12)=CaTiO3_wav_qext
      QE_OPPR(3,1:50,1:50,12)=CaTiO3_wav_qext
      QE_OPPR(4,1:50,1:50,12)=CaTiO3_PlanckMean_qext
      QE_OPPR(5,1:50,1:50,12)=CaTiO3_rosselandMean_qext
      
      QE_OPPR(1,1:50,1:50,13)=Al2O3_wav_qext
      QE_OPPR(2,1:50,1:50,13)=Al2O3_wav_qext
      QE_OPPR(3,1:50,1:50,13)=Al2O3_wav_qext
      QE_OPPR(4,1:50,1:50,13)=Al2O3_PlanckMean_qext
      QE_OPPR(5,1:50,1:50,13)=Al2O3_rosselandMean_qext

      END SUBROUTINE get_cloud_scattering_properties
