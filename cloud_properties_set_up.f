      SUBROUTINE get_cloud_scattering_properties_wrapper
          include 'rcommons.h'
          call get_cloud_scattering_properties(NCLOUDS, NLAYER, NVERT, NIRP, NSOLP, GASCON, METALLICITY)
      END SUBROUTINE get_cloud_scattering_properties_wrapper


      SUBROUTINE get_cloud_scattering_properties(NCLOUDS, NLAYER, NVERT, NIRP, NSOLP, GASCON, METALLICITY)
          implicit none
          integer :: J, L, K, NL, NCLOUDS, NLAYER, NVERT, NIRP, NSOLP
          real :: GAS_CONSTANT_R, GASCON, METALLICITY

          character (len = 40) :: haze_type


          ! Define all the arrays

          ! HAZE ARRAYS ARE DIFFERENT THAN THE OTHER ONES
          real, dimension(50, 100)  :: HAZE_RosselandMean_tau_per_bar, HAZE_RosselandMean_pi0, HAZE_RosselandMean_gg
          real, dimension(50, 100)  :: HAZE_PlanckMean_tau_per_bar, HAZE_PlanckMean_pi0, HAZE_PlanckMean_gg
          real, dimension(500, 100) :: HAZE_wav_tau_per_bar, HAZE_wav_pi0, HAZE_wav_gg
          real, dimension(100)      :: haze_pressure_array_pascals

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

          REAL tcon_1X_MET_KCl(NLAYER)
          REAL tcon_1X_MET_ZnS(NLAYER)
          REAL tcon_1X_MET_Na2S(NLAYER)
          REAL tcon_1X_MET_MnS(NLAYER)
          REAL tcon_1X_MET_Cr(NLAYER)
          REAL tcon_1X_MET_SiO2(NLAYER)
          REAL tcon_1X_MET_Mg2SiO4(NLAYER)
          REAL tcon_1X_MET_VO(NLAYER)
          REAL tcon_1X_MET_Ni(NLAYER)
          REAL tcon_1X_MET_Fe(NLAYER)
          REAL tcon_1X_MET_CaSiO4(NLAYER)
          REAL tcon_1X_MET_CaTiO3(NLAYER)
          REAL tcon_1X_MET_Al2O3(NLAYER)


          REAL tcon_100X_MET_KCl(NLAYER)
          REAL tcon_100X_MET_ZnS(NLAYER)
          REAL tcon_100X_MET_Na2S(NLAYER)
          REAL tcon_100X_MET_MnS(NLAYER)
          REAL tcon_100X_MET_Cr(NLAYER)
          REAL tcon_100X_MET_SiO2(NLAYER)
          REAL tcon_100X_MET_Mg2SiO4(NLAYER)
          REAL tcon_100X_MET_VO(NLAYER)
          REAL tcon_100X_MET_Ni(NLAYER)
          REAL tcon_100X_MET_Fe(NLAYER)
          REAL tcon_100X_MET_CaSiO4(NLAYER)
          REAL tcon_100X_MET_CaTiO3(NLAYER)
          REAL tcon_100X_MET_Al2O3(NLAYER)


          REAL tcon_300X_MET_KCl(NLAYER)
          REAL tcon_300X_MET_ZnS(NLAYER)
          REAL tcon_300X_MET_Na2S(NLAYER)
          REAL tcon_300X_MET_MnS(NLAYER)
          REAL tcon_300X_MET_Cr(NLAYER)
          REAL tcon_300X_MET_SiO2(NLAYER)
          REAL tcon_300X_MET_Mg2SiO4(NLAYER)
          REAL tcon_300X_MET_VO(NLAYER)
          REAL tcon_300X_MET_Ni(NLAYER)
          REAL tcon_300X_MET_Fe(NLAYER)
          REAL tcon_300X_MET_CaSiO4(NLAYER)
          REAL tcon_300X_MET_CaTiO3(NLAYER)
          REAL tcon_300X_MET_Al2O3(NLAYER)

          REAL CORFACT(51)
          REAL TCONDS(3, 51, 13)

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

          haze_type = 'soot'
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

      tcon_1X_MET_Al2O3=(/1685.09,1693.08,1700.92,1708.84,1716.78,1724.79,1732.91,1741.42,1750.37,1759.75,
     &                     1769.31,1779.06,1789.37,1799.94,1810.65,1820.94,1830.99,1841.08,1851.02,1861.02,
     &                     1871.26,1881.74,1892.36,1903.03,1913.69,1924.39,1935.05,1945.95,1957.45,1969.45,
     &                     1980.72,1989.12,1995.54,2002.05,2010.76,2021.69,2033.01,2044.36,2055.67,2066.99,
     &                     2078.33,2089.65,2100.99,2112.62,2124.88,2137.43,2150.32,2163.28,2176.89,2191.32,2198.76/)

      tcon_1X_MET_CaSiO4=(/1508.24,1518.14,1527.48,1536.15,1544.81,1556.36,1565.02,1573.68,1585.23,
     &                       1593.89,1604.95,1614.10,1625.65,1634.31,1645.86,1656.04,1666.07,1677.62,
     &                       1689.17,1700.23,1712.27,1722.50,1735.37,1746.92,1758.47,1772.69,1784.45,
     &                       1796.00,1810.24,1821.99,1835.26,1847.97,1862.41,1876.85,1891.14,1903.65,
     &                       1919.06,1934.48,1949.03,1963.46,1980.70,1995.22,2011.50,2026.98,2044.31,
     &                       2060.60,2078.90,2097.18,2118.35,2139.54,2150.11/)

      tcon_1X_MET_Cr=(/1180,1189,1197,1205,1213,1221,1229,1237,1246,1254,1262,1270,1279,1289,1300,1311,1321,
     & 1332,1342,1353,1364,1374,1385,1396,1406,1417,1428,1439,1453,1466,1480,1493,1507,1520,1534,1548,1561,
     & 1575,1588,1602,1617,1633,1649,1665,1681,1697,1713,1729,1745,1761,1791/)

      tcon_1X_MET_Fe=(/1362.14,1373.43,1384.60,1395.83,1407.00,1418.26,1430.00,1442.48,1455.61,1469.10,1482.58,
     & 1496.08,1509.58,1523.08,1536.57,1550.07,1563.57,1577.13,1590.57,1604.74,1619.69,1635.41,1651.59,1667.75,
     & 1684.03,1700.80,1718.31,1736.60,1755.27,1773.98,1792.93,1812.32,1832.10,1852.28,1872.54,1892.90,1913.24,
     & 1934.27,1956.41,1978.37,2008.05,2030.80,2051.13,2081.89,2103.84,2132.13,2157.98,2190.91,2221.92,2247.77,2263.48/)

      tcon_1X_MET_KCl=(/617.032,621.573,626.038,630.552,635.053,639.555,644.050,648.556,653.049,657.552,662.043,
     & 666.609,671.436,676.532,681.879,687.233,692.462,697.665,702.916,708.306,713.767,719.366,725.024,730.775,
     & 736.460,742.266,748.065,753.932,759.779,765.571,771.346,777.201,783.301,789.715,796.379,803.117,809.863,
     & 816.737,823.798,831.052,838.426,845.980,853.873,862.074,870.494,878.913,887.351,895.768,904.198,912.867,917.385/)

      tcon_1X_MET_Mg2SiO4=(/1370.00,1378.31,1386.62,1395.03,1403.56,1412.29,1421.17,1430.18,1439.17,1448.16,1457.24,
     & 1466.52,1475.94,1485.57,1495.23,1505.09,1515.04,1525.21,1535.33,1545.80,1556.37,1567.12,1578.02,1589.13,1600.29,
     & 1611.67,1623.20,1634.84,1646.61,1658.58,1670.74,1683.03,1695.59,1708.44,1721.29,1734.41,1747.86,1761.64,1775.79,
     & 1789.95,1804.36,1819.11,1834.34,1850.41,1867.40,1885.24,1903.85,1923.11,1943.09,1963.67,1974.17/)

      tcon_1X_MET_MnS=(/1113.43,1120.07,1126.64,1133.21,1139.9,1146.46,1153.15,1159.65,1166.25,1173.01,1180.09,
     & 1187.73,1194.69,1202.12,1209.26,1216.50,1223.92,1231.43,1239.19,1247.17,1255.39,1263.78,1272.32,1281.40,
     & 1289.69,1297.76,1306.00,1315.24,1324.18,1333.27,1342.44,1351.39,1360.54,1369.99,1379.72,1389.42,1399.22,
     & 1409.04,1418.99,1428.77,1438.60,1449.11,1459.19,1469.78,1481.06,1492.70,1504.21,1515.49,
     & 1527.84,1540.17,1545.90/)

      tcon_1X_MET_Na2S=(/776.182,781.957,787.715,793.430,799.135,804.761,810.359,815.868,821.297,
     & 826.731,832.157,837.761,843.767,850.179,856.911,863.688,870.458,877.165,884.006,890.930,
     & 897.948,905.112,912.216,919.374,926.464,933.659,940.930,948.253,955.902,963.944,972.372,
     & 980.941,989.400,998.113,1007.19,1016.45,1025.51,1035.47,1044.52,1054.07,1063.67,1073.49,
     & 1083.32,1093.15,1103.16,1113.41,1123.84,1134.53,1145.92,1158.12,1164.38/)

      tcon_1X_MET_Ni=(/1315.67,1323.99,1333.24,1343.43,1353.31,1363.35,1373.31,1383.34,1393.43,
     & 1403.31,1412.86,1421.18,1432.30,1443.43,1454.55,1465.97,1478.76,1491.77,1504.48,1515.78,
     & 1529.72,1540.84,1554.77,1568.23,1580.96,1593.81,1607.69,1622.14,1635.56,1650.51,1666.23,
     & 1681.67,1697.25,1713.65,1728.42,1746.81,1766.68,1786.23,1800.17,1817.54,1838.73,1857.11,
     & 1878.31,1896.69,1917.47,1936.26,1957.45,1973.04,1994.22,2018.22,2028.81/)

      tcon_1X_MET_SiO2=(/1334.63,1342.58,1350.30,1358.48,1366.64,1374.85,1383.15,1391.59,1400.10,
     & 1408.68,1417.25,1425.87,1434.53,1443.14,1451.71,1460.28,1468.90,1477.44,1486.12,1494.77,
     & 1503.91,1513.72,1524.07,1534.57,1544.98,1555.25,1565.35,1575.40,1585.44,1595.55,1606.04,
     & 1617.06,1628.55,1640.25,1651.98,1664.07,1676.82,1689.97,1703.53,1717.16,1731.36,1746.01,
     & 1761.10,1776.94,1793.70,1811.36,1829.60,1848.37,1867.54,1887.00,1896.46/)

      tcon_1X_MET_CaTiO3=(/1595,1604,1612,1621,1630,1638,1649,1659,1669,1679,1689,1699,1709,1719,
     & 1730,1740,1750,1760,1770,1780,1790,1800,1812,1825,1838,1850,1863,1876,1889,1902,1914,1927,
     & 1940,1953,1966,1979,1994,2009,2023,2038,2053,2068,2083,2098,2113,2128,2143,2158,2173,2187,
     & 2200/)

      tcon_1X_MET_VO=(/1363.07,1371.53,1379.98,1389.36,1399.69,1409.65,1419.40,1429.94,1439.11,1450.21,
     & 1458.82,1470.07,1481.33,1492.58,1503.61,1513.72,1526.34,1536.53,1547.03,1557.31,1569.87,1580.06,
     & 1591.07,1603.14,1616.05,1627.63,1639.06,1652.94,1664.91,1677.84,1690.74,1703.56,1716.58,1729.50,
     & 1742.92,1756.97,1771.03,1783.97,1796.90,1810.39,1824.44,1838.46,1852.55,1867.09,1880.66,1897.51,
     & 1914.24,1929.95,1945.67,1964.18,1972.03/)

      tcon_1X_MET_ZnS=(/708.296,712.010,716.704,719.507,722.309,727.329,730.718,736.323,739.126,744.731,
     & 747.534,753.139,755.942,761.547,764.350,769.955,772.758,778.363,783.969,788.079,792.377,797.982,
     & 803.587,806.390,811.995,815.183,823.206,828.812,834.364,840.022,845.628,851.233,856.839,862.444,
     & 868.049,873.655,879.260,884.865,890.471,896.076,901.682,907.287,915.695,921.897,929.708,936.547,
     & 943.722,949.327,955.732,963.063,966.143/)


      tcon_100X_MET_Al2O3=(/2008.0,2021.553,2040.998,2054.791,2074.059,2088.464,2107.995,2123.472,
     & 2143.832,2160.112,2180.654,2198.012,2219.303,2237.854,2259.873,2278.311,2298.424,2317.606,
     & 2339.447,2360.099,2382.107,2403.775,2426.41,2449.387,2472.649,2496.96,2520.821,2546.476,
     & 2570.889,2597.782,2622.287,2650.216,2675.258,2704.599,2730.131,2760.918,2786.946,2819.38,
     & 2846.748,2882.587,2916.348,2958.405,2976.504,2999.0,2999.0,2999.0,2999.0,2999.0,2999.0,2999.0/)

      tcon_100X_MET_CaSiO4=(/1758.0,1771.951,1791.969,1806.543,1826.945,1842.532,1863.738,
     & 1880.146,1901.607,1919.068,1941.232,1959.913,1982.801,2002.418,2025.486,2046.162,
     & 2069.885,2092.086,2116.975,2140.705,2166.215,2191.475,2218.046,2244.952,2272.09,
     & 2300.214,2327.416,2356.78,2384.948,2416.25,2445.379,2478.527,2508.123,2542.757,
     & 2572.768,2608.922,2639.36,2676.977,2707.387,2746.163,2777.36,2818.139,2849.285,
     & 2891.169,2922.669,2964.971,2979.242,2999.0,2999.0,2999.0/)

      tcon_100X_MET_Cr=(/1442.0,1451.168,1464.322,1474.405,1488.573,1498.903,1512.854,1523.764,
     & 1538.071,1549.606,1564.202,1576.709,1592.145,1605.185,1620.388,1633.932,1649.403,1663.646,
     & 1679.392,1694.574,1711.08,1727.225,1743.955,1760.845,1777.806,1795.63,1813.288,1832.357,
     & 1850.666,1870.873,1889.368,1910.51,1929.633,1952.094,1971.803,1995.614,2015.907,2040.843,
     & 2060.391,2085.607,2107.402,2136.244,2160.656,2194.919,2233.88,2287.967,2331.596,2393.455,
     & 2436.423,2499.0/)

      tcon_100X_MET_Fe=(/1451.0,1464.553,1483.998,1498.543,1518.945,1534.888,1556.652,1573.484,
     & 1595.494,1613.387,1636.092,1655.213,1678.634,1699.547,1724.712,1747.219,1773.005,1797.32,
     & 1824.748,1851.082,1879.593,1907.962,1937.978,1968.712,2000.212,2033.544,2066.95,2103.062,
     & 2137.803,2176.893,2214.344,2257.414,2297.026,2343.883,2385.988,2437.232,2482.228,2538.564,
     & 2587.219,2649.628,2701.766,2770.44,2826.418,2900.294,2942.985,2999.0,2999.0,2999.0,2999.0,
     & 2999.0/)

      tcon_100X_MET_KCl=(/651.0,656.182,663.617,668.891,676.258,681.985,689.798,695.719,703.423,
     & 709.781,717.89,724.805,733.321,740.798,749.71,757.495,766.262,774.44,783.583,792.361,801.864,
     & 811.25,821.091,831.107,841.284,851.942,862.441,873.749,884.547,896.671,908.23,921.28,932.663,
     & 946.242,958.784,973.909,986.702,1002.647,1016.114,1033.263,1046.939,1064.958,1079.689,1099.584,
     & 1115.334,1137.109,1153.419,1176.598,1194.268,1220.0/)

      tcon_100X_MET_Mg2SiO4=(/1599.0,1611.756,1630.057,1643.791,1663.059,1677.108,1696.081,
     & 1711.472,1731.832,1748.112,1768.654,1785.712,1806.471,1824.854,1846.873,1866.633,1889.325,
     & 1910.1,1932.957,1955.017,1979.025,2002.838,2027.932,2053.365,2079.05,2106.035,2132.76,
     & 2161.617,2189.316,2220.25,2249.379,2282.79,2313.296,2349.22,2381.023,2419.523,2452.609,
     & 2493.66,2527.545,2571.011,2607.337,2655.034,2692.915,2744.044,2784.248,2839.829,2881.418,
     & 2940.442,2983.008,3045.0/)

      tcon_100X_MET_MnS=(/1235.0,1245.364,1260.234,1270.781,1285.516,1296.971,1312.596,1324.438,
     & 1339.846,1352.243,1367.92,1380.709,1396.145,1409.749,1426.001,1440.196,1456.183,1470.893,
     & 1487.147,1502.574,1519.08,1535.438,1552.66,1569.845,1586.806,1604.809,1622.944,1642.357,
     & 1660.666,1680.873,1699.368,1720.641,1740.22,1759.621,1765.892,1772.0,1772.0,1773.964,
     & 1783.956,1798.96,1822.892,1854.468,1880.563,1915.918,1944.932,1985.274,2018.709,2066.278,
     & 2104.026,2159.0/)

      tcon_100X_MET_Na2S=(/872.0,879.175,889.47,897.148,907.916,916.478,928.197,937.079,948.634,
     & 958.012,969.905,979.907,992.149,1003.056,1016.162,1027.61,1040.503,1052.66,1066.374,1079.656,
     & 1094.161,1108.375,1123.137,1138.258,1153.766,1170.093,1186.318,1203.868,1220.769,1239.579,
     & 1257.149,1277.379,1296.047,1317.978,1337.239,1360.614,1380.907,1406.184,1427.471,1454.607,
     & 1476.402,1504.909,1526.795,1556.21,1578.177,1608.47,1630.08,1660.748,1682.834,1715.0/)

      tcon_100X_MET_Ni=(/1381.0,1392.56,1409.146,1421.662,1439.231,1452.395,1470.253,1484.124,
     & 1502.282,1517.156,1536.076,1552.111,1571.806,1589.007,1609.453,1627.84,1648.984,1668.606,
     & 1690.447,1711.328,1733.836,1756.2,1779.819,1803.974,1828.689,1854.856,1881.103,1909.291,
     & 1936.051,1966.103,1994.77,2027.659,2057.71,2093.104,2124.459,2162.624,2196.151,2237.917,
     & 2273.104,2318.223,2355.831,2405.314,2445.299,2499.294,2541.985,2601.161,2647.644,2713.686,
     & 2763.481,2836.0/)

      tcon_100X_MET_SiO2=(/1507.0,1517.763,1533.205,1544.534,1560.402,1572.683,1589.424,1602.45,
     & 1619.508,1632.881,1649.639,1663.91,1681.475,1696.878,1715.227,1731.518,1750.084,1767.373,
     & 1786.674,1804.951,1824.458,1843.712,1863.887,1884.41,1905.248,1927.243,1949.195,1972.661,
     & 1994.726,2019.488,2043.069,2069.822,2093.499,2121.484,2146.567,2176.918,2202.946,2235.123,
     & 2261.188,2294.668,2322.873,2359.916,2389.378,2429.169,2460.669,2504.246,2537.273,2584.238,
     & 2620.781,2674.0/)

      tcon_100X_MET_CaTiO3=(/1802.0,1815.553,1834.998,1849.167,1869.002,1884.176,1904.823,1920.809,
     & 1941.719,1958.749,1980.373,1998.312,2020.136,2039.136,2061.679,2081.898,2105.105,2126.346,
     & 2149.711,2172.017,2196.025,2182.225,2120.226,2092.583,2116.329,2141.677,2167.447,2195.78,
     & 2223.948,2255.25,2284.379,2317.659,2347.71,2383.22,2415.023,2453.523,2486.609,2527.66,2561.545,
     & 2604.94,2640.839,2687.978,2725.438,2775.919,2814.88,2868.801,2909.983,2966.483,2951.223,2929.0/)

      tcon_100X_MET_VO=(/1361.0,1370.965,1385.263,1395.781,1410.516,1421.615,1436.682,1448.438,
     & 1463.846,1476.562,1492.78,1506.009,1521.978,1536.032,1552.808,1567.725,1584.743,1600.633,
     & 1618.411,1635.263,1653.269,1671.075,1689.773,1708.823,1728.207,1748.347,1767.913,1789.172,
     & 1809.829,1832.754,1854.023,1878.297,1900.152,1926.484,1951.567,1982.019,2008.487,2041.55,
     & 2069.788,2105.809,2134.869,2173.028,2203.332,2244.127,2275.213,2318.191,2350.402,2396.074,
     & 2427.799,2474.0/)

      tcon_100X_MET_ZnS=(/763.0,767.783,774.646,779.891,787.258,792.629,799.884,805.719,813.423,
     & 819.462,827.031,833.504,841.489,848.516,856.904,864.23,872.482,879.947,888.073,895.902,
     & 904.405,912.825,921.682,930.716,939.924,949.584,959.128,969.423,979.282,990.23,1000.402,
     & 1012.149,1023.076,1035.779,1046.529,1059.708,1071.619,1086.306,1098.035,1112.98,1124.946,
     & 1140.678,1153.305,1170.251,1182.685,1199.887,1212.935,1231.476,1245.531,1266.0/)


      tcon_300X_MET_Al2O3=(/2103.0,2118.546,2140.851,2157.048,2179.717,2196.6,2219.48,2237.495,
     & 2261.156,2279.705,2302.951,2322.213,2345.634,2365.983,2390.099,2411.691,2436.445,2459.333,
     & 2484.73,2508.935,2534.944,2560.688,2587.751,2614.952,2642.09,2670.393,2698.072,2727.78,
     & 2755.948,2787.25,2816.379,2849.659,2879.71,2914.873,2945.331,2977.408,2987.995,2999.0,
     & 2999.0,2999.0,2999.0,2999.0,2999.0,2999.0,2999.0,2999.0,2999.0,2999.0,2999.0,2999.0/)

      tcon_300X_MET_CaSiO4=(/1855.0,1870.546,1892.851,1909.424,1932.66,1949.956,1973.394,1991.832,
     & 2016.043,2035.662,2060.529,2081.414,2106.964,2129.676,2156.938,2181.013,2208.346,2233.32,
     & 2260.748,2287.312,2316.322,2344.538,2373.569,2403.125,2433.171,2464.648,2495.667,2528.921,
     & 2560.376,2595.131,2627.034,2663.315,2695.642,2733.336,2765.586,2804.423,2837.068,2877.404,2909.986,
     & 2949.748,2972.398,2999.0,2999.0,2999.0,2999.0,2999.0,2999.0,2999.0,2999.0,2999.0/)

      tcon_300X_MET_Cr=(/1523.0,1534.161,1550.175,1561.91,1578.345,1591.039,1608.339,1619.764,
     & 1634.071,1646.562,1662.78,1676.61,1693.643,1708.596,1726.421,1742.254,1760.304,1776.88,
     & 1795.165,1812.492,1830.999,1849.5,1869.182,1889.018,1908.888,1929.705,1950.226,1972.335,
     & 1993.461,2016.901,2038.632,2063.428,2085.739,2112.02,2135.312,2163.116,2185.614,2213.867,
     & 2238.629,2270.243,2295.885,2329.524,2356.04,2391.793,2419.563,2458.191,2490.402,2536.21,
     & 2571.95,2624.0/)

      tcon_300X_MET_Fe=(/1494.0,1508.749,1529.91,1545.296,1566.831,1583.6,1606.48,1624.495,
     & 1648.156,1667.343,1691.67,1712.414,1737.964,1760.394,1787.131,1811.541,1839.906,1865.813,
     & 1894.257,1922.23,1953.241,1984.025,2016.501,2049.886,2084.293,2120.619,2156.888,2196.203,
     & 2234.23,2277.067,2318.218,2365.594,2409.305,2461.041,2507.625,2564.436,2614.726,2677.588,
     & 2731.456,2800.264,2856.248,2926.307,2958.716,2999.0,2999.0,2999.0,2999.0,2999.0,2999.0,
     & 2999.0/)

      tcon_300X_MET_KCl=(/667.0,672.581,680.588,686.267,694.201,699.985,707.798,714.393,723.198,
     & 729.781,737.89,745.105,754.154,761.798,770.71,779.023,788.822,797.687,807.337,816.59,
     & 826.594,836.462,846.796,857.303,867.964,879.301,890.754,903.075,914.812,927.817,939.839,
     & 953.805,967.009,982.473,995.911,1012.21,1026.327,1043.989,1059.193,1078.617,1094.429,
     & 1115.182,1131.596,1153.793,1171.615,1196.276,1215.032,1241.666,1261.343,1290.0/)

      tcon_300X_MET_Mg2SiO4=(/1688.0,1702.35,1722.939,1737.92,1758.888,1775.244,1797.566,
     & 1814.821,1837.381,1855.705,1878.951,1898.814,1923.299,1944.265,1968.905,1991.219,
     & 2017.005,2041.073,2067.993,2093.623,2121.133,2148.538,2177.569,2207.125,2237.171,
     & 2268.827,2300.324,2334.41,2367.273,2403.865,2438.08,2477.364,2513.334,2555.609,
     & 2592.787,2637.727,2676.107,2724.112,2765.382,2818.284,2862.302,2919.594,2961.683,
     & 3018.795,3066.459,3132.3,3180.821,3249.673,3299.066,3371.0/)

      tcon_300X_MET_MnS=(/1274.0,1285.56,1302.146,1313.91,1330.345,1343.039,1360.339,1373.787,
     & 1391.395,1405.837,1424.217,1439.21,1457.308,1473.16,1492.034,1508.782,1527.864,1545.373,
     & 1564.674,1583.181,1603.188,1622.712,1642.887,1663.214,1683.568,1704.885,1725.883,1746.053,
     & 1760.607,1772.0,1772.0,1773.181,1777.279,1787.368,1812.003,1841.517,1865.78,1896.294,
     & 1923.228,1957.597,1985.375,2021.972,2051.855,2092.169,2123.669,2167.329,2201.58,2250.319,
     & 2289.272,2346.0/)

      tcon_300X_MET_Na2S=(/909.0,916.972,928.411,937.277,949.744,958.834,971.111,980.753,993.409,
     & 1003.968,1017.483,1028.507,1041.815,1053.903,1068.582,1081.403,1095.843,1109.4,1124.638,
     & 1139.344,1155.35,1171.225,1187.955,1205.041,1222.487,1240.809,1258.944,1278.683,1297.931,
     & 1319.313,1339.195,1362.166,1383.566,1408.789,1431.185,1458.317,1481.697,1510.782,1535.109,
     & 1566.172,1591.387,1624.412,1650.087,1684.585,1710.282,1745.747,1771.435,1807.83,1832.325,1868.0/)

      tcon_300X_MET_Ni=(/1418.0,1430.357,1448.087,1461.039,1479.174,1493.108,1512.081,1526.798,
     & 1546.057,1561.793,1581.795,1598.712,1619.471,1637.572,1659.066,1678.369,1700.544,1721.1,
     & 1743.957,1766.017,1790.025,1813.838,1838.932,1864.756,1891.41,1919.393,1947.072,1977.106,
     & 2006.213,2038.69,2069.207,2104.446,2137.229,2175.799,2209.841,2251.326,2287.941,2333.685,
     & 2372.783,2422.86,2464.314,2518.929,2563.544,2623.795,2671.459,2737.522,2789.305,2862.863,
     & 2917.878,2998.0/)

      tcon_300X_MET_SiO2=(/1582.0,1593.958,1611.116,1624.039,1642.174,1655.751,1674.167,1688.798,
     & 1708.057,1723.474,1742.936,1759.111,1778.806,1796.007,1816.453,1834.84,1855.984,1875.606,
     & 1897.447,1918.558,1941.566,1963.988,1987.114,2010.583,2034.329,2059.318,2084.134,2110.965,
     & 2136.786,2165.516,2192.333,2223.003,2250.777,2283.526,2312.64,2347.721,2377.278,2414.319,
     & 2446.466,2487.516,2520.85,2564.699,2600.054,2647.752,2685.055,2736.69,2776.241,2832.428,
     & 2874.593,2936.0/)

      tcon_300X_MET_CaTiO3=(/1871.0,1885.749,1906.91,1922.672,1944.774,1961.244,1983.566,2001.158,
     & 2024.269,2042.705,2065.951,2085.814,2110.299,2081.592,2013.962,2002.219,2028.005,2052.073,
     & 2078.993,2104.623,2132.133,2159.325,2187.864,2217.125,2247.171,2278.648,2309.667,2343.247,
     & 2375.641,2411.865,2446.08,2485.233,2520.748,2562.494,2599.223,2643.627,2681.565,2729.027,
     & 2769.862,2822.143,2865.306,2909.455,2860.211,2802.21,2834.125,2878.745,2919.112,2975.313,
     & 2984.549,2998.0/)

      tcon_300X_MET_VO=(/1375.0,1385.364,1400.234,1411.158,1426.459,1437.615,1452.682,1464.776,
     & 1480.733,1493.562,1509.78,1523.309,1539.81,1554.314,1571.614,1586.989,1604.524,1620.88,
     & 1639.165,1656.492,1674.999,1693.5,1713.182,1733.018,1752.888,1773.705,1794.226,1816.498,
     & 1838.094,1862.047,1884.241,1909.56,1932.326,1959.252,1983.44,2012.818,2038.404,2070.123,
     & 2096.188,2129.526,2156.877,2192.804,2221.424,2260.002,2289.844,2331.385,2366.451,2416.224,
     & 2452.365,2505.0/)

      tcon_300X_MET_ZnS=(/798.0,803.182,810.617,816.643,825.144,830.985,838.798,845.393,854.198,
     & 860.781,868.89,876.105,885.154,892.798,901.71,909.759,919.042,927.44,936.583,945.361,954.864,
     & 964.25,974.091,984.107,994.284,1004.942,1015.441,1026.912,1038.179,1050.524,1061.621,1074.542,
     & 1086.836,1101.358,1114.347,1130.009,1143.244,1159.903,1174.673,1193.546,1208.931,1229.238,1246.073,
     & 1268.334,1281.597,1299.0,1299.0,1299.0,1299.0,1299.0/)

      tconds(1,1:51,1)=tcon_1X_MET_KCl
      tconds(1,1:51,2)=tcon_1X_MET_ZnS
      tconds(1,1:51,3)=tcon_1X_MET_Na2S
      tconds(1,1:51,4)=tcon_1X_MET_MnS
      tconds(1,1:51,5)=tcon_1X_MET_Cr
      tconds(1,1:51,6)=tcon_1X_MET_SiO2
      tconds(1,1:51,7)=tcon_1X_MET_Mg2SiO4
      tconds(1,1:51,8)=tcon_1X_MET_VO
      tconds(1,1:51,9)=tcon_1X_MET_Ni
      tconds(1,1:51,10)=tcon_1X_MET_Fe
      tconds(1,1:51,11)=tcon_1X_MET_CaSiO4
      tconds(1,1:51,12)=tcon_1X_MET_CaTiO3
      tconds(1,1:51,13)=tcon_1X_MET_Al2O3

      tconds(2,1:51,1)=tcon_100X_MET_KCl
      tconds(2,1:51,2)=tcon_100X_MET_ZnS
      tconds(2,1:51,3)=tcon_100X_MET_Na2S
      tconds(2,1:51,4)=tcon_100X_MET_MnS
      tconds(2,1:51,5)=tcon_100X_MET_Cr
      tconds(2,1:51,6)=tcon_100X_MET_SiO2
      tconds(2,1:51,7)=tcon_100X_MET_Mg2SiO4
      tconds(2,1:51,8)=tcon_100X_MET_VO
      tconds(2,1:51,9)=tcon_100X_MET_Ni
      tconds(2,1:51,10)=tcon_100X_MET_Fe
      tconds(2,1:51,11)=tcon_100X_MET_CaSiO4
      tconds(2,1:51,12)=tcon_100X_MET_CaTiO3
      tconds(2,1:51,13)=tcon_100X_MET_Al2O3

      tconds(3,1:51,1)=tcon_300X_MET_KCl
      tconds(3,1:51,2)=tcon_300X_MET_ZnS
      tconds(3,1:51,3)=tcon_300X_MET_Na2S
      tconds(3,1:51,4)=tcon_300X_MET_MnS
      tconds(3,1:51,5)=tcon_300X_MET_Cr
      tconds(3,1:51,6)=tcon_300X_MET_SiO2
      tconds(3,1:51,7)=tcon_300X_MET_Mg2SiO4
      tconds(3,1:51,8)=tcon_300X_MET_VO
      tconds(3,1:51,9)=tcon_300X_MET_Ni
      tconds(3,1:51,10)=tcon_300X_MET_Fe
      tconds(3,1:51,11)=tcon_300X_MET_CaSiO4
      tconds(3,1:51,12)=tcon_300X_MET_CaTiO3
      tconds(3,1:51,13)=tcon_300X_MET_Al2O3

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
