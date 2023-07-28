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

          REAL tcon_10X_MET_KCl(NLAYER)
          REAL tcon_10X_MET_ZnS(NLAYER)
          REAL tcon_10X_MET_Na2S(NLAYER)
          REAL tcon_10X_MET_MnS(NLAYER)
          REAL tcon_10X_MET_Cr(NLAYER)
          REAL tcon_10X_MET_SiO2(NLAYER)
          REAL tcon_10X_MET_Mg2SiO4(NLAYER)
          REAL tcon_10X_MET_VO(NLAYER)
          REAL tcon_10X_MET_Ni(NLAYER)
          REAL tcon_10X_MET_Fe(NLAYER)
          REAL tcon_10X_MET_CaSiO4(NLAYER)
          REAL tcon_10X_MET_CaTiO3(NLAYER)
          REAL tcon_10X_MET_Al2O3(NLAYER)

          
          REAL tcon_30X_MET_KCl(NLAYER)
          REAL tcon_30X_MET_ZnS(NLAYER)
          REAL tcon_30X_MET_Na2S(NLAYER)
          REAL tcon_30X_MET_MnS(NLAYER)
          REAL tcon_30X_MET_Cr(NLAYER)
          REAL tcon_30X_MET_SiO2(NLAYER)
          REAL tcon_30X_MET_Mg2SiO4(NLAYER)
          REAL tcon_30X_MET_VO(NLAYER)
          REAL tcon_30X_MET_Ni(NLAYER)
          REAL tcon_30X_MET_Fe(NLAYER)
          REAL tcon_30X_MET_CaSiO4(NLAYER)
          REAL tcon_30X_MET_CaTiO3(NLAYER)
          REAL tcon_30X_MET_Al2O3(NLAYER)
          
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
          REAL TCONDS(5, 51, 13)

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

          haze_type = 'soot-2xpi0'
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
              read(2,*) HAZE_wav_gg
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

      CORFACT =   (/0.005,0.018,0.050,0.135,.367,1.000,1.000,1.000,
     &              1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,
     &              1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,
     &              1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,
     &              1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,
     &              1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,
     &              1.000,1.000,1.000/)




      tcon_1X_MET_Al2O3 = (/1690.51,1694.72,1700.551,1707.694,1715.872,1724.545,
     &                     1733.048,1741.593,1750.45,1759.672,1769.316,1779.323,
     &                     1789.557,1799.923,1810.322,1820.674,1830.914,1841.019,
     &                     1851.073,1861.174,1871.418,1881.812,1892.312,1902.882,
     &                     1913.503,1924.302,1935.507,1946.978,1958.25,1968.79,
     &                     1978.382,1987.252,1995.68,2003.995,2012.612,2022.009,
     &                     2032.487,2043.791,2055.401,2066.868,2078.179,2089.546,
     &                     2101.105,2112.879,2124.909,2138.543,2152.32,2165.116,
     &                     2175.892,2184.253/)

      tcon_10X_MET_Al2O3 = (/1793.847,1796.725,1804.404,1816.022,1830.039,1844.913,
     &                     1859.184,1872.394,1886.823,1902.104,1916.936,1931.316,
     &                     1946.993,1961.781,1977.08,1993.23,2008.805,2023.762,
     &                     2039.948,2055.79,2072.042,2089.898,2107.592,2124.689,
     &                     2142.687,2160.609,2178.375,2198.119,2217.721,2236.435,
     &                     2255.62,2274.924,2292.88,2311.937,2330.166,2347.03,
     &                     2364.239,2382.733,2400.32,2420.796,2441.881,2462.458,
     &                     2483.501,2506.392,2528.479,2558.213,2591.015,2623.469,
     &                     2654.412,2684.302/)

      tcon_30X_MET_Al2O3 = (/1900.965,1904.403,1913.326,1925.906,1940.594,1956.196,
     &                     1971.961,1987.714,2003.79,2020.201,2036.869,2053.747,
     &                     2070.846,2088.062,2105.425,2122.931,2140.457,2158.049,
     &                     2175.904,2194.107,2212.799,2232.048,2251.669,2271.566,
     &                     2291.825,2312.459,2333.537,2355.133,2377.025,2399.025,
     &                     2421.169,2443.489,2466.097,2489.319,2513.37,2538.611,
     &                     2565.689,2594.995,2626.307,2658.954,2691.448,2722.115,
     &                     2749.908,2774.42,2796.023,2815.818,2834.279,2851.268,
     &                     2866.565,2879.943/)

      tcon_100X_MET_Al2O3 = (/2015.396,2024.867,2037.947,2053.505,2071.329,2089.609,
     &                     2107.211,2124.495,2142.764,2160.859,2179.998,2199.137,
     &                     2218.919,2238.39,2258.387,2278.191,2298.43,2318.592,
     &                     2339.257,2360.156,2381.804,2403.908,2426.434,2449.377,
     &                     2472.781,2496.81,2521.159,2546.129,2571.251,2597.063,
     &                     2622.901,2649.548,2676.12,2703.619,2730.938,2759.156,
     &                     2787.518,2819.347,2852.073,2886.935,2920.118,2949.667,
     &                     2973.034,2989.913,2999.258,3002.765,3001.943,3000.877,
     &                     2999.007,2998.993/)

      tcon_300X_MET_Al2O3 = (/2111.441,2122.378,2137.581,2155.697,2176.526,2197.931,
     &                     2218.556,2238.731,2259.812,2280.459,2301.934,2323.073,
     &                     2344.991,2366.781,2389.488,2412.288,2435.916,2459.73,
     &                     2484.316,2509.267,2534.916,2560.981,2587.596,2614.686,
     &                     2642.157,2670.146,2698.375,2727.298,2756.518,2786.767,
     &                     2817.643,2850.906,2884.037,2916.276,2944.833,2968.277,
     &                     2985.261,2996.107,3000.954,3001.916,3000.59,2999.918,
     &                     2999.007,2998.999,2999.002,2999.004,2998.996,2999.012,
     &                     2998.983,2999.017/)

      tcon_1X_MET_CaSiO4 = (/1514.925,1519.874,1526.84,1535.157,1544.884,1554.989,
     &                     1564.885,1574.535,1584.584,1594.474,1604.349,1614.656,
     &                     1624.863,1634.945,1645.438,1655.887,1666.822,1677.633,
     &                     1688.911,1700.196,1711.564,1723.283,1735.152,1747.017,
     &                     1759.288,1771.64,1784.282,1796.713,1809.378,1822.226,
     &                     1835.349,1848.777,1862.366,1876.143,1890.452,1904.654,
     &                     1919.293,1933.96,1948.973,1964.195,1979.834,1995.33,
     &                     2011.153,2027.233,2043.794,2062.726,2082.353,2100.669,
     &                     2116.415,2128.722/)

      tcon_10X_MET_CaSiO4 = (/1578.887,1581.908,1589.958,1602.141,1616.839,1632.441,
     &                     1647.416,1661.288,1676.44,1692.498,1708.115,1723.311,
     &                     1739.946,1755.772,1772.494,1790.53,1808.228,1825.391,
     &                     1843.824,1861.795,1880.192,1900.358,1920.322,1939.629,
     &                     1960.022,1980.437,2000.861,2023.812,2046.877,2069.206,
     &                     2092.358,2116.017,2139.04,2165.381,2192.1,2217.901,
     &                     2244.188,2271.577,2297.426,2327.301,2357.729,2386.968,
     &                     2416.316,2447.646,2476.576,2511.886,2548.538,2583.215,
     &                     2615.753,2648.262/)

      tcon_30X_MET_CaSiO4 = (/1665.649,1669.128,1678.41,1691.634,1707.172,1723.749,
     &                     1740.551,1757.368,1774.519,1791.992,1809.698,1827.63,
     &                     1845.893,1864.492,1883.558,1903.138,1923.086,1943.387,
     &                     1964.167,1985.44,2007.305,2029.808,2052.731,2075.958,
     &                     2099.581,2123.644,2148.29,2173.685,2199.655,2226.05,
     &                     2252.954,2280.421,2308.56,2337.566,2367.197,2397.216,
     &                     2427.785,2459.203,2491.891,2526.496,2563.073,2601.312,
     &                     2640.798,2680.703,2719.863,2757.282,2791.357,2820.521,
     &                     2844.208,2862.722/)

      tcon_100X_MET_CaSiO4 = (/1765.53,1775.32,1789.012,1805.419,1824.284,1843.704,
     &                     1862.528,1881.074,1900.675,1920.11,1940.485,1960.734,
     &                     1981.862,2002.922,2024.806,2046.783,2069.578,2092.626,
     &                     2116.544,2140.917,2166.084,2191.787,2218.065,2244.816,
     &                     2271.984,2299.72,2327.759,2356.559,2385.663,2415.764,
     &                     2446.102,2477.571,2509.04,2541.662,2574.054,2607.587,
     &                     2640.693,2675.196,2709.162,2744.634,2779.493,2817.368,
     &                     2854.573,2892.2,2926.156,2955.22,2976.943,2992.121,
     &                     2999.782,3002.659/)

      tcon_300X_MET_CaSiO4 = (/1863.385,1874.403,1889.789,1908.181,1929.337,1951.119,
     &                     1972.249,1993.055,2014.928,2036.575,2059.464,2082.404,
     &                     2106.53,2130.784,2156.018,2181.292,2207.431,2233.747,
     &                     2260.741,2287.961,2315.917,2344.347,2373.473,2403.242,
     &                     2433.515,2464.512,2495.983,2528.344,2560.932,2594.439,
     &                     2627.944,2662.386,2696.634,2731.978,2767.088,2804.049,
     &                     2840.6,2878.624,2914.204,2945.474,2969.986,2987.769,
     &                     2998.013,3002.699,3002.584,3001.226,2999.006,2999.006,
     &                     2998.987,2999.017/)

      tcon_1X_MET_Cr = (/1186.393,1190.691,1196.637,1203.832,1212.116,1220.834,
     &                     1229.202,1237.312,1245.373,1253.461,1261.8,1270.622,
     &                     1279.823,1289.52,1299.701,1310.227,1321.049,1331.783,
     &                     1342.396,1353.021,1363.623,1374.294,1384.916,1395.392,
     &                     1405.976,1416.859,1428.363,1440.269,1452.844,1465.937,
     &                     1479.324,1493.067,1506.804,1520.397,1533.999,1547.437,
     &                     1560.901,1574.467,1588.422,1602.72,1617.539,1632.922,
     &                     1648.664,1664.767,1680.918,1698.345,1715.341,1730.734,
     &                     1743.478,1753.194/)

      tcon_10X_MET_Cr = (/1289.598,1291.635,1297.073,1305.302,1315.221,1325.734,
     &                     1335.797,1345.073,1355.09,1365.597,1375.727,1385.537,
     &                     1396.322,1406.576,1417.38,1428.998,1440.367,1451.37,
     &                     1463.186,1474.689,1486.408,1499.18,1511.775,1523.936,
     &                     1536.82,1549.771,1562.756,1577.351,1592.137,1606.701,
     &                     1622.181,1638.579,1655.826,1677.8,1701.854,1726.41,
     &                     1751.726,1777.734,1803.083,1834.217,1867.245,1899.867,
     &                     1932.619,1966.64,1997.593,2034.048,2070.954,2105.196,
     &                     2137.083,2169.443/)

      tcon_30X_MET_Cr = (/1365.892,1367.841,1373.879,1382.578,1392.812,1403.733,
     &                     1414.797,1425.854,1437.103,1448.534,1460.087,1471.746,
     &                     1483.567,1495.532,1507.708,1520.119,1532.677,1545.391,
     &                     1558.351,1571.564,1585.082,1598.938,1613.015,1627.269,
     &                     1641.785,1656.611,1671.84,1687.581,1703.753,1720.325,
     &                     1737.479,1755.468,1774.689,1795.508,1817.68,1840.73,
     &                     1864.428,1888.78,1914.175,1941.344,1970.538,2001.888,
     &                     2035.845,2072.825,2113.086,2156.666,2201.684,2245.256,
     &                     2285.32,2321.209/)

      tcon_100X_MET_Cr = (/1446.914,1453.505,1462.724,1473.76,1486.446,1499.499,
     &                     1512.11,1524.454,1537.423,1550.291,1563.865,1577.373,
     &                     1591.447,1605.453,1619.914,1634.294,1649.077,1663.935,
     &                     1679.305,1694.903,1710.923,1727.187,1743.8,1760.728,
     &                     1777.949,1795.625,1813.604,1832.14,1850.901,1870.3,
     &                     1889.813,1910.031,1930.324,1951.553,1972.768,1994.706,
     &                     2016.346,2038.867,2061.099,2083.929,2106.918,2133.77,
     &                     2163.185,2199.278,2238.951,2288.761,2341.741,2392.948,
     &                     2436.733,2470.425/)

      tcon_300X_MET_Cr = (/1529.007,1536.829,1547.997,1561.43,1576.55,1591.76,
     &                     1606.077,1619.805,1634.021,1647.888,1662.484,1677.191,
     &                     1693.059,1709.156,1725.857,1742.607,1759.841,1777.081,
     &                     1794.777,1812.676,1831.09,1849.845,1869.087,1888.793,
     &                     1908.919,1929.578,1950.513,1972.046,1993.844,2016.392,
     &                     2039.145,2062.89,2086.72,2111.316,2135.759,2161.281,
     &                     2186.715,2213.542,2240.311,2268.741,2297.16,2327.352,
     &                     2357.255,2389.361,2421.594,2460.463,2500.658,2539.411,
     &                     2572.713,2599.666/)

      tcon_1X_MET_Fe = (/1369.778,1375.758,1383.971,1394.089,1405.725,1418.161,
     &                     1430.424,1442.838,1455.687,1468.912,1482.393,1495.975,
     &                     1509.542,1523.075,1536.522,1549.918,1563.331,1576.909,
     &                     1590.833,1605.209,1620.06,1635.377,1651.149,1667.399,
     &                     1684.077,1701.165,1718.651,1736.578,1754.976,1773.818,
     &                     1793.015,1812.505,1832.187,1852.049,1872.234,1892.124,
     &                     1913.094,1935.23,1957.295,1980.73,2004.599,2029.262,
     &                     2053.9,2079.152,2105.67,2134.941,2165.117,2193.011,
     &                     2215.838,2234.038/)

      tcon_10X_MET_Fe = (/1359.979,1363.319,1372.201,1385.652,1401.883,1419.125,
     &                     1435.695,1451.078,1467.955,1485.918,1503.461,1520.596,
     &                     1539.371,1557.32,1576.504,1597.439,1618.195,1638.488,
     &                     1660.291,1681.663,1703.907,1728.75,1753.721,1778.135,
     &                     1803.904,1829.719,1855.812,1885.513,1915.692,1945.188,
     &                     1975.88,2007.352,2038.439,2074.78,2112.294,2149.07,
     &                     2186.787,2226.163,2264.021,2309.192,2356.358,2402.651,
     &                     2449.585,2499.56,2546.526,2606.248,2669.924,2731.377,
     &                     2789.497,2846.633/)

      tcon_30X_MET_Fe = (/1407.446,1411.683,1421.583,1435.544,1451.982,1469.59,
     &                     1487.499,1505.465,1523.822,1542.571,1561.65,1581.101,
     &                     1601.085,1621.653,1642.992,1665.159,1687.947,1711.281,
     &                     1735.262,1759.93,1785.48,1812.039,1839.38,1867.369,
     &                     1896.112,1925.691,1956.348,1988.336,2021.432,2055.424,
     &                     2090.431,2126.607,2164.265,2203.816,2245.003,2287.485,
     &                     2331.423,2377.21,2425.577,2477.559,2533.016,2590.941,
     &                     2649.888,2707.988,2763.778,2816.693,2865.097,2907.171,
     &                     2942.382,2970.772/)

      tcon_100X_MET_Fe = (/1458.225,1467.779,1481.279,1497.582,1516.486,1536.106,
     &                     1555.324,1574.369,1594.414,1614.262,1635.142,1656.07,
     &                     1678.147,1700.418,1723.838,1747.64,1772.557,1797.939,
     &                     1824.315,1851.295,1879.319,1908.131,1937.908,1968.695,
     &                     2000.428,2033.32,2067.092,2102.342,2138.501,2176.46,
     &                     2215.291,2256.253,2298.08,2342.402,2387.471,2435.465,
     &                     2484.259,2536.583,2589.831,2648.428,2707.993,2773.229,
     &                     2836.482,2893.726,2939.946,2974.59,2995.433,3005.66,
     &                     3006.104,3003.706/)

      tcon_300X_MET_Fe = (/1501.946,1512.242,1526.7,1544.055,1564.149,1584.966,
     &                     1605.31,1625.495,1646.929,1668.293,1690.792,1713.33,
     &                     1737.159,1761.305,1786.616,1812.172,1838.806,1865.908,
     &                     1894.121,1923.011,1953.045,1984.062,2016.364,2049.912,
     &                     2084.457,2120.239,2157.036,2195.512,2235.034,2276.604,
     &                     2319.233,2364.296,2410.409,2459.431,2509.437,2562.874,
     &                     2617.546,2677.823,2738.868,2802.6,2862.179,2914.164,
     &                     2954.436,2983.047,2998.761,3005.309,3004.269,3002.377,
     &                     2998.996,2999.011/)

      tcon_1X_MET_KCl = (/620.112,622.526,625.857,629.946,634.588,639.433,
     &                     644.046,648.521,652.983,657.466,662.037,666.752,
     &                     671.628,676.656,681.805,687.05,692.345,697.657,
     &                     702.984,708.361,713.817,719.374,725.014,730.723,
     &                     736.483,742.28,748.086,753.881,759.661,765.464,
     &                     771.333,777.331,783.48,789.802,796.306,802.988,
     &                     809.819,816.758,823.805,831.006,838.433,846.109,
     &                     854.023,862.136,870.399,879.515,888.529,896.712,
     &                     903.494,908.658/)

      tcon_10X_MET_KCl = (/616.371,617.636,621.007,626.111,632.265,638.794,
     &                     645.053,650.84,657.129,663.764,670.197,676.452,
     &                     683.324,689.881,696.858,704.434,711.905,719.171,
     &                     726.952,734.522,742.251,750.694,759.036,767.099,
     &                     775.634,784.205,792.831,802.599,812.472,822.068,
     &                     832.011,842.15,852.014,863.29,874.728,885.783,
     &                     897.071,908.872,920.083,933.191,946.665,959.714,
     &                     972.862,986.883,999.914,1016.069,1033.011,1049.17,
     &                     1064.371,1079.474/)

      tcon_30X_MET_KCl = (/634.836,636.235,639.892,645.094,651.201,657.708,
     &                     664.292,670.87,677.58,684.428,691.394,698.492,
     &                     705.767,713.21,720.858,728.715,736.704,744.807,
     &                     753.07,761.504,770.164,779.081,788.177,797.415,
     &                     806.847,816.503,826.46,836.79,847.399,858.198,
     &                     869.215,880.48,892.061,904.06,916.391,928.952,
     &                     941.775,954.899,968.383,982.387,996.843,1011.637,
     &                     1026.829,1042.547,1058.988,1076.381,1094.182,1111.452,
     &                     1127.641,1142.84/)

      tcon_100X_MET_KCl = (/653.809,657.406,662.437,668.463,675.364,682.441,
     &                     689.268,695.987,703.105,710.176,717.652,725.173,
     &                     733.103,741.058,749.31,757.571,766.093,774.665,
     &                     783.505,792.484,801.78,811.299,821.075,831.095,
     &                     841.321,851.81,862.502,873.611,884.895,896.547,
     &                     908.319,920.64,933.102,946.178,959.34,973.242,
     &                     987.244,1002.033,1016.74,1032.244,1047.697,1064.283,
     &                     1080.85,1098.632,1116.36,1137.157,1158.206,1178.101,
     &                     1194.859,1208.199/)

      tcon_300X_MET_KCl = (/670.062,673.97,679.324,685.649,692.966,700.522,
     &                     707.78,714.879,722.414,729.915,737.793,745.596,
     &                     753.736,761.953,770.646,779.437,788.515,797.671,
     &                     807.149,816.749,826.548,836.48,846.714,857.273,
     &                     868.132,879.378,890.872,902.746,914.836,927.478,
     &                     940.323,953.769,967.405,981.752,996.258,1011.629,
     &                     1027.101,1043.512,1060.023,1077.603,1095.176,1113.981,
     &                     1132.746,1152.917,1173.043,1196.564,1220.365,1242.897,
     &                     1261.867,1276.893/)

      tcon_1X_MET_Mg2SiO4 = (/1375.557,1380.021,1386.25,1393.957,1402.809,1412.147,
     &                     1421.17,1430.075,1439.065,1448.144,1457.308,1466.578,
     &                     1475.971,1485.533,1495.234,1505.073,1515.041,1525.147,
     &                     1535.394,1545.813,1556.385,1567.122,1578.023,1589.1,
     &                     1600.309,1611.666,1623.171,1634.818,1646.608,1658.59,
     &                     1670.743,1683.073,1695.593,1708.314,1721.281,1734.5,
     &                     1747.967,1761.661,1775.569,1789.737,1804.216,1819.107,
     &                     1834.551,1850.655,1867.561,1887.06,1906.958,1925.523,
     &                     1941.241,1953.447/)

      tcon_10X_MET_Mg2SiO4 = (/1434.816,1437.589,1444.985,1456.177,1469.679,1484.011,
     &                     1497.765,1510.504,1524.414,1539.151,1553.479,1567.42,
     &                     1582.683,1597.203,1612.54,1629.077,1645.306,1661.056,
     &                     1677.992,1694.536,1711.553,1730.309,1748.956,1767.039,
     &                     1786.119,1805.2,1824.294,1845.755,1867.33,1888.227,
     &                     1909.905,1932.074,1953.679,1978.438,2003.606,2027.979,
     &                     2052.888,2078.937,2103.74,2132.809,2162.784,2191.951,
     &                     2221.497,2253.185,2283.029,2321.175,2361.957,2401.434,
     &                     2438.753,2475.441/)

      tcon_30X_MET_Mg2SiO4 = (/1515.139,1518.098,1526.614,1538.803,1553.121,1568.391,
     &                     1583.865,1599.336,1615.094,1631.141,1647.425,1663.959,
     &                     1680.846,1698.095,1715.812,1734.019,1752.553,1771.377,
     &                     1790.592,1810.227,1830.414,1851.239,1872.534,1894.2,
     &                     1916.314,1938.913,1962.139,1986.155,2010.807,2035.965,
     &                     2061.724,2088.169,2115.444,2143.761,2172.89,2202.563,
     &                     2232.84,2263.857,2295.826,2329.155,2363.713,2399.249,
     &                     2435.91,2474.006,2514.014,2556.48,2600.093,2642.635,
     &                     2682.975,2721.958/)

      tcon_100X_MET_Mg2SiO4 = (/1605.879,1614.989,1627.624,1642.701,1660.079,1678.019,
     &                     1695.481,1712.713,1730.817,1748.727,1767.603,1786.405,
     &                     1806.048,1825.757,1846.344,1867.062,1888.573,1910.294,
     &                     1932.708,1955.484,1978.987,2003.016,2027.763,2053.18,
     &                     2079.149,2105.841,2133.029,2161.17,2189.842,2219.755,
     &                     2250.15,2281.992,2314.246,2348.168,2382.306,2418.077,
     &                     2453.85,2491.634,2529.41,2569.526,2609.596,2652.547,
     &                     2695.463,2741.638,2787.465,2840.185,2892.989,2942.628,
     &                     2984.089,3016.856/)

      tcon_300X_MET_Mg2SiO4 = (/1695.727,1705.743,1719.844,1736.78,1756.345,1776.561,
     &                     1796.235,1815.673,1836.288,1856.799,1878.311,1899.746,
     &                     1922.222,1944.754,1968.283,1992.008,2016.626,2041.531,
     &                     2067.398,2093.778,2121.007,2148.847,2177.501,2206.986,
     &                     2237.295,2268.654,2300.702,2333.923,2367.794,2403.135,
     &                     2439.023,2476.559,2514.464,2554.078,2593.881,2635.778,
     &                     2677.928,2723.047,2768.47,2816.309,2863.766,2914.367,
     &                     2964.564,3018.16,3071.19,3132.28,3193.687,3251.822,
     &                     3300.426,3338.432/)

      tcon_1X_MET_MnS = (/1117.921,1121.455,1126.364,1132.372,1139.178,1146.25,
     &                     1152.965,1159.631,1166.38,1173.268,1180.298,1187.428,
     &                     1194.664,1201.941,1209.236,1216.545,1223.924,1231.508,
     &                     1239.214,1247.239,1255.519,1263.943,1272.368,1280.842,
     &                     1289.355,1297.912,1306.562,1315.288,1324.089,1333.085,
     &                     1342.232,1351.456,1360.718,1370.132,1379.691,1389.384,
     &                     1399.142,1408.991,1418.857,1428.73,1438.709,1448.903,
     &                     1459.399,1470.098,1481.11,1493.461,1505.809,1517.23,
     &                     1526.77,1534.083/)

      tcon_10X_MET_MnS = (/1138.599,1140.635,1146.066,1154.285,1164.199,1174.719,
     &                     1184.812,1194.154,1204.344,1215.127,1225.596,1235.763,
     &                     1246.881,1257.425,1268.474,1280.294,1291.817,1302.953,
     &                     1314.951,1326.672,1338.696,1351.911,1365.011,1377.679,
     &                     1391.029,1404.344,1417.571,1432.292,1446.978,1461.124,
     &                     1475.8,1490.832,1505.448,1522.147,1539.071,1555.41,
     &                     1572.071,1589.457,1605.92,1625.123,1644.761,1663.615,
     &                     1682.382,1702.117,1719.887,1740.263,1760.487,1778.959,
     &                     1796.038,1813.618/)

      tcon_30X_MET_MnS = (/1189.149,1191.877,1198.598,1208.01,1218.985,1230.637,
     &                     1242.391,1254.084,1265.919,1277.881,1289.924,1302.062,
     &                     1314.374,1326.851,1339.552,1352.483,1365.532,1378.686,
     &                     1392.034,1405.603,1419.478,1433.719,1448.229,1462.998,
     &                     1478.181,1493.938,1510.517,1528.181,1546.885,1566.325,
     &                     1585.94,1604.847,1622.138,1637.198,1649.68,1659.874,
     &                     1668.932,1678.402,1689.807,1704.35,1722.242,1743.087,
     &                     1766.415,1791.747,1818.762,1847.334,1876.448,1904.681,
     &                     1931.402,1957.449/)

      tcon_100X_MET_MnS = (/1240.618,1247.813,1257.874,1269.927,1283.755,1297.947,
     &                     1311.646,1325.087,1339.148,1352.89,1367.15,1381.229,
     &                     1395.882,1410.418,1425.434,1440.394,1455.794,1471.176,
     &                     1486.923,1502.816,1519.106,1535.58,1552.39,1569.558,
     &                     1587.064,1605.029,1623.221,1641.879,1661.059,1681.696,
     &                     1702.225,1721.993,1739.294,1752.768,1762.245,1768.401,
     &                     1772.619,1778.695,1787.913,1803.621,1824.078,1850.577,
     &                     1880.073,1912.921,1946.912,1987.173,2029.24,2069.897,
     &                     2104.96,2133.258/)

      tcon_300X_MET_MnS = (/1280.282,1288.302,1299.503,1312.884,1328.255,1344.061,
     &                     1359.433,1374.652,1390.741,1406.636,1423.237,1439.697,
     &                     1456.842,1473.901,1491.53,1509.128,1527.353,1545.696,
     &                     1564.538,1583.547,1603.004,1622.688,1642.936,1663.832,
     &                     1685.121,1706.683,1726.716,1743.625,1756.47,1764.573,
     &                     1769.838,1775.594,1783.384,1795.91,1813.332,1837.136,
     &                     1864.296,1894.408,1924.785,1955.707,1986.826,2020.325,
     &                     2053.795,2089.615,2125.467,2168.365,2212.434,2254.572,
     &                     2290.583,2319.594/)

      tcon_1X_MET_Na2S = (/780.091,783.188,787.451,792.657,798.527,804.6,
     &                     810.305,815.777,821.2,826.638,832.205,837.979,
     &                     844.016,850.308,856.825,863.521,870.318,877.173,
     &                     884.061,890.997,897.995,905.053,912.167,919.29,
     &                     926.421,933.59,940.877,948.383,956.115,964.078,
     &                     972.277,980.748,989.378,998.256,1007.267,1016.406,
     &                     1025.674,1035.081,1044.577,1054.127,1063.742,1073.443,
     &                     1083.193,1093.105,1103.094,1114.356,1125.696,1136.223,
     &                     1145.149,1152.113/)

      tcon_10X_MET_Na2S = (/798.505,800.227,804.817,811.766,820.144,829.028,
     &                     837.541,845.403,853.926,862.897,871.579,880.012,
     &                     889.286,898.136,907.539,917.738,927.797,937.593,
     &                     948.119,958.413,969.05,980.836,992.604,1004.051,
     &                     1016.129,1028.212,1040.346,1054.046,1067.868,1081.29,
     &                     1095.213,1109.441,1123.323,1139.262,1155.485,1171.208,
     &                     1187.272,1204.052,1220.02,1238.759,1258.06,1276.769,
     &                     1295.589,1315.587,1334.082,1356.752,1380.344,1402.717,
     &                     1423.719,1444.678/)

      tcon_30X_MET_Na2S = (/835.576,837.486,842.558,849.815,858.384,867.565,
     &                     876.894,886.227,895.717,905.351,915.098,924.986,
     &                     935.095,945.436,956.08,967.051,978.257,989.681,
     &                     1001.387,1013.397,1025.797,1038.64,1051.811,1065.248,
     &                     1079.003,1093.11,1107.659,1122.748,1138.257,1154.081,
     &                     1170.261,1186.842,1203.924,1221.663,1239.946,1258.626,
     &                     1277.747,1297.372,1317.59,1338.587,1360.176,1382.081,
     &                     1404.283,1426.873,1450.05,1474.11,1498.386,1521.759,
     &                     1543.705,1564.65/)

      tcon_100X_MET_Na2S = (/875.816,880.853,888.01,896.677,906.732,917.167,
     &                     927.375,937.469,948.06,958.504,969.457,980.419,
     &                     991.949,1003.506,1015.546,1027.663,1040.278,1053.085,
     &                     1066.346,1079.864,1093.898,1108.298,1123.104,1138.309,
     &                     1153.893,1169.994,1186.478,1203.581,1221.005,1239.19,
     &                     1257.671,1277.009,1296.586,1317.166,1337.942,1359.9,
     &                     1381.968,1405.324,1428.669,1453.335,1477.714,1503.335,
     &                     1528.38,1554.519,1579.81,1608.218,1636.123,1662.113,
     &                     1683.57,1700.519/)

      tcon_300X_MET_Na2S = (/913.263,919.009,927.074,936.734,947.837,959.272,
     &                     970.447,981.5,993.067,1004.505,1016.658,1028.881,
     &                     1041.713,1054.576,1067.958,1081.441,1095.514,1109.783,
     &                     1124.525,1139.559,1155.206,1171.311,1187.909,1204.976,
     &                     1222.501,1240.649,1259.224,1278.488,1298.155,1318.738,
     &                     1339.743,1361.846,1384.309,1408.0,1431.982,1457.329,
     &                     1482.775,1509.642,1536.474,1564.87,1593.019,1622.715,
     &                     1651.864,1682.597,1712.489,1745.748,1778.209,1808.215,
     &                     1832.834,1852.099/)

      tcon_1X_MET_Ni = (/1320.903,1326.055,1333.238,1342.122,1352.302,1363.043,
     &                     1373.467,1383.348,1393.013,1402.703,1412.471,1422.341,
     &                     1432.51,1443.21,1454.651,1466.474,1478.934,1491.219,
     &                     1503.686,1516.392,1529.172,1541.832,1554.555,1567.624,
     &                     1580.881,1594.113,1607.848,1621.712,1636.004,1650.944,
     &                     1665.966,1681.086,1696.724,1713.586,1730.663,1747.664,
     &                     1765.226,1782.989,1801.318,1819.904,1838.572,1857.358,
     &                     1877.201,1897.299,1916.613,1937.652,1958.762,1977.946,
     &                     1994.264,2006.786/)

      tcon_10X_MET_Ni = (/1302.846,1305.727,1313.396,1325.005,1339.016,1353.898,
     &                     1368.199,1381.475,1396.045,1411.555,1426.694,1441.459,
     &                     1457.605,1472.982,1489.28,1506.916,1524.276,1541.162,
     &                     1559.318,1577.075,1595.415,1615.722,1635.988,1655.698,
     &                     1676.501,1697.318,1718.227,1741.835,1765.668,1788.843,
     &                     1812.932,1837.626,1861.887,1890.031,1918.906,1947.069,
     &                     1975.896,2005.983,2034.753,2068.773,2104.045,2138.453,
     &                     2173.225,2210.25,2244.834,2288.209,2334.037,2377.995,
     &                     2419.437,2460.428/)

      tcon_30X_MET_Ni = (/1343.641,1346.865,1355.269,1367.233,1381.308,1396.347,
     &                     1411.613,1426.918,1442.574,1458.593,1474.911,1491.533,
     &                     1508.557,1525.979,1543.912,1562.377,1581.201,1600.343,
     &                     1619.915,1639.952,1660.606,1681.987,1703.934,1726.357,
     &                     1749.35,1772.954,1797.323,1822.627,1848.672,1875.292,
     &                     1902.579,1930.643,1959.704,1990.058,2021.485,2053.699,
     &                     2086.772,2120.858,2156.215,2193.308,2231.926,2271.706,
     &                     2312.783,2355.526,2400.523,2448.431,2497.721,2545.652,
     &                     2590.347,2631.385/)

      tcon_100X_MET_Ni = (/1387.193,1395.426,1406.972,1420.832,1436.784,1453.228,
     &                     1469.207,1484.951,1501.534,1518.001,1535.396,1552.819,
     &                     1571.157,1589.609,1608.89,1628.337,1648.529,1668.94,
     &                     1690.05,1711.517,1733.676,1756.352,1779.777,1804.009,
     &                     1828.957,1854.751,1881.144,1908.573,1936.579,1965.839,
     &                     1995.593,2026.76,2058.428,2091.888,2125.735,2161.506,
     &                     2197.61,2236.068,2274.843,2316.353,2358.066,2402.807,
     &                     2447.699,2496.328,2545.128,2602.729,2661.463,2717.335,
     &                     2764.708,2802.56/)

      tcon_300X_MET_Ni = (/1424.665,1433.321,1445.481,1460.05,1476.842,1494.154,
     &                     1510.989,1527.616,1545.207,1562.681,1581.073,1599.459,
     &                     1618.788,1638.206,1658.455,1678.835,1700.006,1721.446,
     &                     1743.668,1766.316,1789.801,1813.972,1838.993,1864.846,
     &                     1891.463,1919.051,1947.319,1976.672,2006.649,2038.051,
     &                     2070.112,2103.808,2138.071,2174.277,2211.042,2250.138,
     &                     2289.719,2331.928,2374.612,2420.479,2466.707,2516.405,
     &                     2566.356,2620.522,2674.966,2739.169,2804.558,2866.622,
     &                     2919.191,2961.105/)

      tcon_1X_MET_SiO2 = (/1340.013,1344.242,1350.143,1357.436,1365.799,1374.631,
     &                     1383.183,1391.618,1400.101,1408.657,1417.259,1425.875,
     &                     1434.495,1443.1,1451.697,1460.254,1468.764,1477.3,
     &                     1486.002,1494.99,1504.332,1514.045,1524.088,1534.356,
     &                     1544.736,1555.05,1565.217,1575.285,1585.407,1595.731,
     &                     1606.32,1617.198,1628.412,1639.994,1651.976,1664.289,
     &                     1676.926,1689.924,1703.31,1717.095,1731.283,1745.951,
     &                     1761.232,1777.207,1793.96,1813.082,1832.474,1850.476,
     &                     1865.656,1877.358/)

      tcon_10X_MET_SiO2 = (/1362.704,1365.088,1371.455,1381.085,1392.707,1405.044,
     &                     1416.89,1427.868,1439.879,1452.624,1465.027,1477.093,
     &                     1490.279,1502.796,1515.952,1530.069,1543.864,1557.217,
     &                     1571.59,1585.628,1600.035,1615.872,1631.583,1646.792,
     &                     1662.836,1678.871,1694.874,1712.798,1730.766,1748.13,
     &                     1766.136,1784.549,1802.446,1822.872,1843.577,1863.596,
     &                     1884.072,1905.539,1926.026,1950.093,1974.993,1999.331,
     &                     2024.117,2050.843,2076.33,2109.813,2146.215,2181.882,
     &                     2215.743,2248.733/)

      tcon_30X_MET_SiO2 = (/1431.943,1433.934,1441.0,1451.327,1463.52,1476.537,
     &                     1489.713,1502.869,1516.252,1529.865,1543.658,1557.654,
     &                     1571.951,1586.551,1601.539,1616.915,1632.52,1648.316,
     &                     1664.386,1680.75,1697.512,1714.751,1732.346,1750.244,
     &                     1768.529,1787.233,1806.453,1826.298,1846.602,1867.224,
     &                     1888.232,1909.696,1931.739,1954.537,1977.899,2001.618,
     &                     2025.766,2050.467,2075.882,2102.338,2129.741,2157.924,
     &                     2187.076,2217.564,2249.987,2285.026,2321.651,2357.943,
     &                     2392.88,2427.24/)

      tcon_100X_MET_SiO2 = (/1512.783,1520.334,1530.939,1543.685,1558.45,1573.735,
     &                     1588.564,1603.131,1618.399,1633.424,1649.125,1664.73,
     &                     1681.066,1697.434,1714.586,1731.881,1749.775,1767.785,
     &                     1786.311,1805.04,1824.274,1843.822,1863.86,1884.439,
     &                     1905.507,1927.16,1949.223,1972.091,1995.307,2019.324,
     &                     2043.542,2068.755,2094.219,2120.984,2147.822,2175.801,
     &                     2203.733,2233.22,2262.603,2293.692,2324.702,2357.963,
     &                     2391.266,2427.092,2462.748,2504.636,2547.231,2587.752,
     &                     2622.042,2649.471/)

      tcon_300X_MET_SiO2 = (/1588.42,1596.918,1608.833,1623.105,1639.616,1656.698,
     &                     1673.367,1689.824,1707.102,1724.166,1742.04,1759.767,
     &                     1778.235,1796.684,1815.925,1835.293,1855.502,1876.044,
     &                     1897.333,1918.971,1941.219,1963.83,1986.942,2010.524,
     &                     2034.546,2059.237,2084.456,2110.62,2137.242,2164.916,
     &                     2192.945,2222.24,2251.779,2282.565,2313.456,2346.025,
     &                     2378.757,2413.452,2448.258,2485.39,2522.643,2562.603,
     &                     2602.461,2645.164,2687.705,2737.486,2787.873,2835.545,
     &                     2875.773,2907.801/)

      tcon_1X_MET_CaTiO3 = (/1601.37,1605.679,1611.878,1619.741,1628.97,1638.929,
     &                     1648.815,1658.732,1668.75,1678.812,1689.02,1699.041,
     &                     1709.2,1719.397,1729.603,1739.8,1749.96,1759.896,
     &                     1769.706,1779.64,1790.163,1801.136,1812.637,1824.678,
     &                     1837.203,1850.166,1863.18,1876.001,1888.737,1901.521,
     &                     1914.396,1927.033,1939.747,1952.806,1966.166,1979.939,
     &                     1993.987,2008.389,2023.103,2038.033,2053.04,2067.937,
     &                     2082.917,2097.999,2113.085,2129.325,2145.154,2159.427,
     &                     2171.174,2180.081/)

      tcon_10X_MET_CaTiO3 = (/1660.897,1663.954,1672.115,1684.466,1699.349,1715.113,
     &                     1730.184,1744.049,1758.941,1774.484,1789.42,1803.873,
     &                     1819.822,1835.041,1851.188,1868.677,1885.892,1902.614,
     &                     1920.543,1938.001,1955.853,1975.396,1994.717,2013.376,
     &                     2033.068,2052.752,2072.501,2094.856,2117.191,2138.354,
     &                     2159.452,2179.701,2196.421,2210.068,2220.052,2227.41,
     &                     2235.699,2247.399,2261.1,2282.689,2309.028,2337.518,
     &                     2366.778,2395.883,2422.931,2456.523,2491.773,2525.42,
     &                     2557.079,2588.518/)

      tcon_30X_MET_CaTiO3 = (/1733.44,1737.232,1746.591,1759.779,1775.208,1791.608,
     &                     1808.156,1824.668,1841.643,1859.507,1878.744,1899.736,
     &                     1922.56,1946.811,1971.972,1997.374,2021.914,2043.755,
     &                     2060.545,2070.89,2075.543,2076.574,2076.32,2077.335,
     &                     2082.209,2092.886,2110.071,2132.532,2157.541,2183.053,
     &                     2208.372,2233.062,2256.596,2278.753,2299.69,2320.13,
     &                     2341.357,2364.821,2392.003,2424.252,2461.648,2503.234,
     &                     2547.726,2593.555,2638.567,2680.519,2717.084,2746.86,
     &                     2770.239,2789.211/)

      tcon_100X_MET_CaTiO3 = (/1809.311,1818.823,1832.127,1848.077,1866.423,1885.316,
     &                     1903.669,1921.773,1940.849,1959.699,1979.415,1998.989,
     &                     2019.461,2039.926,2061.171,2082.456,2107.622,2136.238,
     &                     2157.753,2164.85,2161.778,2152.324,2140.784,2131.262,
     &                     2128.0,2135.502,2157.771,2191.522,2224.368,2254.326,
     &                     2284.783,2316.623,2348.765,2382.46,2416.405,2452.095,
     &                     2487.857,2525.639,2563.358,2603.332,2643.17,2685.591,
     &                     2727.94,2778.027,2828.509,2872.365,2906.082,2929.883,
     &                     2942.165,2945.745/)

      tcon_300X_MET_CaTiO3 = (/1878.949,1889.406,1904.021,1921.498,1941.641,1962.412,
     &                     1982.51,2002.248,2027.215,2054.498,2071.984,2075.332,
     &                     2069.744,2058.384,2046.366,2036.972,2035.001,2043.967,
     &                     2068.57,2101.867,2131.835,2159.516,2187.982,2217.226,
     &                     2247.233,2278.221,2309.911,2342.825,2376.379,2411.376,
     &                     2446.952,2484.188,2521.774,2561.004,2600.369,2641.767,
     &                     2684.445,2736.391,2787.618,2822.44,2841.857,2851.693,
     &                     2853.831,2859.278,2867.771,2884.337,2909.799,2946.68,
     &                     2978.821,2994.89/)

      tcon_1X_MET_VO = (/1368.564,1373.377,1380.144,1388.754,1398.567,1409.202,
     &                     1419.4,1429.398,1439.388,1449.485,1459.919,1470.367,
     &                     1481.16,1492.267,1503.385,1514.468,1525.433,1536.354,
     &                     1547.177,1557.984,1569.055,1580.219,1591.64,1603.412,
     &                     1615.364,1627.414,1639.905,1652.463,1665.057,1677.717,
     &                     1690.573,1703.516,1716.555,1729.951,1743.339,1756.745,
     &                     1770.21,1783.715,1797.248,1810.885,1824.34,1838.024,
     &                     1852.314,1867.072,1881.985,1898.838,1915.896,1931.733,
     &                     1945.137,1955.561/)

      tcon_10X_MET_VO = (/1325.753,1328.318,1335.133,1345.452,1357.915,1371.177,
     &                     1383.959,1395.885,1409.131,1423.377,1437.403,1451.152,
     &                     1466.134,1480.425,1495.696,1512.344,1528.785,1544.738,
     &                     1561.712,1578.059,1594.338,1611.605,1628.289,1644.223,
     &                     1661.275,1678.644,1696.438,1717.049,1738.222,1759.007,
     &                     1780.476,1802.2,1823.241,1847.144,1871.256,1894.441,
     &                     1918.045,1942.672,1965.868,1992.59,2019.729,2045.734,
     &                     2071.789,2099.585,2125.152,2156.066,2187.947,2217.964,
     &                     2246.073,2274.271/)

      tcon_30X_MET_VO = (/1346.645,1349.366,1356.633,1367.008,1379.225,1392.294,
     &                     1405.592,1418.972,1432.718,1446.823,1461.195,1475.81,
     &                     1490.739,1505.994,1521.715,1537.944,1554.501,1571.295,
     &                     1588.341,1605.569,1623.009,1640.7,1658.53,1676.533,
     &                     1694.939,1713.908,1733.636,1754.286,1775.652,1797.524,
     &                     1819.938,1842.974,1866.824,1891.768,1917.681,1944.379,
     &                     1971.922,2000.34,2029.658,2060.072,2091.326,2123.078,
     &                     2155.357,2188.289,2222.027,2256.794,2291.466,2324.252,
     &                     2354.041,2380.72/)

      tcon_100X_MET_VO = (/1366.375,1373.408,1383.227,1394.981,1408.475,1422.36,
     &                     1435.882,1449.258,1463.346,1477.272,1491.885,1506.416,
     &                     1521.562,1536.65,1552.318,1568.074,1584.466,1601.018,
     &                     1618.068,1635.336,1653.127,1671.281,1689.85,1708.789,
     &                     1728.121,1748.0,1768.223,1789.105,1810.229,1832.001,
     &                     1854.076,1877.378,1901.196,1926.602,1952.689,1980.811,
     &                     2009.563,2040.328,2071.238,2103.959,2136.583,2171.196,
     &                     2205.421,2241.875,2277.802,2318.674,2359.199,2397.022,
     &                     2428.416,2453.156/)

      tcon_300X_MET_VO = (/1380.633,1387.969,1398.145,1410.233,1424.109,1438.359,
     &                     1452.146,1465.706,1479.985,1494.118,1509.034,1523.906,
     &                     1539.407,1554.893,1571.057,1587.34,1604.218,1621.22,
     &                     1638.769,1656.614,1675.033,1693.827,1713.074,1732.743,
     &                     1752.853,1773.558,1794.619,1816.381,1838.492,1861.425,
     &                     1884.595,1908.731,1933.088,1958.718,1984.61,2011.944,
     &                     2039.414,2068.463,2097.459,2128.136,2158.622,2190.746,
     &                     2222.637,2257.425,2292.503,2334.206,2376.939,2417.811,
     &                     2452.648,2480.623/)

      tcon_1X_MET_ZnS = (/710.495,712.638,715.382,718.984,722.905,727.223,
     &                     731.19,735.39,739.67,744.114,748.188,752.535,
     &                     756.546,760.708,764.903,769.431,773.769,778.396,
     &                     783.226,788.153,793.051,797.691,802.186,806.909,
     &                     811.973,817.205,822.547,828.179,834.156,839.964,
     &                     845.819,851.23,856.836,862.449,868.048,873.656,
     &                     879.258,884.632,890.243,895.966,902.024,908.396,
     &                     915.213,922.215,929.212,936.786,944.085,950.425,
     &                     955.764,959.688/)

      tcon_10X_MET_ZnS = (/696.332,697.452,700.455,704.992,710.465,716.265,
     &                     721.821,726.947,732.494,738.322,743.951,749.408,
     &                     755.402,761.102,767.11,773.574,779.904,786.038,
     &                     792.636,799.075,805.686,812.956,820.164,827.132,
     &                     834.465,841.765,848.984,856.973,864.906,872.517,
     &                     880.41,888.494,896.325,905.222,914.198,922.832,
     &                     931.629,940.814,949.482,959.501,969.703,979.504,
     &                     989.335,999.822,1009.49,1021.25,1033.42,1044.921,
     &                     1055.688,1066.487/)

      tcon_30X_MET_ZnS = (/726.669,727.281,730.433,735.157,740.759,746.735,
     &                     752.774,758.787,764.885,771.079,777.356,783.732,
     &                     790.248,796.893,803.691,810.63,817.636,824.698,
     &                     831.864,839.161,846.653,854.385,862.295,870.344,
     &                     878.551,886.913,895.46,904.238,913.189,922.274,
     &                     931.533,940.99,950.691,960.713,970.976,981.373,
     &                     991.905,1002.581,1013.431,1024.556,1035.903,1047.453,
     &                     1059.412,1072.072,1085.715,1100.541,1116.033,1131.259,
     &                     1145.541,1158.861/)

      tcon_100X_MET_ZnS = (/765.566,769.002,773.79,779.517,786.138,792.984,
     &                     799.621,806.131,812.941,819.67,826.791,833.911,
     &                     841.362,848.827,856.58,864.296,872.177,880.016,
     &                     888.018,896.077,904.386,912.888,921.666,930.707,
     &                     939.967,949.509,959.222,969.232,979.405,990.033,
     &                     1000.804,1011.976,1023.237,1035.109,1047.087,1059.68,
     &                     1072.243,1085.381,1098.437,1112.259,1125.886,1140.131,
     &                     1154.12,1168.925,1183.43,1200.17,1216.931,1232.704,
     &                     1245.968,1256.597/)

      tcon_300X_MET_ZnS = (/800.763,804.592,809.913,816.254,823.635,831.303,
     &                     838.72,845.941,853.465,860.92,868.816,876.656,
     &                     884.797,892.952,901.491,910.031,918.758,927.479,
     &                     936.457,945.53,954.837,964.317,974.061,984.046,
     &                     994.28,1004.886,1015.709,1026.853,1038.162,1049.942,
     &                     1061.937,1074.577,1087.409,1100.88,1114.558,1129.166,
     &                     1143.935,1159.611,1175.481,1192.905,1210.642,1230.244,
     &                     1249.322,1266.625,1280.679,1291.323,1297.8,1301.014,
     &                     1301.203,1300.468/)

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

      tconds(2,1:51,1)=tcon_10X_MET_KCl
      tconds(2,1:51,2)=tcon_10X_MET_ZnS
      tconds(2,1:51,3)=tcon_10X_MET_Na2S
      tconds(2,1:51,4)=tcon_10X_MET_MnS
      tconds(2,1:51,5)=tcon_10X_MET_Cr
      tconds(2,1:51,6)=tcon_10X_MET_SiO2
      tconds(2,1:51,7)=tcon_10X_MET_Mg2SiO4
      tconds(2,1:51,8)=tcon_10X_MET_VO
      tconds(2,1:51,9)=tcon_10X_MET_Ni
      tconds(2,1:51,10)=tcon_10X_MET_Fe
      tconds(2,1:51,11)=tcon_10X_MET_CaSiO4
      tconds(2,1:51,12)=tcon_10X_MET_CaTiO3
      tconds(2,1:51,13)=tcon_10X_MET_Al2O3

      tconds(3,1:51,1)=tcon_30X_MET_KCl
      tconds(3,1:51,2)=tcon_30X_MET_ZnS
      tconds(3,1:51,3)=tcon_30X_MET_Na2S
      tconds(3,1:51,4)=tcon_30X_MET_MnS
      tconds(3,1:51,5)=tcon_30X_MET_Cr
      tconds(3,1:51,6)=tcon_30X_MET_SiO2
      tconds(3,1:51,7)=tcon_30X_MET_Mg2SiO4
      tconds(3,1:51,8)=tcon_30X_MET_VO
      tconds(3,1:51,9)=tcon_30X_MET_Ni
      tconds(3,1:51,10)=tcon_30X_MET_Fe
      tconds(3,1:51,11)=tcon_30X_MET_CaSiO4
      tconds(3,1:51,12)=tcon_30X_MET_CaTiO3
      tconds(3,1:51,13)=tcon_30X_MET_Al2O3

      tconds(4,1:51,1)=tcon_100X_MET_KCl
      tconds(4,1:51,2)=tcon_100X_MET_ZnS
      tconds(4,1:51,3)=tcon_100X_MET_Na2S
      tconds(4,1:51,4)=tcon_100X_MET_MnS
      tconds(4,1:51,5)=tcon_100X_MET_Cr
      tconds(4,1:51,6)=tcon_100X_MET_SiO2
      tconds(4,1:51,7)=tcon_100X_MET_Mg2SiO4
      tconds(4,1:51,8)=tcon_100X_MET_VO
      tconds(4,1:51,9)=tcon_100X_MET_Ni
      tconds(4,1:51,10)=tcon_100X_MET_Fe
      tconds(4,1:51,11)=tcon_100X_MET_CaSiO4
      tconds(4,1:51,12)=tcon_100X_MET_CaTiO3
      tconds(4,1:51,13)=tcon_100X_MET_Al2O3

      tconds(5,1:51,1)=tcon_300X_MET_KCl
      tconds(5,1:51,2)=tcon_300X_MET_ZnS
      tconds(5,1:51,3)=tcon_300X_MET_Na2S
      tconds(5,1:51,4)=tcon_300X_MET_MnS
      tconds(5,1:51,5)=tcon_300X_MET_Cr
      tconds(5,1:51,6)=tcon_300X_MET_SiO2
      tconds(5,1:51,7)=tcon_300X_MET_Mg2SiO4
      tconds(5,1:51,8)=tcon_300X_MET_VO
      tconds(5,1:51,9)=tcon_300X_MET_Ni
      tconds(5,1:51,10)=tcon_300X_MET_Fe
      tconds(5,1:51,11)=tcon_300X_MET_CaSiO4
      tconds(5,1:51,12)=tcon_300X_MET_CaTiO3
      tconds(5,1:51,13)=tcon_300X_MET_Al2O3

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
