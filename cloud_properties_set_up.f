      SUBROUTINE get_cloud_scattering_properties_wrapper
          include 'rcommons.h'
          call get_cloud_scattering_properties(NCLOUDS, NLAYER, NVERT, NIRP, NSOLP, GASCON, METALLICITY)
      END SUBROUTINE get_cloud_scattering_properties_wrapper


      SUBROUTINE get_cloud_scattering_properties(NCLOUDS, NLAYER, NVERT, NIRP, NSOLP, GASCON, METALLICITY)
          implicit none
          integer :: I, J, L, K, NL, NCLOUDS, NLAYER, NVERT, NIRP, NSOLP
          real :: GAS_CONSTANT_R, GASCON, METALLICITY
          integer :: MET_INDEX, LOW_MET_INDEX, HIGH_MET_INDEX
          real :: COND_CURVE_METS(5) ! These are log10(Z/Zsun) values for the pre-calculated curves

          character (len = 40) :: haze_type

          ! Define all the arrays

          ! HAZE ARRAYS ARE DIFFERENT THAN THE OTHER ONES
          real, dimension(50, 100)  :: HAZE_RosselandMean_tau_per_bar, HAZE_RosselandMean_pi0, HAZE_RosselandMean_gg
          real, dimension(50, 100)  :: HAZE_PlanckMean_tau_per_bar, HAZE_PlanckMean_pi0, HAZE_PlanckMean_gg
          real, dimension(500, 100) :: HAZE_wav_tau_per_bar, HAZE_wav_pi0, HAZE_wav_gg
          real, dimension(100)      :: haze_pressure_array_pascals

          ! These are 50 by 50 because that's what the data in CLOUD_DATA is
          ! That can change but use to code from Elsie and Isaac
          real, dimension(100,100) :: KCl_rosselandMean_gg
          real, dimension(100,100) :: KCl_rosselandMean_kext
          real, dimension(100,100) :: KCl_rosselandMean_pi0

          real, dimension(100,100) :: KCl_PlanckMean_gg
          real, dimension(100,100) :: KCl_PlanckMean_kext
          real, dimension(100,100) :: KCl_PlanckMean_pi0

          real, dimension(100,100) :: KCl_wav_gg
          real, dimension(100,100) :: KCl_wav_kext
          real, dimension(100,100) :: KCl_wav_pi0


          real, dimension(100,100) :: ZnS_rosselandMean_gg
          real, dimension(100,100) :: ZnS_rosselandMean_kext
          real, dimension(100,100) :: ZnS_rosselandMean_pi0

          real, dimension(100,100) :: ZnS_PlanckMean_gg
          real, dimension(100,100) :: ZnS_PlanckMean_kext
          real, dimension(100,100) :: ZnS_PlanckMean_pi0

          real, dimension(100,100) :: ZnS_wav_gg
          real, dimension(100,100) :: ZnS_wav_kext
          real, dimension(100,100) :: ZnS_wav_pi0


          real, dimension(100,100) :: Na2S_rosselandMean_gg
          real, dimension(100,100) :: Na2S_rosselandMean_kext
          real, dimension(100,100) :: Na2S_rosselandMean_pi0

          real, dimension(100,100) :: Na2S_PlanckMean_gg
          real, dimension(100,100) :: Na2S_PlanckMean_kext
          real, dimension(100,100) :: Na2S_PlanckMean_pi0

          real, dimension(100,100) :: Na2S_wav_gg
          real, dimension(100,100) :: Na2S_wav_kext
          real, dimension(100,100) :: Na2S_wav_pi0


          real, dimension(100,100) :: MnS_rosselandMean_gg
          real, dimension(100,100) :: MnS_rosselandMean_kext
          real, dimension(100,100) :: MnS_rosselandMean_pi0

          real, dimension(100,100) :: MnS_PlanckMean_gg
          real, dimension(100,100) :: MnS_PlanckMean_kext
          real, dimension(100,100) :: MnS_PlanckMean_pi0

          real, dimension(100,100) :: MnS_wav_gg
          real, dimension(100,100) :: MnS_wav_kext
          real, dimension(100,100) :: MnS_wav_pi0


          real, dimension(100,100) :: Cr_rosselandMean_gg
          real, dimension(100,100) :: Cr_rosselandMean_kext
          real, dimension(100,100) :: Cr_rosselandMean_pi0

          real, dimension(100,100) :: Cr_PlanckMean_gg
          real, dimension(100,100) :: Cr_PlanckMean_kext
          real, dimension(100,100) :: Cr_PlanckMean_pi0

          real, dimension(100,100) :: Cr_wav_gg
          real, dimension(100,100) :: Cr_wav_kext
          real, dimension(100,100) :: Cr_wav_pi0


          real, dimension(100,100) :: SiO2_rosselandMean_gg
          real, dimension(100,100) :: SiO2_rosselandMean_kext
          real, dimension(100,100) :: SiO2_rosselandMean_pi0

          real, dimension(100,100) :: SiO2_PlanckMean_gg
          real, dimension(100,100) :: SiO2_PlanckMean_kext
          real, dimension(100,100) :: SiO2_PlanckMean_pi0

          real, dimension(100,100) :: SiO2_wav_gg
          real, dimension(100,100) :: SiO2_wav_kext
          real, dimension(100,100) :: SiO2_wav_pi0


          real, dimension(100,100) :: Mg2SiO4_rosselandMean_gg
          real, dimension(100,100) :: Mg2SiO4_rosselandMean_kext
          real, dimension(100,100) :: Mg2SiO4_rosselandMean_pi0

          real, dimension(100,100) :: Mg2SiO4_PlanckMean_gg
          real, dimension(100,100) :: Mg2SiO4_PlanckMean_kext
          real, dimension(100,100) :: Mg2SiO4_PlanckMean_pi0

          real, dimension(100,100) :: Mg2SiO4_wav_gg
          real, dimension(100,100) :: Mg2SiO4_wav_kext
          real, dimension(100,100) :: Mg2SiO4_wav_pi0


          real, dimension(100,100) :: VO_rosselandMean_gg
          real, dimension(100,100) :: VO_rosselandMean_kext
          real, dimension(100,100) :: VO_rosselandMean_pi0

          real, dimension(100,100) :: VO_PlanckMean_gg
          real, dimension(100,100) :: VO_PlanckMean_kext
          real, dimension(100,100) :: VO_PlanckMean_pi0

          real, dimension(100,100) :: VO_wav_gg
          real, dimension(100,100) :: VO_wav_kext
          real, dimension(100,100) :: VO_wav_pi0


          real, dimension(100,100) :: Ni_rosselandMean_gg
          real, dimension(100,100) :: Ni_rosselandMean_kext
          real, dimension(100,100) :: Ni_rosselandMean_pi0

          real, dimension(100,100) :: Ni_PlanckMean_gg
          real, dimension(100,100) :: Ni_PlanckMean_kext
          real, dimension(100,100) :: Ni_PlanckMean_pi0

          real, dimension(100,100) :: Ni_wav_gg
          real, dimension(100,100) :: Ni_wav_kext
          real, dimension(100,100) :: Ni_wav_pi0


          real, dimension(100,100) :: Fe_rosselandMean_gg
          real, dimension(100,100) :: Fe_rosselandMean_kext
          real, dimension(100,100) :: Fe_rosselandMean_pi0

          real, dimension(100,100) :: Fe_PlanckMean_gg
          real, dimension(100,100) :: Fe_PlanckMean_kext
          real, dimension(100,100) :: Fe_PlanckMean_pi0

          real, dimension(100,100) :: Fe_wav_gg
          real, dimension(100,100) :: Fe_wav_kext
          real, dimension(100,100) :: Fe_wav_pi0


          real, dimension(100,100) :: CaSiO4_rosselandMean_gg
          real, dimension(100,100) :: CaSiO4_rosselandMean_kext
          real, dimension(100,100) :: CaSiO4_rosselandMean_pi0

          real, dimension(100,100) :: CaSiO4_PlanckMean_gg
          real, dimension(100,100) :: CaSiO4_PlanckMean_kext
          real, dimension(100,100) :: CaSiO4_PlanckMean_pi0

          real, dimension(100,100) :: CaSiO4_wav_gg
          real, dimension(100,100) :: CaSiO4_wav_kext
          real, dimension(100,100) :: CaSiO4_wav_pi0


          real, dimension(100,100) :: CaTiO3_rosselandMean_gg
          real, dimension(100,100) :: CaTiO3_rosselandMean_kext
          real, dimension(100,100) :: CaTiO3_rosselandMean_pi0

          real, dimension(100,100) :: CaTiO3_PlanckMean_gg
          real, dimension(100,100) :: CaTiO3_PlanckMean_kext
          real, dimension(100,100) :: CaTiO3_PlanckMean_pi0

          real, dimension(100,100) :: CaTiO3_wav_gg
          real, dimension(100,100) :: CaTiO3_wav_kext
          real, dimension(100,100) :: CaTiO3_wav_pi0


          real, dimension(100,100) :: Al2O3_rosselandMean_gg
          real, dimension(100,100) :: Al2O3_rosselandMean_kext
          real, dimension(100,100) :: Al2O3_rosselandMean_pi0

          real, dimension(100,100) :: Al2O3_PlanckMean_gg
          real, dimension(100,100) :: Al2O3_PlanckMean_kext
          real, dimension(100,100) :: Al2O3_PlanckMean_pi0

          real, dimension(100,100) :: Al2O3_wav_gg
          real, dimension(100,100) :: Al2O3_wav_kext
          real, dimension(100,100) :: Al2O3_wav_pi0

          ! SET UP THE CONDENSATION CURVES
          ! SHOULD BE MET DEPENDENT EVENTUALLY

          REAL CORFACT(80)
          REAL TCONDS(6, 80, 13)
          REAL DUMMY_TCONDS(NCLOUDS+1, 80) ! +1 is pressure axis in the files

          REAL KE_OPPR(5, 100, 100, 13)
          REAL PI0_OPPR(5, 100, 100, 13)
          REAL G0_OPPR(5, 100, 100, 13)

          REAL DENSITY(13)
          REAL FMOLW(13)
          REAL CLOUD_MOLAR_MASSES(13)

          

          real, dimension(100) :: input_temperature_array
          real, dimension(80) :: input_pressure_array_cgs

          real, dimension(100) :: input_particle_size_array_in_meters
          real, dimension(80) :: particle_size_vs_layer_array_in_meters

          REAL, dimension (100)  :: CLOUD_WAV_GRID
          REAL, dimension (500) :: HAZE_WAV_GRID
          REAL :: exp_92_lnsig2_pi, sigma

      COMMON /CLOUD_PROPERTIES/ TCONDS, KE_OPPR, PI0_OPPR, G0_OPPR,
     &                              DENSITY, FMOLW,
     &                              CORFACT,
     &                              input_particle_size_array_in_meters,
     &                              input_temperature_array,
     &                              particle_size_vs_layer_array_in_meters,
     &                              input_pressure_array_cgs,
     &                              HAZE_RosselandMean_tau_per_bar, HAZE_RosselandMean_pi0, HAZE_RosselandMean_gg,
     &                              HAZE_PlanckMean_tau_per_bar,HAZE_PlanckMean_pi0, HAZE_PlanckMean_gg,
     &                              HAZE_wav_tau_per_bar,HAZE_wav_pi0, HAZE_wav_gg,
     &                              haze_pressure_array_pascals, HAZE_WAV_GRID, CLOUD_WAV_GRID, exp_92_lnsig2_pi

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
          open (2, file='../CLOUD_DATA/KCl_rosselandMean_kext.txt')
          open (3, file='../CLOUD_DATA/KCl_rosselandMean_pi0.txt')
          open (4, file='../CLOUD_DATA/KCl_PlanckMean_gg.txt')
          open (5, file='../CLOUD_DATA/KCl_PlanckMean_kext.txt')
          open (6, file='../CLOUD_DATA/KCl_PlanckMean_pi0.txt')
          open (7, file='../CLOUD_DATA/KCl_wav_gg.txt')
          open (8, file='../CLOUD_DATA/KCl_wav_kext.txt')
          open (9, file='../CLOUD_DATA/KCl_wav_pi0.txt')
          
          read(1,*) KCl_rosselandMean_gg
          read(2,*) KCl_rosselandMean_kext
          read(3,*) KCl_rosselandMean_pi0
          read(4,*) KCl_PlanckMean_gg
          read(5,*) KCl_PlanckMean_kext
          read(6,*) KCl_PlanckMean_pi0
          read(7,*) KCl_wav_gg
          read(8,*) KCl_wav_kext
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
          open (2, file='../CLOUD_DATA/ZnS_rosselandMean_kext.txt')
          open (3, file='../CLOUD_DATA/ZnS_rosselandMean_pi0.txt')
          open (4, file='../CLOUD_DATA/ZnS_wav_gg.txt')
          open (5, file='../CLOUD_DATA/ZnS_wav_kext.txt')
          open (6, file='../CLOUD_DATA/ZnS_wav_pi0.txt')
          open (7, file='../CLOUD_DATA/ZnS_PlanckMean_gg.txt')
          open (8, file='../CLOUD_DATA/ZnS_PlanckMean_kext.txt')
          open (9, file='../CLOUD_DATA/ZnS_PlanckMean_pi0.txt')
          
          read(1,*) ZnS_rosselandMean_gg
          read(2,*) ZnS_rosselandMean_kext
          read(3,*) ZnS_rosselandMean_pi0
          read(4,*) ZnS_wav_gg
          read(5,*) ZnS_wav_kext
          read(6,*) ZnS_wav_pi0
          read(7,*) ZnS_PlanckMean_gg
          read(8,*) ZnS_PlanckMean_kext
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
          open (2, file='../CLOUD_DATA/Na2S_rosselandMean_kext.txt')
          open (3, file='../CLOUD_DATA/Na2S_rosselandMean_pi0.txt')
          open (4, file='../CLOUD_DATA/Na2S_wav_gg.txt')
          open (5, file='../CLOUD_DATA/Na2S_wav_kext.txt')
          open (6, file='../CLOUD_DATA/Na2S_wav_pi0.txt')
          open (7, file='../CLOUD_DATA/Na2S_PlanckMean_gg.txt')
          open (8, file='../CLOUD_DATA/Na2S_PlanckMean_kext.txt')
          open (9, file='../CLOUD_DATA/Na2S_PlanckMean_pi0.txt')
          
          read(1,*) Na2S_rosselandMean_gg
          read(2,*) Na2S_rosselandMean_kext
          read(3,*) Na2S_rosselandMean_pi0
          read(4,*) Na2S_wav_gg
          read(5,*) Na2S_wav_kext
          read(6,*) Na2S_wav_pi0
          read(7,*) Na2S_PlanckMean_gg
          read(8,*) Na2S_PlanckMean_kext
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
          open (2, file='../CLOUD_DATA/MnS_rosselandMean_kext.txt')
          open (3, file='../CLOUD_DATA/MnS_rosselandMean_pi0.txt')
          open (4, file='../CLOUD_DATA/MnS_wav_gg.txt')
          open (5, file='../CLOUD_DATA/MnS_wav_kext.txt')
          open (6, file='../CLOUD_DATA/MnS_wav_pi0.txt')
          open (7, file='../CLOUD_DATA/MnS_PlanckMean_gg.txt')
          open (8, file='../CLOUD_DATA/MnS_PlanckMean_kext.txt')
          open (9, file='../CLOUD_DATA/MnS_PlanckMean_pi0.txt')
          
          read(1,*) MnS_rosselandMean_gg
          read(2,*) MnS_rosselandMean_kext
          read(3,*) MnS_rosselandMean_pi0
          read(4,*) MnS_wav_gg
          read(5,*) MnS_wav_kext
          read(6,*) MnS_wav_pi0
          read(7,*) MnS_PlanckMean_gg
          read(8,*) MnS_PlanckMean_kext
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
          open (2, file='../CLOUD_DATA/Cr_rosselandMean_kext.txt')
          open (3, file='../CLOUD_DATA/Cr_rosselandMean_pi0.txt')
          open (4, file='../CLOUD_DATA/Cr_wav_gg.txt')
          open (5, file='../CLOUD_DATA/Cr_wav_kext.txt')
          open (6, file='../CLOUD_DATA/Cr_wav_pi0.txt')
          open (7, file='../CLOUD_DATA/Cr_PlanckMean_gg.txt')
          open (8, file='../CLOUD_DATA/Cr_PlanckMean_kext.txt')
          open (9, file='../CLOUD_DATA/Cr_PlanckMean_pi0.txt')
          
          read(1,*) Cr_rosselandMean_gg
          read(2,*) Cr_rosselandMean_kext
          read(3,*) Cr_rosselandMean_pi0
          read(4,*) Cr_wav_gg
          read(5,*) Cr_wav_kext
          read(6,*) Cr_wav_pi0
          read(7,*) Cr_PlanckMean_gg
          read(8,*) Cr_PlanckMean_kext
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
          open (2, file='../CLOUD_DATA/SiO2_rosselandMean_kext.txt')
          open (3, file='../CLOUD_DATA/SiO2_rosselandMean_pi0.txt')
          open (4, file='../CLOUD_DATA/SiO2_wav_gg.txt')
          open (5, file='../CLOUD_DATA/SiO2_wav_kext.txt')
          open (6, file='../CLOUD_DATA/SiO2_wav_pi0.txt')
          open (7, file='../CLOUD_DATA/SiO2_PlanckMean_gg.txt')
          open (8, file='../CLOUD_DATA/SiO2_PlanckMean_kext.txt')
          open (9, file='../CLOUD_DATA/SiO2_PlanckMean_pi0.txt')
          
          read(1,*) SiO2_rosselandMean_gg
          read(2,*) SiO2_rosselandMean_kext
          read(3,*) SiO2_rosselandMean_pi0
          read(4,*) SiO2_wav_gg
          read(5,*) SiO2_wav_kext
          read(6,*) SiO2_wav_pi0
          read(7,*) SiO2_PlanckMean_gg
          read(8,*) SiO2_PlanckMean_kext
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
          open (2, file='../CLOUD_DATA/Mg2SiO4_rosselandMean_kext.txt')
          open (3, file='../CLOUD_DATA/Mg2SiO4_rosselandMean_pi0.txt')
          open (4, file='../CLOUD_DATA/Mg2SiO4_wav_gg.txt')
          open (5, file='../CLOUD_DATA/Mg2SiO4_wav_kext.txt')
          open (6, file='../CLOUD_DATA/Mg2SiO4_wav_pi0.txt')
          open (7, file='../CLOUD_DATA/Mg2SiO4_PlanckMean_gg.txt')
          open (8, file='../CLOUD_DATA/Mg2SiO4_PlanckMean_kext.txt')
          open (9, file='../CLOUD_DATA/Mg2SiO4_PlanckMean_pi0.txt')
          
          read(1,*) Mg2SiO4_rosselandMean_gg
          read(2,*) Mg2SiO4_rosselandMean_kext
          read(3,*) Mg2SiO4_rosselandMean_pi0
          read(4,*) Mg2SiO4_wav_gg
          read(5,*) Mg2SiO4_wav_kext
          read(6,*) Mg2SiO4_wav_pi0
          read(7,*) Mg2SiO4_PlanckMean_gg
          read(8,*) Mg2SiO4_PlanckMean_kext
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
          open (2, file='../CLOUD_DATA/VO_rosselandMean_kext.txt')
          open (3, file='../CLOUD_DATA/VO_rosselandMean_pi0.txt')
          open (4, file='../CLOUD_DATA/VO_wav_gg.txt')
          open (5, file='../CLOUD_DATA/VO_wav_kext.txt')
          open (6, file='../CLOUD_DATA/VO_wav_pi0.txt')
          open (7, file='../CLOUD_DATA/VO_PlanckMean_gg.txt')
          open (8, file='../CLOUD_DATA/VO_PlanckMean_kext.txt')
          open (9, file='../CLOUD_DATA/VO_PlanckMean_pi0.txt')
          
          read(1,*) VO_rosselandMean_gg
          read(2,*) VO_rosselandMean_kext
          read(3,*) VO_rosselandMean_pi0
          read(4,*) VO_wav_gg
          read(5,*) VO_wav_kext
          read(6,*) VO_wav_pi0
          read(7,*) VO_PlanckMean_gg
          read(8,*) VO_PlanckMean_kext
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
          open (2, file='../CLOUD_DATA/Ni_rosselandMean_kext.txt')
          open (3, file='../CLOUD_DATA/Ni_rosselandMean_pi0.txt')
          open (4, file='../CLOUD_DATA/Ni_wav_gg.txt')
          open (5, file='../CLOUD_DATA/Ni_wav_kext.txt')
          open (6, file='../CLOUD_DATA/Ni_wav_pi0.txt')
          open (7, file='../CLOUD_DATA/Ni_PlanckMean_gg.txt')
          open (8, file='../CLOUD_DATA/Ni_PlanckMean_kext.txt')
          open (9, file='../CLOUD_DATA/Ni_PlanckMean_pi0.txt')
          
          read(1,*) Ni_rosselandMean_gg
          read(2,*) Ni_rosselandMean_kext
          read(3,*) Ni_rosselandMean_pi0
          read(4,*) Ni_wav_gg
          read(5,*) Ni_wav_kext
          read(6,*) Ni_wav_pi0
          read(7,*) Ni_PlanckMean_gg
          read(8,*) Ni_PlanckMean_kext
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
          open (2, file='../CLOUD_DATA/Fe_rosselandMean_kext.txt')
          open (3, file='../CLOUD_DATA/Fe_rosselandMean_pi0.txt')
          open (4, file='../CLOUD_DATA/Fe_wav_gg.txt')
          open (5, file='../CLOUD_DATA/Fe_wav_kext.txt')
          open (6, file='../CLOUD_DATA/Fe_wav_pi0.txt')
          open (7, file='../CLOUD_DATA/Fe_PlanckMean_gg.txt')
          open (8, file='../CLOUD_DATA/Fe_PlanckMean_kext.txt')
          open (9, file='../CLOUD_DATA/Fe_PlanckMean_pi0.txt')
          
          read(1,*) Fe_rosselandMean_gg
          read(2,*) Fe_rosselandMean_kext
          read(3,*) Fe_rosselandMean_pi0
          read(4,*) Fe_wav_gg
          read(5,*) Fe_wav_kext
          read(6,*) Fe_wav_pi0
          read(7,*) Fe_PlanckMean_gg
          read(8,*) Fe_PlanckMean_kext
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
          open (2, file='../CLOUD_DATA/CaSiO4_rosselandMean_kext.txt')
          open (3, file='../CLOUD_DATA/CaSiO4_rosselandMean_pi0.txt')
          open (4, file='../CLOUD_DATA/CaSiO4_wav_gg.txt')
          open (5, file='../CLOUD_DATA/CaSiO4_wav_kext.txt')
          open (6, file='../CLOUD_DATA/CaSiO4_wav_pi0.txt')
          open (7, file='../CLOUD_DATA/CaSiO4_PlanckMean_gg.txt')
          open (8, file='../CLOUD_DATA/CaSiO4_PlanckMean_kext.txt')
          open (9, file='../CLOUD_DATA/CaSiO4_PlanckMean_pi0.txt')
          
          read(1,*) CaSiO4_rosselandMean_gg
          read(2,*) CaSiO4_rosselandMean_kext
          read(3,*) CaSiO4_rosselandMean_pi0
          read(4,*) CaSiO4_wav_gg
          read(5,*) CaSiO4_wav_kext
          read(6,*) CaSiO4_wav_pi0
          read(7,*) CaSiO4_PlanckMean_gg
          read(8,*) CaSiO4_PlanckMean_kext
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
          open (2, file='../CLOUD_DATA/CaTiO3_rosselandMean_kext.txt')
          open (3, file='../CLOUD_DATA/CaTiO3_rosselandMean_pi0.txt')
          open (4, file='../CLOUD_DATA/CaTiO3_wav_gg.txt')
          open (5, file='../CLOUD_DATA/CaTiO3_wav_kext.txt')
          open (6, file='../CLOUD_DATA/CaTiO3_wav_pi0.txt')
          open (7, file='../CLOUD_DATA/CaTiO3_PlanckMean_gg.txt')
          open (8, file='../CLOUD_DATA/CaTiO3_PlanckMean_kext.txt')
          open (9, file='../CLOUD_DATA/CaTiO3_PlanckMean_pi0.txt')
          
          read(1,*) CaTiO3_rosselandMean_gg
          read(2,*) CaTiO3_rosselandMean_kext
          read(3,*) CaTiO3_rosselandMean_pi0
          read(4,*) CaTiO3_wav_gg
          read(5,*) CaTiO3_wav_kext
          read(6,*) CaTiO3_wav_pi0
          read(7,*) CaTiO3_PlanckMean_gg
          read(8,*) CaTiO3_PlanckMean_kext
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
          open (2, file='../CLOUD_DATA/Al2O3_rosselandMean_kext.txt')
          open (3, file='../CLOUD_DATA/Al2O3_rosselandMean_pi0.txt')
          open (4, file='../CLOUD_DATA/Al2O3_wav_gg.txt')
          open (5, file='../CLOUD_DATA/Al2O3_wav_kext.txt')
          open (6, file='../CLOUD_DATA/Al2O3_wav_pi0.txt')
          open (7, file='../CLOUD_DATA/Al2O3_PlanckMean_gg.txt')
          open (8, file='../CLOUD_DATA/Al2O3_PlanckMean_kext.txt')
          open (9, file='../CLOUD_DATA/Al2O3_PlanckMean_pi0.txt')
          
          read(1,*) Al2O3_rosselandMean_gg
          read(2,*) Al2O3_rosselandMean_kext
          read(3,*) Al2O3_rosselandMean_pi0
          read(4,*) Al2O3_wav_gg
          read(5,*) Al2O3_wav_kext
          read(6,*) Al2O3_wav_pi0
          read(7,*) Al2O3_PlanckMean_gg
          read(8,*) Al2O3_PlanckMean_kext
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
          ! Previously in ropprmulti.f
          HAZE_WAV_GRID = (/0.1, 0.101, 0.102, 0.103, 0.104, 0.105, 0.107, 0.108, 0.109, 0.11, 0.111, 0.112, 0.114,
     &                  0.115, 0.116, 0.117, 0.119, 0.12, 0.121, 0.122, 0.124, 0.125, 0.126, 0.128, 0.129, 0.13,
     &                  0.132, 0.133, 0.135, 0.136, 0.138, 0.139, 0.14, 0.142, 0.143, 0.145, 0.147, 0.148, 0.15,
     &                  0.151, 0.153, 0.155, 0.156, 0.158, 0.16, 0.161, 0.163, 0.165, 0.166, 0.168, 0.17, 0.172,
     &                  0.174, 0.176, 0.177, 0.179, 0.181, 0.183, 0.185, 0.187, 0.189, 0.191, 0.193, 0.195, 0.197,
     &                  0.199, 0.202, 0.204, 0.206, 0.208, 0.21, 0.213, 0.215, 0.217, 0.219, 0.222, 0.224, 0.227,
     &                  0.229, 0.231, 0.234, 0.236, 0.239, 0.241, 0.244, 0.247, 0.249, 0.252, 0.255, 0.257, 0.26,
     &                  0.263, 0.266, 0.268, 0.271, 0.274, 0.277, 0.28, 0.283, 0.286, 0.289, 0.292, 0.295, 0.299,
     &                  0.302, 0.305, 0.308, 0.311, 0.315, 0.318, 0.322, 0.325, 0.328, 0.332, 0.335, 0.339, 0.343,
     &                  0.346, 0.35, 0.354, 0.358, 0.361, 0.365, 0.369, 0.373, 0.377, 0.381, 0.385, 0.389, 0.393,
     &                  0.398, 0.402, 0.406, 0.41, 0.415, 0.419, 0.424, 0.428, 0.433, 0.437, 0.442, 0.447, 0.452,
     &                  0.456, 0.461, 0.466, 0.471, 0.476, 0.481, 0.487, 0.492, 0.497, 0.502, 0.508, 0.513, 0.519,
     &                  0.524, 0.53, 0.535, 0.541, 0.547, 0.553, 0.559, 0.564, 0.57, 0.577, 0.583, 0.589, 0.595,
     &                  0.602, 0.608, 0.615, 0.621, 0.628, 0.634, 0.641, 0.648, 0.655, 0.662, 0.669, 0.676, 0.683,
     &                  0.691, 0.698, 0.705, 0.713, 0.721, 0.728, 0.736, 0.744, 0.752, 0.76, 0.768, 0.776, 0.784,
     &                  0.793, 0.801, 0.81, 0.819, 0.827, 0.836, 0.845, 0.854, 0.863, 0.872, 0.882, 0.891, 0.901,
     &                  0.91, 0.92, 0.93, 0.94, 0.95, 0.96, 0.97, 0.98, 0.991, 1.002, 1.012, 1.023, 1.034, 1.045,
     &                  1.056, 1.067, 1.079, 1.09, 1.102, 1.114, 1.126, 1.138, 1.15, 1.162, 1.174, 1.187, 1.2,
     &                  1.212, 1.225, 1.238, 1.252, 1.265, 1.279, 1.292, 1.306, 1.32, 1.334, 1.348, 1.363, 1.377,
     &                  1.392, 1.407, 1.422, 1.437, 1.452, 1.468, 1.483, 1.499, 1.515, 1.531, 1.548, 1.564, 1.581,
     &                  1.598, 1.615, 1.632, 1.65, 1.667, 1.685, 1.703, 1.721, 1.74, 1.758, 1.777, 1.796, 1.815,
     &                  1.834, 1.854, 1.874, 1.894, 1.914, 1.934, 1.955, 1.976, 1.997, 2.018, 2.04, 2.062, 2.084,
     &                  2.106, 2.128, 2.151, 2.174, 2.197, 2.221, 2.244, 2.268, 2.293, 2.317, 2.342, 2.367, 2.392,
     &                  2.418, 2.443, 2.47, 2.496, 2.523, 2.549, 2.577, 2.604, 2.632, 2.66, 2.688, 2.717, 2.746,
     &                  2.775, 2.805, 2.835, 2.865, 2.896, 2.927, 2.958, 2.99, 3.022, 3.054, 3.086, 3.119, 3.153,
     &                  3.186, 3.22, 3.255, 3.289, 3.325, 3.36, 3.396, 3.432, 3.469, 3.506, 3.543, 3.581, 3.619,
     &                  3.658, 3.697, 3.736, 3.776, 3.817, 3.857, 3.899, 3.94, 3.982, 4.025, 4.068, 4.111, 4.155,
     &                  4.199, 4.244, 4.289, 4.335, 4.382, 4.428, 4.476, 4.523, 4.572, 4.62, 4.67, 4.72, 4.77,
     &                  4.821, 4.872, 4.924, 4.977, 5.03, 5.084, 5.138, 5.193, 5.248, 5.304, 5.361, 5.418, 5.476,
     &                  5.534, 5.594, 5.653, 5.714, 5.775, 5.836, 5.898, 5.961, 6.025, 6.089, 6.154, 6.22, 6.286,
     &                  6.354, 6.421, 6.49, 6.559, 6.629, 6.7, 6.772, 6.844, 6.917, 6.991, 7.065, 7.141, 7.217,
     &                  7.294, 7.372, 7.451, 7.53, 7.61, 7.692, 7.774, 7.857, 7.941, 8.025, 8.111, 8.198, 8.285,
     &                  8.374, 8.463, 8.553, 8.645, 8.737, 8.83, 8.924, 9.02, 9.116, 9.213, 9.312, 9.411, 9.511,
     &                  9.613, 9.716, 9.819, 9.924, 10.03, 10.137, 10.245, 10.355, 10.465, 10.577, 10.69, 10.804,
     &                  10.919, 11.036, 11.154, 11.273, 11.393, 11.515, 11.638, 11.762, 11.887, 12.014, 12.143,
     &                  12.272, 12.403, 12.536, 12.669, 12.805, 12.941, 13.079, 13.219, 13.36, 13.503, 13.647,
     &                  13.793, 13.94, 14.089, 14.239, 14.391, 14.545, 14.7, 14.857, 15.015, 15.176, 15.338,
     &                  15.501, 15.667, 15.834, 16.003, 16.174, 16.347, 16.521, 16.697, 16.876, 17.056, 17.238,
     &                  17.422, 17.608, 17.796, 17.986, 18.178, 18.372, 18.568, 18.766, 18.966, 19.169, 19.373,
     &                  19.58, 19.789, 20.0/)
          CLOUD_WAV_GRID = (/ 
     &                  0.1       ,  0.10549764,  0.11129751,  0.11741624,  0.12387136, 
     &                  0.13068136,  0.13786574,  0.1454451 ,  0.15344114,  0.16187678, 
     &                  0.17077617,  0.18016482,  0.19006963,  0.20051896,  0.21154277, 
     &                  0.22317262,  0.23544184,  0.24838557,  0.2620409 ,  0.27644696, 
     &                  0.29164501,  0.30767859,  0.32459363,  0.34243861,  0.36126464, 
     &                  0.38112565,  0.40207855,  0.42418337,  0.44750342,  0.47210553, 
     &                  0.49806018,  0.52544171,  0.55432858,  0.58480355,  0.61695392, 
     &                  0.6508718 ,  0.68665436,  0.72440411,  0.76422921,  0.80624375, 
     &                  0.8505681 ,  0.89732923,  0.94666113,  0.99870511,  1.05361028, 
     &                  1.11153393,  1.17264202,  1.23710961,  1.30512139,  1.37687221, 
     &                  1.45256763,  1.53242451,  1.61667163,  1.70555034,  1.79931529, 
     &                  1.89823509,  2.00259314,  2.11268842,  2.22883634,  2.35136964, 
     &                  2.48063938,  2.6170159 ,  2.7608899 ,  2.91267357,  3.07280176, 
     &                  3.24173321,  3.41995189,  3.60796839,  3.80632136,  4.01557904, 
     &                  4.23634095,  4.46923955,  4.71494206,  4.97415241,  5.24761319, 
     &                  5.53610785,  5.8404629 ,  6.16155028,  6.50028987,  6.85765214, 
     &                  7.23466087,  7.63239618,  8.05199753,  8.49466702,  8.96167288, 
     &                  9.45435302,  9.97411891, 10.52245965, 11.10094616, 11.71123575, 
     &                  12.35507684, 13.03431396, 13.75089307, 14.5068671 , 15.30440182,
     &                  16.14578209, 17.03341839, 17.96985369, 18.9577708 , 20.0
     &                  /)

          input_pressure_array_cgs = (/1.000e+00, 1.263e+00, 1.594e+00, 2.013e+00, 2.541e+00, 
     &    3.209e+00, 4.051e+00, 5.115e+00, 6.458e+00, 8.154e+00, 
     &    1.030e+01, 1.300e+01, 1.641e+01, 2.072e+01, 2.617e+01, 
     &    3.304e+01, 4.171e+01, 5.266e+01, 6.649e+01, 8.396e+01, 
     &    1.060e+02, 1.338e+02, 1.690e+02, 2.134e+02, 2.694e+02, 
     &    3.401e+02, 4.294e+02, 5.422e+02, 6.846e+02, 8.644e+02, 
     &    1.091e+03, 1.378e+03, 1.740e+03, 2.197e+03, 2.774e+03, 
     &    3.502e+03, 4.421e+03, 5.583e+03, 7.049e+03, 8.900e+03, 
     &    1.124e+04, 1.419e+04, 1.791e+04, 2.262e+04, 2.856e+04, 
     &    3.605e+04, 4.552e+04, 5.748e+04, 7.257e+04, 9.163e+04, 
     &    1.157e+05, 1.461e+05, 1.844e+05, 2.329e+05, 2.940e+05, 
     &    3.712e+05, 4.687e+05, 5.918e+05, 7.472e+05, 9.434e+05, 
     &    1.191e+06, 1.504e+06, 1.899e+06, 2.397e+06, 3.027e+06, 
     &    3.822e+06, 4.826e+06, 6.093e+06, 7.693e+06, 9.713e+06, 
     &    1.226e+07, 1.548e+07, 1.955e+07, 2.468e+07, 3.117e+07, 
     &    3.935e+07, 4.968e+07, 6.273e+07, 7.920e+07, 1.000e+08/)

          input_particle_size_array_in_meters = (/
     &    1.00000000e-09, 1.12332403e-09, 1.26185688e-09, 1.41747416e-09,
     &    1.59228279e-09, 1.78864953e-09, 2.00923300e-09, 2.25701972e-09,
     &    2.53536449e-09, 2.84803587e-09, 3.19926714e-09, 3.59381366e-09,
     &    4.03701726e-09, 4.53487851e-09, 5.09413801e-09, 5.72236766e-09,
     &    6.42807312e-09, 7.22080902e-09, 8.11130831e-09, 9.11162756e-09,
     &    1.02353102e-08, 1.14975700e-08, 1.29154967e-08, 1.45082878e-08,
     &    1.62975083e-08, 1.83073828e-08, 2.05651231e-08, 2.31012970e-08,
     &    2.59502421e-08, 2.91505306e-08, 3.27454916e-08, 3.67837977e-08,
     &    4.13201240e-08, 4.64158883e-08, 5.21400829e-08, 5.85702082e-08,
     &    6.57933225e-08, 7.39072203e-08, 8.30217568e-08, 9.32603347e-08,
     &    1.04761575e-07, 1.17681195e-07, 1.32194115e-07, 1.48496826e-07,
     &    1.66810054e-07, 1.87381742e-07, 2.10490414e-07, 2.36448941e-07,
     &    2.65608778e-07, 2.98364724e-07, 3.35160265e-07, 3.76493581e-07,
     &    4.22924287e-07, 4.75081016e-07, 5.33669923e-07, 5.99484250e-07,
     &    6.73415066e-07, 7.56463328e-07, 8.49753436e-07, 9.54548457e-07,
     &    1.07226722e-06, 1.20450354e-06, 1.35304777e-06, 1.51991108e-06,
     &    1.70735265e-06, 1.91791026e-06, 2.15443469e-06, 2.42012826e-06,
     &    2.71858824e-06, 3.05385551e-06, 3.43046929e-06, 3.85352859e-06,
     &    4.32876128e-06, 4.86260158e-06, 5.46227722e-06, 6.13590727e-06,
     &    6.89261210e-06, 7.74263683e-06, 8.69749003e-06, 9.77009957e-06,
     &    1.09749877e-05, 1.23284674e-05, 1.38488637e-05, 1.55567614e-05,
     &    1.74752840e-05, 1.96304065e-05, 2.20513074e-05, 2.47707636e-05,
     &    2.78255940e-05, 3.12571585e-05, 3.51119173e-05, 3.94420606e-05,
     &    4.43062146e-05, 4.97702356e-05, 5.59081018e-05, 6.28029144e-05,
     &    7.05480231e-05, 7.92482898e-05, 8.90215085e-05, 1.00000000e-04
     &    /)

          input_temperature_array = (/
     &    100.0, 139.39393939393938, 178.78787878787878, 218.1818181818182, 257.57575757575756,
     &    296.96969696969694, 336.3636363636364, 375.75757575757575, 415.1515151515151, 454.5454545454545,
     &    493.9393939393939, 533.3333333333333, 572.7272727272727, 612.1212121212121, 651.5151515151515,
     &    690.9090909090909, 730.3030303030303, 769.6969696969696, 809.090909090909, 848.4848484848484,
     &    887.8787878787878, 927.2727272727273, 966.6666666666666, 1006.060606060606, 1045.4545454545455,
     &    1084.8484848484848, 1124.2424242424242, 1163.6363636363635, 1203.030303030303, 1242.4242424242423,
     &    1281.8181818181818, 1321.212121212121, 1360.6060606060605, 1400.0, 1439.3939393939393,
     &    1478.7878787878788, 1518.181818181818, 1557.5757575757575, 1596.9696969696968, 1636.3636363636363,
     &    1675.7575757575755, 1715.151515151515, 1754.5454545454545, 1793.9393939393938, 1833.3333333333333,
     &    1872.7272727272725, 1912.121212121212, 1951.5151515151513, 1990.9090909090908, 2030.3030303030303,
     &    2069.6969696969695, 2109.090909090909, 2148.4848484848485, 2187.8787878787875, 2227.272727272727,
     &    2266.6666666666665, 2306.060606060606, 2345.4545454545455, 2384.8484848484845, 2424.242424242424,
     &    2463.6363636363635, 2503.030303030303, 2542.424242424242, 2581.8181818181815, 2621.212121212121,
     &    2660.6060606060605, 2700.0, 2739.393939393939, 2778.7878787878785, 2818.181818181818,
     &    2857.5757575757575, 2896.9696969696965, 2936.363636363636, 2975.7575757575755, 3015.151515151515,
     &    3054.5454545454545, 3093.9393939393935, 3133.333333333333, 3172.7272727272725, 3212.121212121212,
     &    3251.515151515151, 3290.9090909090905, 3330.30303030303, 3369.6969696969695, 3409.090909090909,
     &    3448.484848484848, 3487.8787878787875, 3527.272727272727, 3566.6666666666665, 3606.060606060606,
     &    3645.454545454545, 3684.8484848484845, 3724.242424242424, 3763.6363636363635, 3803.0303030303025,
     &    3842.424242424242, 3881.8181818181815, 3921.212121212121, 3960.6060606060605, 4000.0
     &    /)

          particle_size_vs_layer_array_in_meters = (/1.000e-07, 1.000e-07, 1.000e-07, 1.000e-07, 1.000e-07, 
     &    1.000e-07, 1.000e-07, 1.000e-07, 1.000e-07, 1.000e-07, 
     &    1.000e-07, 1.000e-07, 1.000e-07, 1.000e-07, 1.000e-07, 
     &    1.000e-07, 1.000e-07, 1.000e-07, 1.000e-07, 1.000e-07, 
     &    1.000e-07, 1.000e-07, 1.000e-07, 1.000e-07, 1.000e-07, 
     &    1.000e-07, 1.000e-07, 1.000e-07, 1.000e-07, 1.000e-07, 
     &    1.000e-07, 1.000e-07, 1.000e-07, 1.000e-07, 1.000e-07, 
     &    1.000e-07, 1.000e-07, 1.000e-07, 1.000e-07, 1.000e-07, 
     &    1.012e-07, 1.042e-07, 1.079e-07, 1.126e-07, 1.185e-07, 
     &    1.260e-07, 1.355e-07, 1.474e-07, 1.625e-07, 1.816e-07, 
     &    2.056e-07, 2.359e-07, 2.743e-07, 3.227e-07, 3.837e-07, 
     &    4.609e-07, 5.583e-07, 6.812e-07, 8.365e-07, 1.033e-06, 
     &    1.280e-06, 1.593e-06, 1.987e-06, 2.485e-06, 3.114e-06, 
     &    3.908e-06, 4.911e-06, 6.177e-06, 7.776e-06, 9.794e-06, 
     &    1.234e-05, 1.556e-05, 1.962e-05, 2.475e-05, 3.123e-05, 
     &    3.940e-05, 4.973e-05, 6.276e-05, 7.922e-05, 1.000e-04/)

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

      CORFACT =   (/1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,
     &              1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,
     &              1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,
     &              1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,
     &              1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,
     &              1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,
     &              1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,
     &              1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,
     &              1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,
     &              1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000/)

      ! The condensation curves get read from file now instead of being hardcoded. in CLOUD_DATA, there should be both the orig (50 pressures, from Eliza Kempton's group) and new (80 pressures, interpolated) files for each metallicity.
      ! The first column is pressure (cgs), then each cloud species gets its own column in the usual order
      open(UNIT=10, FILE='../CLOUD_DATA/condcurves_1xsolar_new.txt')
      read(10,*) dummy_tconds
      close(10)
      ! write(*,*) dummy_tconds(2,:)
      do i=1, 13
        do J=1, 80
            TCONDS(1,J,I) = dummy_tconds(I+1,J)
        end do
      end do
      open(UNIT=11, FILE='../CLOUD_DATA/condcurves_10xsolar_new.txt')
      read(11,*) dummy_tconds
      close(11)
      do i=1, 13
        do J=1, 80
            TCONDS(2,J,I) = dummy_tconds(I+1,J)
        end do
      end do
      open(UNIT=12, FILE='../CLOUD_DATA/condcurves_30xsolar_new.txt')
      read(12,*) dummy_tconds
      close(12)
      do i=1, 13
        do J=1, 80
            TCONDS(3,J,I) = dummy_tconds(I+1,J)
        end do
      end do
      open(UNIT=13, FILE='../CLOUD_DATA/condcurves_100xsolar_new.txt')
      read(13,*) dummy_tconds
      close(13)
      do i=1, 13
        do J=1, 80
            TCONDS(4,J,I) = dummy_tconds(I+1,J)
        end do
      end do
      open(UNIT=14, FILE='../CLOUD_DATA/condcurves_300xsolar_new.txt')
      read(14,*) dummy_tconds
      close(14)
      do i=1, 13
        do J=1, 80
            TCONDS(5,J,I) = dummy_tconds(I+1,J)
        end do
      end do

      G0_OPPR(1,1:100,1:100,1)=KCl_wav_gg
      G0_OPPR(2,1:100,1:100,1)=KCl_wav_gg
      G0_OPPR(3,1:100,1:100,1)=KCl_wav_gg
      G0_OPPR(4,1:100,1:100,1)=KCl_PlanckMean_gg
      G0_OPPR(5,1:100,1:100,1)=KCl_rosselandMean_gg

      G0_OPPR(1,1:100,1:100,2)=ZnS_wav_gg
      G0_OPPR(2,1:100,1:100,2)=ZnS_wav_gg
      G0_OPPR(3,1:100,1:100,2)=ZnS_wav_gg
      G0_OPPR(4,1:100,1:100,2)=ZnS_PlanckMean_gg
      G0_OPPR(5,1:100,1:100,2)=ZnS_rosselandMean_gg

      G0_OPPR(1,1:100,1:100,3)=Na2S_wav_gg
      G0_OPPR(2,1:100,1:100,3)=Na2S_wav_gg
      G0_OPPR(3,1:100,1:100,3)=Na2S_wav_gg
      G0_OPPR(4,1:100,1:100,3)=Na2S_PlanckMean_gg
      G0_OPPR(5,1:100,1:100,3)=Na2S_rosselandMean_gg

      G0_OPPR(1,1:100,1:100,4)=MnS_wav_gg
      G0_OPPR(2,1:100,1:100,4)=MnS_wav_gg
      G0_OPPR(3,1:100,1:100,4)=MnS_wav_gg
      G0_OPPR(4,1:100,1:100,4)=MnS_PlanckMean_gg
      G0_OPPR(5,1:100,1:100,4)=MnS_rosselandMean_gg

      G0_OPPR(1,1:100,1:100,5)=Cr_wav_gg
      G0_OPPR(2,1:100,1:100,5)=Cr_wav_gg
      G0_OPPR(3,1:100,1:100,5)=Cr_wav_gg
      G0_OPPR(4,1:100,1:100,5)=Cr_PlanckMean_gg
      G0_OPPR(5,1:100,1:100,5)=Cr_rosselandMean_gg

      G0_OPPR(1,1:100,1:100,6)=SiO2_wav_gg
      G0_OPPR(2,1:100,1:100,6)=SiO2_wav_gg
      G0_OPPR(3,1:100,1:100,6)=SiO2_wav_gg
      G0_OPPR(4,1:100,1:100,6)=SiO2_PlanckMean_gg
      G0_OPPR(5,1:100,1:100,6)=SiO2_rosselandMean_gg

      G0_OPPR(1,1:100,1:100,7)=Mg2SiO4_wav_gg
      G0_OPPR(2,1:100,1:100,7)=Mg2SiO4_wav_gg
      G0_OPPR(3,1:100,1:100,7)=Mg2SiO4_wav_gg
      G0_OPPR(4,1:100,1:100,7)=Mg2SiO4_PlanckMean_gg
      G0_OPPR(5,1:100,1:100,7)=Mg2SiO4_rosselandMean_gg

      G0_OPPR(1,1:100,1:100,8)=VO_wav_gg
      G0_OPPR(2,1:100,1:100,8)=VO_wav_gg
      G0_OPPR(3,1:100,1:100,8)=VO_wav_gg
      G0_OPPR(4,1:100,1:100,8)=VO_PlanckMean_gg
      G0_OPPR(5,1:100,1:100,8)=VO_rosselandMean_gg
      
      G0_OPPR(1,1:100,1:100,9)=Ni_wav_gg
      G0_OPPR(2,1:100,1:100,9)=Ni_wav_gg
      G0_OPPR(3,1:100,1:100,9)=Ni_wav_gg
      G0_OPPR(4,1:100,1:100,9)=Ni_PlanckMean_gg
      G0_OPPR(5,1:100,1:100,9)=Ni_rosselandMean_gg
      
      G0_OPPR(1,1:100,1:100,10)=Fe_wav_gg
      G0_OPPR(2,1:100,1:100,10)=Fe_wav_gg
      G0_OPPR(3,1:100,1:100,10)=Fe_wav_gg
      G0_OPPR(4,1:100,1:100,10)=Fe_PlanckMean_gg
      G0_OPPR(5,1:100,1:100,10)=Fe_rosselandMean_gg
      
      G0_OPPR(1,1:100,1:100,11)=CaSiO4_wav_gg
      G0_OPPR(2,1:100,1:100,11)=CaSiO4_wav_gg
      G0_OPPR(3,1:100,1:100,11)=CaSiO4_wav_gg
      G0_OPPR(4,1:100,1:100,11)=CaSiO4_PlanckMean_gg
      G0_OPPR(5,1:100,1:100,11)=CaSiO4_rosselandMean_gg
      
      G0_OPPR(1,1:100,1:100,12)=CaTiO3_wav_gg
      G0_OPPR(2,1:100,1:100,12)=CaTiO3_wav_gg
      G0_OPPR(3,1:100,1:100,12)=CaTiO3_wav_gg
      G0_OPPR(4,1:100,1:100,12)=CaTiO3_PlanckMean_gg
      G0_OPPR(5,1:100,1:100,12)=CaTiO3_rosselandMean_gg
      
      G0_OPPR(1,1:100,1:100,13)=Al2O3_wav_gg
      G0_OPPR(2,1:100,1:100,13)=Al2O3_wav_gg
      G0_OPPR(3,1:100,1:100,13)=Al2O3_wav_gg
      G0_OPPR(4,1:100,1:100,13)=Al2O3_PlanckMean_gg
      G0_OPPR(5,1:100,1:100,13)=Al2O3_rosselandMean_gg
      
      PI0_OPPR(1,1:100,1:100,1)=KCl_wav_pi0
      PI0_OPPR(2,1:100,1:100,1)=KCl_wav_pi0
      PI0_OPPR(3,1:100,1:100,1)=KCl_wav_pi0
      PI0_OPPR(4,1:100,1:100,1)=KCl_PlanckMean_pi0
      PI0_OPPR(5,1:100,1:100,1)=KCl_rosselandMean_pi0
      
      PI0_OPPR(1,1:100,1:100,2)=ZnS_wav_pi0
      PI0_OPPR(2,1:100,1:100,2)=ZnS_wav_pi0
      PI0_OPPR(3,1:100,1:100,2)=ZnS_wav_pi0
      PI0_OPPR(4,1:100,1:100,2)=ZnS_PlanckMean_pi0
      PI0_OPPR(5,1:100,1:100,2)=ZnS_rosselandMean_pi0
      
      PI0_OPPR(1,1:100,1:100,3)=Na2S_wav_pi0
      PI0_OPPR(2,1:100,1:100,3)=Na2S_wav_pi0
      PI0_OPPR(3,1:100,1:100,3)=Na2S_wav_pi0
      PI0_OPPR(4,1:100,1:100,3)=Na2S_PlanckMean_pi0
      PI0_OPPR(5,1:100,1:100,3)=Na2S_rosselandMean_pi0
      
      PI0_OPPR(1,1:100,1:100,4)=MnS_wav_pi0
      PI0_OPPR(2,1:100,1:100,4)=MnS_wav_pi0
      PI0_OPPR(3,1:100,1:100,4)=MnS_wav_pi0
      PI0_OPPR(4,1:100,1:100,4)=MnS_PlanckMean_pi0
      PI0_OPPR(5,1:100,1:100,4)=MnS_rosselandMean_pi0
      
      PI0_OPPR(1,1:100,1:100,5)=Cr_wav_pi0
      PI0_OPPR(2,1:100,1:100,5)=Cr_wav_pi0
      PI0_OPPR(3,1:100,1:100,5)=Cr_wav_pi0
      PI0_OPPR(4,1:100,1:100,5)=Cr_PlanckMean_pi0
      PI0_OPPR(5,1:100,1:100,5)=Cr_rosselandMean_pi0

      PI0_OPPR(1,1:100,1:100,6)=SiO2_wav_pi0
      PI0_OPPR(2,1:100,1:100,6)=SiO2_wav_pi0
      PI0_OPPR(3,1:100,1:100,6)=SiO2_wav_pi0
      PI0_OPPR(4,1:100,1:100,6)=SiO2_PlanckMean_pi0
      PI0_OPPR(5,1:100,1:100,6)=SiO2_rosselandMean_pi0

      PI0_OPPR(1,1:100,1:100,7)=Mg2SiO4_wav_pi0
      PI0_OPPR(2,1:100,1:100,7)=Mg2SiO4_wav_pi0
      PI0_OPPR(3,1:100,1:100,7)=Mg2SiO4_wav_pi0
      PI0_OPPR(4,1:100,1:100,7)=Mg2SiO4_PlanckMean_pi0
      PI0_OPPR(5,1:100,1:100,7)=Mg2SiO4_rosselandMean_pi0

      PI0_OPPR(1,1:100,1:100,8)=VO_wav_pi0
      PI0_OPPR(2,1:100,1:100,8)=VO_wav_pi0
      PI0_OPPR(3,1:100,1:100,8)=VO_wav_pi0
      PI0_OPPR(4,1:100,1:100,8)=VO_PlanckMean_pi0
      PI0_OPPR(5,1:100,1:100,8)=VO_rosselandMean_pi0
      
      PI0_OPPR(1,1:100,1:100,9)=Ni_wav_pi0
      PI0_OPPR(2,1:100,1:100,9)=Ni_wav_pi0
      PI0_OPPR(3,1:100,1:100,9)=Ni_wav_pi0
      PI0_OPPR(4,1:100,1:100,9)=Ni_PlanckMean_pi0
      PI0_OPPR(5,1:100,1:100,9)=Ni_rosselandMean_pi0
      
      PI0_OPPR(1,1:100,1:100,10)=Fe_wav_pi0
      PI0_OPPR(2,1:100,1:100,10)=Fe_wav_pi0
      PI0_OPPR(3,1:100,1:100,10)=Fe_wav_pi0
      PI0_OPPR(4,1:100,1:100,10)=Fe_PlanckMean_pi0
      PI0_OPPR(5,1:100,1:100,10)=Fe_rosselandMean_pi0
      
      PI0_OPPR(1,1:100,1:100,11)=CaSiO4_wav_pi0
      PI0_OPPR(2,1:100,1:100,11)=CaSiO4_wav_pi0
      PI0_OPPR(3,1:100,1:100,11)=CaSiO4_wav_pi0
      PI0_OPPR(4,1:100,1:100,11)=CaSiO4_PlanckMean_pi0
      PI0_OPPR(5,1:100,1:100,11)=CaSiO4_rosselandMean_pi0
      
      PI0_OPPR(1,1:100,1:100,12)=CaTiO3_wav_pi0
      PI0_OPPR(2,1:100,1:100,12)=CaTiO3_wav_pi0
      PI0_OPPR(3,1:100,1:100,12)=CaTiO3_wav_pi0
      PI0_OPPR(4,1:100,1:100,12)=CaTiO3_PlanckMean_pi0
      PI0_OPPR(5,1:100,1:100,12)=CaTiO3_rosselandMean_pi0
      
      PI0_OPPR(1,1:100,1:100,13)=Al2O3_wav_pi0
      PI0_OPPR(2,1:100,1:100,13)=Al2O3_wav_pi0
      PI0_OPPR(3,1:100,1:100,13)=Al2O3_wav_pi0
      PI0_OPPR(4,1:100,1:100,13)=Al2O3_PlanckMean_pi0
      PI0_OPPR(5,1:100,1:100,13)=Al2O3_rosselandMean_pi0
      
      KE_OPPR(1,1:100,1:100,1)=KCl_wav_kext
      KE_OPPR(2,1:100,1:100,1)=KCl_wav_kext
      KE_OPPR(3,1:100,1:100,1)=KCl_wav_kext
      KE_OPPR(4,1:100,1:100,1)=KCl_PlanckMean_kext
      KE_OPPR(5,1:100,1:100,1)=KCl_rosselandMean_kext
      
      KE_OPPR(1,1:100,1:100,2)=ZnS_wav_kext
      KE_OPPR(2,1:100,1:100,2)=ZnS_wav_kext
      KE_OPPR(3,1:100,1:100,2)=ZnS_wav_kext
      KE_OPPR(4,1:100,1:100,2)=ZnS_PlanckMean_kext
      KE_OPPR(5,1:100,1:100,2)=ZnS_rosselandMean_kext
      
      KE_OPPR(1,1:100,1:100,3)=Na2S_wav_kext
      KE_OPPR(2,1:100,1:100,3)=Na2S_wav_kext
      KE_OPPR(3,1:100,1:100,3)=Na2S_wav_kext
      KE_OPPR(4,1:100,1:100,3)=Na2S_PlanckMean_kext
      KE_OPPR(5,1:100,1:100,3)=Na2S_rosselandMean_kext
      
      KE_OPPR(1,1:100,1:100,4)=MnS_wav_kext
      KE_OPPR(2,1:100,1:100,4)=MnS_wav_kext
      KE_OPPR(3,1:100,1:100,4)=MnS_wav_kext
      KE_OPPR(4,1:100,1:100,4)=MnS_PlanckMean_kext
      KE_OPPR(5,1:100,1:100,4)=MnS_rosselandMean_kext
      
      KE_OPPR(1,1:100,1:100,5)=Cr_wav_kext
      KE_OPPR(2,1:100,1:100,5)=Cr_wav_kext
      KE_OPPR(3,1:100,1:100,5)=Cr_wav_kext
      KE_OPPR(4,1:100,1:100,5)=Cr_PlanckMean_kext
      KE_OPPR(5,1:100,1:100,5)=Cr_rosselandMean_kext
      
      KE_OPPR(1,1:100,1:100,6)=SiO2_wav_kext
      KE_OPPR(2,1:100,1:100,6)=SiO2_wav_kext
      KE_OPPR(3,1:100,1:100,6)=SiO2_wav_kext
      KE_OPPR(4,1:100,1:100,6)=SiO2_PlanckMean_kext
      KE_OPPR(5,1:100,1:100,6)=SiO2_rosselandMean_kext

      KE_OPPR(1,1:100,1:100,7)=Mg2SiO4_wav_kext
      KE_OPPR(2,1:100,1:100,7)=Mg2SiO4_wav_kext
      KE_OPPR(3,1:100,1:100,7)=Mg2SiO4_wav_kext
      KE_OPPR(4,1:100,1:100,7)=Mg2SiO4_PlanckMean_kext
      KE_OPPR(5,1:100,1:100,7)=Mg2SiO4_rosselandMean_kext
      
      KE_OPPR(1,1:100,1:100,8)=VO_wav_kext
      KE_OPPR(2,1:100,1:100,8)=VO_wav_kext
      KE_OPPR(3,1:100,1:100,8)=VO_wav_kext
      KE_OPPR(4,1:100,1:100,8)=VO_PlanckMean_kext
      KE_OPPR(5,1:100,1:100,8)=VO_rosselandMean_kext
      
      KE_OPPR(1,1:100,1:100,9)=Ni_wav_kext
      KE_OPPR(2,1:100,1:100,9)=Ni_wav_kext
      KE_OPPR(3,1:100,1:100,9)=Ni_wav_kext
      KE_OPPR(4,1:100,1:100,9)=Ni_PlanckMean_kext
      KE_OPPR(5,1:100,1:100,9)=Ni_rosselandMean_kext
      
      KE_OPPR(1,1:100,1:100,10)=Fe_wav_kext
      KE_OPPR(2,1:100,1:100,10)=Fe_wav_kext
      KE_OPPR(3,1:100,1:100,10)=Fe_wav_kext
      KE_OPPR(4,1:100,1:100,10)=Fe_PlanckMean_kext
      KE_OPPR(5,1:100,1:100,10)=Fe_rosselandMean_kext
      
      KE_OPPR(1,1:100,1:100,11)=CaSiO4_wav_kext
      KE_OPPR(2,1:100,1:100,11)=CaSiO4_wav_kext
      KE_OPPR(3,1:100,1:100,11)=CaSiO4_wav_kext
      KE_OPPR(4,1:100,1:100,11)=CaSiO4_PlanckMean_kext
      KE_OPPR(5,1:100,1:100,11)=CaSiO4_rosselandMean_kext
      
      KE_OPPR(1,1:100,1:100,12)=CaTiO3_wav_kext
      KE_OPPR(2,1:100,1:100,12)=CaTiO3_wav_kext
      KE_OPPR(3,1:100,1:100,12)=CaTiO3_wav_kext
      KE_OPPR(4,1:100,1:100,12)=CaTiO3_PlanckMean_kext
      KE_OPPR(5,1:100,1:100,12)=CaTiO3_rosselandMean_kext
      
      KE_OPPR(1,1:100,1:100,13)=Al2O3_wav_kext
      KE_OPPR(2,1:100,1:100,13)=Al2O3_wav_kext
      KE_OPPR(3,1:100,1:100,13)=Al2O3_wav_kext
      KE_OPPR(4,1:100,1:100,13)=Al2O3_PlanckMean_kext
      KE_OPPR(5,1:100,1:100,13)=Al2O3_rosselandMean_kext

      ! Setting up correction factor for mean vs median particle volume:
      sigma = 2.0 ! should eventually not be hardcoded
      exp_92_lnsig2_pi = EXP(-9.0/2.0 * LOG(sigma)*LOG(sigma)) / (4.D0*DATAN(1.D0))
      
      ! Thomas interpolating cloud condensation curves in metallicity
      ! Temporarily define MET_INDEX (this is just to determine whether we need to interpolate)
      WRITE(*,*), 'METALLICITY = ', METALLICITY
      IF (METALLICITY .gt. -0.1 .AND. METALLICITY .lt. 0.1) THEN
        MET_INDEX = 1
      ELSE IF (METALLICITY .gt. 0.9 .AND. METALLICITY .lt. 1.1) THEN
        MET_INDEX = 2
      ELSE IF (METALLICITY .gt. 1.4 .AND. METALLICITY .lt. 1.6) THEN
        MET_INDEX = 3
      ELSE IF (METALLICITY .gt. 1.9 .AND. METALLICITY .lt. 2.1) THEN
        MET_INDEX = 4
      ELSE IF (METALLICITY .gt. 2.37 .AND. METALLICITY .lt. 2.57) THEN
        MET_INDEX = 5
      ELSE
        ! If metallicity isnt close to one of the pre-calculated curves, interpolate between nearest 2
        COND_CURVE_METS = (/0.0, 1.0, 1.5, 2.0, 2.5/) ! these are log10(Z/Zsun)
        If (METALLICITY .gt. 0.01 .AND. METALLICITY .lt. 0.99) THEN
            LOW_MET_INDEX = 1
            HIGH_MET_INDEX = 2
        ELSE IF (METALLICITY .gt. 1.01	.AND. METALLICITY .lt. 1.49) THEN
            LOW_MET_INDEX	= 2
            HIGH_MET_INDEX = 3
        ELSE IF (METALLICITY .gt. 1.51 .AND. METALLICITY .lt. 1.99) THEN
            LOW_MET_INDEX = 3 
            HIGH_MET_INDEX = 4
        ELSE IF (METALLICITY .gt. 2.01 .AND. METALLICITY .lt. 2.49) THEN
            LOW_MET_INDEX = 4
            HIGH_MET_INDEX = 5
        ELSE
            WRITE(*,*) 'something went wrong with the metallicity, look at cloud_properties_set_up.f'
            stop
        END IF
        MET_INDEX = 6
        DO J = 1, 80
            DO I = 1, NCLOUDS
                ! WRITE(*,*) COND_CURVE_METS(LOW_MET_INDEX), COND_CURVE_METS(HIGH_MET_INDEX)
                !TCONDS(6,J,I) = TCONDS(LOW_MET_INDEX,J,I) + ((METALLICITY - (COND_CURVE_METS(LOW_MET_INDEX))) * (TCONDS(HIGH_MET_INDEX,J,I) - TCONDS(LOW_MET_INDEX,J,I)) / (COND_CURVE_METS(HIGH_MET_INDEX) - COND_CURVE_METS(LOW_MET_INDEX)))
                TCONDS(MET_INDEX,J,I) = (METALLICITY - (COND_CURVE_METS(LOW_MET_INDEX)))
                TCONDS(MET_INDEX,J,I) = TCONDS(MET_INDEX,J,I) * (TCONDS(HIGH_MET_INDEX,J,I) - TCONDS(LOW_MET_INDEX,J,I))
                TCONDS(MET_INDEX,J,I) = TCONDS(MET_INDEX,J,I) / (COND_CURVE_METS(HIGH_MET_INDEX) - COND_CURVE_METS(LOW_MET_INDEX))
                TCONDS(MET_INDEX,J,I) = TCONDS(MET_INDEX,J,I) + TCONDS(LOW_MET_INDEX,J,I)
            END DO
        END DO
        
        WRITE(*,*) 'interpolating condensation curves between ', 
     & COND_CURVE_METS(LOW_MET_INDEX), ' and ', COND_CURVE_METS(HIGH_MET_INDEX)
      ENDIF

      END SUBROUTINE get_cloud_scattering_properties
