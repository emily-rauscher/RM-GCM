      SUBROUTINE get_cloud_scattering_properties_wrapper
          include 'rcommons.h'
          call get_cloud_scattering_properties(NCLOUDS, NLAYER, NVERT, NIRP, NSOLP, GASCON, METALLICITY, HAZETYPE)
      END SUBROUTINE get_cloud_scattering_properties_wrapper


      SUBROUTINE get_cloud_scattering_properties(NCLOUDS, NLAYER, NVERT, NIRP, NSOLP, GASCON, METALLICITY, HAZETYPE)
          implicit none
          integer :: J, L, K, NL, NCLOUDS, NLAYER, NVERT, NIRP, NSOLP
          real :: GAS_CONSTANT_R, GASCON, METALLICITY
          character (len = 40) :: HAZETYPE

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
          REAL TCONDS(4, 51, 13)

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

          !haze_type = 'soot'

          !EDIT —— I think changed?
          if (HAZETYPE .eq. 'soot') THEN
              !write(*,*) "Model being run with soot hazes"x
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
          else if (HAZETYPE .eq. 'sulfur') THEN
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
          else if (HAZETYPE .eq. 'tholin') THEN
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
          else if (HAZETYPE .eq. 'soot-2xpi0') THEN
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




      tcon_10X_MET_Al2O3 = (/1793.841,1796.732,1804.399,1816.024,1830.038,1844.912,
     &                     1859.187,1872.39,1886.825,1902.101,1916.941,1931.308,
     &                     1946.999,1961.776,1977.082,1993.228,2008.809,2023.756,
     &                     2039.951,2055.79,2072.04,2089.899,2107.594,2124.687,
     &                     2142.686,2160.612,2178.37,2198.123,2217.72,2236.436,
     &                     2255.616,2274.932,2292.872,2311.943,2330.163,2347.035,
     &                     2364.23,2382.746,2400.307,2420.805,2441.873,2462.468,
     &                     2483.483,2506.415,2528.454,2558.228,2591.005,2623.478,
     &                     2654.402,2684.311/)

      tcon_10X_MET_CaSiO4 = (/1578.882,1581.914,1589.953,1602.143,1616.839,1632.441,
     &                     1647.42,1661.284,1676.442,1692.495,1708.12,1723.303,
     &                     1739.954,1755.767,1772.495,1790.529,1808.232,1825.385,
     &                     1843.828,1861.794,1880.19,1900.359,1920.324,1939.627,
     &                     1960.022,1980.441,2000.855,2023.816,2046.876,2069.208,
     &                     2092.352,2116.027,2139.028,2165.389,2192.095,2217.908,
     &                     2244.176,2271.595,2297.408,2327.314,2357.72,2386.981,
     &                     2416.295,2447.674,2476.549,2511.904,2548.528,2583.226,
     &                     2615.741,2648.271/)

      tcon_10X_MET_Cr = (/1289.595,1291.64,1297.07,1305.303,1315.221,1325.734,
     &                     1335.799,1345.071,1355.092,1365.595,1375.731,1385.532,
     &                     1396.327,1406.573,1417.381,1428.997,1440.37,1451.366,
     &                     1463.188,1474.689,1486.407,1499.18,1511.776,1523.935,
     &                     1536.819,1549.773,1562.752,1577.354,1592.136,1606.702,
     &                     1622.176,1638.587,1655.816,1677.807,1701.849,1726.417,
     &                     1751.716,1777.751,1803.065,1834.231,1867.236,1899.88,
     &                     1932.599,1966.668,1997.567,2034.066,2070.946,2105.206,
     &                     2137.074,2169.45/)

      tcon_10X_MET_Fe = (/1359.974,1363.324,1372.197,1385.654,1401.882,1419.124,
     &                     1435.699,1451.074,1467.958,1485.915,1503.467,1520.587,
     &                     1539.379,1557.315,1576.506,1597.437,1618.199,1638.481,
     &                     1660.296,1681.662,1703.904,1728.751,1753.723,1778.133,
     &                     1803.904,1829.724,1855.804,1885.518,1915.69,1945.191,
     &                     1975.872,2007.366,2038.423,2074.791,2112.287,2149.08,
     &                     2186.769,2226.189,2263.994,2309.212,2356.342,2402.673,
     &                     2449.55,2499.607,2546.478,2606.28,2669.901,2731.404,
     &                     2789.465,2846.662/)

      tcon_10X_MET_KCl = (/616.369,617.639,621.006,626.112,632.265,638.794,
     &                     645.055,650.839,657.13,663.763,670.199,676.448,
     &                     683.327,689.879,696.858,704.433,711.907,719.168,
     &                     726.954,734.522,742.25,750.695,759.037,767.098,
     &                     775.633,784.206,792.828,802.601,812.472,822.069,
     &                     832.008,842.155,852.009,863.294,874.726,885.786,
     &                     897.066,908.88,920.075,933.197,946.66,959.72,
     &                     972.852,986.895,999.901,1016.077,1033.007,1049.174,
     &                     1064.368,1079.477/)

      tcon_10X_MET_Mg2SiO4 = (/1434.81,1437.596,1444.981,1456.179,1469.679,1484.01,
     &                     1497.768,1510.501,1524.416,1539.149,1553.484,1567.412,
     &                     1582.689,1597.198,1612.541,1629.076,1645.31,1661.051,
     &                     1677.996,1694.535,1711.552,1730.31,1748.958,1767.037,
     &                     1786.119,1805.204,1824.288,1845.759,1867.329,1888.228,
     &                     1909.899,1932.084,1953.668,1978.445,2003.601,2027.985,
     &                     2052.876,2078.954,2103.723,2132.821,2162.775,2191.965,
     &                     2221.476,2253.215,2283.001,2321.192,2361.95,2401.44,
     &                     2438.751,2475.443/)

      tcon_10X_MET_MnS = (/1138.595,1140.64,1146.063,1154.286,1164.198,1174.719,
     &                     1184.815,1194.151,1204.345,1215.125,1225.599,1235.757,
     &                     1246.886,1257.421,1268.475,1280.293,1291.82,1302.95,
     &                     1314.954,1326.671,1338.695,1351.911,1365.012,1377.678,
     &                     1391.028,1404.347,1417.567,1432.295,1446.977,1461.125,
     &                     1475.796,1490.838,1505.441,1522.152,1539.068,1555.414,
     &                     1572.063,1589.469,1605.908,1625.132,1644.755,1663.624,
     &                     1682.37,1702.133,1719.873,1740.272,1760.484,1778.961,
     &                     1796.04,1813.616/)

      tcon_10X_MET_Na2S = (/798.502,800.23,804.815,811.767,820.143,829.028,
     &                     837.543,845.401,853.927,862.896,871.582,880.008,
     &                     889.291,898.133,907.54,917.737,927.799,937.59,
     &                     948.121,958.412,969.049,980.836,992.604,1004.05,
     &                     1016.129,1028.214,1040.342,1054.049,1067.868,1081.291,
     &                     1095.21,1109.447,1123.316,1139.267,1155.482,1171.213,
     &                     1187.265,1204.063,1220.009,1238.767,1258.054,1276.778,
     &                     1295.576,1315.605,1334.066,1356.763,1380.339,1402.723,
     &                     1423.715,1444.682/)

      tcon_10X_MET_Ni = (/1302.841,1305.732,1313.392,1325.007,1339.015,1353.897,
     &                     1368.202,1381.471,1396.048,1411.552,1426.699,1441.452,
     &                     1457.612,1472.977,1489.282,1506.914,1524.28,1541.156,
     &                     1559.321,1577.074,1595.413,1615.723,1635.99,1655.696,
     &                     1676.5,1697.322,1718.22,1741.839,1765.666,1788.845,
     &                     1812.926,1837.637,1861.875,1890.039,1918.9,1947.077,
     &                     1975.883,2006.003,2034.733,2068.788,2104.034,2138.469,
     &                     2173.2,2210.284,2244.801,2288.229,2334.026,2378.007,
     &                     2419.428,2460.437/)

      tcon_10X_MET_SiO2 = (/1362.697,1365.095,1371.45,1381.087,1392.706,1405.044,
     &                     1416.892,1427.864,1439.88,1452.622,1465.032,1477.087,
     &                     1490.285,1502.792,1515.954,1530.067,1543.868,1557.212,
     &                     1571.594,1585.627,1600.033,1615.873,1631.585,1646.79,
     &                     1662.836,1678.874,1694.869,1712.801,1730.765,1748.132,
     &                     1766.132,1784.557,1802.437,1822.879,1843.574,1863.602,
     &                     1884.063,1905.553,1926.012,1950.104,1974.985,1999.343,
     &                     2024.098,2050.869,2076.304,2109.828,2146.208,2181.887,
     &                     2215.743,2248.735/)

      tcon_10X_MET_CaTiO3 = (/1660.892,1663.96,1672.111,1684.468,1699.349,1715.112,
     &                     1730.187,1744.046,1758.943,1774.481,1789.426,1803.865,
     &                     1819.83,1835.036,1851.189,1868.675,1885.895,1902.608,
     &                     1920.547,1938.0,1955.851,1975.397,1994.719,2013.374,
     &                     2033.068,2052.757,2072.495,2094.86,2117.191,2138.356,
     &                     2159.45,2179.707,2196.418,2210.073,2220.05,2227.415,
     &                     2235.688,2247.415,2261.083,2282.7,2309.02,2337.527,
     &                     2366.764,2395.904,2422.908,2456.537,2491.765,2525.427,
     &                     2557.072,2588.523/)

      tcon_10X_MET_VO = (/1325.748,1328.323,1335.129,1345.453,1357.915,1371.176,
     &                     1383.962,1395.882,1409.133,1423.375,1437.408,1451.145,
     &                     1466.141,1480.421,1495.697,1512.343,1528.789,1544.733,
     &                     1561.716,1578.058,1594.336,1611.606,1628.291,1644.221,
     &                     1661.275,1678.648,1696.432,1717.053,1738.221,1759.009,
     &                     1780.471,1802.209,1823.23,1847.151,1871.252,1894.448,
     &                     1918.034,1942.688,1965.852,1992.602,2019.72,2045.746,
     &                     2071.77,2099.61,2125.129,2156.081,2187.939,2217.973,
     &                     2246.065,2274.278/)

      tcon_10X_MET_ZnS = (/696.328,697.457,700.452,704.993,710.464,716.265,
     &                     721.822,726.945,732.494,738.321,743.953,749.405,
     &                     755.404,761.1,767.11,773.573,779.905,786.036,
     &                     792.637,799.075,805.685,812.956,820.165,827.131,
     &                     834.465,841.766,848.982,856.975,864.905,872.518,
     &                     880.408,888.497,896.321,905.224,914.196,922.835,
     &                     931.624,940.82,949.476,959.505,969.7,979.508,
     &                     989.328,999.831,1009.482,1021.254,1033.418,1044.921,
     &                     1055.691,1066.483/)














      tcon_100X_MET_Al2O3 = (/2015.394,2024.871,2037.947,2053.505,2071.331,2089.608,
     &                     2107.214,2124.489,2142.771,2160.852,2180.005,2199.131,
     &                     2218.926,2238.385,2258.392,2278.187,2298.434,2318.589,
     &                     2339.259,2360.153,2381.806,2403.908,2426.434,2449.378,
     &                     2472.779,2496.812,2521.156,2546.132,2571.246,2597.069,
     &                     2622.894,2649.555,2676.111,2703.627,2730.926,2759.169,
     &                     2787.502,2819.36,2852.059,2886.946,2920.112,2949.677,
     &                     2973.03,2989.923,2999.253,3002.777,3001.935,3000.888,
     &                     2999.0,2999.0/)

      tcon_100X_MET_CaSiO4 = (/1765.527,1775.325,1789.011,1805.418,1824.286,1843.703,
     &                     1862.532,1881.069,1900.682,1920.103,1940.492,1960.727,
     &                     1981.869,2002.916,2024.812,2046.779,2069.583,2092.623,
     &                     2116.547,2140.915,2166.085,2191.786,2218.065,2244.817,
     &                     2271.982,2299.723,2327.755,2356.564,2385.657,2415.77,
     &                     2446.094,2477.579,2509.03,2541.673,2574.042,2607.6,
     &                     2640.678,2675.211,2709.144,2744.652,2779.472,2817.386,
     &                     2854.556,2892.214,2926.147,2955.232,2976.937,2992.13,
     &                     2999.779,3002.665/)

      tcon_100X_MET_Cr = (/1446.911,1453.51,1462.723,1473.76,1486.447,1499.499,
     &                     1512.112,1524.45,1537.427,1550.286,1563.87,1577.368,
     &                     1591.452,1605.449,1619.918,1634.291,1649.08,1663.932,
     &                     1679.307,1694.901,1710.924,1727.187,1743.8,1760.729,
     &                     1777.948,1795.627,1813.602,1832.143,1850.898,1870.304,
     &                     1889.808,1910.036,1930.318,1951.56,1972.761,1994.714,
     &                     2016.337,2038.876,2061.087,2083.943,2106.898,2133.79,
     &                     2163.159,2199.304,2238.924,2288.776,2341.732,2392.952,
     &                     2436.738,2470.422/)

      tcon_100X_MET_Fe = (/1458.223,1467.782,1481.279,1497.581,1516.488,1536.105,
     &                     1555.328,1574.363,1594.421,1614.255,1635.149,1656.063,
     &                     1678.154,1700.411,1723.844,1747.635,1772.562,1797.936,
     &                     1824.318,1851.293,1879.32,1908.131,1937.908,1968.696,
     &                     2000.426,2033.323,2067.087,2102.348,2138.494,2176.469,
     &                     2215.281,2256.265,2298.066,2342.417,2387.453,2435.484,
     &                     2484.236,2536.609,2589.802,2648.457,2707.963,2773.25,
     &                     2836.466,2893.741,2939.935,2974.603,2995.426,3005.667,
     &                     3006.108,3003.701/)

      tcon_100X_MET_KCl = (/653.809,657.408,662.437,668.463,675.365,682.44,
     &                     689.27,695.985,703.108,710.173,717.655,725.17,
     &                     733.105,741.056,749.312,757.569,766.095,774.664,
     &                     783.506,792.483,801.781,811.298,821.075,831.096,
     &                     841.321,851.811,862.5,873.613,884.893,896.549,
     &                     908.316,920.644,933.098,946.182,959.335,973.247,
     &                     987.238,1002.039,1016.733,1032.252,1047.688,1064.293,
     &                     1080.839,1098.644,1116.348,1137.164,1158.202,1178.102,
     &                     1194.862,1208.196/)

      tcon_100X_MET_Mg2SiO4 = (/1605.876,1614.994,1627.624,1642.701,1660.081,1678.018,
     &                     1695.484,1712.708,1730.824,1748.72,1767.61,1786.398,
     &                     1806.055,1825.751,1846.35,1867.057,1888.577,1910.29,
     &                     1932.711,1955.482,1978.989,2003.016,2027.763,2053.181,
     &                     2079.148,2105.844,2133.025,2161.175,2189.836,2219.762,
     &                     2250.142,2282.001,2314.236,2348.179,2382.293,2418.091,
     &                     2453.835,2491.651,2529.391,2569.546,2609.572,2652.572,
     &                     2695.434,2741.668,2787.436,2840.203,2892.978,2942.629,
     &                     2984.097,3016.844/)

      tcon_100X_MET_MnS = (/1240.617,1247.816,1257.874,1269.926,1283.757,1297.946,
     &                     1311.649,1325.084,1339.154,1352.885,1367.156,1381.224,
     &                     1395.887,1410.414,1425.438,1440.391,1455.797,1471.174,
     &                     1486.925,1502.814,1519.107,1535.58,1552.39,1569.558,
     &                     1587.063,1605.03,1623.218,1641.882,1661.055,1681.7,
     &                     1702.223,1721.996,1739.294,1752.772,1762.243,1768.408,
     &                     1772.611,1778.705,1787.899,1803.634,1824.062,1850.59,
     &                     1880.056,1912.939,1946.894,1987.184,2029.234,2069.898,
     &                     2104.965,2133.253/)

      tcon_100X_MET_Na2S = (/875.815,880.855,888.01,896.676,906.733,917.166,
     &                     927.377,937.466,948.064,958.5,969.461,980.415,
     &                     991.953,1003.502,1015.55,1027.66,1040.281,1053.083,
     &                     1066.348,1079.863,1093.899,1108.298,1123.104,1138.309,
     &                     1153.892,1169.996,1186.476,1203.584,1221.002,1239.194,
     &                     1257.666,1277.015,1296.58,1317.173,1337.934,1359.909,
     &                     1381.958,1405.334,1428.658,1453.347,1477.7,1503.349,
     &                     1528.365,1554.536,1579.793,1608.228,1636.116,1662.114,
     &                     1683.574,1700.513/)

      tcon_100X_MET_Ni = (/1387.191,1395.43,1406.971,1420.831,1436.785,1453.227,
     &                     1469.211,1484.947,1501.54,1517.995,1535.403,1552.813,
     &                     1571.163,1589.604,1608.895,1628.334,1648.533,1668.937,
     &                     1690.053,1711.516,1733.678,1756.351,1779.777,1804.01,
     &                     1828.955,1854.753,1881.141,1908.578,1936.574,1965.845,
     &                     1995.585,2026.769,2058.418,2091.9,2125.722,2161.521,
     &                     2197.594,2236.085,2274.824,2316.374,2358.042,2402.834,
     &                     2447.669,2496.361,2545.095,2602.751,2661.45,2717.34,
     &                     2764.712,2802.555/)

      tcon_100X_MET_SiO2 = (/1512.78,1520.339,1530.938,1543.684,1558.451,1573.735,
     &                     1588.567,1603.127,1618.404,1633.418,1649.131,1664.724,
     &                     1681.071,1697.429,1714.591,1731.877,1749.779,1767.783,
     &                     1786.313,1805.039,1824.276,1843.821,1863.86,1884.44,
     &                     1905.505,1927.162,1949.22,1972.095,1995.303,2019.329,
     &                     2043.536,2068.762,2094.211,2120.993,2147.812,2175.812,
     &                     2203.721,2233.233,2262.589,2293.708,2324.684,2357.982,
     &                     2391.244,2427.116,2462.725,2504.65,2547.223,2587.751,
     &                     2622.051,2649.46/)

      tcon_100X_MET_CaTiO3 = (/1809.309,1818.827,1832.127,1848.077,1866.425,1885.315,
     &                     1903.672,1921.768,1940.856,1959.691,1979.422,1998.982,
     &                     2019.468,2039.92,2061.18,2082.454,2107.621,2136.242,
     &                     2157.766,2164.847,2161.784,2152.323,2140.785,2131.261,
     &                     2127.997,2135.501,2157.759,2191.532,2224.364,2254.327,
     &                     2284.777,2316.629,2348.755,2382.47,2416.392,2452.108,
     &                     2487.841,2525.655,2563.336,2603.351,2643.144,2685.618,
     &                     2727.914,2778.043,2828.503,2872.376,2906.081,2929.888,
     &                     2942.169,2945.743/)

      tcon_100X_MET_VO = (/1366.373,1373.412,1383.227,1394.98,1408.477,1422.359,
     &                     1435.885,1449.254,1463.351,1477.267,1491.891,1506.411,
     &                     1521.567,1536.646,1552.323,1568.071,1584.47,1601.016,
     &                     1618.071,1635.334,1653.128,1671.281,1689.85,1708.79,
     &                     1728.119,1748.003,1768.22,1789.108,1810.225,1832.006,
     &                     1854.07,1877.385,1901.188,1926.611,1952.679,1980.822,
     &                     2009.55,2040.342,2071.224,2103.975,2136.565,2171.216,
     &                     2205.399,2241.9,2277.78,2318.69,2359.19,2397.028,
     &                     2428.417,2453.156/)

      tcon_100X_MET_ZnS = (/765.565,769.004,773.79,779.517,786.139,792.984,
     &                     799.623,806.129,812.944,819.668,826.793,833.909,
     &                     841.364,848.825,856.582,864.295,872.179,880.015,
     &                     888.019,896.076,904.386,912.888,921.666,930.707,
     &                     939.967,949.51,959.22,969.233,979.403,990.036,
     &                     1000.801,1011.979,1023.233,1035.113,1047.082,1059.685,
     &                     1072.238,1085.387,1098.43,1112.266,1125.878,1140.139,
     &                     1154.111,1168.935,1183.419,1200.177,1216.925,1232.709,
     &                     1245.965,1256.601/)














      tcon_300X_MET_Al2O3 = (/2111.43,2122.387,2137.576,2155.697,2176.527,2197.93,
     &                     2218.559,2238.725,2259.819,2280.451,2301.942,2323.066,
     &                     2344.998,2366.775,2389.494,2412.282,2435.921,2459.726,
     &                     2484.319,2509.265,2534.918,2560.981,2587.596,2614.687,
     &                     2642.156,2670.149,2698.371,2727.303,2756.513,2786.775,
     &                     2817.634,2850.914,2884.031,2916.283,2944.83,2968.283,
     &                     2985.258,2996.113,3000.949,3001.922,3000.581,2999.924,
     &                     2999.0,2999.0,2999.0,2999.0,2999.0,2999.0,
     &                     2999.0,2999.0/)

      tcon_300X_MET_CaSiO4 = (/1863.375,1874.412,1889.784,1908.181,1929.338,1951.118,
     &                     1972.253,1993.049,2014.935,2036.566,2059.472,2082.396,
     &                     2106.537,2130.777,2156.024,2181.286,2207.435,2233.743,
     &                     2260.745,2287.959,2315.919,2344.346,2373.472,2403.243,
     &                     2433.513,2464.515,2495.979,2528.349,2560.926,2594.446,
     &                     2627.936,2662.396,2696.623,2731.991,2767.075,2804.064,
     &                     2840.587,2878.635,2914.197,2945.483,2969.98,2987.778,
     &                     2998.006,3002.708,3002.575,3001.232,2999.0,2999.0,
     &                     2999.0,2999.0/)

      tcon_300X_MET_Cr = (/1528.998,1536.837,1547.992,1561.431,1576.55,1591.759,
     &                     1606.079,1619.801,1634.026,1647.883,1662.49,1677.185,
     &                     1693.064,1709.152,1725.861,1742.603,1759.845,1777.079,
     &                     1794.779,1812.675,1831.091,1849.845,1869.087,1888.794,
     &                     1908.918,1929.58,1950.51,1972.049,1993.84,2016.397,
     &                     2039.139,2062.897,2086.713,2111.324,2135.75,2161.291,
     &                     2186.704,2213.555,2240.298,2268.756,2297.145,2327.371,
     &                     2357.235,2389.386,2421.572,2460.48,2500.649,2539.415,
     &                     2572.715,2599.663/)

      tcon_300X_MET_Fe = (/1501.934,1512.252,1526.695,1544.055,1564.149,1584.964,
     &                     1605.313,1625.488,1646.936,1668.285,1690.8,1713.322,
     &                     1737.166,1761.298,1786.622,1812.166,1838.81,1865.904,
     &                     1894.124,1923.008,1953.046,1984.061,2016.364,2049.914,
     &                     2084.454,2120.242,2157.031,2195.518,2235.026,2276.613,
     &                     2319.222,2364.309,2410.394,2459.449,2509.417,2562.898,
     &                     2617.52,2677.848,2738.846,2802.619,2862.167,2914.179,
     &                     2954.425,2983.061,2998.747,3005.322,3004.255,3002.38,
     &                     2999.0,2999.0/)

      tcon_300X_MET_KCl = (/670.057,673.975,679.321,685.649,692.966,700.522,
     &                     707.781,714.877,722.416,729.912,737.796,745.593,
     &                     753.738,761.951,770.648,779.435,788.516,797.67,
     &                     807.15,816.748,826.549,836.48,846.714,857.274,
     &                     868.131,879.38,890.871,902.748,914.834,927.481,
     &                     940.32,953.772,967.401,981.757,996.252,1011.635,
     &                     1027.094,1043.519,1060.015,1077.612,1095.166,1113.992,
     &                     1132.734,1152.931,1173.03,1196.574,1220.359,1242.901,
     &                     1261.867,1276.894/)

      tcon_300X_MET_Mg2SiO4 = (/1695.715,1705.754,1719.839,1736.781,1756.345,1776.56,
     &                     1796.238,1815.667,1836.295,1856.791,1878.318,1899.738,
     &                     1922.229,1944.747,1968.289,1992.003,2016.63,2041.527,
     &                     2067.401,2093.776,2121.009,2148.846,2177.5,2206.987,
     &                     2237.293,2268.658,2300.698,2333.929,2367.787,2403.143,
     &                     2439.014,2476.57,2514.452,2554.091,2593.866,2635.795,
     &                     2677.91,2723.067,2768.448,2816.333,2863.739,2914.397,
     &                     2964.532,3018.196,3071.155,3132.308,3193.669,3251.839,
     &                     3300.417,3338.444/)

      tcon_300X_MET_MnS = (/1280.272,1288.31,1299.498,1312.884,1328.255,1344.059,
     &                     1359.436,1374.647,1390.746,1406.63,1423.243,1439.691,
     &                     1456.848,1473.896,1491.534,1509.124,1527.357,1545.693,
     &                     1564.54,1583.545,1603.005,1622.688,1642.936,1663.833,
     &                     1685.119,1706.684,1726.717,1743.626,1756.471,1764.576,
     &                     1769.833,1775.599,1783.375,1795.917,1813.32,1837.145,
     &                     1864.285,1894.419,1924.776,1955.72,1986.811,2020.344,
     &                     2053.774,2089.641,2125.442,2168.384,2212.422,2254.582,
     &                     2290.579,2319.6/)

      tcon_300X_MET_Na2S = (/913.256,919.015,927.071,936.735,947.838,959.271,
     &                     970.448,981.497,993.071,1004.5,1016.662,1028.876,
     &                     1041.717,1054.572,1067.962,1081.438,1095.517,1109.781,
     &                     1124.527,1139.558,1155.207,1171.311,1187.909,1204.976,
     &                     1222.5,1240.651,1259.222,1278.492,1298.151,1318.743,
     &                     1339.737,1361.853,1384.302,1408.008,1431.973,1457.339,
     &                     1482.764,1509.654,1536.461,1564.884,1593.004,1622.732,
     &                     1651.846,1682.617,1712.47,1745.763,1778.2,1808.223,
     &                     1832.83,1852.103/)

      tcon_300X_MET_Ni = (/1424.654,1433.332,1445.475,1460.051,1476.842,1494.153,
     &                     1510.991,1527.611,1545.212,1562.674,1581.079,1599.453,
     &                     1618.793,1638.201,1658.46,1678.831,1700.009,1721.442,
     &                     1743.67,1766.314,1789.802,1813.971,1838.993,1864.847,
     &                     1891.461,1919.053,1947.315,1976.677,2006.643,2038.058,
     &                     2070.104,2103.818,2138.06,2174.29,2211.028,2250.153,
     &                     2289.702,2331.947,2374.591,2420.503,2466.681,2516.435,
     &                     2566.323,2620.56,2674.929,2739.195,2804.541,2866.634,
     &                     2919.187,2961.109/)

      tcon_300X_MET_SiO2 = (/1588.41,1596.928,1608.827,1623.105,1639.616,1656.697,
     &                     1673.369,1689.819,1707.108,1724.159,1742.046,1759.76,
     &                     1778.24,1796.678,1815.93,1835.289,1855.506,1876.04,
     &                     1897.335,1918.969,1941.22,1963.83,1986.942,2010.525,
     &                     2034.544,2059.24,2084.453,2110.625,2137.237,2164.922,
     &                     2192.938,2222.248,2251.77,2282.576,2313.444,2346.038,
     &                     2378.743,2413.468,2448.24,2485.409,2522.622,2562.627,
     &                     2602.435,2645.194,2687.677,2737.508,2787.859,2835.557,
     &                     2875.767,2907.809/)

      tcon_300X_MET_CaTiO3 = (/1878.938,1889.416,1904.016,1921.497,1941.642,1962.408,
     &                     1982.516,2002.243,2027.212,2054.502,2071.992,2075.326,
     &                     2069.752,2058.376,2046.372,2036.964,2035.005,2043.956,
     &                     2068.562,2101.874,2131.829,2159.515,2187.982,2217.227,
     &                     2247.232,2278.226,2309.908,2342.831,2376.372,2411.385,
     &                     2446.944,2484.199,2521.763,2561.019,2600.353,2641.788,
     &                     2684.426,2736.396,2787.625,2822.448,2841.846,2851.716,
     &                     2853.807,2859.299,2867.748,2884.352,2909.779,2946.68,
     &                     2978.84,2994.876/)

      tcon_300X_MET_VO = (/1380.622,1387.979,1398.139,1410.234,1424.108,1438.358,
     &                     1452.148,1465.702,1479.99,1494.113,1509.039,1523.9,
     &                     1539.412,1554.889,1571.061,1587.336,1604.221,1621.217,
     &                     1638.771,1656.613,1675.034,1693.826,1713.074,1732.744,
     &                     1752.852,1773.56,1794.617,1816.384,1838.488,1861.43,
     &                     1884.589,1908.738,1933.08,1958.726,1984.601,2011.955,
     &                     2039.402,2068.475,2097.445,2128.151,2158.605,2190.766,
     &                     2222.614,2257.45,2292.479,2334.221,2376.929,2417.813,
     &                     2452.652,2480.618/)

      tcon_300X_MET_ZnS = (/800.758,804.596,809.91,816.255,823.635,831.303,
     &                     838.721,845.939,853.467,860.917,868.818,876.653,
     &                     884.799,892.949,901.493,910.029,918.76,927.477,
     &                     936.458,945.529,954.838,964.317,974.061,984.046,
     &                     994.279,1004.887,1015.708,1026.855,1038.16,1049.945,
     &                     1061.934,1074.581,1087.405,1100.884,1114.552,1129.172,
     &                     1143.929,1159.619,1175.472,1192.914,1210.634,1230.251,
     &                     1249.318,1266.63,1280.677,1291.327,1297.8,1301.013,
     &                     1301.208,1300.46/)

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


      tconds(3,1:51,1)=tcon_100X_MET_KCl
      tconds(3,1:51,2)=tcon_100X_MET_ZnS
      tconds(3,1:51,3)=tcon_100X_MET_Na2S
      tconds(3,1:51,4)=tcon_100X_MET_MnS
      tconds(3,1:51,5)=tcon_100X_MET_Cr
      tconds(3,1:51,6)=tcon_100X_MET_SiO2
      tconds(3,1:51,7)=tcon_100X_MET_Mg2SiO4
      tconds(3,1:51,8)=tcon_100X_MET_VO
      tconds(3,1:51,9)=tcon_100X_MET_Ni
      tconds(3,1:51,10)=tcon_100X_MET_Fe
      tconds(3,1:51,11)=tcon_100X_MET_CaSiO4
      tconds(3,1:51,12)=tcon_100X_MET_CaTiO3
      tconds(3,1:51,13)=tcon_100X_MET_Al2O3

      tconds(4,1:51,1)=tcon_300X_MET_KCl
      tconds(4,1:51,2)=tcon_300X_MET_ZnS
      tconds(4,1:51,3)=tcon_300X_MET_Na2S
      tconds(4,1:51,4)=tcon_300X_MET_MnS
      tconds(4,1:51,5)=tcon_300X_MET_Cr
      tconds(4,1:51,6)=tcon_300X_MET_SiO2
      tconds(4,1:51,7)=tcon_300X_MET_Mg2SiO4
      tconds(4,1:51,8)=tcon_300X_MET_VO
      tconds(4,1:51,9)=tcon_300X_MET_Ni
      tconds(4,1:51,10)=tcon_300X_MET_Fe
      tconds(4,1:51,11)=tcon_300X_MET_CaSiO4
      tconds(4,1:51,12)=tcon_300X_MET_CaTiO3
      tconds(4,1:51,13)=tcon_300X_MET_Al2O3

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
