      SUBROUTINE OPPRMULTI(TAURAY,TAUL,TAUGAS,TAUAER,solar_calculation_indexer, DPG,
     &             LLA, LLS, JDBLE, JDBLEDBLE, JN, JN2, iblackbody_above, ISL, IR, IRS,EMISIR,
     &             EPSILON, HEATI, HEATS, HEAT, SOLNET,TPI, SQ3, SBK,AM, AVG, ALOS,
     &  SCDAY,RGAS,GANGLE,GWEIGHT,GRATIO,EMIS,RSFX,NPROB,SOL,RAYPERBAR,WEIGHT,
     &  GOL,WOL,WAVE,TT,Y3,U0,FDEGDAY,
     &  WOT,GOT,PTEMPG,PTEMPT,G0,OPD,PTEMP,
     &  uG0,uTAUL,W0,uW0,uopd,U1S,
     &  U1I,TOON_AK,B1,B2,EE1,EM1,
     &  EM2,EL1,EL2,GAMI,AF,
     &  BF,EF,SFCS,B3,CK1,CK2,
     &  CP,CPB,CM,CMB,DIRECT,EE3,
     &  EL3,FNET,TMI,AS,DF,
     &  DS,XK,DIREC,DIRECTU,DINTENT,
     &  UINTENT,TMID,TMIU,tslu,total_downwelling,alb_tot,
     &  tiru,firu,fird,fsLu,fsLd,fsLn,alb_toa,fupbs,
     &  fdownbs,fnetbs,fdownbs2,fupbi,fdownbi,fnetbi,
     &  qrad,alb_tomi,alb_toai, p_pass,
     &  PI0_TEMP, G0_TEMP, tauaer_temp,j1,denom,kount)
!
!     **************************************************************
!     *  Purpose             :  CaLculates optical properties      *
!     *                         such as single scattering albedo,  *
!     *                         asymmetry parameter, etc.          *
!     *                         This routine is case dependent and *
!     *                         wiLL have to be repLaced by the    *
!     *                         user.                              *
!     *  Subroutines Called  :  None                               *
!     *  Output              :  TAUL, W0, G0, OPD, Y3              *
!     * ************************************************************
!
      include 'rcommons.h'

      INTEGER LLA, LLS, JDBLE, JDBLEDBLE, JN, JN2, iblackbody_above, ISL, IR, IRS, j1,kount, MET_INDEX
      REAL EMISIR, EPSILON, HEATI(NLAYER), HEATS(NLAYER), HEAT(NLAYER), SOLNET
      REAL TPI, SQ3, SBK,AM, AVG, ALOS
      REAL SCDAY, RGAS, GANGLE(3), GWEIGHT(3), GRATIO(3), EMIS(5), RSFX(5),NPROB(5), SOL(5),RAYPERBAR(5),WEIGHT(5)
      REAL GOL(5,2*NL+2), WOL(5,2*NL+2), WAVE(5+1), TT(NL+1), Y3(5,3,2*NL+2), U0, FDEGDAY
      REAL WOT, GOT, PTEMPG(5), PTEMPT(5), G0(5,2*NL+2), OPD( 5,2*NL+2), PTEMP(5,2*NL+2)
      REAL uG0(5,2*NL+2), uTAUL(5,2*NL+2), W0(5,2*NL+2), uW0(5,2*NL+2), uopd(5,2*NL+2),  U1S( 5)
      REAL U1I(5), TOON_AK(5,2*NL+2), B1(5,2*NL+2), B2(  5,2*NL+2), EE1( 5,2*NL+2), EM1(5,2*NL+2)
      REAL EM2(5,2*NL+2), EL1( 5,2*NL+2), EL2(5,2*NL+2), GAMI(5,2*NL+2), AF(5,4*NL+4)
      REAL BF(5,4*NL+4), EF(5,4*NL+4), SFCS(5), B3(5,2*NL+2), CK1(5,2*NL+2), CK2(5,2*NL+2)
      REAL CP(5,2*NL+2), CPB(5,2*NL+2), CM(5,2*NL+2), CMB(5,2*NL+2), DIRECT(5,2*NL+2), EE3(5,2*NL+2)
      REAL EL3(5,2*NL+2), FNET(5,2*NL+2), TMI(5,2*NL+2), AS(5,4*NL+4), DF(5,4*NL+4)
      REAL DS(5,4*NL+4), XK(5,4*NL+4), DIREC(5,2*NL+2), DIRECTU(5,2*NL+2), DINTENT(5,3,2*NL+2)
      REAL UINTENT(5,3,2*NL+2), TMID(5,2*NL+2), TMIU(5,2*NL+2), tslu,total_downwelling,alb_tot
      REAL tiru,firu(2),fird(2),fsLu(3), fsLd(3),fsLn(3),alb_toa(3), fupbs(NL+1)
      REAL fdownbs(NL+1),fnetbs(NL+1),fdownbs2(NL+1), fupbi(NL+1),fdownbi(NL+1),fnetbi(NL+1)
      REAL qrad(NL+1),alb_tomi,alb_toais

      REAL DENOM
      REAL DPG(NLAYER), p_pass(NLAYER), layer_pressure_bar(NLAYER)
      REAL CONDFACT(NLAYER,NCLOUDS)

      REAL PI0_TEMP(NSOL + NIR, NVERT, NCLOUDS)
      REAL G0_TEMP(NSOL + NIR, NVERT, NCLOUDS)
      REAL tauaer_temp(NTOTAL, NLAYER, NCLOUDS)

      REAL CLOUDLOC(NL+1,NCLOUDS)
      INTEGER BASELEV
      INTEGER TOPLEV(NCLOUDS)

      real, dimension(NIR+NSOL,2*NL+2) :: TAURAY,TAUL,TAUGAS,TAUAER
      real, dimension(NIR+NSOL,NL+1) :: TAU_HAZE

      ! These are hardcoded to 50 but they are just lookup tables
      ! Don't worry about expanding the GCM to more levels
      real, dimension(50) :: input_temperature_array
      real, dimension(50) :: input_pressure_array_cgs

      real, dimension(50) :: input_particle_size_array_in_meters
      real, dimension(50) :: particle_size_vs_layer_array_in_meters

      REAL QE_OPPR(NSOL + NIR, 50, 50, NCLOUDS)
      REAL PI0_OPPR(NSOL + NIR, 50, 50, NCLOUDS)
      REAL G0_OPPR(NSOL + NIR, 50, 50, NCLOUDS)

      ! HAZE ARRAYS ARE DIFFERENT THAN THE OTHER ONES
      real, dimension(50, 100)  :: HAZE_RosselandMean_tau_per_bar, HAZE_RosselandMean_pi0, HAZE_RosselandMean_gg
      real, dimension(50, 100)  :: HAZE_PlanckMean_tau_per_bar, HAZE_PlanckMean_pi0, HAZE_PlanckMean_gg
      real, dimension(500, 100) :: HAZE_wav_tau_per_bar, HAZE_wav_pi0, HAZE_wav_gg
      real, dimension(100)      :: haze_pressure_array_pascals

      REAL TCONDS(5,51,NCLOUDS)
      REAL CORFACT(51)

      REAL DENSITY(NCLOUDS)
      REAL FMOLW(NCLOUDS)
      REAL MOLEF(NCLOUDS)

      INTEGER CLOUD_WAVELENGTH_INDEXES(NSOL + NIR)
      INTEGER HAZE_WAVELENGTH_INDEXES(NSOL + NIR)
      INTEGER WAV_LOC

      INTEGER K,J,L, iradgas
      INTEGER size_loc, temp_loc, solar_calculation_indexer, layer_index, haze_layer_index
      real particle_size

      REAL, dimension (500) :: HAZE_WAV_GRID
      REAL, dimension (50)  :: CLOUD_WAV_GRID

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

      ! This really should be moved to cloud_properties_set_up
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


      CLOUD_WAV_GRID = (/0.1, 0.111, 0.124, 0.138, 0.154, 0.172, 0.191, 0.213, 0.238, 0.265, 0.295, 0.329,
     &                   0.366, 0.408, 0.454, 0.506, 0.564, 0.629, 0.7, 0.78, 0.869, 0.969, 1.079, 1.202,
     &                   1.34, 1.493, 1.663, 1.853, 2.065, 2.301, 2.563, 2.856, 3.182, 3.546, 3.95, 4.401,
     &                   4.904, 5.464, 6.088, 6.783, 7.558, 8.421, 9.383, 10.454, 11.648, 12.978, 14.46,
     &                   16.111, 17.951, 20.0/)

      ! THE THREE Condensation curve sets are for 1X, 100X, and 300X Met
      ! Sorry that this is bad code
      ! Malsky
      IF (METALLICITY .gt. -0.1 .AND. METALLICITY .lt. 0.1) THEN
          MET_INDEX = 1
      ELSE IF (METALLICITY .gt. 0.9 .AND. METALLICITY .lt. 1.1) THEN
          MET_INDEX = 2
      ELSE IF (METALLICITY .gt. 1.9 .AND. METALLICITY .lt. 2.1) THEN
          MET_INDEX = 3
      ELSE IF (METALLICITY .gt. 2.37 .AND. METALLICITY .lt. 2.57) THEN
          MET_INDEX = 4
      ELSE
          write(*,*) 'Something is wrong with your metallicity'
          write(*,*) 'Check ropprrmulti'
          write(*,*) 'THE THREE Condensation curve sets are for 1X, 100X, and 300X Met'
          stop
      END IF


      DO J  = 2,NLAYER
          layer_pressure_bar(J)  = (p_pass(J)-p_pass(J-1)) * 1e-5
      END DO

      layer_pressure_bar(1)=(10.0**(LOG10(layer_pressure_bar(2))-(LOG10(layer_pressure_bar(3))-
     &                       LOG10(layer_pressure_bar(2)))))

      ! Find the index of the haze wavelengths in the visible
      HAZE_WAVELENGTH_INDEXES(1) = MINLOC(ABS((HAZE_WAV_GRID) - (RAYSCATLAM(1))),1)
      HAZE_WAVELENGTH_INDEXES(2) = MINLOC(ABS((HAZE_WAV_GRID) - (RAYSCATLAM(2))),1)
      HAZE_WAVELENGTH_INDEXES(3) = MINLOC(ABS((HAZE_WAV_GRID) - (RAYSCATLAM(3))),1)

      ! Find the index of the haze wavelengths at 5.0 microns
      HAZE_WAVELENGTH_INDEXES(4) = MINLOC(ABS((HAZE_WAV_GRID) - (5.0)),1)
      HAZE_WAVELENGTH_INDEXES(5) = MINLOC(ABS((HAZE_WAV_GRID) - (5.0)),1)

      ! Find the index of the cloud wavelengths in the visible
      CLOUD_WAVELENGTH_INDEXES(1) = MINLOC(ABS((CLOUD_WAV_GRID) - (RAYSCATLAM(1))),1)
      CLOUD_WAVELENGTH_INDEXES(2) = MINLOC(ABS((CLOUD_WAV_GRID) - (RAYSCATLAM(2))),1)
      CLOUD_WAVELENGTH_INDEXES(3) = MINLOC(ABS((CLOUD_WAV_GRID) - (RAYSCATLAM(3))),1)

      ! Find the index of the haze wavelengths at 5.0 microns
      CLOUD_WAVELENGTH_INDEXES(4) = MINLOC(ABS((CLOUD_WAV_GRID) - (5.0)),1)
      CLOUD_WAVELENGTH_INDEXES(5) = MINLOC(ABS((CLOUD_WAV_GRID) - (5.0)),1)

      Y3(:,:,:) = 0.0


      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!         GET THE HAZE DATA FIRST       !!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      IF (HAZES) THEN
          DO J = 1, NLAYER
              haze_layer_index = MINLOC(ABS((haze_pressure_array_pascals) - (p_pass(J))),1)  ! Both of these are in PA

              ! This grabs the optical depth per bar, then multiply it by the pressure in bars
              IF (PICKET_FENCE_CLOUDS .eq. .False.) THEN
                  DO L = solar_calculation_indexer,NSOLP
                      WAV_LOC = HAZE_WAVELENGTH_INDEXES(2) !THIS IS THE DOUBLE GRAY VERSION
                      TAU_HAZE(L,J) = HAZE_wav_tau_per_bar(WAV_LOC, haze_layer_index) * layer_pressure_bar(J)
                  END DO
              ELSE
                  DO L = solar_calculation_indexer,NSOLP
                      WAV_LOC = HAZE_WAVELENGTH_INDEXES(L)
                      TAU_HAZE(L,J) = HAZE_wav_tau_per_bar(WAV_LOC, haze_layer_index) * layer_pressure_bar(J)
                  END DO
              END IF
          END DO

          DO J = 1, NLAYER
              haze_layer_index = MINLOC(ABS((haze_pressure_array_pascals) - (p_pass(J))),1)  ! Both of these are in PA
              temp_loc         = MINLOC(ABS(input_temperature_array - (TT(J))),1) ! Not needed for the stellar calc

              ! This grabs the optical depth per bar, then multiply it by the pressure in bars
              IF (PICKET_FENCE_CLOUDS .eq. .False.) THEN
                  DO L = NSOLP+1,NTOTAL
                      WAV_LOC = HAZE_WAVELENGTH_INDEXES(4)
                      TAU_HAZE(L,J) = HAZE_wav_tau_per_bar(WAV_LOC, haze_layer_index) * layer_pressure_bar(J)
                  END DO
              ELSE
                  TAU_HAZE(NSOLP+1,J)=HAZE_PlanckMean_tau_per_bar(temp_loc, haze_layer_index)*layer_pressure_bar(J)
                  TAU_HAZE(NSOLP+2,J)=HAZE_RosselandMean_tau_per_bar(temp_loc, haze_layer_index)*layer_pressure_bar(J)
              END IF
          END DO
      ELSE
           TAU_HAZE = 0.0
      END IF


      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!         CLOUD SCATTERING PROPERTIES       !!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      DO J = 1,NLAYER - 1
          layer_index   = MINLOC(ABS(input_pressure_array_cgs - (p_pass(J) * 10.0)),1)
          temp_loc      = MINLOC(ABS(input_temperature_array - (TT(J))),1)
          particle_size = particle_size_vs_layer_array_in_meters(layer_index) ! Convert to CGS
          size_loc      = MINLOC(ABS(input_particle_size_array_in_meters - (particle_size)), 1)

          DO I = 1,NCLOUDS
              ! GET THE SCATTERING PROPERTIES
              IF (PICKET_FENCE_CLOUDS .eq. .False.) THEN
                  DO L = solar_calculation_indexer,NSOLP
                      WAV_LOC = CLOUD_WAVELENGTH_INDEXES(2)
                      PI0_TEMP(L,J,I) = PI0_OPPR(L,WAV_LOC,size_loc,I)
                      G0_TEMP(L,J,I)  = G0_OPPR(L,WAV_LOC,size_loc,I)
                  END DO

                  DO L = NSOLP+1,NTOTAL
                      WAV_LOC = CLOUD_WAVELENGTH_INDEXES(4)
                      PI0_TEMP(L,J,I) = PI0_OPPR(L,WAV_LOC,size_loc,I)
                      G0_TEMP(L,J,I)  = G0_OPPR(L,WAV_LOC,size_loc,I)
                  END DO
              ELSE
                  DO L = solar_calculation_indexer,NSOLP
                      WAV_LOC = CLOUD_WAVELENGTH_INDEXES(L)
                      PI0_TEMP(L,J,I) = PI0_OPPR(L,WAV_LOC,size_loc,I)
                      G0_TEMP(L,J,I)  = G0_OPPR(L,WAV_LOC,size_loc,I)
                  END DO

                  DO L = NSOLP+1,NTOTAL
                      PI0_TEMP(L,J,I) = PI0_OPPR(L,temp_loc,size_loc,I)
                      G0_TEMP(L,J,I)  = G0_OPPR(L,temp_loc,size_loc,I)
                  END DO
              END IF

              CONDFACT(J,I) = min(max((Tconds(MET_INDEX,layer_index,I)-TT(J))/10.,0.0),1.0)

              CLOUDLOC(J,I) = NINT(CONDFACT(J,I))*J
              BASELEV = MAXVAL(CLOUDLOC(1:50,I),1)
              TOPLEV(I)  = max(BASELEV-AERLAYERS,0)

              ! DPG is CGS before that 10x
              IF (PICKET_FENCE_CLOUDS .eq. .False.) THEN
                  DO L = solar_calculation_indexer,NSOLP
                      WAV_LOC = CLOUD_WAVELENGTH_INDEXES(2)

                      tauaer_temp(L,J,I) = (DPG(J)*10.0)*molef(I)*3./4./particle_size/density(I)*fmolw(I)*
     &                              CONDFACT(J,I)*MTLX*CORFACT(layer_index)*QE_OPPR(L,WAV_LOC,size_loc,I)
                  END DO
                  DO L = NSOLP+1,NTOTAL
                      WAV_LOC = CLOUD_WAVELENGTH_INDEXES(4)
                      tauaer_temp(L,J,I) = (DPG(J)*10.0)*molef(I)*3./4./particle_size/density(I)*fmolw(I)*
     &                              CONDFACT(J,I)*MTLX*CORFACT(layer_index)*QE_OPPR(L,WAV_LOC,size_loc,I)
                  END DO
              ELSE
                  DO L = solar_calculation_indexer,NSOLP
                      WAV_LOC = CLOUD_WAVELENGTH_INDEXES(L)
                      tauaer_temp(L,J,I) = (DPG(J)*10.0)*molef(I)*3./4./particle_size/density(I)*fmolw(I)*
     &                              CONDFACT(J,I)*MTLX*CORFACT(layer_index)*QE_OPPR(L,WAV_LOC,size_loc,I)
                  END DO

                  DO L = NSOLP+1,NTOTAL
                      tauaer_temp(L,J,I) = (DPG(J)*10.0)*molef(I)*3./4./particle_size/density(I)*fmolw(I)*
     &                                  CONDFACT(J,I)*MTLX*CORFACT(layer_index)*QE_OPPR(L,temp_loc,size_loc,I)
                  END DO
              END IF
          END DO
      END DO

      ! Uncomment for compact clouds I think
      DO I = 1,NCLOUDS
          DO J = 1, TOPLEV(I)
              tauaer_temp(:,J,I) = 0.0
          END DO

          !if (TOPLEV(I) + 4 .le. NLAYER) THEN
          !    tauaer_temp(:,TOPLEV(I)+4,I) = tauaer_temp(:,TOPLEV(I)+4,I)*0.01 !i.e.e(-3)
          !    tauaer_temp(:,TOPLEV(I)+3,I) = tauaer_temp(:,TOPLEV(I)+3,I)*0.03 !i.e.e(-3)
          !    tauaer_temp(:,TOPLEV(I)+2,I) = tauaer_temp(:,TOPLEV(I)+2,I)*0.1 !i.e.e(-2)
          !    tauaer_temp(:,TOPLEV(I)+1,I) = tauaer_temp(:,TOPLEV(I)+1,I)*0.3 !i.e.e(-1)
          !ENDIF
      END DO


      IF (PICKET_FENCE_CLOUDS .eq. .FALSE.) THEN
          !     SW AT STANDARD VERTICAL RESOLUTION
          DO J = 1,NLAYER
              haze_layer_index = MINLOC(ABS((haze_pressure_array_pascals) - (p_pass(J))),1) ! Pascals

              DO L = solar_calculation_indexer,NSOLP
                  WAV_LOC = CLOUD_WAVELENGTH_INDEXES(2)
                  TAUAER(L,J) = SUM(tauaer_temp(L,J,1:NCLOUDS)) + TAU_HAZE(L,J)
                  WOL(L,J)    = SUM(tauaer_temp(L,J,1:NCLOUDS)/(TAUAER(L,J)+1e-8) * PI0_TEMP(L,J,1:NCLOUDS))
     &                    + (TAU_HAZE(L,J) * HAZE_wav_pi0(WAV_LOC, haze_layer_index) / (TAUAER(L,J) + 1e-8))
                  GOL(L,J)    = SUM(tauaer_temp(L,J,1:NCLOUDS)/(TAUAER(L,J)+1e-8) * G0_TEMP(L,J,1:NCLOUDS))
     &                    + (TAU_HAZE(L,J) * HAZE_wav_gg(WAV_LOC, haze_layer_index)  / (TAUAER(L,J) + 1e-8))
              END DO
          END DO

    !     LW AT 2X VERTICAL RESOLUTION (FOR PERFORMANCE).
          k = 1
          DO J = 1,NDBL,2
              haze_layer_index = MINLOC(ABS((haze_pressure_array_pascals) - (p_pass(K))),1)  ! Both of these are in pa
              temp_loc         = MINLOC(ABS(input_temperature_array - (TT(K))),1) ! Not needed for the stellar calc

              JJ = J

              DO L = NSOLP+1,NTOTAL
                  ! GREP CHECK THIS
                  WAV_LOC = CLOUD_WAVELENGTH_INDEXES(4)
                  TAUAER(L,JJ) = SUM(tauaer_temp(L,K,1:NCLOUDS)) + TAU_HAZE(L,K)
                  WOL(L,JJ)    = SUM(tauaer_temp(L,K,1:NCLOUDS)/(TAUAER(L,JJ)+1e-8)*PI0_TEMP(L,K,1:NCLOUDS)) +
     &                              (TAU_HAZE(L,K) * HAZE_wav_pi0(WAV_LOC, haze_layer_index) / (TAUAER(L,JJ) + 1e-8))
                  GOL(L,JJ)    = SUM(tauaer_temp(L,K,1:NCLOUDS)/(TAUAER(L,JJ)+1e-8)*G0_TEMP(L,K,1:NCLOUDS))  +
     &                              (TAU_HAZE(L,K) * HAZE_wav_gg(WAV_LOC, haze_layer_index)  / (TAUAER(L,JJ) + 1e-8))
              END DO
              JJ = J+1
              DO L = NSOLP+1,NTOTAL
                  TAUAER(L,JJ) = TAUAER(L,JJ-1)
                  WOL(L,JJ)    = WOL(L,JJ-1)
                  GOL(L,JJ)    = GOL(L,JJ-1)
              END DO
              k = k+1
          END DO
      ELSE
    !     SW AT STANDARD VERTICAL RESOLUTION
          DO J = 1,NLAYER
              haze_layer_index = MINLOC(ABS((haze_pressure_array_pascals) - (p_pass(J))),1) ! Pascals
              DO L = solar_calculation_indexer,NSOLP
                  WAV_LOC = CLOUD_WAVELENGTH_INDEXES(L)
                  TAUAER(L,J) = SUM(tauaer_temp(L,J,1:NCLOUDS)) + TAU_HAZE(L,J)
                  WOL(L,J)    = SUM(tauaer_temp(L,J,1:NCLOUDS)/(TAUAER(L,J)+1e-8) * PI0_TEMP(L,J,1:NCLOUDS))
     &                    + (TAU_HAZE(L,J) * HAZE_wav_pi0(WAV_LOC, haze_layer_index) / (TAUAER(L,J) + 1e-8))
                  GOL(L,J)    = SUM(tauaer_temp(L,J,1:NCLOUDS)/(TAUAER(L,J)+1e-8) * G0_TEMP(L,J,1:NCLOUDS))
     &                    + (TAU_HAZE(L,J) * HAZE_wav_gg(WAV_LOC, haze_layer_index)  / (TAUAER(L,J) + 1e-8))
              END DO
          END DO

    !     LW AT 2X VERTICAL RESOLUTION (FOR PERFORMANCE).
          k = 1
          DO J = 1,NDBL,2
              haze_layer_index = MINLOC(ABS((haze_pressure_array_pascals) - (p_pass(K))),1) ! Both of these are in pa
              temp_loc         = MINLOC(ABS(input_temperature_array - (TT(K))),1) ! Not needed for the stellar calc

              JJ = J

              TAUAER(NSOLP+1,JJ) = SUM(tauaer_temp(NSOLP+1,K,1:NCLOUDS)) + TAU_HAZE(NSOLP+1,K)
              WOL(NSOLP+1,JJ) =
     &        SUM(tauaer_temp(NSOLP+1,K,1:NCLOUDS)/(TAUAER(NSOLP+1,JJ)+1e-8)*PI0_TEMP(NSOLP+1,K,1:NCLOUDS))
     &        + (TAU_HAZE(NSOLP+1,K) * HAZE_PlanckMean_pi0(temp_loc, haze_layer_index) / (TAUAER(NSOLP+1,JJ) + 1e-8))
              GOL(NSOLP+1,JJ) =
     &        SUM(tauaer_temp(NSOLP+1,K,1:NCLOUDS)/(TAUAER(NSOLP+1,JJ)+1e-8)*G0_TEMP(NSOLP+1,K,1:NCLOUDS))
     &        + (TAU_HAZE(NSOLP+1,K) * HAZE_PlanckMean_gg(temp_loc, haze_layer_index) / (TAUAER(NSOLP+1,JJ) + 1e-8))


              TAUAER(NSOLP+2,JJ) = SUM(tauaer_temp(NSOLP+2,K,1:NCLOUDS)) + TAU_HAZE(NSOLP+2,K)
              WOL(NSOLP+2,JJ) =
     &        SUM(tauaer_temp(NSOLP+2,K,1:NCLOUDS)/(TAUAER(NSOLP+2,JJ)+1e-8)*PI0_TEMP(NSOLP+2,K,1:NCLOUDS))
     &        + (TAU_HAZE(NSOLP+2,K) * HAZE_RosselandMean_pi0(temp_loc, haze_layer_index) / (TAUAER(NSOLP+2,JJ) + 1e-8))
              GOL(NSOLP+2,JJ) =
     &        SUM(tauaer_temp(NSOLP+2,K,1:NCLOUDS)/(TAUAER(NSOLP+2,JJ)+1e-8)*G0_TEMP(NSOLP+2,K,1:NCLOUDS))
     &        + (TAU_HAZE(NSOLP+2,K) * HAZE_RosselandMean_gg(temp_loc, haze_layer_index) / (TAUAER(NSOLP+2,JJ) + 1e-8))

              JJ = J+1
              DO L = NSOLP+1,NTOTAL
                  TAUAER(L,JJ) = TAUAER(L,JJ-1)
                  WOL(L,JJ)    = WOL(L,JJ-1)
                  GOL(L,JJ)    = GOL(L,JJ-1)
              END DO
              k = k+1
          END DO
      END IF
      
      ! Smooth out the cloud properties after doubling
      DO L = NSOLP+1,NTOTAL
          DO J = 2, NDBL-1, 2
              TAUAER(L,J) = (TAUAER(L,J+1) + TAUAER(L,J-1)) / 2.0
              WOL(L,J) = (WOL(L,J+1) + WOL(L,J-1)) / 2.0
              GOL(L,J) = (GOL(L,J+1) + GOL(L,J-1)) / 2.0
          END DO
      END DO


      iradgas = 1
      DO J = 1,NLAYER
          j1 = max(1, j-1)

!         First the solar at standard resolution
          DO L = solar_calculation_indexer,NSOLP
              TAUL(L,J) = TAUGAS(L,J)+TAURAY(L,J)+TAUAER(L,J)

              if(TAUL(L,J) .lt. 1d-6 ) then
                  TAUL(L,J) = 1d-6
              endif

              utauL(L,j)  = TAUL(L,J)
              WOT = (TAURAY(L,J)+TAUAER(L,J)*WOL(L,J))/TAUL(L,J)

              if (iradgas.eq.0) then
                  wot = woL(L,j)
              endif

              WOT       = min(1.0 - 1d-6,WOT)
              uw0(L,j)  = WOT
              DENOM     = (TAURAY(L,J) + TAUAER(L,J) * WOL(L,J))


              if( DENOM .LE. 1d-6 ) then
                  DENOM = 1d-6
              endif

              if( DENOM .GT. 1d-6 ) then
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


              IF (.False.) THEN
                  TAUL(L,J)   = TAUL(L,J) * (1.-WOT*(GOT*GOT))
                  W0(L,J)     = (1.-(GOT*GOT))*WOT/(1.-WOT*(GOT*GOT))
                  G0(L,J)     = GOT/(1.+GOT)
                  OPD(L,J)    = 0.0
                  OPD(L,J)    = OPD(L,J1)+TAUL(L,J)
              ELSE
                  W0(L,J)= uw0(L,J)
                  G0(L,J)= ug0(L,J)
                  TAUL(L,J)= utaul(L,J)
                  OPD(L,J)= uOPD(L,J)
              ENDIF


!             HERE'S WHERE YOU CAN HARDWIRE VALUES
              if( taul(L,j) .lt. 0.) then
                  write(*,*) 'ERROR! The VISIBLE layer optical depth is less than 0:', taul(L,j)
                  stop
              endif
          END DO
      END DO

!     NOW AGAIN FOR THE IR
      DO J = 1,NDBL
          j1 = max( 1, j-1 )
          DO L = NSOLP+1,NTOTAL
              TAUL(L,J) = TAUGAS(L,J)+TAURAY(L,J)+TAUAER(L,J)

              !if (iradgas.eq.0) then
              !    tauL(L,j) = tauaer(L,j)
              !endif

              !if( TAUL(L,J) .lt. 1d-6 ) then
              !    TAUL(L,J) = 1d-6
              !endif

              utauL(L,j)  = TAUL(L,J)
              WOT         = (TAURAY(L,J)+TAUAER(L,J)*WOL(L,J))/TAUL(L,J)
              if (iradgas.eq.0) then
                  wot = woL(L,j)
              endif

              WOT         = min(1.-1d-6,WOT)
              uw0(L,j)    = WOT
              DENOM       = (TAURAY(L,J)+ TAUAER(L,J)*WOL(L,J))

              if( DENOM .LE. 1d-6 ) then
                  DENOM = 1d-6
              endif

              if( DENOM .GT. 1d-6 ) then
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

              IF (.false.) THEN
                  TAUL(L,J)   = TAUL(L,J) * (1.-WOT*(GOT*GOT))
                  W0(L,J)     = (1.-(GOT*GOT))*WOT/(1.-WOT*(GOT*GOT))
                  G0(L,J)     = GOT/(1.+GOT)
                  OPD(L,J)    = 0.0
                  OPD(L,J)    = OPD(L,J1)+TAUL(L,J)
              ELSE
                  W0(L,J)= uw0(L,J)
                  G0(L,J)= ug0(L,J)
                  TAUL(L,J)= utaul(L,J)
                  OPD(L,J)= uOPD(L,J)
             END IF


             if(taul(L,j) .lt. 0.) then
                 write(*,*) 'ERROR! The IR layer optical depth is less than 0:', taul(L,j)
                 stop
             endif

          END DO

          DO I = 1,NGAUSS
              DO L = NSOLP+1,NTOTAL
                  Y3(L,I,J) =   EXP(-TAUL(L,J)/GANGLE(I))
              END DO
          END DO
      END DO

      RETURN
      END

