C**********************************************************
C             MODULES CLOUDSETTINGS AND CLOUDTABLES
C**********************************************************
C     The cloud state that used to be shared through COMMON/CLOUDY/ and
C     COMMON/CLOUD_PROPERTIES/, each re-declared by hand in every
C     routine that used it (the copies had drifted apart in type and
C     layout). Both modules are PRIVATE apart from the listed variables,
C     so a bare "use" never exports the params.i constants (NL, MG, ...)
C     into a routine that includes params.i itself.

      MODULE cloudsettings
C     Cloud settings: the fort.7 &INCLOUDY keys (read in INISET) and the
C     prescribed-cloud fields derived from them (MAKECLOUDS).
      include 'params.i'
      SAVE
      PRIVATE
      PUBLIC :: AEROSOLMODEL, AERTOTTAU, CLOUDBASE, CLOUDTOP, CLDFRCT,
     &          AERHFRAC, PI0AERSW, ASYMSW, EXTFACTLW, PI0AERLW, ASYMLW,
     &          DELTASCALE, SIG_AREA, PHI_LON, TAUAEROSOL, AEROPROF,
     &          MAXTAU, MAXTAULOC, TCON, AEROSOLCOMP, MTLX, METALLICITY,
     &          HAZES, PICKET_FENCE_CLOUDS, MOLEF, AERLAYERS, GRAYCLDV,
     &          C_TO_O, CLDRAMP, HAZETYPE

      CHARACTER(30) :: AEROSOLMODEL, AEROSOLCOMP
      REAL AERTOTTAU, CLOUDBASE, CLOUDTOP, CLDFRCT, AERHFRAC
      REAL PI0AERSW, ASYMSW, EXTFACTLW, PI0AERLW, ASYMLW
      REAL SIG_AREA, PHI_LON
      LOGICAL DELTASCALE, HAZES, PICKET_FENCE_CLOUDS, GRAYCLDV
      REAL TAUAEROSOL(NL+1,MG,2,JG), AEROPROF(NL+1), TCON(NL+1)
      REAL MAXTAU
      INTEGER MAXTAULOC, AERLAYERS
C     Cloud species mixing ratios (fort.7 MOLEF), one per species.
      REAL MOLEF(NCLOUDS)
C     Metallicity (log10 solar) and C/O, which also pick the corr-k
C     opacity tables.
      REAL MTLX, METALLICITY, C_TO_O
C
C     Cloud ramp: every cloud and haze optical depth, prescribed or dye,
C     is scaled by model day/CLDRAMP until day CLDRAMP. Counted from the
C     start of the run (KOUNT is restored on restart). 0 (default) = no
C     ramp.
      REAL CLDRAMP
C
C     Haze optical properties: CLOUD_DATA/haze_<HAZETYPE>_*.txt, e.g.
C     soot (default), sulfur, tholin, soot-2xpi0.
      CHARACTER(30) :: HAZETYPE
      END MODULE cloudsettings


      MODULE cloudtables
C     Cloud and haze optical-property tables and condensation curves,
C     read by get_cloud_scattering_properties (cloud_properties_set_up.f)
C     and used by OPPRMULTI and OPPRMULTI_CORRK.
      include 'params.i'
      SAVE
      PRIVATE
      PUBLIC :: TCONDS, KE_OPPR, PI0_OPPR, G0_OPPR, DENSITY, FMOLW,
     &          CORFACT, input_particle_size_array_in_meters,
     &          input_temperature_array,
     &          particle_size_vs_layer_array_in_meters,
     &          input_pressure_array_cgs,
     &          HAZE_RosselandMean_tau_per_bar, HAZE_RosselandMean_pi0,
     &          HAZE_RosselandMean_gg, HAZE_PlanckMean_tau_per_bar,
     &          HAZE_PlanckMean_pi0, HAZE_PlanckMean_gg,
     &          HAZE_wav_tau_per_bar, HAZE_wav_pi0, HAZE_wav_gg,
     &          haze_pressure_array_pascals, HAZE_WAV_GRID,
     &          CLOUD_WAV_GRID, exp_92_lnsig2_pi

C     Condensation curves: (metallicity slot, pressure, species).
      REAL TCONDS(7,80,NCLOUDS)
C     Mie optics: (band, particle size, temperature, species).
      REAL KE_OPPR(5,100,100,NCLOUDS)
      REAL PI0_OPPR(5,100,100,NCLOUDS)
      REAL G0_OPPR(5,100,100,NCLOUDS)
      REAL DENSITY(NCLOUDS), FMOLW(NCLOUDS)
      REAL CORFACT(80)
      REAL input_particle_size_array_in_meters(100)
      REAL input_temperature_array(100)
      REAL particle_size_vs_layer_array_in_meters(80)
      REAL input_pressure_array_cgs(80)
      REAL, dimension(50,100) :: HAZE_RosselandMean_tau_per_bar,
     &     HAZE_RosselandMean_pi0, HAZE_RosselandMean_gg
      REAL, dimension(50,100) :: HAZE_PlanckMean_tau_per_bar,
     &     HAZE_PlanckMean_pi0, HAZE_PlanckMean_gg
      REAL, dimension(500,100) :: HAZE_wav_tau_per_bar, HAZE_wav_pi0,
     &     HAZE_wav_gg
      REAL haze_pressure_array_pascals(100)
      REAL HAZE_WAV_GRID(500)
      REAL CLOUD_WAV_GRID(100)
C     Lognormal correction exp(-9/2 ln^2 sigma)/pi to the cloud optics.
      REAL exp_92_lnsig2_pi
      END MODULE cloudtables
