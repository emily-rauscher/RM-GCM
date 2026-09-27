      SUBROUTINE get_gas_opacity_corrk_wrapper(METALLICITY, C_TO_O, FBASEFLUX, GASCON, with_TiO_and_VO)
          use corrkmodule, only : corrk_setup, TS_CORRK, PS_CORRK, TS_LOG_CORRK, 
     &           PS_LOG_CORRK, WGTS_CORRK, WNO_EDGES, WNO_CTRS, STEL_SPEC, INT_SPEC, TAURAY_PER_DPG,
     &           OPAC_CORRK, PLANCK_INTS, PLANCK_TS, NWNO, MINWNOSTEL, CLOUD_KEXT, CLOUD_A, CLOUD_G
          ! include 'rcommons.h'
          ! include 'nwno.inc'
          REAL :: METALLICITY, C_TO_O, FBASEFLUX, GASCON
          REAL :: with_TiO_and_VO

          ! REAL :: OPAC_CORRK(73, 20, NWNO, 8)
          ! REAL :: TS_CORRK(73), PS_CORRK(20), TS_LOG_CORRK(73), PS_LOG_CORRK(20), WGTS_CORRK(8) 
          ! REAL :: WNO_EDGES(NWNO+1), WNO_CTRS(NWNO), STEL_SPEC(NWNO), INT_SPEC(NWNO)
          ! REAL :: OPAC_CIA(73, 20, NWNO), TAURAY_PER_DPG(73, 20, NWNO)
          ! REAL :: PLANCK_INTS(NWNO, 3925), PLANCK_TS(3925)
          
          WRITE(*,*) "In get_gas_opacity_corrk_wrapper"
          call corrk_setup(METALLICITY,C_TO_O, FBASEFLUX, GASCON, with_TiO_and_VO,TS_CORRK, PS_CORRK, TS_LOG_CORRK,
     &           PS_LOG_CORRK, WGTS_CORRK, WNO_EDGES, WNO_CTRS, STEL_SPEC, INT_SPEC, TAURAY_PER_DPG,
     &           OPAC_CORRK, PLANCK_INTS, PLANCK_TS, NWNO, MINWNOSTEL, CLOUD_KEXT, CLOUD_A, CLOUD_G)
          ! write(*,*) "After call to corrk_setup"
          ! write(*,*) nwno
          ! write(*,*) TS_CORRK(1), PS_CORRK(1), TS_LOG_CORRK(1), PS_LOG_CORRK(1)
          ! write(*,*) WGTS_CORRK(1), WNO_EDGES(1), WNO_CTRS(1)
          ! write(*,*) "stel_spec:", STEL_SPEC 
          ! write(*,*) INT_SPEC
          ! write(*,*) TAURAY_PER_DPG(1,1,1), OPAC_CORRK(1,1,1,1)
          ! ! write(*,*) PLANCK_INTS(1,1), PLANCK_TS, PLANCK_TS(1)
          ! write(*,*) "OPAC_CORRK(30,12,6,:):", OPAC_CORRK(30,12,6,:)
          
          
          ! stop
      END SUBROUTINE get_gas_opacity_corrk_wrapper