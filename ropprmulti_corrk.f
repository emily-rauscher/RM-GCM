      SUBROUTINE OPPRMULTI_CORRK(TAURAY,TAUL,TAUGAS,TAUAER,solar_calculation_indexer, DPG,
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
     &  PI0_TEMP, G0_TEMP, tauaer_temp,j1,denom,kount, itspd)
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
      use corrkmodule, only : MINWNOSTEL, CLOUD_A, CLOUD_G, CLOUD_KEXT, NWNO
      include 'rcommons.h'

      INTEGER LLA, LLS, JDBLE, JDBLEDBLE, JN, JN2, iblackbody_above, ISL, IR, IRS, j1,kount, itspd
      INTEGER, AUTOMATIC :: MET_INDEX
      REAL EMISIR, EPSILON, HEATI(NLAYER), HEATS(NLAYER), HEAT(NLAYER), SOLNET
      REAL TPI, SQ3, SBK,AM, AVG, ALOS
      REAL SCDAY, RGAS, GANGLE(3), GWEIGHT(3), GRATIO(3), EMIS(NTOTAL), RSFX(NTOTAL),NPROB(NTOTAL), SOL(NTOTAL)
      REAL RAYPERBAR(NTOTAL),WEIGHT(NTOTAL)
      REAL GOL(NTOTAL,2*NL+2), WOL(NTOTAL,2*NL+2), WAVE(5+1), TT(NL+1), Y3(NTOTAL,3,2*NL+2), U0, FDEGDAY
      REAL WOT, GOT, PTEMPG(NTOTAL), PTEMPT(NTOTAL), G0(NTOTAL,2*NL+2), OPD(NTOTAL,2*NL+2), PTEMP(NTOTAL,2*NL+2)
      REAL uG0(NTOTAL,2*NL+2), uTAUL(NTOTAL,2*NL+2), W0(NTOTAL,2*NL+2), uW0(NTOTAL,2*NL+2), uopd(NTOTAL,2*NL+2),  U1S(NTOTAL)
      REAL U1I(NTOTAL), TOON_AK(NTOTAL,2*NL+2), B1(NTOTAL,2*NL+2), B2(  5,2*NL+2), EE1(NTOTAL,2*NL+2), EM1(NTOTAL,2*NL+2)
      REAL EM2(NTOTAL,2*NL+2), EL1(NTOTAL,2*NL+2), EL2(NTOTAL,2*NL+2), GAMI(NTOTAL,2*NL+2), AF(NTOTAL,4*NL+4)
      REAL BF(NTOTAL,4*NL+4), EF(NTOTAL,4*NL+4), SFCS(NTOTAL), B3(NTOTAL,2*NL+2), CK1(NTOTAL,2*NL+2), CK2(NTOTAL,2*NL+2)
      REAL CP(NTOTAL,2*NL+2), CPB(NTOTAL,2*NL+2), CM(NTOTAL,2*NL+2), CMB(NTOTAL,2*NL+2), DIRECT(NTOTAL,2*NL+2), EE3(NTOTAL,2*NL+2)
      REAL EL3(NTOTAL,2*NL+2), FNET(NTOTAL,2*NL+2), TMI(NTOTAL,2*NL+2), AS(NTOTAL,4*NL+4), DF(NTOTAL,4*NL+4)
      REAL DS(NTOTAL,4*NL+4), XK(NTOTAL,4*NL+4), DIREC(NTOTAL,2*NL+2), DIRECTU(NTOTAL,2*NL+2), DINTENT(NTOTAL,3,2*NL+2)
      REAL UINTENT(NTOTAL,3,2*NL+2), TMID(NTOTAL,2*NL+2), TMIU(NTOTAL,2*NL+2), tslu,total_downwelling,alb_tot
      REAL tiru,firu(NIR),fird(NIR),fsLu(NSOL), fsLd(NSOL),fsLn(NSOL),alb_toa(NSOL), fupbs(NL+1)
      REAL fdownbs(NL+1),fnetbs(NL+1),fdownbs2(NL+1), fupbi(NL+1),fdownbi(NL+1),fnetbi(NL+1)
      REAL qrad(NL+1),alb_tomi,alb_toais

      REAL, AUTOMATIC :: FACTOR, RAMP
      REAL DENOM

      REAL DPG(NLAYER), p_pass(NLAYER)
      real, automatic :: layer_pressure_bar(NLAYER)
      real, automatic :: CONDFACT(NLAYER,NCLOUDS)

      REAL PI0_TEMP(NSOL + NIR, NVERT, NCLOUDS)
      REAL G0_TEMP(NSOL + NIR, NVERT, NCLOUDS)
      REAL tauaer_temp(NTOTAL, NLAYER, NCLOUDS)

      real, automatic :: CLOUDLOC(NL+1,NCLOUDS)
      integer, automatic :: BASELEV
      integer, automatic :: TOPLEV(NCLOUDS)

      real, dimension(NIR+NSOL,2*NL+2) :: TAURAY,TAUL,TAUGAS,TAUAER
      real, automatic, dimension(NIR+NSOL,NL+1) :: TAU_HAZE

      ! These are hardcoded to 50 but they are just lookup tables
      ! Don't worry about expanding the GCM to more levels
      real, dimension(100) :: input_temperature_array
      real, dimension(50) :: input_pressure_array_cgs

      real, dimension(100) :: input_particle_size_array_in_meters
      real, dimension(50) :: particle_size_vs_layer_array_in_meters

      REAL KE_OPPR(5, 100, 100, NCLOUDS) ! hardcoded for PF and DG!
      REAL PI0_OPPR(5, 100, 100, NCLOUDS) ! hardcoded for PF and DG!
      REAL G0_OPPR(5, 100, 100, NCLOUDS) ! hardcoded for PF and DG!

      ! HAZE ARRAYS ARE DIFFERENT THAN THE OTHER ONES
      real, dimension(50, 100)  :: HAZE_RosselandMean_tau_per_bar, HAZE_RosselandMean_pi0, HAZE_RosselandMean_gg
      real, dimension(50, 100)  :: HAZE_PlanckMean_tau_per_bar, HAZE_PlanckMean_pi0, HAZE_PlanckMean_gg
      real, dimension(500, 100) :: HAZE_wav_tau_per_bar, HAZE_wav_pi0, HAZE_wav_gg
      real, dimension(100)      :: haze_pressure_array_pascals

      REAL TCONDS(6,51,NCLOUDS)
      REAL CORFACT(51)

      REAL DENSITY(NCLOUDS)
      REAL FMOLW(NCLOUDS)
    !   REAL MOLEF(NCLOUDS)

      INTEGER, AUTOMATIC :: CLOUD_WAVELENGTH_INDEXES(5)
      INTEGER, AUTOMATIC :: HAZE_WAVELENGTH_INDEXES(5)
      INTEGER, AUTOMATIC :: WAV_LOC, chan_idx, stel_idx

      INTEGER, AUTOMATIC :: K, J, L, iradgas
      INTEGER, AUTOMATIC :: size_loc, temp_loc, layer_index, haze_layer_index
      INTEGER solar_calculation_indexer
      REAL, AUTOMATIC :: particle_size
      REAL, AUTOMATIC :: p_weight, size_weight, tconds_interp, corfact_interp
      REAL, AUTOMATIC :: tconds_lo, tconds_hi, kext_interp

      REAL, dimension (500) :: HAZE_WAV_GRID
      REAL, dimension (100)  :: CLOUD_WAV_GRID
      REAL :: exp_92_lnsig2_pi

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

      ! TAUAER/WOL/GOL are only filled for a subset of (L,J) index ranges below
      ! (e.g. solar L only up to J=NLAYER, IR L for the full doubled grid).
      ! Zero them fully here so unfilled elements are deterministic instead of
      ! whatever happened to be in the caller's buffer / uninitialized memory
      ! (which differed between OMP threads and inflated SUM(TAUAER) wildly).
      TAUAER = 0.0
      WOL    = 0.0
      GOL    = 0.0

      ! THE THREE Condensation curve sets are for 1X, 100X, and 300X Met
      ! Sorry that this is bad code
      ! Malsky
    !   write(*,*) 'made it to ropprmulti.f'
    !   write(*,*) 'input_pressure_array_cgs:', input_pressure_array_cgs
    !   write(*,*) NWNO
    !   flush(6)
    

      IF (aerosolcomp .eq. 'standard') THEN
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
            !   write(*,*) 'Something is wrong with your metallicity'
            !   write(*,*) 'Check ropprrmulti'
            !   write(*,*) 'THE THREE Condensation curve sets are for 1X, 100X, and 300X Met'
            !   stop
            ! If not in the list of real curves, use the interpolated curve
            MET_INDEX = 6
        END IF
        ! write(*,*) 'METALLICITY INDEX:', MET_INDEX
        ! flush(6)


        DO J  = 2,NLAYER
            layer_pressure_bar(J)  = (p_pass(J)-p_pass(J-1)) * 1e-5
        END DO

        layer_pressure_bar(1)=(10.0**(LOG10(layer_pressure_bar(2))-(LOG10(layer_pressure_bar(3))-
     &                       LOG10(layer_pressure_bar(2)))))

        !   WRITE(*,*) 'haze wl grid:', HAZE_WAV_GRID
        !   WRITE(*,*) 'haze_wl_indices:', HAZE_WAVELENGTH_INDEXES
        !   WRITE(*,*) 'RAYSCATLAM:', RAYSCATLAM
        ! Find the index of the haze wavelengths in the visible
        ! HAZE_WAVELENGTH_INDEXES(1) = MINLOC(ABS((HAZE_WAV_GRID) - (RAYSCATLAM(1))),1)
        ! HAZE_WAVELENGTH_INDEXES(2) = MINLOC(ABS((HAZE_WAV_GRID) - (RAYSCATLAM(2))),1)
        ! HAZE_WAVELENGTH_INDEXES(3) = MINLOC(ABS((HAZE_WAV_GRID) - (RAYSCATLAM(3))),1)

        ! ! Find the index of the haze wavelengths at 5.0 microns
        ! HAZE_WAVELENGTH_INDEXES(4) = MINLOC(ABS((HAZE_WAV_GRID) - (5.0)),1)
        ! HAZE_WAVELENGTH_INDEXES(5) = MINLOC(ABS((HAZE_WAV_GRID) - (5.0)),1)

        ! ! Find the index of the cloud wavelengths in the visible
        ! CLOUD_WAVELENGTH_INDEXES(1) = MINLOC(ABS((CLOUD_WAV_GRID) - (RAYSCATLAM(1))),1)
        ! CLOUD_WAVELENGTH_INDEXES(2) = MINLOC(ABS((CLOUD_WAV_GRID) - (RAYSCATLAM(2))),1)
        ! CLOUD_WAVELENGTH_INDEXES(3) = MINLOC(ABS((CLOUD_WAV_GRID) - (RAYSCATLAM(3))),1)

        ! ! Find the index of the haze wavelengths at 5.0 microns
        ! CLOUD_WAVELENGTH_INDEXES(4) = MINLOC(ABS((CLOUD_WAV_GRID) - (5.0)),1)
        ! CLOUD_WAVELENGTH_INDEXES(5) = MINLOC(ABS((CLOUD_WAV_GRID) - (5.0)),1)

        Y3(:,:,:) = 0.0

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!!!!         GET THE HAZE DATA FIRST       !!!!!!!!!!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! write(*,*) "Getting haze data"
        ! flush(6)
        IF (HAZES) THEN
            write(*,*) "ERROR: you are trying to run a correlated-k model with hazes"
            stop
            DO J = 1, NLAYER
                haze_layer_index = NEAREST_INDEX(haze_pressure_array_pascals, 100, p_pass(J))  ! Both of these are in PA
                ! This grabs the optical depth per bar, then multiply it by the pressure in bars
                IF (PICKET_FENCE_CLOUDS .eqv. .FALSE.) THEN
                    DO L = MAX(solar_calculation_indexer,MINWNOSTEL*8),NSOL
                        WAV_LOC = HAZE_WAVELENGTH_INDEXES(2) !THIS IS THE DOUBLE GRAY VERSION
                        TAU_HAZE(L,J) = HAZE_wav_tau_per_bar(WAV_LOC, haze_layer_index) * layer_pressure_bar(J)
                    END DO
                ELSE
                    DO L = MAX(solar_calculation_indexer,MINWNOSTEL*8),NSOL
                        WAV_LOC = HAZE_WAVELENGTH_INDEXES(L)
                        TAU_HAZE(L,J) = HAZE_wav_tau_per_bar(WAV_LOC, haze_layer_index) * layer_pressure_bar(J)
                    END DO
                END IF
            END DO

            DO J = 1, NLAYER
                haze_layer_index = NEAREST_INDEX(haze_pressure_array_pascals, 100, p_pass(J))  ! Both of these are in PA
                temp_loc         = NEAREST_INDEX(input_temperature_array, 100, TT(J)) ! Not needed for the stellar calc

                ! This grabs the optical depth per bar, then multiply it by the pressure in bars
                IF (PICKET_FENCE_CLOUDS .eqv. .FALSE.) THEN
                    DO L = NSOL+1,NTOTAL
                        WAV_LOC = HAZE_WAVELENGTH_INDEXES(4)
                        TAU_HAZE(L,J) = HAZE_wav_tau_per_bar(WAV_LOC, haze_layer_index) * layer_pressure_bar(J)
                    END DO
                ELSE
                    TAU_HAZE(NSOL+1,J)=HAZE_PlanckMean_tau_per_bar(temp_loc, haze_layer_index)*layer_pressure_bar(J)
                    TAU_HAZE(NSOL+2,J)=HAZE_RosselandMean_tau_per_bar(temp_loc, haze_layer_index)*layer_pressure_bar(J)
                END IF
            END DO
        ELSE
            TAU_HAZE = 0.0
        END IF
        ! write(*,*) "HAZE check over"
        ! flush(6)
        ! write(*,*) "kext checks:", CLOUD_KEXT(1,4,4), CLOUD_KEXT(1,4,5), CLOUD_KEXT(1,4,6)


        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!!!!         CLOUD SCATTERING PROPERTIES       !!!!!!!!!!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! Add these checks before the main loops
        ! if (MINWNOSTEL < 0) then
        !     write(*,*) "Error: MINWNOSTEL not properly initialized:", MINWNOSTEL
        !     stop
        ! endif
        ! if (NWNO <= 0) then
        !     write(*,*) "Error: NWNO not properly initialized:", NWNO
        !     stop
        ! endif
        ! write(*,*) NLAYER
        ! write(*,*) 'Size of input_pressure_array_cgs = ', size(input_pressure_array_cgs)
        ! write(*,*) 'Size of p_pass = ', size(p_pass)
        ! flush(6)
        ! write(*,*) "MTLX:", MTLX
        DO J = 1,NLAYER - 1
            ! write(*,*) "I am here!"
            ! flush(6)
            ! write(*,*) particle_size_vs_layer_array_in_meters
            ! write(*,*) input_pressure_array_cgs

            ! flush(6)

            layer_index   = NEAREST_INDEX(input_pressure_array_cgs, 50, p_pass(J) * 10.0)
            ! write(*,*) input_pressure_array_cgs, p_pass(J) * 10.0, layer_index
            ! if (layer_index < 1 .or. layer_index > size(input_pressure_array_cgs)) then
            !     write(*,*) "Error: invalid layer_index:", layer_index
            !     flush(6)
            !     stop
            ! endif
            ! write(*,*) "layer_index:", layer_index
            ! write(*,*) "particle_size_vs_layer_array_in_meters:", particle_size_vs_layer_array_in_meters
            ! write(*,*) "input_particle_size_array_in_meters:", input_particle_size_array_in_meters
            ! write(*,*) "input_temperature_array:", input_temperature_array
            ! flush(6)

            CALL LOG_INTERP_WEIGHTS(input_pressure_array_cgs, 50, p_pass(J) * 10.0,
     &                              layer_index, p_weight)
            particle_size = LERP(particle_size_vs_layer_array_in_meters(layer_index),
     &                            particle_size_vs_layer_array_in_meters(layer_index+1), p_weight)
            CALL LOG_INTERP_WEIGHTS(input_particle_size_array_in_meters, 100, particle_size,
     &                              size_loc, size_weight)
            !corr-k version:
            DO I = 1, NCLOUDS
                DO L = MAX(solar_calculation_indexer,MINWNOSTEL*8+1),NSOL
                    chan_idx = MODULO(L-1, 8) + 1
                    stel_idx = MODULO((L-chan_idx)/8,NWNO) + 1

                    PI0_TEMP(L,J,I) = LERP(CLOUD_A(I,size_loc,stel_idx), CLOUD_A(I,size_loc+1,stel_idx), size_weight)
                    G0_TEMP(L,J,I)  = LERP(CLOUD_G(I,size_loc,stel_idx), CLOUD_G(I,size_loc+1,stel_idx), size_weight)
                END DO
                ! since spectral, should be able to copy IR to vis at the end, but for consistency we'll do it here
                DO L = NSOL+1,NTOTAL
                    chan_idx = MODULO(L-1, 8) + 1
                    stel_idx = MODULO((L-chan_idx)/8,NWNO) + 1
                    PI0_TEMP(L,J,I) = LERP(CLOUD_A(I,size_loc,stel_idx), CLOUD_A(I,size_loc+1,stel_idx), size_weight)
                    G0_TEMP(L,J,I)  = LERP(CLOUD_G(I,size_loc,stel_idx), CLOUD_G(I,size_loc+1,stel_idx), size_weight)
                END DO
                tconds_lo = TCONDS(MET_INDEX,layer_index,I)
                tconds_hi = TCONDS(MET_INDEX,layer_index+1,I)
                tconds_interp = tconds_lo + p_weight * (tconds_hi - tconds_lo)
                kext_interp = LERP(CLOUD_KEXT(I,size_loc,stel_idx), CLOUD_KEXT(I,size_loc+1,stel_idx), size_weight)
                CONDFACT(J,I) = min(max((tconds_interp-TT(J))/10.,0.0),1.0)


                CLOUDLOC(J,I) = NINT(CONDFACT(J,I))*J
                BASELEV = MAXVAL(CLOUDLOC(1:50,I),1)
                TOPLEV(I)  = max(BASELEV-AERLAYERS,0)

                corfact_interp = LERP(CORFACT(layer_index), CORFACT(layer_index+1), p_weight)

                DO L = MAX(solar_calculation_indexer,MINWNOSTEL*8+1),NSOL
                    chan_idx = MODULO(L-1, 8) + 1
                    stel_idx = MODULO((L-chan_idx)/8,NWNO) + 1
                    tauaer_temp(L,J,I) = (DPG(J)*10.0)*molef(I)*3./4./particle_size/particle_size/particle_size/density(I)*
     &                              fmolw(I)*CONDFACT(J,I)*MTLX*corfact_interp*
     &                              LERP(CLOUD_KEXT(I,size_loc,stel_idx), CLOUD_KEXT(I,size_loc+1,stel_idx), size_weight) / 1.0e4 ! convert k from cm^2 to m^2
     &                              * exp_92_lnsig2_pi ! correction factor for mean vs median radius, divided by pi
                END DO
                DO L = NSOL+1,NTOTAL
                    chan_idx = MODULO(L-1, 8) + 1
                    stel_idx = MODULO((L-chan_idx)/8,NWNO) + 1
                    tauaer_temp(L,J,I) = (DPG(J)*10.0)*molef(I)*3./4./particle_size/particle_size/particle_size/density(I)*
     &                              fmolw(I)*CONDFACT(J,I)*MTLX*corfact_interp*
     &                              LERP(CLOUD_KEXT(I,size_loc,stel_idx), CLOUD_KEXT(I,size_loc+1,stel_idx), size_weight) / 1.0e4 ! convert k from cm^2 to m^2
     &                              * exp_92_lnsig2_pi ! correction factor for mean vs median radius, divided by pi
                END DO
            END DO
        END DO
        ! write(*,*) "CONDFACT:", CONDFACT(:,6)
        ! write(*,*) "CONDCURVE:", Tconds(MET_INDEX,:,6)
        ! stop

            ! DO I = 1,NCLOUDS
            !     ! GET THE SCATTERING PROPERTIES
            !     IF (PICKET_FENCE_CLOUDS .eqv. .False.) THEN
            !         DO L = MAX(solar_calculation_indexer,MINWNOSTEL*8),NSOL
            !             WAV_LOC = CLOUD_WAVELENGTH_INDEXES(2)
            !             PI0_TEMP(L,J,I) = PI0_OPPR(1,WAV_LOC,size_loc,I)
            !             G0_TEMP(L,J,I)  = G0_OPPR(1,WAV_LOC,size_loc,I)
            !         END DO

            !         DO L = NSOL+1,NTOTAL
            !             WAV_LOC = CLOUD_WAVELENGTH_INDEXES(4)
            !             PI0_TEMP(L,J,I) = PI0_OPPR(1,WAV_LOC,size_loc,I)
            !             G0_TEMP(L,J,I)  = G0_OPPR(1,WAV_LOC,size_loc,I)
            !         END DO
            !     ELSE
            !         DO L = MAX(solar_calculation_indexer,MINWNOSTEL*8),NSOL
            !             WAV_LOC = CLOUD_WAVELENGTH_INDEXES(L)
            !             PI0_TEMP(L,J,I) = PI0_OPPR(L,WAV_LOC,size_loc,I)
            !             G0_TEMP(L,J,I)  = G0_OPPR(L,WAV_LOC,size_loc,I)
            !         END DO

            !         DO L = NSOL+1,NTOTAL
            !             PI0_TEMP(L,J,I) = PI0_OPPR(L,temp_loc,size_loc,I)
            !             G0_TEMP(L,J,I)  = G0_OPPR(L,temp_loc,size_loc,I)
            !         END DO
            !     END IF

                ! CONDFACT(J,I) = min(max((Tconds(MET_INDEX,layer_index,I)-TT(J))/10.,0.0),1.0)

                ! CLOUDLOC(J,I) = NINT(CONDFACT(J,I))*J
                ! BASELEV = MAXVAL(CLOUDLOC(1:50,I),1)
                ! TOPLEV(I)  = max(BASELEV-AERLAYERS,0)

                ! DPG is CGS before that 10x
    !             IF (PICKET_FENCE_CLOUDS .eqv. .FALSE.) THEN
    !                 DO L = MAX(solar_calculation_indexer,MINWNOSTEL*8),NSOL
    !                     WAV_LOC = CLOUD_WAVELENGTH_INDEXES(2)

    !                     tauaer_temp(L,J,I) = (DPG(J)*10.0)*molef(I)*3./4./particle_size/particle_size/particle_size/density(I)*
    !  &                              fmolw(I)*CONDFACT(J,I)*MTLX*CORFACT(layer_index)*KE_OPPR(1,WAV_LOC,size_loc,I) / 1.0e4 ! convert k from cm^2 to m^2
    !  &                              * exp_92_lnsig2_pi ! correction factor for mean vs median radius, divided by pi
    !                 END DO
    !                 DO L = NSOL+1,NTOTAL
    !                     WAV_LOC = CLOUD_WAVELENGTH_INDEXES(4)
    !                     tauaer_temp(L,J,I) = (DPG(J)*10.0)*molef(I)*3./4./particle_size/particle_size/particle_size/density(I)*
    !  &                              fmolw(I)*CONDFACT(J,I)*MTLX*CORFACT(layer_index)*KE_OPPR(1,WAV_LOC,size_loc,I) / 1.0e4 ! convert k from cm^2 to m^2
    !  &                              * exp_92_lnsig2_pi ! correction factor for mean vs median radius, divided by pi
    !                 END DO
    !             ELSE
    !                 DO L = MAX(solar_calculation_indexer,MINWNOSTEL*8),NSOL
    !                     WAV_LOC = CLOUD_WAVELENGTH_INDEXES(L)
    !                     tauaer_temp(L,J,I) = (DPG(J)*10.0)*molef(I)*3./4./particle_size/particle_size/particle_size/density(I)*
    !  &                              fmolw(I)*CONDFACT(J,I)*MTLX*CORFACT(layer_index)*KE_OPPR(L,WAV_LOC,size_loc,I) / 1.0e4 ! convert k from cm^2 to m^2
    !  &                              * exp_92_lnsig2_pi ! correction factor for mean vs median radius, divided by pi
    !                 END DO

    !                 DO L = NSOL+1,NTOTAL
    !                     tauaer_temp(L,J,I) = (DPG(J)*10.0)*molef(I)*3./4./particle_size/particle_size/particle_size/density(I)*
    !  &                              fmolw(I)*CONDFACT(J,I)*MTLX*CORFACT(layer_index)*KE_OPPR(L,WAV_LOC,size_loc,I) / 1.0e4 ! convert k from cm^2 to m^2
    !  &                              * exp_92_lnsig2_pi ! correction factor for mean vs median radius, divided by pi
    !                 END DO
    !             END IF
        !     END DO
        ! END DO


        ! Uncomment for compact clouds I think (all the same for corrk)
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

        ! tauaer_temp is dimensioned to NLAYER but only computed up to NLAYER-1 (=NVERT) above;
        ! extend the last computed layer downward so the read at J=NLAYER below isn't uninitialized.
        tauaer_temp(:,NLAYER,:) = tauaer_temp(:,NLAYER-1,:)

        ! Corrk version:
        DO J=1, NLAYER
            haze_layer_index = NEAREST_INDEX(haze_pressure_array_pascals, 100, p_pass(J)) ! Pascals

            DO L = MAX(solar_calculation_indexer,MINWNOSTEL*8+1),NSOL
                WAV_LOC = CLOUD_WAVELENGTH_INDEXES(2) ! hazes should be not used, so doesn't matter that it's double gray
                TAUAER(L,J) = SUM(tauaer_temp(L,J,1:NCLOUDS)) + TAU_HAZE(L,J)
                WOL(L,J)    = SUM(tauaer_temp(L,J,1:NCLOUDS)/(TAUAER(L,J)+1e-8) * PI0_TEMP(L,min(J,NVERT),1:NCLOUDS))
     &                    + (TAU_HAZE(L,J) * HAZE_wav_pi0(WAV_LOC, haze_layer_index) / (TAUAER(L,J) + 1e-8))
                GOL(L,J)    = SUM(tauaer_temp(L,J,1:NCLOUDS)/(TAUAER(L,J)+1e-8) * G0_TEMP(L,min(J,NVERT),1:NCLOUDS))
     &                    + (TAU_HAZE(L,J) * HAZE_wav_gg(WAV_LOC, haze_layer_index)  / (TAUAER(L,J) + 1e-8))
            END DO
        END DO


        ! IF (PICKET_FENCE_CLOUDS .eqv. .FALSE.) THEN
            !     SW AT STANDARD VERTICAL RESOLUTION
        DO J = 1,NLAYER
            haze_layer_index = NEAREST_INDEX(haze_pressure_array_pascals, 100, p_pass(J)) ! Pascals

            DO L = MAX(solar_calculation_indexer,MINWNOSTEL*8+1),NSOL
                WAV_LOC = CLOUD_WAVELENGTH_INDEXES(2)
                TAUAER(L,J) = SUM(tauaer_temp(L,J,1:NCLOUDS)) + TAU_HAZE(L,J)
                WOL(L,J)    = SUM(tauaer_temp(L,J,1:NCLOUDS)/(TAUAER(L,J)+1e-8) * PI0_TEMP(L,min(J,NVERT),1:NCLOUDS))
     &                    + (TAU_HAZE(L,J) * HAZE_wav_pi0(WAV_LOC, haze_layer_index) / (TAUAER(L,J) + 1e-8))
                GOL(L,J)    = SUM(tauaer_temp(L,J,1:NCLOUDS)/(TAUAER(L,J)+1e-8) * G0_TEMP(L,min(J,NVERT),1:NCLOUDS))
     &                    + (TAU_HAZE(L,J) * HAZE_wav_gg(WAV_LOC, haze_layer_index)  / (TAUAER(L,J) + 1e-8))
            END DO
        END DO
        ! write(*,*) "TAUAER:", TAUAER
    !     LW AT 2X VERTICAL RESOLUTION (FOR PERFORMANCE).
        k = 1
        DO J = 1,NDBL,2
            haze_layer_index = NEAREST_INDEX(haze_pressure_array_pascals, 100, p_pass(K))  ! Both of these are in pa
            temp_loc         = NEAREST_INDEX(input_temperature_array, 100, TT(K)) ! Not needed for the stellar calc

            JJ = J

            DO L = NSOL+1,NTOTAL
                ! GREP CHECK THIS
                WAV_LOC = CLOUD_WAVELENGTH_INDEXES(4)
                TAUAER(L,JJ) = SUM(tauaer_temp(L,K,1:NCLOUDS)) + TAU_HAZE(L,K)
                WOL(L,JJ)    = SUM(tauaer_temp(L,K,1:NCLOUDS)/(TAUAER(L,JJ)+1e-8)*PI0_TEMP(L,min(K,NVERT),1:NCLOUDS)) !+
    !  &                              (TAU_HAZE(L,K) * HAZE_wav_pi0(WAV_LOC, haze_layer_index) / (TAUAER(L,JJ) + 1e-8))
                GOL(L,JJ)    = SUM(tauaer_temp(L,K,1:NCLOUDS)/(TAUAER(L,JJ)+1e-8)*G0_TEMP(L,min(K,NVERT),1:NCLOUDS))  !+
    !  &                              (TAU_HAZE(L,K) * HAZE_wav_gg(WAV_LOC, haze_layer_index)  / (TAUAER(L,JJ) + 1e-8))
            END DO
            JJ = J+1
            DO L = NSOL+1,NTOTAL
                TAUAER(L,JJ) = TAUAER(L,JJ-1)
                WOL(L,JJ)    = WOL(L,JJ-1)
                GOL(L,JJ)    = GOL(L,JJ-1)
            END DO
            k = k+1
        END DO
        ! ELSE
    !     !     SW AT STANDARD VERTICAL RESOLUTION
    !         DO J = 1,NLAYER
    !             haze_layer_index = MINLOC(ABS((haze_pressure_array_pascals) - (p_pass(J))),1) ! Pascals
    !             DO L = MAX(solar_calculation_indexer,MINWNOSTEL*8),NSOL
    !                 WAV_LOC = CLOUD_WAVELENGTH_INDEXES(L)
    !                 TAUAER(L,J) = SUM(tauaer_temp(L,J,1:NCLOUDS)) + TAU_HAZE(L,J)
    !                 WOL(L,J)    = SUM(tauaer_temp(L,J,1:NCLOUDS)/(TAUAER(L,J)+1e-8) * PI0_TEMP(L,J,1:NCLOUDS))
    !  &                    + (TAU_HAZE(L,J) * HAZE_wav_pi0(WAV_LOC, haze_layer_index) / (TAUAER(L,J) + 1e-8))
    !                 GOL(L,J)    = SUM(tauaer_temp(L,J,1:NCLOUDS)/(TAUAER(L,J)+1e-8) * G0_TEMP(L,J,1:NCLOUDS))
    !  &                    + (TAU_HAZE(L,J) * HAZE_wav_gg(WAV_LOC, haze_layer_index)  / (TAUAER(L,J) + 1e-8))
    !             END DO
    !         END DO

    !     !     LW AT 2X VERTICAL RESOLUTION (FOR PERFORMANCE).
    !         k = 1
    !         DO J = 1,NDBL,2
    !             haze_layer_index = MINLOC(ABS((haze_pressure_array_pascals) - (p_pass(K))),1) ! Both of these are in pa
    !             temp_loc         = MINLOC(ABS(input_temperature_array - (TT(K))),1) ! Not needed for the stellar calc

    !             JJ = J

    !             TAUAER(NSOL+1,JJ) = SUM(tauaer_temp(NSOL+1,K,1:NCLOUDS)) + TAU_HAZE(NSOL+1,K)
    !             WOL(NSOL+1,JJ) =
    !  &        SUM(tauaer_temp(NSOL+1,K,1:NCLOUDS)/(TAUAER(NSOL+1,JJ)+1e-8)*PI0_TEMP(NSOL+1,K,1:NCLOUDS))
    !  &        + (TAU_HAZE(NSOL+1,K) * HAZE_PlanckMean_pi0(temp_loc, haze_layer_index) / (TAUAER(NSOL+1,JJ) + 1e-8))
    !             GOL(NSOL+1,JJ) =
    !  &        SUM(tauaer_temp(NSOL+1,K,1:NCLOUDS)/(TAUAER(NSOL+1,JJ)+1e-8)*G0_TEMP(NSOL+1,K,1:NCLOUDS))
    !  &        + (TAU_HAZE(NSOL+1,K) * HAZE_PlanckMean_gg(temp_loc, haze_layer_index) / (TAUAER(NSOL+1,JJ) + 1e-8))


    !             TAUAER(NSOL+2,JJ) = SUM(tauaer_temp(NSOL+2,K,1:NCLOUDS)) + TAU_HAZE(NSOL+2,K)
    !             WOL(NSOL+2,JJ) =
    !  &        SUM(tauaer_temp(NSOL+2,K,1:NCLOUDS)/(TAUAER(NSOL+2,JJ)+1e-8)*PI0_TEMP(NSOL+2,K,1:NCLOUDS))
    !  &        + (TAU_HAZE(NSOL+2,K) * HAZE_RosselandMean_pi0(temp_loc, haze_layer_index) / (TAUAER(NSOL+2,JJ) + 1e-8))
    !             GOL(NSOL+2,JJ) =
    !  &        SUM(tauaer_temp(NSOL+2,K,1:NCLOUDS)/(TAUAER(NSOL+2,JJ)+1e-8)*G0_TEMP(NSOL+2,K,1:NCLOUDS))
    !  &        + (TAU_HAZE(NSOL+2,K) * HAZE_RosselandMean_gg(temp_loc, haze_layer_index) / (TAUAER(NSOL+2,JJ) + 1e-8))

    !             JJ = J+1
    !             DO L = NSOL+1,NTOTAL
    !                 TAUAER(L,JJ) = TAUAER(L,JJ-1)
    !                 WOL(L,JJ)    = WOL(L,JJ-1)
    !                 GOL(L,JJ)    = GOL(L,JJ-1)
    !             END DO
    !             k = k+1
    !         END DO
    !     END IF


        ! Smooth out the cloud properties after doubling
        DO L = NSOL+1,NTOTAL
            DO J = 2, NDBL-1, 2
                TAUAER(L,J) = (TAUAER(L,J+1) + TAUAER(L,J-1)) / 2.0
                WOL(L,J) = (WOL(L,J+1) + WOL(L,J-1)) / 2.0
                GOL(L,J) = (GOL(L,J+1) + GOL(L,J-1)) / 2.0
            END DO
        END DO

        ! ramp = 5.0  ! Set an appropriate value for ramping up clouds linearly
        ! ! Apply a ramp to the cloud properties
        ! IF (KOUNT/ITSPD .LT. ramp) THEN
        !   factor = (KOUNT/ramp)/ITSPD
        !   ! write(*,*) 'Ramping up the cloud properties by a factor of:', factor
        !   write(*,*) 'TAUAER before ramp:', TAUAER
        !   TAUAER = TAUAER * factor
        !   write(*,*) 'TAUAER after ramp:', TAUAER
        ! ENDIF
      END IF
    !   write(*,*) TAUAER ! should be 0 if aerosoltype is none
      iradgas = 1
      DO J = 1,NLAYER
          j1 = max(1, j-1)

            !         First the solar at standard resolution
          DO L = MAX(solar_calculation_indexer,MINWNOSTEL*8+1),NSOL
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
          DO L = NSOL+1,NTOTAL
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
              DO L = NSOL+1,NTOTAL
                  Y3(L,I,J) =   EXP(-TAUL(L,J)/GANGLE(I))
              END DO
          END DO
      END DO
      RETURN

      CONTAINS

      INTEGER FUNCTION NEAREST_INDEX(ARR, N, VAL)
          ! Explicit-loop nearest-value lookup. Replaces MINLOC(ABS(ARR-VAL),1),
          ! which forces the compiler to materialize a temporary array; under
          ! -recursive these temporaries can land in shared scratch storage and
          ! be clobbered by other OpenMP threads, corrupting the result.
          INTEGER, INTENT(IN) :: N
          REAL, INTENT(IN) :: ARR(N), VAL
          INTEGER :: KK
          REAL :: BESTDIFF, DIFFVAL
          NEAREST_INDEX = 1
          BESTDIFF = ABS(ARR(1) - VAL)
          DO KK = 2, N
              DIFFVAL = ABS(ARR(KK) - VAL)
              IF (DIFFVAL .LT. BESTDIFF) THEN
                  BESTDIFF = DIFFVAL
                  NEAREST_INDEX = KK
              ENDIF
          END DO
      END FUNCTION NEAREST_INDEX

      SUBROUTINE LOG_INTERP_WEIGHTS(ARR, N, VAL, IDX_LO, W)
          ! Locate VAL within the monotonically increasing, strictly
          ! positive array ARR(1:N) and return the lower bracketing index
          ! IDX_LO (1 <= IDX_LO <= N-1) and interpolation weight W in [0,1]
          ! such that, for any array Y defined on the same grid,
          ! Y(VAL) ~= (1-W)*Y(IDX_LO) + W*Y(IDX_LO+1), with W computed from
          ! the position of VAL between ARR(IDX_LO) and ARR(IDX_LO+1) in
          ! log-space. VAL outside [ARR(1),ARR(N)] is clamped to the nearest
          ! endpoint (W=0 or W=1, no extrapolation).
          INTEGER, INTENT(IN)  :: N
          REAL,    INTENT(IN)  :: ARR(N), VAL
          INTEGER, INTENT(OUT) :: IDX_LO
          REAL,    INTENT(OUT) :: W
          INTEGER :: KK

          IF (VAL .LE. ARR(1)) THEN
              IDX_LO = 1
              W = 0.0
              RETURN
          ELSE IF (VAL .GE. ARR(N)) THEN
              IDX_LO = N - 1
              W = 1.0
              RETURN
          ENDIF

          IDX_LO = N - 1
          DO KK = 1, N - 1
              IF (VAL .LT. ARR(KK+1)) THEN
                  IDX_LO = KK
                  EXIT
              ENDIF
          END DO

          W = (LOG(VAL) - LOG(ARR(IDX_LO))) / (LOG(ARR(IDX_LO+1)) - LOG(ARR(IDX_LO)))
      END SUBROUTINE LOG_INTERP_WEIGHTS

      REAL FUNCTION LERP(Y_LO, Y_HI, W)
          REAL, INTENT(IN) :: Y_LO, Y_HI, W
          LERP = Y_LO + W * (Y_HI - Y_LO)
      END FUNCTION LERP

      END

