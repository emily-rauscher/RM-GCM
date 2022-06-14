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
     &  PI0_TEMP, G0_TEMP, tauaer_temp,j1,denom)
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

      INTEGER LLA, LLS, JDBLE, JDBLEDBLE, JN, JN2, iblackbody_above, ISL, IR, IRS, j1
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
      REAL DPG(NLAYER), p_pass(NLAYER)
      REAL CONDFACT(NLAYER,NCLOUDS)

      REAL PI0_TEMP(NSOL + NIR, NVERT, NCLOUDS)
      REAL G0_TEMP(NSOL + NIR, NVERT, NCLOUDS)
      REAL tauaer_temp(NTOTAL, NLAYER, NCLOUDS)

      REAL CLOUDLOC(NL+1,NCLOUDS)
      INTEGER BASELEV
      INTEGER TOPLEV(NCLOUDS)

      real, dimension(NIR+NSOL,2*NL+2) :: TAURAY,TAUL,TAUGAS,TAUAER
      real, dimension(NIR+NSOL,2*NL+2) :: TAU_HAZE

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
      real, dimension(50, 50) :: HAZE_Rosseland_tau_per_bar, HAZE_Rosseland_pi0, HAZE_Rosseland_gg
      real, dimension(50) :: HAZE_500nm_tau_per_bar
      real, dimension(50) :: HAZE_500nm_pi0, HAZE_500nm_gg
      real, dimension(50) :: HAZE_650nm_tau_per_bar, HAZE_650nm_pi0, HAZE_650nm_gg
      real, dimension(50) :: HAZE_800nm_tau_per_bar, HAZE_800nm_pi0, HAZE_800nm_gg
      real, dimension(50) :: HAZE_5000nm_tau_per_bar, HAZE_5000nm_pi0, HAZE_5000nm_gg
      real, dimension(50) :: haze_pressure_array_pascals

      REAL TCONDS(51,NCLOUDS)
      REAL CORFACT(51)

      REAL DENSITY(NCLOUDS)
      REAL FMOLW(NCLOUDS)
      REAL MOLEF(NCLOUDS)

      INTEGER WAVELENGTH_INDEXES(NSOL + NIR)
      INTEGER WAV_LOC

      INTEGER K,J,L, iradgas
      INTEGER size_loc, temp_loc, solar_calculation_indexer, layer_index, haze_layer_index
      real particle_size


      COMMON /CLOUD_PROPERTIES/ TCONDS, QE_OPPR, PI0_OPPR, G0_OPPR,
     &                              DENSITY, FMOLW,
     &                              CORFACT,
     &                              input_particle_size_array_in_meters,
     &                              input_temperature_array,
     &                              particle_size_vs_layer_array_in_meters,
     &                              input_pressure_array_cgs,
     &                              HAZE_Rosseland_tau_per_bar,HAZE_Rosseland_pi0, HAZE_Rosseland_gg,
     &                              HAZE_500nm_tau_per_bar,
     &                              HAZE_500nm_pi0, HAZE_500nm_gg,
     &                              HAZE_650nm_tau_per_bar, HAZE_650nm_pi0, HAZE_650nm_gg,
     &                              HAZE_800nm_tau_per_bar, HAZE_800nm_pi0, HAZE_800nm_gg,
     &                              HAZE_5000nm_tau_per_bar, HAZE_5000nm_pi0, HAZE_5000nm_gg,
     &                              haze_pressure_array_pascals


      Y3(:,:,:) = 0.0

      ! These three correspond to the wavelengths in the wavelength dependent scattering parameters
      ! 1 is 0.51 microns
      ! 2 is 0.64 microns
      ! 3 is 0.80 microns
      ! 4 is 4.89 microns
      WAVELENGTH_INDEXES = (/8, 11, 14, 38, 38/)

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!         GET THE HAZE DATA FIRST       !!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      write(*,*) INCLUDE_HAZES

      ! Do the starlight at 1x resolution
      IF (INCLUDE_HAZES) THEN
          DO J = 1, NLAYER
              haze_layer_index = MINLOC(ABS((haze_pressure_array_pascals) - (p_pass(J))),1)  ! Both of these are in PA

              ! This grabs the optical depth per bar, then multiply it by the pressure in bars
              TAU_HAZE(1,J) = HAZE_500nm_tau_per_bar(haze_layer_index) * (p_pass(J) * 1e-5)
              TAU_HAZE(2,J) = HAZE_650nm_tau_per_bar(haze_layer_index) * (p_pass(J) * 1e-5)
              TAU_HAZE(3,J) = HAZE_800nm_tau_per_bar(haze_layer_index) * (p_pass(J) * 1e-5)
          END DO

          ! Do the thermal at 2x resolution
          DO J = 1, NLAYER
              haze_layer_index = MINLOC(ABS((haze_pressure_array_pascals) - (p_pass(J))),1)  ! Both of these are in PA
              temp_loc         = MINLOC(ABS(input_temperature_array - (TT(J))),1) ! Not needed for the stellar calc

              ! This grabs the optical depth per bar, then multiply it by the pressure in bars
              DO L = NSOLP+1,NTOTAL
                  TAU_HAZE(L,2*J-1) = HAZE_Rosseland_tau_per_bar(temp_loc, haze_layer_index) * (p_pass(J) * 1e-5)
                  TAU_HAZE(L,2*J)   = HAZE_Rosseland_tau_per_bar(temp_loc, haze_layer_index) * (p_pass(J) * 1e-5)
              END DO
          END DO
      ELSE
           TAU_HAZE = 0.0
      END IF

      write(*,*) TAU_HAZE(1,:)
      stop





      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!         CLOUD SCATTERING PROPERTIES       !!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      DO J = 1,NLAYER
          layer_index   = MINLOC(ABS(input_pressure_array_cgs - (p_pass(J) * 10.0)),1)
          temp_loc      = MINLOC(ABS(input_temperature_array - (TT(J))),1)
          particle_size = particle_size_vs_layer_array_in_meters(layer_index) ! Convert to CGS
          size_loc      = MINLOC(ABS(input_particle_size_array_in_meters - (particle_size)), 1)

          DO I = 1,NCLOUDS
              DO L = solar_calculation_indexer,NSOLP
                  WAV_LOC = WAVELENGTH_INDEXES(L)
                  PI0_TEMP(L,J,I) = PI0_OPPR(L,WAV_LOC,size_loc,I)
                  G0_TEMP(L,J,I)  = G0_OPPR(L,WAV_LOC,size_loc,I)
              END DO

              DO L = NSOLP+1,NTOTAL
                  PI0_TEMP(L,J,I) = PI0_OPPR(L,temp_loc,size_loc,I)
                  G0_TEMP(L,J,I)  = G0_OPPR(L,temp_loc,size_loc,I)
              END DO

              CONDFACT(J,I) = min(max((Tconds(layer_index,I)-TT(J))/10.,0.0),1.0)

              ! STOP some weird behaviour, I don't know if this should be taken out. Probably
              !IF (J .gt. 5) THEN
              !    IF ((CONDFACT(J-1,I) .eq. 0) .AND. (CONDFACT(J-2,I) .eq. 0) .AND. (CONDFACT(J-3,I) .eq. 0)) THEN
              !        CONDFACT(J,I) = 0.0
              !    END IF
              !END IF

              CLOUDLOC(J,I) = NINT(CONDFACT(J,I))*J
              BASELEV = MAXVAL(CLOUDLOC(1:50,I),1)
              TOPLEV(I)  = max(BASELEV-AERLAYERS,0)

              ! DPG is CGS before that 10x
              DO L = solar_calculation_indexer,NSOLP
                  WAV_LOC = WAVELENGTH_INDEXES(L)
                  tauaer_temp(L,J,I) = (DPG(J)*10.0)*molef(I)*3./4./particle_size/density(I)*fmolw(I)*
     &                              CONDFACT(J,I)*MTLX*CORFACT(layer_index)*QE_OPPR(L,WAV_LOC,size_loc,I)
              END DO

              DO L = NSOLP+1,NTOTAL
                  tauaer_temp(L,J,I) = (DPG(J)*10.0)*molef(I)*3./4./particle_size/density(I)*fmolw(I)*
     &                                  CONDFACT(J,I)*MTLX*CORFACT(layer_index)*QE_OPPR(L,temp_loc,size_loc,I)
              END DO
          END DO
      END DO

      ! Uncomment for compact clouds I think
      DO I = 1,NCLOUDS
          DO J = 1, TOPLEV(I)
              tauaer_temp(:,J,I) = 0.0
          END DO
      END DO


!     SW AT STANDARD VERTICAL RESOLUTION
      DO J = 1,NLAYER
          haze_layer_index = MINLOC(ABS((haze_pressure_array_pascals) - (p_pass(J))),1)  ! Both of these are in pascals

          TAUAER(1,J) = SUM(tauaer_temp(1,J,1:NCLOUDS)) + TAU_HAZE(1,J)
          WOL(1,J)    = SUM(tauaer_temp(1,J,1:NCLOUDS)/(TAUAER(1,J)+1e-8) * PI0_TEMP(1,J,1:NCLOUDS))
     &                    + (TAU_HAZE(1,J) * HAZE_500nm_pi0(haze_layer_index) / TAUAER(1,J))
          GOL(1,J)    = SUM(tauaer_temp(1,J,1:NCLOUDS)/(TAUAER(1,J)+1e-8) * G0_TEMP(1,J,1:NCLOUDS))
     &                    + (TAU_HAZE(1,J) * HAZE_500nm_gg(haze_layer_index)  / TAUAER(1,J))

          TAUAER(2,J) = SUM(tauaer_temp(2,J,1:NCLOUDS)) + TAU_HAZE(2,J)
          WOL(2,J)    = SUM(tauaer_temp(2,J,1:NCLOUDS)/(TAUAER(2,J)+1e-8) * PI0_TEMP(2,J,1:NCLOUDS))
     &                   + (TAU_HAZE(2,J) * HAZE_650nm_pi0(haze_layer_index) / TAUAER(2,J))
          GOL(2,J)    = SUM(tauaer_temp(2,J,1:NCLOUDS)/(TAUAER(2,J)+1e-8) * G0_TEMP(2,J,1:NCLOUDS))
     &                   + (TAU_HAZE(2,J) * HAZE_650nm_gg(haze_layer_index)  / TAUAER(2,J))

          TAUAER(3,J) = SUM(tauaer_temp(3,J,1:NCLOUDS)) + TAU_HAZE(3,J)
          WOL(3,J)    = SUM(tauaer_temp(3,J,1:NCLOUDS)/(TAUAER(3,J)+1e-8) * PI0_TEMP(3,J,1:NCLOUDS))
     &                   + (TAU_HAZE(3,J) * HAZE_800nm_pi0(haze_layer_index) / TAUAER(3,J))
          GOL(3,J)    = SUM(tauaer_temp(3,J,1:NCLOUDS)/(TAUAER(3,J)+1e-8) * G0_TEMP(3,J,1:NCLOUDS))
     &                   + (TAU_HAZE(3,J) * HAZE_800nm_gg(haze_layer_index)  / TAUAER(3,J))
      END DO


!     LW AT 2X VERTICAL RESOLUTION (FOR PERFORMANCE).
      k = 1
      DO J = 1,NDBL,2
          JJ = J

          haze_layer_index = MINLOC(ABS((haze_pressure_array_pascals) - (p_pass(J))),1)  ! Both of these are in pascals
          temp_loc         = MINLOC(ABS(input_temperature_array - (TT(J))),1) ! Not needed for the stellar calc

          DO L = NSOLP+1,NTOTAL
              TAUAER(L,J) = SUM(tauaer_temp(L,J,1:NCLOUDS)) + TAU_HAZE(L,J)
              WOL(L,J)    = SUM(tauaer_temp(L,J,1:NCLOUDS)/(TAUAER(1,J)+1e-8) * PI0_TEMP(L,J,1:NCLOUDS))
     &                       + (TAU_HAZE(L,J) * HAZE_Rosseland_pi0(temp_loc, haze_layer_index) / TAUAER(L,J))
              GOL(L,J)    = SUM(tauaer_temp(L,J,1:NCLOUDS)/(TAUAER(1,J)+1e-8) * G0_TEMP(L,J,1:NCLOUDS))
     &                       + (TAU_HAZE(L,J) * HAZE_Rosseland_gg(temp_loc, haze_layer_index)  / TAUAER(L,J))
          END DO
          JJ = J+1
          DO L = NSOLP+1,NTOTAL
              TAUAER(L,JJ) = TAUAER(L,JJ-1)
              WOL(L,JJ)    = WOL(L,JJ-1)
              GOL(L,JJ)    = GOL(L,JJ-1)
          END DO
          k = k+1
      END DO

      ! Smooth out the cloud properties after doubling
      DO L = NSOLP+1,NTOTAL
          DO J = 2, NDBL, 2
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

              WOT       = min(1.-1d-6,WOT)
              uw0(L,j)  = WOT
              DENOM     = (TAURAY(L,J) + TAUAER(L,J) * WOL(L,J))

              if(DENOM .LE. 1d-6 ) then
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
              uOPD(L,J)   = uOPD(L,J1)+uTAUL(L,J)

              IF (.TRUE.) THEN
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

              if (iradgas.eq.0) then
                  tauL(L,j) = tauaer(L,j)
              endif

              if( TAUL(L,J) .lt. 1d-6 ) then
                  TAUL(L,J) = 1d-6
              endif

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

              IF (.TRUE.) THEN
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

