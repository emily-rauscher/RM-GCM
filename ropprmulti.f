      SUBROUTINE OPPRMULTI(TAURAY,TAUL,TAUGAS,TAUAER,solar_calculation_indexer)
!
!     **************************************************************
!     *  Purpose             :  CaLculates optical properties      *
!     *                         such as single scattering albedo,  *
!     *                         asymmetry parameter, etc.          *
!     *                         This routine is case dependent and *
!     *                         wiLL have to be repLaced by the    *
!     *                         user.                              *
!     *  Subroutines Called  :  None                               *
!     *  Input               :  PAH2O, RDH2O, CO2, O3, ETC         *
!     *  Output              :  TAUL, W0, G0, OPD, Y3              *
!     * ************************************************************
!
      include 'rcommons.h'

      REAL TAUFACT

      REAL CONDFACT(NLAYER,NCLOUDS)
      REAL CLOUDLOC(NLAYER,NCLOUDS)

      REAL QE_TEMP(NSOL + NIR, NVERT, NCLOUDS)
      REAL PI0_TEMP(NSOL + NIR, NVERT, NCLOUDS)
      REAL G0_TEMP(NSOL + NIR, NVERT, NCLOUDS)
      REAL TAUAER_OPPR(NTOTAL, NLAYER, NCLOUDS)

      real, dimension(NIR+NSOL,2*NL+2) :: TAURAY,TAUL,TAUGAS,TAUAER

      ! These are hardcoded to 50 but they are just lookup tables
      ! Don't worry about expanding the GCM to more levels
      real, dimension(50) :: input_temperature_array
      real, dimension(50) :: input_pressure_array_cgs

      real, dimension(50) :: input_particle_size_array_in_meters
      real, dimension(50) :: particle_size_vs_layer_array_in_meters

      REAL QE_OPPR(NSOL + NIR, 50, 50, NCLOUDS)
      REAL PI0_OPPR(NSOL + NIR, 50, 50, NCLOUDS)
      REAL G0_OPPR(NSOL + NIR, 50, 50, NCLOUDS)

      REAL TCONDS(51,NCLOUDS)
      REAL CORFACT(51)

      REAL DENSITY(NCLOUDS)
      REAL FMOLW(NCLOUDS)
      REAL MOLEF(NCLOUDS)

      INTEGER K,J,BASELEV,TOPLEV,L
      INTEGER size_loc, temp_loc, pressure_loc,solar_calculation_indexer
      real particle_size

      COMMON /CLOUD_PROPERTIES/ TCONDS, QE_OPPR, PI0_OPPR, G0_OPPR,
     &                           DENSITY, FMOLW, MOLEF,
     &                           CORFACT,
     &                           input_particle_size_array_in_meters,
     &                           input_temperature_array,
     &                           particle_size_vs_layer_array_in_meters,
     &                           input_pressure_array_cgs

      DO J = 1,NLAYER-1
          ! Get the index of the closest pressure
          pressure_loc = MINLOC(ABS(input_pressure_array_cgs - (PRESS(J))),1)

          ! Using the pressure, find the particle size
          particle_size = particle_size_vs_layer_array_in_meters(pressure_loc)

          ! Find the location of the particle size
          ! Is this necessary, can I just use the pressure?
          size_loc = MINLOC(ABS(input_particle_size_array_in_meters - (particle_size)), 1)

          ! Get the array index of the closest temperature
          temp_loc = MINLOC(ABS(input_temperature_array - (TT(J))),1)

          DO I = 1,NCLOUDS
              DO L = 1,NTOTAL
                  QE_TEMP(L, J, I) = QE_OPPR(L, temp_loc,size_loc,I)
                  PI0_TEMP(L,J, I) = PI0_OPPR(L,temp_loc,size_loc,I)
                  G0_TEMP(L,J, I)  = G0_OPPR(L, temp_loc,size_loc,I)
              END DO
          END DO
      END DO

      DO I = 1,NCLOUDS
          DO J = 1,NLAYER-1
              pressure_loc = MINLOC(ABS(input_pressure_array_cgs - (PRESS(J))),1)
              temp_loc     = MINLOC(ABS(input_temperature_array - (TT(J))),1)

              particle_size = particle_size_vs_layer_array_in_meters(pressure_loc)
              size_loc = MINLOC(ABS(input_particle_size_array_in_meters - (particle_size)), 1)

              CONDFACT(J,I) = min(max((Tconds(pressure_loc,I)-TT(J)) / 10.0, 0.0), 1.0)

              TAUFACT = DPG(J)*10.*molef(I)*3./4./particle_size/density(I)*fmolw(I)*
     &                  CONDFACT(J,I)*MTLX*CORFACT(pressure_loc)

              DO L = 1,NTOTAL
                  TAUAER_OPPR(L,J,I) = TAUFACT*QE_OPPR(L,temp_loc,size_loc,I)
              END DO

              CLOUDLOC(J,I) = NINT(CONDFACT(J,I))*J
          END DO

          ! uncomment this section for compact cloud
          BASELEV = MAXVAL(CLOUDLOC(1:NVERT,I),1)
          TOPLEV  = max(BASELEV-AERLAYERS,0)  !changed from 1 to 0

          DO J = 1,TOPLEV
              DO L = 1,NTOTAL
                  TAUAER_OPPR(L,J,I) = 0.0
              END DO
          END DO

          ! MALSKY CODE
          DO L = 1,NTOTAL
              TAUAER_OPPR(L,TOPLEV+2,I) = TAUAER_OPPR(L,TOPLEV+2,I) * 0.135335
              TAUAER_OPPR(L,TOPLEV+1,I) = TAUAER_OPPR(L,TOPLEV+1,I) * 0.367879
          END DO
      END DO


!     SW AT STANDARD VERTICAL RESOLUTION
      DO L = solar_calculation_indexer,NSOLP
          DO J = 1,NLAYER
              TAUAER(L,J) = SUM(TAUAER_OPPR(L,J,1:NCLOUDS))
              WOL(L,J)    = SUM(TAUAER_OPPR(L,J,1:NCLOUDS)/(TAUAER(L,J)+1e-8) * PI0_TEMP(L,J,1:NCLOUDS))
              GOL(L,J)    = SUM(TAUAER_OPPR(L,J,1:NCLOUDS)/(TAUAER(L,J)+1e-8) * G0_TEMP(L, J,1:NCLOUDS))
          END DO
      END DO

!     LW AT 2X VERTICAL RESOLUTION (FOR PERFORMANCE).
      k = 1
      DO J = 1,NDBL,2
          JJ = J
          DO L = NSOLP+1,NTOTAL
              TAUAER(L,JJ) = SUM(TAUAER_OPPR(L,K,1:NCLOUDS))
              WOL(L,JJ) = SUM(TAUAER_OPPR(L,K,1:NCLOUDS)/(TAUAER(L,JJ)+1e-8)*PI0_TEMP(L,K,1:NCLOUDS))
              GOL(L,JJ) = SUM(TAUAER_OPPR(L,K,1:NCLOUDS)/(TAUAER(L,JJ)+1e-8)*G0_TEMP(L,K,1:NCLOUDS))
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

              if(TAUL(L,J) .lt. EPSILON ) then
                  TAUL(L,J) = EPSILON
              endif

              utauL(L,j)  = TAUL(L,J)
              WOT = (TAURAY(L,J)+TAUAER(L,J)*WOL(L,J))/TAUL(L,J)
              if (iradgas.eq.0) then
                  wot = woL(L,j)
              endif

              WOT       = min(1.-EPSILON,WOT)
              uw0(L,j)  = WOT
              DENOM     = (TAURAY(L,J) + TAUAER(L,J) * WOL(L,J))

              if( DENOM .LE. EPSILON ) then
                  DENOM = EPSILON
              endif

              if( DENOM .GT. EPSILON ) then
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



              IF (deltascale) THEN
                  FO          = GOT*GOT
                  DEN         = 1.-WOT*FO
                  TAUL(L,J)   = TAUL(L,J) * DEN
                  W0(L,J)     = (1.-FO)*WOT/DEN
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

              if( TAUL(L,J) .lt. EPSILON ) then
                  TAUL(L,J) = EPSILON
              endif

              utauL(L,j)  = TAUL(L,J)
              WOT         = (TAURAY(L,J)+TAUAER(L,J)*WOL(L,J))/TAUL(L,J)
              if (iradgas.eq.0) then
                  wot = woL(L,j)
              endif

              WOT         = min(1.-EPSILON,WOT)
              uw0(L,j)    = WOT
              DENOM       = (TAURAY(L,J)+ TAUAER(L,J)*WOL(L,J))

              if( DENOM .LE. EPSILON ) then
                  DENOM = EPSILON
              endif

              if( DENOM .GT. EPSILON ) then
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

              IF (deltascale) THEN
                  FO          = GOT*GOT
                  DEN         = 1.-WOT*FO
                  TAUL(L,J)   = TAUL(L,J) * DEN
                  W0(L,J)     = (1.-FO)*WOT/DEN
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

