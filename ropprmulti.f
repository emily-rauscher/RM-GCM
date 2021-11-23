      SUBROUTINE OPPRMULTI
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
      REAL TAUAER_OPPR(NTOTAL, NLAYER, NCLOUDS)

      REAL CONDFACT(NLAYER,NCLOUDS)
      REAL CLOUDLOC(NLAYER,NCLOUDS)

      REAL TCONDS(NLAYER,NCLOUDS)

      REAL QE_OPPR(NSOL + NIR, NVERT, NVERT, NCLOUDS)
      REAL PI0_OPPR(NSOL + NIR, NVERT, NVERT, NCLOUDS)
      REAL G0_OPPR(NSOL + NIR, NVERT, NVERT, NCLOUDS)

      REAL QE_TEMP(NSOL + NIR, NVERT, NCLOUDS)
      REAL PI0_TEMP(NSOL + NIR, NVERT, NCLOUDS)
      REAL G0_TEMP(NSOL + NIR, NVERT, NCLOUDS)

      REAL DENSITY(NCLOUDS)
      REAL FMOLW(NCLOUDS)
      REAL MOLEF(NCLOUDS)

      REAL CORFACT(NLAYER)

      real, dimension(NVERT) :: input_particle_size_array_in_meters
      real, dimension(NVERT) :: input_temperature_array
      real, dimension(NVERT) :: particle_size_vs_layer_array_in_meters

      INTEGER K,J,BASELEV,TOPLEV

!     The mass of these layers is this pressure thickness divide by G.
!     (NOTE - THE TOP LAYER IS FROM PTOP TO 0, SO AVERAGE = PTOP/2)
!     P_aerad=PRESSURE AT EDGES OF LAYERS (PASCALS) 
!     PRESS - PRESSURE AT EDGE OF LAYER (dyne/cm^2)
!     DPG   - MASS OF LAYER (G / CM**2)!    D PBAR 
!     PBARS - THICKNESS OF LAYER IN PRESSURE (BARS)
!     PRESSMID- PRESSURE AT CENTER OF LAYER (dyne/cm^2

      COMMON /CLOUD_PROPERTIES/ TCONDS, QE_OPPR, PI0_OPPR, G0_OPPR,
     &                              DENSITY, FMOLW, MOLEF,
     &                              CORFACT,
     &                              input_particle_size_array_in_meters,
     &                              input_temperature_array,
     &                              particle_size_vs_layer_array_in_meters

      DO J = 1,NLAYER -1
          size_loc = MINLOC(ABS(input_particle_size_array_in_meters -
     &                      (particle_size_vs_layer_array_in_meters(J))), 1)

          ! Get the array index of the closest temperature
          temp_loc = MINLOC(ABS(input_temperature_array - (TT(J))),1)

          DO I = 1,NCLOUDS
              DO L = 1,NTOTAL
                  QE_TEMP(L, J, I) = QE_OPPR(L,temp_loc,size_loc,I)
                  PI0_TEMP(L,J, I) = PI0_OPPR(L,temp_loc,size_loc,I)
                  G0_TEMP(L, J, I) = G0_OPPR(L,temp_loc,size_loc, I)
              END DO
          END DO
      END DO


      DO I = 1,NCLOUDS
          DO J = 1,NLAYER-1
              CONDFACT(J,I) = min(max((Tconds(J,I)-TT(J))/10.,0.0),1.0)
              TAUFACT = DPG(J)*10.*molef(I)*3./4./ particle_size_vs_layer_array_in_meters(J)
     &                  /density(I)*fmolw(I)*CONDFACT(J,I)*MTLX*CORFACT(J)

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
              TAUAER_OPPR(L,TOPLEV+2,I) = TAUAER_OPPR(L,TOPLEV+2,I) * 0.367879
              TAUAER_OPPR(L,TOPLEV+1,I) = TAUAER_OPPR(L,TOPLEV+1,I) * 0.135335
          END DO
      END DO

      DO L =1, NTOTAL
          DO I =1, NCLOUDS
              TAUAER_OPPR(L,NLAYER,I) = 0.0
          END DO
      END DO

!     SW AT STANDARD VERTICAL RESOLUTION
      DO L = LLS,NSOLP
          DO J = 1,NLAYER
              TAUAER(L,J) = SUM(TAUAER_OPPR(L,J,1:NCLOUDS))

              WOL(L,J) = SUM(TAUAER_OPPR(L,J,1:NCLOUDS)/(SUM(TAUAER_OPPR(L,J,1:NCLOUDS))+1e-8) *
     &                       PI0_TEMP(L,J,1:NCLOUDS))
              GOL(L,J) = SUM(TAUAER_OPPR(L,J,1:NCLOUDS)/(SUM(TAUAER_OPPR(L,J,1:NCLOUDS))+1e-8) *
     &                       G0_TEMP(L, J,1:NCLOUDS))

          END DO
      END DO


!     LW AT 2X VERTICAL RESOLUTION (FOR PERFORMANCE).
      k = 1
      DO J = 1,NDBL,2
          JJ = J
          DO L = NSOLP+1,NTOTAL
              TAUAER(L,JJ) = SUM(TAUAER_OPPR(L,K,1:NCLOUDS))
              WOL(L,JJ)    = SUM(TAUAER_OPPR(L,K,1:NCLOUDS)/(SUM(TAUAER_OPPR(L,K,1:NCLOUDS))+1e-8) *
     &                        PI0_TEMP(L,K,1:NCLOUDS))
              GOL(L,JJ)    = SUM(TAUAER_OPPR(L,K,1:NCLOUDS)/(SUM(TAUAER_OPPR(L,K,1:NCLOUDS))+1e-8) *
     &                        G0_TEMP(L,K,1:NCLOUDS))
          END DO
          JJ = J+1
          DO L = NSOLP+1,NTOTAL
              TAUAER(L,JJ) = TAUAER(L,JJ-1)
              WOL(L,JJ)    = WOL(L,JJ-1)
              GOL(L,JJ)    = GOL(L,JJ-1)
          END DO
          k = k+1
      END DO

      iradgas = 1
      DO J = 1,NLAYER
          j1 = max(1, j-1)

!         First the solar at standard resolution
          DO L = LLS,NSOLP
              TAUL(L,J) = TAUGAS(L,J)+TAURAY(L,J)+TAUAER(L,J)

              if( TAUL(L,J) .lt. EPSILON ) then
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

          IF(IR .EQ. 1) THEN
              DO I = 1,NGAUSS
                  DO L = LLS,LLA
                      Y3(L,I,J) =   EXP(-TAUL(L,J)/GANGLE(I))
                  END DO
              END DO
          END IF
      END DO

      RETURN
      END

