      SUBROUTINE OPPR1
!
!     **********************************************************
!     *  Purpose             :  Calculate Planck Function and  *
!     *                         and its derivative at ground   *
!     *                         and at all altitudes.          *
!     *  Subroutines Called  :  None                           *
!     *  Input               :  TGRND, NLOW, WEIGHT            *
!     *  Output              :  PTEMP, PTEMPG, SLOPE           *
!     * ********************************************************
!
      include 'rcommons.h'
      real  ITP, ITG, IT1, SBK, SBKoverPI,g11
      DIMENSION  PTEMP2(NTOTAL-NSOLP),PLTEMP1(NTOTAL-NSOLP)
      DIMENSION  T(NLAYER)
      integer kindex
!     **************************************
!     * CALCULATE PTEMP AND SLOPE          *
!     **************************************
!
!     CALCULATE THE WAVELENGTH DEPENDENT PLANK FUNCTION AT THE GROUND.

      SBK=5.6704E-8
      SBKoverPI=SBK/PI
      T=t_aerad
!     NOTE TO ERIN:
!     Here is where the ground temperature is treated seperately, currently
!     comented out since it is not defined seperately.
!     ITG --  temperature at the ground, and L is the index over wavelength bins
!     (ntotal = 2 bins, Nsolp = number of solar bins i.e.1)
!     PTEMPG -- Plank temperature of the ground layer (i.e. sigma T^4 /pi)

!    THE CODE BELOW IS A MESS. IT DEALS WITH THE
!     if( iblackbody_above .ne. 0 )then
!
!      CALCULATE THE WAVELENGTH DEPENDENT PLANK FUNCTION AT THE TOP
!      OF THE MODEL, BUT I SIMPLY SET IT UP SO THE LOOPS BELOW ARE
!      SUFFICIENT.

      DO J = 1,NDBL ! NLAYER ! MTR
          kindex = max( 1, j-1 )
          IT1 = TTsub(J)*TTsub(J)*TTsub(J)*TTsub(J)*SBKoverPI

!         KINDEX MAKES AS DEFINED ABOVE MAKES THE TOP LAYER ISOTHERMAL;
!         BELOW THE TOP SLOPE IS REDIFINED USING THE EXTRAPOLATED
!         VALUE. FIND PLANK FUNCTION AT BOTTOM OF EACH LAYER.
!         NOTE: IF YOU FORCE SLOPE=0, THEN YOU HAVE ISOTHERMAL
!         LAYERS WITH TT(J) CORRESPONDING TO AVERAGE TEMPERATURE
!         OF LAYER AND TT(NLAYER) SHOULD BE SET TO TGRND.

          DO L = NSOLP+1,NTOTAL
              PTEMP(L,J)=IT1
              SLOPE(L,J)   = (PTEMP(L,J)-PTEMP(L,KINDEX))/TAUL(L,J)
              IF( TAUL(L,J) .le. 1.0E-6 ) THEN
                  SLOPE(L,J) = 0.
              END IF
          END DO
      END DO

      RETURN
      END

