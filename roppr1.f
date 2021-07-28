      SUBROUTINE OPPR1(Beta_IR)
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
      real, dimension(2) :: Beta_IR
      DIMENSION  PTEMP2(NTOTAL-NSOLP),PLTEMP1(NTOTAL-NSOLP)
      DIMENSION  T(NLAYER)
      integer kindex

!     CALCULATE THE WAVELENGTH DEPENDENT PLANK FUNCTION AT THE GROUND.

      SBK=5.6704E-8
      SBKoverPI=SBK/PI

      T=t_aerad
      DO 300 J            =   1,NDBL
          kindex          = max( 1, j-1 )

          Beta_IR(1) = 1
          IT1 = TTsub(J)*TTsub(J)*TTsub(J)*TTsub(J)*SBKoverPI * Beta_IR(1)

          DO 200 L = NSOLP+1,NTOTAL
              PTEMP(L,J)  = IT1
              SLOPE(L,J)  = (PTEMP(L,J)-PTEMP(L,KINDEX))/TAUL(L,J)
              if( TAUL(L,J) .le. 1.0E-6 ) SLOPE(L,J) = 0.
 200      CONTINUE
 300  CONTINUE
      RETURN
      END

