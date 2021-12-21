      SUBROUTINE OPPR1(TAUL, SLOPE)
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
      real, dimension(5,2*NL+2) :: TAUL
      real, dimension(NTOTAL,NDBL) :: SLOPE

      integer kindex
!     **************************************
!     * CALCULATE PTEMP AND SLOPE          *
!     **************************************

      SBK=5.6704E-8
      SBKoverPI=SBK/PI
      T=t_aerad

      DO 300 J            =   1,NDBL ! NLAYER ! MTR
          kindex          = max( 1, j-1 )
          IT1 = TTsub(J)*TTsub(J)*TTsub(J)*TTsub(J)*SBKoverPI

          DO 200 L        = NSOLP+1,NTOTAL
              PTEMP(L,J)=IT1
              SLOPE(L,J)   = (PTEMP(L,J)-PTEMP(L,KINDEX))/TAUL(L,J)

              if( TAUL(L,J) .le. 1.0E-6 ) SLOPE(L,J) = 0.
 200      CONTINUE
 300  CONTINUE

      RETURN
      END

