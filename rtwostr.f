      SUBROUTINE TWOSTR(TAUL,solar_calculation_indexer)
!
!    ******************************************************************
!    *  Purpose             :  Defines matrix properties and sets up  *
!    *                         matrix coefficients that do not depend *
!    *                         on zenith angle or temperature.        *
!    *  Subroutines Called  :  None                                   *
!    *  Input               :  W0, G0                                 *
!    * ****************************************************************
!
      include 'rcommons.h'

      real, dimension(5,2*NL+2) :: TAUL
      integer solar_calculation_indexer

       DO 10 L    =  solar_calculation_indexer,LLA
          if( L .LE. NSOLP )then
            U1I(L) = SQ3  !2.d0 !SQ3
          else
            U1I(L) = 2.d0
          endif
  10      U1S(L)  =  TPI/U1I(L)

!      HERE WE DEFINE LAYER PROPERTIES FOLLOWING GENERAL SCHEME
!      OF MEADOR AND WEAVOR. THEN WE SET UP LAYER PROPERTIES
!      NEEDED FOR MATRIX.
!
       DO 14 J          =  1,NLAYER
          DO 14 L       =  solar_calculation_indexer,NSOLP
!            THESE ARE FOR TWO STREAM AND HEMISPHERIC MEANS
             B1(L,J)    =  0.5*U1I(L)*(2. - W0(L,J)*(1. + G0(L,J)))
             B2(L,J)    =  0.5*U1I(L)*W0(L,J)*(1. - G0(L,J))
             TOON_AK(L,J)    = SQRT(ABS(B1(L,J)*B1(L,J) - B2(L,J)*B2(L,J)))
             GAMI(L,J)  =  B2(L,J)/(B1(L,J) + TOON_AK(L,J))
             EE1(L,J)   =  EXP(-TOON_AK(L,J)*TAUL(L,J))
             EL1(L,J)   =  1.0 + GAMI(L,J) *EE1(L,J)  !e1
             EM1(L,J)   =  1.0 - GAMI(L,J) * EE1(L,J) !e2
             EL2(L,J)   =  GAMI(L,J) + EE1(L,J)       !e3
             EM2(L,J)   =  GAMI(L,J) - EE1(L,J)       !e4
  14  CONTINUE

       DO 15 J          =  1,NDBL
          DO 15 L       =  NSOLP+1,NTOTAL
!            THESE ARE FOR TWO STREAM AND HEMISPHERIC MEANS
             B1(L,J)    =  0.5*U1I(L)*(2. - W0(L,J)*(1. + G0(L,J)))
             B2(L,J)    =  0.5*U1I(L)*W0(L,J)*(1. - G0(L,J))
             TOON_AK(L,J)    = SQRT(ABS(B1(L,J)*B1(L,J) - B2(L,J)*B2(L,J)))
             GAMI(L,J)  =  B2(L,J)/(B1(L,J) + TOON_AK(L,J))
             EE1(L,J)   =  EXP(-TOON_AK(L,J)*TAUL(L,J))

             EL1(L,J)   =  1.0 + GAMI(L,J) *EE1(L,J)  !e1
             EM1(L,J)   =  1.0 - GAMI(L,J) * EE1(L,J) !e2
             EL2(L,J)   =  GAMI(L,J) + EE1(L,J)       !e3
             EM2(L,J)   =  GAMI(L,J) - EE1(L,J)       !e4
  15  CONTINUE

!
!     WE SEEK TO SOLVE AX(L-1)+BX(L)+EX(L+1) = D.
!     L=2N FOR EVEN L, L=N+1 FOR ODD L. THE MEAN INTENSITY (TMI/4PI)
!     AND THE NET FLUX (FNET) ARE RELATED TO X'S AS NOTED IN ADD.
!     FIRST WE SET UP THE COEFFICIENTS THAT ARE INDEPENDENT OF SOLAR
!     ANGLE OR TEMPARATURE: A(I),B(I),E(I). D(I) IS DEFINED IN ADD.
!
      J                 =  0
      DO 18 JD          =  2,JN,2
         J              =  J + 1
         DO 18 L        =  solar_calculation_indexer,NSOLP
!          HERE ARE THE EVEN MATRIX ELEMENTS
             AF(L,JD)   =  EM1(L,J+1)*EL1(L,J)-EM2(L,J+1)*EL2(L,J)
             BF(L,JD)   =  EM1(L,J+1)* EM1(L,J)-EM2(L,J+1)*EM2(L,J)
             EF(L,JD)  = EL1(L,J+1)*EM2(L,J+1) - EL2(L,J+1)*EM1(L,J+1)
!          HERE ARE THE ODD MATRIX ELEMENTS EXCEPT FOR THE TOP.
             AF(L,JD+1) =  EM1(L,J)*EL2(L,J)-EL1(L,J)*EM2(L,J)
             BF(L,JD+1) =  EL1(L,J+1)*EL1(L,J) - EL2(L,J+1)*EL2(L,J)
             EF(L,JD+1) =  EL2(L,J)*EM2(L,J+1)-EL1(L,J)*EM1(L,J+1)
  18  CONTINUE
      J                 =  0
      DO 19 JD          =  2,JN2,2
         J              =  J + 1
         DO 19 L        =  NSOLP+1,NTOTAL
!          HERE ARE THE EVEN MATRIX ELEMENTS
             AF(L,JD)   =  EM1(L,J+1)*EL1(L,J)-EM2(L,J+1)*EL2(L,J)
             BF(L,JD)   =  EM1(L,J+1)* EM1(L,J)-EM2(L,J+1)*EM2(L,J)
             EF(L,JD)  = EL1(L,J+1)*EM2(L,J+1) - EL2(L,J+1)*EM1(L,J+1)
!          HERE ARE THE ODD MATRIX ELEMENTS EXCEPT FOR THE TOP.
             AF(L,JD+1) =  EM1(L,J)*EL2(L,J)-EL1(L,J)*EM2(L,J)
             BF(L,JD+1) =  EL1(L,J+1)*EL1(L,J) - EL2(L,J+1)*EL2(L,J)
             EF(L,JD+1) =  EL2(L,J)*EM2(L,J+1)-EL1(L,J)*EM1(L,J+1)
  19  CONTINUE

!
!     HERE ARE THE TOP AND BOTTOM BOUNDARY CONDITIONS AS WELL AS THE
!     BEGINNING OF THE TRIDIAGONAL SOLUTION DEFINITIONS. I ASSUME
!     NO DIFFUSE RADIATION IS INCIDENT AT UPPER BOUNDARY.
!
      DO 20 L        = solar_calculation_indexer,NSOLP
         AF(L,1)     = 0.0
         BF(L,1)     = EL1(L,1)
         EF(L,1)     = -EM1(L,1)
         AF(L,JDBLE) = EL1(L,NLAYER)-RSFX(L)*EL2(L,NLAYER)
         BF(L,JDBLE) = EM1(L,NLAYER)-RSFX(L)*EM2(L,NLAYER)
         EF(L,JDBLE) = 0.0
  20  CONTINUE
      DO 21 L        = NSOLP+1,NTOTAL
         AF(L,1)     = 0.0
         BF(L,1)     = EL1(L,1)
         EF(L,1)     = -EM1(L,1)
         AF(L,JDBLEDBLE) = EL1(L,NDBL)-RSFX(L)*EL2(L,NDBL)
         BF(L,JDBLEDBLE) = EM1(L,NDBL)-RSFX(L)*EM2(L,NDBL)
         EF(L,JDBLEDBLE) = 0.0
  21  CONTINUE
      RETURN
      END
