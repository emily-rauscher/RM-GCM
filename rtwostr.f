      SUBROUTINE TWOSTR(TAUL,solar_calculation_indexer,
     &   NPROB, SOL, RAYPERBAR, WEIGHT, GOL, WOL, TAUCONST,
     &   WAVE, TT, Y3, PTEMPG, PTEMPT, G0, OPD,
     &   PTEMP, uG0, uTAUL, W0, uW0,
     &   uopd, U1S, U1I, TOON_AK, B1, B2, EE1,
     &   EM1, EM2, EL1, EL2, GAMI, AF,
     &   BF, EF, SFCS, B3, CK1, CK2,
     &   CP, CPB, CM, CMB, DIRECT,
     &   FNET, EE3,EL3, TMI, AS,
     &   DF, DS, XK, DIREC, DIRECTU,
     &   DINTENT, UINTENT, TMID, TMIU,
     &   firu, fird, fsLu, fsLd, fsLn, alb_toa, fupbs, fdownbs,
     &   fnetbs, fdownbs2, fupbi, fdownbi, fnetbi)
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
      INTEGER LLA, LLS, JDBLE, JDBLEDBLE, JN, JN2, iblackbody_above, ISL, IR, IRS
      REAL EPSILON, SOLNET, EMISIR, TPI, SQ3, SBK, AM, AVG, ALOS, SCDAY, RGAS
      REAL U0, FDEGDAY, WOT, GOT, tslu, total_downwelling, alb_tot,tiru, alb_tomi, alb_toai
      REAL HEATI(NL+1), HEATS(NL+1), HEAT(NL+1), GANGLE(3), GWEIGHT(3), GRATIO(3), EMIS(5), RSFX(5)
      REAL NPROB(5), SOL(5), RAYPERBAR(5), WEIGHT(5), GOL(5,2*NL+2), WOL(5,2*NL+2),TAUCONST(5)
      REAL WAVE(5+1), TT(NL+1), Y3(5,3,2*NL+2), PTEMPG(5), PTEMPT(5), G0(5,2*NL+2), OPD( 5,2*NL+2)
      REAL PTEMP(5,2*NL+2), uG0(5,2*NL+2), uTAUL(5,2*NL+2), W0(5,2*NL+2), uW0(5,2*NL+2)
      REAL uopd(5,2*NL+2), U1S(5), U1I(5), TOON_AK(5,2*NL+2), B1(5,2*NL+2), B2(5,2*NL+2), EE1(5,2*NL+2)
      REAL EM1(5,2*NL+2), EM2(5,2*NL+2), EL1(5,2*NL+2), EL2(5,2*NL+2), GAMI(5,2*NL+2), AF(5,4*NL+4)
      REAL BF(5,4*NL+4), EF(5,4*NL+4), SFCS(5), B3(5,2*NL+2), CK1(5,2*NL+2), CK2(5,2*NL+2)
      REAL CP(5,2*NL+2), CPB(5,2*NL+2), CM(5,2*NL+2), CMB(5,2*NL+2), DIRECT(5,2*NL+2)
      REAL FNET(5,2*NL+2), EE3(5,2*NL+2),EL3(5,2*NL+2), TMI(5,2*NL+2), AS(5,4*NL+4)
      REAL DF(5,4*NL+4), DS(5,4*NL+4), XK(5,4*NL+4), DIREC(5,2*NL+2), DIRECTU(5,2*NL+2)
      REAL DINTENT(5,3,2*NL+2), UINTENT(5,3,2*NL+2), TMID(5,2*NL+2), TMIU(5,2*NL+2)
      REAL firu(2), fird(2), fsLu(3), fsLd(3), fsLn(3), alb_toa(3), fupbs(NL+1), fdownbs(NL+1)
      REAL fnetbs(NL+1), fdownbs2(NL+1), fupbi(NL+1), fdownbi(NL+1), fnetbi(NL+1)


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
