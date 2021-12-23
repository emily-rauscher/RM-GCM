      SUBROUTINE OPPR1(TAUL, SLOPE, TTsub, t_pass)
!     &   NPROB, SOL, RAYPERBAR, WEIGHT, GOL, WOL, TAUCONST,
!     &   WAVE, TT, Y3, PTEMPG, PTEMPT, G0, OPD,
!     &   PTEMP, uG0, uTAUL, W0, uW0,
!     &   uopd, U1S, U1I, TOON_AK, B1, B2, EE1,
!     &   EM1, EM2, EL1, EL2, GAMI, AF,
!     &   BF, EF, SFCS, B3, CK1, CK2,
!     &   CP, CPB, CM, CMB, DIRECT,
!     &   FNET, EE3,EL3, TMI, AS,
!     &   DF, DS, XK, DIREC, DIRECTU,
!     &   DINTENT, UINTENT, TMID, TMIU,
!     &   firu, fird, fsLu, fsLd, fsLn, alb_toa, fupbs, fdownbs,
!     &   fnetbs, fdownbs2, fupbi, fdownbi, fnetbi)

!     **********************************************************
!     *  Purpose             :  Calculate Planck Function and  *
!     *                         and its derivative at ground   *
!     *                         and at all altitudes.          *
!     *  Subroutines Called  :  None                           *
!     *  Input               :  NLOW, WEIGHT            *
!     *  Output              :  PTEMP, PTEMPG, SLOPE           *
!     * ********************************************************
!
      include 'rcommons.h'
      real  ITP, ITG, IT1, SBKoverPI,g11,SBK
      DIMENSION  T(NLAYER), T_pass(NLAYER)
      real, dimension(5,2*NL+2) :: TAUL
      real, dimension(NTOTAL,NDBL) :: SLOPE
      real, dimension(NDBL) :: TTsub

!      INTEGER LLA, LLS, JDBLE, JDBLEDBLE, JN, JN2, iblackbody_above, ISL, IR, IRS
!      REAL EPSILON, SOLNET, EMISIR, TPI, SQ3, AM, AVG, ALOS, SCDAY, RGAS
!      REAL U0, FDEGDAY, WOT, GOT, tslu, total_downwelling, alb_tot,tiru, alb_tomi, alb_toai
!      REAL HEATI(NL+1), HEATS(NL+1), HEAT(NL+1), GANGLE(3), GWEIGHT(3), GRATIO(3), EMIS(5), RSFX(5)
!      REAL NPROB(5), SOL(5), RAYPERBAR(5), WEIGHT(5), GOL(5,2*NL+2), WOL(5,2*NL+2),TAUCONST(5)
!      REAL WAVE(5+1), TT(NL+1), Y3(5,3,2*NL+2), PTEMPG(5), PTEMPT(5), G0(5,2*NL+2), OPD( 5,2*NL+2)
!      REAL PTEMP(5,2*NL+2), uG0(5,2*NL+2), uTAUL(5,2*NL+2), W0(5,2*NL+2), uW0(5,2*NL+2)
!      REAL uopd(5,2*NL+2), U1S(5), U1I(5), TOON_AK(5,2*NL+2), B1(5,2*NL+2), B2(5,2*NL+2), EE1(5,2*NL+2)
!      REAL EM1(5,2*NL+2), EM2(5,2*NL+2), EL1(5,2*NL+2), EL2(5,2*NL+2), GAMI(5,2*NL+2), AF(5,4*NL+4)
!      REAL BF(5,4*NL+4), EF(5,4*NL+4), SFCS(5), B3(5,2*NL+2), CK1(5,2*NL+2), CK2(5,2*NL+2)
!      REAL CP(5,2*NL+2), CPB(5,2*NL+2), CM(5,2*NL+2), CMB(5,2*NL+2), DIRECT(5,2*NL+2)
!      REAL FNET(5,2*NL+2), EE3(5,2*NL+2),EL3(5,2*NL+2), TMI(5,2*NL+2), AS(5,4*NL+4)
!      REAL DF(5,4*NL+4), DS(5,4*NL+4), XK(5,4*NL+4), DIREC(5,2*NL+2), DIRECTU(5,2*NL+2)
!      REAL DINTENT(5,3,2*NL+2), UINTENT(5,3,2*NL+2), TMID(5,2*NL+2), TMIU(5,2*NL+2)
!      REAL firu(2), fird(2), fsLu(3), fsLd(3), fsLn(3), alb_toa(3), fupbs(NL+1), fdownbs(NL+1)
!      REAL fnetbs(NL+1), fdownbs2(NL+1), fupbi(NL+1), fdownbi(NL+1), fnetbi(NL+1)

      integer kindex
!     **************************************
!     * CALCULATE PTEMP AND SLOPE          *
!     **************************************

      SBK=5.6704E-8
      SBKoverPI=SBK/PI
      T=t_pass

      DO 300 J            =   1,NDBL
          kindex          = max(1,j-1)


          IT1 = TTsub(J)*TTsub(J)*TTsub(J)*TTsub(J)*SBKoverPI

          DO 200 L        = NSOLP+1,NTOTAL
              PTEMP(L,J)=IT1
              SLOPE(L,J)   = (PTEMP(L,J)-PTEMP(L,KINDEX))/TAUL(L,J)

              if( TAUL(L,J) .le. 1.0E-6 ) SLOPE(L,J) = 0.
 200      CONTINUE
 300  CONTINUE

      RETURN
      END

