      SUBROUTINE OPPR1(TAUL, SLOPE, t,
     &             LLA, LLS, JDBLE, JDBLEDBLE, JN, JN2, iblackbody_above, ISL, IR, IRS, EMISIR,
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
     &  qrad,alb_tomi,alb_toai, num_layers)
!
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

      integer kindex, J, L, num_layers, K, index_num

      INTEGER LLA, LLS, JDBLE, JDBLEDBLE, JN, JN2, iblackbody_above, ISL, IR, IRS
      REAL EMISIR, EPSILON, HEATI(NLAYER), HEATS(NLAYER), HEAT(NLAYER), SOLNET
      REAL TPI, SQ3, SBK,AM, AVG, ALOS
      REAL SCDAY, RGAS, GANGLE(3), GWEIGHT(3), GRATIO(3), EMIS(5), RSFX(5),NPROB(5), SOL(5),RAYPERBAR(5),WEIGHT(5)
      REAL GOL(5,2*NL+2), WOL(5,2*NL+2), WAVE(5+1), TT(NL+1), Y3(5,3,2*NL+2), U0, FDEGDAY
      REAL WOT, GOT, PTEMPG(5), PTEMPT(5), G0(5,2*NL+2), OPD( 5,2*NL+2), PTEMP(5,2*NL+2)
      REAL uG0(5,2*NL+2), uTAUL(5,2*NL+2), W0(5,2*NL+2), uW0(5,2*NL+2), uopd(5,2*NL+2),  U1S( 5)
      REAL U1I(5), TOON_AK(5,2*NL+2), B1(5,2*NL+2), B2(5,2*NL+2), EE1( 5,2*NL+2), EM1(5,2*NL+2)
      REAL EM2(5,2*NL+2), EL1( 5,2*NL+2), EL2(5,2*NL+2), GAMI(5,2*NL+2), AF(5,4*NL+4)
      REAL BF(5,4*NL+4), EF(5,4*NL+4), SFCS(5), B3(5,2*NL+2), CK1(5,2*NL+2), CK2(5,2*NL+2)
      REAL CP(5,2*NL+2), CPB(5,2*NL+2), CM(5,2*NL+2), CMB(5,2*NL+2), DIRECT(5,2*NL+2), EE3(5,2*NL+2)
      REAL EL3(5,2*NL+2), FNET(5,2*NL+2), TMI(5,2*NL+2), AS(5,4*NL+4), DF(5,4*NL+4)
      REAL DS(5,4*NL+4), XK(5,4*NL+4), DIREC(5,2*NL+2), DIRECTU(5,2*NL+2), DINTENT(5,3,2*NL+2)
      REAL UINTENT(5,3,2*NL+2), TMID(5,2*NL+2), TMIU(5,2*NL+2), tslu,total_downwelling,alb_tot
      REAL tiru,firu(2),fird(2),fsLu(3), fsLd(3),fsLn(3),alb_toa(3), fupbs(NL+1)
      REAL fdownbs(NL+1),fnetbs(NL+1),fdownbs2(NL+1), fupbi(NL+1),fdownbi(NL+1),fnetbi(NL+1)
      REAL qrad(NL+1),alb_tomi,alb_toai

      real  ITP, ITG, IT1, SBKoverPI,g11
      real, DIMENSION(NLAYER) :: T
      real, dimension(5,2*NL+2) :: TAUL
      real, dimension(NTOTAL,NDBL) :: SLOPE
      real, dimension(2*NL+2) :: ttsub


!     **************************************
!     * CALCULATE PTEMP AND SLOPE          *
!     **************************************

      !K  =  1
      !DO J  = 1, (2*NL+2)-1,2
      !    L  =  J
      !    TTsub(L) = tt(K)
      !    L  =  L+1
      !    TTsub(L) = t(K)
      !    K  =  K+1
      !END DO


      SBK=5.6704E-8
      SBKoverPI=SBK/PI

      DO 300 J            =   1,NDBL

          !IT1 = TTsub(J)*TTsub(J)*TTsub(J)*TTsub(J)*SBKoverPI

          if (MOD(J, 2) .eq. 0) THEN
              index_num = J / 2
              IT1 = T(index_num)*T(index_num)*T(index_num)*T(index_num)*SBKoverPI
          ELSE
              index_num = (J / 2) + 1
              IT1 = TT(index_num)*TT(index_num)*TT(index_num)*TT(index_num)*SBKoverPI
          END IF

          DO 200 L        = NSOLP+1,NTOTAL



              kindex          = max(1,j-1)
              PTEMP(L,J)=IT1
              SLOPE(L,J)   = (PTEMP(L,J)-PTEMP(L,KINDEX)) / TAUL(L,J)




              if( TAUL(L,J) .le. 1.0E-6 ) THEN
                  SLOPE(L,J) = 0.
              END IF

 200      CONTINUE
 300  CONTINUE



      RETURN
      END

