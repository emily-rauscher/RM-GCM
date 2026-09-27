      SUBROUTINE OPPR1_CORRK(TAUL, SLOPE, t,
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
     &  qrad,alb_tomi,alb_toai, num_layers, iband)
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
      use corrkmodule, only : PLANCK_INTS, NWNO
      include 'rcommons.h'
      
      integer num_layers, iband
      INTEGER, AUTOMATIC :: kindex, J, L, index_num

      INTEGER LLA, LLS, JDBLE, JDBLEDBLE, JN, JN2, iblackbody_above, ISL, IR, IRS
      REAL EMISIR, EPSILON, HEATI(NLAYER), HEATS(NLAYER), HEAT(NLAYER), SOLNET
      REAL TPI, SQ3, SBK,AM, AVG, ALOS
      REAL SCDAY, RGAS, GANGLE(3), GWEIGHT(3), GRATIO(3), EMIS(NBATCH), RSFX(NBATCH),NPROB(NBATCH), SOL(NBATCH)
      REAL RAYPERBAR(NBATCH),WEIGHT(NBATCH)
      REAL GOL(NBATCH,2*NL+2), WOL(NBATCH,2*NL+2), WAVE(5+1), TT(NL+1), Y3(NBATCH,3,2*NL+2), U0, FDEGDAY
      REAL WOT, GOT, PTEMPG(NBATCH), PTEMPT(NBATCH), G0(NBATCH,2*NL+2), OPD( NBATCH,2*NL+2), PTEMP(NBATCH,2*NL+2)
      REAL uG0(NBATCH,2*NL+2), uTAUL(NBATCH,2*NL+2), W0(NBATCH,2*NL+2), uW0(NBATCH,2*NL+2), uopd(NBATCH,2*NL+2),  U1S( NBATCH)
      REAL U1I(NBATCH), TOON_AK(NBATCH,2*NL+2), B1(NBATCH,2*NL+2), B2(NBATCH,2*NL+2), EE1( NBATCH,2*NL+2), EM1(NBATCH,2*NL+2)
      REAL EM2(NBATCH,2*NL+2), EL1( NBATCH,2*NL+2), EL2(NBATCH,2*NL+2), GAMI(NBATCH,2*NL+2), AF(NBATCH,4*NL+4)
      REAL BF(NBATCH,4*NL+4), EF(NBATCH,4*NL+4), SFCS(NBATCH), B3(NBATCH,2*NL+2), CK1(NBATCH,2*NL+2), CK2(NBATCH,2*NL+2)
      REAL CP(NBATCH,2*NL+2), CPB(NBATCH,2*NL+2), CM(NBATCH,2*NL+2), CMB(NBATCH,2*NL+2), DIRECT(NBATCH,2*NL+2), EE3(NBATCH,2*NL+2)
      REAL EL3(NBATCH,2*NL+2), FNET(NBATCH,2*NL+2), TMI(NBATCH,2*NL+2), AS(NBATCH,4*NL+4), DF(NBATCH,4*NL+4)
      REAL DS(NBATCH,4*NL+4), XK(NBATCH,4*NL+4), DIREC(NBATCH,2*NL+2), DIRECTU(NBATCH,2*NL+2), DINTENT(NBATCH,3,2*NL+2)
      REAL UINTENT(NBATCH,3,2*NL+2), TMID(NBATCH,2*NL+2), TMIU(NBATCH,2*NL+2), tslu,total_downwelling,alb_tot
      REAL tiru,firu(NKGAUSS),fird(NKGAUSS),fsLu(NKGAUSS), fsLd(NKGAUSS),fsLn(NKGAUSS),alb_toa(NKGAUSS), fupbs(NL+1)
      REAL fdownbs(NL+1),fnetbs(NL+1),fdownbs2(NL+1), fupbi(NL+1),fdownbi(NL+1),fnetbi(NL+1)
      REAL qrad(NL+1),alb_tomi,alb_toai

      real, automatic :: IT1
      real, DIMENSION(NLAYER) :: T
      real, dimension(NBATCH,2*NL+2) :: TAUL
      real, dimension(NBATCH,NDBL) :: SLOPE
      real, automatic :: localT
      INTEGER, AUTOMATIC :: temp_idx

      ! Thomas, make these data entries instead of regular vars.
    !   data PLANCK_C_1 /1.4724444e-50/ ! pre-factor for the planck function (2h/c^2)
    !   data PLANCK_C_2 /4.8014493e-11/ ! exponential factor for the planck function (h/kb)
    !   SBK=5.6704E-8
    !   SBKoverPI=SBK/PI
  
      ! PLANCK_TS is uniformly spaced 75..3999 K by 1 K (3925 entries).
      ! Index is computed directly instead of via a 3925-element linear scan.
      DO J = 1, NDBL
          IF (MOD(J, 2) .EQ. 0) THEN
              index_num = J / 2
              localT = T(index_num)
          ELSE
              index_num = (J / 2) + 1
              localT = TT(index_num)
          END IF
          kindex = MAX(1, J-1)
          IF (localT .GE. 75.) THEN
              temp_idx = INT(localT - 75.0) + 1
              temp_idx = MAX(1, MIN(temp_idx, 3924))
              IT1 = PLANCK_INTS(iband, temp_idx) +
     &              (PLANCK_INTS(iband, temp_idx+1) -
     &               PLANCK_INTS(iband, temp_idx)) *
     &              (localT - (74.0 + REAL(temp_idx)))
              DO L = NKGAUSS+1, NBATCH
                  PTEMP(L,J) = IT1
                  SLOPE(L,J) = (IT1 - PTEMP(L,kindex)) / TAUL(L,J)
                  IF (TAUL(L,J) .LE. 1.0E-6) SLOPE(L,J) = 0.
              END DO
          ELSE
              DO L = NKGAUSS+1, NBATCH
                  PTEMP(L,J) = 0.0
                  SLOPE(L,J) = -PTEMP(L,kindex) / TAUL(L,J)
                  IF (TAUL(L,J) .LE. 1.0E-6) SLOPE(L,J) = 0.
              END DO
          END IF
      END DO
    !   write(*,*) SUM(PTEMP(1:NTOTAL,1:NDBL), dim=1)/8



      RETURN
      END

