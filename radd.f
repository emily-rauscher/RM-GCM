      SUBROUTINE ADD(TAUL,solar_calculation_indexer, SLOPE,
     &               LLA, LLS, JDBLE, JDBLEDBLE, JN, JN2, iblackbody_above, ISL, IR, IRS, EMISIR,
     &               EPSILON, HEATI, HEATS, HEAT, SOLNET,TPI, SQ3, SBK,AM, AVG, ALOS,
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
!     ***************************************************************
!     *  Purpose             :  Defines source terms, form matrix   *
!     *                         for multiple layers and solve tri-  *
!     *                         diagnol equations to obtain mean    *
!     *                         intensity and net flux.             *
!     *  Subroutines Called  :  None                                *
!     *  Input               :  G0, U0, RSFX, TAUL, OPD             *
!     *  Output              :  B3, EE3, DIRECT, SFCS, CPB, CMB, DS *
!     * *************************************************************
!
      include 'rcommons.h'


      INTEGER LLA, LLS, JDBLE, JDBLEDBLE, JN, JN2, iblackbody_above, ISL, IR, IRS, kindex, J, K, L
      REAL EMISIR, EPSILON, HEATI(NLAYER), HEATS(NLAYER), HEAT(NLAYER), SOLNET
      REAL TPI, SQ3, SBK,AM, AVG, ALOS
      REAL SCDAY, RGAS, GANGLE(3), GWEIGHT(3), GRATIO(3), EMIS(5), RSFX(5),NPROB(5), SOL(5),RAYPERBAR(5),WEIGHT(5)
      REAL GOL(5,2*NL+2), WOL(5,2*NL+2), WAVE(5+1), TT(NL+1), Y3(5,3,2*NL+2), U0, FDEGDAY
      REAL WOT, GOT, PTEMPG(5), PTEMPT(5), G0(5,2*NL+2), OPD( 5,2*NL+2), PTEMP(5,2*NL+2)
      REAL uG0(5,2*NL+2), uTAUL(5,2*NL+2), W0(5,2*NL+2), uW0(5,2*NL+2), uopd(5,2*NL+2),  U1S( 5)
      REAL U1I(5), TOON_AK(5,2*NL+2), B1(5,2*NL+2), B2(  5,2*NL+2), EE1( 5,2*NL+2), EM1(5,2*NL+2)
      REAL EM2(5,2*NL+2), EL1( 5,2*NL+2), EL2(5,2*NL+2), GAMI(5,2*NL+2), AF(5,4*NL+4)
      REAL BF(5,4*NL+4), EF(5,4*NL+4), SFCS(5), B3(5,2*NL+2), CK1(5,2*NL+2), CK2(5,2*NL+2)
      REAL CP(5,2*NL+2), CPB(5,2*NL+2), CM(5,2*NL+2), CMB(5,2*NL+2), DIRECT(5,2*NL+2), EE3(5,2*NL+2)
      REAL EL3(5,2*NL+2), FNET(5,2*NL+2), TMI(5,2*NL+2), AS(5,4*NL+4), DF(5,4*NL+4)
      REAL DS(5,4*NL+4), XK(5,4*NL+4), DIREC(5,2*NL+2), DIRECTU(5,2*NL+2), DINTENT(5,3,2*NL+2)
      REAL UINTENT(5,3,2*NL+2), TMID(5,2*NL+2), TMIU(5,2*NL+2), tslu,total_downwelling,alb_tot
      REAL tiru,firu(2),fird(2),fsLu(3), fsLd(3),fsLn(3),alb_toa(3), fupbs(NL+1)
      REAL fdownbs(NL+1),fnetbs(NL+1),fdownbs2(NL+1), fupbi(NL+1),fdownbi(NL+1),fnetbi(NL+1)
      REAL qrad(NL+1),alb_tomi,alb_toai, x4_add


      real C2_VAR, C1_VAR


!     THIS SUBROUTINE FORMS THE MATRIX FOR THE MULTIPLE LAYERS AND
!     USES A TRIDIAGONAL ROUTINE TO FIND RADIATION IN THE ENTIRE
!     ATMOSPHERE.

      real, dimension(5,2*NL+2) :: TAUL
      integer solar_calculation_indexer
      real, dimension(NTOTAL,NDBL) :: SLOPE

!     ******************************
!     *   CALCULATIONS FOR SOLAR   *
!     ******************************

      IF(ISL .NE. 0)  THEN
        DU0                =  1./U0
        DO 10 J            =  1,NLAYER
            j1 = max( 1, j-1 )
            DO 10 L        =  1,NSOLP
               B3(L,J)     =  0.5*(1.-SQ3*G0(L,J)*U0)
               B4          =  1. - B3(L,J)
               X2          =  TAUL(L,J)*DU0
               EE3(L,J)    =  EXP(-X2)
               X3          =  OPD(L,J)*DU0
               EL3(L,J)    =  EXP(-X3)*SOL(L)

               DIRECT(L,J) = U0*EL3(L,J)
               C1_VAR          =  B1(L,J) - DU0

               if( ABS(C1_VAR).lt. 1e-6 ) THEN
               C1_VAR = SIGN(1e-6,C1_VAR)
               endif

               C2_VAR          =  TOON_AK(L,J)*TOON_AK(L,J) - DU0*DU0

               if( ABS(C2_VAR) .le. 1e-6 ) then
               C2_VAR = 1e-6
               endif

               CP1         =  W0(L,J)*(B3(L,J)*C1_VAR+B4*B2(L,J))/C2_VAR
               CPB(L,J)    =  CP1 * EL3(L,J)

               if( j .ne. 1 ) then
                 x4_add = eL3(L,j1)
               else
                 x4_add = soL(L)
               endif

               CP(L,J)     =  CP1 * x4_add
               CM1         =  ( CP1*B2(L,J) + W0(L,J)*B4 )/C1_VAR
               CMB(L,J)    =  CM1 * EL3(L,J)
               CM(L,J)     =  CM1 * x4_add

  10  CONTINUE
        DO 20 L            =  1,NSOLP
          SFCS(L)         =  DIRECT(L,NLAYER) * RSFX(L)
  20  CONTINUE
      END IF


!     ******************************
!     * CALCULATIONS FOR INFRARED. *
!     ******************************
!
      IF(IRS .NE. 0)  THEN
        DO 30 J           =   1,NDBL
           KINDEX         = max(1,J-1)
           DO 30 L        = NSOLP+1,NTOTAL
              B3(L,J)     = 1.0/(B1(L,J)+B2(L,J))
              CP(L,J)     = (PTEMP(L,KINDEX)+SLOPE(L,J)*B3(L,J))*U1S(L)
              CPB(L,J)    = CP(L,J) + SLOPE(L,J)*TAUL(L,J)*U1S(L)
              CM(L,J)     = (PTEMP(L,KINDEX)-SLOPE(L,J)*B3(L,J))*U1S(L)
              CMB(L,J)    = CM(L,J) + SLOPE(L,J)*TAUL(L,J)*U1S(L)
              EL3(L,J)    = 0.0
              DIRECT(L,J) = 0.0
              EE3(L,J)    = 0.0
  30  CONTINUE


      DO 40 L             = NSOLP+1,NTOTAL
  40     SFCS(L)          = EMIS(L)*PTEMPG(L)*3.141592653589
      END IF
!
      J                =  0
      DO 42 JD         =  2,JN,2
         J             =  J + 1
         DO 42 L       =  solar_calculation_indexer,NSOLP
!           HERE ARE THE EVEN MATRIX ELEMENTS
            DF(L,JD) = (CP(L,J+1) - CPB(L,J))*EM1(L,J+1) -
     &                  (CM(L,J+1) - CMB(L,J))*EM2(L,J+1)
!           HERE ARE THE ODD MATRIX ELEMENTS EXCEPT FOR THE TOP.
            DF(L,JD+1) =  EL2(L,J) * (CP(L,J+1)-CPB(L,J)) +
     &                    EL1(L,J) * (CMB(L,J) - CM(L,J+1))
  42  CONTINUE


!     AGAIN FOR THE IR
      J                =  0
      DO 43 JD         =  2,JN2,2
         J             =  J + 1
         DO 43 L       =  NSOLP+1,LLA

!           HERE ARE THE EVEN MATRIX ELEMENTS
            DF(L,JD) = (CP(L,J+1) - CPB(L,J))*EM1(L,J+1) -
     &                  (CM(L,J+1) - CMB(L,J))*EM2(L,J+1)

!           HERE ARE THE ODD MATRIX ELEMENTS EXCEPT FOR THE TOP.
            DF(L,JD+1) =  EL2(L,J) * (CP(L,J+1)-CPB(L,J)) +
     &                    EL1(L,J) * (CMB(L,J) - CM(L,J+1))

   43  CONTINUE
!     HERE ARE THE TOP AND BOTTOM BOUNDARY CONDITIONS AS WELL AS THE
!     BEGINNING OF THE TRIDIAGONAL SOLUTION DEFINITIONS. I ASSUME NO
!     DIFFUSE RADIATION IS INCIDENT AT THE TOP.
!
!VIS
      DO 44 L        = solar_calculation_indexer,NSOLP
         DF(L,1)     = -CM(L,1)
         DF(L,JDBLE) = SFCS(L)+RSFX(L)*CMB(L,NLAYER)-CPB(L,NLAYER)
         DS(L,JDBLE) = DF(L,JDBLE)/BF(L,JDBLE)
  44     AS(L,JDBLE) = AF(L,JDBLE)/BF(L,JDBLE)
!IR
      DO 45 L        = NSOLP+1,LLA
         DF(L,1)     = -CM(L,1)
         DF(L,JDBLEDBLE) = SFCS(L)+RSFX(L)*CMB(L,NDBL)-CPB(L,NDBL)
         DS(L,JDBLEDBLE) = DF(L,JDBLEDBLE)/BF(L,JDBLEDBLE)
  45     AS(L,JDBLEDBLE) = AF(L,JDBLEDBLE)/BF(L,JDBLEDBLE)

!
!     (Where the magic happens...)
!
!     ********************************************
!     *     WE SOLVE THE TRIDIAGONAL EQUATIONS   *
!     ********************************************

      DO 46 J               = 2, JDBLE
         DO 46 L            = solar_calculation_indexer,NSOLP
            X               = 1./(BF(L,JDBLE+1-J) - EF(L,JDBLE+1-J)*AS(L,JDBLE+2-J))
            AS(L,JDBLE+1-J) = AF(L,JDBLE+1-J)*X
            DS(L,JDBLE+1-J) = (DF(L,JDBLE+1-J) - EF(L,JDBLE+1-J) *DS(L,JDBLE+2-J))*X
  46  CONTINUE


!   NOW IR
      DO 47 J               = 2, JDBLEDBLE
         DO 47 L            = NSOLP+1,LLA
            X               = 1./(BF(L,JDBLEDBLE+1-J) - EF(L,JDBLEDBLE+1-J)*AS(L,JDBLEDBLE+2-J))
            AS(L,JDBLEDBLE+1-J) = AF(L,JDBLEDBLE+1-J)*X
            DS(L,JDBLEDBLE+1-J) = (DF(L,JDBLEDBLE+1-J) - EF(L,JDBLEDBLE+1-J)*DS(L,JDBLEDBLE+2-J))*X
  47  CONTINUE


      DO 48 L       = solar_calculation_indexer,NTOTAL
  48     XK(L,1)    = DS(L,1)

      DO 50 J       = 2, JDBLE
         DO 50 L    = solar_calculation_indexer,NTOTAL
            XK(L,J) = DS(L,J) - AS(L,J)*XK(L,J-1)
  50  CONTINUE

      DO 51 J       = 2, JDBLEDBLE
         DO 51 L    = NSOLP+1,NTOTAL
            XK(L,J) = DS(L,J) - AS(L,J)*XK(L,J-1)
  51  CONTINUE


!  ***************************************************************
!     CALCULATE LAYER COEFFICIENTS, NET FLUX AND MEAN INTENSITY
!  ***************************************************************

      do J = 1,NLAYER
        do L = solar_calculation_indexer,NSOLP
          CK1(L,J)   = XK(L,2*J-1)
          CK2(L,J)   = XK(L,2*J)

          FNET(L,J)  = CK1(L,J)  *( EL1(L,J) -EL2(L,J))   +
     &                 CK2(L,J) *( EM1(L,J)-EM2(L,J) ) + CPB(L,J) -
     &                  CMB(L,J) - DIRECT(L,J)
!
          TMI(L,J)   =  EL3(L,J) + U1I(L) *( CK1(L,J)  *
     &                  ( EL1(L,J) + EL2(L,J))   +
     &                   CK2(L,J) *( EM1(L,J)+EM2(L,J) ) +
     &                   CPB(L,J) + CMB(L,J) )
        enddo
      enddo

!  AND AGAIN FOR IR

      do J = 1,NDBL
        do L = NSOLP+1,NTOTAL
          CK1(L,J)   = XK(L,2*J-1)
          CK2(L,J)   = XK(L,2*J)

          FNET(L,J)  = CK1(L,J) * (EL1(L,J) -EL2(L,J)) +
     &                 CK2(L,J) * (EM1(L,J)-EM2(L,J) ) + CPB(L,J) -
     &                 CMB(L,J) - DIRECT(L,J)

          TMI(L,J)   =  EL3(L,J) + U1I(L) * (CK1(L,J) *
     &                  (EL1(L,J) + EL2(L,J)) +
     &                  CK2(L,J) * (EM1(L,J) + EM2(L,J)) +
     &                  CPB(L,J) + CMB(L,J))

        enddo
      enddo

      RETURN
      END
