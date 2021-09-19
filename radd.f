      SUBROUTINE ADD

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
!     THIS SUBROUTINE FORMS THE MATRIX FOR THE MULTIPLE LAYERS AND
!     USES A TRIDIAGONAL ROUTINE TO FIND RADIATION IN THE ENTIRE
!     ATMOSPHERE.
!
      
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
               C1          =  B1(L,J) - DU0
               if( ABS(C1).lt.EPSILON ) THEN
               c1 = SIGN(EPSILON,C1)
               endif
               C2          =  AK(L,J)*AK(L,J) - DU0*DU0
               if( ABS(C2) .le. EPSILON ) then
               c2 = EPSILON
               endif
               CP1         =  W0(L,J)*(B3(L,J)*C1+B4*B2(L,J))/C2
               CPB(L,J)    =  CP1 * EL3(L,J)
               if( j .ne. 1 ) then
                 x4 = eL3(L,j1)
               else
                 x4 = soL(L)
               endif
               CP(L,J)     =  CP1 * X4
               CM1         =  ( CP1*B2(L,J) + W0(L,J)*B4 )/C1
               CMB(L,J)    =  CM1 * EL3(L,J)
               CM(L,J)     =  CM1 * X4

  10  CONTINUE

!
!       CALCULATE SFCS, THE SOURCE AT THE BOTTOM.
!
        DO 20 L            =  1,NSOLP
          SFCS(L)         =  DIRECT(L,NLAYER) * RSFX(L)
  20  CONTINUE
      END IF
!     ******************************
!     * CALCULATIONS FOR INFRARED. *
!     ******************************
!
      IF(IRS .NE. 0)  THEN
        DO 30 J           =   1,NLAYER
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
  40     SFCS(L)          = EMIS(L)*PTEMPG(L)*PI
      END IF
!
      J                =  0
      DO 42 JD         =  2,JN,2
         J             =  J + 1
         DO 42 L       =  LLS,NSOLP
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
      DO 44 L        = LLS,NSOLP
         DF(L,1)     = -CM(L,1)
         DF(L,JDBLE) = SFCS(L)+RSFX(L)*CMB(L,NLAYER)-CPB(L,NLAYER)
         DS(L,JDBLE) = DF(L,JDBLE)/BF(L,JDBLE)
  44     AS(L,JDBLE) = AF(L,JDBLE)/BF(L,JDBLE)
!IR
      DO 45 L        = NSOLP+1,LLA
         DF(L,1)     = -CM(L,1)
         DF(L,JDBLE) = SFCS(L)+RSFX(L)*CMB(L,NLAYER)-CPB(L,NLAYER)
         DS(L,JDBLE) = DF(L,JDBLE)/BF(L,JDBLE)
  45     AS(L,JDBLE) = AF(L,JDBLE)/BF(L,JDBLE)



!
!     (Where the magic happens...)
!
!     ********************************************
!     *     WE SOLVE THE TRIDIAGONAL EQUATIONS   *
!     ********************************************


      DO 46 J               = 2, JDBLE
         DO 46 L            = LLS,NSOLP
            X               = 1./(BF(L,JDBLE+1-J) -
     &                         EF(L,JDBLE+1-J)*AS(L,JDBLE+2-J))
            AS(L,JDBLE+1-J) = AF(L,JDBLE+1-J)*X
            DS(L,JDBLE+1-J) = (DF(L,JDBLE+1-J) - EF(L,JDBLE+1-J)
     &                         *DS(L,JDBLE+2-J))*X
  46  CONTINUE


!   NOW IR
      DO 47 J               = 2, JDBLE
         DO 47 L            = NSOLP+1,LLA
            X               = 1./(BF(L,JDBLE+1-J) -
     &                         EF(L,JDBLE+1-J)*AS(L,JDBLE+2-J))
            AS(L,JDBLE+1-J) = AF(L,JDBLE+1-J)*X
        DS(L,JDBLE+1-J) = (DF(L,JDBLE+1-J) - EF(L,JDBLE+1-J)
     &                         *DS(L,JDBLE+2-J))*X
  47  CONTINUE

      DO 48 L       = LLS,LLA
  48     XK(L,1)    = DS(L,1)
!
      DO 50 J       = 2, JDBLE
         DO 50 L    = LLS,LLA
            XK(L,J) = DS(L,J) - AS(L,J)*XK(L,J-1)
  50  CONTINUE

      DO 51 J       = 2, JDBLE
         DO 51 L    = NSOLP+1,LLA
            XK(L,J) = DS(L,J) - AS(L,J)*XK(L,J-1)
  51  CONTINUE


!
!  ***************************************************************
!     CALCULATE LAYER COEFFICIENTS, NET FLUX AND MEAN INTENSITY
!  ***************************************************************
!
      do J = 1,NLAYER
        do L = LLS,NSOLP
          CK1(L,J)   = XK(L,2*J-1)
          CK2(L,J)   = XK(L,2*J)

          FNET(L,J)  = CK1(L,J)  *( EL1(L,J) -EL2(L,J))   +
     &                 CK2(L,J) *( EM1(L,J)-EM2(L,J) ) + CPB(L,J) -
     &                  CMB(L,J) - DIRECT(L,J)
!
          TMI(L,J)   =  EL3(L,J) + U1I(L) *(CK1(L,J)  *
     &                  (EL1(L,J) + EL2(L,J))   +
     &                   CK2(L,J) *( EM1(L,J)+EM2(L,J) ) +
     &                   CPB(L,J) + CMB(L,J) )
        enddo
      enddo

!  AND AGAIN FOR IR

      do J = 1,NLAYER
        do L = NSOLP+1,NTOTAL

          CK1(L,J)   = XK(L,2*J-1)                                         
          CK2(L,J)   = XK(L,2*J)

          FNET(L,J)  = CK1(L,J) * (EL1(L,J) -EL2(L,J))   +
     &                 CK2(L,J) * (EM1(L,J)-EM2(L,J) ) + CPB(L,J) -
     &                  CMB(L,J) - DIRECT(L,J)
!                                                            
          TMI(L,J)   =  EL3(L,J) + U1I(L) *( CK1(L,J)  *
     &                  ( EL1(L,J) + EL2(L,J))   +
     &                   CK2(L,J) * (EM1(L,J)+EM2(L,J) ) +
     &                   CPB(L,J) + CMB(L,J) )                             
        enddo
      enddo



      RETURN
      END

