      SUBROUTINE NEWFLUX1(TAUL, SLOPE,
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
     &  qrad,alb_tomi,alb_toai, num_layers, Y1, Y2, Y4, Y8, A1, A2, A3, A4, A5, A7, Y5)
!
!     **************************************************************
!     *  Purpose             :  Calculate upward and downward      *
!     *                         intensities and fluxes using Gauss *
!     *                         Quadrature angles and weights.     *
!     *  Subroutines Called  :  None                               *
!     *  Input               :  PTEMP, SLOPE, Y3, B3, EE1, EE2     *
!     *  Output              :  DINTENT, UINTENT, DIREC, DIRECTU   *
!     * ************************************************************
!
      include 'rcommons.h'

      INTEGER LLA, LLS, JDBLE, JDBLEDBLE, JN, JN2, iblackbody_above, ISL, IR, IRS, M, I, L, kindex, J
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
      REAL qrad(NL+1),alb_tomi,alb_toai

      real, dimension(5,2*NL+2) :: TAUL
      real, dimension(NTOTAL,NDBL) :: SLOPE
!
!     LOCAL DIMENSIONS
      REAL, DIMENSION(NTOTAL,NGAUSS,NDBL) :: Y1, Y2, Y4, Y8
      REAL, DIMENSION(NTOTAL,NDBL)        :: A1, A2, A3, A4, A5, A7, Y5

      A3(:,:) = 0.0
      A7(:,:) = 0.0
      DO 200 J           =  1,NDBL
          kindex         = max( 1, j-1 )
          DO 100  L      =  NSOLP+1,NTOTAL
!            HERE WE DO NO SCATTERING COEFFICIENTS
             A3(L,J)     =  PTEMP(L,KINDEX)*TPI
             A4(L,J)     =  TPI*SLOPE(L,J)
             A7(L,J)     =  A3(L,J)
             Y5(L,J)     =  A4(L,J)*TAUL(L,J)
 100      CONTINUE


!         HERE WE DO SCATTERING
          IF(IRS .NE. 0) THEN
              DO 50 L    =  NSOLP+1,NTOTAL
                A1(L,J)  =  U1I(L) - TOON_AK(L,J)
                A2(L,J)  =  GAMI(L,J)*(TOON_AK(L,J)+U1I(L))
                A3(L,J)  =  A3(L,J)+(SLOPE(L,J)*(TPI*B3(L,J)-U1S(L)))
                A7(L,J)  =  A7(L,J)-(SLOPE(L,J)*(TPI*B3(L,J)-U1S(L)))
 50           CONTINUE
          ENDIF
  200 CONTINUE


!     CALCULATIONS FOR ALL GAUSS POINTS. HERE WE DO NO SCATTERING COEFFI
!
      DO 400       J         =  1,NDBL
         DO 350    I         =  1,NGAUSS
            DO 300 L         =  NSOLP+1,NTOTAL
               Y1(L,I,J)  =  0.0
               Y2(L,I,J)  =  0.0
               Y4(L,I,J)  =  A7(L,J) - A4(L,J)*GANGLE(I)
               Y8(L,I,J)  =  A3(L,J)+A4(L,J)*GANGLE(I)
 300        CONTINUE
!
!           HERE WE DO SCATTERING
            IF(IRS .NE. 0) THEN
              DO 325 L    =  NSOLP+1,NTOTAL
                 YA        =  A1(L,J)*(Y3(L,I,J)-EE1(L,J))/
     &                             (TOON_AK(L,J)*GANGLE(I)-1.)
                 YB        =  A2(L,J)*(1.- EE1(L,J)*Y3(L,I,J))/
     &                             (TOON_AK(L,J)*GANGLE(I)+1.)
                 CKP= CK1(L,J)+CK2(L,J)
                 CKM= CK1(L,J) -CK2(L,J)
                 Y1(L,I,J) =  CKP*YB+CKM*YA
                 Y2(L,I,J) =  CKP*YA+CKM*YB
 325          CONTINUE
            ENDIF
 350     CONTINUE
 400  CONTINUE




!
      DO 450 J             =  1,NDBL
         DO 425  L         =  NSOLP+1,NTOTAL
            TMID(L,J) = 0.0
            TMIU(L,J) = 0.0
            DIREC(L,J)     =  0.0
            DIRECTU(L,J)   =  0.0
 425     CONTINUE
 450  CONTINUE


!     DIREC IS DOWNWARD FLUX. DIRECTU IS UPWARD FLUX.
!     CALCULATE DINTENT THE DOWNWARD INTENSITY AND DIREC THE DOWNWARD FL

       DO 500 I             = 1,NGAUSS
          DO 475 L          = NSOLP+1,NTOTAL
             if( iblackbody_above .eq. 1 )then
               DINTENT(L,I,1) = PTEMPT(L)*Y3(L,I,1)*TPI +Y1(L,I,1)+(1.-Y3(L,I,1))*Y4(L,I,1)
             else
               DINTENT(L,I,1) = (1.-Y3(L,I,1))*Y4(L,I,1)+Y1(L,I,1)
             endif

             TMID(L,1)      = TMID(L,1)+DINTENT(L,I,1)*GRATIO(I)
             DIREC(L,1)     = DIREC(L,1)+DINTENT(L,I,1)*GWEIGHT(I)
 475      CONTINUE
 500   CONTINUE


!      DINTENT IS DOWNWARD INTENSITY * TPI. DIREC IS THE DOWNWARD FLUX.
       DO 530        J           = 2,NDBL
           DO 520    I           = 1,NGAUSS
              DO 510 L           = NSOLP+1,NTOTAL
                 DINTENT(L,I,J)  = DINTENT(L,I,J-1)*Y3(L,I,J)
     &                              +Y1(L,I,J)+Y5(L,J)+
     &                              (1.-Y3(L,I,J))*Y4(L,I,J)
                 TMID(L,J)       = TMID(L,J)+DINTENT(L,I,J)*GRATIO(I)
                 DIREC(L,J)      = DIREC(L,J)+DINTENT(L,I,J)*
     &                              GWEIGHT(I)
 510          CONTINUE
 520       CONTINUE
 530   CONTINUE

!
!     UINTENT IS THE UPWARD INTENSITY * TPI. DIRECTU IS THE UPWARD FLUX.
!     ASSUME THAT THE REFLECTIVITY IS LAMBERT.


       DO 570     I               =  1,NGAUSS
          DO 560  L               =  NSOLP+1,NTOTAL
             UINTENT(L,I,NDBL)  =  PTEMPG(L)*EMIS(L)
     &                               *TPI+2.*RSFX(L)*DIREC(L,NDBL)
             TMIU(L,NDBL)       =  TMIU(L,NDBL)+
     &                               UINTENT(L,I,NDBL)*GRATIO(I)
             DIRECTU(L,NDBL)    =  DIRECTU(L,NDBL)+
     &                               UINTENT(L,I,NDBL)*GWEIGHT(I)
 560      CONTINUE
 570   CONTINUE


!
      DO 650        M              = 2,NDBL
          J                        = NDBL-M+1
          DO 640    I              = 1,NGAUSS
             DO 630 L              = NSOLP+1,NTOTAL
                  UINTENT(L,I,J)    = (UINTENT(L,I,J+1)-Y5(L,J+1))*Y3(L,I,J+1)+Y2(L,I,J+1)+(1.-Y3(L,I,J+1))*Y8(L,I,J+1)
                  TMIU(L,J)        = TMIU(L,J)+UINTENT(L,I,J)*GRATIO(I)
                  DIRECTU(L,J)     = DIRECTU(L,J) + GWEIGHT(I)*UINTENT(L,I,J)
 630         CONTINUE
 640      CONTINUE
 650  CONTINUE





      RETURN
      END