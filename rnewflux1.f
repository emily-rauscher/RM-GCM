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
      REAL SCDAY, RGAS, GANGLE(3), GWEIGHT(3), GRATIO(3), EMIS(NBATCH), RSFX(NBATCH),NPROB(NBATCH), SOL(NBATCH)
      REAL RAYPERBAR(NBATCH),WEIGHT(NBATCH)
      REAL GOL(NBATCH,2*NL+2), WOL(NBATCH,2*NL+2), WAVE(5+1), TT(NL+1), Y3(NBATCH,3,2*NL+2), U0, FDEGDAY
      REAL WOT, GOT, PTEMPG(NBATCH), PTEMPT(NBATCH), G0(NBATCH,2*NL+2), OPD( NBATCH,2*NL+2), PTEMP(NBATCH,2*NL+2)
      REAL uG0(NBATCH,2*NL+2), uTAUL(NBATCH,2*NL+2), W0(NBATCH,2*NL+2), uW0(NBATCH,2*NL+2), uopd(NBATCH,2*NL+2),  U1S( NBATCH)
      REAL U1I(NBATCH), TOON_AK(NBATCH,2*NL+2), B1(NBATCH,2*NL+2), B2(  NBATCH,2*NL+2), EE1( NBATCH,2*NL+2), EM1(NBATCH,2*NL+2)
      REAL EM2(NBATCH,2*NL+2), EL1( NBATCH,2*NL+2), EL2(NBATCH,2*NL+2), GAMI(NBATCH,2*NL+2), AF(NBATCH,4*NL+4)
      REAL BF(NBATCH,4*NL+4), EF(NBATCH,4*NL+4), SFCS(NBATCH), B3(NBATCH,2*NL+2), CK1(NBATCH,2*NL+2), CK2(NBATCH,2*NL+2)
      REAL CP(NBATCH,2*NL+2), CPB(NBATCH,2*NL+2), CM(NBATCH,2*NL+2), CMB(NBATCH,2*NL+2), DIRECT(NBATCH,2*NL+2), EE3(NBATCH,2*NL+2)
      REAL EL3(NBATCH,2*NL+2), FNET(NBATCH,2*NL+2), TMI(NBATCH,2*NL+2), AS(NBATCH,4*NL+4), DF(NBATCH,4*NL+4)
      REAL DS(NBATCH,4*NL+4), XK(NBATCH,4*NL+4), DIREC(NBATCH,2*NL+2), DIRECTU(NBATCH,2*NL+2), DINTENT(NBATCH,3,2*NL+2)
      REAL UINTENT(NBATCH,3,2*NL+2), TMID(NBATCH,2*NL+2), TMIU(NBATCH,2*NL+2), tslu,total_downwelling,alb_tot
      REAL tiru,firu(NKGAUSS),fird(NKGAUSS),fsLu(NKGAUSS), fsLd(NKGAUSS),fsLn(NKGAUSS),alb_toa(NKGAUSS), fupbs(NL+1)
      REAL fdownbs(NL+1),fnetbs(NL+1),fdownbs2(NL+1), fupbi(NL+1),fdownbi(NL+1),fnetbi(NL+1)
      REAL qrad(NL+1),alb_tomi,alb_toai

      real, dimension(NBATCH,2*NL+2) :: TAUL
      real, dimension(NBATCH,NDBL) :: SLOPE
!
!     LOCAL DIMENSIONS
      REAL, DIMENSION(NBATCH,NGAUSS,NDBL) :: Y1, Y2, Y4, Y8
      REAL, DIMENSION(NBATCH,NDBL)        :: A1, A2, A3, A4, A5, A7, Y5

      A3(:,:) = 0.0
      A7(:,:) = 0.0
      DO 200 J           =  1,NDBL
          kindex         = max( 1, j-1 )
          DO 100  L      =  NKGAUSS+1,NBATCH
!            HERE WE DO NO SCATTERING COEFFICIENTS
             A3(L,J)     =  PTEMP(L,KINDEX)*TPI
             A4(L,J)     =  TPI*SLOPE(L,J)
             A7(L,J)     =  A3(L,J)
             Y5(L,J)     =  A4(L,J)*TAUL(L,J)
 100      CONTINUE


!         HERE WE DO SCATTERING
          IF(IRS .NE. 0) THEN
              DO 50 L    =  NKGAUSS+1,NBATCH
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
            DO 300 L         =  NKGAUSS+1,NBATCH
               Y1(L,I,J)  =  0.0
               Y2(L,I,J)  =  0.0
               Y4(L,I,J)  =  A7(L,J) - A4(L,J)*GANGLE(I)
               Y8(L,I,J)  =  A3(L,J)+A4(L,J)*GANGLE(I)
 300        CONTINUE
!
!           HERE WE DO SCATTERING
            IF(IRS .NE. 0) THEN
              DO 325 L    =  NKGAUSS+1,NBATCH
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
         DO 425  L         =  NKGAUSS+1,NBATCH
            TMID(L,J) = 0.0
            TMIU(L,J) = 0.0
            DIREC(L,J)     =  0.0
            DIRECTU(L,J)   =  0.0
 425     CONTINUE
 450  CONTINUE


!     DIREC IS DOWNWARD FLUX. DIRECTU IS UPWARD FLUX.
!     CALCULATE DINTENT THE DOWNWARD INTENSITY AND DIREC THE DOWNWARD FL

       DO 500 I             = 1,NGAUSS
          DO 475 L          = NKGAUSS+1,NBATCH
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
              DO 510 L           = NKGAUSS+1,NBATCH
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
          DO 560  L               =  NKGAUSS+1,NBATCH
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
             DO 630 L              = NKGAUSS+1,NBATCH
                  UINTENT(L,I,J)    = (UINTENT(L,I,J+1)-Y5(L,J+1))*Y3(L,I,J+1)+Y2(L,I,J+1)+(1.-Y3(L,I,J+1))*Y8(L,I,J+1)
                  TMIU(L,J)        = TMIU(L,J)+UINTENT(L,I,J)*GRATIO(I)
                  DIRECTU(L,J)     = DIRECTU(L,J) + GWEIGHT(I)*UINTENT(L,I,J)
 630         CONTINUE
 640      CONTINUE
 650  CONTINUE





      RETURN
      END