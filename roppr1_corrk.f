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
      use corrkmodule, only : PLANCK_INTS, PLANCK_TS, NWNO, NEAREST_INDEX
      include 'rcommons.h'
      
      integer kindex, J, L, num_layers, K, index_num

      INTEGER LLA, LLS, JDBLE, JDBLEDBLE, JN, JN2, iblackbody_above, ISL, IR, IRS
      REAL EMISIR, EPSILON, HEATI(NLAYER), HEATS(NLAYER), HEAT(NLAYER), SOLNET
      REAL TPI, SQ3, SBK,AM, AVG, ALOS
      REAL SCDAY, RGAS, GANGLE(3), GWEIGHT(3), GRATIO(3), EMIS(NTOTAL), RSFX(NTOTAL),NPROB(NTOTAL), SOL(NTOTAL)
      REAL RAYPERBAR(NTOTAL),WEIGHT(NTOTAL)
      REAL GOL(NTOTAL,2*NL+2), WOL(NTOTAL,2*NL+2), WAVE(5+1), TT(NL+1), Y3(NTOTAL,3,2*NL+2), U0, FDEGDAY
      REAL WOT, GOT, PTEMPG(NTOTAL), PTEMPT(NTOTAL), G0(NTOTAL,2*NL+2), OPD( NTOTAL,2*NL+2), PTEMP(NTOTAL,2*NL+2)
      REAL uG0(NTOTAL,2*NL+2), uTAUL(NTOTAL,2*NL+2), W0(NTOTAL,2*NL+2), uW0(NTOTAL,2*NL+2), uopd(NTOTAL,2*NL+2),  U1S( NTOTAL)
      REAL U1I(NTOTAL), TOON_AK(NTOTAL,2*NL+2), B1(NTOTAL,2*NL+2), B2(NTOTAL,2*NL+2), EE1( NTOTAL,2*NL+2), EM1(NTOTAL,2*NL+2)
      REAL EM2(NTOTAL,2*NL+2), EL1( NTOTAL,2*NL+2), EL2(NTOTAL,2*NL+2), GAMI(NTOTAL,2*NL+2), AF(NTOTAL,4*NL+4)
      REAL BF(NTOTAL,4*NL+4), EF(NTOTAL,4*NL+4), SFCS(NTOTAL), B3(NTOTAL,2*NL+2), CK1(NTOTAL,2*NL+2), CK2(NTOTAL,2*NL+2)
      REAL CP(NTOTAL,2*NL+2), CPB(NTOTAL,2*NL+2), CM(NTOTAL,2*NL+2), CMB(NTOTAL,2*NL+2), DIRECT(NTOTAL,2*NL+2), EE3(NTOTAL,2*NL+2)
      REAL EL3(NTOTAL,2*NL+2), FNET(NTOTAL,2*NL+2), TMI(NTOTAL,2*NL+2), AS(NTOTAL,4*NL+4), DF(NTOTAL,4*NL+4)
      REAL DS(NTOTAL,4*NL+4), XK(NTOTAL,4*NL+4), DIREC(NTOTAL,2*NL+2), DIRECTU(NTOTAL,2*NL+2), DINTENT(NTOTAL,3,2*NL+2)
      REAL UINTENT(NTOTAL,3,2*NL+2), TMID(NTOTAL,2*NL+2), TMIU(NTOTAL,2*NL+2), tslu,total_downwelling,alb_tot
      REAL tiru,firu(NIR),fird(NIR),fsLu(NSOL), fsLd(NSOL),fsLn(NSOL),alb_toa(NSOL), fupbs(NL+1)
      REAL fdownbs(NL+1),fnetbs(NL+1),fdownbs2(NL+1), fupbi(NL+1),fdownbi(NL+1),fnetbi(NL+1)
      REAL qrad(NL+1),alb_tomi,alb_toai

      real  ITP, ITG, IT1, SBKoverPI,g11
      real, DIMENSION(NLAYER) :: T
      real, dimension(NTOTAL,2*NL+2) :: TAUL
      real, dimension(NTOTAL,NDBL) :: SLOPE
      real, dimension(2*NL+2) :: ttsub
      real :: localT


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
      ! COMMON/PLANCK_INT/PLANCK_INTS, PLANCK_TS
      ! REAL :: PLANCK_INTS(NWNO, 3925), PLANCK_TS(3925)
      INTEGER :: temp_idx, gauss_idx, stel_idx
      logical :: lo_temp_flag

      ! Thomas, make these data entries instead of regular vars.
    !   data PLANCK_C_1 /1.4724444e-50/ ! pre-factor for the planck function (2h/c^2)
    !   data PLANCK_C_2 /4.8014493e-11/ ! exponential factor for the planck function (h/kb)
    !   SBK=5.6704E-8
    !   SBKoverPI=SBK/PI
  
      DO J            =   1,NDBL
          lo_temp_flag = .FALSE.

          !IT1 = TTsub(J)*TTsub(J)*TTsub(J)*TTsub(J)*SBKoverPI
          ! I think for corr-k this needs to be replaced with an integral of the Planck Function over the spectral bins
          ! I don't know how to rapidly integrate the planck function, there is definitely a faster method than, e.g., trapz
          ! Trying a look-up table for now
          if (MOD(J, 2) .eq. 0) THEN
              index_num = J / 2
              temp_idx = NEAREST_INDEX(PLANCK_TS, 3925, T(index_num))
              if (T(index_num) .LT. PLANCK_TS(temp_idx)) THEN
                temp_idx = temp_idx - 1
              END IF
              ! clamp so temp_idx and temp_idx+1 stay within PLANCK_TS/PLANCK_INTS bounds (1:3925)
              temp_idx = max(1, min(temp_idx, 3924))
              if (T(index_num) .GE. 75.) THEN
                  lo_temp_flag = .TRUE.
              END IF
              localT = T(index_num)
             
              ! IT1 = T(index_num)*T(index_num)*T(index_num)*T(index_num)*SBKoverPI
          ELSE
              index_num = (J / 2) + 1
              temp_idx = NEAREST_INDEX(PLANCK_TS, 3925, TT(index_num))
              if (TT(index_num) .LT. PLANCK_TS(temp_idx)) THEN
                temp_idx = temp_idx - 1
              END IF
              ! clamp so temp_idx and temp_idx+1 stay within PLANCK_TS/PLANCK_INTS bounds (1:3925)
              temp_idx = max(1, min(temp_idx, 3924))
              if (TT(index_num) .GE. 75.) THEN
                lo_temp_flag = .TRUE.
              END IF
              localT = TT(index_num)
            !   IT1 = TT(index_num)*TT(index_num)*TT(index_num)*TT(index_num)*SBKoverPI
          END IF
          IF (lo_temp_flag) THEN
              
            DO L        = NSOL+1,NTOTAL
                ! We've now figured out the temp index for the layer, so now we grab the planck integral
                ! Modulo is used here to figure out which bin we're in from the big array
                !   write(*,*) "L: ", L, "J: ", J, "temp_idx: ", temp_idx, "MODULO: ", MODULO(L-1,8)+1
                gauss_idx = MODULO(L-1,8)+1 ! This is the index of the gauss point
                ! stel_idx = (L - gauss_idx)/8 + 1 - NWNO! This is the index of the wavenumber bin
                stel_idx = (L-gauss_idx-NSOL)/8 + 1
                ! IT1 = PLANCK_INTS(stel_idx, temp_idx) ! Nearest-neighbor interpolation in T
                IT1 = PLANCK_INTS(stel_idx, temp_idx) + (PLANCK_INTS(stel_idx, temp_idx+1) - 
     &                PLANCK_INTS(stel_idx, temp_idx)) * 
     &                (localT - PLANCK_TS(temp_idx)) / (PLANCK_TS(temp_idx+1) - PLANCK_TS(temp_idx)) ! linear T interpolation
    !  &                (T(index_num) - PLANCK_TS(temp_idx)) / (PLANCK_TS(temp_idx+1) - PLANCK_TS(temp_idx)) ! linear T interpolation
                


                kindex     = max(1,j-1)
                PTEMP(L,J) = IT1
                SLOPE(L,J) = (PTEMP(L,J)-PTEMP(L,KINDEX)) / TAUL(L,J)


                if( TAUL(L,J) .le. 1.0E-6 ) THEN
                    SLOPE(L,J) = 0.
                END IF

            END DO
          ELSE
            DO L        = NSOL+1,NTOTAL
                IT1 = 0.0
                kindex = max(1,j-1)
                PTEMP(L,J)=IT1
                SLOPE(L,J)   = (PTEMP(L,J)-PTEMP(L,KINDEX)) / TAUL(L,J)

                if( TAUL(L,J) .le. 1.0E-6 ) THEN
                    SLOPE(L,J) = 0.
                END IF
            END DO
          END IF
        END DO
    !   write(*,*) SUM(PTEMP(1:NTOTAL,1:NDBL), dim=1)/8



      RETURN
      END

