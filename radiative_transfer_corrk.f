      subroutine opacity_wrapper_corrk(t, p_pass, tau_IRe, tau_Ve,
     &                           Beta_V, Beta_IR, gravity_SI, incident_starlight_fraction,
     &             LLA, LLS, JDBLE, JDBLEDBLE, JN, JN2, iblackbody_above, ISL, IR, IRS,EMISIR,
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
     &  qrad,alb_tomi,alb_toai, num_layers,
     &  dpe, Pl, Tl, pe,
     &  k_IR, k_lowP, k_hiP, Tin, Pin, Freedman_met,
     &  Freedman_T, Freedman_P, Tl10, Pl10, temperature_val, pressure_val, k_IRl, k_Vl, tau_ray_temp)
          use corrkmodule, only : corrk_setup, TS_CORRK, PS_CORRK, TS_LOG_CORRK, 
     &           PS_LOG_CORRK, WGTS_CORRK, WNO_EDGES, WNO_CTRS, STEL_SPEC, INT_SPEC, TAURAY_PER_DPG,
     &           OPAC_CORRK, PLANCK_INTS, PLANCK_TS
          include 'rcommons.h'
          include 'nwno.inc'
          ! integer :: NWNO

          ! Thomas inclusions from corr-k common block
    !       COMMON/CORRKGAS/OPAC_CORRK, TS_CORRK, PS_CORRK, TS_LOG_CORRK, PS_LOG_CORRK, WGTS_CORRK, WNO_EDGES, WNO_CTRS, STEL_SPEC,
    !  &                    INT_SPEC, TAURAY_PER_DPG
    !       REAL :: OPAC_CORRK(NTGRID, 20, NWNO, 8)
    !       REAL :: TS_CORRK(NTGRID), PS_CORRK(20), TS_LOG_CORRK(NTGRID), PS_LOG_CORRK(20), WGTS_CORRK(8) 
    !       REAL :: WNO_EDGES(NWNO+1), WNO_CTRS(NWNO), STEL_SPEC(NWNO), INT_SPEC(NWNO)
    !       REAL :: TAURAY_PER_DPG(NTGRID, 20, NWNO)
          ! End Thomas inclusions
          REAL :: GASCON

          INTEGER LLA, LLS, JDBLE, JDBLEDBLE, JN, JN2, iblackbody_above, ISL, IR, IRS
          REAL EMISIR, EPSILON, HEATI(NL+1), HEATS(NL+1), HEAT(NL+1), SOLNET
          REAL TPI, SQ3, SBK,AM, AVG, ALOS
          REAL SCDAY, RGAS, GANGLE(3), GWEIGHT(3), GRATIO(3), EMIS(NTOTAL), RSFX(NTOTAL),NPROB(NTOTAL), SOL(NTOTAL)
          REAL RAYPERBAR(NTOTAL),WEIGHT(NTOTAL)
          REAL GOL(NTOTAL,2*NL+2), WOL(NTOTAL,2*NL+2), WAVE(NTOTAL+1), TT(NL+1), Y3(NTOTAL,3,2*NL+2), U0, FDEGDAY
          REAL WOT, GOT, PTEMPG(NTOTAL), PTEMPT(NTOTAL), G0(NTOTAL,2*NL+2), OPD(NTOTAL,2*NL+2), PTEMP(NTOTAL,2*NL+2)
          REAL uG0(NTOTAL,2*NL+2), uTAUL(NTOTAL,2*NL+2), W0(NTOTAL,2*NL+2), uW0(NTOTAL,2*NL+2), uopd(NTOTAL,2*NL+2),  U1S(NTOTAL)
          REAL U1I(NTOTAL), TOON_AK(NTOTAL,2*NL+2), B1(NTOTAL,2*NL+2), B2(  5,2*NL+2), EE1(NTOTAL,2*NL+2), EM1(NTOTAL,2*NL+2)
          REAL EM2(NTOTAL,2*NL+2), EL1(NTOTAL,2*NL+2), EL2(NTOTAL,2*NL+2), GAMI(NTOTAL,2*NL+2), AF(NTOTAL,4*NL+4)
          REAL BF(NTOTAL,4*NL+4), EF(NTOTAL,4*NL+4), SFCS(NTOTAL), B3(NTOTAL,2*NL+2), CK1(NTOTAL,2*NL+2), CK2(NTOTAL,2*NL+2)
          REAL CP(NTOTAL,2*NL+2), CPB(NTOTAL,2*NL+2), CM(NTOTAL,2*NL+2), CMB(NTOTAL,2*NL+2)
          REAL DIRECT(NTOTAL,2*NL+2), EE3(NTOTAL,2*NL+2)
          REAL EL3(NTOTAL,2*NL+2), FNET(NTOTAL,2*NL+2), TMI(NTOTAL,2*NL+2), AS(NTOTAL,4*NL+4), DF(NTOTAL,4*NL+4)
          REAL DS(NTOTAL,4*NL+4), XK(NTOTAL,4*NL+4), DIREC(NTOTAL,2*NL+2), DIRECTU(NTOTAL,2*NL+2), DINTENT(NTOTAL,3,2*NL+2)
          REAL UINTENT(NTOTAL,3,2*NL+2), TMID(NTOTAL,2*NL+2), TMIU(NTOTAL,2*NL+2), tslu,total_downwelling,alb_tot
          REAL tiru,firu(NIR),fird(NIR),fsLu(NSOL), fsLd(NSOL),fsLn(NSOL),alb_toa(NSOL), fupbs(NL+1)
          REAL fdownbs(NL+1),fnetbs(NL+1),fdownbs2(NL+1), fupbi(NL+1),fdownbi(NL+1),fnetbi(NL+1)
          REAL qrad(NL+1),alb_tomi,alb_toai

          real, dimension(NIR, NL+1) :: k_IRl
          real :: tau_ray_temp(NSOL, NL+1)
          real, dimension(NSOL, NL+1) :: k_Vl

          integer :: NLAYER, J, k
          real :: Tirr, Tint, gravity_SI, incident_starlight_fraction

          real, dimension(NIR) :: Beta_IR
          real, dimension(NSOL) :: Beta_V

          real, dimension(NIR,NL+2) :: tau_IRe
          real, dimension(NSOL,NL+2) :: tau_Ve

          real, dimension(NL+1) :: dpe, Pl, Tl, pe, p_pass, t
          real :: k_IR, k_lowP, k_hiP, Tin, Pin, Freedman_met
          real :: Freedman_T, Freedman_P, Tl10, Pl10, temperature_val, pressure_val

          ! dummy vars
          real :: total_layer_taus(NLAYER), dummy_weights(NSOL)
          integer :: gauss_idx, stel_idx

          ! WRITE(*,*) "OPAC_CORRK: ", OPAC_CORRK(1,1,1,1), OPAC_CORRK(1,1,2,1), OPAC_CORRK(1,1,3,1)

          ! write(*,*) 'gravity_SI: ', gravity_SI
          do J = 1, NL+1
             pe(J) = p_pass(J)
          end do

          DO J = 1, NL
              dpe(J) = pe(J+1) - pe(J)
              pl(J) = dpe(J) / log(pe(J+1)/pe(J))
          END DO

          if (tt(1) .ge. 1.0) then
              DO J = 1, NLAYER
                  Tl(J) = tt(J)
              END DO
          else
              write(*,*) 'Something went wrong with the temperature profile'
          end if

          dpe(NLAYER) = 10.0 ** (LOG10(dpe(NLAYER-1)) + (LOG10(dpe(NLAYER-1)) - LOG10(dpe(NLAYER-2))))
          pl(NLAYER)  = 10.0 ** (LOG10(pl(NLAYER-1))  + (LOG10(pl(NLAYER-1))  - LOG10(pl(NLAYER-2))))
          Tl(NLAYER)  = Tl(NLAYER-1) + ABS(Tl(NLAYER-1) - Tl(NLAYER-2)) / 2.0
          CALL calculate_opacities_corrk(NLAYER, NSOL, NIR, Tl, Pl, dpe, tau_IRe, tau_Ve, gravity_SI, k_IRl, k_Vl,
     &                                   OPAC_CORRK, TS_CORRK, PS_CORRK, TS_LOG_CORRK, PS_LOG_CORRK,
     &                                   tau_ray_temp, TAURAY_PER_DPG,NWNO,NTGRID,NPGRID)
          
          ! DO L = 1, NSOL
          !   ! write(*,*)  'L: ', L
          !   gauss_idx = MODULO(L-1,8)+1
          !   stel_idx = (L - gauss_idx) / 8 + 1
          !   ! write(*,*) 'gauss_idx: ', gauss_idx
          !   ! stel_idx = (L - chan_idx) / 8 + 1
          !   dummy_weights(L) = WGTS_CORRK(gauss_idx) * STEL_SPEC(stel_idx)
          !   ! dummy_star = STEL_SPEC(stel_idx)
          ! END DO
          ! DO J = 1, NLAYER
          !     total_layer_taus(J) = sum(tau_IRe(:,J) * dummy_weights)
          ! END DO

      end subroutine opacity_wrapper_corrk

      subroutine calculate_opacities_corrk(NLAYER, NSOL, NIR, Tl, Pl, dpe, tau_IRe, tau_Ve, gravity_SI, k_IRl, k_Vl,
     &                                   OPAC_CORRK, TS_CORRK, PS_CORRK, TS_LOG_CORRK, PS_LOG_CORRK,
     &                                   tau_ray_temp, TAURAY_PER_DPG, NWNO, NTGRID, NPGRID)

        implicit none
        integer :: k, NLAYER, J, NSOL, NIR, i, NWNO, NTGRID, NPGRID

        real :: R, gravity_SI

        real, dimension(NLAYER) :: dpe, Pl, Tl, pe
        real, dimension(NIR, NLAYER) :: k_IRl
        real, dimension(NSOL,NLAYER) :: k_Vl
        real, dimension(NIR,NLAYER+1) :: tau_IRe
        real, dimension(NSOL,NLAYER+1) :: tau_Ve, tau_ray_temp
        ! corrk commons:
        REAL :: OPAC_CORRK(NTGRID, NPGRID, NWNO, 8), TAURAY_PER_DPG(NTGRID, NPGRID, NWNO)
        REAL :: TS_CORRK(NTGRID), PS_CORRK(NPGRID), TS_LOG_CORRK(NTGRID), PS_LOG_CORRK(NPGRID)
        ! WRITE(*,*) "OPAC_CORRK: ", OPAC_CORRK(1,1,1,1), OPAC_CORRK(1,1,2,1), OPAC_CORRK(1,1,3,1)
        tau_Ve(:,1) = 0.0
        tau_IRe(:,1) = 0.0
        ! write(*,*) 'loc(tau_ray_temp)1: ', LOC(tau_ray_temp)
        do k = 1, NLAYER
          call local_opacities_corrk(Tl(k), Pl(k), k_IRl, k_Vl, OPAC_CORRK, TS_CORRK, PS_CORRK, TS_LOG_CORRK, PS_LOG_CORRK, k,
     &                                 NLAYER, NIR, NSOL, tau_ray_temp, TAURAY_PER_DPG, NWNO, NTGRID, NPGRID)
          ! double-gray overwrite:
          ! k_IRl(:,k) = 1.e-3 ! should be m^2/kg (= cm^2/g / 10)
          ! k_Vl(:,k) = 1.e-4
          ! DPG bug fix adds:
          k_Vl(:,k) = k_IRl(:,k) ! spectral, so these are interchangeable

          ! DPG bug fix removes:
          ! tau_IRe(:,k) = ((k_IRl(:,k) * dpe(k)) / gravity_SI) ! k, dpe, and gravity_SI are all in SI units
          ! tau_Ve(:,k)  = tau_IRe(:,k) ! spectral, so these are the same 
          ! tau_Ve(:,k) = ((k_Vl(:,k) * dpe(k)) / gravity_SI) ! k, dpe, and gravity_SI are all in SI units
          ! tau_ray_temp(:,k) = tau_ray_temp(:,k) !* dpe(k) / gravity_SI

          
        end do
      end subroutine calculate_opacities_corrk

      subroutine local_opacities_corrk(Tin, Pin, k_IRl, k_Vl, OPAC_CORRK, TS_CORRK, PS_CORRK, TS_LOG_CORRK, PS_LOG_CORRK, k,
     &                                 NLAYER, NIR, NSOL, tau_ray_temp, TAURAY_PER_DPG, NWNO, NTGRID, NPGRID)
        use corrkmodule, only : NEAREST_INDEX
        implicit none
        real :: Tin, Pin
        integer :: NLAYER, NIR, NSOL, I, NWNO, NTGRID, NPGRID
        real, dimension(NIR, NLAYER) :: k_IRl
        real, dimension(NSOL, NLAYER) :: k_Vl, tau_ray_temp
        ! For now, interpolating molceular opacities (cm^2/molecule) linearly in both pressure and temperature
        ! Since we're explicitly spectral now, I think K_IRl and K_Vl are the same thing
        ! chan_idx is the overall idx in the flat GCM stuff, so should always be wave_idx * 8 + gauss_idx
        integer :: T_idx, P_idx, wave_idx, gauss_idx, chan_idx, k
        ! corrk commons:
        REAL :: OPAC_CORRK(NTGRID, NPGRID, NWNO, 8), TAURAY_PER_DPG(NTGRID, NPGRID, NWNO)
        REAL :: TS_CORRK(NTGRID), PS_CORRK(NPGRID), TS_LOG_CORRK(NTGRID), PS_LOG_CORRK(NPGRID)
        ! WRITE(*,*) "OPAC_CORRK: ", OPAC_CORRK(1,1,1,1), OPAC_CORRK(1,1,2,1), OPAC_CORRK(1,1,3,1)

        ! NTGRID = 73
        ! NPGRID = 20

        T_idx = NEAREST_INDEX(TS_CORRK, NTGRID, Tin)
        P_idx = NEAREST_INDEX(PS_CORRK, NPGRID, Pin)

        ! Nearest-neighbor the rayleigh scattering optical depth (faster than bilinear and accurate enough)
        do chan_idx = 1, NSOL

          gauss_idx = modulo(chan_idx - 1, 8) + 1
          wave_idx = MODULO((chan_idx - gauss_idx)/8,NWNO) + 1

          tau_ray_temp(chan_idx, k) = TAURAY_PER_DPG(T_idx, P_idx, wave_idx)
        end do
        ! Set up indices for bilinear interpolation
        if (PS_CORRK(P_idx) .gt. Pin) then
          P_idx = P_idx - 1
        end if
        if (TS_CORRK(T_idx) .gt. Tin) then
          T_idx = T_idx - 1
        end if

        ! Clamp to valid range so T_idx+1/P_idx+1 stay in bounds; Tin/Pin outside
        ! the table get the edge cell (i.e. extrapolation by clamping)
        T_idx = max(1, min(T_idx, NTGRID-1))
        P_idx = max(1, min(P_idx, NPGRID-1))


        do chan_idx = 1, NSOL
          gauss_idx = modulo(chan_idx - 1, 8) + 1
          wave_idx = MODULO((chan_idx - gauss_idx)/8,NWNO) + 1
          ! k_IRl(chan_idx, k) = OPAC_CORRK(T_idx, P_idx, wave_idx, gauss_idx)
          ! interpolate molecular line opacities

          call bilinear_interpolation(LOG10(Tin), log10(Pin), TS_LOG_CORRK(T_idx), TS_LOG_CORRK(T_idx+1), PS_LOG_CORRK(P_idx), 
     &                                  PS_LOG_CORRK(P_idx+1), 
     &                                  OPAC_CORRK(T_idx, P_idx, wave_idx, gauss_idx), OPAC_CORRK(T_idx+1, P_idx, 
     &                                  wave_idx, gauss_idx), OPAC_CORRK(T_idx, P_idx+1, wave_idx, gauss_idx), 
     &                                  OPAC_CORRK(T_idx+1, P_idx+1, wave_idx, gauss_idx), k_IRl(chan_idx, k))
          ! interpolate CIA opacities (now rolled into molecular opacities)
    !       call bilinear_interpolation(LOG10(Tin), log10(Pin), TS_LOG_CORRK(T_idx), TS_LOG_CORRK(T_idx+1), PS_LOG_CORRK(P_idx), 
    !  &                                  PS_LOG_CORRK(P_idx+1), 
    !  &                                  OPAC_CIA(T_idx, P_idx, wave_idx), OPAC_CIA(T_idx+1, P_idx, wave_idx),
    !  &                                  OPAC_CIA(T_idx, P_idx+1, wave_idx), OPAC_CIA(T_idx+1, P_idx+1, wave_idx), 
    !  &                                  k_CIA(chan_idx, k))

          ! interpolate Rayleigh scattering opacities (This is now nearest-neighbor, see above)
    !       call bilinear_interpolation(LOG10(Tin), log10(Pin), TS_LOG_CORRK(T_idx), TS_LOG_CORRK(T_idx+1), PS_LOG_CORRK(P_idx), 
    !  &                       PS_LOG_CORRK(P_idx+1), 
    !  &                       TAURAY_PER_DPG(T_idx, P_idx, wave_idx), TAURAY_PER_DPG(T_idx+1, P_idx, wave_idx),
    !  &                       TAURAY_PER_DPG(T_idx, P_idx+1, wave_idx), TAURAY_PER_DPG(T_idx+1, P_idx+1, wave_idx),
    !  &                       tau_ray_temp(chan_idx, k))

          k_IRl(chan_idx, k) = 10**k_IRl(chan_idx, k) ! take out of logspace
          ! write(*,*), 'Pin, Tin, k_IRl: ', Pin, Tin, k_IRl(chan_idx, :)
          tau_ray_temp(chan_idx, k) = 10**tau_ray_temp(chan_idx, k) ! should be SI, m^2/kg

        end do


      end subroutine local_opacities_corrk


      subroutine bilinear_interpolation(x, y, x1, x2, y1, y2, Q11, Q12, Q21, Q22, result)
        implicit none
        real, intent(in) :: x, y
        real, intent(in) :: x1, x2, y1, y2
        real, intent(in) :: Q11, Q12, Q21, Q22
        real, intent(out) :: result
        real :: R1, R2
    
        ! Perform the bilinear interpolation
        R1 = ((x2 - x) / (x2 - x1)) * Q11 + ((x - x1) / (x2 - x1)) * Q21
        R2 = ((x2 - x) / (x2 - x1)) * Q12 + ((x - x1) / (x2 - x1)) * Q22
        result = ((y2 - y) / (y2 - y1)) * R1 + ((y - y1) / (y2 - y1)) * R2
      end subroutine bilinear_interpolation