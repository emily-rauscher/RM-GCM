      SUBROUTINE TWO_STREAM_WRAPPER(MEAN_HEMISPHERIC_INTENSITY_UP, MEAN_HEMISPHERIC_INTENSITY_DOWN,
     &                              MEAN_HEMISPHERIC_FLUX_UP, MEAN_HEMISPHERIC_FLUX_DOWN,
     &                              QUADRATURE_FLUX, QUADRATURE_INTENSITY, Beta_V, Beta_IR, term1,
     &                              total_heat_ir, total_heat_vi, HEMISPHERIC_INTENSITY_DOWN, HEMISPHERIC_INTENSITY_UP)
          include 'rcommons.h'

          integer :: NVERT, J, I
          real, dimension(NVERT) :: SINGLE_W0S
          real, dimension(NVERT) :: SINGLE_G0S
          real, dimension(NVERT) :: SINGLE_TAULS
          real, dimension(NVERT) :: SINGLE_TEMPS

          real, dimension(NVERT) :: MEAN_HEMISPHERIC_INTENSITY_UP, MEAN_HEMISPHERIC_INTENSITY_DOWN
          real, dimension(NVERT) :: MEAN_HEMISPHERIC_FLUX_UP, MEAN_HEMISPHERIC_FLUX_DOWN
          real, dimension(NVERT) :: QUADRATURE_FLUX, QUADRATURE_INTENSITY
          real, dimension(NVERT) :: HEMISPHERIC_INTENSITY_DOWN, HEMISPHERIC_INTENSITY_UP

          real :: mu_0, NU, FLUX_SURFACE_QUADRATURE, EMIS_VAL, RSFX_VAL, Beta_IR_VAL, Beta_V_VAL, term1
          real, dimension(2) :: Beta_IR
          real, dimension(3) :: Beta_V

          real, dimension(NVERT) :: total_heat_ir
          real, dimension(NVERT) :: total_heat_vi


          DO J = 1, NVERT
            SINGLE_W0S(J) = W0(2,J)
            SINGLE_G0S(J) = G0(2,J)

            SINGLE_TAULS(J) = TAUL(1,J)
            SINGLE_TEMPS(J) = TT(J)
          END DO

          ! This is to calculate the incident fraction of the starlight
          mu_0 = COS(ALAT * PI / 180.0) * COS(ALON * PI / 180.0)

          IF (mu_0 .le. 0) THEN
            mu_0 = 0
          END IF

          NU = 1e13
          FLUX_SURFACE_QUADRATURE = 1.0

          total_heat_ir(1) = 0.0
          DO I = 1, 1
              NU   = 1e13
              EMIS_VAL = EMIS(2)
              RSFX_VAL = RSFX(2)
              !Beta_IR_VAL = Beta_IR(I)
              Beta_IR_VAL = 1.0
              CALL GET_IR_INTENSITY(MEAN_HEMISPHERIC_INTENSITY_UP, MEAN_HEMISPHERIC_INTENSITY_DOWN,
     &                              MEAN_HEMISPHERIC_FLUX_UP, MEAN_HEMISPHERIC_FLUX_DOWN,
     &                              NU, EMIS_VAL, RSFX_VAL, NVERT, SINGLE_TEMPS, SINGLE_TAULS, SINGLE_W0S, SINGLE_G0S,
     &                              Beta_IR_VAL, HEMISPHERIC_INTENSITY_DOWN, HEMISPHERIC_INTENSITY_UP)
              DO J = 1, NVERT
                   total_heat_ir(J) = total_heat_ir(J)  +
     &                        ((MEAN_HEMISPHERIC_INTENSITY_UP(J+1)/2 - MEAN_HEMISPHERIC_INTENSITY_DOWN(J+1)/2) -
     &                         (MEAN_HEMISPHERIC_INTENSITY_UP(J)/2   - MEAN_HEMISPHERIC_INTENSITY_DOWN(J)/2)) * TERM1
              END DO
          END DO



 !         DO I = 1, 3
 !             EMIS_VAL = EMIS(1)
 !             RSFX_VAL = RSFX(1)
 !             Beta_V_VAL = Beta_V(I)
 !             CALL GET_SOLAR_INTENSITY(QUADRATURE_FLUX, QUADRATURE_INTENSITY,
 !    &                                NU, FLUX_SURFACE_QUADRATURE,
 !    &                                EMIS(2), RSFX(2), NVERT, TAULS, W0, G0, mu_0, Beta_V_VAL)
!
!              DO J = 1, NVERT
!                  total_heat_vi(J) = total_heat_vi(J) + (QUADRATURE_FLUX(J+1) - QUADRATURE_FLUX(J))*TERM1
!              END DO!
!
 !         END DO
!
                 
      END SUBROUTINE TWO_STREAM_WRAPPER
      




      
      SUBROUTINE GET_IR_INTENSITY(MEAN_HEMISPHERIC_INTENSITY_UP, MEAN_HEMISPHERIC_INTENSITY_DOWN,
     &                            MEAN_HEMISPHERIC_FLUX_UP, MEAN_HEMISPHERIC_FLUX_DOWN,   
     &                            NU, EMIS, RSFX, NVERT, TEMPS, TAULS, W0, G0, Beta_IR,
     &                            HEMISPHERIC_INTENSITY_DOWN, HEMISPHERIC_INTENSITY_UP)

          integer :: NVERT, J, L, Z, NGAUSS
          real :: NU, mu_1, mu, temp_val_2, BB_TOP_OF_ATM, BB_BOTTOM_OF_ATM, SFCS_HEMISPHERIC
          real(16) :: temp_val_1
          real :: e_con, pi, bolz_constant, clight, h_constant

          real, dimension(NVERT) :: TAULS, W0, G0, TEMPS

          real ::  EMIS, RSFX
          real :: Beta_IR
          real, dimension(NVERT) :: y1, y2
          real, dimension(NVERT) :: LAMBDAS, GAMMA
          real, dimension(NVERT) :: temp_e_val
          real, dimension(NVERT) :: e1, e2, e3, e4
          real, dimension(NVERT) :: temp_gamma_val
          real, dimension(NVERT) :: CP, CPB, CM, CMB

          real, dimension(2*NVERT) :: A, B, D, E
          real, dimension(2*NVERT) :: AS, DS, X, Y

          real, dimension(NVERT) :: B0, B1
          real, dimension(NVERT) :: SOURCE_G, SOURCE_H, SOURCE_J, SOURCE_K
          real, dimension(NVERT) :: ALPHA_1, ALPHA_2, SIGMA_1, SIGMA_2
          real, dimension(NVERT) :: SOURCE_Y1, SOURCE_Y2, source_temp
          real, dimension(NVERT) :: HEMISPHERIC_INTENSITY_DOWN, HEMISPHERIC_INTENSITY_UP

          real, dimension(NVERT) :: MEAN_HEMISPHERIC_INTENSITY_UP, MEAN_HEMISPHERIC_INTENSITY_DOWN
          real, dimension(NVERT) :: MEAN_HEMISPHERIC_FLUX_UP, MEAN_HEMISPHERIC_FLUX_DOWN

          real, dimension(3) :: GAUSS_ANGLES, GUASS_WEIGHTS, GAUSS_RATIOS
          integer :: I

          DO J = 1, NVERT
            MEAN_HEMISPHERIC_INTENSITY_UP(J) = 0
            MEAN_HEMISPHERIC_INTENSITY_DOWN(J) = 0
            MEAN_HEMISPHERIC_FLUX_UP(J) = 0
            MEAN_HEMISPHERIC_FLUX_DOWN(J) = 0
          END DO

          GAUSS_ANGLES  = (/0.2123405382, 0.5905331356, 0.9114120405/)
          GUASS_WEIGHTS = (/0.0698269799, 0.2292411064, 0.2009319137/)
          GAUSS_RATIOS  = (/0.4679139346, 0.3607615730, 0.1713244924/)

          mu   = 1.0
          NGAUSS = 3

          bolz_constant = 1.380649e-23
          h_constant    = 6.62607015e-34
          clight        = 3e8
          pi            = 3.141592653589
          e_con         = 2.718281828459

          !temp_val_1 = 2.0 * NU * NU / (CLIGHT * CLIGHT)
          !temp_val_1 = temp_val_1 * (h_constant)
          !temp_val_1 = temp_val_1 * NU
          !temp_val_2 = e_con ** (h_constant * NU / (bolz_constant * TEMPS(1))) - 1.0
          !BB_TOP_OF_ATM = REAL(temp_val_1 * (1.0 / temp_val_2))

          !temp_val_1 = 2.0 * NU * NU / (CLIGHT * CLIGHT)
          !temp_val_1 = temp_val_1 * (h_constant)
          !temp_val_1 = temp_val_1 * NU
          !temp_val_2 = e_con ** (h_constant * NU / (bolz_constant * TEMPS(NVERT))) - 1.0
          !BB_BOTTOM_OF_ATM = REAL(temp_val_1 * (1.0 / temp_val_2))

          BB_TOP_OF_ATM    = TEMPS(1) * TEMPS(1) * TEMPS(1) * TEMPS(1) * 1.8049459031e-8 * Beta_IR
          BB_BOTTOM_OF_ATM = TEMPS(NVERT) * TEMPS(NVERT) * TEMPS(NVERT) * TEMPS(NVERT) *
     &                       1.8049459031e-8 * Beta_IR

          SFCS_HEMISPHERIC = EMIS * PI * BB_BOTTOM_OF_ATM
          mu_1 = 0.5

          DO J = 1, NVERT
              y1(J) = (2.0 - (W0(J) * (1.0 + G0(J))));
              y2(J) = (W0(J) * (1.0 - G0(J)));

              LAMBDAS(J)    =  (y1(J)*y1(J) - y2(J)*y2(J) ) ** 0.5
              GAMMA(J)      =  y2(J) / (y1(J) + LAMBDAS(J));
              temp_e_val(J) =  e_con ** (-LAMBDAS(J) * TAULS(J));

              e1(J) = 1.0 + GAMMA(J) * temp_e_val(J);
              e2(J) = 1.0 - GAMMA(J) * temp_e_val(J);
              e3(J) = GAMMA(J) + temp_e_val(J);
              e4(J) = GAMMA(J) - temp_e_val(J);
          END DO

          J = 1
          DO L = 2, 2*NVERT-1, 2
              ! HERE ARE THE EVEN MATRIX ELEMENTS
              A(L)   =  e2(J+1) * e1(J)   - e3(J) * e4(J+1)
              B(L)   =  e2(J+1) * e2(J)   - e4(J+1) * e4(J)
              D(L)   =  e1(J+1) * e4(J+1) - e3(J+1) * e2(J+1)

              ! HERE ARE THE ODD MATRIX ELEMENTS EXCEPT FOR THE TOP.
              A(L+1) =  e2(J)   * e3(J)   - e1(J)   * e4(J)
              B(L+1) =  e1(J+1) * e1(J)   - e3(J+1) * e3(J)
              D(L+1) =  e3(J)   * e4(J+1) - e1(J)   * e2(J+1)
              J = J + 1
          END DO

          ! HERE ARE THE TOP AND BOTTOM BOUNDARY CONDITIONS AS WELL AS THE
          ! BEGINNING OF THE TRIDIAGONAL SOLUTION DDINITIONS. I ASSUME
          ! NO DIFFUSE RADIATION IS INCIDENT AT UPPER BOUNDARY.
          A(1) = 0.0;
          B(1) = e1(1);
          D(1) = -e2(1);

          A(2*NVERT) = e1(NVERT) - RSFX * e3(NVERT);
          B(2*NVERT) = e2(NVERT) - RSFX * e4(NVERT);
          D(2*NVERT) = 0.0;

          ! This is the part of the code that solves for the blackbody stuff
          DO J=1, NVERT
              IF (1 .GE. J-1) THEN
                  KINDEX = 1
              ELSE
                  KINDEX = J-1
              END IF

              ! STUPID HACK
              ! PREVENTS NU ** 3 from being too big
              ! OMG FORTRAN
              !temp_val_1 = 2.0 * NU * NU / (CLIGHT * CLIGHT)
              !temp_val_1 = temp_val_1 * (h_constant)
              !temp_val_1 = temp_val_1 * NU
              !temp_val_2 = e_con ** (h_constant * NU / (bolz_constant * TEMPS(J))) - 1.0
              !B0(J) = REAL(temp_val_1 * (1.0 / temp_val_2))

              ! This is if you want double gray, you'll also have to fix boundary stuff
              B0(J) = TEMPS(J) * TEMPS(J) * TEMPS(J) * TEMPS(J) * 1.8049459031e-8 * Beta_IR
              B1(J) = (B0(J) - B0(KINDEX)) / TAULS(J)
          END DO

          DO J=1, NVERT
              IF (1 .GE. J-1) THEN
                  KINDEX = 1
              ELSE
                  KINDEX = J-1
              END IF

              temp_gamma_val(J)   = 1.0 / (y1(J) + y2(J));

              CP(J)  = (B0(KINDEX) + B1(J) * temp_gamma_val(J)) * 2.0 * PI * mu_1
              CPB(J) = CP(J) + B1(J) * TAULS(J) * 2.0 * PI * mu_1

              CM(J)  = (B0(KINDEX) - B1(J) * temp_gamma_val(J)) * 2.0 * PI * mu_1
              CMB(J) = CM(J) + B1(J) * TAULS(J) * 2.0 * PI * mu_1
          END DO

          J = 1
          DO L = 2, 2*NVERT-1, 2
              E(L)   = (CP(J+1) - CPB(J)) * e2(J+1) + (CMB(J) - CM(J+1))*e4(J+1)
              E(L+1) = e3(J) * (CP(J+1) - CPB(J)) + e1(J) * (CMB(J) - CM(J+1))
              J = J + 1
          END DO

          ! HERE ARE THE TOP AND BOTTOM BOUNDARY CONDITIONS AS WELL AS THE
          ! BEGINNING OF THE TRIDIAGONAL SOLUTION. I ASSUME NO
          ! DIFFUSE RADIATION IS INCIDENT AT THE TOP.
          E(1) = 0.0 - CM(1);
          E(2*NVERT)  = SFCS_HEMISPHERIC - CPB(NVERT) + RSFX * CMB(NVERT)

          DS(2*NVERT) = E(2*NVERT) / B(2*NVERT)
          AS(2*NVERT) = A(2*NVERT) / B(2*NVERT)

          DO L = 1, 2*NVERT - 1
              X(2*NVERT-L)  = 1.0 / (B(2*NVERT-L) - D(2*NVERT-L) * AS(2*NVERT-L+1))
              AS(2*NVERT-L) = A(2*NVERT-L) * X(2*NVERT-L)
              DS(2*NVERT-L) = (E(2*NVERT-L) - (D(2*NVERT-L) * DS(2*NVERT-L+1))) * X(2*NVERT-L)
          END DO

          Y(1) = DS(1)

          DO L = 2, 2*NVERT
              Y(L) = DS(L) - AS(L) * Y(L-1)
          END DO

          DO J = 2, NVERT+1
              SOURCE_Y1(J-1) = Y(2*J-3)
              SOURCE_Y2(J-1) = Y(2*J-2)
          END DO

          DO J=1, NVERT
              SOURCE_G(J) = (SOURCE_Y1(J) + SOURCE_Y2(J)) * (1.0/mu_1 - LAMBDAS(J))
              SOURCE_H(J) = (SOURCE_Y1(J) - SOURCE_Y2(J)) * GAMMA(J) * (1.0/mu_1 + LAMBDAS(J))
              SOURCE_J(J) = (SOURCE_Y1(J) + SOURCE_Y2(J)) * GAMMA(J) * (1.0/mu_1 + LAMBDAS(J))
              SOURCE_K(J) = (SOURCE_Y1(J) - SOURCE_Y2(J)) * (1.0/mu_1 - LAMBDAS(J))

              source_temp(J) = (1.0 / (y1(J) + y2(J))) - mu_1

              ALPHA_1(J)     = 2 * pi * (B0(J) + (B1(J) * source_temp(J)))
              ALPHA_2(J)     = 2 * pi * (B1(J))

              SIGMA_1(J) = 2.0 * pi * (B0(J) - (B1(J) * source_temp(J)))
              SIGMA_2(J) = 2.0 * pi * B1(J)
          END DO

          DO I=1, NGAUSS
            mu = GAUSS_ANGLES(I)

            !BB_TOP_OF_ATM
            HEMISPHERIC_INTENSITY_DOWN(1) = BB_TOP_OF_ATM * e_con ** (-TAULS(1)/mu) +
     &                      SOURCE_J(1)/(LAMBDAS(1)*mu + 1.0) * (1.0 - e_con ** (-TAULS(1)*(LAMBDAS(1)+1.0/mu))) + 
     &                      SOURCE_K(1)/(LAMBDAS(1)*mu - 1.0) * (e_con ** (-TAULS(1)/mu)-e_con**(-TAULS(1)*LAMBDAS(1)))+
     &                      SIGMA_1(1) * (1.0 - e_con ** (-TAULS(1)/mu)) + 
     &                      SIGMA_2(1) * (mu * e_con ** (-TAULS(1)/mu) + TAULS(1) - mu);

            MEAN_HEMISPHERIC_INTENSITY_DOWN(1) = MEAN_HEMISPHERIC_INTENSITY_DOWN(1) + 
     &                                           GAUSS_RATIOS(I) * HEMISPHERIC_INTENSITY_DOWN(1)
            MEAN_HEMISPHERIC_FLUX_DOWN(1)      = MEAN_HEMISPHERIC_INTENSITY_DOWN(1) + 
     &                                           GUASS_WEIGHTS(I) * HEMISPHERIC_INTENSITY_DOWN(1)

            DO J = 2, NVERT
              HEMISPHERIC_INTENSITY_DOWN(J) = HEMISPHERIC_INTENSITY_DOWN(J-1) * e_con ** (-TAULS(J)/mu) + 
     &                       SOURCE_J(J)/(LAMBDAS(J)*mu + 1.0) * (1.0 - e_con ** (-TAULS(J)*(LAMBDAS(J)+1.0/mu))) + 
     &                       SOURCE_K(J)/(LAMBDAS(J)*mu-1.0) * (e_con ** (-TAULS(J)/mu)-e_con**(-TAULS(J)*LAMBDAS(J)))+
     &                       SIGMA_1(J) * (1.0 - e_con ** (-TAULS(J)/mu)) + 
     &                       SIGMA_2(J) * (mu * e_con ** (-TAULS(J)/mu) + TAULS(J) - mu)

              MEAN_HEMISPHERIC_INTENSITY_DOWN(J) = MEAN_HEMISPHERIC_INTENSITY_DOWN(J) +
     &                                             GAUSS_RATIOS(I) * HEMISPHERIC_INTENSITY_DOWN(J)
              MEAN_HEMISPHERIC_FLUX_DOWN(J)      = MEAN_HEMISPHERIC_FLUX_DOWN(J) + GUASS_WEIGHTS(I) *
     &                                             HEMISPHERIC_INTENSITY_DOWN(J)

            END DO

            HEMISPHERIC_INTENSITY_UP(NVERT)      = 2.0 * BB_BOTTOM_OF_ATM * EMIS * pi
            MEAN_HEMISPHERIC_INTENSITY_UP(NVERT) = MEAN_HEMISPHERIC_INTENSITY_UP(NVERT) + 
     &                                              GAUSS_RATIOS(I) * HEMISPHERIC_INTENSITY_UP(NVERT)
            MEAN_HEMISPHERIC_FLUX_UP(NVERT)      = MEAN_HEMISPHERIC_FLUX_UP(NVERT) +
     &                                              GUASS_WEIGHTS(I) * HEMISPHERIC_INTENSITY_UP(NVERT)

            !Calculate the upward intensity next
            DO Z = 1, NVERT - 1
              J = NVERT - Z
              HEMISPHERIC_INTENSITY_UP(J) = HEMISPHERIC_INTENSITY_UP(J+1) * e_con ** (-TAULS(J+1)/mu) + 
     &                           SOURCE_G(J)/(LAMBDAS(J)*mu-1.0)*(e_con ** (-TAULS(J+1)/mu)-e_con **
     &                           (-TAULS(J+1)*(LAMBDAS(J)))) +
     &                           SOURCE_H(J)/(LAMBDAS(J)*mu+1.0)*(1.0 - e_con ** (-TAULS(J+1)*(LAMBDAS(J) + 1.0/mu)))+
     &                           ALPHA_1(J) * (1.0 - e_con ** (-TAULS(J+1)/mu)) + 
     &                           ALPHA_2(J) * (mu - ((TAULS(J+1) + mu) * (e_con ** (-TAULS(J+1)/mu))))

              MEAN_HEMISPHERIC_INTENSITY_UP(J) = MEAN_HEMISPHERIC_INTENSITY_UP(J) +
     &                                           GAUSS_RATIOS(I) * HEMISPHERIC_INTENSITY_UP(J)
              MEAN_HEMISPHERIC_FLUX_UP(J)      = MEAN_HEMISPHERIC_FLUX_UP(J) +
     &                                           GUASS_WEIGHTS(I) * HEMISPHERIC_INTENSITY_UP(J)
            END DO
          END DO
      END SUBROUTINE GET_IR_INTENSITY



      SUBROUTINE GET_SOLAR_INTENSITY(QUADRATURE_FLUX, QUADRATURE_INTENSITY,
     &                               NU, FLUX_SURFACE_QUADRATURE_TOT,
     &                               EMIS, RSFX, NVERT, TAULS, W0, G0, mu_0, Beta_V)
          integer :: NVERT, J, L
          real :: RSFX, mu_0, mu_1, NU, FLUX_SURFACE_QUADRATURE, FLUX_SURFACE_QUADRATURE_TOT

          real :: e_con, pi
          real :: Beta_V
          real, dimension(NVERT) :: QUADRATURE_FLUX, QUADRATURE_INTENSITY, TAULS, W0, G0
          real, dimension(NVERT+1) :: TAUCS

          real, dimension(NVERT) :: DIRECT_QUADRATURE
          real, dimension(NVERT) :: y1, y2, y3, y4
          real, dimension(NVERT) :: LAMBDAS, GAMMA
          real, dimension(NVERT) :: temp_e_val
          real, dimension(NVERT) :: e1, e2, e3, e4
          real, dimension(NVERT) :: temp_gamma_val
          real, dimension(NVERT) :: CP, CPB, CM, CMB

          real, dimension(2*NVERT) :: A, B, D, E
          real, dimension(2*NVERT) :: AS, DS, X, Y

          e_con = 2.718281828459
          pi    = 3.141592653589

          FLUX_SURFACE_QUADRATURE = FLUX_SURFACE_QUADRATURE_TOT * Beta_V

          TAUCS(1) = 0.0
          DO N = 1, NVERT
              TAUCS(N+1) = TAUCS(N) + TAULS(N)
          END DO

          DO J = 1, NVERT
              DIRECT_QUADRATURE(J) = mu_0 * pi * FLUX_SURFACE_QUADRATURE*e_con ** (-1.0 * (TAUCS(J) + TAULS(J)) / mu_0)
          END DO

          ! Boundary conditions
          mu_1 = 0.577350;
          SFCS_QUADRATURE = RSFX * mu_0 * e_con ** (-(TAUCS(NVERT + 1))/mu_0) * pi * FLUX_SURFACE_QUADRATURE

          ! HERE WE FIND LAYER PROPERTIES FOLLOWING GENERAL SCHEME
          ! OF MEADOR AND WEAVOR. THEN WE SET UP LAYER PROPERTIES
          ! NEEDED FOR MATRIX

          DO J = 1, NVERT
              y1(J) = 0.86602540378 * (2.0 - (W0(J) * (1.0 + G0(J))))
              y2(J) = 0.86602540378 * W0(J) * (1.0 - G0(J))
              y3(J) = 0.5 * (1.0 - 1.73205080756 * mu_0 * G0(J))
              y4(J) = 1.0 - (y3(J))

              LAMBDAS(J)    =  (y1(J) ** 2.0 - y2(J) ** 2.0) ** 0.5
              GAMMA(J)      =  y2(J) / (y1(J) + LAMBDAS(J))
              temp_e_val(J) =  e_con ** (-LAMBDAS(J) * TAULS(J))

              e1(J)   =  1.0 + GAMMA(J) * temp_e_val(J)
              e2(J)   =  1.0 - GAMMA(J) * temp_e_val(J)
              e3(J)   =  GAMMA(J) + temp_e_val(J)
              e4(J)   =  GAMMA(J) - temp_e_val(J)
          END DO

          J = 1
          DO L = 2, 2*NVERT-1, 2
              ! HERE ARE THE EVEN MATRIX ELEMENTS
              A(L)   =  e2(J+1) * e1(J)   - e3(J) * e4(J+1)
              B(L)   =  e2(J+1) * e2(J)   - e4(J+1) * e4(J)
              D(L)   =  e1(J+1) * e4(J+1) - e3(J+1) * e2(J+1)

              ! HERE ARE THE ODD MATRIX ELEMENTS EXCEPT FOR THE TOP.
              A(L+1) =  e2(J)   * e3(J)   - e1(J)   * e4(J)
              B(L+1) =  e1(J+1) * e1(J)   - e3(J+1) * e3(J)
              D(L+1) =  e3(J)   * e4(J+1) - e1(J)   * e2(J+1)
              J = J + 1
          END DO

          ! HERE ARE THE TOP AND BOTTOM BOUNDARY CONDITIONS AS WELL AS THE
          ! BEGINNING OF THE TRIDIAGONAL SOLUTION DDINITIONS. I ASSUME
          ! NO DIFFUSE RADIATION IS INCIDENT AT UPPER BOUNDARY.
          A(1) = 0.0;
          B(1) = e1(1);
          D(1) = -e2(1);

          A(2*NVERT) = e1(NVERT) - RSFX * e3(NVERT);
          B(2*NVERT) = e2(NVERT) - RSFX * e4(NVERT);
          D(2*NVERT) = 0.0;

          DO J = 1, NVERT
              temp_gamma_val(J)   = 1.0 / (y1(J) + y2(J))

              CP(J) =  W0(J) * PI * FLUX_SURFACE_QUADRATURE * 
     &               e_con ** (-(TAUCS(J)) / mu_0) * 
     &               (((y1(J) - 1.0 / mu_0) * y3(J)) + (y4(J) * y2(J))) / 
     &               (LAMBDAS(J) ** 2.0 - (1.0 / (mu_0 ** 2.0)))


              CPB(J) = W0(J) * PI * FLUX_SURFACE_QUADRATURE * 
     &               e_con ** (-(TAUCS(J) + TAULS(J)) / mu_0) *   
     &               (((y1(J) - 1.0 / mu_0) * y3(J)) + (y4(J) * y2(J))) / 
     &               (LAMBDAS(J) ** 2.0 - (1.0 / (mu_0 ** 2.0)))


              CM(J) = W0(J) * PI * FLUX_SURFACE_QUADRATURE * 
     &             e_con ** (-(TAUCS(J)) / mu_0) *   
     &             (((y1(J) + 1.0 / mu_0) * y4(J)) + (y2(J) * y3(J))) / 
     &             ((LAMBDAS(J) ** 2.0) - (1.0 / (mu_0 ** 2.0)))


              CMB(J) = W0(J) * PI * FLUX_SURFACE_QUADRATURE * 
     &              e_con ** (-(TAUCS(J) + TAULS(J)) / mu_0) *   
     &               (((y1(J) + 1.0 / mu_0) * y4(J)) + (y2(J) * y3(J))) / 
     &               ((LAMBDAS(J) ** 2.0) - (1.0 / (mu_0 ** 2.0)))
          END DO

          J = 1
          DO L = 2, 2*NVERT-1, 2
              E(L)   = (CP(J+1) - CPB(J)) * e2(J+1) + (CMB(J) - CM(J+1))*e4(J+1)
              E(L+1) = e3(J) * (CP(J+1) - CPB(J)) + e1(J) * (CMB(J) - CM(J+1))
              J = J + 1
          END DO

          ! HERE ARE THE TOP AND BOTTOM BOUNDARY CONDITIONS AS WELL AS THE
          ! BEGINNING OF THE TRIDIAGONAL SOLUTION. I ASSUME NO
          ! DIFFUSE RADIATION IS INCIDENT AT THE TOP.
          E(1) = 0.0 - CM(1);
          E(2*NVERT)  = SFCS_QUADRATURE - CPB(NVERT) + RSFX * CMB(NVERT)

          DS(2*NVERT) = E(2*NVERT) / B(2*NVERT)
          AS(2*NVERT) = A(2*NVERT) / B(2*NVERT)

          DO L = 1, 2*NVERT - 1
              X(2*NVERT-L)  = 1.0 / (B(2*NVERT-L) - D(2*NVERT-L) * AS(2*NVERT-L+1))
              AS(2*NVERT-L) = A(2*NVERT-L) * X(2*NVERT-L)
              DS(2*NVERT-L) = (E(2*NVERT-L) - (D(2*NVERT-L) * DS(2*NVERT-L+1))) * X(2*NVERT-L)
          END DO

          Y(1) = DS(1)

          DO L = 2, 2*NVERT
              Y(L) = DS(L) - AS(L) * Y(L-1)
          END DO


          DO J = 1, NVERT
              QUADRATURE_FLUX(J) = Y(2*J-1)*(e1(J)-e3(J))+Y(2*J)*(e2(J)-e4(J)) + CPB(J) - CMB(J)-DIRECT_QUADRATURE(J)
              QUADRATURE_INTENSITY(J) = (1.0 / mu_1) * (Y(2*J-1)*(e1(J) + e3(J)) + 
     &                   (Y(2*J) * (e2(J) + e4(J))) + CPB(J) + CMB(J))+DIRECT_QUADRATURE(J)/mu_0
          END DO
      RETURN
      END SUBROUTINE GET_SOLAR_INTENSITY