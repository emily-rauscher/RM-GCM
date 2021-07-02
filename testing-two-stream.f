      SUBROUTINE TWO_STREAM_WRAPPER
          include 'rcommons.h'

          integer :: NLAYER, J
          real, dimension(NLAYER) :: SINGLE_W0S
          real, dimension(NLAYER) :: SINGLE_G0S
          real, dimension(NLAYER) :: SINGLE_TAULS
          real, dimension(NLAYER) :: SINGLE_TEMPS

          real :: NU

          NU = 1e13
        
          DO J = 1, NLAYER
            SINGLE_W0S(J) = W0(1,J)
            SINGLE_G0S(J) = G0(1,J)

            SINGLE_TAULS(J) = TAUL(1,J)
            SINGLE_TEMPS(J) = TT(J)
          END DO

          ! I also created a variable for NUM_NUS
          CALL GET_IR_INTENSITY(NUS(1), EMIS(1), RSFX(1), NLAYER, SINGLE_TEMPS, SINGLE_TAULS, SINGLE_W0S, SINGLE_G0S)
          CALL GET_SOLAR_INTENSITY(NUS(2), EMIS(1), RSFX(1), NLAYER, SINGLE_TAULS, SINGLE_W0S, SINGLE_G0S)

      END SUBROUTINE TWO_STREAM_WRAPPER
      
      
      SUBROUTINE GET_IR_INTENSITY(NU, EMIS, RSFX, NLAYER, TEMPS, TAULS, W0, G0)
          integer :: NLAYER, J, L, Z
          real :: NU, mu_1, mu, temp_val_2, BB_TOP_OF_ATM, BB_BOTTOM_OF_ATM, SFCS_HEMISPHERIC
          real(16) :: temp_val_1
          real :: e_con, pi, bolz_constant, clight, h_constant

          real, dimension(50) :: TAULS, W0, G0, TEMPS

          real ::  EMIS, RSFX
          real, dimension(NLAYER) :: y1, y2
          real, dimension(NLAYER) :: LAMBDAS, GAMMA
          real, dimension(NLAYER) :: temp_e_val
          real, dimension(NLAYER) :: e1, e2, e3, e4
          real, dimension(NLAYER) :: temp_gamma_val
          real, dimension(NLAYER) :: CP, CPB, CM, CMB

          real, dimension(2*NLAYER) :: A, B, D, E
          real, dimension(2*NLAYER) :: AS, DS, X, Y

          real, dimension(NLAYER) :: B0, B1
          real, dimension(NLAYER) :: SOURCE_G, SOURCE_H, SOURCE_J, SOURCE_K
          real, dimension(NLAYER) :: ALPHA_1, ALPHA_2, SIGMA_1, SIGMA_2
          real, dimension(NLAYER) :: SOURCE_Y1, SOURCE_Y2, source_temp
          real, dimension(NLAYER) :: HEMISPHERIC_INTENSITY_DOWN, HEMISPHERIC_INTENSITY_UP

          mu   = 1.0

          bolz_constant = 1.380649e-23
          h_constant    = 6.62607015e-34
          clight        = 3e8
          pi            = 3.141592653589
          e_con         = 2.718281828459

          temp_val_1 = 2.0 * NU * NU / (CLIGHT * CLIGHT)
          temp_val_1 = temp_val_1 * (h_constant)
          temp_val_1 = temp_val_1 * NU
          temp_val_2 = e_con ** (h_constant * NU / (bolz_constant * TEMPS(NLAYER))) - 1.0
          BB_TOP_OF_ATM = REAL(temp_val_1 * (1.0 / temp_val_2))

          temp_val_1 = 2.0 * NU * NU / (CLIGHT * CLIGHT)
          temp_val_1 = temp_val_1 * (h_constant)
          temp_val_1 = temp_val_1 * NU
          temp_val_2 = e_con ** (h_constant * NU / (bolz_constant * TEMPS(NLAYER))) - 1.0
          BB_BOTTOM_OF_ATM = REAL(temp_val_1 * (1.0 / temp_val_2))

          SFCS_HEMISPHERIC = EMIS * PI * BB_BOTTOM_OF_ATM
          mu_1 = 0.5

          DO J = 1, NLAYER
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
          DO L = 2, 2*NLAYER-1, 2
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

          A(2*NLAYER) = e1(NLAYER) - RSFX * e3(NLAYER);
          B(2*NLAYER) = e2(NLAYER) - RSFX * e4(NLAYER);
          D(2*NLAYER) = 0.0;

          ! This is the part of the code that solves for the blackbody stuff
          DO J=1, NLAYER
              IF (1 .GE. J-1) THEN
                  KINDEX = 1
              ELSE
                  KINDEX = J-1
              END IF

              ! STUPID HACK
              ! PREVENTS NU ** 3 from being too big
              ! OMG FORTRAN
              temp_val_1 = 2.0 * NU * NU / (CLIGHT * CLIGHT)
              temp_val_1 = temp_val_1 * (h_constant)
              temp_val_1 = temp_val_1 * NU
              temp_val_2 = e_con ** (h_constant * NU / (bolz_constant * TEMPS(J))) - 1.0

              B0(J) = REAL(temp_val_1 * (1.0 / temp_val_2))

              ! This is if you want double gray, you'll also have to fix boundary stuff
              ! B0(J) = TEMPS(J) * TEMPS(J) * TEMPS(J) * TEMPS(J) * 1.8049459031e-8

              B1(J) = (B0(J) - B0(KINDEX)) / TAULS(J)
          END DO

          DO J=1, NLAYER
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
          DO L = 2, 2*NLAYER-1, 2
              E(L)   = (CP(J+1) - CPB(J)) * e2(J+1) + (CMB(J) - CM(J+1))*e4(J+1)
              E(L+1) = e3(J) * (CP(J+1) - CPB(J)) + e1(J) * (CMB(J) - CM(J+1))
              J = J + 1
          END DO

          ! HERE ARE THE TOP AND BOTTOM BOUNDARY CONDITIONS AS WELL AS THE
          ! BEGINNING OF THE TRIDIAGONAL SOLUTION. I ASSUME NO
          ! DIFFUSE RADIATION IS INCIDENT AT THE TOP.
          E(1) = 0.0 - CM(1);
          E(2*NLAYER)  = SFCS_HEMISPHERIC - CPB(NLAYER) + RSFX * CMB(NLAYER)

          DS(2*NLAYER) = E(2*NLAYER) / B(2*NLAYER)
          AS(2*NLAYER) = A(2*NLAYER) / B(2*NLAYER)

          DO L = 1, 2*NLAYER - 1
              X(2*NLAYER-L)  = 1.0 / (B(2*NLAYER-L) - D(2*NLAYER-L) * AS(2*NLAYER-L+1))
              AS(2*NLAYER-L) = A(2*NLAYER-L) * X(2*NLAYER-L)
              DS(2*NLAYER-L) = (E(2*NLAYER-L) - (D(2*NLAYER-L) * DS(2*NLAYER-L+1))) * X(2*NLAYER-L)
          END DO

          Y(1) = DS(1)

          DO L = 2, 2*NLAYER
              Y(L) = DS(L) - AS(L) * Y(L-1)
          END DO

          DO J = 2, NLAYER+1
              SOURCE_Y1(J-1) = Y(2*J-3)
              SOURCE_Y2(J-1) = Y(2*J-2)
          END DO

          DO J=1, NLAYER
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

          HEMISPHERIC_INTENSITY_DOWN(1) = BB_TOP_OF_ATM * e_con ** (-TAULS(1)/mu) + 
     &                    SOURCE_J(1)/(LAMBDAS(1)*mu + 1.0) * (1.0 - e_con ** (-TAULS(1)*(LAMBDAS(1)+1.0/mu))) + 
     &                    SOURCE_K(1)/(LAMBDAS(1)*mu - 1.0) * (e_con ** (-TAULS(1)/mu) - e_con ** (-TAULS(1)*LAMBDAS(1)))+ 
     &                    SIGMA_1(1) * (1.0 - e_con ** (-TAULS(1)/mu)) + 
     &                    SIGMA_2(1) * (mu * e_con ** (-TAULS(1)/mu) + TAULS(1) - mu);

          DO J = 2, NLAYER
              HEMISPHERIC_INTENSITY_DOWN(J) = HEMISPHERIC_INTENSITY_DOWN(J-1) * e_con ** (-TAULS(J)/mu) + 
     &                     SOURCE_J(J)/(LAMBDAS(J)*mu + 1.0) * (1.0 - e_con ** (-TAULS(J)*(LAMBDAS(J)+1.0/mu))) + 
     &                     SOURCE_K(J)/(LAMBDAS(J)*mu - 1.0) * (e_con ** (-TAULS(J)/mu) - e_con ** (-TAULS(J)*LAMBDAS(J)))+ 
     &                     SIGMA_1(J) * (1.0 - e_con ** (-TAULS(J)/mu)) + 
     &                     SIGMA_2(J) * (mu * e_con ** (-TAULS(J)/mu) + TAULS(J) - mu)
          END DO

          HEMISPHERIC_INTENSITY_UP(NLAYER) = 2.0 * BB_BOTTOM_OF_ATM * EMIS * pi

          !Calculate the upward intensity next
          DO Z = 1, NLAYER - 1
              J = NLAYER - Z
              HEMISPHERIC_INTENSITY_UP(J) = HEMISPHERIC_INTENSITY_UP(J+1) * e_con ** (-TAULS(J+1)/mu) + 
     &                         SOURCE_G(J)/(LAMBDAS(J)*mu-1.0)*(e_con ** (-TAULS(J+1)/mu)-e_con ** (-TAULS(J+1)*(LAMBDAS(J)))) + 
     &                         SOURCE_H(J)/(LAMBDAS(J)*mu+1.0) * (1.0 - e_con ** (-TAULS(J+1) * (LAMBDAS(J) + 1.0/mu))) + 
     &                         ALPHA_1(J) * (1.0 - e_con ** (-TAULS(J+1)/mu)) + 
     &                         ALPHA_2(J) * (mu - ((TAULS(J+1) + mu) * (e_con ** (-TAULS(J+1)/mu))))
          END DO

      END SUBROUTINE GET_IR_INTENSITY



      SUBROUTINE GET_SOLAR_INTENSITY(NU, EMIS, RSFX, NLAYER, TAULS, W0, G0)
          integer :: NLAYER, J, L
          real :: RSFX, FLUX_SURFACE_QUADRATURE, mu_0, mu_1

          real :: e_con, pi

          real, dimension(NLAYER) :: QUADRATURE_FLUX, QUADRATURE_INTENSITY, TAULS, W0, G0
          real, dimension(NLAYER+1) :: TAUCS

          real, dimension(NLAYER) :: DIRECT_QUADRATURE
          real, dimension(NLAYER) :: y1, y2, y3, y4
          real, dimension(NLAYER) :: LAMBDAS, GAMMA
          real, dimension(NLAYER) :: temp_e_val
          real, dimension(NLAYER) :: e1, e2, e3, e4
          real, dimension(NLAYER) :: temp_gamma_val
          real, dimension(NLAYER) :: CP, CPB, CM, CMB

          real, dimension(2*NLAYER) :: A, B, D, E
          real, dimension(2*NLAYER) :: AS, DS, X, Y

          e_con = 2.718281828459
          pi    = 3.141592653589

          mu_0 = 1.0
          FLUX_SURFACE_QUADRATURE = 1.0
          

          TAUCS(1) = 0.0
          DO N = 1, NLAYER
              TAUCS(N+1) = TAUCS(N) + TAULS(N)
          END DO

          DO J = 1, NLAYER
              DIRECT_QUADRATURE(J) = mu_0 * pi * FLUX_SURFACE_QUADRATURE*e_con ** (-1.0 * (TAUCS(J) + TAULS(J)) / mu_0)
          END DO

          ! Boundary conditions
          mu_1 = 0.577350;
          SFCS_QUADRATURE = RSFX * mu_0 * e_con ** (-(TAUCS(NLAYER + 1))/mu_0) * pi * FLUX_SURFACE_QUADRATURE

          ! HERE WE FIND LAYER PROPERTIES FOLLOWING GENERAL SCHEME
          ! OF MEADOR AND WEAVOR. THEN WE SET UP LAYER PROPERTIES
          ! NEEDED FOR MATRIX

          DO J = 1, NLAYER
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
          DO L = 2, 2*NLAYER-1, 2
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

          A(2*NLAYER) = e1(NLAYER) - RSFX * e3(NLAYER);
          B(2*NLAYER) = e2(NLAYER) - RSFX * e4(NLAYER);
          D(2*NLAYER) = 0.0;

          DO J = 1, NLAYER
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
          DO L = 2, 2*NLAYER-1, 2
              E(L)   = (CP(J+1) - CPB(J)) * e2(J+1) + (CMB(J) - CM(J+1))*e4(J+1)
              E(L+1) = e3(J) * (CP(J+1) - CPB(J)) + e1(J) * (CMB(J) - CM(J+1))
              J = J + 1
          END DO

          ! HERE ARE THE TOP AND BOTTOM BOUNDARY CONDITIONS AS WELL AS THE
          ! BEGINNING OF THE TRIDIAGONAL SOLUTION. I ASSUME NO
          ! DIFFUSE RADIATION IS INCIDENT AT THE TOP.
          E(1) = 0.0 - CM(1);
          E(2*NLAYER)  = SFCS_QUADRATURE - CPB(NLAYER) + RSFX * CMB(NLAYER)

          DS(2*NLAYER) = E(2*NLAYER) / B(2*NLAYER)
          AS(2*NLAYER) = A(2*NLAYER) / B(2*NLAYER)

          DO L = 1, 2*NLAYER - 1
              X(2*NLAYER-L)  = 1.0 / (B(2*NLAYER-L) - D(2*NLAYER-L) * AS(2*NLAYER-L+1))
              AS(2*NLAYER-L) = A(2*NLAYER-L) * X(2*NLAYER-L)
              DS(2*NLAYER-L) = (E(2*NLAYER-L) - (D(2*NLAYER-L) * DS(2*NLAYER-L+1))) * X(2*NLAYER-L)
          END DO

          Y(1) = DS(1)

          DO L = 2, 2*NLAYER
              Y(L) = DS(L) - AS(L) * Y(L-1)
          END DO


          DO J = 1, NLAYER
              QUADRATURE_FLUX(J) = Y(2*J-1)*(e1(J)-e3(J)) + Y(2*J)*(e2(J)-e4(J)) + CPB(J) - CMB(J) - DIRECT_QUADRATURE(J)

              QUADRATURE_INTENSITY(J) = (1.0 / mu_1) * (Y(2*J-1)*(e1(J) + e3(J)) + 
     &                   (Y(2*J) * (e2(J) + e4(J))) + CPB(J) + CMB(J))+DIRECT_QUADRATURE(J)/mu_0
          END DO
      RETURN
      END SUBROUTINE GET_SOLAR_INTENSITY