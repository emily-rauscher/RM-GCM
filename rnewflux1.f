      SUBROUTINE NEWFLUX1
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
!
!     LOCAL DIMENSIONS
      DIMENSION Y1(NTOTAL,NGAUSS,NLAYER),     
     &           Y2(NTOTAL,NGAUSS,NLAYER),     
     &           Y4(NTOTAL,NGAUSS,NLAYER),     
     &           Y8(NTOTAL,NGAUSS,NLAYER),     
     &           A1(NTOTAL,NLAYER),A2(NTOTAL,NLAYER),     
     &           A3(NTOTAL,NLAYER),A4(NTOTAL,NLAYER),     
     &           A7(NTOTAL,NLAYER),Y5(NTOTAL,NLAYER)
!

      DO 200 J           =  1,NLAYER
          kindex         = max( 1, j-1 )
          DO 100  L      =  NSOLP+1,NTOTAL
!            HERE WE DO NO SCATTERING COEFFICIENTS
             A3(L,J)     =  PTEMP(L,KINDEX)*TPI
             A4(L,J)     =  TPI*SLOPE(L,J)
             A7(L,J)     =  A3(L,J)
             Y5(L,J)     =  A4(L,J)*TAUL(L,J)
 100      CONTINUE
!
!         HERE WE DO SCATTERING
          IF(IRS .NE. 0) THEN
              DO 50 L    =  NSOLP+1,NTOTAL
                X4       =  SLOPE(L,J)*(TPI*B3(L,J)-U1S(L))
                A1(L,J)  =  U1I(L) - AK(L,J)
                A2(L,J)  =  GAMI(L,J)*(AK(L,J)+U1I(L))
                A3(L,J)  =  A3(L,J)+X4
                A7(L,J)  =  A7(L,J)-X4
 50           CONTINUE
          ENDIF
  200 CONTINUE
!
!     CALCULATIONS FOR ALL GAUSS POINTS. HERE WE DO NO SCATTERING COEFFI
!
      DO 400       J         =  1,NLAYER
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
     &                             (AK(L,J)*GANGLE(I)-1.)            
                 YB        =  A2(L,J)*(1.- EE1(L,J)*Y3(L,I,J))/
     &                             (AK(L,J)*GANGLE(I)+1.)            
                 CKP= CK1(L,J)+CK2(L,J)
                 CKM= CK1(L,J) -CK2(L,J)
                 Y1(L,I,J) =  CKP*YB+CKM*YA                    
                 Y2(L,I,J) =  CKP*YA+ CKM*YB                   
 325          CONTINUE
            ENDIF
 350     CONTINUE
 400  CONTINUE
!
      DO 450 J             =  1,NLAYER
         DO 425  L         =  NSOLP+1,NTOTAL
            TMID(L,J) = 0.0
            TMIU(L,J) = 0.0
            DIREC(L,J)     =  0.0
            DIRECTU(L,J)   =  0.0
 425     CONTINUE
 450  CONTINUE
!
!     DIREC IS DOWNWARD FLUX. DIRECTU IS UPWARD FLUX.
!     CALCULATE DINTENT THE DOWNWARD INTENSITY AND DIREC THE DOWNWARD FL
!
       DO 500 I             = 1,NGAUSS
          DO 475 L          = NSOLP+1,NTOTAL
             if( iblackbody_above .eq. 1 )then
               DINTENT(L,I,1) = PTEMPT(L)*Y3(L,I,1)*TPI     
     &                           +Y1(L,I,1)+     
     &                           (1.-Y3(L,I,1))*Y4(L,I,1)
             else
               DINTENT(L,I,1) = (1.-Y3(L,I,1))*Y4(L,I,1)     
     &                           +Y1(L,I,1)
             endif
             TMID(L,1)      = TMID(L,1)+DINTENT(L,I,1)*GRATIO(I)
             DIREC(L,1)     = DIREC(L,1)+DINTENT(L,I,1)*
     &                         GWEIGHT(I)
 475      CONTINUE
 500   CONTINUE
!
!      DINTENT IS DOWNWARD INTENSITY * TPI. DIREC IS THE DOWNWARD FLUX.
!
       DO 530        J           = 2,NLAYER
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
             UINTENT(L,I,NLAYER)  =  PTEMPG(L)*EMIS(L)     
     &                               *TPI+2.*RSFX(L)*DIREC(L,NLAYER)
             TMIU(L,NLAYER)       =  TMIU(L,NLAYER)+     
     &                               UINTENT(L,I,NLAYER)*GRATIO(I)
             DIRECTU(L,NLAYER)    =  DIRECTU(L,NLAYER)+     
     &                               UINTENT(L,I,NLAYER)*GWEIGHT(I)
 560      CONTINUE
 570   CONTINUE
!
      DO 650        M              = 2,NLAYER
          J                        = NLAYER-M+1
          DO 640    I              = 1,NGAUSS
             DO 630 L              = NSOLP+1,NTOTAL
                 UINTENT(L,I,J)    = (UINTENT(L,I,J+1)-Y5(L,J+1))     
     &                               *Y3(L,I,J+1)+Y2(L,I,J+1)+     
     &                               (1.-Y3(L,I,J+1))*Y8(L,I,J+1)
                  TMIU(L,J)        = TMIU(L,J)+UINTENT(L,I,J)*GRATIO(I)
                  DIRECTU(L,J)     = DIRECTU(L,J) +     
     &                               GWEIGHT(I)*UINTENT(L,I,J)
 630         CONTINUE
 640      CONTINUE
 650  CONTINUE

      RETURN
      END

