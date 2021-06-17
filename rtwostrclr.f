      SUBROUTINE TWOSTRCLR
!
!    ******************************************************************
!    *  Purpose             :  Defines matrix properties and sets up  *
!    *                         matrix coefficients that do not depend *
!    *                         on zenith angle or temperature.        *
!    *  Subroutines Called  :  None                                   *
!    *  Input               :  W0, G0                                 *
!    *  Output              :  B1, B2, GAMI, ACON, EL1, AF, ETC       *
!    * ****************************************************************
!
      include 'rcommons.h'
       DO 10 L    =  LLS,LLA
          if( L .LE. NSOLP )then
            U1I(L) = SQ3  !2.d0 !SQ3
          else
            U1I(L) = 2.d0
          endif
  10      U1S(L)  =  TPI/U1I(L)
!       write(*,*) 'U1S',U1S
!       write(*,*) 'U1I',U1I 
!       write(*,*) 'W0',W0
!       write(*,*) 'U1I',U1I
!       write(*,*)'U1S',U1S  
!       stop

!       DO J = 1, NLAYER
!        DO L= 1,2 
!        W0(L,J) = 0.9
!        G0(L,J) = 0.848
!        ENDDO
!       ENDDO
!      HERE WE DEFINE LAYER PROPERTIES FOLLOWING GENERAL SCHEME
!      OF MEADOR AND WEAVOR. THEN WE SET UP LAYER PROPERTIES
!      NEEDED FOR MATRIX.
!
       DO 14 J          =  1,NLAYER
          DO 14 L       =  LLS,LLA
!            THESE ARE FOR TWO STREAM AND HEMISPHERIC MEANS
             B1CLR(L,J)    =  0.5*U1I(L)*(2. - W0CLR(L,J)*(1. + G0CLR(L,J)))
             B2CLR(L,J)    =  0.5*U1I(L)*W0CLR(L,J)*(1. - G0CLR(L,J))
             AKCLR(L,J)    = SQRT(ABS(B1CLR(L,J)*B1CLR(L,J) - B2CLR(L,J)*B2CLR(L,J)))
             GAMICLR(L,J)  =  B2CLR(L,J)/(B1CLR(L,J) + AKCLR(L,J))
             EE1CLR(L,J)   =  EXP(-AKCLR(L,J)*TAULCLR(L,J))
!             write(*,*) 'EE1',J,EE1(L,J)
!             write(*,*) 'GAMI',J,GAMI(L,J)
             EL1CLR(L,J)   =  1.0 + GAMICLR(L,J) *EE1CLR(L,J)  !e1                          
             EM1CLR(L,J)   =  1.0 - GAMICLR(L,J) * EE1CLR(L,J) !e2                      
             EL2CLR(L,J)   =  GAMICLR(L,J) + EE1CLR(L,J)       !e3                        
             EM2CLR(L,J)   =  GAMICLR(L,J) - EE1CLR(L,J)       !e4

  14  CONTINUE
!
!     WE SEEK TO SOLVE AX(L-1)+BX(L)+EX(L+1) = D.
!     L=2N FOR EVEN L, L=N+1 FOR ODD L. THE MEAN INTENSITY (TMI/4PI)
!     AND THE NET FLUX (FNET) ARE RELATED TO X'S AS NOTED IN ADD.
!     FIRST WE SET UP THE COEFFICIENTS THAT ARE INDEPENDENT OF SOLAR
!     ANGLE OR TEMPARATURE: A(I),B(I),E(I). D(I) IS DEFINED IN ADD.
!
      J                 =  0
      DO 18 JD          =  2,JN,2
         J              =  J + 1
         DO 18 L        =  LLS,LLA
!          HERE ARE THE EVEN MATRIX ELEMENTS
             AFCLR(L,JD)   =  EM1CLR(L,J+1)*EL1CLR(L,J)-
     &                         EM2CLR(L,J+1)*EL2CLR(L,J)
             BFCLR(L,JD)   =  EM1CLR(L,J+1)* EM1CLR(L,J)-
     &                         EM2CLR(L,J+1)*EM2CLR(L,J)
             EFCLR(L,JD)   =  EL1CLR(L,J+1)*EM2CLR(L,J+1) - 
     &                         EL2CLR(L,J+1)*EM1CLR(L,J+1)
!          HERE ARE THE ODD MATRIX ELEMENTS EXCEPT FOR THE TOP. 
             AFCLR(L,JD+1) =  EM1CLR(L,J)*EL2CLR(L,J)-
     &                         EL1CLR(L,J)*EM2CLR(L,J)
             BFCLR(L,JD+1) =  EL1CLR(L,J+1)*EL1CLR(L,J)-
     &                         EL2CLR(L,J+1)*EL2CLR(L,J)      
             EFCLR(L,JD+1) =  EL2CLR(L,J)*EM2CLR(L,J+1)-
     &                         EL1CLR(L,J)*EM1CLR(L,J+1)
  18  CONTINUE
!
!     HERE ARE THE TOP AND BOTTOM BOUNDARY CONDITIONS AS WELL AS THE
!     BEGINNING OF THE TRIDIAGONAL SOLUTION DEFINITIONS. I ASSUME
!     NO DIFFUSE RADIATION IS INCIDENT AT UPPER BOUNDARY.
!
      DO 20 L        = LLS,LLA
         AFCLR(L,1)     = 0.0
         BFCLR(L,1)     = EL1CLR(L,1)
         EFCLR(L,1)     = -EM1CLR(L,1)
         AFCLR(L,JDBLE) = EL1CLR(L,NLAYER)-RSFX(L)*EL2CLR(L,NLAYER)
         BFCLR(L,JDBLE) = EM1CLR(L,NLAYER)-RSFX(L)*EM2CLR(L,NLAYER)
         EFCLR(L,JDBLE) = 0.0
  20  CONTINUE
!          write(*,*)'RSFX in twostream',RSFX
!         write(*,*)'AF',AF
!         write(*,*)'BF',BF
!         write(*,*)'EF',EF
!      write(*,*) 'w0',w0
!      write(*,*) 'g0',w0

      RETURN
      END

