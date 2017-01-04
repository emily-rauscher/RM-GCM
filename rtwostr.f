      SUBROUTINE TWOSTR
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
!
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
!       write(*,*) 'G0',G0
!       write(*,*) 'U1I',U1I
!       write(*,*)'U1S',U1S  
!    
!      HERE WE DEFINE LAYER PROPERTIES FOLLOWING GENERAL SCHEME
!      OF MEADOR AND WEAVOR. THEN WE SET UP LAYER PROPERTIES
!      NEEDED FOR MATRIX.
!
       DO 14 J          =  1,NLAYER
          DO 14 L       =  LLS,LLA
!            THESE ARE FOR TWO STREAM AND HEMISPHERIC MEANS
             B1(L,J)    =  0.5*U1I(L)*(2. - W0(L,J)*(1. + G0(L,J)))
             B2(L,J)    =  0.5*U1I(L)*W0(L,J)*(1. - G0(L,J))
             AK(L,J)    = SQRT(ABS(B1(L,J)*B1(L,J) - B2(L,J)*B2(L,J)))
             GAMI(L,J)  =  B2(L,J)/(B1(L,J) + AK(L,J))
             EE1(L,J)   =  EXP(-AK(L,J)*TAUL(L,J))
!             write(*,*) 'EE1',J,EE1(L,J)
!             write(*,*) 'GAMI',J,GAMI(L,J)
             EL1(L,J)   =  1.0 + GAMI(L,J) *EE1(L,J)  !e1                          
             EM1(L,J)   =  1.0 - GAMI(L,J) * EE1(L,J) !e2                      
             EL2(L,J)   =  GAMI(L,J) + EE1(L,J)       !e3                        
             EM2(L,J)   =  GAMI(L,J) - EE1(L,J)       !e4

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
             AF(L,JD)   =  EM1(L,J+1)*EL1(L,J)-EM2(L,J+1)*EL2(L,J)
             BF(L,JD)   =  EM1(L,J+1)* EM1(L,J)-EM2(L,J+1)*EM2(L,J)
             EF(L,JD)  = EL1(L,J+1)*EM2(L,J+1) - EL2(L,J+1)*EM1(L,J+1)
!          HERE ARE THE ODD MATRIX ELEMENTS EXCEPT FOR THE TOP. 
             AF(L,JD+1) =  EM1(L,J)*EL2(L,J)-EL1(L,J)*EM2(L,J)
             BF(L,JD+1) =  EL1(L,J+1)*EL1(L,J) - EL2(L,J+1)*EL2(L,J)      
             EF(L,JD+1) =  EL2(L,J)*EM2(L,J+1)-EL1(L,J)*EM1(L,J+1)
  18  CONTINUE
!
!     HERE ARE THE TOP AND BOTTOM BOUNDARY CONDITIONS AS WELL AS THE
!     BEGINNING OF THE TRIDIAGONAL SOLUTION DEFINITIONS. I ASSUME
!     NO DIFFUSE RADIATION IS INCIDENT AT UPPER BOUNDARY.
!
      DO 20 L        = LLS,LLA
         AF(L,1)     = 0.0
         BF(L,1)     = EL1(L,1)
         EF(L,1)     = -EM1(L,1)
         AF(L,JDBLE) = EL1(L,NLAYER)-RSFX(L)*EL2(L,NLAYER)
         BF(L,JDBLE) = EM1(L,NLAYER)-RSFX(L)*EM2(L,NLAYER)
         EF(L,JDBLE) = 0.0
  20  CONTINUE
!          write(*,*)'RSFX in twostream',RSFX
!         write(*,*)'AF',AF
!         write(*,*)'BF',BF
!         write(*,*)'EF',EF
      RETURN
      END

