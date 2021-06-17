      SUBROUTINE OPPRCLR
!
!     **************************************************************
!     *  Purpose             :  CaLculates optical properties      *
!     *                         such as single scattering albedo,  *
!     *                         asymmetry parameter, etc.          *
!     *                         This routine is case dependent and *
!     *                         wiLL have to be repLaced by the    *
!     *                         user.                              *
!     *  Subroutines Called  :  None                               *
!     *  Input               :  PAH2O, RDH2O, CO2, O3, ETC         *
!     *  Output              :  TAUL, W0, G0, OPD, Y3              *
!     * ************************************************************
!
      include 'rcommons.h'

!     W0(NWAVE,NLAYER) : SINGLE SCATTERING ALBEDO *** delta scaled ***
!     G0(NWAVE,NLAYER) : ASYMMETRY PARAMETER *** delta scaled ***
!     OPD(NWAVE,NLAYER): cumulative OPTICAL DEPTH *** delta scaled ***
!     SFL(NWAVE)       : SOLAR FLUX
!    uW0(NWAVE,NLAYER)  : unscaled SINGLE SCATTERING ALBEDO 
!    uG0(NWAVE,NLAYER)  : unscaled ASYMMETRY PARAMETER 
!    uTAUL(NWAVE,NLAYER): unscaled OPTICAL DEPTH of layer
!
!     ASSUME THAT P IS SAME ON ALL SIGMA LEVELS. IF PSURFACE
!     VARIES A LOT, THEN WE MAY NEED TO CALCULATE TAUO3,
!     TAUCO2, TAUO2 FOR EACH POINT.
!
!     NOTE : THE TOP LAYER IS A DUMMY. IT CONTAINS A DEFINED GAS
!            AMOUNT. DIFFERENT MODELS WILL REQUIRE DIFFERENT
!            TREATMENT OF THIS.
!     CALCULATE TOTAL OPTICAL DEPTH INCLUDING GASES. THEN
!     GIVEN THE AEROSOL OPTICAL DEPTHS AND CLOUD OPTICAL DEPTHS,
!     CALCULATE FINAL OPTICAL PROPERTIES. WE USE A DELTA
!     TWO STREAM APPROACH TO FIND W0, SINGLE SCATTERING ALBEDO,
!     G0, ASYMMMETRY PARAMETER, TAUL, LAYER OPTICAL DEPTH,
!     OPD, CUMULATIVE OPTICAL DEPTH TO BASE OF LAYER.
!     open( 23, file='tau.dat', form='formatted', status='unknown' )   
!     write(23,*) nlayer, ntotal
!
!   NOTE: THIS IS DOUBLE GRAY SPECIFIC; WOULD REQUIRE GENERALIZATION
         
!  110   CONTINUE
      DO 200  J          =   1,NLAYER
!       DO 130 L          = 1,NTOTAL
!    VIS
       TAUAER(1,J)       = 0.0
       WOL(1,J)          = PI0AERSW 
       GOL(1,J)          = ASYMSW
       TAUS(1,J)         = PI0AERSW*0.0
       
!           IF (TAUS(2,J) .NE. 0) THEN 
!            TTAS         = TAUS(1,J)
!           ELSE
!            TTAS         = 1.
!           ENDIF

!       GOL(1,J)          = G01(1,J)/TTAS

!    IR
       TAUAER(2,J)       = 0.0
       WOL(2,J)          = PI0AERLW
       GOL(2,J)          = ASYMLW
       TAUS(2,J)         = 0.0

!         if (TAUS(2,J).ne.0) then
!            TTAS         = TAUS(2,J)
!         else
!            TTAS         = 1.
!         ENDIF

!       GOL(2,J)          = G01(1,J)/TTAS

 200  CONTINUE


!            do j  = 1,nlayer
!            write(*,*)'j,Tauray/Taugas',j,TAURAY(1,J)/TAUGAS(1,J)
!            enddo

!     iradgas = 0: no gas in radiative xfer!
      iradgas = 1     

      DO 500 J           = 1,NLAYER
          j1             = max( 1, j-1 )
!
          DO 400 L       = LLS,LLA
 
!MTR              TAUL(L,J)   = TAUH2O(L,J)+TAUGAS(L,J)+   &
!MTR                               PARAY(L,J)+TAUAER(L,J)+TAUCLD(L,J)
!              THE CODE ORIGINALLY INCLUDED TAUCLD AS WELL AS TAUAER. I
!              DO NOT SEE THE POINT OF HAVING BOTH OPTICAL DEPTHS IN OUR
!              MODEL SO I AM OMMITTING TAUCLD, W0CLD, GCLD, ETC

              TAUL(L,J)   = TAUGAS(L,J)+TAURAY(L,J)+TAUAER(L,J)!+TAUCLD(L,J)
             
             if (iradgas.eq.0) then
             tauL(L,j) = tauaer(L,j)
             endif

             if( TAUL(L,J) .lt. EPSILON ) then
             TAUL(L,J) = EPSILON
             endif

             utauL(L,j)  = TAUL(L,J)
             WOT         = (TAURAY(L,J)+TAUAER(L,J)*WOL(L,J))/TAUL(L,J)
             if (iradgas.eq.0) then
              wot = woL(L,j)
             endif
 
             WOT         = min(1.-EPSILON,WOT)
             uw0(L,j)    = WOT
!             write(*,*) 'WOT',WOT
             DENOM       = (TAURAY(L,J)+ TAUAER(L,J)*WOL(L,J))
             if( DENOM .LE. EPSILON ) then 
             DENOM = EPSILON then
             endif
             if( DENOM .GT. EPSILON ) then
               GOT = ( GOL(L,J)* WOL(L,J)*TAUAER(L,J) ) / DENOM
!       print*, j, L, GCLD(l,j), gol(l,j),wol(l,j),got
             else
               GOT = 0.
             endif
             if (iradgas.eq.0) then
             GOT = goL(L,j)
             endif
             ug0(L,j)    = GOT
             uOPD(L,J)   = 0.0
             uOPD(L,J)   = uOPD(L,J1)+uTAUL(L,J)
             IF (deltascale) THEN
             FO          = GOT*GOT
             DEN         = 1.-WOT*FO
             TAUL(L,J)   = TAUL(L,J) * DEN
             W0(L,J)     = (1.-FO)*WOT/DEN
             G0(L,J)     = GOT/(1.+GOT)
             OPD(L,J)    = 0.0
             OPD(L,J)    = OPD(L,J1)+TAUL(L,J)             
             ELSE 
                  W0(L,J)= uw0(L,J)
                  G0(L,J)= ug0(L,J)
                TAUL(L,J)= utaul(L,J)
                 OPD(L,J)= uOPD(L,J)
             ENDIF 
!!!!!!!!!!!!!!!!HERE'S WHERE YOU CAN HARDWIRE VALUES!!!!!!!!!
 
             if( taul(L,j).lt.0. ) then
       write(*,*) 'taul lt 0'
       stop
       endif
                       
 400        CONTINUE
           
          IF(IR .EQ. 1) THEN
             DO 450 I        =   1,NGAUSS
                DO 425 L     =   LLS,LLA
                   Y3(L,I,J) =   EXP(-TAUL(L,J)/GANGLE(I))
 425            CONTINUE
 450         CONTINUE
          ENDIF
 500  CONTINUE

!      DO J = 1,NLAYER 
!      write(*,*) 'W0',J,W0(1,J)
!      ENDDO 
!      DO J = 1,NLAYER
!      write(*,*) 'G0',J,G0(1,J)
!      ENDDO
!      DO J = 1,NLAYER
!      write(*,*) 'WOL',J,WOL(1,J)
!      ENDDO
!      DO J = 1,NLAYER
!      write(*,*) 'GOL',J,GOL(1,J)
!      ENDDO
!      DO J = 1,NLAYER
!      write(*,*) 'TAURAY',J,TAURAY(1,J)
!      ENDDO
!      DO J = 1,NLAYER
!      write(*,*) 'TAUAER',J,TAUAER(2,J)
!      ENDDO
!      DO J = 1,NLAYER
!      write(*,*) 'TAUGAS',J,TAUGAS(1,J)
!      ENDDO
!      DO J = 1,NLAYER
!      write(*,*) 'TAUL',J,TAUL(1,J)
!      ENDDO
!      DO J = 1,NLAYER
!      write(*,*) 'OPD',J,OPD(1,J) 
!      ENDDO     
!         stop
      RETURN
      END

