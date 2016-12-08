      SUBROUTINE OPPR
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
!
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
!i
!      DO 120 J 
!      DO 110 L           =   1,NWAVE
!...ej hack: ssa properties calculated in setuprad.f
!          TAUA(1,16)       =   1.0
!          TAUA(2,16)       =   .5
!          TAUS(1,16)       =   .98
!          TAUS(2,16)       =   0.4
!          G01(1,16)        =   0.7
!          G01(2,16)        =  0.99
!          TAUA(1,17)       =   2.0
!          TAUA(2,17)       =   1.0
!          TAUS(1,17)       =   1.96
!          TAUS(2,17)       =   0.8
!          G01(1,17)        =   0.7
!          G01(2,17)        =  0.99
!          TAUA(1,18)       =   4.0
!          TAUA(2,18)       =   2.0  
!          TAUS(1,18)       =   3.92
!          TAUS(2,18)       =   1.6
!          G01(1,18)        =   0.7
!          G01(2,18)        =  0.99
!          TAUA(1,19)       =   8.0
!         TAUA(2,19)       =   4.0 
!          TAUS(1,19)       =   7.84
!          TAUS(2,19)       =   3.2
!          G01(1,19)        =   0.7
!          G01(2,19)        =  0.99
!          TAUA(1,20)       =   10.0
!          TAUA(2,20)       =   5.0   
!          TAUS(1,20)       =   9.8
!          TAUS(2,20)       =   4.0
!          G01(1,20)        =   0.7
!          G01(2,20)        =  0.99




!  110   CONTINUE


      DO 200  J          =   1,NLAYER
!
!MTR      DO 110 L           =   1,NWAVE
!...ej hack: ssa properties calculated in setuprad.f
!          TAUA(L,J)       =   0.0
!          TAUS(L,J)       =   0.0
!          G01(L,J)        =   0.0
!  110   CONTINUE
!
!       do ig            =   1,NGROUP
!         DO 120 I       =   1,NRAD
!             R2Z        =   XSECTA(I,ig)*CAER(I,J,ig)
!          DO 115 L      =   1,NWAVE
!...ej hack: ssa properties calculated in setuprad.f
!             TAUA(L,J)  =   TAUA(L,J)+RDQEXT(I,ig,L)*R2Z
!             TAUS(L,J)  =   TAUS(L,J)+QSCAT(I,ig,L)*R2Z
!             G01(L,J)   =   G01(L,J) +QBRQS(I,ig,L)*R2Z
!...dbg:
!       if(L == 9 .and. CAER(I,J,ig) > 0.) then
!MTR       if( i == NRAD ) then
!        print*, j, L, taua(l,j), taus(l,j), g01(l,j)
!MTR       endif
!MTR 115      CONTINUE
!MTR  120     CONTINUE
!MTR       enddo
!
       DO 130 L          = 1,NTOTAL
       TAUAER(L,J)       = MAX(TAUA(NPROB(L),J),EPSILON)
       WOL(L,J)          = TAUS(NPROB(L),J)/TAUAER(L,J)
       TAUAER(L,J)       = TAUA(NPROB(L),J)
       if( WOL(L,J) .ne. 0. ) then
         TTAS = TAUS(NPROB(L),J)
       else
         TTAS = 1.
       endif
!...ej hack:
!       GOL(L,J)          = G01(NPROB(L),J)/TTAS
       GOL(L,J)          = G01(NPROB(L),J)
!       print*, j, l, G01(NPROB(L),J), ttas, gol(l,j)
       if( gol(l,j) > 300. ) stop
 130   CONTINUE
!
 200  CONTINUE
!stop
!
!     iradgas = 0: no gas in radiative xfer
!
      iradgas = 1     

      DO 500 J           = 1,NLAYER
          j1             = max( 1, j-1 )
!
!     Bergstrom water vapor continuum fix:
!
!      <qcorr> is layer average water vapor mixing ratio
!      <pcorr> is layer average pressure
!
!     For layer 0, calculate mixing ratio [g/g] from vapor column
!     [g/cm^2]
!     and average pressure [dyne/cm^2]
!
!MTR          if( j .eq. 1 )then
!MTR             qcorr = rdh2o(1) * g / ptop
!MTR             pcorr = p(1) / 2.
!MTR           else
!MTR             qcorr = q(j1)
!MTR             pcorr = p(j1)
!MTR           endif
!
!MTR           cco = exp(1800./t(j1))*(qcorr*pcorr/2.87 + pcorr/4610.)
!
          DO 400 L       = LLS,LLA
!MTR              TAUH2O(L,J) = RDH2O(J)*PAH2O(L,J)
!
!     Bergstrom water vapor continuum fix (next two statements)
!
!MTR              if( L .GT. NSOLP+30 .AND. L .LE. NSOLP+36 ) then
!MTR                TAUH2O(L,J) = RDH2O(J)*PAH2O(L,J)*cco
!MTR              else
!MTR                TAUH2O(L,J) = TAUH2O(L,J)
!MTR              endif
!MTR              if (L.gt.nsolp+36)    &
!MTR              tauh2o(L,j) = tauh2o(L,j) + cco*rdh2o(j)*contnm(L-nsolp)
 
!MTR              TAUL(L,J)   = TAUH2O(L,J)+TAUGAS(L,J)+   &
!MTR                               PARAY(L,J)+TAUAER(L,J)+TAUCLD(L,J)
              TAUL(L,J)   = TAUGAS(L,J)+   
     &                              TAURAY(L,J)+TAUAER(L,J)+TAUCLD(L,J)
!              write(*,*)'TAUL',J, TAUL(L,J),TAUGAS(L,J)
!              write(*,*)'TAUGAS',  TAUGAS(L,J)
!              write(*,*)'TAURAY',  TAURAY(L,J)
!              write(*,*)'TAUAER',  TAUAER(L,J)
!              write(*,*)'TAUCLD',  TAUCLD(L,J)              
!...dbg:
             if (iradgas.eq.0) tauL(L,j) = tauaer(L,j)
             if( TAUL(L,J) .lt. EPSILON ) TAUL(L,J) = EPSILON
!...dbg:
!          if( L == 80 .and. J < NLAYER ) print*, 'J tauh2o para taua
!          tauc: ', j,   &
!            tauh2o(l,j), paray(L,J), tauaer(L,J), taucld(L,J)
!MTR           if( J < NLAYER )  &
!MTR          write(23,*) J, L, p(j)/1.e3, q(j), tauh2o(l,j), taugas(L,J),
!MTR paray(L,J),  &
!MTR                       tauaer(L,J), taucld(L,J), taul(L,J)
!MTR           if( J < NLAYER )  
!     &    write(23,*) J, L, p(j)/1.e3, taugas(L,J), paray(L,J), 
!     &                    tauaer(L,J), taucld(L,J), taul(L,J)



             utauL(L,j)  = TAUL(L,J)
             WOT         = (TAURAY(L,J)+TAUAER(L,J)*WOL(L,J)+   
     &                           TAUCLD(L,J)*WCLD(L,J))/TAUL(L,J)
!             write(*,*) 'WOT',L,J,WOT
             if (iradgas.eq.0) wot = woL(L,j)

             WOT         = min(1.-EPSILON,WOT)
             uw0(L,j)    = WOT
             DENOM       = (TAURAY(L,J)+TAUCLD(L,J)*WCLD(L,J)+   
     &                          TAUAER(L,J)*WOL(L,J))
             if( DENOM .GT. EPSILON ) then
               GOT = ( WCLD(L,J)*GCLD(L,J)*TAUCLD(L,J) + GOL(L,J)*   
     &                WOL(L,J)*TAUAER(L,J) ) / DENOM
!       print*, j, L, GCLD(l,j), gol(l,j),wol(l,j),got
             else
               GOT = 0.
             endif
             if (iradgas.eq.0) GOT = goL(L,j)
             ug0(L,j)    = GOT
             FO          = GOT**2
             DEN         = 1.-WOT*FO
             TAUL(L,J)   = TAUL(L,J) * DEN
!       print*, j, L, GOT, DEN, TAUL(L,J)
       if( taul(L,j) < 0. ) stop
             W0(L,J)     = (1.-FO)*WOT/DEN
             G0(L,J)     = GOT/(1.+GOT)
             OPD(L,J)    = 0.0
             OPD(L,J)    = OPD(L,J1)+TAUL(L,J)
!MTR             uOPD(L,J)   = 0.0
!MTR             uOPD(L,J)   = uOPD(L,J1)+uTAUL(L,J)
 400        CONTINUE
!
          IF(IR .EQ. 1) THEN
             DO 450 I        =   1,NGAUSS
                DO 425 L     =   LLS,LLA
                   Y3(L,I,J) =   EXP(-TAUL(L,J)/GANGLE(I))
 425            CONTINUE
 450         CONTINUE
          ENDIF
 500  CONTINUE

!       DO J=1,NLAYER
!         write(*,*) ,P_AERAD(J),TAUA(1,J)
!         write(*,*) 'W0',W0 
!       ENDDO
      RETURN
      END

