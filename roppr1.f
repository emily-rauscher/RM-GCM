      SUBROUTINE OPPR1
!
!     **********************************************************
!     *  Purpose             :  Calculate Planck Function and  *
!     *                         and its derivative at ground   *
!     *                         and at all altitudes.          *
!     *  Subroutines Called  :  None                           *
!     *  Input               :  TGRND, NLOW, WEIGHT            *
!     *  Output              :  PTEMP, PTEMPG, SLOPE           *
!     * ********************************************************
!
      include 'rcommons.h'
      real  ITP, ITG, IT1, SBK, SBKoverPI,g11
      DIMENSION  PTEMP2(NTOTAL-NSOLP),PLTEMP1(NTOTAL-NSOLP)
      integer kindex
!     **************************************
!     * CALCULATE PTEMP AND SLOPE          *
!     **************************************
!
!     CALCULATE THE WAVELENGTH DEPENDENT PLANK FUNCTION AT THE GROUND.
!       ITG                 = ANINT(100.*TGRND) - NLOW
!       write(*,*) 'TAUL top of OPPR1',TAUL
!
!     

     
      g11=g0(1,1)
      SBK=5.6704E-8
      SBKoverPI=SBK/PI
!       write(*,*)'ITG,TGRND',ITG,TGRND
!MTR        CALL GATHER(NIRP,PLTEMP1,PLANK(1,ITG),LTEMP)
!NOTE TO ERIN: 
! Here is where the ground temperature is treated seperately, currently
! comented out since it is not defined seperately. 
! ITG --  temperature at the ground, and L is the index over wavelength bins
! (ntotal = 2 bins, Nsolp = number of solar bins i.e.1)
! PTEMPG -- Plank temperature of the ground layer (i.e. sigma T^4 /pi)

!        DO 100 L            =   NSOLP+1,NTOTAL
!           PTEMPG(L)        =  ITG*ITG*ITG*ITG*SBKoverPI! PLTEMP1(L-NSOLP)*WEIGHT(L)
!  100   CONTINUE
!          PTEMPG=ITG*ITG*ITG*ITG*SBKoverPI
!

! THE CODE BELOW IS A MESS. IT DEALS WITH THE 
!      if( iblackbody_above .ne. 0 )then
!
!       CALCULATE THE WAVELENGTH DEPENDENT PLANK FUNCTION AT THE TOP 
!       OF THE MODEL, BUT I SIMPLY SET IT UP SO THE LOOPS BELOW ARE
!       SUFFICIENT.
!MTR         ITP                 = ANINT(100.*t_above) - NLOW
!
!         ITP =  ANINT(T_ABOVE)
!MTR         CALL GATHER(NIRP,PLTEMP1,PLANK(1,ITP),LTEMP)
!         write(*,*),'NSOLP+1,NTOTAL',NSOLP+1,NTOTAL
!         DO 400 L            =   NSOLP+1,NTOTAL
!            PTEMPT(L)        =   PLTEMP1(L-NSOLP)*WEIGHT(L)
!             PTEMPT(L)= ITP*ITP*ITP*ITP*SBKoverPI
!             PTEMP(L,1)=PTEMPT(L)
!             PTEMPG(L)=ITG*ITG*ITG*ITG*SBKoverPI
! 400     CONTINUE
!
!            write(*,*) 'PTEMP(L,1)',PTEMP(L,1)
!            write(*,*) 'PTEMP(L,NLAYER)',PTEMP(L,NLAYER)  
!      endif
!
!          write(*,*) 'TAUL line52',TAUL
!          write(*,*) 'PTEMP',PTEMP
!          write(*,*) 'IT1',IT1
!          PTEMP(L,1)=PTEMPT(L)
!          write(*,*) 'NLAYER:',NLAYER
!             PTEMPT(L)= ITP*ITP*ITP*ITP*SBKoverPI
!             PTEMP(L,1)=PTEMPT(L)
!MTR             PTEMP(L,NLAYER)=ITG*ITG*ITG*ITG*5.67E-8
!          write(*,*) 'TAUL line60',TAUL            
!MTRXXX             write(*,*) 'PTEMPG',PTEMPG

        DO 300 J            =   1,NLAYER ! MTR
!                kindex = j-1
           kindex          = max( 1, j-1 )
!MTR             IT1             = ANINT(100.*TT(J)) - NLOW
!MTR             CALL GATHER(NIRP,PTEMP2,PLANK(1,IT1),LTEMP)
!
!
            IT1 = TT(J)*TT(J)*TT(J)*TT(J)*SBKoverPI
!MTRXXX            write(*,*) 'IT1',IT1
!           KINDEX MAKES AS DEFINED ABOVE MAKES THE TOP LAYER ISOTHERMAL;
!           BELOW THE TOP SLOPE IS REDIFINED USING THE EXTRAPOLATED
!           VALUE. FIND PLANK FUNCTION AT BOTTOM OF EACH LAYER.
!           NOTE: IF YOU FORCE SLOPE=0, THEN YOU HAVE ISOTHERMAL
!           LAYERS WITH TT(J) CORRESPONDING TO AVERAGE TEMPERATURE
!           OF LAYER AND TT(NLAYER) SHOULD BE SET TO TGRND.
!            write(*,*) 'KINDEX',KINDEX
!            write(*,*) 'J',J
            DO 200 L        = NSOLP+1,NTOTAL
!               PTEMP(L,J)   = PTEMP2(L-NSOLP)*WEIGHT(L)
               PTEMP(L,J)=IT1
!               write(*,*)'IT1',IT1
!               write(*,*)'PTEMP(L,J)',PTEMP(L,J)
!               write(*,*)'TAUL(L,J)',TAUL(L,J)
!               write(*,*)'indices',J,KINDEX
              SLOPE(L,J)   = (PTEMP(L,J)-PTEMP(L,KINDEX))/TAUL(L,J)
!               write(*,*)' SLOPE(L,J)',Slope(L,J)
               if( TAUL(L,J) .le. 1.0E-6 ) SLOPE(L,J) = 0.
 200        CONTINUE
 300     CONTINUE
!                   write(*,*) 'TAUL line88',TAUL 
!         DO 450 L            =   NSOLP+1,NTOTAL
!               SLOPE(L,1)   = 0
!               SLOPE(L,NLAYER)=0
!          if( TAUL(L,NLAYER) .gt. 1.0E-6 )
!     &          SLOPE(L,NLAYER)=(PTEMP(L,NLAYER)-PTEMP(L,NLAYER-1))
!     &                         /TAUL(L,NLAYER)
!450     CONTINUE
!      
!         DO J=1,NLAYER
!         write(*,*)'SLOPE',J, SLOPE(2,J)
!         ENDDO
!         DO J=1,NLAYER
!         write(*,*)'PTEMP',J,PTEMP(2,J)
!         ENDDO
!         DO J=1,NLAYER
!        write(*,*)'TAUL',J,TAUL(2,J)
!         ENDDO
          G0(1,1)=g11
!         write(*,*) 'TAU bottom of OPPR1',TAUL
      RETURN
      END

