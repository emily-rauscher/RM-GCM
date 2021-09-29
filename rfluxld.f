      SUBROUTINE FLUXLD
!
!     *******************************************************************
!     *  Purpose             :  Evaluate the IR flux to see if it meets * 
!     *                         the criteria for Flux Limited Diffusion *
!     *                         (see Rauscher & Menou,2012), and        *
!     *                         replaces it if so to compute heating    * 
!     *                         The equations were copied from Cnikos.f *
!     *                         following coding by Emily R. ~MTR       * 
!     *  Subroutines Called  :  None                                    *
!     *  Input               :  TGRND, NLOW, WEIGHT                     *
!     *  Output              :  PTEMP, PTEMPG, SLOPE                    *
!     * ****************************************************************
!
      include 'rcommons.h'
      real  EXPCORR,GRAD,FNETDIFF,TFAC,TAULIMIT
      DIMENSION PTEMP2(NTOTAL-NSOLP),PLTEMP1(NTOTAL-NSOLP)
      DIMENSION T(NLAYER)
      integer kindex
!     **************************************     
       T=T_aerad
       DO 100 L     =  NSOLP+1,NTOTAL  ! LOOP OVER IR BINS IF ANY
         DO 200 ILEV  =  1,NLAYER  ! LOOP OVER LAYERS
!          IF THIS CONDITIONAL, BASED ON RESOLUTION DEFINED IN INVARPARAM, IS MET...
           IF (OPD(L,ILEV).GT.TAULIMIT) THEN

!      CHOOSE THE OPACITY POWER LAW CASE; IN PRACTICE, THIS HAS FOUND TO
!      REDUCE THE COMPUTING TIME SURPRISINGLY
!      FOR THE OPACITY POWER LAW = 0 (CONSTANT IN HEIGHT) CASE:

       IF (OPACIR_POWERLAW.EQ.0) THEN
           IF (ILEV.EQ.1) THEN !IF at the top of the atmosphere, use boundaries
               GRAD=(TT(1)-T(1))/((P_AERAD(1)-PLAYER(1))/100.)
               EXPCORR=TAUCONST(L)*1.E6*0.5*P_AERAD(ILEV)
           ELSEIF (ILEV.EQ.NLAYER) THEN !IF AT THE BOTTOM, USE BOTTOM BOUNDARY
               GRAD=(T(ILEV-1)-T(ILEV))/(PLAYER(ILEV-1)-PLAYER(ILEV))*100.
               EXPCORR=TAUCONST(L)*1.E6*0.5*(2*P_AERAD(ILEV)-PLAYER(ILEV)-PLAYER(ILEV-1))
           ELSE
               GRAD=(T(ILEV-1)-T(ILEV))/(PLAYER(ILEV-1)-PLAYER(ILEV))*100.
               EXPCORR=TAUCONST(L)*1.E6*0.5*(PLAYER(ILEV)-PLAYER(ILEV-1))
           ENDIF
           FNETDIFF=3.0242E-10*GRAD*(TT(ILEV)*TT(ILEV)*TT(ILEV))/TAUCONST(L)

!          FOR AN OPACITY POWERLAW OF UNITY
!          111111111111111111111111111111111111111111111111111111111111111111111
           ELSEIF (OPACIR_POWERLAW.EQ.1) THEN
               IF (ILEV.EQ.1) THEN !IF at the top of the atmosphere, use boundaries
                   GRAD=(TT(1)-T(1))/((P_AERAD(1)-PLAYER(1))/100.)
                   EXPCORR=TAUCONST(L)*1.E6*0.5*P_AERAD(ILEV)*(1.E2*P_AERAD(ILEV)/OPACIR_REFPRES)
               ELSEIF (ILEV.EQ.NLAYER) THEN !IF AT THE BOTTOM, USE BOTTOM BOUNDARY
                   GRAD=(T(ILEV-1)-T(ILEV))/(PLAYER(ILEV-1)-PLAYER(ILEV))*100.
                   EXPCORR=TAUCONST(L)*1.E6*0.5
     &             *(2*P_AERAD(ILEV)-PLAYER(ILEV)-PLAYER(ILEV-1))
     &              *(1.E2*P_AERAD(ILEV)/OPACIR_REFPRES)
               ELSE
                   GRAD=(T(ILEV-1)-T(ILEV))/(PLAYER(ILEV-1)-PLAYER(ILEV))*100.
                   EXPCORR=TAUCONST(L)*1.E6*0.5*(PLAYER(ILEV)-PLAYER(ILEV-1))
     &              *(1.E2*P_AERAD(ILEV)/OPACIR_REFPRES)
           ENDIF
           FNETDIFF=3.0242E-10*GRAD*(TT(ILEV)*TT(ILEV)*TT(ILEV))/TAUCONST(L)*(OPACIR_REFPRES/P_AERAD(ILEV)/1.E2)
!
!            NOW FOR AN OPACITY POWERLAW EQUAL TO TWO, JUST ANOTHER FACTOR 
!         22222222222222222222222222222222222222222222222222222222222222222222222
          ELSEIF (OPACIR_POWERLAW.EQ.2) THEN
            IF (ILEV.EQ.1) THEN !IF at the top of the atmosphere, use boundaries
               GRAD=(TT(1)-T(1))/((P_AERAD(1)-PLAYER(1))/100.)
               EXPCORR=TAUCONST(L)*1.E6*0.5*P_AERAD(ILEV)
     &              *(1.E2*P_AERAD(ILEV)/OPACIR_REFPRES)
     &              *(1.E2*P_AERAD(ILEV)/OPACIR_REFPRES)

            ELSEIF (ILEV.EQ.NLAYER) THEN !IF AT THE BOTTOM, USE BOTTOM BOUNDARY
               GRAD=(T(ILEV-1)-T(ILEV))/(PLAYER(ILEV-1)-PLAYER(ILEV))*100.
               !NB. *100. CONVERTS DENOMINATOR FROM PASCALS TO MBAR FOLLOWING ER
               EXPCORR=TAUCONST(L)*1.E6*0.5
     &              *(2*P_AERAD(ILEV)-PLAYER(ILEV)-PLAYER(ILEV-1))
     &              *(1.E2*P_AERAD(ILEV)/OPACIR_REFPRES)
     &              *(1.E2*P_AERAD(ILEV)/OPACIR_REFPRES)
            ELSE
               GRAD=(T(ILEV-1)-T(ILEV))/(PLAYER(ILEV-1)-PLAYER(ILEV))*100.
              EXPCORR=TAUCONST(L)*1.E6*0.5*(PLAYER(ILEV)-PLAYER(ILEV-1))
     &              *(1.E2*P_AERAD(ILEV)/OPACIR_REFPRES)
     &              *(1.E2*P_AERAD(ILEV)/OPACIR_REFPRES)
            ENDIF
       FNETDIFF=3.0242E-10*GRAD*(TT(ILEV)*TT(ILEV)*TT(ILEV))/TAUCONST(L)
     &           *(OPACIR_REFPRES/P_AERAD(ILEV)/1.E2)
     &           *(OPACIR_REFPRES/P_AERAD(ILEV)/1.E2)


!          NOW FOR AN OPACITY POWERLAW EQUAL TO THREE, KEEP ADDING FACTORS...
!         333333333333333333333333333333333333333333333333333333333333333333333333
          ELSEIF (OPACIR_POWERLAW.EQ.3) THEN
            IF (ILEV.EQ.1) THEN !IF at the top of the atmosphere, use boundaries
               GRAD=(TT(1)-T(1))/((P_AERAD(1)-PLAYER(1))/100.)
               EXPCORR=TAUCONST(L)*1.E6*0.5*P_AERAD(ILEV)
     &              *(1.E2*P_AERAD(ILEV)/OPACIR_REFPRES)
     &              *(1.E2*P_AERAD(ILEV)/OPACIR_REFPRES)
     &              *(1.E2*P_AERAD(ILEV)/OPACIR_REFPRES)
            ELSEIF (ILEV.EQ.NLAYER) THEN !IF AT THE BOTTOM, USE BOTTOM BOUNDARY
               GRAD=(T(ILEV-1)-T(ILEV))/(PLAYER(ILEV-1)-PLAYER(ILEV))*100.
               !NB. *100. CONVERTS DENOMINATOR FROM PASCALS TO MBAR FOLLOWING ER
               EXPCORR=TAUCONST(L)*1.E6*0.5
     &              *(2*P_AERAD(ILEV)-PLAYER(ILEV)-PLAYER(ILEV-1))
     &              *(1.E2*P_AERAD(ILEV)/OPACIR_REFPRES)
     &              *(1.E2*P_AERAD(ILEV)/OPACIR_REFPRES)
     &              *(1.E2*P_AERAD(ILEV)/OPACIR_REFPRES)
            ELSE
               GRAD=(T(ILEV-1)-T(ILEV))/(PLAYER(ILEV-1)-PLAYER(ILEV))*100.
              EXPCORR=TAUCONST(L)*1.E6*0.5*(PLAYER(ILEV)-PLAYER(ILEV-1))
     &              *(1.E2*P_AERAD(ILEV)/OPACIR_REFPRES)
     &              *(1.E2*P_AERAD(ILEV)/OPACIR_REFPRES)
     &              *(1.E2*P_AERAD(ILEV)/OPACIR_REFPRES)
            ENDIF
       FNETDIFF=3.0242E-10*GRAD*(TT(ILEV)*TT(ILEV)*TT(ILEV))/TAUCONST(L)
     &           *(OPACIR_REFPRES/P_AERAD(ILEV)/1.E2)
     &           *(OPACIR_REFPRES/P_AERAD(ILEV)/1.E2)
     &           *(OPACIR_REFPRES/P_AERAD(ILEV)/1.E2)

!       OR ELSE, AN ARBITRARY POWERLAW VALUE WHERE THE EXPONENTIATION IS ACTUALLY COMPUTED
          ELSE
            IF (ILEV.EQ.1) THEN !IF at the top of the atmosphere, use boundaries
               GRAD=(TT(1)-T(1))/((P_AERAD(1)-PLAYER(1))/100.)
               EXPCORR=TAUCONST(L)*1.E6*0.5*P_AERAD(ILEV)
     &             *(1.E2*P_AERAD(ILEV)/OPACIR_REFPRES)**OPACIR_POWERLAW
            ELSEIF (ILEV.EQ.NLAYER) THEN !IF AT THE BOTTOM, USE BOTTOM BOUNDARY
               GRAD=(T(ILEV-1)-T(ILEV))/(PLAYER(ILEV-1)-PLAYER(ILEV))*100.
               !NB. *100. CONVERTS DENOMINATOR FROM PASCALS TO MBAR FOLLOWING ER
               EXPCORR=TAUCONST(L)*1.E6*0.5
     &              *(2*P_AERAD(ILEV)-PLAYER(ILEV)-PLAYER(ILEV-1))
     &             *(1.E2*P_AERAD(ILEV)/OPACIR_REFPRES)**OPACIR_POWERLAW
            ELSE
               GRAD=(T(ILEV-1)-T(ILEV))/(PLAYER(ILEV-1)-PLAYER(ILEV))*100.
              EXPCORR=TAUCONST(L)*1.E6*0.5*(PLAYER(ILEV)-PLAYER(ILEV-1))
     &             *(1.E2*P_AERAD(ILEV)/OPACIR_REFPRES)**OPACIR_POWERLAW
            ENDIF
       FNETDIFF=3.0242E-10*GRAD*(TT(ILEV)*TT(ILEV)*TT(ILEV))/TAUCONST(L)
     &           *(OPACIR_REFPRES/P_AERAD(ILEV)/1.E2)**OPACIR_POWERLAW
          ENDIF

!       NOW, REGARDLESS OF WHAT POWERLAW WAS USED, COMPLETE THE CALCULATION
            TFAC=1.-EXP(-1.66*EXPCORR/43.7)
            FNET(L,ILEV)=(1.0-TFAC)*FNET(L,ILEV) + TFAC*FNETDIFF
       ENDIF

        
200   CONTINUE      
100   CONTINUE
      RETURN
      END

