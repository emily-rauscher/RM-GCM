      SUBROUTINE CLOUDFRACTER

!     ***************************************************************
!     Purpose:    Recomputes fluxes for cases with cloud fractions <1.  
!
!     Input:      Cloud fraction and scattering paramters taken from
!                 taken from interface common block.
!
!     Output:     Fluxes and intensities for the weighted fraction of 
!                 clear and cloudy skies assuming a maximum overlap.
!     **************************************************************

      include 'rcommons.h'
      REAL FNETO(NTOTAL,NLAYER),DIRECO(NTOTAL,NLAYER)
     &     ,DIRECTUO(NTOTAL,NLAYER),AEROPROFO(NLAYER)
     &    ,UINTENTO(NTOTAL,NGAUSS,NLAYER),DINTENTO(NTOTAL,NGAUSS,NLAYER)
     &     ,TMIDO(NTOTAL,NLAYER),TMIUO(NTOTAL,NLAYER)
     &     ,W0O(NTOTAL,NLAYER),G0O(NTOTAL,NLAYER)
     &     ,OPDO(NTOTAL,NLAYER),TAULO(NTOTAL,NLAYER)
     &     ,TMIO(NTOTAL,NLAYER),CK1O(NTOTAL,NLAYER)
     &     ,EL1O(NTOTAL,NLAYER),CK2O(NTOTAL,NLAYER)
     &     ,EM1O(NTOTAL,NLAYER),CPBO(NTOTAL,NLAYER)
     &     ,EL2O(NTOTAL,NLAYER),EM2O(NTOTAL,NLAYER)
     &     ,CMBO(NTOTAL,NLAYER)
!     POTENTIALLY OVERCAST RESULTS WERE COMPUTED PREVIOUSLY, ROUTINELY.     
!     FIRST, SAVE THE PREVIOUS OVERCAST RESULTS SO AS NOT TO OVERWRITE
!     THEM WHEN COMPUTING NEW FLUXES AND INTENSITIES. 
      FNETO     =  FNET
      DIRECO    =  DIREC
      DIRECTUO  =  DIRECTU
      TMIUO     =  TMIU
      TMIDO     =  TMID
      TMIO      =  TMI
      W0O       =  W0
      G0O       =  G0
      OPDO      =  OPD
      TAULO     =  TAUL
      UINTENTO  =  UINTENT
      DINTENTO  =  DINTENT
      CK1O      =  CK1
      EL1O      =  EL1
      CK2O      =  CK2
      EM1O      =  EM1
      CPBO      =  CPB
      EL2O      =  EL2
      EM2O      =  EM2
      CMBO      =  CMBO
      AEROPROFO =  AEROPROF

!     REMOVE THOSE CLOUDS--WE WANT CLEAR SKY FLUXES
         DO J  = 1, NLAYER
           AEROPROF(J) = 0.0
         ENDDO
!     NOW CALL ALL THOSE PROGRAMS AGAIN
      CALL OPPR
        IF(IR .NE. 0) THEN
      CALL OPPR1
        ENDIF
!     IF NO INFRARED SCATTERING THEN SET INDEX TO NUMBER OF
!     SOLAR INTERVALS
!
        IF(IRS .EQ. 0) THEN
          LLA  =  NSOLP
        ENDIF
!
!     IF EITHER SOLAR OR INFRARED SCATTERING CALCULATIONS ARE REQUIRED
!     CALL THE TWO STREAM CODE AND FIND THE SOLUTION
        IF(ISL .NE. 0 .OR. IRS .NE. 0 ) THEN
          CALL TWOSTR
          CALL ADD
        ENDIF
!     IF INFRARED CALCULATIONS ARE REQUIRED THEN CALL NEWFLUX1 FOR
!     A MORE ACCURATE SOLUTION
        IF(IR .NE. 0) THEN
         CALL NEWFLUX1
        ENDIF
!     NOW WITH EVERYTHING RECOMPUTED, WE HAVE TO COMPUTE THE EFFECTIVE
!     NET FLUXES BASED ON THE CLOUD COVER ASSUMING A WEIGHTED AVERAGE
!     CONSISTENT WITH A MAXIMUM OVERLAP SCENARIO. FOR NON-AVERAGED THAT
!     ARE SUBSEQUENTLY NOT CONSEQUENTIAL, ONLY DIAGNOSTIC, JUST SET TO 
!     THE CLOUDY VALUE

        DO 556 J  = 1, NLAYER
         DO 555 L = 1, NTOTAL
           FNET(L,J) = (1.-CLDFRCT)*FNET(L,J)+CLDFRCT*FNETO(L,J)
           DIREC(L,J) = (1.-CLDFRCT)*DIREC(L,J)+CLDFRCT*DIRECO(L,J)
           DIRECTU(L,J)=(1.-CLDFRCT)*DIRECTU(L,J)+CLDFRCT*DIRECTUO(L,J)
           TMIU(L,J)  = (1.-CLDFRCT)*TMIU(L,J)+CLDFRCT*TMIUO(L,J)
           TMID(L,J)  = (1.-CLDFRCT)*TMID(L,J)+CLDFRCT*TMIDO(L,J)
           OPD(L,J)  =  OPDO(L,J)
           W0(L,J)   =  W0O(L,J)
           G0(L,J)   =  G0O(L,J)
           TAUL(L,J) =  TAULO(L,J)
!          DINTENT(L,J)=(1.-CLDFRCT)*DINTENT(L,J)+CLDFRCT*DINTENTO(L,J)
!          UINTENT(L,J)=(1.-CLDFRCT)*UINTENT(L,J)+CLDFRCT*UINTENTO(L,J)
 555     CONTINUE
          AEROPROF(J) = AEROPROFO(J)
 556    CONTINUE
    
       END
