! *********************************************************************
!
!     COMPUTE PLANCK FUNCTION TABLE. WAVE IS IN UNITS OF MICRONS.
!
! **********************************************************************
!
!
      include 'rcommons.h'
!     Set <iblackbody_above> = 1 to include a source of radiation
!     at the top of the radiative transfer model domain
!    
      iblackbody_above = ir_above_aerad
      t_above = tabove_aerad
!
!     Set <ibeyond_spectrum> = 1 to include blackbody radiation at
!     wavelengths longer than WAVE(NWAVE+1) in PLANK(NWAVE+1-NSOL) and
!     at wavelengths shorter than WAVE(NSOL+1) in PLANK(1)
!
      ibeyond_spectrum = 0
!
!  ;COMMENTING ALL THIS OUT TO SIMPLY COMPUTE SIGMA T^4 FOR BLACKBODY
!  EMMISION
!   DIRECTLY LATER, RATHER THAN PRODUCING A TABLE OF VALUES USING THE
!   NLOW TO
!   NHIGH VALUES IN THE GLOBRAD.H

      DO 380 J  =   1,NCOUNT
         JJ     =   NLOW+J
         T1     =   0.01 * FLOAT(JJ)

         if( ibeyond_spectrum .eq. 1 )then

           plank(nwave+1-nsol,j) = t1**4
           DO I =   NSOL+2,NWAVE
              K =   I-NSOL
              V =   1.438E4  /  WAVE(I)
            CALL PLNK(V,T1,PLANK(K,J))
           ENDDO

         else

           DO I =   NSOL+1,NWAVE+1
              K =   I-NSOL
              V =   1.438E4  /  WAVE(I)
              CALL PLNK(V,T1,PLANK(K,J))
           ENDDO

         endif

 380  CONTINUE
      DO 410 J   =   1,NCOUNT

         if( ibeyond_spectrum .eq. 1 )then

           plank(1,j) = plank(2,j)*sbk/pi
           DO L  =   NSOL+2,NWAVE
              K  =   L-NSOL
              PLANK(K,J) = (PLANK(K+1,J)-PLANK(K,J))*SBK/PI
           ENDDO

         else

           DO L  =   NSOL+1,NWAVE
              K  =   L-NSOL
              PLANK(K,J) = (PLANK(K+1,J)-PLANK(K,J))*SBK/PI
           ENDDO

         endif

 410  CONTINUE

      RETURN
      END


