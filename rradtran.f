      SUBROUTINE RADTRAN
!
!     **************************************************************
!     Purpose:    Driver routine for radiative transfer model.  
!
!     Input:      Temperature, vapor, and aerosol profiles and solar
!                 zenith angle are taken from interface common block.
!
!     Output:     Profiles of radiative fluxes, heating
!                 rates for air and particles; vertically 
!                 integrated optical depths; and albedos (which
!                 are all loaded into interface common block).
!     **************************************************************
!
      include 'rcommons.h'

      integer, parameter :: nwave_alb = NTOTAL
      real wavea(nwave_alb),albedoa(nwave_alb),t(NZ),p(NZ)
      real maxopd(nwave_alb)

      integer jflip

!     Reset flag for computation of solar fluxes

      ISL        = isl_aerad

!     Get atmospheric profiles from interface common block

      do k = 1,nvert
        t(k) = t_aerad(k)
        p(k)=player(k)*10.  ! to convert from Pa to Dyne cm^2
      enddo

!
!     INTERPOLATE TEMPERATURE FROM LAYER CENTER (T) TO LAYER EDGE (TT)
!
!      TT(1) = tabove_aerad
      DO 12 J = 2, NVERT
         TT(J) = T(J-1) * (PRESS(J)/P(J-1)) **    
     &              (log(T(J)/T(J-1))/log(P(J)/P(J-1)))
   12 CONTINUE

!     SINCE WE DON'T HAVE A GROUND (YET,...ERIN)
!     EXTRAPOLATE FOR THE BOUNDARY TEMPERATURES
!     We need a top temperature boundary.  Zero degrees at zero pressure
!     does not work since the slope of such a layer in log p cannot be
!     computed.!     We instead introduce a pressure at 0.5 * P(1)
!     (sigma level). 
!     Since this is not logarithmically spaced, the extrapolation was
!     treated differently at the top.
!     TOP
      TT(1)=((T(1)-TT(2))/log(P(1)/PRESS(2)))*log(PRESS(1)/P(1))+T(1)
!     BOTTOM
      TT(NLAYER)=T(NVERT) * (PRESS(NLAYER)/P(NVERT)) **
     &              (log(T(NVERT)/T(NVERT-1))/log(P(NVERT)/P(NVERT-1)))

!     HERE, INSTEAD OF SPECIFYING THE GROUND AND TOP TEMPERATURES, WE
!     USE THE EXTRAPOLATED VALUES TO DEFINE THESE; ALTERNATIVELY THEY
!     MAY BE DEFINED HERE, WHERE THEY WILL BE PASSED TO THE MODEL
      TGRND=TT(NLAYER)
      TABOVE_AERAD=TT(1)
!      write(*,*)'TGRND',TGRND   
!     WATER VAPOR (G / CM**2)
!    
!     create a T array for the double resolution IR by combinging the
!     layer center and edge temperatures 
      
                 K  =  1
      DO 46 J  = 1, NDBL-1,2
                 L  =  J
         TTsub(L) = TT(K)
                 L  =  L+1
         TTsub(L) = T(K)
                 K  =  K+1
 46    CONTINUE

!     Solar zenith angle 
      u0 = u0_aerad
!
!     SURFACE REFLECTIVITY AND EMISSIVITY
!
!...Hack: use spectrally dependent surface albedo

      DO 20 L =  1,NSOLP
         RSFX(L) = ALBSW  !ALBEDO_SFC
         EMIS(L) =  1.0 - RSFX(L)
 20   CONTINUE

!...Hack: specify EMIS based on RSFX rather than visa versa

      DO 30 L =  NSOLP+1,NTOTAL
         EMIS(L) =  EMISIR
         RSFX(L) = 1.0 - EMIS(L)
        if( wave(nprob(L)).gt.wavea(nwave_alb) ) then
         rsfx(L) = albedoa(nwave_alb)
        endif
         EMIS(L) = 1.0 - RSFX(L)
 30   CONTINUE

!     SET WAVELENGTH LIMITS LLA AND LLS BASED ON VALUES OF ISL AND IR

        LLA                   =  NTOTAL
        LLS                   =  1
        IF(ISL  .EQ. 0) THEN
          LLS   =  NSOLP+1
        ENDIF
!
        IF(IR   .EQ. 0) THEN
          LLA  =  NSOLP
        ENDIF

!     CALCULATE THE OPTICAL PROPERTIES

        IF(AEROSOLCOMP .EQ. 'All') THEN
          CALL OPPRMULTI
        ELSE 
          CALL OPPR
        ENDIF

!     IF INFRARED CALCULATIONS ARE REQUIRED THEN CALCULATE
!     THE PLANK FUNCTION
        
        IF(IR .NE. 0) THEN
           CALL OPPR1
        ENDIF

!     IF NO INFRARED SCATTERING THEN SET INDEX TO NUMBER OF
!     SOLAR INTERVALS

        IF(IRS .EQ. 0) THEN
          LLA  =  NSOLP
        ENDIF

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
        
         !CLOUD FRACTION     
!        NOW, IF WE ARE INCLUDING AEROSOLS, AND WE WOULD LIKE A  CLOUD
!        FRACTION LESS THAN UNITY, THEN RECOMPUTE THESE FLUXES FOR A
!        CLEAR SKY AND COMBINE IN A WEIGHTED AVERAGE ASSUMING MAXIMUM OVERLAP. 

!        HERE WE TAKE THE DOUBLE RESOLUTION FLUXES AND EXTRACT JUST
!        THOSE THAT CORRESPOND TO THE STANDARD (VISIBLE) PRESSURE
!        LEVELS. THESE VALUES SHOULD BE SUPERIOR TO THOSE COMPUTED
!        WITHOUT DOUBLING
                      K     =  1
            DO        J     =  1,NLAYER
             FNET(2,J)      =  DIRECTU(2,k)-DIREC(2,k)  
             DIRECTU(2,J)   =  DIRECTU(2,K)
             DIREC(2,J)     =  DIREC(2,K)
             OPD(2,J)       =  OPD(2,K)
             TAUL(2,J)      =  TAUL(2,K)
             W0(2,J)        =  W0(2,K)
             G0(2,J)        =  G0(2,K)
             ug0(2,J)       =  ug0(2,K)
             uOPD(2,J)      =  uOPD(2,K)
             uTAUL(2,j)     =  uTAUL(2,K)
             uW0(2,J)       =  uW0(2,k)
             TMIU(2,J)      =  TMIU(2,k)
             TMID(2,J)      =  TMID(2,k)
                      K     =  K+2.  
            ENDDO

          IF(FLXLIMDIF) THEN 
            IF(OPD(2,NLAYER) .GT. TAULIMIT) THEN !ASSUMES DOUBLE GRAY,REMOVE THIS LINE OTHERWISE!!
            CALL FLUXLD
            ENDIF
          ENDIF
        ENDIF
          
!     ATTENTION! THE FOLLOWING IS A MODEL-SPECIFIC
!     MODIFICATION:
!     HERE WE PRESCRIBE THE BOTTOM BOUNDARY CONDITION NET FLUX IN THE IR.
!     BE AWARE: IT ALSO AFFECTS THE UPWARD FLUX AT THE BASE IN NEWFLUX.

         FNET(2,NLAYER)=FBASEFLUX  

!     SINCE WE ARE PLACING A CONSTRAINT ON THE NET FLUX AT THE BOTTOM OF
!     THE MODEL (BY DEFINING  NET FLUX AT THE BASE TO EQUAL FBASEFLUX)
!     TO BE SELF CONSISTENT, WE CAN REDIFINE THE UPWARD FLUX FROM THE
!     SURFACE TO BE NETFLUX MINUS DOWNWARD FLUX

!      HERE WE DERIVE THE UPWARD FLUX FROM THE NET FLUX, SELF CONSISTENT
!      WITH BOTTOM BOUNDARY CONDITION

         DIRECTU(2,NLAYER)    =  FBASEFLUX+DIREC(2,NLAYER)
  
!     CALCULATE INFRAFRED AND SOLAR HEATING RATES (DEG/DAY),

           DO 500 J      =  1,NVERT
           HEATS(J)   =  0.0
           HEATI(J)   =  0.0
           TERM1      =  FDEGDAY/(DPG(J+1)*G)
           IF(ISL .NE. 0) THEN
             DO 480 L     =  1,NSOLP
               HEATS(J)   =  HEATS(J)+(FNET(L,J+1)-FNET(L,J))*TERM1
 480         CONTINUE
           ENDIF
           IF(IR .NE. 0) THEN
             DO 490 L     =  NSOLP+1,NTOTAL
               HEATI(J)   =  HEATI(J)+(FNET(L,J+1)-FNET(L,J))*TERM1   
 490         CONTINUE
           ENDIF
            HEAT(J)        =  HEATS(J)+HEATI(J)


!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     !!!!!!!!!!!!!!!         MALSKY CODE          !!!!!!!!!!!!!!!!
!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            CALL TWO_STREAM_WRAPPER

!     Load heating rates [deg_K/s] into interface common block

           heats_aerad(j) =  heats(j)/scday
           heati_aerad(j) =  heati(j)/scday

 500    CONTINUE
 
!     Here we Calculate (4 * pi * mean_intensity) for the IR.
!
      IF (IR .NE. 0) THEN
        DO J = 1, NVERT
          DO L = NSOLP+1, NTOTAL
            TMI(L,J) = TMIU(L,J)+TMID(L,J)
          end do
        end do
      ENDIF


!   Calculate some diagnostic quantities (formerly done in radout.f) and
!   load them into the interface common block.  None of these presently 
!   influence any microphysical processes -- hence, the following code
!   only
!   needs to be executed before the aerosol model writes its output.  
!   Not all of the calculated quantities are presently being
!   loaded into the interface common block.  
!
!   Load optical depths into interface common block

        do i = 1, nwave
          opd_aerad(i) = uopd(i,nlayer)
        enddo
!
!     <tsLu> and <tsLd> are total upwelling and downwelling solar
!     fluxes at top-of-atmosphere
!
        tsLu = 0.
        tsLd = 0.
!
!     <fupbs>, <fdownbs>, and <fnetbs> are total upwelling, downwelling,
!     and net
!     solar fluxes at grid boundaries
!
        do 507 j = 1, nlayer
          fupbs(j) = 0.
          fdownbs(j) = 0.
          fnetbs(j) = 0.
          fdownbs2(j)= 0.
 507    continue
!
!     <fsLu> and <fsLd> are upwelling, downwelling, and net
!     solar fluxes at top-of-atmosphere (spectrally-resolved)
!
!     <alb_toa> is albedo at top-of-atmosphere (spectrally-resolved)
!
        do 509 i = 1, nsoL
          fsLu(i) = 0.
          fsLd(i) = 0.
          alb_toa(i) = 0.
 509    continue
!
!      <alb_tomi> and <alb_toai> are total solar albedos at top-of-model 
!      and top-of-atmosphere
!
        alb_tomi = 0.
        alb_toai = 0.
!
!     CALCULATE SOLAR ABSORBED BY GROUND, SOLNET, AND UPWARD AND
!     DOWNWARD
!     LONGWAVE FLUXES AT SURFACE, XIRUP AND XIRDOWN (WATTS/M**2)
!
        SOLNET  = 0.0
        IF (ISL .NE. 0) THEN
          DO 510 L       =  1,NSOLP
            SOLNET = SOLNET - FNET(L,NLAYER)
            fp = ck1(L,1)*eL2(L,1) - ck2(L,1)*em2(L,1) + cp(L,1)
             fsLu(L) = fsLu(L)+fp
            do 510 j = 1, NLAYER !nlayer
              fp  =  ck1(L,j)* eL1(L,j) + ck2(L,j)*em1(L,j) + cpb(L,j)
              fm  =  ck1(L,j)* eL2(L,j) + ck2(L,j)*em2(L,j) + cmb(L,j)
              fupbs(j) = fupbs(j) + fp
              fdownbs2(j) = fdownbs2(J) + fm
              fnetbs(j) = fnetbs(j) + fnet(L,j)
              if (L.eq.nsolp) then 
              fdownbs(J) = fupbs(j) - fnetbs(j)
              endif

 510      CONTINUE
          do 508 i = 1, nsoL
            fsLd(i) = psol_aerad*u0_aerad !u0*solfx(i)
            alb_toa(i) = fsLu(i)/fsLd(i)
            tsLu = tsLu + fsLu(i)
            tsLd = tsLd + fsLd(i)
 508      continue
          alb_tomi = fupbs(1)/fdownbs(1)
          alb_toai = tsLu/tsLd

!      Load albedos into interface common block
!
          alb_toai_aerad = alb_toai
          alb_tomi_aerad = alb_tomi
          do i = 1, NSOL
            alb_toa_aerad(i) = alb_toa(i)
          enddo
!
!      Load fluxes into interface common block
!
          do j = 1, nlayer
            fsl_up_aerad(j) = fupbs(j)
            fsl_dn_aerad(j) = fdownbs(j)
          enddo

        ENDIF
!
!     <tiru> is total upwelling infrared flux at top-of-atmosphere;
!     <fupbi>, <fdownbi>, and <fnetbi> are total upwelling, downwelling,
!     and net
!     infrared fluxes at grid boundaries
!
        tiru = 0.
        do 606 j = 1, nlayer
          fupbi(j)   =  0.0
          fdownbi(j)   =  0.0
          fnetbi(j)   =  0.0
 606    continue
!
!     <firu> is upwelling infrared flux at top-of-atmosphere
!     (spectrally-resolved)
!
        do 609 i = 1, nir
          firu(i) = 0.
 609    continue

        XIRDOWN = 0.0
        XIRUP   = 0.0

        IF (IR .NE. 0) THEN
          DO 520 L        =  NSOLP+1,NTOTAL
             XIRDOWN = XIRDOWN + DIREC  (L,NLAYER)
             XIRUP   = XIRUP   + DIRECTU(L,NLAYER)
             firu(L-nsol ) = firu( L-nsol ) +    
     &                                directu(L,1)
             do 520 j = 1, nlayer
               fupbi(j) = fupbi(j) + directu(L,j)
               fdownbi(j) = fdownbi(j) + direc  (L,j)
               fnetbi(j) = fnetbi(j) + directu(L,j) - direc(L,j)
 520      CONTINUE

          do 529 i = 1, nir
            tiru = tiru + firu(i)
 529      continue
!
!      Load fluxes into interface common block
!
          do j = 1, nlayer
               jflip=nlayer+1-j
            fir_up_aerad(j) = fupbi(jflip)
            fir_dn_aerad(j) = fdownbi(jflip)
            fir_net_aerad(j)= fnetbi(jflip)
            fsl_up_aerad(j) = fupbs(jflip)
            fsl_dn_aerad(j) = fdownbs(jflip)
            fsl_net_aerad(j)= fnetbs(jflip)
            
!            write(*,*)'fir_up_aerad',j,fir_up_aerad(j)
!            write(*,*)'fir_dn_aerad',j,fir_dn_aerad(j)
          enddo
    
        ENDIF
C     RFLUXES  Array to hold fluxes at top and bottom of atmosphere           
C     1st index - flux 1=SW, 2=LW                                         
C     2nd index - Direction 1=DN, 2=UP                                    
C     3rd index - Where 1=TOP, 2=SURFACE   
      RFLUXES_aerad(1,1,1)=fsl_dn_aerad(NLAYER)   ! SW down top                            
      RFLUXES_aerad(1,1,2)=fsl_dn_aerad(1)/(1.0-ALBSW)   ! SW down bottom                
      RFLUXES_aerad(1,2,1)=fsl_up_aerad(NLAYER)  ! SW up top             
      RFLUXES_aerad(1,2,2)=RFLUXES_aerad(1,1,2)*ALBSW   ! SW up bottom                                                                       
      RFLUXES_aerad(2,1,1)=fir_dn_aerad(NLAYER)   ! LW down top                           
      RFLUXES_aerad(2,1,2)=fir_dn_aerad(1)       ! LW down bottom                       
      RFLUXES_aerad(2,2,1)=fir_up_aerad(NLAYER)       ! LW up top                          
      RFLUXES_aerad(2,2,2)=fir_up_aerad(1)   ! LW up bottom  
!      write(*,*)'RFLUXES_aerad',rfluxes_aerad
      return
      END

