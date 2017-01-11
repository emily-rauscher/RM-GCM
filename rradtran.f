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
!      integer, parameter :: nwave_alb = 142
      real wavea(nwave_alb),albedoa(nwave_alb),t(NZ),p(NZ)
      real maxopd(nwave_alb)
!      real RFLUXES_aerad(2,2,2) 
      integer jflip
!  Reset flag for computation of solar fluxes
!
      ISL        = isl_aerad
!
!     Get atmospheric profiles from interface common block
!
      do k = 1,nvert
        t(k) = t_aerad(k)
!       qv(k) = qv_aerad(k)
        p(k)=player(k)*10.  ! to convert from Pa to Dyne cm^2
!        p(k)=player(k)
      enddo
!        t(NLAYER)=tgrnd !grnd temp from upwelling base flux +downwelling
!     write(*,*),'p_aerad',p_aerad
!      write(*,*)'P',P
!     write(*,*)'PRESS',PRESS
      
!...dbg:
!      close(21)
!      open(unit=21,file='testprof.dat',form='formatted',status='unknown')
!      write(21,*) p, t
!      close(21)
!      stop
!
!     INTERPOLATE TEMPERATURE FROM LAYER CENTER (T) TO LAYER EDGE (TT)
!
!      TT(1) = tabove_aerad
      DO 12 J = 2, NVERT
         TT(J) = T(J-1) * (PRESS(J)/P(J-1)) **    
     &              (log(T(J)/T(J-1))/log(P(J)/P(J-1)))
   12 CONTINUE
!      TT(NLAYER)=TGRND
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
!      write(*,*) 'T',T
!      write(*,*) 'TT',TT


!M      DO 10 J = 2, NLAYER
!M         RDH2O(J)   = Q(J-1) * DPG(J-1)
!M   10 CONTINUE
!
!     AEROSOL CONCENTRATIONS (# / CM**2)
!
!MTR      do ig = 1, NGROUP
!MTR        DO 15 J = 2, NVERT
!MTR           DO 15 I = 1, NRAD
!              CAER(I,J,ig)  = pc_aerad(J-1,I,ig) 
!MTR   15   CONTINUE
!MTR      enddo
!
!     Solar zenith angle 
!
      u0 = u0_aerad
!
!     SURFACE REFLECTIVITY AND EMISSIVITY
!
!...Hack: use spectrally dependent surface albedo
!
!M      open(21,file='arctas_albedo.dat',status='unknown',form='formatted')
!      open(21,file='crystal-FACE_open-ocean_albedo.dat',status='unknown',
!      &
!           form='formatted')
!M      do i = 1, nwave_alb
!M        read(21,*) wavea(i), albedoa(i)
!        print*, i, wavea(i), albedoa(i)
!M      enddo
!M      close(21)
!      stop
!M      wavea = wavea / 1.e3
      DO 20 L =  1,NSOLP
         RSFX(L) = ALBSW  !ALBEDO_SFC
!         EMIS(L) =  1.0
!         write(*,*) 'EMIS', EMIS(L)
!      write(*,*) ' calling interpol'
!        call interpol( albedoa, wavea, wave(nprob(L)),
!     &                  nwave_alb, 1, rsfx(L) )
!         if( wave(nprob(L)) > wavea(nwave_alb) ) 
!     &           rsfx(L) = albedoa(nwave_alb)
         EMIS(L) =  1.0 - RSFX(L)
!         print*, L, wave(nprob(L)), rsfx(L)
!      write(*,*) 'L, wave(nprob(L)), rsfx(L),albedoa(nwave_alb)'
!      write(*,*) L, wave(nprob(L)), rsfx(L),albedoa(nwave_alb)
 20   CONTINUE
!...Hack: specify EMIS based on RSFX rather than visa versa
      DO 30 L =  NSOLP+1,NTOTAL
         EMIS(L) =  EMISIR
         RSFX(L) = 1.0 - EMIS(L)

!         call interpol( albedoa, wavea, wave(nprob(L)), 
!     &    nwave_alb, 1, rsfx(L) )
        if( wave(nprob(L)).gt.wavea(nwave_alb) ) then
         rsfx(L) = albedoa(nwave_alb)
        endif
         EMIS(L) = 1.0 - RSFX(L)
!         print*, L, wave(nprob(L)), rsfx(L), emis(l)
 30   CONTINUE
!      write(*,*) 'just interpoled'
!
!     SET WAVELENGTH LIMITS LLA AND LLS BASED ON VALUES OF ISL AND IR
!
        LLA                   =  NTOTAL
        LLS                   =  1
        IF(ISL  .EQ. 0) THEN
          LLS   =  NSOLP+1
        ENDIF
!
        IF(IR   .EQ. 0) THEN
          LLA  =  NSOLP
        ENDIF
!
!     CALCULATE THE OPTICAL PROPERTIES
!
!        write(*,*) ' calling OPPR'
        CALL OPPR
!        write(*,*) 'called OPPR'
!     IF INFRARED CALCULATIONS ARE REQUIRED THEN CALCULATE
!     THE PLANK FUNCTION
        
!       write(*,*)'TUAL pre OPPR1',TAUL
!        write(*,*) ' calling OPPR1'
        IF(IR .NE. 0) THEN
           CALL OPPR1
        ENDIF
!       write(*,*)'TAUL AFTER OPPR1',TAUL
!        write(*,*) 'PTEMP',PTEMP
!        write(*,*) 'SLOPE', SLOPE
!     IF NO INFRARED SCATTERING THEN SET INDEX TO NUMBER OF
!     SOLAR INTERVALS
!
        IF(IRS .EQ. 0) THEN
          LLA  =  NSOLP
        ENDIF
!
!     IF EITHER SOLAR OR INFRARED SCATTERING CALCULATIONS ARE REQUIRED
!     CALL THE TWO STREAM CODE AND FIND THE SOLUTION
!        write(*,*) 'RSFX in RADTRAN',RSFX
        IF(ISL .NE. 0 .OR. IRS .NE. 0 ) THEN
          CALL TWOSTR
          CALL ADD
        ENDIF
!     IF INFRARED CALCULATIONS ARE REQUIRED THEN CALL NEWFLUX1 FOR
!     A MORE ACCURATE SOLUTION
!        write(*,*) 'FNET(radtran)',FNET       
        IF(IR .NE. 0) THEN
!        write(*,*) 'HEATI',HEATI
         CALL NEWFLUX1
!        write(*,*) 'HEATI',HEATI
          
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
!          write(*,*), NSOLP
!MTR          write(*,*) 'NLAYER',NLAYER
!        DO IJ=1,NLAYER 
!           write(*,*) 'FNET',IJ,P_aerad(IJ),FNET(2,IJ)
!        ENDDO 
!           write(*,*) 'FBASEFLUX',FBASEFLUX
!        DO IJ=1,NLAYER
!           write(*,*) 'DPG*G',IJ,DPG(IJ)*G
!        ENDDO
!         write(*,*)'PSOL_aerad*u0_aerad', PSOL_aerad*u0_aerad
!         write(*,*)'Index,Mass factor, Flux Divergence'
           DO 500 J      =  1,NVERT !!MTRNVERT
           HEATS(J)   =  0.0
           HEATI(J)   =  0.0
!           write(*,*)  'DPG' 
           TERM1      =  FDEGDAY/(DPG(J+1)*G)
!           write(*,*) J,TERM1         
!           write(*,*)'FNET',FNET
           IF(ISL .NE. 0) THEN
             DO 480 L     =  1,NSOLP
               HEATS(J)   =  HEATS(J)+(FNET(L,J+1)-FNET(L,J))*TERM1
!               write(*,*)'FNET',(FNET(L,J+1),FNET(L,J))
 480         CONTINUE
           ENDIF
!           write(*,*) 'IR= ',IR
           IF(IR .NE. 0) THEN
!             write(*,*) 'IR.NE.0'
             DO 490 L     =  NSOLP+1,NTOTAL
                HEATI(J)  =  HEATI(J)+( DIRECTU(L,J+1)-DIREC(L,J+1)    
     &                        -(DIRECTU(L,J)-DIREC(L,J)) )*TERM1
   
!MTR      HERE'S WHERE YOU WOULD CALL FLUX LIMITED DIFFUSION CODE
!MTR                write(*,*) 'HEATI(J)',HEATI(J)
!MTR                 write(*,*) ' DIRECTU(L,J+1)', DIRECTU(L,J+1)
!MTR                 write(*,*) ' DIRECTU(L,J)', DIRECTU(L,J)
!MTR                 write(*,*) ' DIREC(L,J+1)', DIREC(L,J+1)
!MTR                 write(*,*) ' DIREC(L,J)', DIREC(L,J)

 490         CONTINUE
!MTR                    write(*,*) 'PRESS', PRESS
!MTR                    write(*,*)'Player',PLAYER
!MTR                    write(*,*) 'OPD',OPD
!MTR                    write(*,*) 'T',T,TT
!MTR                    write(*,*) 'OOM_IN',OOM_IN
!MTR                    write(*,*) 'NL',NL
!MTR                    write(*,*) 'OPACIR_POWERLAW',OPACIR_POWERLAW
!MTR                    write(*,*) 'OPACIR_REFPRES',OPACIR_REFPRES
!MTR                    write(*,*) 'stop'
!MTR             stop
           ENDIF
           HEAT(J)        =  HEATS(J)+HEATI(J)
!
!     Load heating rates [deg_K/s] into interface common block
!          
!            write(*,*) 'heats', heats
!            write(*,*) 'scday', scday
           heats_aerad(j) =  heats(j)/scday
           heati_aerad(j) =  heati(j)/scday
!            write(*,*) 'heats_aerad_in radtran',heats_aerad(j)
!           write(*,*) 'j,bar,solar,IR',j,heats_aerad(j),heati_aerad(j)
 500    CONTINUE
!           DO 502 IJ=1,NVERT !NLAYER
!           write(*,*) 'LWUP,DN',IJ,P_aerad(IJ)
!           write(*,*) '     ',DIRECTU(2,IJ),DIREC(2,IJ)
!            write(*,*) 'FNETS(IJ+1,IJ)',IJ,FNET(1,IJ+1),FNET(1,IJ)
! 502   CONTINUE        
 
!     Here we Calculate (4 * pi * mean_intensity) for the IR.
!
      IF (IR .NE. 0) THEN
        DO J = 1, NVERT
          DO L = NSOLP+1, NTOTAL
            TMI(L,J) = TMIU(L,J)+TMID(L,J)
          end do
        end do
      ENDIF

!     Here we compute the heating rates for droplets
!     (C11 converts W/m^2 to erg/cm^2)
!
!MTR      C11 = 1000.
!MTR      do ig = 1, NGROUP
!MTR       DO I = 1, NRAD
!MTR        DO J = 1, NVERT
! 
!MTR         QRAD(j) = 0.
!
!MTR         IF (IR .NE. 0) THEN
!MTR           DO L = NSOLP+1, NTOTAL
!MTR             X = TMI(L,J)-4.0*PI*PTEMP(L,J)
!MTR             if( abs(X/TMI(L,J)) .lt. EPSILON ) x = 0.
!MTR             QRAD(j) = QRAD(j) + X*C11*XSECTA(I,ig) *    
!MTR      &      (RDQEXT(I,ig,NPROB(L))-QSCAT(I,ig,NPROB(L)))
!MTR           end do
!MTR         ENDIF
!        
!MTR         IF (ISL .NE. 0) THEN
!MTR           DO L = 1, NSOLP
!MTR             QRAD(j) = QRAD(j) + TMI(L,J)*C11*XSECTA(I,ig) *    
!MTR      &      (RDQEXT(I,ig,NPROB(L))-QSCAT(I,ig,NPROB(L)))
!MTR           end do
!MTR         ENDIF
!
!     Load layer averages of droplet heating rates into interface common
!     block
!
!        if (j.eq.nvert) then
!          qrad_aerad(i,j,ig) = qrad(j)
!        else if (j.gt.1) then
!          print*, i, j, ig
!          qrad_aerad(i,j-1,ig) = 0.5 * ( qrad(j) + qrad(j-1) )
!        endif
!
!MTR         end do     ! j=1,nvert
!MTR        end do      ! i=1,nrad
!MTR      end do       ! ig=1,ngroup
!
!
!   Calculate some diagnostic quantities (formerly done in radout.f) and
!   load them into the interface common block.  None of these presently 
!   influence any microphysical processes -- hence, the following code
!   only
!   needs to be executed before the aerosol model writes its output.  
!   Not all of the calculated quantities are presently being
!   loaded into the interface common block.  
!
!   Load optical depths into interface common block
!
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
!        write(*,*) 'SOLAR: DOWN, UP, NET'
        IF (ISL .NE. 0) THEN
          DO 510 L       =  1,NSOLP
            SOLNET = SOLNET - FNET(L,NLAYER)
            fp = ck1(L,1)*eL2(L,1) - ck2(L,1)*em2(L,1) + cp(L,1)
!            fsLu( nprob(L) ) = fsLu( nprob(L) ) + fp
             fsLu(L) = fsLu(L)+fp
!             write(*,*) 'solnet',SOLNET
!            write(*,*),'fsLu(L)',fsLu(L)
!            write(*,*),'ck1(L,1)',ck1(L,1)
!            write(*,*),'el2(L,1),',eL2(L,1)
!            write(*,*),'ck2(L,1)',ck2(L,1)
!            write(*,*),'em2(L,1),',em2(L,1)
!            write(*,*),'cp', cp(L,1)
!            write(*,*),'product fp=',ck1(L,1)*eL2(L,1) -
!     &  ck2(L,1)*em2(L,1) + cp(L,1)
            do 510 j = 1, NLAYER !nlayer
              fp  =  ck1(L,j)* eL1(L,j) + ck2(L,j)*em1(L,j) + cpb(L,j)
              fm  =  ck1(L,j)* eL2(L,j) + ck2(L,j)*em2(L,j) + cmb(L,j)
              fupbs(j) = fupbs(j) + fp
              fdownbs2(j) = fdownbs2(J) + fm
              fnetbs(j) = fnetbs(j) + fnet(L,j)
              if (L.eq.nsolp) then 
              fdownbs(J) = fupbs(j) - fnetbs(j)
              endif
!              write(*,*),j,fupbs(j),fdownbs2(j)
!              write(*,*),j,fnetbs(j),DIRECT(1,J)
 510      CONTINUE
          do 508 i = 1, nsoL
            fsLd(i) = psol_aerad*u0_aerad !u0*solfx(i)
            alb_toa(i) = fsLu(i)/fsLd(i)
            tsLu = tsLu + fsLu(i)
            tsLd = tsLd + fsLd(i)
 508      continue
          alb_tomi = fupbs(1)/fdownbs(1)
          alb_toai = tsLu/tsLd
!          write(*,*) 'alb_toai',alb_toai
!          if ((alb_toai.gt.1.0) .or. (alb_toai.lt.0.0)) then
!          write(*,*) 'alb_toai crazy!!!',alb_toai
!          write(*,*) 'aeroprof',aeroprof
!          write(*,*) 'w0',w0
!          write(*,*)  'G0',G0
!          write(*,*)  'TAUL',TAUL
!          write(*,*) 'OPD',OPD 
!          write(*,*) 'PSOL_aerad',PSOL_aerad
!         write(*,*) 'u0_aerad',u0_aerad
!          write(*,*) 'alat',alat
!          write(*,*) 'alon',alon
!          write(*,*) 'fnet',fnet
!          CALL ADD
!          write(*,*) 'fnet',fnet          
!          stop
!          endif

!          write(*,*) 'psol_aerad raddtran',psol_aerad
!          write(*,*) 'u0_aerad raddtran',u0_aerad
!         write(*,*) 'alb_tomi',alb_tomi
!          write(*,*) 'alb_toai',alb_toai
!          write(*,*) 'fups(1)',fupbs(1)
!          write(*,*) 'fdownbs2(1)',fdownbs2(1)
!          write(*,*) 'fsLu(1)',fsLu(1)
!          write(*,*) 'fsLd(1)',fsLd(1)
!          write(*,*) 'diff=',fsLu(1)-fsLd(1)
!          write(*,*) 'fnet(1,1)',fnet(1,1)
!          write(*,*) 'fup-fdown',fupbs(1)-fdownbs2(1)
!          write(*,*) 'fnet(1,2)',fnet(1,2)
!          write(*,*) 'fup-fdown',fupbs(2)-fdownbs2(2)
!          stop
!          write(*,*) 'tsLu',tsLu
!          write(*,*) 'tsLd',tsLd
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
!               write(*,*)'IR: DOWN,UP,NET'
          DO 520 L        =  NSOLP+1,NTOTAL
             XIRDOWN = XIRDOWN + DIREC  (L,NLAYER)
             XIRUP   = XIRUP   + DIRECTU(L,NLAYER)
             firu(L-nsol ) = firu( L-nsol ) +    
     &                                directu(L,1)
             do 520 j = 1, nlayer
               fupbi(j) = fupbi(j) + directu(L,j)
               fdownbi(j) = fdownbi(j) + direc  (L,j)
               fnetbi(j) = fnetbi(j) + directu(L,j) - direc(L,j)
!               write(*,*),fdownbi(j),fupbi(j),fnetbi(j)
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

