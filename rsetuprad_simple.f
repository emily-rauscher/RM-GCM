      SUBROUTINE SETUPRAD_SIMPLE
!
!     *********************************************************
!     *  Purpose            :  Defines all constants, and     *
!     *                        calculates pressure averaged   *
!     *                        absorption coefficients.       *
!     * *******************************************************
!
      include 'rcommons.h'
!
! **********************************************************************
!
!           LOCAL DECLARATIONS
!
! **********************************************************************
!     
      REAL bar2kamg,p(NZ),ABSCOEFF(NWAVE),G,WVO,AM
!MTR      DIMENSION AKO3(4,6), AKCO2(6,6), AKH2O(54,6), PJ(6)
      dimension rup_1(NGROUP)
      dimension rhoi(NRAD), dbnds(NRAD+1)
      dimension zbnds(6), pbnds(6), rn2ds(NRAD,6)
!      dimension tauem(NWAVE,6), ssam(NWAVE,6), asmm(NWAVE,6)
      dimension tauem(5,NWAVE), ssam(5,NWAVE), asmm(5,NWAVE)
      dimension temparr(6,NWAVE)
      dimension pbndsm(6)
      integer i1, i2, indorder(5)
!
      logical all_ok
      DATA AVG    / 6.02252E+23  /
      DATA PI    /3.14159265359/

!
! **********************************************************************
!
!            DEFINE CONSTANTS 
!
! **********************************************************************
!
!
!     UNITS ARE (CM**2)/GM
!MTR      DATA (AH2O(I),I=1,77)  /      14*0.0, 0.0000, 0.1965, 9.2460,
!MTR     &       0.0000, 0.1765, 9.2460, 0.0000, 0.0000, 0.2939, 0.5311,
!MTR     &     113.40, 0.0000, 0.0000, 0.3055, 5.2180, 113.00, 0.0000,
!MTR     &       0.3355, 5.5090, 124.90, 3*0.00, 0.0000, 0.3420, 7.1190,
!MTR     &       95.740, 4*0.00,   0.00, 4*0.00, 4*.3012,4*5.021,4*63.17, 
!MTR     &        4*699.1,0.0000, 6.3321, 5.336,  123.4,  7*0.0    /
!
!MTR      DATA (PSH2O(I),I=1,77) /     
!MTR     &       14*0.0, 3*0.54, 3*0.54, 0.0, 4*0.54, 0.0, 4*0.52,     
!MTR     &       4*0.44, 3*0.00, 4*0.62, 5*0.0, 20*0.60,   4*0.60,     
!MTR     &       7*0.0                     /
!
!     ***********************
!     *  DATA FOR INFRARED  *
!     ***********************
!
!MTR      DATA (PSH2O(I),I=78,148)  /   71*0.0                  /
!
!
!     GAUSS ANGLES AND GAUSS WEIGHTS FOR GAUSSIAN INTEGRATION
!     MOMENTS (USE FIRST MOMENT VALUES) N=3
!
      DATA GANGLE  /  0.2123405382, 0.5905331356,     
     &               0.9114120405                       /
      DATA GRATIO/0.4679139346, 0.3607615730, 0.1713244924/
      DATA GWEIGHT /  0.0698269799, 0.2292411064,     
     &                 0.2009319137                       /
!
!
!     ALOS   - LOCSHMIDT'S NUMBER (#/CM**3)
!     AM     - MOLECULAR WEIGHT OF AIR (G/MOL)
!     AVG    - AVAGODROS' NUMBER (#/MOL)
!     G      - GRAVITY (CM/S**2)
!     PI     - PI
!     RGAS   - UNIVERSAL GAS CONSTANT (ERG / MOL K)
!     SCDAY  - NUMBER OF SECONDS IN ONE DAY (S)
!
!      DATA AM     / 28.966       /
      AM= RGAS/R_AIR
      DATA ALOS   / 2.68719E19   /
!MTR      DATA AVG    / 6.02252E+23  /
!      DATA G      / 980.6        /
      G= GA*100.
!MTR      DATA PI     / 3.1415926536 /
      DATA RGAS   / 8.31430E+07  /
      DATA SBK    / 5.6697E-8    /
      DATA SCDAY  / 86400.0      /


!     EPSILON - roundoff error precision
!
      DATA EPSILON / ALMOST_ZERO  /
!
!     EXPMAX - LARGEST (NEGATIVE) EXP ARGUMENT
!
!MTR      DATA EXPMAX  / POWMAX       /
!
!     THIS ROUTINE ASSUMES THAT PRESSURE DOES NOT VARY WITH
!     LOCATION OR TIME.
!
!     CO2MOL - MOLECULAR WEIGHT OF CO2 (G/MOL)
!     O3MOL  - MOLECULAR WEIGHT OF O3 (G/MOL)
!     O2MOL  - MOLECULAR WEIGHT OF O2 (G/MOL)
!
      DATA O2MOL   /  32.         /
!
!
!MTR      DATA (AKH2O(1,I),I=1,6)  / 0.02080, 0.01550, 0.01040, 0.00510,
!MTR     &                            0.00200, 0.00062       /
!
!     CORERAD - RADIUS OF CORE OF AEROSOL PARTICLES
!     COREREAL- REAL PART OF REFRACTIVE INDEX OF CORES
!     COREIMAG- IMAGINARY PART OF REFRACTIVE INDEX OF CORES
!
      DATA CORERAD  / 0.0        /
      DATA COREREAL / 1.25       /
      DATA COREIMAG / 0.5        /
!
!     
!      write(*,*) 'top of pradsetup,, P',P
!      write(*,*) 'top of pradsetup,PR',Pr
!      write(*,*) 'top of pradsetup,Pfull',p_full
!      write(*,*) 'top of pradsetup,paerad',p_aerad 
      AM= RGAS/R_AIR
!      write(*,*)'ABSLW1',ABSLW1
      ABSCOEFF(1)=ABSSW1
      ABSCOEFF(2)=ABSLW1
!      write(*,*)'ABSCOEFF', ABSCOEFF
!      write(*,*) 'Im in rsetuprad_simple'
!     DERIVED PARAMETERS
!
      SQ3     =   SQRT(3.)
      JDBLE   =   2*NLAYER
      JN      =   JDBLE-1
      TPI     =   2.*PI
!      CPCON   =   1.006
      CPCON   =   GASCON/AKAP/1000.  ! Cp in J/gm/K 
      FDEGDAY =   1.0E-4*G*SCDAY/CPCON
!  Open output print files
!
      prtofil = 'carma.p'
      radofil = 'rad.p'


!MTR      open(unit=LUNOPRT,file=prtofil,status='unknown')
!MTR      open(unit=LUNORAD,file=radofil,status='unknown')
!
!     Initialize the aerosol concentrations
!
!MTR ayayay...
!MTR      open(unit=18,file='dndz_0808_lawson.dat',status='unknown',form='formatted')
!MTR      read(18,*) i1, i2
!MTR      read(18,*) dbnds
!      print*, dbnds*1.e4
!MTR      read(18,*) r
!MTR      r = r / 2.
!MTR      rup(:,1) = dbnds(2:NBIN+1) / 2.
!      print*, r*1.e4
!      print*, rup*1.e4
!      read(18,*) zbnds
!!      print*, zbnds
!      read(18,*) pbnds
!      pbnds = pbnds * 1.e3
!      print*, pbnds
!      read(18,*) rn2ds
!!      print*, rn2ds(:,1)
!MTR      close(18)

!      open(unit=19,file='20080418_bext_w0_g.txt',status='unknown',form='formatted')
!      open(unit=19,file='20080418_bext_w0_g_liquid_only.txt',status='unknown',form='formatted')
!      open(unit=19,file='20080427_bext_w0_g.txt',status='unknown',form='formatted')
!      open(unit=19,file='20080426_bext_w0_g_v2.txt',status='unknown',form='formatted')
!      open(unit=19,file='20080418_bext_w0_g_v2.txt',status='unknown',form='formatted')
!      open(unit=19,file='RF31_bext_w0_g_no_ice.txt',status='unknown',form='formatted')
!      open(unit=19,file='RF23_bext_w0_g_no_ice.txt',status='unknown',form='formatted')
!MTR      open(unit=19,file='Simulation1and8_RF31_bext_w0_g.txt',status='unknown',form='formatted')
!      open(unit=19,file='Simulation2_RF31_noIce_bext_w0_g.txt',status='unknown',form='formatted')
!      open(unit=19,file='Simulation3_RF31_noAerosol_bext_w0_g.txt',status='unknown',form='formatted')
!      open(unit=19,file='Simulation4and5_RF23_bext_w0_g.txt',status='unknown',form='formatted')
!      open(unit=19,file='Simulation6_RF23_noIce_bext_w0_g.txt',status='unknown',form='formatted')
!      open(unit=19,file='Simulation7_RF23_noAerosol_bext_w0_g.txt',status='unknown',form='formatted')
!MTR      tauem = 0.
!MTR      ssam = 0.
!MTR      asmm = 0.
!MTR      read(19,*) tauem, ssam, asmm
!...Swap order of variables
      
!      pbndsm = (/ 651.162, 645.037, 640.145, 636.170, 630.406 /)
!      pbndsm = (/ 946.548, 936.719, 929.154, 924.076, 918.341 /)
!      pbndsm = (/ 1946.548, 936.719, 929.154, 924.076, 918.341 /)
!...For RF31
!      pbndsm = (/ 1018.5, 941., 936., 929., 924., 918. /)
!...Hack: use shifted cloud boundaries
!MTR      pbndsm = (/ 783.771, 721.362, 717.383, 711.811, 707.820, 703.022
!MTR/)
!MTR      pbndsm = (/ 1018.5, 721.362, 717.383, 711.811, 707.820, 703.022 /)
!...For RF23
!      pbndsm = (/ 1029.53, 651., 645., 640.5, 636., 632. /)

! change order of vertical levels

!      temparr = tauem
!      do iwave = 1, NWAVE
!        tauem(1,iwave) = temparr(2,iwave)
!        tauem(2,iwave) = temparr(3,iwave)
!        tauem(3,iwave) = temparr(5,iwave)
!        tauem(5,iwave) = temparr(1,iwave)
!      enddo

!      temparr = ssam
!      do iwave = 1, NWAVE
!        ssam(1,iwave) = temparr(2,iwave)
!        ssam(2,iwave) = temparr(3,iwave)
!        ssam(3,iwave) = temparr(5,iwave)
!        ssam(5,iwave) = temparr(1,iwave)
!      enddo

!      temparr = asmm
!      do iwave = 1, NWAVE
!        asmm(1,iwave) = temparr(2,iwave)
!        asmm(2,iwave) = temparr(3,iwave)
!        asmm(3,iwave) = temparr(5,iwave)
!        asmm(5,iwave) = temparr(1,iwave)
!      enddo

      do iwave = 1, NWAVE
        do iz = 1, NVERT
          taua(iwave,iz) = 0.
          taus(iwave,iz) = 0.
          g01(iwave,iz) = 0.
        enddo
        do iz = 2, NVERT
          plev = p_aerad(iz)/1.e3
          do ip = 1, 4
            if( plev < pbndsm(ip) .and. plev >= pbndsm(ip+1) ) then
!MTR              taua(iwave,iz) = tauem(ip,iwave)*dz(iz)
!MTR              taus(iwave,iz) = taua(iwave,iz) * ssam(ip,iwave)
              g01(iwave,iz) = asmm(ip,iwave)
!...dbg:
!      if( iwave == 9 ) print*, iz, ip, plev, tauem(ip,iwave),
!      taua(iwave,iz)
!      if( iz == 13 ) print*, iwave, plev, taus(iwave,iz),
!      taua(iwave,iz), g01(iwave,iz)
!      if( iz >= 1 ) print*, iwave, ip, plev, tauem(ip,iwave),
!      taus(iwave,iz)
            endif
          enddo
        enddo
      enddo
!      print*, sum( taua(9,:) )
!      stop


!MTR      do ig = 1, NGROUP
!MTR        do j = 1, NLAYER
!MTR          do i = 1, NRAD
!MTR            caer(i,j,ig) = 0.
!MTR          enddo
!MTR        enddo
!MTR      enddo

!MTR      ig = 1
!MTR      ibinm = 1
!      do j = 1, NVERT
!        plev = p_aerad(j)
!        do ip = 1, 6
!          if( plev < pbnds(ip) .and. plev >= pbnds(ip+1) ) then
!            do i = ibinm, NRAD
!              caer(i,j,ig) = rn2ds(i,ip)
!            enddo
!!...dbg:
!!     print*, j, ip, plev/1.e3, caer(3,j,ig)
!          endif
!        enddo
!!        print*, j, plev/1.e3, caer(5,j,ig)
!      enddo

!MTR      riwp = 0.
!MTR      do ig = 1, NGROUP
!MTR        do j = 1, NLAYER
!MTR          do i = 1, NRAD
!...increase small crystals
!            if( i <= 5 ) caer(i,j,ig) = caer(i,j,ig)*10.
!            caer(i,j,ig) = 0.
!            caer(i,j,ig) = caer(i,j,ig) *scalef
!MTR            riwp = riwp + 4./3.*PI*rhoi(i)*r(i,ig)**3*caer(i,j,ig)*2.e4
!MTR            caer(i,j,ig) = caer(i,j,ig) * 2.e4
!      if( caer(i,j,ig) > 0. ) then
!        print*, j,i,r(i,ig)*1.e4, caer(i,j,ig), riwp
!      endif
!MTR          enddo
!MTR        enddo
!MTR      enddo
!      print*, 'iwp = ', riwp*1.e4
!
!     Inverse of nprob matrix:  NPROBI(i,1) is first probability
!     interval
!     for wavelength interval i, and NPROBI(i,2) is number of
!     probability
!     intervals in the wavelength interval
!!
!      do i = 1, nwave
!        kount = 0
!        do j = 1, ntotal
!          if (nprob(j).eq.i) then
!            if (kount.eq.0) nprobi(i,1) = j
!            kount = kount + 1
!          endif
!        enddo
!        nprobi(i,2) = kount
!              PLANK(K,J)=T1*T1*T1*T1
!      enddo
!
!     Load wavelengths into interface common block
!
!      do i = 1, NWAVE+1
!        wave_aerad(i) = wave(i)
!      enddo
!
!     Get scalars from interface common block:
!
!     ISL        - do solar calculations when = 1
!     IR         - do infrared calculations when = 1
!     IRS        - do infrared scattering when = 1
!     FLXLIMDIF  - do flux limited diffusion correction = 1 ~MTR
!     UO         - SOLAR ZENITH ANGLE
!     EMISIR     - SURFACE IR EMISSIVITY
!     PTOP       - PRESSURE AT TOP OF MODEL (DYNES/CM**2)
!     PBOT       - PRESSURE AT BOTTOM OF MODEL (DYNES/CM**2)
!     SFC_WIND   - wind speed at 10 m altitude (m/s)
!     SFC_ALB    - surface albedo when fixed
!
      ISL          = 1
      IR           = 1
      IRS          = 1
      FLXLIMDIF    = 1
      EMISIR       = 1.
!      PTOP         = 53.6e3
!      PBOT         = 1000.e3
      PTOP         =p_aerad(1)*10. !P(1)
      PBOT         =p_aerad(NL+1)*10.
!       WRITE(*,*) 'PTOP',PTOP
!       write(*,*) 'PBOT',PBOT
!       write(*,*) 'ISL', ISL
!       write(*,*) 'PLANK',PLANK
     
!      write(*,*) 'p_aerad',p_aerad
    
!...RF31
!      PBOT         = 1018.5e3
!...RF23
!      PBOT         = 1024.5e3

!      if( ifix_sfc_alb .eq. 1 )then
!        SFC_ALB  = sfc_alb_aerad
!      else
!        SFC_WIND = sfc_wind_aerad
!      endif
!      SFC_ALB = 0.07
!      SFC_ALB = 0.9
      ALBEDO_SFC = SWALB !SFC_ALB MTR
!
!     Get atmospheric pressure profile from interface common block
!     [ dyne / cm^2 ]
!
      do k = 1,NZ
!        p(k) = p_aerad(k)
        press(k)=p_aerad(k)*10.
!        write(*,*),'k,p(k)',k,p(k)
      enddo
!      write(*,*) 'line 355'
 

!     The layer thickness in pressure should be the difference between
!     the top and bottom edge of the layer.  Eg. So for layer 3, the layer thickness
!     is the pressure at the boundaries 4 minus at the boundary 3.  
!     These boundary pressures are given in p_full (NZ+1 of them). 
!     The mass of these layers is this pressure thickness divide by G.
!    
!     (NOTE - THE TOP LAYER IS FROM PTOP TO 0, SO AVERAGE = PTOP/2)
!     P_aerad=PRESSURE AT EDGES OF LAYERS (PASCALS) 
!     PRESS - PRESSURE AT EDGE OF LAYER (dyne/cm^2)
!     DPG   - MASS OF LAYER (G / CM**2)!    D PBAR 
!     PBARS - THICKNESS OF LAYER IN PRESSURE (BARS)

!      write(*,*) 'P_aerad in rsetuprad',P_aerad
!      write(*,*) 'G in rsetuprad',G
!      write(*,*) 'NLAYER',NLAYER
!      write(*,*) 'NL',NL
      PRESS=P_aerad*10.
      PBAR(1)  = P_AERAD(2)*1e-6
      DPG(1)= PRESS(1)/G
      DO 45 J  = 2,NLAYER
         PBAR(J)  = (p_aerad(J)-p_aerad(J-1))*1e-5
         DPG(J) = (PRESS(J)-PRESS(J-1)) / G
!         write(*,*),PRESS(J),press(j-1),DPG(j)  
 45    CONTINUE
!        Here we are saying that the top layer has a thickness that
!        extends from the top of the atmosphere (pressure = 0) down to
!        the interface ABOVE the top sigma levels. The atmosphere
!        above the upper most flux boundary (Press(1)) is not explicitly
!        heated or emitting; it simply attenuates the incoming and
!        outgoing radiation by a very small amount.
!        In practice, the small amount of mass is non very signicant,
!        and the absorption coefficient at such rarified pressures and
!        temperatures is probably not very valid. 
!        DPG and PBAR include this mass above the model in index.  Since
!        there are N vertical levels and N+1 layers (i.e. boundaries,
!        interfaces), the arrays of DPG and PBAR have N+1 too account
!        for this cap above the model.  The attenuation through it must
!        still be calculated using these values.

!         DPG(2)=PRESS(1)/G

!       
!         PRESS(NZ) = p_aerad(nz)*10.
!         write(*,*),'PBAR(J)',PBAR
!         write(*,*),'PRESS(J)',PRESS
!         write(*,*),'DPG(J)',DPG
!         write(*,*),'G', G
             
!         print*, 'j p dpg: ', j, press(j)/1.e3, dpg(j-1)
! 45   CONTINUE
!      PBAR(NLAYER)  = P(NVERT)/1.0E6
!      PRESS(NLAYER) = PBOT
!      DPG(NVERT)    = (PRESS(NLAYER)-PRESS(NVERT)) / G
!

!     SKIN TEMPERATURE
!
!  For RF31
!MTR don't understand      tgrnd = 264.673
!  For RF23
!      TT(NLAYER) = TGRND
!      write(*,*)'TGRND',TGRND
!      tgrnd = t_aerad(NZ)
!     AMOUNT OF WATER VAPOR ABOVE MODEL DOMAIN (GM / CM**2)
!     From 1976 U.S. Standard Atmosphere, mid-latitude sounding:
!
!     For z_top = 5 km
!     For z_top = 5 km
!     RDH2O(1) = .13
!
!     For z_top = 2.6 km
!     RDH2O(1) = .64
!
!     For z_top = 19 km
      RDH2O(1)   = 1.6e-4
      RDH2O(1)   = 1.6e-3
!
!     CALCULATE CLOUD AND AEROSOL OPACITY.
!
      DO 100  J             = 1,NLAYER
          DO 50 L           = 1,NTOTAL
               TAUAER(L,J)  = 0.
               TAUCLD(L,J)  = 0.
               WCLD(L,J)    = 0.
               WOL(L,J)     = 0.
               GOL(L,J)     = 0.
               GCLD(L,J)    = 0.
               TAURAY(L,J)    = 0.
 50        CONTINUE
100   CONTINUE
!
!      DO 120 L           =   NSOLP+1,NTOTAL
!         LTEMP(L-NSOLP)  =   NPROB(L) - NSOL
! 120  CONTINUE
!
      X                  =   ALOS/AVG
!
!
!     CALCULATE ABSORPTION COEFFICIENTS
!     HERE WE FIND TAUGAS. IT IS TAUCO2+TAUO2+TAUO3.
!
!MTR      PM             =   PTOP/G

!MTR       write(*,*),'PM',PM
!MTR       
!       DO 305 L       =   1,NTOTAL
!          IF (OPACIR_POWERLAW.eq.0) THEN !MTR Modif to avoid exponent
!                        TAUGAS(L,1)=ABSCOEFF(L)*PM
!                      ELSE IF (OPACIR_POWERLAW.eq.1) THEN
!                        TAUGAS(L,1)=ABSCOEFF(L)*PM
!      &                  *MAX(1E-6,(PRESS(1)/10./OPACIR_REFPRES))
!                      ELSE IF (OPACIR_POWERLAW.eq.2) THEN
!                        TAUGAS(L,1)=ABSCOEFF(L)*PM
!      &                 *MAX(1E-6,(PRESS(1)/10./OPACIR_REFPRES))
!      &                  *MAX(1E-6,(PRESS(1)/10./OPACIR_REFPRES))
!                      ELSE IF (OPACIR_POWERLAW.eq.3) THEN
!                        TAUGAS(L,1)=ABSCOEFF(L)*PM
!      &                  *MAX(1E-6,(PRESS(1)/10./OPACIR_REFPRES))
!      &                  *MAX(1E-6,(PRESS(1)/10./OPACIR_REFPRES))
!      &                  *MAX(1E-6,(PRESS(1)/10./OPACIR_REFPRES))
!                      ELSE
!                        TAUGAS(L,1)=ABSCOEFF(L)*PM
!      &                 *MAX(1E-6,(PRESS(1)/10./OPACIR_REFPRES)**OPACIR_POWERLAW)
!           ENDIF     
!  305  CONTINUE
!
      DO 308   J     =   1,NLAYER
         PM          =   DPG(J) 
!         write(*,*) 'PM', PM 
!         write(*,*) 'ABSCOEFF',ABSCOEFF
!         write(*,*) 'OPACIR_REFPRES',OPACIR_REFPRES
!         write(*,*) 'PRESS(J)',PRESS(J)   
         DO 306 L    =   1,NTOTAL
         IF (OPACIR_POWERLAW.eq.0) THEN !MTR Modif to avoid exponent
                       TAUGAS(L,J)=ABSCOEFF(L)*PM
                     ELSE IF (OPACIR_POWERLAW.eq.1) THEN
                       TAUGAS(L,J)=ABSCOEFF(L)*PM
     &                  *MAX(1E-6,(PRESS(J)/10./OPACIR_REFPRES))
                     ELSE IF (OPACIR_POWERLAW.eq.2) THEN
                       TAUGAS(L,J)=ABSCOEFF(L)*PM
     &                  *MAX(1E-6,(PRESS(J)/10./OPACIR_REFPRES))
     &                  *MAX(1E-6,(PRESS(J)/10./OPACIR_REFPRES))
                     ELSE IF (OPACIR_POWERLAW.eq.3) THEN
                       TAUGAS(L,J)=ABSCOEFF(L)*PM
     &                  *MAX(1E-6,(PRESS(J)/10./OPACIR_REFPRES))
     &                  *MAX(1E-6,(PRESS(J)/10./OPACIR_REFPRES))
     &                  *MAX(1E-6,(PRESS(J)/10./OPACIR_REFPRES))
                     ELSE
                       TAUGAS(L,J)=ABSCOEFF(L)*PM
     &                  *MAX(1E-6,(PRESS(J)/10./OPACIR_REFPRES)**OPACIR_POWERLAW)
          ENDIF
!      write(*,*)'L,J,TAUGAS',L,J,TAUGAS(L,J)
 306  CONTINUE
 308  CONTINUE

      DO  L           =   NSOLP+1,NTOTAL
          TAUCONST(L)=ABSCOEFF(L)/GA/100.
      ENDDO

  
!     WAVE MUST BE IN MICRONS
!     CALCULATE RAYLEIGH OPTICAL DEPTH PARAMETERS.
!
!      write(*,*) 'G', G
!      write(*,*) 'AM',AM
!      write(*,*) 'AVG',AVG
      bar2kamg=AVG*1.e5/(2.686763e25*G/100.*AM*1e-3)*1e-3
!      bar2kamg= 2.686763e25*G*AM/(1e5*1e-3*AVG)
!        write(*,*) 'bar2kamg',bar2kamg
    
!      write(*,*) 'kmamg2bar',kmamg2bar
!it is derived by the relation:
!km-amagats= Pressure * Avogadro's# * mole fraction/  (Loshchmidt's# *
!gravity* molecular weight)
!where pressure is in Pa, Avogadro's=6.0221367*10^23 mol^-1
!Lo=2.6867630*10^25 m^-3, gravity= 24.40 m s^-2
!molec. weight for mole fraction .86 H2 and .136 He (von Zahn)
! 2.27 *10^-3 kg/mole

      DO 310 L      =   1,NTOTAL
!
!MTR          WVO       =    WAVE(NPROB(L))
!MTR          TAURAY(L) =   (8.46E-9/WVO**4) *     &
!MTR                        ( 1.+0.0113/WVO**2+0.00013/WVO**4 )
              WVO = 0.6  !Hardwire for H2!! MTR
              RAYPERBAR(L) =   (bar2kamg*0.000219/WVO/WVO/WVO/WVO) *     
     &           ( 1.+0.0157248/WVO/WVO+0.0001978/WVO/WVO/WVO/WVO)*0.
              
!              write(*,*) 'RAYPERBAR(L)', RAYPERBAR(L)
!              write(*,*) 'PBAR',PBAR
 310  CONTINUE
!     WE DO NOT INCLUDE RAYLEIGH SCATTERING IN INFRARED
!
      DO 330 J          =   1,NLAYER
         DO 320 L       =   1,NTOTAL
           if( L .LE. NSOLP )then
!             PARAY(L,J+1) = TAURAY(L)*DPG(J)*G ; 
              TAURAY(L,J) = RAYPERBAR(L)*PBAR(J) !PER BAR X LAYER THICKNESS IN BAR
           else
             TAURAY(L,J) = 0.0
           endif
 320     CONTINUE
!
!         DO 325 L       =   1,NLAYER
!            if( L .LE. NSOLP )then
!              TAURAY(L,1) = RAYPERBAR(L)*PTOP
!              TAURAY(L,1) = RAYPERBAR(L)*PBAR(1)
!            else
!              TAURAY(L,1) = 0.0
!            endif
! 325     CONTINUE
!
 330  CONTINUE
!
! **********************************************************************
!
!            CALCULATE THE AEROSOL EXTINCTION CROSS SECTIONS
!
! **********************************************************************
!
!     Get <is_grp_ice> and radius grid from interface common block
!     and calculate cross-sectional area for each bin.
!
!
!
!  extra absorption of solar radiation 
!  choices are I_NO_SOOT, I_CORE_FIX, I_CORE_FRAC, I_SHELL_FRAC
!
!
!      CALCULATE THE CENTER OF THE WAVELENGTH INTERVAL OF AN IR INTERVAL
!
!
!        RDQEXT(I,ig,L)   =   0.0
!        QSCAT(I,ig,L)    =   0.0
!        QBRQS(I,ig,L)    =   0.0
!
!        DO 104 J         =   1,6
!
!  Core no bigger than particle
!
!         corerad_safe = min( rr, corerad )
! 
!         CALL MIESS(RR,REAL,TMAG,thetd,n_thetd,QEXTD,QSCATD,CTBRQS,  &
!                    CORERAD_safe,COREREAL,COREIMAG,WVNO)
! 
!         RDQEXT(I,ig,L)     =   RDQEXT(I,ig,L)+QEXTD/6.
!         QSCAT(I,ig,L)      =   QSCAT(I,ig,L)+QSCATD/6.
!         QBRQS(I,ig,L)      =   QBRQS(I,ig,L)+CTBRQS/6.
!         RR                 =   RR+DDR
! 
! 104    CONTINUE
!
! 107   CONTINUE
! 110  CONTINUE
!
!      enddo      ! ig=1,NGROUP
!
!      endif   
!
!      if ( do_mie ) then

!     Write extinction and scattering coefficients to data file
!
!        open(LUNMIE,file='mie.data',form='formatted')
!
!        write(LUNMIE,*) NWAVE,NRAD,NGROUP
!        write(LUNMIE,*) (r(1,ig),ig=1,NGROUP)
!        write(LUNMIE,*) (rup(1,ig),ig=1,NGROUP)

!        do ig = 1,NGROUP
!         do i = 1,NRAD
!          do L = 1,NWAVE
!            write(LUNMIE,*) rdqext(i,ig,L), qscat(i,ig,L), qbrqs(i,ig,L)
!          enddo
!         enddo
!        enddo
!
!c        close(LUNMIE)
!
!      endif   
!
!     Write some values to print file
!
!      if( myproc .eq. 0 )then
!        do ig = 1, NGROUP
!          WRITE (LUNORAD,500) ig
!          DO I = 1, NRAD
!            sizparm6 = 2.*pi*r(i,ig)/(wave(6)*1.e-4)
!            sizparm24 = 2.*pi*r(i,ig)/(wave(24)*1.e-4)
!            WRITE(LUNORAD,505) I,R(I,ig),rdqext(i,ig,6),sizparm6,
!     1         rdqext(i,ig,24),sizparm24
!          ENDDO
!        enddo
!      endif
! 500  FORMAT(/," SETUPRAD: igroup = ",i4,//,     &
!             "   i     r(cm)   rdqext(6)      x6",     &
!             "   rdqext(24)      x24",/)
! 505  FORMAT(I4,5(1PE11.2))
!
! *********************************************************************
!
!                              CHECK SUM OF WEIGHTS
!
! **********************************************************************
!

!MTR      SUM0           =   0.0
!MTR       SUM1          =   0.0
!MTR       SUM2          =   0.0
!MTR       DO 340 L      =   1,NSOLP
!MTR          SUM0        =   SUM0+WEIGHT(L)
!MTR  340  CONTINUE
!MTR       DO 350 L      =   NSOLP+1,NTOTAL
!MTR          SUM1       =   SUM1+WEIGHT(L)
!MTR  350  CONTINUE
!MTR       SUM2          =   SUM0+SUM1
!
!      IF ( ABS(NWAVE-SUM2) .GT. 1.E-3 ) WRITE(LUNORAD,355)
!      SUM,SUM1,SUM2
!
!MTR  355  FORMAT(//,"SETUPRAD: Error in weights ",/,     &
!MTR               " Sum of weights for solar = ",1PE15.5,/,     &
!MTR               " sum of weights for ir = ",1PE15.5,/,     &
!MTR               " total sum = ",1PE15.5)
!
       DO 360 L   =   1,NSOLP
          SOL(L)  = PSOL_aerad !  SOLC_IN  !SOLFX(NPROB(L)) * WEIGHT(L)
  360  CONTINUE
! *********************************************************************
!
!     COMPUTE PLANCK FUNCTION TABLE. WAVE IS IN UNITS OF MICRONS.
!
! **********************************************************************
!
!
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
!  ;COMMENTING ALL THIS OUT TO SIMPLY COMPUTE SIGMA T^4 FOR BLACKBODY EMMISION
!   DIRECTLY LATER, RATHER THAN PRODUCING A TABLE OF VALUES USING THE NLOW TO
!   NHIGH VALUES IN THE GLOBRAD.H

!MTR      DO 380 J  =   1,NCOUNT
!MTR         JJ     =   NLOW+J
!MTR         T1     =   0.01 * FLOAT(JJ)

!MTR         if( ibeyond_spectrum .eq. 1 )then

!MTR           plank(nwave+1-nsol,j) = t1**4
!MTR           DO I =   NSOL+2,NWAVE
!MTR              K =   I-NSOL
!MTR              V =   1.438E4  /  WAVE(I)
!MTR            CALL PLNK(V,T1,PLANK(K,J))
!MTR           ENDDO

!MTR         else

!MTR           DO I =   NSOL+1,NWAVE+1
!MTR              K =   I-NSOL
!MTR              V =   1.438E4  /  WAVE(I)
!MTR!              CALL PLNK(V,T1,PLANK(K,J))
!MTR           ENDDO

!MTR         endif

!MTR 380  CONTINUE
!
!MTR      DO 410 J   =   1,NCOUNT

!MTR         if( ibeyond_spectrum .eq. 1 )then

!MTR           plank(1,j) = plank(2,j)*sbk/pi
!MTR           DO L  =   NSOL+2,NWAVE
!MTR              K  =   L-NSOL
!MTR              PLANK(K,J) = (PLANK(K+1,J)-PLANK(K,J))*SBK/PI
!MTR           ENDDO

!MTR         else

!MTR           DO L  =   NSOL+1,NWAVE
!MTR              K  =   L-NSOL
!MTR              PLANK(K,J) = (PLANK(K+1,J)-PLANK(K,J))*SBK/PI
!MTR           ENDDO

!MTR         endif

!MTR 410  CONTINUE

      RETURN
      END
