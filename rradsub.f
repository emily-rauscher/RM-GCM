      subroutine radsub(iffirst, deltaz,pr, p_pass,t_pass,qh2o_pass, 
     &                radheat,htlw,htsw,rfluxes,alat,alon,KOUNT,ITSPD)

!     iffirst is just the indicator for numbering and runs the setup
!     deltaz--the layer thickness in meters
!     p_pass--the layer boundary pressures in pascal (NL+1)
!     t_pass--mid-layer temperatures in K and one bottom boundary temp
!     both p_ and t_ pass begin at the top and go down.

      include 'rcommons.h'
      PARAMETER(PI2=2.0*3.14159265359) 
      character*16 clearheat_file_prefix
      character*4 clearheat_file_suffix
      character*22 clearheat_file
      character*2 char_ihour
      integer iffirst
      real t_pass(NZ)
      real p_pass(NZ)
      real deltaz(NZ)
      real qh2o_pass(NZ)
      real radheat(NZ)
      real heats_aerad_tot(NZ), heati_aerad_tot(NZ), radheat_tot(NZ)
      real wave_pass(1) !wave_pass(45)
      real cheats(NZ), cheati(NZ), cheat(NZ)
      real htlw(NZ),htsw(NZ)
      real pbot_pass, ptop_pass,PSOL,PSOL_aerad
      real rfluxes(nwave,2,2)
      REAL SSLON,SSLAT  ! ER:
      REAL DLENGTH  ! ER: half-length of solar dayinteger ibinmin
      real PI2
 582       FORMAT(I4,5(F12.3))
      clearheat_file_prefix = 'clear_heat_rf31_'
      clearheat_file_prefix = 'clear_heat_rf23_'
      clearheat_file_suffix = '.dat'
      ibinm = ibinmin
      ifsetup = 0
      if( iffirst == 1 ) ifsetup = 1
!      GRAV = 980.6d+0
!      RGAS = 8.31430d+07
!      WTMOL_AIR = 28.966d+0
!      R_AIR = RGAS / WTMOL_AIR
!MTR      use phyiscal_constants

!@ Keep an Eye on this, Mike
      if_diurnal = 0

      heats_aerad_tot = 0.
      heati_aerad_tot = 0.
!
!  Reverse the vertical index for radiation code
!  MTR NOT SURE IF VERTICAL INVERSION IS NECESSARY...

!MTR     do iz = 1, NZ

!MTR     k = NZ + 1 - iz
  
!MTR        t_aerad(iz) = t_pass(k)
!MTR        p_aerad(iz) = p_pass(k)
!        qv_aerad(iz) = qh2o_pass(k)
!MTR        dz(iz) = deltaz(k)
!        print*, iz, k, t_aerad(iz), p_aerad(iz)/1.e3, &
!             deltaz(iz)/1.e5, qv_aerad(iz)

!      enddo

       t_aerad=t_pass
       p_aerad=p_pass
      
       dz=deltaz

!      u0_aerad = 0.78
!      u0_aerad = 0.318
!      u0_aerad = 0.448
!      u0_aerad = 0.349  ! Barrow Alaska, 080427, 0z
!      u0_aerad = 0.4798 ! RF23, (72:21:37.44, 154:51:39.6, 18 Apr, 14:52 local)
!      u0_aerad = 0.4481 ! RF31, (72:21:37.44, 154:51:39.6, 26 Apr, 17:05 local)
!!      u0_aerad = 0.01
!      u0_aerad= 
      ir_above_aerad = 0
      tabove_aerad = 0


!...Parameters for dirunal variation of u0

!MTR      if( if_diurnal == 1 ) then
!MTR
!MTR        iday = 117  ! RF31
!MTR        iday = 109  ! RF23
!MTR        rlat = 72.3604
!MTR        rad_start = 0. * 3600.
!MTR
!MTR        saz = 2. * PI / 365. * iday
!MTR
!MTR        declin = 0.006918 - 0.399912*cos(saz)    +0.070257*sin(saz)  &
!MTR                          - 0.006758*cos(2.*saz) +0.000907*sin(2.*saz)
!MTR&
!MTR                          - 0.002697*cos(3.*saz) +0.001480*sin(3.*saz)
!MTR
!MTR        zsin = sin(declin) * sin( rlat * PI/180. )
!MTR       zcos = cos(declin) * cos( rlat * PI/180. )

!MTR        time = 1.*SCDAY/24.
!MTR        sun_angle = PI + ( time + rad_start )*2.*PI/SCDAY
!MTR        u0_aerad = zsin + zcos*cos(sun_angle)
!MTR        u0_aerad = max( 0.*ONE, u0 )
!        print*, itime, time/SCDAY*24., sun_angle, u0

!      endif

!@ The following lines of code are taken from cnikos and may require adjustment
C ER modif for non-synchronous orbit

      IF (PORB.NE.0) THEN
!         write(*,*) 'PORB.NE.0'
!         write(*,*) 'PORB',PORB
!         write(*,*) 'KOUNT',KOUNT
!         write(*,*) 'ITSPD',ITSPD
         SSLON=(1./PORB-1.)*KOUNT*360./ITSPD
C         SSLON=ALON-SSLON
         SSLON=MOD(SSLON,360.)
      ELSE
         SSLON=0.  ! substellar longitude
C         SSLON=ALON  !local longitudinal angle to star, in degrees
      ENDIF
C ER modif for non-zero obliquity
      IF (OBLIQ.EQ.0) THEN
         SSLAT=0.  ! substellar latitude
C        SSLAT=ALAT  !local latitudinal angle to star, in degrees
         DLENGTH=PI/2.
      ELSE
         SSLAT=ASIN(SIN(OBLIQ*PI/180.)
     +        *SIN(PI2*KOUNT/ITSPD/PORB))*180./PI
C        SSLAT=ALAT-SSLAT
         IF (SSLAT.GT.0) THEN
            IF (ALAT.GT.90.-SSLAT) THEN
               DLENGTH=PI
            ELSEIF (ALAT.LT.-90.+SSLAT) THEN
               DLENGTH=0.
            ELSE
               DLENGTH=ACOS(-1.*TAN(ALAT/360.*PI2)*TAN(SSLAT/360.*PI2))
            ENDIF
         ELSEIF (ALAT.LT.-90.-SSLAT) THEN
            DLENGTH=PI
         ELSEIF (ALAT.GT.90+SSLAT) THEN
            DLENGTH=0.
         ELSE
            DLENGTH=ACOS(-1.*TAN(ALAT/360.*PI2)*TAN(SSLAT/360.*PI2))
         ENDIF
      ENDIF

C Setup SW code                                                           
      IF (LBIN) THEN
        call BinaryFlux(SOLC,KOUNT,ITSPD)
        SOLC=SOLC*(1.0-TOAALB)
!        WRITE(88,*) SOLC
      ELSE
        SOLC=SOLC_IN*(1.0-TOAALB)
      ENDIF
!      LDIUR=.TRUE.   !Diurnally averaged if false                        
C      LDIUR=.FALSE.   !Diurnally averaged if false
C                                                                         
C       YCLOCK=3.14159      !time of day in radians                        
!       YCLOCK= (DOY-REAL(INT(DOY))*2.*3.14159                            
C                                                                         
!      CALL SOLANG(LDIUR,DOY,YCLOCK,ALAT,ALON,AMU0,RDAYL,CDISSEM)          
C         
C     globally averaged solar constant, vertical rays
      AMU0=1.0
      PSOL=SOLC/4.
      IF(.NOT.L1DZENITH) THEN
         DDAY=FORCE1DDAYS
CC       DDAY is rampup time from uniform to full heating
C        (1D for DDAYs, then linear to full in another DDAYs)
         IF(DAY.GT.DDAY) THEN
            DFAC=MIN(1.0,(DAY - DDAY)/DDAY)
            IF(.NOT.LDIUR) THEN
CC Modif for hemispheric forcing of Hot Jupiters
!              write(*,*) 'AMU0',AMU0
!              write(*,*) 'DFAC',DFAC
!              write(*,*) 'ALAT',ALAT
!              write(*,*) 'SSLAT', SSLAT
!              write(*,*) 'PI2', PI2
!              write(*,*) 'ALON',ALON
!              write(*,*) 'PI2',PI2
!              write(*,*) 'SSLON',SSLON
!              write(*,*) (1.0-DFAC)*AMU0
!              write(*,*) DFAC *SIN(ALAT/360.*PI2)*SIN(SSLAT/360.*PI2)
!              write(*,*) COS(ALAT/360.*PI2)*COS(SSLAT/360.*PI2)
!              write(*,*) COS((ALON-SSLON)/360.*PI2)
              
               AMU0=(1.0-DFAC)*AMU0
     &              +DFAC*MAX(0.0,SIN(ALAT/360.*PI2)*SIN(SSLAT/360.*PI2)
     &                           +COS(ALAT/360.*PI2)*COS(SSLAT/360.*PI2)
     &                           *COS((ALON-SSLON)/360.*PI2))
               PSOL=(1.0-DFAC)*PSOL + DFAC*SOLC
            ELSE
CC ER modif for diurnal forcing (a la Liu & Schneider 2010)
C               PSOL=(1.0-DFAC)*PSOL + DFAC*SOLC/PI*COS(ALAT/360.*PI2)
C ER modif for non-zero obliquity
               PSOL=(1.0-DFAC)*PSOL
     &              +DFAC*SOLC/PI*
     &              (SIN(ALAT/360.*PI2)*SIN(SSLAT/360.*PI2)*DLENGTH
     &              +COS(ALAT/360.*PI2)*COS(SSLAT/360.*PI2)*SIN(DLENGTH)
     &              )
            ENDIF
         ENDIF
      ENDIF    
!End CNIKOS lines
      u0_aerad = max( 0.*ONE, AMU0 )  
      PSOL_aerad=PSOL
      do_mie_aerad = .false.
         IF (AMU0.GT.0) THEN
           isl_aerad=1
         ELSE 
           isl_aerad=0
         ENDIF
      ir_aerad = 1
      ntime = 1
      if( if_diurnal == 1 ) ntime = 24

      itime1 = 12
      do itime = 1, ntime
!      do itime = itime1, itime1
!        write(char_ihour,'(i2.2)') itime
!        clearheat_file = clearheat_file_prefix // char_ihour //clearheat_file_suffix
!        print*, clearheat_file

!MTR        if( if_diurnal == 1 ) then
!MTR          time = (itime-1.)*SCDAY/24.
!MTR          sun_angle = PI + ( time + rad_start )*2.*PI/SCDAY
!MTR          u0_aerad = zsin + zcos*cos(sun_angle)
!MTR          u0_aerad = max( 1.e-5*ONE, u0_aerad )
!print*, time, sun_angle, u0_aerad
!stop
!MTR        endif
        call setuprad_simple
!        write(*,*) 'called setuprad'
        pc_aerad = 0.
!          write(*,*) 'not calling radtran'
        call radtran
!          write(*,*) 'called radtran'
!        write(*,*) 'RSFX',RSFX
!      print*, 'TOA fsolu fsold fsoln cfsol = ', fupbs(1), fdownbs(1),
!      &
!             fdownbs(1)-fupbs(1), fdownbs(1)-fupbs(1)-370.48325
!      print*, 'TOA firu fird firn cfir = ', fupbi(1), fdownbi(1),   &
!             fdownbi(1)-fupbi(1), fdownbi(1)-fupbi(1)+301.20316
!      print*, 'opd = ', uopd(9,nlayer)
!      print*, 'iwp tau cfnet = ', riwp*1.e4, uopd(9,nlayer)-0.1111,   &
!      print*, riwp*1.e4, uopd(9,nlayer)-0.1111,   &
!              fdownbs(1)-fupbs(1)-370.62216,
!              fdownbi(1)-fupbi(1)+259.67572,  &
!              fdownbs(1)-fupbs(1)-370.62216 +  &
!              fdownbi(1)-fupbi(1)+259.67572

!        call radout

!        WRITE(6,560)
! 560    FORMAT(" RADOUT:",/,      &
!               " j   p(j)    press(j)    t(j) ",      &
!               "     tt(j)  rdh2o(j)  ctot       firu     fird  ",
!               &
!               "    fsLu     fsLd   ")
!!
!        DO 565 J = 1, NVERT
!           ctot = 0.
!           do ig = 1, NGROUP
!             do i = 1, NRAD
!               ctot = ctot + caer(i,j,ig)
!             enddo
!           enddo
!           WRITE(6,562) J,P(J),PRESS(J),T(J),TT(J),      &
!                        RDH2O(J),Ctot,fupbi(j),fdownbi(j),      &
!                        fupbs(j),fdownbs(j)
! 562       FORMAT(I3,11(1PE9.2))
! 565    CONTINUE
!        stop
!
!...Read in clear-sky radiative heating rates
!
        cheats = 0.
        cheati = 0.
        cheat = 0.
!      open(unit=19,file='clear_heat_rf23_1600.dat',status='unknown',form='formatted')
!      open(unit=19,file='clear_heat_rf23_noice.dat',status='unknown',form='formatted')
!      open(unit=19,file='clear_heat_rf31_sza.dat',status='unknown',form='formatted')
!      open(unit=19,file='clear_heat_rf31_alb.dat',status='unknown',form='formatted')
!      open(unit=19,file='clear_heat_rf23_alb.dat',status='unknown',form='formatted')
!      open(unit=19,file='clear_heat_rf31_constalb.dat',status='unknown',
!      &
!           form='formatted')
!MTR        if( if_diurnal == 0 ) then
!      open(unit=19,file='clear_heat_rf31_varalb.dat',status='unknown',
!      &
!           form='formatted')
!MTR      open(unit=19,file='clearheat_rf31.dat',status='unknown', 
!MTR     &      form='formatted')
!          open(unit=19,file='clearheat_rf23.dat',status='unknown', &
!               form='formatted')
!          open(unit=19,file='clearheat_rf23_oceanalbedo.dat',status='unknown',
!          &
!               form='formatted')
!          open(unit=19,file='clearheat_rf31_oceanalbedo.dat',status='unknown',
!          &
!               form='formatted')
!MTR        endif

!MTR
!MTR        if( if_diurnal == 1 ) then
!MTR          open(unit=19,file=clearheat_file,status='unknown', 
!MTR     &          form='formatted')
!MTR        endif
!MTR        do i = 1, NZ
!MTR          read(19,*) i1, r1, r2, r3, r4, r5
!MTR          cheats(i) = r3
!MTR          cheati(i) = r4
!MTR          cheat(i) = r5
!          print*, r3, r4, r5
!MTR        enddo
!MTR        close(19)
!
!  Calculate radiative heating rate
!
!      print*, 'iffirst = ', iffirst
!      print*, 'iz   p   T   hs   hi   hnet'
!print*, clearheat_file
!print*, cheati(1:5)

!        write(*,*) 'heats_aerad',heats_aerad
!        write(*,*) 'heati_aerad',heati_aerad
        do iz = 1,NZ
          jz = NZ + 1 - iz
!          print*, heati_aerad(jz)*SCDAY
!if( iz < 6 ) print*, iz, jz, heati_aerad(jz)*SCDAY - cheati(iz)
          radheat(iz) = heats_aerad(jz) + heati_aerad(jz)
          heats_aerad_tot(iz) = heats_aerad_tot(iz) +   
     &                           heats_aerad(jz)*SCDAY - cheats(iz)
          heati_aerad_tot(iz) = heati_aerad_tot(iz) +   
     &                           heati_aerad(jz)*SCDAY - cheati(iz)
          radheat_tot(iz) = radheat_tot(iz) +   
     &                       heats_aerad(jz)*SCDAY - cheats(iz) +  
     &                       heati_aerad(jz)*SCDAY - cheati(iz)
!       print*,'iz jz p T hs hi heat = ', iz, jz, p_pass(iz)/1.e3,  &
!             t_pass(iz), heats_aerad(jz)*24.*3600.,  &
!             heati_aerad(jz)*24.*3600.,  &
!             radheat(iz)*24.*3600.
!       print*,'iz p T hi heat = ', iz, p_pass(iz)/1.e3,  &
!        WRITE(6,582)  iz, p_pass(iz)/1.e3,  &
!        WRITE(19,582)  iz, p_pass(iz)/1.e3,  &
!             t_pass(iz),  &
!             heats_aerad(jz)*24.*3600.-cheats(iz),  &
!             heati_aerad(jz)*24.*3600.-cheati(iz),  &
!             radheat(iz)*24.*3600.-cheat(iz)
!             heats_aerad(jz)*24.*3600.,  &
!             heati_aerad(jz)*24.*3600.,  &
!             radheat(iz)*24.*3600.
        enddo
!         write(*,*)'heats_SW'
!         do iz=1,NZ
!         write(*,*),heats_aerad(iz)
!         enddo
!         write(*,*)'heats_LW'
!         do iz=1,NZ
!         write(*,*),heati_aerad(iz)
!         enddo
!MTR        close(19)

      enddo

      if( if_diurnal == 1 ) then
        heats_aerad_tot = heats_aerad_tot / ntime
        heati_aerad_tot = heati_aerad_tot / ntime
        radheat_tot = radheat_tot / ntime
      endif

!MTR      do iz = 1,NZ
!MTR        WRITE(6,582)  iz, p_pass(iz) !/1.e3,  
!MTR     &        t_pass(iz),  
!MTR     &        heats_aerad_tot(iz),  
!MTR     &        heati_aerad_tot(iz),  
!MTR     &        radheat_tot(iz)
!MTR      enddo
      htlw=heati_aerad_tot
      htsw=heats_aerad_tot
      rfluxes=rfluxes_aerad
      return

      end

