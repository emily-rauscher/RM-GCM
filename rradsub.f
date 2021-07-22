      subroutine radsub(iffirst,pr, p_pass,t_pass,qh2o_pass, 
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
!      real deltaz(NZ)
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
      if( iffirst.eq. 1 ) ifsetup = 1


      if_diurnal = 0
      heats_aerad_tot = 0.
      heati_aerad_tot = 0.

       t_aerad=t_pass
       p_aerad=p_pass


      ir_above_aerad = 0
      tabove_aerad = 0


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
!      write(*,*) 'solc raddsub',SOLC        
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
      if ((AMU0.gt.0) .and. (AMU0.lt.1e-6)) THEN
       AMU0 = 0.0
      endif
      

!     write(*,*) 'AMU0',AMU0
!      u0_aerad = max( 0.*ONE, AMU0 )  
      u0_aerad = max(0*ONE, AMU0 ) 
      PSOL_aerad=PSOL




      do_mie_aerad = .false.
!         DAY/NIGHT SW CONDITIONAL
!         IF (AMU0.GT.0) THEN
         IF (AMU0.GT.EPSILON) THEN
           isl_aerad=1
         ELSE 
           isl_aerad=0
         ENDIF
      ir_aerad = 1
      ntime = 1
      if( if_diurnal.eq.1 ) ntime = 24

      itime1 = 12
      do itime = 1, ntime
        call setuprad_simple
        pc_aerad = 0.
        call radtran

!
!...Read in clear-sky radiative heating rates
!
        cheats = 0.
        cheati = 0.
        cheat = 0.


        do iz = 1,NZ
          jz = NZ + 1 - iz

          radheat(iz) = heats_aerad(jz) + heati_aerad(jz)
          heats_aerad_tot(iz) = heats_aerad_tot(iz) +   
     &                           heats_aerad(jz)*SCDAY - cheats(iz)
          heati_aerad_tot(iz) = heati_aerad_tot(iz) +   
     &                           heati_aerad(jz)*SCDAY - cheati(iz)
          radheat_tot(iz) = radheat_tot(iz) +   
     &                       heats_aerad(jz)*SCDAY - cheats(iz) +  
     &                       heati_aerad(jz)*SCDAY - cheati(iz)

        enddo

      enddo

      if( if_diurnal.eq. 1 ) then
        heats_aerad_tot = heats_aerad_tot / ntime
        heati_aerad_tot = heati_aerad_tot / ntime
        radheat_tot = radheat_tot / ntime
      endif


      htlw=heati_aerad_tot
      htsw=heats_aerad_tot
      rfluxes=rfluxes_aerad
      return

      end

