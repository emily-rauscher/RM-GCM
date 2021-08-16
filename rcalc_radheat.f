!***********************************************************************  
!*                         SUBROUTINE CALC_RADHEAT                     *  
!*********************************************************************** 
      SUBROUTINE CALC_RADHEAT(pr,t,prflux,alat1,alon,htlw,htsw,
     $             DOY,cf,ic,rfluxes,swalb,kount,itspd)

!...Calculate radiative heating rate profiles and corresponding vertical
!...wind speed.

!.. pr is an array of NL+1, in pascals, starting at the center of the top layer
! and going down to the center of the bottom layer, with an extra
! element for the surface

! T is an array of Temperature, NL+1, in Kelvin, starting at the top
! mid-layer down to the bottom mid layer with an extra index for the
! base boundary temperature.
!  PR and T are passed to P_full and T_full
! NZ is NL+1

!Prflux is an array of NL+1, in pascals of pressurs at the edge of
!layers, starting at the top of the atmosphere and working down. The
!first is zero and the bottom is the bottom pressure with the same value
!as the bottom pressure in Pr (i.e. P0 adjusted for dynamics).

!      use physical_constants

       include 'rcommons.h'

       real, parameter :: BK = 1.38054e-16
       real, parameter :: L      = 2.5e10
       
      REAL PR(NZ),T(NZ),Cpd,PRFLUX(NZ)
      real, dimension(NZ) :: p_full, t_full, radheat
      real, dimension(NZ) :: z, htsw,htlw
      real, dimension(45) :: wave_pass
      real, dimension(NWAVE,2,2) :: rfluxes

      t_full=t
      tgrnd=t(NZ)
      p_full=prflux  ! PRFLUX is the layer edge pressures
      rfluxes=fluxes 
      iffirst = 1
      R_AIR = GASCON*10000.
      Cpd = GASCON/AKAP !1.004e+7
      PREF = P0
      GRAV = GA *100.
      eps= Rd/Rv
      RdCp = Rd/Cpd
      player=pr

      call radsub(iffirst, pr,p_full,t_full,qh2o_full, radheat,htlw,htsw,rfluxes,alat1,alon,KOUNT,ITSPD)

      iffirst = 0
      

      IF ((LFLUXDIAG).AND.(KOUNTP-KOUTP.LT.NTSTEP_IN)) THEN !-1)) THEN
         WRITE(63,*) 'LATITUDE, LONGITUDE:',ALAT1,ALON
         WRITE(63,2010)'1)P(Bars)',
     $        'FLUXES(Wm-2): 2)LW UP','3)LW DOWN',         
     $        '4)LW NET',
     $        '5)SW UP','6)SW DOWN',
     $        '7)SW NET'
 2010    FORMAT(1X,A9,1X,A21,4X,A9,4X,A8,4x,A7,4x,A9,4x,A8)                               

C     ER Modif: output pressures in bar instead of mbar
         DO ILAY=NLAYER,1,-1                                                  
          WRITE(63,2013),P_FULL(NLAYER-ILAY+1)*1e-5,FIR_UP_AERAD(ILAY),                            
     $    FIR_DN_AERAD(ILAY),FIR_NET_AERAD(ILAY),
     $    FSL_UP_AERAD(ILAY),FSL_DN_AERAD(ILAY),FSL_NET_AERAD(ILAY)                         
     $
 2013       FORMAT(2X,F12.6,3X,E12.5,3X,E12.5,3X,E12.5,3X,E12.5         
     $             ,3X,E12.5,3X,E12.5)

      END DO
!
         WRITE(63,2023)'PRESSURE (Bars)','HEATING RATES: SW (K/DAY)'          
     $        ,'LW (K/DAY)','TOTAL (K/DAY)'                                                 
 2023    FORMAT(2X,A15,2X,A25,4x,A10,4x,A13)                                        
         DO LHT=NLAYER,1,-1                                                  
         WRITE(63,2020) PR(NLAYER-LHT+1)*1e-5,HTSW(LHT),HTLW(LHT),
     $          HTSW(LHT)+HTLW(LHT)
     $  
 2020       FORMAT(2X,F12.6,19X,E12.5,3X,E12.5,3x,E12.5)   
         END DO
         WRITE(63,*)

!  HERE WE WRITE TO FILE 62 ADDTIONAL RADIATIVE TRANSFER BY PRODUCTS
         WRITE(62,*)'LATITUDE, LONGITUDE:',ALAT1,ALON
         WRITE(62,*)'  Cosine of the incidencd angle, mu0:',U0
         write(62,*)' Top of atmosphere albedo: ',alb_toai
         write(62,*)'  Top of model albedo: ',alb_tomi
         write(62,*)'  SW Flux absorbed at bottom boundary (Wm-2): '
     &                 ,SOLNET
         WRITE(62,*)'  Upwelling LW Flux at top-of-atmopshere: '
     &                 ,tiru
         Write(62,*)'  Up- & Downward SW flux at TOA: '
         WRITE(62,*)   tslu,tsld
         WRITE(62,2031)'1)P(Bars)',
     &        '2)TAULS (SW & LW)', 
     &        '3)CUMMULATIVE TAUS',
     &        '4)PI0s',
     &        '5)G0s',
     &        '6)DIRECT SW (d-scaled)',
     &        '7)4PI*|INTENSITIES|(W/M^2)'
 2031    FORMAT(3X,A9,6X,A17,9X,A18,12X,A6,11x,A5,7x,A20,3x,A26)            
         
          DO IL=1,NLAYER                                                  
          WRITE(62,2033),P_FULL(IL)*1e-5,uTAUL(1,IL),
     $    uTAUL(2,IL),uOPD(1,IL),uOPD(2,IL),uW0(1,IL),uW0(2,IL),uG0(1,IL),
     $    uG0(2,IL),DIRECT(1,IL),TMI(1,IL),TMI(2,IL)                   
 2033       FORMAT(F12.6,2X,E12.5,1X,E12.5,2X,E12.5,1x,E12.5,3X,
     $             F7.4,1x,F7.4,2X,F7.4,1X,F7.4,2X,E12.5,2X,
     $             E12.5,1X,E12.5) 
          END DO
          WRITE(62,*) ''
      ENDIF
      end

