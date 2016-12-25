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
!      integer, dimension(20) :: ibinmin

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

!...Calculate z and dz based on hydrostatic equilibrium

!pmean is the average pressure for the layer.i.e. for each layer,
! the base and top are added and divided by 2. It's in pascals. 

!      z(1) = 0.
!      do iz = 2, NZ
!        scale = R_AIR*0.5*(t_full(iz)+t_full(iz-1)) / GRAV
!        pmean = 0.5*(prflux(iz)+prflux(iz-1))
!        z(iz) = z(iz-1) - scale * (prflux(iz)-prflux(iz-1)) / pmean
!        write(*,*) 'scale',scale
!        write(*,*) 'z(iz)',scale * (prflux(iz)-prflux(iz-1)) / pmean
!        write(*,*)'pmean',pmean
!  print*, iz, pmean/1.e3, t_full(iz), scale/1.e5, z(iz)/1.e5
!      enddo
!        dz(1)=z(2)*(-1)
!      do iz = 2,NZ-1
!        dz(iz) = -1*(z(iz+1)-z(iz))
!        print*, iz, z(iz)/1.e5, dz(iz)/1.e5
!      enddo
!        write(*,*)'dz',dz

!      do iz = 1, NZ
!print*, 'p_full = ', p_full(iz)/1.e3
!      enddo
!print*, 't_full = ', t_full
!print*, 'qh2o_full = ', qh2o_full

! open( 9, file='tropical_profile.dat', form='formatted',
! status='unknown' )
! do i = 1, 208
!   read(9,*) ii, p_full(i), t_full(i), qh2o_full(i)
! enddo


!...Call radiative transfer code

!      do i = 1,2 !20
!        ibinmin(i) = i
!      enddo

!...for case with small crystal artifacts removed


!do ii = 1, 20
!      do ii = 1, 1

!        ibinm = ibinmin(ii)
!         write(*,*) 'calling radsub'
!         write(*,*) 'iffirst',iffirst
!         write(*,*) 'dz', dz
!         write(*,*) 'p_full', p_full
!         write(*,*) 't_full', t_full
!         write(*,*) 'qh2ofull',qh2o_full
!         write(*,*) 'radheat',radheat
!         write(*,*) 'wave_pass',wave_pass
!         write(*,*) 'ibinm',ibinm        
         call radsub(iffirst, pr,p_full,t_full,qh2o_full,
     &               radheat,htlw,htsw,rfluxes,alat1,alon,KOUNT,ITSPD)

        iffirst = 0
      
!         write(*,*) 'radheat',radheat
!         write(*,*) 'htlw',htlw
!         write(*,*) 'htsw',htsw
!         write(*,*) 'fluxes calcrad',rfluxes 
!      enddo
!        write(*,*) 'fsl_up_aerad',fsl_up_aerad
!        write(*,*) 'fsl_dn_aerad',fsl_dn_aerad
!        write(*,*) 'fsl_net_aerad',fsl_net_aerad
!        write(*,*) 'fir_up_aerad',fir_up_aerad
!        write(*,*) 'fir_dn_aerad',fir_dn_aerad
!        write(*,*) 'fir_net_aerad',fir_net_aerad
!        write(*,*) 'htlw',htlw
!        write(*,*) 'htsw',htsw
!        write(*,*) 'pr',pr
!!      LFLUXDIAG=.TRUE. !! Switch on or off diagnostics in fort.63
!below
!          write(*,*)'PSOL',alat1,alon
!          write(*,*) 'FNET',FNET(1,2)
!         DO ILAY=NLAYER,1,-1                                                  
!
!          WRITE(*,2014),P_FULL(NLAYER-ILAY+1)*1e-5,FIR_UP_AERAD(ILAY),                            
!     $    FIR_DN_AERAD(ILAY),FIR_NET_AERAD(ILAY),
!     $    FSL_UP_AERAD(ILAY),FSL_DN_AERAD(ILAY),FSL_NET_AERAD(ILAY)                         
!     $
! 2014       FORMAT(2X,F12.6,3X,E12.5,3X,E12.5,3X,E12.5,3X,E12.5
!     $             ,3X,E12.5,3X,E12.5)
!         END DO  
!         stop         
C ER Modif: only write fort.63 every kountp timesteps
C     (kountp-1 b/c in cmorc nikos called when mod(kount,ntstep).eq.1)
      ! write(*,*) 'SOLNET',SOLNET
      ! write(*,*) 'alb_toai_aerad',alb_toai_aerad
      ! write(*,*)' alb_tomi_aerad',alb_tomi_aerad
      ! write(*,*) 'TMI',TMI
      IF ((LFLUXDIAG).AND.(KOUTP.EQ.KOUNTP-1)) THEN
         WRITE(63,*) 'LATITUDE, LONGITUDE:',ALAT1,ALON
         WRITE(63,2010)'1)P(Bars)',
     $        'FLUXES(Wm-2): 2)LW UP','3)LW DOWN',         
     $        '4)LW NET',
     $        '5)SW UP','6)SW DOWN',
     $        '7)SW NET'
 2010    FORMAT(1X,A16,1X,A21,4X,A9,4X,A8,4x,A7,4x,A9,4x,A8)                               

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
      ENDIF  
      end

