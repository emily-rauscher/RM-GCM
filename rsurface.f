      SUBROUTINE SURFACETEMP(tboa,pboa)
!
!     ***************************************************************************************
!     Takes boundary fluxes and calculates a surface temperature
!     based on the 1D heat equation, assumes no horizontal heat transport
!
!     done in rradtran, using current fluxes. Is recalculated before 
!     fluxes are redone in rnewflux1. 
!
!     relevant fluxes are added to arrays in rradtran at end (reflected and re-emitted fluxes)
!     ****************************************************************************************
!
!
!     NEED to include, Temp array (in rradtran, called TT array, NL+1 dimension
!                      Winds array (in rradiation, UG (IGC,NL) dimensions) -> long,level
!                      Fluxes at BOA (fluxes are passed back to rrflux array which stores all fluxes at all lat long)
!                              (long,lat,6) figure out how to pass current lat/long to this array. 
!                              1--SW,DN,SURF **
!                              2--SW,UP,SURF
!                              3--LW,DN,SURF **
!                              4--LW,UP,SURF 
!                              5--SW,DN-UP,TOA 
!                              6--LW,UP,SURF 
!                      Pressures for calculating density

!      include 'params.i'
      include 'rcommons.h'

!     JH holds current lat. Need to pull prev timestep fluxes and TGRND
      COMMON/LEGAU/ ALPJ(MJP),DALPJ(MJP)                                                                                                                  
     +              ,ALP(NWJ2,2,JGL),DALP(NWJ2,2,JGL)                                                                                                         
     +              ,RLP(NWJ2,2,JGL),RDLP(NWJ2,2,JGL)                                                                                                          
     +              ,SI(JGG),CS(JGG),SISQ(JGG),CSSQ(JGG),SECSQ(JGG)                                                                                            
     +              ,ALAT(JGG),GWT(JGG),AW(JGG),JH,JL,JINC  

!     Boltzman in erg/K, stefboltz in cgs
      PARAMETER(bltz=1.380658E-16,stbltz=5.670367E-5)

      REAL :: SHEAT            !hold relevant fluxes
      REAL :: mass_molar, dens_boa              !mass and density of this block
      REAL :: trans_c                         ! transfer coefficient

      trans_c=10E-3
      mass_molar=28.97

!     size of timestep in SECONDS
!      tday_sec=PI2/WW

!     Mass and Density of bottom layer. Using T(NL) because this is the mid level temperature
!     DPG is calculated as PRESS(J)-PRESS(J-1)/G in rsetup_simple.
!      rcm=RADEA*100.    !to cgs
!      csurf=csurf*1000. !to cgs

!      SEC=ALAT(JH)/57.29578 
!      mass_boa=DPG(NL)*rcm*rcm*(2.*PI/MG)*(PI/2./JG)*cos(SEC)

!     P=(rho/m)k*T  ==> using pressure and temperature of mid level with mass of entire level
      dens_boa=(pboa/tboa)/(GASCON*10E4)
      
!     Flux emitted from surface LW. Store to common block here

      FLWE=SURF_EMIS*stbltz*(TGRND0*TGRND0*TGRND0*TGRND0)
      SHEAT=dens_boa*(CPD*10E4)*trans_c*BOAWIND*(tboa-TGRND0)

!      write(*,*) '-------------> RSURFACE'
!      write(*,*) FLWE,SHEAT,FLWD,FWSD
!      write(*,*) 'MASS=', mass_boa
!      write(*,*) 'Density=', dens_boa
!      write(*,*) 'Pressure=', pboa
!      write(*,*) 'FLWE=', FLWE
!      write(*,*) 'SHEAT=', SHEAT
!      write(*,*) 'FLWD=', FLWD
!      write(*,*) 'FSWD=', FSWD
!      write(*,*) 'Fbase=', FBASEFLUX*1E3
!      write(*,*) '-- POSITIVE FLUXES=', FLWD+FSWD+(FBASEFLUX)*(1E3)
!      write(*,*) '-- NEGATIVE FLUXES=', FLWE+SHEAT
!      write(*,*) tboa, TGRND0, trans_c, gascon

!      IF (ISNAN(FLWE)) THEN
!         FLWE=0.
!      ENDIF
!      IF (ISNAN(SHEAT)) THEN
!         SHEAT=0.
!      ENDIF
      IF (ISNAN(FLWD)) THEN
         FLWD=0.
      ENDIF
!      IF (FSWD.eq.0) THEN 
!         FSWD=0.
!         FLWD=0.
!         FLWE=0.
!         SHEAT=0.
!      ENDIF

!      write(*,*) FLWE,SHEAT,FLWD,FWSD
      
!     FBASEFLUX*10E3
      Flux_sum=(FBASEFLUX*1E3)+FSWD+FLWD-FLWE-SHEAT
!      write(*,*) 'FLUX_SUM=', Flux_sum

      deltaTs=DELTAT*(Flux_sum)/(1000.*CSURF*RHOSURF*GRNDZ)
!      write(*,*) 'CONSTATS=', CSURF*RHOSURF*GRNDZ
!      write(*,*) 'DELTA_Ts=', deltaTs

      TGRND0=TGRND0+deltaTs
!      write(*,*) 'Updated TGRND0=', TGRND0
!      write(*,*) Flux_sum, TT(NL+1),deltaTs,TGRND0
 
      END SUBROUTINE
