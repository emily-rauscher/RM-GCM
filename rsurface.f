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


!     P=(rho/m)k*T  ==> using pressure and temperature of mid level with mass of entire level
      dens_boa=(pboa/tboa)/(GASCON*10E4)
      
!     Flux emitted from surface LW. Store to common block here

      FLWE=SURF_EMIS*stbltz*(TGRND0*TGRND0*TGRND0*TGRND0)
      SHEAT=dens_boa*(CPD*10E4)*trans_c*BOAWIND*(tboa-TGRND0)


      IF (ISNAN(FLWD)) THEN
         FLWD=0.
      ENDIF

      Flux_sum=(FBASEFLUX*1E3)+FSWD+FLWD-FLWE-SHEAT
      SURFES=Flux_sum   
!     above saves Flux Sum as W/m2 in Surface to be writen to array later
!     used to calculate the energy stored in surface

      deltaTs=DELTAT*(Flux_sum)/(1000.*CSURF*RHOSURF*GRNDZ)
      TGRND0=TGRND0+deltaTs
 
      END SUBROUTINE
