      subroutine ec_anom(e,P,x,i)
      real :: x0,x1,x,e,a,b
      real, parameter :: PI = 3.1415

      integer :: i, P

      x0=0
      x=x0
      do
        a=(x0-e*sin(x0)-(2.0*PI*i/P))
        b=1.0-e*cos(x0)
        x1=x0-(a/b)
        if (abs(x0-e*sin(x0)-(2.0*PI*i/P)) <10E-10) exit
        x0=x1
      end do
      x=x1
      end subroutine


      subroutine BinaryFlux(i,output)
      IMPLICIT NONE

      real, parameter         :: SIGMA = 5.67E-8        !StefanBoltzman
      real, parameter         :: RSUN      = 6.95E8         !Radius of Sun
      real, parameter         :: PI        = 3.1415
      real, parameter         :: AU        = 1.496E11

      real                    :: output

      real                    :: Ep                         !Ecc Anomoly Planet
      real                    :: Es                         !Ecc Anomoly Stars

      integer                 :: i                          ! input current timestep


! Period of Planet = PORB (orbital period in units of rotation period) * 2pi/WW
! for Eccentric anomaly equation - need time in units of timesteps

! => P = Torb*1000 = units of timesteps
!These are Binary values from Fort.7
      COMMON/BINVAL/LBIN,PORBST,ECCPL,ECCST,SMAPL,SMAST,STMASS1,
     & STMASS2,STRAD1,STRAD2,STTTEMP1,STTEMP2

      LOGICAL LBIN


! need this for PORB (planet orbit in planet days)
       COMMON/VARPARAM/OOM_IN, LPLOTMAP,NLPLOTMAP_IN,RFCOEFF_IN,
     & NTSTEP_IN, NSKIP_IN, BOTRELAXTIME, FBASEFLUX, FORCE1DDAYS,
     & OPACIR_POWERLAW, OPACIR_REFPRES, SOLC_IN, TOAALB,
     & PORB, OBLIQ, ECCEN

      integer, parameter         :: Pplanet        = Torbplanet*1000  !*TSPD
      integer, parameter         :: Pstars         = Torbstars*1000   !*TSPD


!below are parameters necessary to define so that calculations can be broken up into multiple steps for ease

      real :: AngleP                   !Theta, Planet
      real :: AngleS                   !Theta, Star
      real :: Rp                       !Planet Orbital Distance
      real :: Rs1                      !Orbital Distance - S1
      real :: Rs2                      !Orbital Distance - S2

      real :: Rpx                      !Orbital Distance - P,x
      real :: Rpy                      !Orbital Distance - P,y
      real :: Rs1x                     !Orbital Distance - S1,x
      real :: Rs1y                     !Orbital Distance - S1,y
      real :: Rs2x                     !Orbital Distance - S2,x
      real :: Rs2y                     !Orbital Distance - S2,y
      real :: RP1x                     !Orbital Distance - P to S1,x
      real :: RP1y                     !Orbital Distance - P to S1,y
      real :: RP2x                     !Orbital Distance - P to S2,x
      real :: RP2y                     !Orbital Distance - P to S2,y

      real :: Rp1                      !Distance planet to star1
      real :: Rp2                      !Distance planet to star2
      real :: Fs1                      !Flux recieved from S1
      real :: Fs2                      !Flux recieved from S2
      real :: FTot                     !Total Flux

      real :: Ls1
      real :: Ls2



!at each timestep - need to calcuate Eccentric Annomaly of Planet and Stars
! and then Theta of Planet and Stars
! and then distance from center of mass

! SEE KEPLERS EQUATIONS
!(2pi/P)t=E-esinE
!tan(theta/2)=sqrt((1+e)/(1-e))tan(E/2)
!R=z*a(1-e^2)/(1+ecos(theta)
!where z=1 planet, z=m2/(m1+m2) for star1 and z= m1/(m1+m2) for star2

!Flux is dependent on Distance from STAR not center of mass -> need to caclucate this

!FIRST - given timestep, calculate Ec_p (Ep) and Ec_s (Ec)
! i is the given timestep to be calcuated for

      call ec_anom(eccp,Pplanet,Ep,i)
      call ec_anom(eccs,Pstars,Es,i)

!SECOND - for the calcuated Ep and Es - corresponding angle on the orbit

      AngleP=2*atan(sqrt((1+eccp)/(1-eccp))*tan(Ep/2))
      AngleS=2*atan(sqrt((1+eccs)/(1-eccs))*tan(Es/2))

!THIRD - for calculated angle - Distance to COM for planet and star1 and star2
! R2 is conventionaly defined to be negative

      Rp=ap*(1-eccp**2)/(1+eccp*cos(AngleP))
      Rs1=(mstar2/(mstar1+mstar2))*as*(1-eccs**2)/(1+eccs*cos(AngleS))
      Rs2=-(mstar1/(mstar1+mstar2))*as*(1-eccs**2)/(1+eccs*cos(AngleS))

! FOURTH - distance from planet to each star

      Rpx = Rp*cos(AngleP)
      Rpy = Rp*sin(AngleP)

      Rs1x = Rs1*cos(AngleS)
      Rs1y = Rs1*sin(AngleS)

      Rs2x = Rs2*cos(AngleS)
      Rs2y = Rs2*sin(AngleS)

      RP1x = Rpx - Rs1x
      RP1y = Rpy - Rs1y

      RP2x = Rpx - Rs2x
      RP2y = Rpy - Rs2y

      Rp1=sqrt((RP1x**2)+(RP1y**2))
      Rp2=sqrt((RP2x**2)+(RP2y**2))


! LAST - flux (in units of W/m^2)

      Ls1=4*PI*SIGMA*(Rstar1)**2*(Tstar1**4)
      Ls2=4*PI*SIGMA*(Rstar2)**2*(Tstar2**4)

      Fs1=Ls1/(4*PI*(Rp1**2))
      Fs2=Ls2/(4*PI*(Rp2**2))
      Ftot=Fs1+Fs2
      output=Ftot

      end subroutine


