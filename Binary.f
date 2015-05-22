!#################################
!# Calculates Flux for two stars #
!#################################

!      subroutine BinaryFlux(count,output)
      subroutine BinaryFlux(output,KOUNT,ITSPD)

      include 'params.i'

      real, parameter         :: STEFBOL = 5.67E-8        !StefanBoltzman
      real, parameter         :: RSUN      = 6.95E8         !Radius of Sun
      real, parameter         :: ASTROU        = 1.496E11

      real                    :: output

      real                    :: Ep                         !Ecc Anomoly Planet
      real                    :: Es                         !Ecc Anomoly Stars

!      integer                 :: count                          ! input current timestep


! Period of Planet = PORB (orbital period in units of rotation period) * 2pi/WW
! for Eccentric anomaly equation - need time in units of timesteps

! => P = Torb*1000 = units of timesteps
!These are Binary values from Fort.7
      COMMON/BINVAL/LBIN,PORBST,ECCPL,ECCST,SMAPL,SMAST,STMASS1,
     & STMASS2,STRAD1,STRAD2,STTEMP1,STTEMP2

      LOGICAL LBIN


! need this for PORB (planet orbit in planet days)
       COMMON/VARPARAM/OOM_IN, LPLOTMAP,NLPLOTMAP_IN,RFCOEFF_IN,
     & NTSTEP_IN, NSKIP_IN, BOTRELAXTIME, FBASEFLUX, FORCE1DDAYS,
     & OPACIR_POWERLAW, OPACIR_REFPRES, SOLC_IN, TOAALB,
     & PORB, OBLIQ, ECCEN

       LOGICAL LPLOTMAP

! meed this for pi and for common block with kount/tspd in it

      PARAMETER(MH=2,PI=3.14159265359,PI2=2.0*PI                          
     +,NNP=NN+1,MGPP=MG+2,JGP=JG+1,JGG=JG*NHEM,JGGP=JGG+1,MJP=NWJ2+NWJ2   
     +,NLM=NL-1,NLP=NL+1,NLPP=NL+2,NLA=NL+3,NLB=NL+4,NL2=NL*NL            
     +,IDA=(MG+MG+MG)/2+1,IDB=NWJ2*NL,IDC=IDB+IDB,IDD=MGPP*NL             
     +,IDE=NL2*NN,IDF=NCRAY*(MG+1),IDG=JG*NL,IDH=JG*MG                    
     +,IDI=NNP/2,IDJ=IDI*IDI,IDK=NL*IDI,IDL=MGPP/2,IDM=NNP/2,IDN=IDM*NL   
     +,NWW=1+(MM-1)/MOCT)                                                 



!below are parameters necessary to define so that calculations can be broken up into multiple steps ***CHANGE THIS

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

!      call ec_anom(ECCPL,PORB,Ep,count)
!      call ec_anom(ECCST,PORBST,Es,count)
      call ec_anom(ECCPL,PORB,Ep,KOUNT,ITSPD)
      call ec_anom(ECCST,PORBST,Es,KOUNT,ITSPD)
!SECOND - for the calcuated Ep and Es - corresponding angle on the orbit

      AngleP=2*atan(sqrt((1+ECCPL)/(1-ECCPL))*tan(Ep/2))
      AngleS=2*atan(sqrt((1+ECCST)/(1-ECCST))*tan(Es/2))

!THIRD - for calculated angle - Distance to COM for planet and star1 and star2
! R2 is conventionaly defined to be negative

      Rp=SMAPL*(1-ECCPL**2)/(1+ECCPL*cos(AngleP))
      Rs1=(STMASS2/(STMASS1+STMASS2))*SMAST*(1-ECCST**2)
     &/(1+ECCST*cos(AngleS))
      Rs2=-(STMASS1/(STMASS1+STMASS2))*SMAST*(1-ECCST**2)
     &/(1+ECCST*cos(AngleS))

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

      Ls1=4*PI*STEFBOL*(STRAD1*RSUN)**2*(STTEMP1**4)
      Ls2=4*PI*STEFBOL*(STRAD2*RSUN)**2*(STTEMP2**4)

      Fs1=Ls1/(4*PI*((Rp1*ASTROU)**2))
      Fs2=Ls2/(4*PI*((Rp2*ASTROU)**2))
      Ftot=Fs1+Fs2
      output=Ftot

      end subroutine

!#################################
! solves first kepler equation 
!################################

!      subroutine ec_anom(bine,binP,binx,count)                                                                                   
      subroutine ec_anom(bine,binP,binx,KOUNT,ITSPD)

      include 'params.i'

      real :: binx0,binx1,binx,bine,bina,binb,binP
!      integer :: count                                                                                                           
      PARAMETER(MH=2,PI=3.14159265359,PI2=2.0*PI                                                                           
     +,NNP=NN+1,MGPP=MG+2,JGP=JG+1,JGG=JG*NHEM,JGGP=JGG+1,MJP=NWJ2+NWJ2                                                    
     +,NLM=NL-1,NLP=NL+1,NLPP=NL+2,NLA=NL+3,NLB=NL+4,NL2=NL*NL                                                            
     +,IDA=(MG+MG+MG)/2+1,IDB=NWJ2*NL,IDC=IDB+IDB,IDD=MGPP*NL                                                              
     +,IDE=NL2*NN,IDF=NCRAY*(MG+1),IDG=JG*NL,IDH=JG*MG                                                                     
     +,IDI=NNP/2,IDJ=IDI*IDI,IDK=NL*IDI,IDL=MGPP/2,IDM=NNP/2,IDN=IDM*NL                                                   
     +,NWW=1+(MM-1)/MOCT)                                                                                                   

      binx0=0
      binx=binx0
      do
        bina=(binx0-bine*sin(binx0)-(PI2*KOUNT/(binP*ITSPD)))
        binb=1.0-bine*cos(binx0)
        binx1=binx0-(bina/binb)
        if (abs(binx0-bine*sin(binx0)-(PI2*KOUNT/(binP*ITSPD)))
     & <10E-10) exit
        binx0=binx1
      end do
      binx=binx1
      end subroutine


