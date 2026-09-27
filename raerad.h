!  @(#) aerad.h  McKie  Oct-1995
!  This is the include file for the interface between the aerosol
!  microphysics and radiative transfer components of CARMA.
!
!  Global symbolic constants are defined and common blocks are
!  declared.
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
!  Start of user-defined symbolic constants 
!
!
!  Define # grid pts in x, y, z directions
      parameter( NZ= NL+1)
!
!
!  Define # particle radius bins
!
      parameter( NBIN = 1)
!
!
!  Define # particle elements 
!
      parameter( NELEM = 1 )
!
!
!  Define # particle groups
!
      parameter( NGROUP = 1 )
!
!
!  Define # solutes
!
      parameter( NSOLUTE = 1 )
!
!
!  Define # gases
!
      parameter( NGAS = 1 )
!
!
!  Define # solar wavelength bins
!  Define # infrared wavelength bins
!      parameter( NSOL = 3)
!      parameter( NIR  = 2)
!
!  Define # layers in rad xfer model domain underlying aerosol model domain
!
      parameter( NZ_BELOW = 0 )
!
!
!  End of user-defined symbolic constants
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
!  The remaining symbolic constants will need no attention from most
!  users
!
!
!  Define # x, y,z direction grid box boundaries
!
      parameter( NZP1 = NZ + 1 )
!
!
!  Define total # wavelength bins
!
!      parameter( NWAVE = NSOL + NIR )
!
!
!  Define # layers in radiation model
!
      parameter( NZ_RAD = NZ + NZ_BELOW )
!
!
!  Define logical unit number for aerosol output print file
!
      parameter( LUNOPRT =  9 )
!
!
!  Define logical unit number for time step info output
!
      parameter( LUNOSTEP = 11 )
!
!
!  Define logical unit number for input restart file
!
      parameter( LUNIRES = 12 )
!
!
!  Define logical unit number for output restart file
!
      parameter( LUNORES = 13 )
!
!
!  Define logical unit number for output history file
!
      parameter( LUNOHIS = 14 )
!
!
!  Define logical unit number for input and output of Mie coefficients
!
      parameter( LUNMIE = 15 )
!
!
!  Define logical unit number for print output from radiation submodel
!
      parameter( LUNORAD = 16 )
!
!
!  Define logical unit number for time substep info output
!
      parameter( LUNOSUB = 17 )
!
!

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
!  Declare common blocks for dharma
!


