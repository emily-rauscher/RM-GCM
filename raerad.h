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
      parameter( NSOL = 3)
      parameter( NIR  = 2)
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
      parameter( NWAVE = NSOL + NIR )
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
!  Declare common blocks for input to radiative transfer code
!
!   r_aerad        radius mid-pts from aerosol grids [cm]
!   rup_aerad      upper radii from aerosol grids [cm]
!   p_aerad        pressure [dyne/cm^2]
!   t_aerad        temperature [K]
!   pc_aerad       particle concentration [#/cm^3]
!   qv_aerad       water vapor mixing ratio [g/g]
!   tabove_aerad   blackbody temperature for downwelling IR flux into model [K]
!   ptop_aerad     pressure at top of aerosol model domain [dyne/cm^2]
!   pbot_aerad     pressure at bottom of aerosol model domain [dyne/cm^2]
!   u0_aerad       cosine of solar zenith angle [dimensionless]
!   emisir_aerad   surface IR emissivity [dimensionless]
!   tsfc_aerad     surface temperature [K]
!   h2ocol_aerad   water vapor column above model domain [g/cm^2]
!   isl_aerad      =1 means do solar calculations
!   ir_aerad       =1 means do infrared calculations
!   irs_aerad      =1 means do infrared scattering
!   ir_above_aerad =1 means include downwelling flux into top of model domain
!   psol_aerad     = incoming flux adjusted for latituded and daylength

      common / aerad1 /     
     &   isl_aerad, ir_aerad, irs_aerad, ir_above_aerad,
     &   ifsetup, ibinm,  
     &   riwp,rfluxes_aerad(2,2,2),
     &   r_aerad(NBIN,NGROUP),      
     &   rup_aerad(NBIN,NGROUP),     
     &   p_aerad(NZ_RAD), t_aerad(NZ_RAD),PRESS(NZ_RAD),      
     &   pc_aerad(NZ_RAD,NBIN,NGROUP),     
     &   qv_aerad(NZ_RAD), tabove_aerad,     
     &   ptop_aerad, pbot_aerad, u0_aerad,      
     &   emisir_aerad, tsfc_aerad, h2ocol_aerad,      
     &   iaerad1,psol_aerad
!
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
!  Declare common blocks for output from radiative transfer code
!
!   heati_aerad    infrared heating rates [K/s]
!   heats_aerad    solar heating rates [K/s]
!   qrad_aerad     particle radiative heating rates [K/s]
!   alb_tomi_aerad spectrally-integrated albedo at top-of-model
!   alb_toai_aerad spectrally-integrated albedo at top-of-atmosphere
!   alb_toa_aerad  spectrally-resolved albedo at top-of-atmosphere
!   opd_aerad      spectrally-resolved optical depth
!   opd_tot_aerad  spectrally-integrated optical depth
!   fsl_up_aerad   solar upwelling flux [W m^-2]
!   fsl_dn_aerad   solar downwelling flux [W m^-2]
!   fir_up_aerad   infrared upwelling flux [W m^-2]
!   fir_dn_aerad   infrared downwelling flux [W m^-2]
!   fsl_net_aerad   solar net flux [W m^-2]
!   fir_net_aerad   infrared net flux [W m^-2]!

      common / aerad2 /     
     &   wave_aerad(NWAVE+1),     
     &   heati_aerad(NZ_RAD), heats_aerad(NZ_RAD),     
     &   qrad_aerad(NZ_RAD,NBIN,NGROUP),     
     &   alb_tomi_aerad, alb_toai_aerad,     
     &   alb_toa_aerad(NSOL), opd_aerad(NWAVE), opd_tot_aerad,     
!     &   fsl_up_aerad(NZ_RAD+1), fsl_dn_aerad(NZ_RAD+1),     
!     &   fir_up_aerad(NZ_RAD+1), fir_dn_aerad(NZ_RAD+1),     
!     &   fir_net_aerad(NZ_RAD+1),fsl_net_aerad(NZ_RAD+1),
     &   fsl_up_aerad(NZ_RAD), fsl_dn_aerad(NZ_RAD),
     &   fir_up_aerad(NZ_RAD), fir_dn_aerad(NZ_RAD),
     &   fir_net_aerad(NZ_RAD),fsl_net_aerad(NZ_RAD),
     &   iaerad2
!
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
!  Declare common blocks for dharma
!
!   is_grp_ice_aerad =T means group is ice crystals
!   do_mie_aerad     =T means calculate Mie coefficients
!   ienconc_aerad    particle number conc. element for group
!   sfc_alb_aerad    surface albedo when fixed [dimensionless]
!   sfc_wind_aerad   surface wind [m/s]
!   dr_aerad         aerosol grid bin widths [cm]

      logical is_grp_ice_aerad, do_mie_aerad

      common / aerad3 /     
     &   is_grp_ice_aerad(NGROUP), do_mie_aerad,     
     &   ienconc_aerad(NGROUP),     
     &   sfc_alb_aerad, sfc_wind_aerad, dr_aerad(NBIN,NGROUP),     
     &   iaerad3

