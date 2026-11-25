      module physical_constants

!...Define physical constants

       real, parameter :: PI = 3.14159
       real, parameter :: BK = 1.38054e-16
       real, parameter :: RGAS = 8.31430e+07
       real, parameter :: AVG = 6.02252e23
       real, parameter :: WTMOL_AIR = 28.966e+0
       real, parameter :: WTMOL_H2O = 18.
       real, parameter :: WTMOL_HNO3 = 63.
       real, parameter :: GASCON *10000. ! converting from J/kg/K to erg/K/g R_AIR = RGAS / WTMOL_AIR
       real, parameter :: R_AIR= GASCON
       real, parameter :: GRAV = GA *100. ! converting to cm/s^2 for Rad !980.6e+0
       real, parameter :: PREF = P0 !1013.e3
       real, parameter :: Cp = GASCON/AKAP !1.004e+7
       real, parameter :: RHOSOL = 1.38
       real, parameter :: RHOW = 1.
       real, parameter :: RHOI = 0.93
       real, parameter :: RHONAT = 1.62
       real, parameter :: PRENUC = 2.075e33 * RHOW / RHOI
       real, parameter :: BAL = 6.1121e3
       real, parameter :: BBL = 18.729
       real, parameter :: BCL = 257.87
       real, parameter :: BDL = 227.3
       real, parameter :: BAI = 6.1115e3
       real, parameter :: BBI = 23.036
       real, parameter :: BCI = 279.82
       real, parameter :: BDI = 333.7

       real, parameter :: RLHE_AVE = 2.5e10

       real, parameter :: Rv     = 461.e4
       real, parameter :: Rd     = 287.e4
       real, parameter :: L      = 2.5e10

       real, parameter :: eps  = Rd/Rv
       real, parameter :: RdCp = Rd/Cp

      end module physical_constants

