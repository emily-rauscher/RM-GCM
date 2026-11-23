!     Number of cloud species, indexed as MOLEF is (1 = KCl ... 7 =
!     Mg2SiO4 ... 13 = Al2O3, 14 = MgSiO3, which has corr-k tables but
!     no CLOUD_DATA optics, so is for corr-k only). Every per-species array is sized by it;
!     adding a species means raising it, adding the species' name
!     (corrk.f CLOUD_SPECIES), density and molar mass
!     (cloud_properties_set_up.f) and its optics block there, and
!     supplying its data files (a condensation-curve column, CLOUD_DATA
!     optics tables, corr-k tables in each band folder).
!     No type statement: like the PARAMETERs in params.i, NCLOUDS is an
!     integer by its initial letter, and a type statement here would
!     come before the IMPLICIT in rprecision.h, which Fortran forbids.
!     Scopes with IMPLICIT NONE declare INTEGER NCLOUDS themselves.
      PARAMETER(NCLOUDS=14)
