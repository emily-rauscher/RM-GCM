# legacy/

Sources that are no longer part of the model build. They are kept for
reference and are not compiled by the `makefile`.

- **Pre-rradiation.f radiation chain:** `cmorc.f`, `cmorc_old.f`,
  `cnikos.f`, `cradsw.f`, `csw.f`, `cswtt.f`, `cirrad.f`, `cdedd.f`,
  `ccldtrn.f`, `ccolamt.f`, `co3interp.f`. It was replaced by
  `rradiation.f` and the `r*.f` Toon two-stream code.
- **Other unreachable routines:** `chanal1.f`, `cinisurf.f`, `cinterp.f`,
  `rinterpol.f`, `rplnk.f`, `rsurface.f`, `ranf.f`, `cbinary.f`,
  `cmynetcdf.f`, `double-gray-ropprmulti.f`, and the unused
  `physical_constants.f` module.
- **`fortran_readfort7.f90`:** an old standalone fort.7 reader. It is out
  of date; for example, it assumes 13 cloud species.
- **`compile_nopg`:** the old ifort build script, replaced by the
  `makefile`. It does not build the current model, because it lacks
  `cloudmods.f`, `tracerclds.f`, `dyeadv.f` and `calc_v_fall.f`.

Several of these still contain hand-written COMMON blocks (for example
`/CLOUDY/`) from before those blocks became modules. Bring them up to
date before moving any of them back into the build.
