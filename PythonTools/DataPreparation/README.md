

* salinity_interp.py
  * input: salinity.nc
  * output (appends): source.nc
  * purpose: interpolates salinity fields onto the
             same time points as the other fields
  * reason: the source fields used have salinity at a 
            courser time grid than the other fields

* subsetter.py
  * input: full_domain.nc
  * output (overwritten): subset.nc
  * purpose: given a large domain, select a contiguous
             subset and save into a separate file

* relabel.py
  * input: interp.nc
  * output: input.nc
  * purpose: given an output from interpolator.x,
             relabel fields appropriately for the 
             coarse_grain routine

How To:
1. Interpolate salinity using salinity_interp.py. Then rename source.nc as full_domain.nc.
2. If desired, call subsetter.py to produce subset.nc
3. Rename subset.nc as source.nc
4. Run the interpolator (interpolator.x)
5. Relabel the fields (relabel.py)
6. Move input.nc to the appropriate directory and run coarse_grain.x
