# Unit Tests {#unittests1}

# Unit Testing
[TOC]

---

**NOTE**

Unit-testing is currently out-dated and broadly non-functional.
The results of the previous unit testings are includied below,
which tested that the differentiation schemes achieved the stated order of accuracy.

---

## Spherical Differentiation

As of 12 March 2019, the results of the testing suite are below.
The measured convergence rate  matches the intended.
Note that in the case of sixth order differentiation the reported
convergence rate appears lower because of near convergence at higher
resolutions.

### Second order

```
Beginning tests for differentiation routines.
At highest resolution, 6.55% of tiles were land.
2-norm
Mean convergence rates: (lon, lat) = (-2, -1.99)
  Lon ( Points  Error  Ord )  :  Lat ( Points  Error  Ord )   ( log2 of pts and err )
      ( 005    0.09319  -1.91 )  :      ( 005    -1.449  -1.95 )
      ( 006    -1.82  -2.01 )  :      ( 006    -3.395  -1.99 )
      ( 007    -3.826  -2.01 )  :      ( 007    -5.387  -2 )
      ( 008    -5.834  -2.01 )  :      ( 008    -7.383  -2 )
      ( 009    -7.843  -2 )  :      ( 009    -9.383  -2 )
      ( 010    -9.847  -2 )  :      ( 010    -11.38  -2 )
      ( 011    -11.85  -2 )  :      ( 011    -13.38  -2 )
      ( 012    -13.85  -2 )  :      ( 012    -15.38  -2 )
inf-norm
Mean convergence rates: (lon, lat) = (-2.01, -1.98)
  Lon ( Points  Error  Ord )  :  Lat ( Points  Error  Ord )   ( log2 of pts and err )
      ( 005    2.456  -1.83 )  :      ( 005    -0.06331  -1.88 )
      ( 006    0.621  -2.04 )  :      ( 006    -1.947  -1.97 )
      ( 007    -1.417  -2.05 )  :      ( 007    -3.918  -1.99 )
      ( 008    -3.466  -2.03 )  :      ( 008    -5.911  -2 )
      ( 009    -5.496  -2.02 )  :      ( 009    -7.909  -2 )
      ( 010    -7.514  -2.01 )  :      ( 010    -9.909  -2 )
      ( 011    -9.524  -2 )  :      ( 011    -11.91  -2 )
      ( 012    -11.53  -2 )  :      ( 012    -13.91  -2 )
```

### Fourth order

```
Beginning tests for differentiation routines.
At highest resolution, 6.55% of tiles were land.
2-norm
Mean convergence rates: (lon, lat) = (-4.07, -3.98)
  Lon ( Points  Error  Ord )  :  Lat ( Points  Error  Ord )   ( log2 of pts and err )
      ( 005    -0.4418  -3.92 )  :      ( 005    -4.058  -3.86 )
      ( 006    -4.366  -4.16 )  :      ( 006    -7.92  -3.97 )
      ( 007    -8.522  -4.14 )  :      ( 007    -11.89  -3.99 )
      ( 008    -12.66  -4.09 )  :      ( 008    -15.88  -4 )
      ( 009    -16.76  -4.05 )  :      ( 009    -19.88  -4 )
      ( 010    -20.81  -4.03 )  :      ( 010    -23.88  -4 )
      ( 011    -24.84  -4.02 )  :      ( 011    -27.88  -4 )
      ( 012    -28.86  -4.02 )  :      ( 012    -31.87  -4 )
inf-norm
Mean convergence rates: (lon, lat) = (-3.97, -3.88)
  Lon ( Points  Error  Ord )  :  Lat ( Points  Error  Ord )   ( log2 of pts and err )
      ( 005    2.503  -3.37 )  :      ( 005    -2.677  -3.81 )
      ( 006    -0.862  -3.94 )  :      ( 006    -6.487  -3.86 )
      ( 007    -4.806  -4.05 )  :      ( 007    -10.35  -3.71 )
      ( 008    -8.86  -4.05 )  :      ( 008    -14.06  -3.85 )
      ( 009    -12.91  -4.03 )  :      ( 009    -17.91  -4.15 )
      ( 010    -16.93  -4.02 )  :      ( 010    -22.05  -3.8 )
      ( 011    -20.95  -4.01 )  :      ( 011    -25.86  -4 )
      ( 012    -24.96  -4.01 )  :      ( 012    -29.86  -4 )
```

### Sixth order

```
Beginning tests for differentiation routines.
At highest resolution, 6.55% of tiles were land.
2-norm
Mean convergence rates: (lon, lat) = (-6, -5.22)
  Lon ( Points  Error  Ord )  :  Lat ( Points  Error  Ord )   ( log2 of pts and err )
      ( 005    -0.1769  -5.8 )  :      ( 005    -6.27  -5.77 )
      ( 006    -5.982  -6.29 )  :      ( 006    -12.04  -5.93 )
      ( 007    -12.28  -6.38 )  :      ( 007    -17.97  -5.99 )
      ( 008    -18.66  -6.36 )  :      ( 008    -23.96  -6 )
      ( 009    -25.01  -6.29 )  :      ( 009    -29.96  -5.99 )
      ( 010    -31.3  -6.23 )  :      ( 010    -35.96  -5.01 )
      ( 011    -37.53  -3.08 )  :      ( 011    -40.96  1.23 )
      ( 012    -40.61  -3.08 )  :      ( 012    -39.73  1.23 )
inf-norm
Mean convergence rates: (lon, lat) = (-5.75, -5.02)
  Lon ( Points  Error  Ord )  :  Lat ( Points  Error  Ord )   ( log2 of pts and err )
      ( 005    3.163  -5.37 )  :      ( 005    -4.147  -6.48 )
      ( 006    -2.203  -5.76 )  :      ( 006    -10.63  -4.38 )
      ( 007    -7.968  -6.04 )  :      ( 007    -15.01  -6.11 )
      ( 008    -14.01  -6.06 )  :      ( 008    -21.12  -5.8 )
      ( 009    -20.06  -6.04 )  :      ( 009    -26.92  -6.03 )
      ( 010    -26.1  -5.99 )  :      ( 010    -32.95  -4.85 )
      ( 011    -32.09  -3.76 )  :      ( 011    -37.81  1.3 )
      ( 012    -35.85  -3.76 )  :      ( 012    -36.5  1.3 )
```

---

## Vorticity

As of 8 March 2019, the results for the vorticity tests are below.
The accuracy scales with resolution of course (since it relies on differentiation),
but with reasonable resolution to avoid errors in the derivatives themselves, this shows
that the equations for vorticity are correct.

Note that the norms are normalized with the triangle equality to give a [0,1] metric.

```
Beginning vorticity tests.
  (normalized) 2-norm   error: 4.467e-07
  (normalized) inf-norm error: 1.473e-06
```
