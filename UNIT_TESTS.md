[TOC]
# Unit Tests

## Spherical Differentiation

As of 7 March 2019, the results of the testing suite are below.
The measured convergence rate (~4) matches the intended.

```
Beginning tests for differentiation routines.
2-norm
Mean convergence rates: (lon, lat) = (-4.05, -3.99)
  Lon ( Points  Error )  :  Lat ( Points  Error )   (log2)
      ( 005    -0.9438 )  :      ( 005    -8.705 )
      ( 006    -4.081 )  :      ( 006    -12.64 )
      ( 007    -8.456 )  :      ( 007    -16.63 )
      ( 008    -12.73 )  :      ( 008    -20.62 )
      ( 009    -16.81 )  :      ( 009    -24.62 )
      ( 010    -20.83 )  :      ( 010    -28.62 )
      ( 011    -24.84 )  :      ( 011    -32.62 )
      ( 012    -28.85 )  :      ( 012    -36.6 )
inf-norm
Mean convergence rates: (lon, lat) = (-4.13, -3.97)
  Lon ( Points  Error )  :  Lat ( Points  Error )   (log2)
      ( 005    1.492 )  :      ( 005    -7.469 )
      ( 006    -0.8676 )  :      ( 006    -11.32 )
      ( 007    -5.067 )  :      ( 007    -15.29 )
      ( 008    -9.582 )  :      ( 008    -19.29 )
      ( 009    -13.99 )  :      ( 009    -23.28 )
      ( 010    -18.23 )  :      ( 010    -27.28 )
      ( 011    -22.37 )  :      ( 011    -31.28 )
      ( 012    -26.44 )  :      ( 012    -35.14 )
```

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
