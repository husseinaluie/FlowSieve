[TOC]
# Methods

This file will outline some of the computational methodologies.

## Coarse-Grain Filtering

### Region Selection for Filtering

## Across-scale Energy Transfer

## Baroclinic Energy Transfer


\f[
\Lambda^m = \frac{\overline{\mathbf{\omega}}\cdot \left( \nabla \overline{\rho} \times \nabla \overline{P} \right)}{\overline{\rho}}
\f]

In the case that \f$ \vec{\omega}=\left(\omega_r,0,0\right) \f$, then this reduces to
\f[
\Lambda^m = \frac{\overline{\omega_r}}{\overline{\rho}r^2\cos(\phi)}\left( \partial_{\lambda}\overline{\rho}\partial_{\phi}\overline{p} - \partial_{\phi}\overline{\rho}\partial_{\lambda}\overline{p} \right)
\f]

\f{eqnarray}{

\f}

## Differentiation

The spherical differentiation methods, [lat-derv], [lon-deriv], use a land-avoiding five-points stencil for fourth-order differentiation.
Land avoiding is achieved by simply using a non-centred stencil when appropriate. 
This is done to avoid having to use zero-velocity land cells, as this may introduce artificially steep velocity gradients that would confound derivatives.

### Cartesian derivatives
The secondary differentation tools ([x-deriv], [y-deriv], [z-deriv]) simply apply the chain rule on the spherical differentiation methods.

\f{eqnarray}{
\arraycolsep=2.5pt\def\arraystretch{2.5}
\left[ \begin{array}{c}
{\displaystyle \frac{\partial}{\partial x} } \\
{\displaystyle \frac{\partial}{\partial y} } \\
{\displaystyle \frac{\partial}{\partial z} }
\end{array} \right]
&=
\arraycolsep=2.5pt\def\arraystretch{2.5}
\left[ \begin{array}{ccc}
{\displaystyle \cos(\lambda)\cos(\phi) } 
  & {\displaystyle - \frac{\sin(\lambda)}{r\cos(\phi)} } 
  & {\displaystyle - \frac{\cos(\lambda)\sin(\phi)}{r} } \\
{\displaystyle \sin(\lambda)\cos(\phi) }
  & {\displaystyle   \frac{\cos(\lambda)}{r\cos(\phi)} }
  & {\displaystyle - \frac{\sin(\lambda)\sin(\phi)}{r} }\\
{\displaystyle \sin(\phi) }
  & {\displaystyle   0 }
  & {\displaystyle   \frac{\cos(\phi)}{r} }
\end{array} \right]
\arraycolsep=2.5pt\def\arraystretch{2.5}
\left[ \begin{array}{c}
{\displaystyle \frac{\partial}{\partial r} } \\
{\displaystyle \frac{\partial}{\partial \lambda} } \\
{\displaystyle \frac{\partial}{\partial \phi} }
\end{array} \right]
\f}

In the case that we're restricing ourselves to a spherical shell,
then \f$ \partial/\partial r\equiv0  \f$, so this reduced down to

\f{eqnarray}{
\arraycolsep=2.5pt\def\arraystretch{2.5}
\left[ \begin{array}{c}
{\displaystyle \frac{\partial}{\partial x} } \\
{\displaystyle \frac{\partial}{\partial y} } \\
{\displaystyle \frac{\partial}{\partial z} }
\end{array} \right]
=
\arraycolsep=2.5pt\def\arraystretch{2.5}
\left[ \begin{array}{ccc}
{\displaystyle - \frac{\sin(\lambda)}{r\cos(\phi)} } 
  & {\displaystyle - \frac{\cos(\lambda)\sin(\phi)}{r} } \\
{\displaystyle   \frac{\cos(\lambda)}{r\cos(\phi)} }
  & {\displaystyle - \frac{\sin(\lambda)\sin(\phi)}{r} }\\
{\displaystyle   0 }
  & {\displaystyle   \frac{\cos(\phi)}{r} }
\end{array} \right]
\arraycolsep=2.5pt\def\arraystretch{2.5}
\left[ \begin{array}{c}
{\displaystyle \frac{\partial}{\partial \lambda} } \\
{\displaystyle \frac{\partial}{\partial \phi} }
\end{array} \right]
\f}

### Testing

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

[lat-deriv]: @ref latitude_derivative_at_point "latitude derivative"
[lon-deriv]: @ref longitude_derivative_at_point "longitude derivative"
[x-deriv]: @ref x_derivative_at_point "x derivative"
[y-deriv]: @ref y_derivative_at_point "y derivative"
[z-deriv]: @ref z_derivative_at_point "z derivative"
