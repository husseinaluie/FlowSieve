[TOC]
# Methods

This file will outline some of the computational methodologies.

## Coarse-Grain Filtering

### Region Selection for Filtering

## Across-scale Energy Transfer

## Baroclinic Energy Transfer


\f[
\Lambda^m = \frac{\overline{\omega}\cdot \left( \nabla \overline{\rho} \times \nabla \overline{P} \right)}{\overline{\rho}}
\f]

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
For reasons unknown, the order of convergence (8) is twice the intended (4).

```
Beginning tests for differentiation routines.
Mean convergence rates: (lon, lat) = (-8.08, -7.98)
  Lon ( Points  Error )  :  Lat ( Points  Error )   (log2)
      ( 005    -1.876 )  :      ( 005    -17.29 )
      ( 006    -8.264 )  :      ( 006    -25.17 )
      ( 007    -16.95 )  :      ( 007    -33.14 )
      ( 008    -25.43 )  :      ( 008    -41.13 )
      ( 009    -33.56 )  :      ( 009    -49.13 )
      ( 010    -41.59 )  :      ( 010    -57.13 )
      ( 011    -49.6 )  :      ( 011    -65.13 )
      ( 012    -57.6 )  :      ( 012    -73.08 )
```

[lat-deriv]: @ref latitude_derivative_at_point "latitude derivative"
[lon-deriv]: @ref longitude_derivative_at_point "longitude derivative"
[x-deriv]: @ref x_derivative_at_point "x derivative"
[y-deriv]: @ref y_derivative_at_point "y derivative"
[z-deriv]: @ref z_derivative_at_point "z derivative"
