[TOC]
# Methods

This file will outline some of the computational methodologies.

## Coarse-Grain Filtering

### Region Selection for Filtering

Suppose that we are computing the filter at longitude \f$\lambda_0\f$ and latitude \f$\phi_0\f$ (at indices \f$I_{\lambda_0}\f$ and \f$I_{\phi_0}\f$ respectively) over length scale \f$L\f$. 
Further, suppose that we have uniform (spherical) grid spacing \f$\Delta\lambda\f$ and \f$\Delta\phi\f$. 
Let \f$R_E\f$ be the mean radius of the earth. 

The spacing between latitude bands, in *metres*, is then \f$\Delta\phi_m=\Delta\phi R_E\f$.
The number of points that we'll include in either latitude side is then \f$\Delta\phi_N=\mathrm{ceil}\left(1.1 (L/2) / \Delta\phi_m\right)\f$
Note that the 1.1 factor is to take a region slightly larger than the kernel (assuming a compact support kernel).

The outer loop in the integral is then from \f$I_{\phi_0} - \Delta\phi_N\f$ to \f$I_{\phi_0}+\Delta\phi_N\f$.
Next, at each latitude \f$\phi\f$ within that loop, the same process is applied to compute the bounds for the longitude loop, with
the spacing between longitude band, in *metres*, given by \f$\Delta\lambda_m=\Delta\lambda R_E \cos(\phi)\f$.

## Evolution of Kinetic Energy of Filter Velocities

\f[
\frac{\partial}{\partial t}\left( \frac{\rho_0}{2}\left| \overline{\vec{u}} \right|^2 \right)
+ \nabla\cdot J_{\mathrm{transport}} 
=
-\Pi 
- \rho_0\nu\left| \nabla\overline{\vec{u}}\right|^2
+ \overline{\rho}\vec{g}\cdot\overline{\vec{u}}
+\rho_0\overline{F}_{\mathrm{forcing}}\cdot\overline{\vec{u}}
\f]

### Across-scale Kinetic Energy Transfer

\f[ \Pi = -\rho_0 S_{ij}\tau_{ij} = -\frac{\rho_0}{2}\left(\overline{u}_{i,j}+\overline{u}_{j,i}\right)\left(\overline{u_iu_j} - \overline{u_i}\cdot\overline{u_j}\right) \f]

Unit-wise,
\f[ 
\left[\rho_0 \right]=\frac{\mathrm{kg}}{\mathrm{m}^3}, \qquad
\left[S_{ij} \right]= \frac{1}{\mathrm{s}} , \qquad
\left[\tau_{ij} \right]= \frac{\mathrm{m}^2}{\mathrm{s}^2} , \qquad
\left[ W \right] = \frac{\mathrm{kg}\cdot\mathrm{m}^2}{\mathrm{s}^3}, \qquad
\Rightarrow \qquad
\left[ \Pi\right] = \frac{\mathrm{W}}{\mathrm{m}^3}
\f]

### Advection of Kinetic Energy

\f{eqnarray*}{
J_{\mathrm{transport}} =\qquad\qquad\qquad
\frac{\rho_0}{2} \left| \overline{\vec{u}} \right|^2 \overline{\vec{u}} &\mapsto& \mathrm{Advection}\,\mathrm{by}\,\mathrm{coarse}\mathrm{-}\mathrm{scale}\,\mathrm{flow} \\
+ \overline{P} \overline{u} &\mapsto& \mathrm{Pressure}\mathrm{-}\mathrm{induced}\,\mathrm{transport} \\
- \nu \frac{\rho}{2}  \nabla\left( \left| \overline{\vec{u}} \right|^2 \right) &\mapsto& \mathrm{Diffusion} \\
+ \rho_0  \overline{\vec{u}}\,\,\overline{\tau}(\vec{u}, \vec{u}) &\mapsto& \mathrm{Advection}\,\mathrm{by}\,\mathrm{fine}\mathrm{-}\mathrm{scale}\,\mathrm{flow}
\f}

The following assume incompressibility: \f[\nabla\cdot\overline{\vec{u}}=\overline{u}_{j,j}=0\f]

#### Advection by coarse-scale Flow
\f{eqnarray}{
\nabla\cdot\left(\frac{\rho_0}{2}  \left| \overline{\vec{u}} \right|^2  \overline{\vec{u}}\right) &=&  \frac{\rho_0}{2}  [ (\overline{u}_i \overline{u}_i)  \overline{u}_j ]_{,j} \\
&=&  \frac{\rho_0}{2} (\overline{u}_i \overline{u}_i)_{,j}  \overline{u}_j \\
&=& \rho_0  \overline{u}_i  \overline{u}_{i,j}  \overline{u}_j
\f}

#### Pressure-induced transport

\f{eqnarray}{
\nabla\cdot\left(\overline{P}\,\overline{u}\right) &=& (\overline{p}\,\overline{u}_j)_{,j} \\
&=& \overline{p}_{,j}\, \overline{u}_j
\f}

#### Diffusion

Not implemented.

#### Advection by fine-scale flow

\f{eqnarray}{
\nabla\cdot\left(  \rho_0  \overline{\vec{u}} \,\,\overline{\tau}(\vec{u}, \vec{u})  \right)
&=& \rho_0  [ \overline{u}_i \tau_{ij} ]_{,j} \\
&=& \rho_0 \left(
\overline{u}_{i,j} \left[ \overline{(u_iu_j)}  - \overline{u}_i\,\overline{u}_j \right]
+ \overline{u}_i  \left[ \overline{(u_iu_j)}_{,j} - \overline{u}_{i,j}\,\overline{u}_j \right]
\right)
\f}


## Baroclinic Energy Transfer


\f[
\Lambda^m = \frac{\overline{\mathbf{\omega}}\cdot \left( \nabla \overline{\rho} \times \nabla \overline{p} \right)}{\overline{\rho}}
\f]

Unit-wise, this gives 
\f[ 
\left[\omega \right]=\frac{1}{\mathrm{s}}, \qquad
\left[\nabla \right]= \frac{1}{\mathrm{m}} , \qquad
\left[ p \right]=\left[g\rho z\right] =  \frac{\mathrm{kg}}{\mathrm{m}\cdot\mathrm{s}^2} , \qquad
\left[ W \right] = \frac{\mathrm{kg}\cdot\mathrm{m}^2}{\mathrm{s}^3}, \qquad
\Rightarrow \qquad
\left[ \Lambda^m \right] = \frac{\mathrm{W}}{\mathrm{m}^5} 
\f]

In the case that \f$ \vec{\omega}=\left(\omega_r,0,0\right) \f$, then this reduces to
\f[
\Lambda^m = \frac{\overline{\omega_r}}{\overline{\rho}r^2\cos(\phi)}\left( \partial_{\lambda}\overline{\rho}\partial_{\phi}\overline{p} - \partial_{\phi}\overline{\rho}\partial_{\lambda}\overline{p} \right)
\f]


## Differentiation

The spherical differentiation methods ([spher-deriv]) use a land-avoiding stencil.
Land avoiding is achieved by simply using a non-centred stencil when appropriate. 
This is done to avoid having to specify field values at land cells, as this may introduce artificially steep velocity gradients that would confound derivatives, particularly with pressure and density, where there's no clear extension to land.


### Cartesian derivatives
The secondary differentation tools ([Cart-deriv]) simply apply the chain rule on the spherical differentiation methods.

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

[spher-deriv]: @ref spher_derivative_at_point "spherical derivatives"
[cart-deriv]: @ref Cart_derivative_at_point "Cartesian derivatives"
