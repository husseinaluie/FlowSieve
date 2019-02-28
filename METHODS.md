[TOC]
# Methods

This file will outline some of the computational methodologies.

## Cartesian derivatives

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


