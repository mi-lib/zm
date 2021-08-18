/* ZM - Z's Mathematics Toolbox
 * (C)Copyright Zhidao since 1998, all rights are reserved.
 */

/*!
 * \mainpage

 ZM is a collection of numerical computation and optimization
 tools including:
 - random number generator (Mersenne twister)
 - basic statistics tools
 - special functions (Gauss's error, gamma, beta, Bessel)
 - complex number (arithmetics, quadratic/cubic equation solvers)
 - vectors and matrices of real numbers and complex numbers
 - polynomial expression and equation solver
 - linear equation solvers (generalized, tridiagonal, Lyapunov)
 - matrix decompositions (LU, Choresky, LQ/QR)
 - eigenvalue analysis
 - finite Fourier series and Fast Fourier Transform
 - graph search and rapidly-explored random tree.
 - kd-tree
 - vector sequence and interpolations (linear, Lagrange, spline,
   Akima, PCHIP, polynomial fitting, NURBS)
 - numerical integrator
 - optimization solvers (line search, linear programming, linear
   complementary programming, quadratic programming, nonlinear
   programming based on descent/non-descent methods)
 - nonlinear equation solvers
 - first/second order ordinary differential equation quadratures
   (Runge-Kutta family, embedded Runge-Kutta family, Adams, Gear,
   backward Euler, Butcher-Kuntzmann family, deferred correction,
   Leapflog)
 - multiple classification analysis (K-means, GMM)
 - nonlinear oscillator
 - parser of mathematical expression
 */

#ifndef __ZM_H__
#define __ZM_H__

#include <zm/zm_rand.h>
#include <zm/zm_sf.h>
#include <zm/zm_mca.h>
#include <zm/zm_pex.h>
#include <zm/zm_fourier.h>
#include <zm/zm_nle.h>
#include <zm/zm_ip.h>
#include <zm/zm_nurbs.h>
#include <zm/zm_data.h>
#include <zm/zm_opt.h>
#include <zm/zm_ode.h>
#include <zm/zm_intg.h>
#include <zm/zm_oscil.h>

#endif /* __ZM_H__ */
