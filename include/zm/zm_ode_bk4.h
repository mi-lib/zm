/* ZM - Z's Mathematics Toolbox
 * Copyright (C) 1998 Tomomichi Sugihara (Zhidao)
 *
 * zm_ode_bk4 - ordinary differential equation quadrature:
 * Butcher-Kuntzmann method (Gauss method and Radau method).
 */

#ifndef __ZM_ODE_BK4_H__
#define __ZM_ODE_BK4_H__

/* NOTE: never include this header file in user programs. */

__BEGIN_DECLS

/* Butcher-Kuntzmann method is a two-step fourth-dimension
 * implicit Runge-Kutta method.
 * In implicit Runge-Kutta method, it is known that dimension
 * is no more than twice steps, and thus, BK4 is thought
 * to be the most practical.
 */
__EXPORT zODE *zODEInit_BK4(zODE *ode, int dim, int dummy, zVec (* f)(double,zVec,void*,zVec));
__EXPORT void zODEDestroy_BK4(zODE *ode);

/* 'BK4' is a dummy keyword, since Butcher-Kuntzmann method
 * is a general framework. As body implementations of solvers,
 * 'Gauss' and 'Radau' are prepared.
 */
__EXPORT zVec zODEUpdate_Gauss(zODE *ode, double t, zVec x, double dt, void *util);
__EXPORT zVec zODEUpdate_Radau(zODE *ode, double t, zVec x, double dt, void *util);

#define zODEInit_Gauss    zODEInit_BK4
#define zODEDestroy_Gauss zODEDestroy_BK4
#define zODEInit_Radau    zODEInit_BK4
#define zODEDestroy_Radau zODEDestroy_BK4

__END_DECLS

#endif /* __ZM_ODE_BK4_H__ */
