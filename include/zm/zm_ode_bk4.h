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

/* Butcher-Kuntzmann method is a two-step fourth-dimension implicit Runge-Kutta method.
 * In implicit Runge-Kutta method, it is known that dimension is no more than twice steps, and thus,
 * BK4 is thought to be the most practical.
 */
__ZM_EXPORT zODE *zODECreateBK4(zODE *ode, int dim, int dummy, zVec (* f)(double,zVec,void*,zVec));
__ZM_EXPORT void zODEDestroyBK4(zODE *ode);

/* 'BK4' is a dummy keyword, since Butcher-Kuntzmann method is a general framework.
 * As body implementations of solvers, 'Gauss' and 'Radau' are prepared. */

/* Gauss method */
#define zODECreateGauss  zODECreateBK4
#define zODEDestroyGauss zODEDestroyBK4
__ZM_EXPORT zVec zODEUpdateGauss(zODE *ode, double t, zVec x, double dt, void *util);

/* Radau method */
#define zODECreateRadau  zODECreateBK4
#define zODEDestroyRadau zODEDestroyBK4
__ZM_EXPORT zVec zODEUpdateRadau(zODE *ode, double t, zVec x, double dt, void *util);

__END_DECLS

#endif /* __ZM_ODE_BK4_H__ */
