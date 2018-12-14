/* ZM - Z's Mathematics Toolbox
 * Copyright (C) 1998 Tomomichi Sugihara (Zhidao)
 *
 * zm_ode_rk4 - ordinary differential equation quadrature:
 * classical Runge-Kutta method.
 */

#ifndef __ZM_ODE_RK4_H__
#define __ZM_ODE_RK4_H__

/* NOTE: never include this header file in user programs. */

__BEGIN_DECLS

/* classical Runge-Kutta method is also known as
 * forth-dimension Runge-Kutta method.
 * It assures fifth-order error.
 * Four is the maximum step number which coincides with
 * the dimention of error, so that it is thought to be
 * the most practical, while some objectives exist.
 */
__EXPORT zODE *zODEInit_RK4(zODE *ode, int dim, int dummy, zVec (* f)(double,zVec,void*,zVec));
__EXPORT void zODEDestroy_RK4(zODE *ode);
__EXPORT zVec zODEUpdate_RK4(zODE *ode, double t, zVec x, double dt, void *util);

__END_DECLS

#endif /* __ZM_ODE_RK4_H__ */
