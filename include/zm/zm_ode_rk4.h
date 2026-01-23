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

/* classical Runge-Kutta method is also known as forth-dimension explicit Runge-Kutta method.
 * It assures fifth-order error.
 * Four is the maximum step number which coincides with the dimention of error, so that it is thought to be
 * the most practical (while some claims exist).
 */
__ZM_EXPORT zODE *zODECreateRK4(zODE *ode, int dim, int dummy, zVec (* f)(double,zVec,void*,zVec));
__ZM_EXPORT void zODEDestroyRK4(zODE *ode);
__ZM_EXPORT zVec zODEUpdateRK4(zODE *ode, double t, zVec x, double dt, void *util);

__END_DECLS

#endif /* __ZM_ODE_RK4_H__ */
