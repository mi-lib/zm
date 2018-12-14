/* ZM - Z's Mathematics Toolbox
 * Copyright (C) 1998 Tomomichi Sugihara (Zhidao)
 *
 * zm_ode_heun - ordinary differential equation quadrature:
 * Heun method.
 */

#ifndef __ZM_ODE_HEUN_H__
#define __ZM_ODE_HEUN_H__

/* NOTE: never include this header file in user programs. */

__BEGIN_DECLS

/* Here, the third-dimension Runge-Kutta method is named Heun method.
 * To be strict, it is one of the modified Euler methods by Heun.
 */
__EXPORT zODE *zODEInit_Heun(zODE *ode, int dim, int dummy, zVec (* f)(double,zVec,void*,zVec));
__EXPORT void zODEDestroy_Heun(zODE *ode);
__EXPORT zVec zODEUpdate_Heun(zODE *ode, double t, zVec x, double dt, void *util);

__END_DECLS

#endif /* __ZM_ODE_HEUN_H__ */
