/* ZM - Z's Mathematics Toolbox
 * Copyright (C) 1998 Tomomichi Sugihara (Zhidao)
 *
 * zm_ode_rkg - ordinary differential equation quadrature:
 * Runge-Kutta-Gill method.
 */

#ifndef __ZM_ODE_RKG_H__
#define __ZM_ODE_RKG_H__

/* NOTE: never include this header file in user programs. */

__BEGIN_DECLS

/* Runge-Kutta-Gill method is a variation of calssical
 * Runge-Kutta method. It is advantageous where less
 * workspace is required and rounding error is well absorbed.
 */
__EXPORT zODE *zODEInit_RKG(zODE *ode, int dim, int dummy, zVec (* f)(double,zVec,void*,zVec));
__EXPORT void zODEDestroy_RKG(zODE *ode);
__EXPORT zVec zODEUpdate_RKG(zODE *ode, double t, zVec x, double dt, void *util);

__END_DECLS

#endif /* __ZM_ODE_RKG_H__ */
