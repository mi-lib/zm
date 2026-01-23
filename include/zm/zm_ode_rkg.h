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

/* Runge-Kutta-Gill method is a variation of calssical Runge-Kutta method.
 * It is advantageous where less workspace is required and rounding error is well absorbed.
 */
__ZM_EXPORT zODE *zODECreateRKG(zODE *ode, int dim, int dummy, zVec (* f)(double,zVec,void*,zVec));
__ZM_EXPORT void zODEDestroyRKG(zODE *ode);
__ZM_EXPORT zVec zODEUpdateRKG(zODE *ode, double t, zVec x, double dt, void *util);

__END_DECLS

#endif /* __ZM_ODE_RKG_H__ */
