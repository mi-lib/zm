/* ZM - Z's Mathematics Toolbox
 * Copyright (C) 1998 Tomomichi Sugihara (Zhidao)
 *
 * zm_ode_beuler - ordinary differential equation quadrature:
 * backward Euler method.
 */

#ifndef __ZM_ODE_BEULER_H__
#define __ZM_ODE_BEULER_H__

/* NOTE: never include this header file in user programs. */

__BEGIN_DECLS

/* backward Euler method - easiest way of implicit solutions.
 */
__EXPORT zODE *zODEInit_BEuler(zODE *ode, int dim, int dummy, zVec (* f)(double,zVec,void*,zVec));
__EXPORT void zODEDestroy_BEuler(zODE *ode);
__EXPORT zVec zODEUpdate_BEuler(zODE *ode, double t, zVec x, double dt, void *util);

/* trapezoidal formula method - A-stable implicit solutions.
 */
#define zODEInit_TR zODEInit_BEuler
#define zODEDestroy_TR zODEDestroy_BEuler
__EXPORT zVec zODEUpdate_TR(zODE *ode, double t, zVec x, double dt, void *util);

__END_DECLS

#endif /* __ZM_ODE_BEULER_H__ */
