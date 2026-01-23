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

/* backward Euler method - easiest way of implicit solutions. */
__ZM_EXPORT zODE *zODECreateBEuler(zODE *ode, int dim, int dummy, zVec (* f)(double,zVec,void*,zVec));
__ZM_EXPORT void zODEDestroyBEuler(zODE *ode);
__ZM_EXPORT zVec zODEUpdateBEuler(zODE *ode, double t, zVec x, double dt, void *util);

/* trapezoidal formula method - A-stable implicit solutions. */
#define zODECreateTR zODECreateBEuler
#define zODEDestroyTR zODEDestroyBEuler
__ZM_EXPORT zVec zODEUpdateTR(zODE *ode, double t, zVec x, double dt, void *util);

__END_DECLS

#endif /* __ZM_ODE_BEULER_H__ */
