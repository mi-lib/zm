/* ZM - Z's Mathematics Toolbox
 * Copyright (C) 1998 Tomomichi Sugihara (Zhidao)
 *
 * zm_ode_euler - ordinary differential equation quadrature:
 * Euler method.
 */

#ifndef __ZM_ODE_EULER_H__
#define __ZM_ODE_EULER_H__

/* NOTE: never include this header file in user programs. */

__BEGIN_DECLS

/* Euler method is too known to describe:-)
 */
__EXPORT zODE *zODEInit_Euler(zODE *ode, int dim, int dummy, zVec (* f)(double,zVec,void*,zVec));
__EXPORT void zODEDestroy_Euler(zODE *ode);
__EXPORT zVec zODEUpdate_Euler(zODE *ode, double t, zVec x, double dt, void *util);

__END_DECLS

#endif /* __ZM_ODE_EULER_H__ */
