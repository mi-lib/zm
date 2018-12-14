/* ZM - Z's Mathematics Toolbox
 * Copyright (C) 1998 Tomomichi Sugihara (Zhidao)
 *
 * zm_ode_adams - ordinary differential equation quadrature:
 * Adams's linear multistep method.
 */

#ifndef __ZM_ODE_ADAMS_H__
#define __ZM_ODE_ADAMS_H__

/* NOTE: never include this header file in user programs. */

__BEGIN_DECLS

/* Predictor-Corrector method with a combination of
 * Adams=Bashforth / Adams=Moulton formulae.
 */
__EXPORT zODE *zODEInit_Adams(zODE *ode, int dim, int step, zVec (* f)(double,zVec,void*,zVec));
__EXPORT void zODEDestroy_Adams(zODE *ode);
__EXPORT zVec zODEUpdate_Adams(zODE *ode, double t, zVec x, double dt, void *util);

__END_DECLS

#endif /* __ZM_ODE_ADAMS_H__ */
