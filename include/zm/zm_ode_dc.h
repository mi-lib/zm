/* ZM - Z's Mathematics Toolbox
 * Copyright (C) 1998 Tomomichi Sugihara (Zhidao)
 *
 * zm_ode_dc - ordinary differential equation quadrature:
 * deferred correction with Richardson's extrapolation.
 */

#ifndef __ZM_ODE_DC_H__
#define __ZM_ODE_DC_H__

__BEGIN_DECLS

/* Deferred correction with Richardson's extrapolation.
 * zODEUpdateDC() returns the null pointer without updating x if dt is zero, as zODEUpdate().
 * The multistep methods "Adams" and "Gear" cannot be combined with the deferred correction, since the
 * history is mixed up with the full time step and the half time steps of the trials. The embedded
 * Runge-Kutta methods "RKF45", "CK45" and "DP45" can be combined with it, but are not recommended,
 * since they control the step size by themselves, and the deferred correction only raises the order
 * of accuracy by one at the cost of about three times as many evaluations. */
__ZM_EXPORT zODE *zODECreateDC(zODE *ode, int dim, int step, zVec (* f)(double,zVec,void*,zVec));
__ZM_EXPORT void zODEDestroyDC(zODE *ode);
__ZM_EXPORT zVec zODEUpdateDC(zODE *ode, double t, zVec x, double dt, void *util);

__END_DECLS

#endif /* __ZM_ODE_DC_H__ */
