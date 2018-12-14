/* ZM - Z's Mathematics Toolbox
 * Copyright (C) 1998 Tomomichi Sugihara (Zhidao)
 *
 * zm_ode_erk - ordinary differential equation quadrature:
 * embedded Runge-Kutta method.
 */

#ifndef __ZM_ODE_ERK_H__
#define __ZM_ODE_ERK_H__

/* NOTE: never include this header file in user programs. */

__BEGIN_DECLS

/* Some solvers in the class of embedded Runge-Kutta methods,
 * which automatically modify quantized integration time dt in
 * accordance with the error estimated, are provided.
 */

/* Runge-Kutta-Fehlberg method is forth-dimension-embedded
 * six-step fifth-dimension Runge-Kutta method.
 */
__EXPORT zODE *zODEInit_RKF45(zODE *ode, int dim, int dummy, zVec (* f)(double,zVec,void*,zVec));
__EXPORT void zODEDestroy_RKF45(zODE *ode);
__EXPORT zVec zODEUpdate_RKF45(zODE *ode, double t, zVec x, double dt, void *util);

/* Cash-Karp method is forth-dimension-embedded
 * six-step fifth-dimension Runge-Kutta method.
 */
__EXPORT zODE *zODEInit_CK45(zODE *ode, int dim, int dummy, zVec (* f)(double,zVec,void*,zVec));
__EXPORT void zODEDestroy_CK45(zODE *ode);
__EXPORT zVec zODEUpdate_CK45(zODE *ode, double t, zVec x, double dt, void *util);

/* Dormand-Prince method is forth-dimension-embedded
 * seven-step fifth-dimension Runge-Kutta method.
 */
__EXPORT zODE *zODEInit_DP45(zODE *ode, int dim, int dummy, zVec (* f)(double,zVec,void*,zVec));
__EXPORT void zODEDestroy_DP45(zODE *ode);
__EXPORT zVec zODEUpdate_DP45(zODE *ode, double t, zVec x, double dt, void *util);

__END_DECLS

#endif /* __ZM_ODE_ERK_H__ */
