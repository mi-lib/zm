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

/* Some solvers in the class of embedded Runge-Kutta methods, which automatically modify quantized
 * integration time dt in accordance with the error estimated, are provided. */

/* Runge-Kutta-Fehlberg method is forth-dimension-embedded six-step fifth-dimension Runge-Kutta method. */
__ZM_EXPORT zODE *zODECreateRKF45(zODE *ode, int dim, int dummy, zVec (* f)(double,zVec,void*,zVec));
__ZM_EXPORT void zODEDestroyRKF45(zODE *ode);
__ZM_EXPORT zVec zODEUpdateRKF45(zODE *ode, double t, zVec x, double dt, void *util);

/* Cash-Karp method is forth-dimension-embedded six-step fifth-dimension Runge-Kutta method. */
__ZM_EXPORT zODE *zODECreateCK45(zODE *ode, int dim, int dummy, zVec (* f)(double,zVec,void*,zVec));
__ZM_EXPORT void zODEDestroyCK45(zODE *ode);
__ZM_EXPORT zVec zODEUpdateCK45(zODE *ode, double t, zVec x, double dt, void *util);

/* Dormand-Prince method is forth-dimension-embedded seven-step fifth-dimension Runge-Kutta method. */
__ZM_EXPORT zODE *zODECreateDP45(zODE *ode, int dim, int dummy, zVec (* f)(double,zVec,void*,zVec));
__ZM_EXPORT void zODEDestroyDP45(zODE *ode);
__ZM_EXPORT zVec zODEUpdateDP45(zODE *ode, double t, zVec x, double dt, void *util);

__END_DECLS

#endif /* __ZM_ODE_ERK_H__ */
