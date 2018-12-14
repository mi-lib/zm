/* ZM - Z's Mathematics Toolbox
 * Copyright (C) 1998 Tomomichi Sugihara (Zhidao)
 *
 * zm_ode_dc - ordinary differential equation quadrature:
 * deferred correction with Richardson's extrapolation.
 */

#ifndef __ZM_ODE_DC_H__
#define __ZM_ODE_DC_H__

__BEGIN_DECLS

/*!
 */
__EXPORT zODE *zODEInitDC(zODE *ode, int dim, int step, zVec (* f)(double,zVec,void*,zVec));
__EXPORT void zODEDestroyDC(zODE *ode);
__EXPORT zVec zODEUpdateDC(zODE *ode, double t, zVec x, double dt, void *util);

__END_DECLS

#endif /* __ZM_ODE_DC_H__ */
