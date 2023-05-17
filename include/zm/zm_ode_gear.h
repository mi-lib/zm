/* ZM - Z's Mathematics Toolbox
 * Copyright (C) 1998 Tomomichi Sugihara (Zhidao)
 *
 * zm_ode_gear - ordinary differential equation quadrature:
 * Gear method.
 */

#ifndef __ZM_ODE_GEAR_H__
#define __ZM_ODE_GEAR_H__

/* NOTE: never include this header file in user programs. */

__BEGIN_DECLS

/* Gear method is based on the Backward Differentiation
 * Formula(BDF), one of the implicit linear multistep methods.
 * It has a stiffly stable property against stiff systems.
 *
 * Note that the initial series of variables should be set
 * properly. They are set for zeros in default, which may
 * cause stable but largely different numerical behaviors.
 */
__ZM_EXPORT zODE *zODEInit_Gear(zODE *ode, int dim, int dummy, zVec (* f)(double,zVec,void*,zVec));
__ZM_EXPORT void zODEDestroy_Gear(zODE *ode);
__ZM_EXPORT zVec zODEUpdate_Gear(zODE *ode, double t, zVec x, double dt, void *util);

/*! \brief initialize history with a given vector.
 *
 * zODEInitHist_Gear() is particularly prepared for Gear method.
 * Since Gear method is affected by the choice of initial series
 * of variable, it perhaps works better to set them for the initial
 * value than nothing.
 * Call this function, and it fills the initial history with \a x.
 * \return
 * zODEInitHist_Gear() returns no values.
 * \notes
 * Since it is specialized to Gear method, \a ode should be assigned
 * for Gear and initialized in advance.
 */
__ZM_EXPORT void zODEInitHist_Gear(zODE *ode, zVec x);

__END_DECLS

#endif /* __ZM_ODE_GEAR_H__ */
