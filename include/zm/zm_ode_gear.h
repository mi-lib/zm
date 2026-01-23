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

/* Gear method is based on the Backward Differentiation Formula(BDF), one of the implicit linear multistep
 * methods. It has a stiffly stable property against stiff systems.
 *
 * Note that the initial series of variables should be set properly. They are set for zeros in default,
 * which might cause stable but largely different numerical behaviors.
 */
__ZM_EXPORT zODE *zODECreateGear(zODE *ode, int dim, int dummy, zVec (* f)(double,zVec,void*,zVec));
__ZM_EXPORT void zODEDestroyGear(zODE *ode);
__ZM_EXPORT zVec zODEUpdateGear(zODE *ode, double t, zVec x, double dt, void *util);

/*! \brief initialize history with a given vector.
 *
 * zODEInitHistoryGear() is particularly prepared for Gear method.
 * Since Gear method is affected by the choice of initial series of variable, it perhaps works better to
 * set them for the initial value than nothing.
 * This fills the initial history with \a x.
 * \return
 * zODEInitHistoryGear() returns no values.
 * \notes
 * Since it is specialized to Gear method, \a ode should be assigned for Gear and initialized in advance.
 */
__ZM_EXPORT void zODEInitHistoryGear(zODE *ode, zVec x);

__END_DECLS

#endif /* __ZM_ODE_GEAR_H__ */
