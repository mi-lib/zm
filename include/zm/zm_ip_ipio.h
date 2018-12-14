/* ZM - Z's Mathematics Toolbox
 * Copyright (C) 1998 Tomomichi Sugihara (Zhidao)
 *
 * zm_ip_ipio - interpolation: interpolation-in-order.
 */

#ifndef __ZM_IP_IPIO_H__
#define __ZM_IP_IPIO_H__

/* NOTE: never include this header file in user programs. */

__BEGIN_DECLS

/* ********************************************************** */
/* CLASS: zIPIO
 * interpolator-in-order class
 * this class realizes interpolation-in-order.
 * ********************************************************** */

typedef struct{
  double p_prev, p_cur, p_next;
  double v_prev, v_cur;
  double dt, dt_next;
  double c[4]; /* coefficients */
} zIPIO;

#define zIPIODT(i) (i)->dt

/* METHOD:
 * zIPIOInit, zIPIOCreate, zIPIOSetNextValue, zIPIOUpdate
 * zIPIOLinValue, zIPIOSplineValue
 * - interpolation-in-order.
 *
 * 'zIPIO' realizes the interpolation-in-order; the interpolation
 * is done with values at three times - previous, current and
 * next steps. Updating the inner parameters of the interpolator
 * by setting new value at the next step, the segment between the
 * two control points, previous and current point, is interpolated.
 * #
 * 'zIPIO' has three control points and two time steps.
 * 'p_prev', 'p_cur' and 'p_next' are the previous, current and
 * next point respectively. And the time step 'dt' is a time
 * between 'p_prev' and 'p_cur', while 'dt_next' is a time
 * between 'p_cur' and 'p_next'.
 * Setting these 5 parameters and updating the inner state of
 * the instance of 'zIPIO', the value at the time from the
 * previous point can be calculated.
 * #
 * The procedure to use 'zIPIO' is as follows.
 *  1. Create the new instance 'ip' of zIPIO by calling 'zIPIOCreate()'.
 *     'p0' is the initial value for the interpolation.
 *  2. When the next value comes, set it with 'zIPIOSetNextValue()'.
 *     'p' is the new value and 'dt' is a time step between the
 *     current and next step.
 *  3. Call 'zIOIPUpdate()', and the inner parameters of 'ip'
 *     will be updated; throwing away the previous value,
 *     and sliding back the current and next value to the previous
 *     and current value respectively.
 *  4. The two methods for interpolation are available.
 *     One is the linear interpolation and the other is the spline
 *     interpolation. Each is called by 'zIPIOLinValue()' and
 *     'zIPIOSplineValue()' respectively.
 *     't' is a time width from the previous time.
 *  5. The instance of 'zIPIO', 'ip', can be intialized by using
 *     'zIPIOInit()'.
 * [RETURN VALUE]
 * 'zIPIOInit()', 'zIPIOCreate()', 'zIPIOSetNextValue()' and
 * 'zIPIOUpdate()' return no value.
 * #
 * 'zIPIOLinValue()' and 'zIPIOSplineValue()' return the value
 * interpolated.
 */
__EXPORT void zIPIOInit(zIPIO *ip);
__EXPORT void zIPIOCreate(zIPIO *ip, double p0);
__EXPORT void zIPIOSetNextValue(zIPIO *ip, double p, double dt);
__EXPORT void zIPIOUpdate(zIPIO *ip);
__EXPORT double zIPIOSplineValue(zIPIO *ip, double t);
__EXPORT double zIPIOLinValue(zIPIO *ip, double t);
#define zIPIOValue(i,t)  zIPIOSplineValue(i,t)

/* ********************************************************** */
/* CLASS: zFerguson
 * Ferguson's curve
 * ********************************************************** */

typedef struct{
  double c[4];
} zFerguson;

/* METHOD:
 * zFergusonCreate, zFergusonValue
 * - Ferguson s interpolation curve.
 *
 * 'zFerguson' realizes Ferguson s curve, which is defined as:
 *   x(t)=x(0)*H_00(t/T)+x(T)*H_01(t/T)+x'(0)*H_10(t/T)+x'(T)*H_11(t/T)
 * where
 *   H_00(s) = 2 s^3 - 3 s^2 + 1 = (s-1)^2*(2s+1)
 *   H_01(s) =-2 s^3 + 3 s^2     =- s^2*(2s-3)
 *   H_10(s) =   s^3 - 2 s^2 + s =  s(s-1)^2
 *   H_11(s) =   s^3 -   s^2     =  s^2(s-1)
 * #
 * An advantage of it is that one can specify the boundary condition
 * in both position and velocity.
 * #
 * 'zFergusonCreate()' calculates the coefficients of the curve
 * which satisfies a given set of 'x1', 'v1', 'x2' and 'v2'.
 * 'term' is the terminal time(the time when the value has to be
 * 'x2').
 * #
 * 'zFergusonValue()' calculates the value at each moment 't'.
 * [RETURN VALUE]
 * 'zFergusonCreate()' returns no value.
 * 'zFergusonValue()' returns the value computed.
 */
__EXPORT void zFergusonCreate(zFerguson *ferg, double term, double x1, double v1, double x2, double v2);
__EXPORT double zFergusonValue(zFerguson *ferg, double t);

__END_DECLS

#endif /* __ZM_IP_IPIO_H__ */
