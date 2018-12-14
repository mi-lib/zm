/* ZM - Z's Mathematics Toolbox
 * Copyright (C) 1998 Tomomichi Sugihara (Zhidao)
 *
 * zm_ip_ipio - interpolation: interpolation-in-order.
 */

#include <zm/zm_ip.h>

/* ********************************************************** */
/* CLASS: zIPIO
 * interpolator-in-order class
 * this class realizes interpolation-in-order.
 * ********************************************************** */

/* zIPIOInit
 * - initialization of interpolation-in-order.
 */
void zIPIOInit(zIPIO *ip)
{
  ip->p_prev = ip->p_cur = ip->p_next = 0;
  ip->v_prev = ip->v_cur = 0;
  ip->dt = ip->dt_next = 1.0; /* dummy time step */
  ip->c[0] = ip->c[1] = ip->c[2] = ip->c[3] = 0;
}

/* zIPIOCreate
 * - creation of interpolation-in-order.
 */
void zIPIOCreate(zIPIO *ip, double p0)
{
  zIPIOInit( ip );
  ip->p_prev = ip->p_cur = ip->p_next = p0;
}

/* zIPIOSetNextValue
 * - setting of the next value on interpolation-in-order.
 */
void zIPIOSetNextValue(zIPIO *ip, double p, double dt)
{
  ip->p_prev = ip->p_cur;
  ip->p_cur = ip->p_next;
  ip->v_prev = ip->v_cur;
  ip->dt = ip->dt_next;

  ip->p_next = p;
  ip->dt_next = dt;

  ip->v_cur = 0.5 * ( (ip->p_next-ip->p_cur)/ip->dt_next
                    + (ip->p_cur-ip->p_prev)/ip->dt );
}

/* zIPIOUpdate
 * - updating of interpolation-in-order.
 */
void zIPIOUpdate(zIPIO *ip)
{
  double dp;

  ip->c[0] = ip->p_prev;
  ip->c[1] = ip->v_prev;
  dp = ip->p_cur - ip->p_prev;
  ip->c[2] = (3*dp/ip->dt-(ip->v_cur+2*ip->v_prev))/ip->dt;
  ip->c[3] = (-2*dp/ip->dt+(ip->v_cur+ip->v_prev))/(ip->dt*ip->dt);
}

/* zIPIOLinValue
 * - value on linear interpolation-in-order.
 */
double zIPIOLinValue(zIPIO *ip, double t)
{
  return ip->p_prev + ( ip->p_cur - ip->p_prev ) * t/ip->dt;
}

/* zIPIOSplineValue
 * - value on spline interpolation-in-order.
 */
double zIPIOSplineValue(zIPIO *ip, double t)
{
  return ip->c[0]+(ip->c[1]+(ip->c[2]+ip->c[3]*t)*t)*t;
}

/* ********************************************************** */
/* CLASS: zFerguson
 * Ferguson's curve
 * ********************************************************** */

/* zFergusonCreate
 * - creation of Ferguson s interpolation curve.
 */
void zFergusonCreate(zFerguson *ferg, double term, double x1, double v1, double x2, double v2)
{
  ferg->c[0] = x1;
  ferg->c[1] = v1 / term;
  ferg->c[2] = (3*(x2-x1)-2*v1+v2)/(term*term);
  ferg->c[3] = (2*(x1-x2)+v1-v2)/(term*term*term);
}

/* zFergusonValue
 * - value calculation of Ferguson s interpolation curve.
 */
double zFergusonValue(zFerguson *ferg, double t)
{
  return ( ( ferg->c[3]*t + ferg->c[2] )*t + ferg->c[1] )*t + ferg->c[0];
}
