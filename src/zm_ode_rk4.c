/* ZM - Z's Mathematics Toolbox
 * Copyright (C) 1998 Tomomichi Sugihara (Zhidao)
 *
 * zm_ode_rk4 - ordinary differential equation quadrature:
 * classical Runge-Kutta method.
 */

#include <zm/zm_ode.h>

typedef struct{
  zVec x, k[4];
} _zODE_RK4;

/* zODEInit_RK4
 * - initialize ODE solver based on classical Runge-Kutta method.
 */
zODE* zODEInit_RK4(zODE *ode, int dim, int dummy, zVec (* f)(double,zVec,void*,zVec))
{
  _zODE_RK4 *ws;

  if( !( ws = zAlloc( _zODE_RK4, 1 ) ) ||
      !( ws->x = zVecAlloc( dim ) ) /* incremental vector */ ||
      !( ws->k[0] = zVecAlloc( dim ) ) /* step-1 vector */ ||
      !( ws->k[1] = zVecAlloc( dim ) ) /* step-2 vector */ ||
      !( ws->k[2] = zVecAlloc( dim ) ) /* step-3 vector */ ||
      !( ws->k[3] = zVecAlloc( dim ) ) /* step-4 vector */ ){
    ZALLOCERROR();
    return NULL;
  }
  ode->f = f;
  ode->_ws = ws;
  return ode;
}

/* zODEDestroy_RK4
 * - destroy ODE solver.
 */
void zODEDestroy_RK4(zODE *ode)
{
  _zODE_RK4 *ws;

  ws = ode->_ws;
  zVecFreeAO( 5, ws->x, ws->k[0], ws->k[1], ws->k[2], ws->k[3] );
  zFree( ode->_ws );
  ode->f = NULL;
}

/* zODEUpdate_RK4
 * - directly integrate variable by ODE based on classical Runge-Kutta method.
 */
zVec zODEUpdate_RK4(zODE *ode, double t, zVec x, double dt, void *util)
{
  _zODE_RK4 *ws;
  double dt1, dt2, dt3;

  ws = ode->_ws;
  dt1 = dt * 0.5;
  dt2 = dt / 6;
  dt3 = dt2 * 2;
  ode->f( t, x, util, ws->k[0] );
  ode->cat( x, dt1, ws->k[0], ws->x, util );
  ode->f( t+dt1, ws->x, util, ws->k[1] );
  ode->cat( x, dt1, ws->k[1], ws->x, util );
  ode->f( t+dt1, ws->x, util, ws->k[2] );
  ode->cat( x, dt, ws->k[2], ws->x, util );
  ode->f( t+dt, ws->x, util, ws->k[3] );

  ode->cat( x, dt2, ws->k[0], x, util );
  ode->cat( x, dt3, ws->k[1], x, util );
  ode->cat( x, dt3, ws->k[2], x, util );
  ode->cat( x, dt2, ws->k[3], x, util );
  return x;
}
