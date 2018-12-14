/* ZM - Z's Mathematics Toolbox
 * Copyright (C) 1998 Tomomichi Sugihara (Zhidao)
 *
 * zm_ode_heun - ordinary differential equation quadrature:
 * Heun method.
 */

#include <zm/zm_ode.h>

typedef struct{
  zVec x, k[3];
} _zODE_Heun;

/* zODEInit_Heun
 * - initialize ODE solver based on Heun method.
 */
zODE* zODEInit_Heun(zODE *ode, int dim, int dummy, zVec (* f)(double,zVec,void*,zVec))
{
  _zODE_Heun *ws;

  if( !( ws = zAlloc( _zODE_Heun, 1 ) ) ||
      !( ws->x = zVecAlloc( dim ) ) /* incremental vector */ ||
      !( ws->k[0] = zVecAlloc( dim ) ) /* step-1 vector */ ||
      !( ws->k[1] = zVecAlloc( dim ) ) /* step-2 vector */ ||
      !( ws->k[2] = zVecAlloc( dim ) ) /* step-3 vector */ ){
    ZALLOCERROR();
    return NULL;
  }
  ode->f = f;
  ode->_ws = ws;
  return ode;
}

/* zODEDestroy_Heun
 * - destroy ODE solver.
 */
void zODEDestroy_Heun(zODE *ode)
{
  _zODE_Heun *ws;

  ws = ode->_ws;
  zVecFreeAO( 4, ws->x, ws->k[0], ws->k[1], ws->k[2] );
  zFree( ode->_ws );
  ode->f = NULL;
}

/* zODEUpdate_Heun
 * - directly integrate variable by ODE based on classical Runge-Kutta method.
 */
zVec zODEUpdate_Heun(zODE *ode, double t, zVec x, double dt, void *util)
{
  _zODE_Heun *ws;

  ws = ode->_ws;
  ode->f( t, x, util, ws->k[0] );
  ode->cat( x, dt/3, ws->k[0], ws->x, util );
  ode->f( t+dt/3, ws->x, util, ws->k[1] );
  ode->cat( x, dt*2/3, ws->k[1], ws->x, util );
  ode->f( t+dt*2/3, ws->x, util, ws->k[2] );

  ode->cat( x, 0.25*dt, ws->k[0], x, util );
  ode->cat( x, 0.75*dt, ws->k[2], x, util );
  return x;
}
