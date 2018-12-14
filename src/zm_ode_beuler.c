/* ZM - Z's Mathematics Toolbox
 * Copyright (C) 1998 Tomomichi Sugihara (Zhidao)
 *
 * zm_ode_beuler - ordinary differential equation quadrature:
 * backward Euler method & trapezoidal formula method.
 */

#include <zm/zm_ode.h>
#include <zm/zm_nle.h>

typedef struct{
  zVec v, x;
  double t, dt;
  void *util;
  int iter;
  zNLE nle;
} _zODE_BEuler;

/* backward Euler method */

static zVec _zODE_BEuler_Func(zVec x, zVec y, void *util);

/* (static)
 * _zODE_BEuler_Func
 * - algebraic equation for backward Euler method.
 */
zVec _zODE_BEuler_Func(zVec x, zVec y, void *util)
{
  zODE *ode = util;
  _zODE_BEuler *ws = ode->_ws;

  ode->sub( x, ws->x, y, ws->util );
  ode->f( ws->t, x, ws->util, ws->v );
  return ode->cat( y, -ws->dt, ws->v, y, ws->util );
}

/* zODEInit_BEuler
 * - initialize ODE solver based on backward Euler method.
 */
zODE* zODEInit_BEuler(zODE *ode, int dim, int iter, zVec (* f)(double,zVec,void*,zVec))
{
  _zODE_BEuler *ws;

  if( !( ws = zAlloc( _zODE_BEuler, 1 ) ) ||
      !( ws->v = zVecAlloc( dim ) ) ||
      !( ws->x = zVecAlloc( dim ) ) ||
      !zNLECreate( &ws->nle, dim, dim, 0, _zODE_BEuler_Func, NULL ) ){
    ZALLOCERROR();
    return NULL;
  }
  ode->f = f;
  ode->_ws = ws;
  ws->iter = iter;
  zNLEAssignNR( &ws->nle );
  return ode;
}

/* zODEDestroy_BEuler
 * - destroy ODE solver.
 */
void zODEDestroy_BEuler(zODE *ode)
{
  _zODE_BEuler *ws = ode->_ws;

  zNLEDestroy( &ws->nle );
  zVecFreeAO( 2, ws->v, ws->x );
  zFree( ws );
  ode->f = NULL;
}

/* zODEUpdate_BEuler
 * - directly integrate variable by ODE based on backward Euler method.
 */
zVec zODEUpdate_BEuler(zODE *ode, double t, zVec x, double dt, void *util)
{
  _zODE_BEuler *ws = ode->_ws;

  ws->t = t + dt;
  ws->dt = dt;
  ws->util = util;
  zVecCopyNC( x, ws->x );
  zNLESolve( &ws->nle, x, ode, zTOL, ws->iter, NULL );
  return x;
}

/* trapezoidal formula method */

/* zODEUpdate_TR
 * - directly integrate variable by ODE based on backward Euler method.
 */
zVec zODEUpdate_TR(zODE *ode, double t, zVec x, double dt, void *util)
{
  _zODE_BEuler *ws = ode->_ws;

  ws->t = t + dt;
  ws->dt = 0.5*dt;
  ws->util = util;
  zVecCopyNC( x, ws->x );
  ode->f( t, ws->x, ws->util, ws->v );  /* NOTE: ws->t != t */
  ode->cat( ws->x, ws->dt, ws->v, ws->x, ws->util ); /* NOTE: ws->dt != dt */
  zNLESolve( &ws->nle, x, ode, zTOL, ws->iter, NULL );
  return x;
}
