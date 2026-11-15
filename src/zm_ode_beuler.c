/* ZM - Z's Mathematics Toolbox
 * Copyright (C) 1998 Tomomichi Sugihara (Zhidao)
 *
 * zm_ode_beuler - ordinary differential equation quadrature:
 * backward Euler method & trapezoidal formula method.
 */

#include <zm/zm_ode.h>
#include <zm/zm_nle.h>

typedef struct{
  zVec v, w, x, x0;
  double t, dt;
  void *util;
  int iter;
  zNLE nle;
} _zODE_BEuler;

/* backward Euler method */

/* algebraic equation for backward Euler method. */
static zVec _zODE_BEuler_Func(zVec x, zVec y, void *util)
{
  _zODE_BEuler *ws;

  ws = (_zODE_BEuler *)((zODE *)util)->_ws;
  ((zODE *)util)->sub( x, ws->x, y, ws->util );
  ((zODE *)util)->f( ws->t, x, ws->util, ws->v );
  return ((zODE *)util)->cat( y, -ws->dt, ws->v, y, ws->util );
}

/* create an ODE solver based on backward Euler method. */
zODE* zODECreateBEuler(zODE *ode, int dim, int iter, zVec (* f)(double,zVec,void*,zVec))
{
  _zODE_BEuler *ws;

  if( !( ws = zAlloc( _zODE_BEuler, 1 ) ) ||
      !( ws->v = zVecAlloc( dim ) ) ||
      !( ws->w = zVecAlloc( dim ) ) ||
      !( ws->x = zVecAlloc( dim ) ) ||
      !( ws->x0 = zVecAlloc( dim ) ) ||
      !zNLECreate( &ws->nle, dim, dim, _zODE_BEuler_Func, NULL ) ){
    ZALLOCERROR();
    return NULL;
  }
  ode->f = f;
  ode->_ws = ws;
  ws->iter = iter;
  zNLEAssignNR( &ws->nle );
  return ode;
}

/* destroy an ODE solver. */
void zODEDestroyBEuler(zODE *ode)
{
  _zODE_BEuler *ws;

  ws = (_zODE_BEuler *)ode->_ws;
  zNLEDestroy( &ws->nle );
  zVecFreeAtOnce( 4, ws->v, ws->w, ws->x, ws->x0 );
  zFree( ws );
  ode->f = NULL;
}

/* solve the nonlinear equation of backward Euler method or trapezoidal formula method from the initial
 * guess x, where the residual is to be relatively small to the scale s of the state and its change in a
 * time step. x is restored to ws->x0 if it fails. an explicit predictor such as x + dt f(t+dt,x) is not
 * used for the initial guess, since it can lead the iteration to a wrong solution for a stiff equation. */
static zVec _zODESolveBEuler(zODE *ode, zVec x, double s)
{
  _zODE_BEuler *ws;

  ws = (_zODE_BEuler *)ode->_ws;
  if( !zIsFinite( s ) || zNLESolve( &ws->nle, x, ode, zTOL, zTOL*s, ws->iter, NULL ) < 0 ){
    zVecCopyNC( ws->x0, x );
    return NULL;
  }
  return x;
}

/* directly integrate variable by ODE based on backward Euler method. */
zVec zODEUpdateBEuler(zODE *ode, double t, zVec x, double dt, void *util)
{
  _zODE_BEuler *ws;
  double s;

  ws = (_zODE_BEuler *)ode->_ws;
  ws->t = t + dt;
  ws->dt = dt;
  ws->util = util;
  zVecCopyNC( x, ws->x );
  zVecCopyNC( x, ws->x0 );
  s = _zODEImplicitScale( ode, t, x, dt, ws->v, ws->w, util );
  return _zODESolveBEuler( ode, x, s );
}

/* trapezoidal formula method */

/* directly integrate variable by ODE based on backward Euler method. */
zVec zODEUpdateTR(zODE *ode, double t, zVec x, double dt, void *util)
{
  _zODE_BEuler *ws;
  double s;

  ws = (_zODE_BEuler *)ode->_ws;
  ws->t = t + dt;
  ws->dt = 0.5*dt;
  ws->util = util;
  zVecCopyNC( x, ws->x );
  zVecCopyNC( x, ws->x0 );
  s = _zODEImplicitScale( ode, t, x, dt, ws->v, ws->w, util ); /* NOTE: ws->t != t */
  ode->cat( ws->x, ws->dt, ws->v, ws->x, ws->util ); /* NOTE: ws->dt != dt */
  return _zODESolveBEuler( ode, x, s );
}
