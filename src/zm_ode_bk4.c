/* ZM - Z's Mathematics Toolbox
 * Copyright (C) 1998 Tomomichi Sugihara (Zhidao)
 *
 * zm_ode_bk4 - ordinary differential equation quadrature:
 * Butcher-Kuntzmann method (Gauss method and Radau method).
 */

#include <zm/zm_ode.h>
#include <zm/zm_nle.h>

typedef struct{
  zVec v, x, xt, k;
  zVecStruct k1, k2;
  double t1, t2;
  double a11, a12, a21, a22;
  void *util;
  int iter;
  zNLE nle;
} _zODE_BK4;

static zVec _zODE_BK4_Func(zVec x, zVec y, void *util);

/* (static)
 * _zODE_BK4_Func
 * - algebraic equation for Butcher-Kuntzmann method.
 */
zVec _zODE_BK4_Func(zVec x, zVec y, void *util)
{
  zODE *ode = util;
  _zODE_BK4 *ws = ode->_ws;
  zVecStruct y1, y2;

  /* divide y vector and assign them to half-size vectors */
  zVecSetSize( &y1, zVecSizeNC(ws->x) );
  zVecBuf(&y1) = zVecBuf(y);
  zVecSetSize( &y2, zVecSizeNC(ws->x) );
  zVecBuf(&y2) = zVecBuf(y) + zVecSizeNC(ws->x);

  /* k1 */
  ode->cat( ws->x, ws->a11, &ws->k1, ws->xt, ws->util );
  ode->cat( ws->xt, ws->a12, &ws->k2, ws->xt, ws->util );
  ode->f( ws->t1, ws->xt, ws->util, ws->v );
  ode->sub( &ws->k1, ws->v, &y1, ws->util );
  /* k2 */
  ode->cat( ws->x, ws->a21, &ws->k1, ws->xt, ws->util );
  ode->cat( ws->xt, ws->a22, &ws->k2, ws->xt, ws->util );
  ode->f( ws->t2, ws->xt, ws->util, ws->v );
  ode->sub( &ws->k2, ws->v, &y2, ws->util );
  return y;
}

/* zODEInit_BK4
 * - initialize ODE solver based on Butcher-Kuntzmann method.
 */
zODE* zODEInit_BK4(zODE *ode, int dim, int iter, zVec (* f)(double,zVec,void*,zVec))
{
  _zODE_BK4 *ws;

  if( !( ws = zAlloc( _zODE_BK4, 1 ) ) ||
      !( ws->v = zVecAlloc( dim ) ) ||
      !( ws->x = zVecAlloc( dim ) ) ||
      !( ws->xt = zVecAlloc( dim ) ) ||
      !( ws->k = zVecAlloc( dim*2 ) ) ||
      !zNLECreate( &ws->nle, dim*2, dim*2, 0, _zODE_BK4_Func, NULL ) ){
    ZALLOCERROR();
    return NULL;
  }
  zVecSetSize( &ws->k1, dim );
  zVecBuf(&ws->k1) = zVecBuf(ws->k);
  zVecSetSize( &ws->k2, dim );
  zVecBuf(&ws->k2) = zVecBuf(ws->k) + dim;
  ws->iter = iter;
  ode->f = f;
  ode->_ws = ws;
  zNLEAssignNR( &ws->nle );
  return ode;
}

/* zODEDestroy_BK4
 * - destroy ODE solver.
 */
void zODEDestroy_BK4(zODE *ode)
{
  _zODE_BK4 *ws = ode->_ws;

  zNLEDestroy( &ws->nle );
  zVecFreeAO( 4, ws->v, ws->x, ws->xt, ws->k );
  zFree( ws );
  ode->f = NULL;
}

/* zODEUpdate_Gauss
 * - directly integrate variable by ODE based on Gauss method.
 */
zVec zODEUpdate_Gauss(zODE *ode, double t, zVec x, double dt, void *util)
{
  _zODE_BK4 *ws = ode->_ws;
  const double bk4_c = sqrt(3)/6;

  /* Butcher's array */
  ws->t1 = t + (0.5+bk4_c)*dt;
  ws->t2 = t + (0.5-bk4_c)*dt;
  ws->a11 = ws->a22 = 0.25 *dt;
  ws->a12 = (0.25+bk4_c)*dt;
  ws->a21 = (0.25-bk4_c)*dt;
  dt *= 0.5;

  ws->util = util;
  zVecCopyNC( x, ws->x );
  zVecClear( ws->k );
  zNLESolve( &ws->nle, ws->k, ode, zTOL, ws->iter, NULL );
  ode->cat( x, dt, &ws->k1, x, util );
  ode->cat( x, dt, &ws->k2, x, util );
  return x;
}

/* zODEUpdate_Radau
 * - directly integrate variable by ODE based on Radau method.
 */
zVec zODEUpdate_Radau(zODE *ode, double t, zVec x, double dt, void *util)
{
  _zODE_BK4 *ws = ode->_ws;

  /* Butcher's array */
  ws->t1 = t + dt/3.0;
  ws->t2 = t + dt;
  ws->a11 = -5.0 * ( ws->a12 = -1.0/12 *dt );
  ws->a21 =  3.0 * ( ws->a22 =  0.25   *dt );

  ws->util = util;
  zVecCopyNC( x, ws->x );
  zVecClear( ws->k );
  zNLESolve( &ws->nle, ws->k, ode, zTOL, ws->iter, NULL );
  ode->cat( x, ws->a21, &ws->k1, x, util );
  ode->cat( x, ws->a22, &ws->k2, x, util );
  return x;
}
