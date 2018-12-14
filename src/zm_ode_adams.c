/* ZM - Z's Mathematics Toolbox
 * Copyright (C) 1998 Tomomichi Sugihara (Zhidao)
 *
 * zm_ode_adams - ordinary differential equation quadrature:
 * Adams's linear multistep method.
 */

#include <zm/zm_ode.h>

/* Adams's Linear Multistep method */
typedef struct{
  int step; /* multistep number */
  zVec dx; /* incremental vector */
  zVec wb, wm; /* weight vector */
  zVec x1, x2; /* working memory for PC method */
  zVecRing hist; /* history of differential values */
} _zODE_Adams;

static zODE *_zODEAlloc_Adams(zODE *ode, int dim, int step, zVec (* f)(double,zVec,void*,zVec));
static void _zODESetWeight_Adams(zVec w, int step, double g0);
static zVec _zODEUpdate_Adams(zODE *ode, double t, zVec x, zVec xorg, zVec xnew, zVec w, double dt, void *util);

/* _zODEAlloc_Adams
 * - allocate working space for zODE_Adams.
 *  (static)
 */
zODE *_zODEAlloc_Adams(zODE *ode, int dim, int step, zVec (* f)(double,zVec,void*,zVec))
{
  _zODE_Adams *ws;

  if( !( ws = zAlloc( _zODE_Adams, 1 ) ) ||
      !( ws->dx = zVecAlloc( dim ) ) ||
      !( ws->wb = zVecAlloc( step ) ) ||
      !( ws->wm = zVecAlloc( step ) ) ||
      !( ws->x1 = zVecAlloc( dim ) ) ||
      !( ws->x2 = zVecAlloc( dim ) ) ||
      !zVecRingAlloc( &ws->hist, dim, step ) ){
    ZALLOCERROR();
    return NULL;
  }
  ws->step = step;
  ode->f = f;
  ode->_ws = ws;
  return ode;
}

/* _zODESetWeight_Adams
 * - set weighting coefficients for linear multistep method.
 *   (static)
 */
void _zODESetWeight_Adams(zVec w, int step, double g0)
{
  register int i, j;
  zVec g;

  if( !( g = zVecAlloc( step ) ) ){
    ZALLOCERROR();
    return;
  }
  for( zVecSetElem( g, 0, 1 ), i=1; i<step; i++ ){
    zVecSetElem( g, i, g0 );
    for( j=1; j<=i; j++ )
      zVecElem(g,i) -= zVecElem(g,i-j) / (j+1);
  }
  for( i=0; i<step; i++ ){
    for( j=i; j<step; j++ )
      zVecElem(w,i) += zCombi(j,i) * zVecElem(g,j);
    if( zIsOdd(i) ) zVecElem(w,i) = -zVecElem(w,i);
  }
  zVecFree( g );
}

/* (static)
 * _zODEUpdate_Adams
 * - inner computation of Adams's method.
 */
zVec _zODEUpdate_Adams(zODE *ode, double t, zVec x, zVec xorg, zVec xnew, zVec w, double dt, void *util)
{
  _zODE_Adams *ws;

  ws = ode->_ws;
  ode->f( t, xorg, util, *zRingHead(&ws->hist) );
  zVecRingLS( ws->dx, w, &ws->hist );
  return ode->cat( x, dt, ws->dx, xnew, util );
}

/* zODEInit_Adams
 * - initialize ODE solver based on Predictor-Corrector method
 *   with Adams=Bashforth / Adams=Moulton's formulae.
 */
zODE *zODEInit_Adams(zODE *ode, int dim, int step, zVec (* f)(double,zVec,void*,zVec))
{
  _zODE_Adams *ws;

  if( !( ode = _zODEAlloc_Adams( ode, dim, step, f ) ) )
    return NULL;
  ws = ode->_ws;
  _zODESetWeight_Adams( ws->wb, ws->step, 1 );   /* AB formula */
  _zODESetWeight_Adams( ws->wm, ws->step, 0 ); /* AM formula */
  return ode;
}

/* zODEDestroy_Adams
 * - destroy ODE solver.
 */
void zODEDestroy_Adams(zODE *ode)
{
  _zODE_Adams *ws;

  ws = ode->_ws;
  ws->step = 0;
  zVecFreeAO( 5, ws->dx, ws->wb, ws->wm, ws->x1, ws->x2 );
  zVecRingFree( &ws->hist );
  ode->f = NULL;
  zFree( ode->_ws );
}

/* zODEUpdate_Adams
 * - directly integrate variable by ODE based on
 *   Predictor-Corrector method.
 */
zVec zODEUpdate_Adams(zODE *ode, double t, zVec x, double dt, void *util)
{
  _zODE_Adams *ws;
  double t2;
  zVec xnew, xold;
  register int i, iter = 0;

  ws = ode->_ws;
  t2 = t + dt;
  xold = ws->x1;
  xnew = ws->x2;
  /* predictor */
  zRingDecHead( &ws->hist ); /* increment ring */
  _zODEUpdate_Adams( ode, t,  x, x, xold, ws->wb, dt, util );
  /* corrector by successive substitution method */
  zRingDecHead( &ws->hist ); /* increment ring */
  ZITERINIT( iter );
  for( i=0; i<iter; i++ ){
    _zODEUpdate_Adams( ode, t2, x, xold, xnew, ws->wm, dt, util );
    if( zIsTiny( zVecDist( xold, xnew ) ) ) goto UPDATE;
    zSwap( zVec, xold, xnew );
  }
  ZITERWARN( iter );
 UPDATE:
  zRingIncHead( &ws->hist ); /* reset ring */
  zVecCopy( xnew, x );
  return x;
}
