/* ZM - Z's Mathematics Toolbox
 * Copyright (C) 1998 Tomomichi Sugihara (Zhidao)
 *
 * zm_ode_gear - ordinary differential equation quadrature:
 * Gear method.
 */

#include <zm/zm_ode.h>
#include <zm/zm_nle.h>

typedef struct{
  double t, dt;
  zVec v, vp, a;
  double b;
  zVecRing hist;
  void *util;
  zNLE nle;
} _zODE_Gear;

/* algebraic equation for Gear method. */
static zVec _zODE_Gear_Func(zVec x, zVec y, void *util)
{
  zODE *ode;
  _zODE_Gear *ws;
  int i;

  ode = (zODE *)util;
  ws = (_zODE_Gear *)ode->_ws;
  zVecCopyNC( x, y );
  for( i=0; i<zVecSizeNC(ws->a); i++ )
    ode->cat( y, zVecElem(ws->a,i), *zRingElem(&ws->hist,i), y, ws->util );
  ode->f( ws->t+ws->dt, x, ws->util, ws->v );
  ode->cat( y, -ws->b*ws->dt, ws->v, y, ws->util );
  return y;
}

/* Jacobian matrix of the algebraic equation for Gear method, namely, I - b dt df/dx, where df/dx is
 * computed by the central difference of f. the difference of the residual itself is not used, since
 * the history in the residual can be much larger than the state, e.g., for a stiff equation, and the
 * difference of the residual is lost in the round-off error of the history. */
static zMat _zODE_Gear_Jacobi(const zVec x, zMat j, void *util)
{
  zODE *ode;
  _zODE_Gear *ws;
  int i, k;
  double org, h, xmax;

  ode = (zODE *)util;
  ws = (_zODE_Gear *)ode->_ws;
  xmax = zVecElemAbsMax( x, NULL );
  for( i=0; i<zVecSizeNC(x); i++ ){
    h = zOptDMDiffStep( x, i, xmax, ZM_OPT_DM_EPS, ZM_OPT_DM_DIFF_SCALE_RATIO );
    org = zVecElemNC(x,i);
    zVecSetElemNC( x, i, org + h );
    ode->f( ws->t+ws->dt, x, ws->util, ws->v );
    zVecSetElemNC( x, i, org - h );
    ode->f( ws->t+ws->dt, x, ws->util, ws->vp );
    zVecSetElemNC( x, i, org );
    for( k=0; k<zVecSizeNC(x); k++ )
      zMatSetElemNC( j, k, i, -0.5*ws->b*ws->dt*( zVecElemNC(ws->v,k) - zVecElemNC(ws->vp,k) ) / h );
    zMatElemNC(j,i,i) += 1;
  }
  return j;
}

/* create an ODE solver based on Gear method. */
zODE* zODECreateGear(zODE *ode, int dim, int step, zVec (* f)(double,zVec,void*,zVec))
{
  _zODE_Gear *ws;
  double a[][6] = {
    { -1, 0, 0, 0, 0, 0 },
    { -4.0/3, 1.0/3, 0, 0, 0, 0 },
    { -18.0/11, 9.0/11, -2.0/11, 0, 0, 0 },
    { -48.0/25, 36.0/25, -16.0/25, 3.0/25, 0, 0 },
    { -300.0/137, 300.0/137, -200.0/137, 75.0/137, -12.0/137, 0 },
    { -120.0/49, 150.0/49, -400.0/147, 75.0/49, -24.0/49, 10.0/147 },
  };
  double b[] = {
    1, 2.0/3, 6.0/11, 12.0/25, 60.0/137, 20.0/49,
  };

  if( !( ws = zAlloc( _zODE_Gear, 1 ) ) ){
    ZALLOCERROR();
    return NULL;
  }
  if( step < 1 ){
    ZRUNWARN( ZM_WARN_ODE_GEAR1, step );
    step = 1;
  } else if( step > 6 ){
    ZRUNWARN( ZM_WARN_ODE_GEAR2, step );
    step = 6;
  }
  if( !( ws->v = zVecAlloc( dim ) ) ||
      !( ws->vp = zVecAlloc( dim ) ) ||
      !( ws->a = zVecCloneArray( a[step-1], step ) ) ||
      !zVecRingAlloc( &ws->hist, dim, step ) ||
      !zNLECreate( &ws->nle, dim, dim, _zODE_Gear_Func, _zODE_Gear_Jacobi ) ){
    ZALLOCERROR();
    return NULL;
  }
  ws->b = b[step-1];
  ode->order = step;
  ode->f = f;
  ode->_ws = ws;
  zNLEAssignNR( &ws->nle );
  return ode;
}

/* destroy an ODE solver. */
void zODEDestroyGear(zODE *ode)
{
  _zODE_Gear *ws;

  ws = (_zODE_Gear *)ode->_ws;
  zNLEDestroy( &ws->nle );
  zVecFreeAtOnce( 3, ws->v, ws->vp, ws->a );
  zVecRingFree( &ws->hist );
  ode->f = NULL;
  zFree( ode->_ws );
}

/* initialize the history of states as if the state had stayed at x. */
void zODEInitHistoryGear(zODE *ode, double t, zVec x, void *util)
{
  zVecRingFill( &((_zODE_Gear*)ode->_ws)->hist, x );
}

/* put a state into the history, where zODEUpdateGear() puts the current state at the head of the ring
 * and refers to the previous states after it. */
static void _zODEPushHistoryGear(zODE *ode, double t, zVec x, void *util)
{
  _zODE_Gear *ws;

  ws = (_zODE_Gear *)ode->_ws;
  zVecCopyNC( x, *zRingHead(&ws->hist) );
  zRingDecHead( &ws->hist );
}

/* create the history of states by an implicit one-step method of order not lower than k-1, which is
 * A-stable not to spoil the stability for stiff equations. a more accurate method than necessary is
 * chosen, since the cost for k-1 steps is negligible. */
zVec zODECreateHistoryGear(zODE *ode, double *t, zVec x, double dt, void *util)
{
  zODE starter;

  zODEInit( &starter );
  if( ode->order <= 4 )
    zODEAssign( &starter, Radau, ode->cat, ode->sub );
  else
    zODEAssign( &starter, Gauss, ode->cat, ode->sub );
  /* the deferred correction raises the order of accuracy of Gauss method to five for six steps */
  return ode->order >= 6 ?
    _zODECreateHistoryDC( ode, &starter, _zODEPushHistoryGear, t, x, dt, util ) :
    _zODECreateHistory( ode, &starter, _zODEPushHistoryGear, t, x, dt, util );
}

/* directly integrate variable by ODE based on Gear method. */
zVec zODEUpdateGear(zODE *ode, double t, zVec x, double dt, void *util)
{
  _zODE_Gear *ws;

  ws = (_zODE_Gear *)ode->_ws;
  ws->t = t;
  ws->dt = dt;
  ws->util = util;
  zVecCopyNC( x, *zRingHead(&ws->hist) );
  /* the residual is to be relatively small to the state. the default maximum iteration number is used. */
  zNLESolve( &ws->nle, x, ode, zTOL, zTOL*zVecElemAbsMax(x,NULL), 0, NULL );
  zRingDecHead( &ws->hist );
  return x;
}
