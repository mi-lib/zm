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
  zVec v, a;
  double b;
  zVecRing hist;
  void *util;
  zNLE nle;
} _zODE_Gear;

static zVec _zODE_Gear_Func(zVec x, zVec y, void *util);

/* (static)
 * _zODE_Gear_Func
 * - algebraic equation for Gear method.
 */
zVec _zODE_Gear_Func(zVec x, zVec y, void *util)
{
  zODE *ode = util;
  _zODE_Gear *ws = ode->_ws;
  register int i;

  zVecCopyNC( x, y );
  for( i=0; i<zVecSizeNC(ws->a); i++ )
    ode->cat( y, zVecElem(ws->a,i), *zRingElem(&ws->hist,i), y, util );
  ode->f( ws->t+ws->dt, x, ws->util, ws->v );
  ode->cat( y, -ws->b*ws->dt, ws->v, y, ws->util );
  return y;
}

/* zODEInit_Gear
 * - initialize ODE solver based on Gear method.
 */
zODE* zODEInit_Gear(zODE *ode, int dim, int step, zVec (* f)(double,zVec,void*,zVec))
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
      !( ws->a = zVecCloneArray( a[step-1], step ) ) ||
      !zVecRingAlloc( &ws->hist, dim, step ) ||
      !zNLECreate( &ws->nle, dim, dim, 0, _zODE_Gear_Func, NULL ) ){
    ZALLOCERROR();
    return NULL;
  }
  ws->b = b[step-1];
  ode->f = f;
  ode->_ws = ws;
  zNLEAssignNR( &ws->nle );
  return ode;
}

/* zODEDestroy_Gear
 * - destroy ODE solver.
 */
void zODEDestroy_Gear(zODE *ode)
{
  _zODE_Gear *ws = ode->_ws;

  zNLEDestroy( &ws->nle );
  zVecFreeAO( 2, ws->v, ws->a );
  zVecRingFree( &ws->hist );
  ode->f = NULL;
  zFree( ode->_ws );
}

/* zODEInitHist_Gear
 * - initialize history with a given vector.
 */
void zODEInitHist_Gear(zODE *ode, zVec x)
{
  zVecRingFill( &((_zODE_Gear*)ode->_ws)->hist, x );
}

/* zODEUpdate_Gear
 * - directly integrate variable by ODE based on Gear method.
 */
zVec zODEUpdate_Gear(zODE *ode, double t, zVec x, double dt, void *util)
{
  _zODE_Gear *ws = ode->_ws;

  ws->t = t;
  ws->dt = dt;
  ws->util = util;
  zVecCopyNC( x, *zRingHead(&ws->hist) );
  zNLESolve( &ws->nle, x, ode, zTOL, 0, NULL ); /* use default maximum iteration number */
  zRingDecHead( &ws->hist );
  return x;
}
