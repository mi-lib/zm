/* ZM - Z's Mathematics Toolbox
 * Copyright (C) 1998 Tomomichi Sugihara (Zhidao)
 *
 * zm_ode_dc - ordinary differential equation quadrature:
 * deferred correction with Richardson's extrapolation.
 */

#include <zm/zm_ode.h>

/* create an ODE solver enabling deferred correction. */
zODE *zODECreateDC(zODE *ode, int dim, int step, zVec (* f)(double,zVec,void*,zVec))
{
  ode->_x1 = zVecAlloc( dim );
  ode->_x2 = zVecAlloc( dim );
  if( !ode->_x1 || !ode->_x2 ){
    ZALLOCERROR();
    return NULL;
  }
  if( !zODECreate( ode, dim, step, f ) ) return NULL;
  ode->_h = 0;
  return ode;
}

/* destroy an ODE solver enabling deferred correction. */
void zODEDestroyDC(zODE *ode)
{
  zVecFree( ode->_x1 );
  zVecFree( ode->_x2 );
  zODEDestroy( ode );
}

/* order of accuracy of the method for the deferred correction (one if unknown) */
#define _zODEOrderDC(ode) ( (ode)->order > 0 ? (ode)->order : 1 )

/* try a time step with deferred correction, namely, compute the full-step and the half-step-twice
 * estimations of the state vector, and return the ratio of the estimated error to the tolerance. */
static double _zODETrialDC(zODE *ode, double t, zVec x, double dt, void *util, void *dummy)
{
  double hdt;

  /* full-step update */
  zVecCopyNC( x, ode->_x1 );
  zODEUpdate( ode, t, ode->_x1, dt, util );
  /* half-step-twice update */
  zVecCopyNC( x, ode->_x2 );
  hdt = 0.5 * dt;
  zODEUpdate( ode, t,     ode->_x2, hdt, util );
  zODEUpdate( ode, t+hdt, ode->_x2, hdt, util );
  /* for a method of order p, the error of x2 is estimated as ( x2 - x1 ) / ( 2^p - 1 ). */
  return zODEStepErr( x, ode->_x1, ode->_x2, 1.0 / ( ( 1 << _zODEOrderDC(ode) ) - 1 ) );
}

/* accept the last trial with Richardson's extrapolation x = x2 + ( x2 - x1 ) / ( 2^p - 1 ). */
static void _zODEAcceptDC(zODE *ode, zVec x, void *util, void *dummy)
{
  ode->sub( ode->_x2, ode->_x1, ode->_x1, util );
  ode->cat( ode->_x2, 1.0 / ( ( 1 << _zODEOrderDC(ode) ) - 1 ), ode->_x1, x, util );
}

/* update state vector by solving ODE enabling deferred correction. */
zVec zODEUpdateDC(zODE *ode, double t, zVec x, double dt, void *util)
{
  if( dt == 0 ) return NULL; /* see zODEUpdate */
  /* the error of x2 is in the order of dt^(p+1) */
  return zODEUpdateAdaptive( ode, t, x, dt, _zODEOrderDC(ode), _zODETrialDC, _zODEAcceptDC, util, NULL );
}
