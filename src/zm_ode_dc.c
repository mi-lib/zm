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
  return zODECreate( ode, dim, step, f );
}

/* destroy an ODE solver enabling deferred correction. */
void zODEDestroyDC(zODE *ode)
{
  zVecFree( ode->_x1 );
  zVecFree( ode->_x2 );
  zODEDestroy( ode );
}

/* update state vector by solving ODE enabling deferred correction. */
zVec zODEUpdateDC(zODE *ode, double t, zVec x, double dt, void *util)
{
  double hdt;

  if( dt == 0 ) return NULL; /* see zODEUpdate */
  /* full-step update */
  zVecCopyNC( x, ode->_x1 );
  zODEUpdate( ode, t, ode->_x1, dt, util );
  /* half-step-twice update */
  zVecCopyNC( x, ode->_x2 );
  hdt = 0.5 * dt;
  zODEUpdate( ode, t,     ode->_x2, hdt, util );
  zODEUpdate( ode, t+hdt, ode->_x2, hdt, util );
  /* evaluation (the step is not divided any more if it is already small enough) */
  if( dt > ZM_ODE_DT_TOL && !zODEStepIsTol( x, ode->_x1, ode->_x2 ) ){
    zODEUpdateDC( ode, t,     x, hdt, util );
    zODEUpdateDC( ode, t+hdt, x, hdt, util );
  } else{ /* Richardson's extrapolation x = x2 + ( x2 - x1 ) / ( 2^p - 1 ) for a method of order p */
    ode->sub( ode->_x2, ode->_x1, ode->_x1, util );
    ode->cat( ode->_x2, 1.0 / ( ( 1 << ( ode->order > 0 ? ode->order : 1 ) ) - 1 ), ode->_x1, x, util );
  }
  return x;
}
