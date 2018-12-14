/* ZM - Z's Mathematics Toolbox
 * Copyright (C) 1998 Tomomichi Sugihara (Zhidao)
 *
 * zm_ode_dc - ordinary differential equation quadrature:
 * deferred correction with Richardson's extrapolation.
 */

#include <zm/zm_ode.h>

/* zODEInitDC
 * - initialize ODE solver enabling deferred correction.
 */
zODE *zODEInitDC(zODE *ode, int dim, int step, zVec (* f)(double,zVec,void*,zVec))
{
  ode->_x1 = zVecAlloc( dim );
  ode->_x2 = zVecAlloc( dim );
  if( !ode->_x1 || !ode->_x2 ){
    ZALLOCERROR();
    return NULL;
  }
  return zODEInit( ode, dim, step, f );
}

/* zODEDestroyDC
 * - destroy ODE solver enabling deferred correction.
 */
void zODEDestroyDC(zODE *ode)
{
  zVecFree( ode->_x1 );
  zVecFree( ode->_x2 );
  zODEDestroy( ode );
}

#define ZODE_DC_TOL ( 1.0e-6 )
/* zODEUpdateDC
 * - update state vector by solving ODE enabling deferred correction.
 */
zVec zODEUpdateDC(zODE *ode, double t, zVec x, double dt, void *util)
{
  /* full-step update */
  zVecCopyNC( x, ode->_x1 );
  zODEUpdate( ode, t, ode->_x1, dt, util );
  /* half-step-twice update */
  zVecCopyNC( x, ode->_x2 );
  dt *= 0.5;
  zODEUpdate( ode, t,    ode->_x2, dt, util );
  zODEUpdate( ode, t+dt, ode->_x2, dt, util );
  /* evaluation */
  if( zVecDist( ode->_x1, ode->_x2 ) > ZODE_DC_TOL ){
    zODEUpdateDC( ode, t,    x, dt, util );
    zODEUpdateDC( ode, t+dt, x, dt, util );
  } else{ /* Richardson's extrapolation */
    zVecMul( ode->_x2, 2, x );
    ode->sub( x, ode->_x1, x, util );
  }
  return x;
}
