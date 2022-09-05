/* ZM - Z's Mathematics Toolbox
 * Copyright (C) 1998 Tomomichi Sugihara (Zhidao)
 *
 * zm_ode_euler - ordinary differential equation quadrature:
 * Euler method.
 */

#include <zm/zm_ode.h>

/* initialize ODE solver based on Euler method. */
zODE* zODEInit_Euler(zODE *ode, int dim, int dummy, zVec (* f)(double,zVec,void*,zVec))
{
  if( !( ode->_ws = zVecAlloc( dim ) ) ){
    ZALLOCERROR();
    return NULL;
  }
  ode->f = f;
  return ode;
}

/* destroy ODE solver. */
void zODEDestroy_Euler(zODE *ode)
{
  zVecFree( (zVec)ode->_ws );
  ode->f = NULL;
}

/* directly integrate variable by ODE based on Euler method. */
zVec zODEUpdate_Euler(zODE *ode, double t, zVec x, double dt, void *util)
{
  ode->f( t, x, util, (zVec)ode->_ws );
  return ode->cat( x, dt, (zVec)ode->_ws, x, util );
}
