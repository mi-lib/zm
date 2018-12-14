/* ZM - Z's Mathematics Toolbox
 * Copyright (C) 1998 Tomomichi Sugihara (Zhidao)
 *
 * zm_ode_rkg - ordinary differential equation quadrature:
 * Runge-Kutta-Gill method.
 */

#include <zm/zm_ode.h>

/* Runge-Kutta-Gill method */
typedef struct{
  zVec u, v; /* working space */
} _zODE_RKG;

/* zODEInit_RKG
 * - initialize ODE solver based on Runge-Kutta-Gill method.
 */
zODE* zODEInit_RKG(zODE *ode, int dim, int dummy, zVec (* f)(double,zVec,void*,zVec))
{
  _zODE_RKG *ws;

  if( !( ws = zAlloc( _zODE_RKG, 1 ) ) ||
      !( ws->u = zVecAlloc( dim ) ) /* workspace 1 */ ||
      !( ws->v = zVecAlloc( dim ) ) /* workspace 2 */ ){
    ZALLOCERROR();
    return NULL;
  }
  ode->f = f;
  ode->_ws = ws;
  return ode;
}

/* zODEDestroy_RKG
 * - destroy ODE solver.
 */
void zODEDestroy_RKG(zODE *ode)
{
  _zODE_RKG *ws;

  ws = ode->_ws;
  zVecFreeAO( 2, ws->u, ws->v );
  zFree( ode->_ws );
  ode->f = NULL;
}

/* zODEUpdate_RKG
 * - directly integrate variable by ODE based on Runge-Kutta-Gill method.
 */
zVec zODEUpdate_RKG(zODE *ode, double t, zVec x, double dt, void *util)
{
  _zODE_RKG *ws;
  const double r2 = sqrt(2);
  const double c1 = 1 - 0.5*r2;
  const double c2 = 1 + 0.5*r2;
  double dt1, dt2, dt3;

  ws = ode->_ws;
  dt1 = 0.5 * dt;
  dt2 = c1 * dt;
  dt3 = c2 * dt;
  /* first step */
  ode->f( t, x, util, ws->u );
  ode->cat( x, dt1, ws->u, x, util );
  zVecCopyNC( ws->u, ws->v );
  /* second step */
  ode->f( t+dt1, x, util, ws->u );
  ode->cat( x, dt2, ws->u, x, util );
  ode->cat( x,-dt2, ws->v, x, util );
  zVecMulNCDRC( ws->v, 0.5*(3*r2-4) );
  zVecCatNCDRC( ws->v, 2*c1, ws->u );
  /* third step */
  ode->f( t+dt1, x, util, ws->u );
  ode->cat( x, dt3, ws->u, x, util );
  ode->cat( x,-dt3, ws->v, x, util );
  zVecMulNCDRC( ws->v,-0.5*(3*r2+4) );
  zVecCatNCDRC( ws->v, 2*c2, ws->u );
  /* fourth step */
  ode->f( t+dt, x, util, ws->u );
  ode->cat( x, dt/6, ws->u, x, util );
  ode->cat( x,-dt/3, ws->v, x, util );
  return x;
}
