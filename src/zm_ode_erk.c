/* ZM - Z's Mathematics Toolbox
 * Copyright (C) 1998 Tomomichi Sugihara (Zhidao)
 *
 * zm_ode_erk - ordinary differential equation quadrature:
 * embedded Runge-Kutta method.
 */

#include <zm/zm_ode.h>

#define ZODE_ERK_TOL  ( 1.0e-6 )
#define ZODE_ERK_DT_TOL  ( 1.0e-8 )
typedef struct{
  int stepsize;
  zVec _xc, _xf; /* course and fine estimations */
  zVec *k;
} _zODE_ERK;

/* _zODEInit_ERK
 * - initialize ODE solver based on embedded Runge-Kutta method.
 */
zODE* _zODEInit_ERK(zODE *ode, int dim, int stepsize, zVec (* f)(double,zVec,void*,zVec))
{
  register int i;
  _zODE_ERK *ws;
  bool check = true;

  if( !( ws = zAlloc( _zODE_ERK, 1 ) ) ){
    ZALLOCERROR();
    return NULL;
  }
  if( !( ws->_xc = zVecAlloc( dim ) ) ||
      !( ws->_xf = zVecAlloc( dim ) ) ||
      !( ws->k = zAlloc( zVec, stepsize ) ) ) check = false;
  for( i=0; i<stepsize; i++ )
    if( !( ws->k[i] = zVecAlloc( dim ) ) ) check = false;
  ws->stepsize = stepsize;
  if( !check ){
    ZALLOCERROR();
    return NULL;
  }
  ode->f = f;
  ode->_ws = ws;
  return ode;
}

/* _zODEDestroy_ERK
 * - destroy ODE solver.
 */
void _zODEDestroy_ERK(zODE *ode)
{
  register int i;
  _zODE_ERK *ws;

  ws = ode->_ws;
  zVecFreeAO( 2, ws->_xc, ws->_xf );
  for( i=0; i<ws->stepsize; i++ )
    zVecFree( ws->k[i] );
  zFree( ws->k );
  zFree( ode->_ws );
  ode->f = NULL;
}

/* _zODEUpdate_ERK
 * - directly integrate variable by ODE based on embedded Runge-Kutta method.
 */
zVec _zODEUpdate_ERK(zODE *ode, double t, zVec x, double dt, double a[], double bc[], double bf[], double c[], void *util)
{
  register int i, j, l;
  _zODE_ERK *ws;

  ws = ode->_ws;
  ode->f( t, x, util, ws->k[0] );
  for( l=0, i=1; i<ws->stepsize; i++ ){
    zVecCopyNC( x, ws->_xc );
    for( j=0; j<i; j++ )
      ode->cat( ws->_xc, a[l++]*dt, ws->k[j], ws->_xc, util );
    ode->f( t+c[i-1]*dt, ws->_xc, util, ws->k[i] );
  }
  zVecCopyNC( x, ws->_xc );
  zVecCopyNC( x, ws->_xf );
  for( i=0; i<ws->stepsize; i++ ){
    ode->cat( ws->_xc, bc[i]*dt, ws->k[i], ws->_xc, util );
    ode->cat( ws->_xf, bf[i]*dt, ws->k[i], ws->_xf, util );
  }
  if( dt > ZODE_ERK_DT_TOL && zVecDist( ws->_xc, ws->_xf ) > ZODE_ERK_TOL ){
    dt *= 0.5;
    _zODEUpdate_ERK( ode, t, x, dt, a, bc, bf, c, util );
    _zODEUpdate_ERK( ode, t+dt, x, dt, a, bc, bf, c, util );
  } else
    zVecCopyNC( ws->_xf, x );
  return x;
}

/* Runge-Kutta-Fehlberg method */

zODE* zODEInit_RKF45(zODE *ode, int dim, int dummy, zVec (* f)(double,zVec,void*,zVec))
{
  return _zODEInit_ERK( ode, dim, 6, f );
}

void zODEDestroy_RKF45(zODE *ode)
{
  _zODEDestroy_ERK( ode );
}

zVec zODEUpdate_RKF45(zODE *ode, double t, zVec x, double dt, void *util)
{
  static double a[] = {
    1.0/4,
    3.0/32, 9.0/32,
    1932.0/2197, -7200.0/2197, 7296.0/2197,
    439.0/216, -8.0, 3680.0/513, -845.0/4104,
    -8.0/27, 2.0, -3544.0/2565, 1859.0/4104, -11.0/40,
  };
  static double b4[] = { 25.0/216, 0, 1408.0/2565, 2197.0/4104, -1.0/5, 0 };
  static double b5[] = { 16.0/135, 0, 6656.0/12825, 28561.0/56430, -9.0/50, 2.0/55 };
  static double c[] = { 1.0/4, 3.0/8, 12.0/13, 1.0, 1.0/2 };

  return _zODEUpdate_ERK( ode, t, x, dt, a, b4, b5, c, util );
}

/* Cash-Karp method */

zODE* zODEInit_CK45(zODE *ode, int dim, int dummy, zVec (* f)(double,zVec,void*,zVec))
{
  return _zODEInit_ERK( ode, dim, 6, f );
}

void zODEDestroy_CK45(zODE *ode)
{
  _zODEDestroy_ERK( ode );
}

zVec zODEUpdate_CK45(zODE *ode, double t, zVec x, double dt, void *util)
{
  static double a[] = {
    1.0/5,
    3.0/40, 9.0/40,
    3.0/10,-9.0/10, 6.0/5,
  -11.0/54, 5.0/2, -70.0/27, 35.0/27,
 1631.0/55296, 175.0/512, 575.0/13824, 44275.0/110592, 253.0/4096,
  };
  static double b4[] = { 37.0/378, 0.0, 250.0/621, 125.0/594, 0.0, 512.0/1771 };
  static double b5[] = { 2825.0/27648, 0.0, 18575.0/48384, 13525.0/55296, 277.0/14336, 1.0/4 };
  static double c[] = { 1.0/5, 3.0/10, 3.0/5, 1.0, 7.0/8.0 };

  return _zODEUpdate_ERK( ode, t, x, dt, a, b4, b5, c, util );
}

/* Dormand-Prince method */

zODE* zODEInit_DP45(zODE *ode, int dim, int dummy, zVec (* f)(double,zVec,void*,zVec))
{
  return _zODEInit_ERK( ode, dim, 7, f );
}

void zODEDestroy_DP45(zODE *ode)
{
  _zODEDestroy_ERK( ode );
}

zVec zODEUpdate_DP45(zODE *ode, double t, zVec x, double dt, void *util)
{
  static double a[] = {
    1.0/5,
    3.0/40, 9.0/40,
   44.0/45, -56.0/15, 32.0/9,
19372.0/6561,-25360.0/2187, 64448.0/6561, -212.0/729,
 9017.0/3168, -355.0/33, 46732.0/5247, 49.0/176, -5103.0/18656,
   35.0/384, 0.0, 500.0/1113, 125.0/192, -2187.0/6784, 11.0/84,
  };
  static double b4[] = {
    5179.0/57600, 0.0, 7571.0/16695, 393.0/640,-92097.0/339200, 187.0/2100, 1.0/40 };
  static double b5[] = {
    35.0/384, 0.0, 500.0/1113, 125.0/192, -2187.0/6784, 11.0/84, 0 };
  static double c[] = { 1.0/5, 3.0/10, 4.0/5, 8.0/9, 1.0, 1.0 };

  return _zODEUpdate_ERK( ode, t, x, dt, a, b4, b5, c, util );
}
