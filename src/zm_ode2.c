/* ZM - Z's Mathematics Toolbox
 * Copyright (C) 1998 Tomomichi Sugihara (Zhidao)
 *
 * zm_ode2 - ordinary differential equation quadrature:
 * second-order differential equation solver.
 */

#include <zm/zm_ode.h>

/* *** regular update - a simple extension of first-order ODE *** */

static zVec _zODE2RegularFunc(double t, zVec y, void *ode, zVec dy);

void zODE2AssignFunc(zODE2 *ode, zVec (* catf1)(zVec,double,zVec,zVec,void*), zVec (* catf2)(zVec,double,zVec,zVec,void*), zVec (* subf1)(zVec,zVec,zVec,void*), zVec (* subf2)(zVec,zVec,zVec,void*))
{
  ode->cat_dis = catf1 ? catf1 : _zODECatDefault;
  ode->cat_vel = catf2 ? catf2 : _zODECatDefault;
  ode->sub_dis = subf1 ? subf1 : _zODESubDefault;
  ode->sub_vel = subf2 ? subf2 : _zODESubDefault;
}

/* initialize */
zODE2 *zODE2InitRegular(zODE2 *ode, int dim, int step, zVec (* f)(double,zVec,zVec,void*,zVec))
{
  ode->f = f;
  if( !zODEInit( &ode->_ode, dim*2, step, _zODE2RegularFunc ) )
    return NULL;
  if( !( ode->_x = zVecAlloc(dim*2) ) ){
    ZALLOCERROR();
    return NULL;
  }
  return ode;
}

/* destroy */
void zODE2DestroyRegular(zODE2 *ode)
{
  ode->f = NULL;
  ode->util = NULL;
  zVecFree( ode->_x );
  zODEDestroy( &ode->_ode );
}

/* simply transform second-order ODE to first-order ODE. */
zVec _zODE2RegularFunc(double t, zVec y, void *ode, zVec dy)
{
  zVecStruct x, v, a;
  int size;

  size = zVecSizeNC(y) / 2;
  zVecSetSize( &x, size ); zVecBuf(&x) = zVecBuf(y);
  zVecSetSize( &v, size ); zVecBuf(&v) = zVecBuf(y)  + size;
  zVecSetSize( &a, size ); zVecBuf(&a) = zVecBuf(dy) + size;
  zRawVecCopy( zVecBuf(&v), zVecBuf(dy), size );
  ((zODE2 *)ode)->f( t, &x, &v, ((zODE2 *)ode)->util, &a );
  return dy;
}

/* concatenate state variable vector */
zVec _zODE2Cat(zVec y, double dt, zVec dy, zVec yn, void *ode)
{
  zVecStruct x, xn, v, vn, vf, af;
  int size;

  size = zVecSizeNC(y) / 2;
  zVecSetSize( &x,  size ); zVecBuf(&x)  = zVecBuf(y);
  zVecSetSize( &v,  size ); zVecBuf(&v)  = zVecBuf(y)  + size;
  zVecSetSize( &xn, size ); zVecBuf(&xn) = zVecBuf(yn);
  zVecSetSize( &vn, size ); zVecBuf(&vn) = zVecBuf(yn) + size;
  zVecSetSize( &vf, size ); zVecBuf(&vf) = zVecBuf(dy);
  zVecSetSize( &af, size ); zVecBuf(&af) = zVecBuf(dy) + size;
  ((zODE2 *)ode)->cat_dis( &x, dt, &vf, &xn, ((zODE2 *)ode)->util );
  ((zODE2 *)ode)->cat_vel( &v, dt, &af, &vn, ((zODE2 *)ode)->util );
  return yn;
}

/* subtract a state variable vector from another */
zVec _zODE2Sub(zVec y1, zVec y2, zVec dy, void *ode)
{
  zVecStruct x1, x2, dx, v1, v2, dv;
  int size;

  size = zVecSizeNC(y1) / 2;
  zVecSetSize( &x1, size ); zVecBuf(&x1) = zVecBuf(y1);
  zVecSetSize( &v1, size ); zVecBuf(&v1) = zVecBuf(y1)  + size;
  zVecSetSize( &x2, size ); zVecBuf(&x2) = zVecBuf(y2);
  zVecSetSize( &v2, size ); zVecBuf(&v2) = zVecBuf(y2) + size;
  zVecSetSize( &dx, size ); zVecBuf(&dx) = zVecBuf(dy);
  zVecSetSize( &dv, size ); zVecBuf(&dv) = zVecBuf(dy) + size;
  ((zODE2 *)ode)->sub_dis( &x1, &x2, &dx, ((zODE2 *)ode)->util );
  ((zODE2 *)ode)->sub_vel( &v1, &v2, &dv, ((zODE2 *)ode)->util );
  return dy;
}

/* update state destructively. */
void zODE2UpdateRegular(zODE2 *ode, double t, zVec x, zVec v, double dt, void *util)
{
  zRawVecCopy( zVecBuf(x), zVecBuf(ode->_x), zVecSizeNC(x) );
  zRawVecCopy( zVecBuf(v), zVecBuf(ode->_x)+zVecSizeNC(v), zVecSizeNC(v) );
  ode->util = util;
  zODEUpdate( &ode->_ode, t, ode->_x, dt, ode );
  zRawVecCopy( zVecBuf(ode->_x), zVecBuf(x), zVecSizeNC(x) );
  zRawVecCopy( zVecBuf(ode->_x)+zVecSizeNC(v), zVecBuf(v), zVecSizeNC(v) );
}

/* *** simplest symplectic method *** */

/* initialize */
zODE2 *zODE2InitSympl(zODE2 *ode, int dim, int step, zVec (* f)(double,zVec,zVec,void*,zVec))
{
  ode->f = f;
  if( !( ode->_a = zVecAlloc(dim) ) ){
    ZALLOCERROR();
    return NULL;
  }
  return ode;
}

/* destroy */
void zODE2DestroySympl(zODE2 *ode)
{
  ode->f = NULL;
  zVecFree( ode->_a );
}

/* update state destructively. */
void zODE2UpdateSympl(zODE2 *ode, double t, zVec x, zVec v, double dt, void *util)
{
  ode->f( t, x, v, util, ode->_a );
  ode->cat_vel( v, dt, ode->_a, v, util );
  ode->cat_dis( x, dt, v, x, util );
}

/* *** leapfrog method *** */

/* initialize */
zODE2 *zODE2InitLeapfrog(zODE2 *ode, int dim, int step, zVec (* f)(double,zVec,zVec,void*,zVec))
{
  ode->f = f;
  ode->_x = zVecAlloc(dim);
  ode->_v = zVecAlloc(dim);
  ode->_a = zVecAlloc(dim);
  if( !ode->_x || !ode->_v || !ode->_a ){
    ZALLOCERROR();
    zODE2DestroyLeapfrog( ode );
    return NULL;
  }
  return ode;
}

/* initialize startup history */
void zODE2InitHistLeapfrog(zODE2 *ode, zVec x, zVec v, double dt)
{
  zVecCopy( v, ode->_v ); /* v_-0.5 = v_0 */
}

/* destroy */
void zODE2DestroyLeapfrog(zODE2 *ode)
{
  ode->f = NULL;
  zVecFreeAO( 3, ode->_x, ode->_v, ode->_a );
}

/* update state destructively */
void zODE2UpdateLeapfrog(zODE2 *ode, double t, zVec x, zVec v, double dt, void *util)
{
  ode->f( t, x, v, util, ode->_a );
  /* v_(i+1/2) = v_(i-1/2) + dt * a_(i) */
  ode->cat_vel( ode->_v, dt, ode->_a, ode->_v, util );
  /* x_(i+1) = x_(i) + dt * v_(i+1/2) */
  ode->cat_dis( x, dt, ode->_v, x, util );
  /* v_(i+1) = v_(i) + 0.5dt * a_(i+1/2) */
  ode->f( t, x, ode->_v, util, ode->_a );
  ode->cat_vel( ode->_v, 0.5*dt, ode->_a, v, util );
}
