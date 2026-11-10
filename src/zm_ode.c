/* ZM - Z's Mathematics Toolbox
 * Copyright (C) 1998 Tomomichi Sugihara (Zhidao)
 *
 * zm_ode - ordinary differential equation quadrature.
 */

#include <zm/zm_ode.h>

/* concatenate deviation in finite time step to the current variable vector. */
zVec _zODECatDefault(zVec x, double dt, zVec dx, zVec xn, void *dummy)
{
  return zVecCatNC( x, dt, dx, xn );
}

/* subtract a variable vector from another. */
zVec _zODESubDefault(zVec x1, zVec x2, zVec dx, void *dummy)
{
  return zVecSubNC( x1, x2, dx );
}

/* assign concatenation and subtraction functions to an ODE solver. */
void zODEAssignFunc(zODE *ode, zVec (* catf)(zVec,double,zVec,zVec,void*), zVec (* subf)(zVec,zVec,zVec,void*))
{
  ode->cat = catf ? catf : _zODECatDefault;
  ode->sub = subf ? subf : _zODESubDefault;
}

/* ratio of the estimated error between two estimations of the state vector to the tolerance. */
double zODEStepErr(const zVec x0, const zVec x1, const zVec x2, double r)
{
  int i;
  double e, err = 0;

  for( i=0; i<zVecSizeNC(x0); i++ ){
    e = r * fabs( zVecElemNC(x1,i) - zVecElemNC(x2,i) ) / ( ZM_ODE_TOL * ( 1 + zMax( fabs( zVecElemNC(x0,i) ), fabs( zVecElemNC(x2,i) ) ) ) );
    if( !( e <= err ) ) err = e; /* NaN is propagated */
  }
  return err;
}
