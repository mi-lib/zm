/* ZM - Z's Mathematics Toolbox
 * Copyright (C) 1998 Tomomichi Sugihara (Zhidao)
 *
 * zm_ode - ordinary differential equation quadrature.
 */

#include <zm/zm_ode.h>

/* _zODECatDefault
 * - concatenate deviation in finite time step to the current
 *   variable vector.
 */
zVec _zODECatDefault(zVec x, double dt, zVec dx, zVec xn, void *dummy)
{
  return zVecCatNC( x, dt, dx, xn );
}

/* _zODESubDefault
 * - subtract a variable vector from another.
 */
zVec _zODESubDefault(zVec x1, zVec x2, zVec dx, void *dummy)
{
  return zVecSubNC( x1, x2, dx );
}

void zODEAssignFunc(zODE *ode, zVec (* catf)(zVec,double,zVec,zVec,void*), zVec (* subf)(zVec,zVec,zVec,void*))
{
  ode->cat = catf ? catf : _zODECatDefault;
  ode->sub = subf ? subf : _zODESubDefault;
}
