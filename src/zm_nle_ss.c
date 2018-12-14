/* ZM - Z's Mathematics Toolbox
 * Copyright (C) 1998 Tomomichi Sugihara (Zhidao)
 *
 * zm_nle_ss - nonlinear equation: successive substitution method.
 */

#include <zm/zm_nle.h>

/* successive substitution method to solve equation x=f(x).
 */
zVec zSSSolve(zVec (* f)(zVec,zVec,void*), zVec x, void *util, int iter)
{
  register int i;
  zVec y, e;

  y = zVecAlloc( zVecSizeNC(x) );
  e = zVecAlloc( zVecSizeNC(x) );
  if( !y || !e ){
    x = NULL;
    goto TERMINATE;
  }
  ZITERINIT( iter );
  for( i=0; i<iter; i++ ){
    f( x, y, util );
    zVecSubNC( y, x, e );
    if( zVecIsTiny( e ) ) goto TERMINATE;
    zVecCopyNC( y, x );
  }
  ZITERWARN( iter );
 TERMINATE:
  zVecFree( y );
  zVecFree( e );
  return x;
}
