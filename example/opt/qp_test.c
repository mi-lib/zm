#include <zm/zm_opt.h>

int main(void)
{
  double qarray[] = {
    1, 2,-4, 5
  };
  double aarray[] = {
    1, 2
  };
  zMat q, a;
  zVec z, b, x;
  double cost;

  q = zMatAllocSqr( 2 );
  a = zMatAlloc( 1, 2 );
  z = zVecAlloc( 2 );
  b = zVecAlloc( 1 );
  x = zVecAlloc( 2 );

  zMatCopyArray( qarray, 2, 2, q );
  zMatCopyArray( aarray, 1, 2, a );
  zVecSetElem( b, 0, 3 );
  zMatWrite( q );
  zMatWrite( a );
  zVecWrite( b );
  zQPSolve( q, z, a, b, x, &cost );
  zVecWrite( x );

  zMatFree( q );
  zMatFree( a );
  zVecFree( b );
  zVecFree( x );
  return 0;
}
