/* linear equation solution by
 * 1. LU decomposition method
 * 2. Gauss's method
 */
#include <zm/zm_le.h>

int main(void)
{
#if TEST == 1
  double marray[] = {
    4, 0,-2, 0,
    0, 3, 0, 1,
    8, 0, 1, 1,
    0, 1, 1,-3 };
  double varray[] = {
    2, 4, 10, -1
  };
  unsigned s = 4;
#else
  double marray[] = {
    4, 5,-2, 4,
    4, 3,-3, 1,
    8, 2, 1, 1,
    2, 1, 1, 3 };
  double varray[] = {
    11, 5, 12, 7
  };
  unsigned s = 4;
#endif
  zMat a, l, u;
  zVec b, x;
  zIndex index;

  a = zMatCloneArray( marray, s, s );
  l = zMatAllocSqr( s );
  u = zMatAllocSqr( s );
  b = zVecCloneArray( varray, s );
  x = zVecAlloc( s );
  index = zIndexCreate( s );

  zMatWrite( a );
  zVecWrite( b );

  printf( "LU decomposition method.\n" );
  zLUDecomp( a, l, u, index );
  zMatWrite( l );
  zMatWrite( u );
  zLESolveLU( a, b, x );
  zVecWrite( x );

  printf( "Gauss's elimination method.\n" );
  zLESolveGauss( a, b, x );
  zVecWrite( x );

  zMatFree( a );
  zMatFree( l );
  zMatFree( u );
  zVecFree( b );
  zVecFree( x );
  zIndexFree( index );
  return 0;
}
