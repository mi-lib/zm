#include <zm/zm_le.h>

#define TEST 0

int main(void)
{
#if TEST == 1
  /* regular case */
  double marray[] = {
    1, 2,
   -2,-4 };
  unsigned r = 2, c = 2;
#elif TEST == 2
  /* regular case (row excess) */
  double marray[] = {
    2, 1, 1,
    1, 1, 0 };
  unsigned r = 2, c = 3;
#elif TEST == 3
  /* regular case (row excess) */
  double marray[] = {
    4, 5,-2, 4,
    5, 3, 4, 3,
   -3, 1, 2,-4,
    8, 2, 1, 1,
   -1,-3, 2, 1,
    1, 3, 0,-2 };
  unsigned r = 6, c = 4;
#else
  /* singular case (column excess) */
  double marray[] = {
    4, 5,-2, 4, 5, 3,
    4, 3,-3, 1, 2,-4,
    8, 2, 1, 1,-1,-3,
    4, 5,-2, 4, 5, 3,
    0, 2, 1, 3, 3, 7 };
  unsigned r = 5, c = 6;
#endif
  zMat a, ai, b, d, e;

  a = zMatCloneArray( marray, r, c );
  ai = zMatAlloc( c, r );
  zMPInv( a, ai );
  printf( "A: " );
  zMatPrint( a );
  printf( "A^+: " );
  zMatPrint( ai );

  b = zMatAllocSqr( r );
  d = zMatAlloc( r, c );
  e = zMatAlloc( r, c );
  zMulMatMat( a, ai, b );
  zMulMatMat( b, a, d );
  zMatSub( a, d, e );
  printf( "|| A - AA^+A || = %.10f\n", zMatNorm(e) );
  zMatFree( b );
  zMatFree( d );
  zMatFree( e );

  b = zMatAllocSqr( c );
  d = zMatAlloc( c, r );
  e = zMatAlloc( c, r );
  zMulMatMat( ai, a, b );
  zMulMatMat( b, ai, d );
  zMatSub( ai, d, e );
  printf( "|| A^+ - A^+AA^+ || = %.10f\n", zMatNorm(e) );
  zMatFree( b );
  zMatFree( d );
  zMatFree( e );

  zMatFree( a );
  zMatFree( ai );
  return 0;
}
