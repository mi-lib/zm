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
  zMat a, ai, b;

  a = zMatCloneArray( marray, r, c );
  ai = zMatAlloc( c, r );
  zMPInv( a, ai );
  printf( "A: " );
  zMatWrite( a );
  printf( "A^+: " );
  zMatWrite( ai );

  b = zMatAllocSqr( r );
  zMulMatMat( a, ai, b );
  zMulMatMatDRC( b, a );
  printf( "AA^+A: " );
  zMatWrite( a );
  zMatFree( b );

  b = zMatAllocSqr( c );
  zMulMatMat( ai, a, b );
  zMulMatMatDRC( b, ai );
  printf( "A^+AA^+: " );
  zMatWrite( ai );
  zMatFree( b );

  zMatFree( a );
  zMatFree( ai );
  return 0;
}
