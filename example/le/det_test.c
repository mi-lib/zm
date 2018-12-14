#include <zm/zm_le.h>

#define TEST 1

int main(void)
{
#if TEST == 1
  double marray[] = { /* 66 */
    3,-1, 2, 4,
    2, 1, 1, 3,
   -2, 0, 3,-1,
    0,-2, 1, 3 };
  unsigned s = 4;
#elif TEST == 2
  double marray[] = { /* -4 */
    1, 3,
    2, 2, };
  unsigned s = 2;
#elif TEST == 3
  double marray[] = { /* 6 */
    2, 0, 1,-1,
   -4, 2, 3, 1,
    1, 0,-1, 2,
   -2, 1, 0, 2 };
  unsigned s = 4;
#elif TEST == 4
  double marray[] = { /* 2 */
    2, 1,-1,-1,
   -2,-3, 1, 2,
   -2,-1, 3, 1,
   -1,-4, 1, 2 };
  unsigned s = 4;
#elif TEST == 5
  double marray[] = { /* 74 */
   -3,-2, 2,-5,
   -2, 3,-2, 3,
    4,-1, 2,-2,
   -3, 2,-4, 3 };
  unsigned s = 4;
#elif TEST == 6
  double marray[] = { /* 0 */
    3, 1, 1,-1,
    7, 2,-4, 1,
   19, 6, 0,-3,
   22, 7,-1,-3 };
  unsigned s = 4;
#elif TEST == 7
  double marray[] = { /* 104 */
    1, 3, 2, 5,
    3, 4, 9,11,
    2, 2, 1, 9,
    1, 4, 3, 9 };
  unsigned s = 4;
#elif TEST == 8
  double marray[] = { /* 0 */
    0, 0, 0, 0,
    0, 0, 0, 0,
    0, 0, 0, 0,
    0, 0, 0, 0 };
  unsigned s = 4;
#else
  double marray[] = { /* 0.0001 */
    -3.9,-1, 1, 1,
    18,5.1,-3,-4,
    2, 1, 1.1, 0,
    12, 4, 0,-1.9 };
  unsigned s = 4;
#endif
  zMat m, l, u;
  zIndex idx;

  m = zMatCloneArray( marray, s, s );
  printf( "matrix: " ); zMatWrite( m );
  printf( "determinant = %f\n", zMatDet( m ) );

  /* confirmation */
  l = zMatAllocSqr( s );
  u = zMatAllocSqr( s );
  idx = zIndexCreate( s );
  zLUDecomp( m, l, u, idx );
  printf( "L: " ); zMatWrite( l );
  zMatFree( l );
  zMatFree( u );
  zIndexFree( idx );

  zMatFree( m );
  return 0;
}
