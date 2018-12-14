#include <zm/zm_le.h>

#define TEST 4

int main(void)
{
#if TEST == 1
  /* regular case (column excess) */
  double marray[] = {
    4, 5,-2, 4, 5, 3,
    4, 3,-3, 1, 2,-4,
    8, 2, 1, 1,-1,-3,
    2, 1, 1, 3, 0,-2 };
  unsigned r = 4, c = 6;
#elif TEST == 2
  /* regular case (row excess) */
  double marray[] = {
    4, 5,-2, 4,
    5, 3, 4, 3,
   -3, 1, 2,-4,
    8, 2, 1, 1,
   -1,-3, 2, 1,
    1, 3, 0,-2 };
  unsigned r = 6, c = 4;
#elif TEST == 3
  /* singular case (column excess) */
  double marray[] = {
    4, 5,-2, 4, 5, 3,
    4, 3,-3, 1, 2,-4,
    8, 2, 1, 1,-1,-3,
    4, 5,-2, 4, 5, 3,
    0, 2, 1, 3, 3, 7 };
  unsigned r = 5, c = 6;
#elif TEST == 4
  /* singular case (row excess) */
  double marray[] = {
    4, 5,-2, 4,
    5, 3, 4, 3,
    4, 5,-2, 4,
   -3, 1, 2,-4,
    1,-2, 6,-1,
    5, 3, 4, 3 };
  unsigned r = 6, c = 4;
#else /* TEST == 0 */
  /* regular case */
  double marray[] = {
    4, 5,-2, 4,
    4, 3,-3, 1,
    8, 2, 1, 1,
    2, 1, 1, 3 };
  unsigned r = 4, c = 4;
#endif

  zMat a, l, u;
  zIndex index;
  int rank;

  a = zMatCloneArray( marray, r, c );
  rank = zLUDecompAlloc( a, &l, &u, &index );
  printf( "L: " ); zMatWrite( l );
  printf( "U: " ); zMatWrite( u );
  printf( "A: " ); zMatWrite( a );
  zMulMatMat( l, u, a );
  printf( "LxU: " ); zMatWrite( a );
  printf( "(rank = %d)\n", rank );

  zMatFree( a );
  zMatFree( l );
  zMatFree( u );
  zIndexFree( index );
  return 0;
}
