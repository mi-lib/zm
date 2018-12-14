#include <zm/zm_le.h>

#define TEST 2

int main(void)
{
  /* singular case (column excess) */
#if TEST == 1
  double a_arr[] = {
    2, 1, 1,
    1, 2, 3,
    1, 2, 2,
    2, 2, 2 };
  double b_arr[] = {
    1, 7, 5, 4 };
  unsigned r = 4, c = 3;
#elif TEST == 2
  double a_arr[] = {
    2, 1, 1,
    1, 2, 1 };
  double b_arr[] = {
    3, 2 };
  unsigned r = 2, c = 3;
#else /* TEST == 0 */
  double a_arr[] = {
    4, 5,-2, 4, 5, 3,
    4, 3,-3, 1, 2,-4,
    8, 2, 1, 1,-1,-3,
    4, 5,-2, 4, 5, 3,
    0, 2, 1, 3, 3, 7 };
  double b_arr[] = {
    5, 4, 3, 4, 1 };
  unsigned r = 5, c = 6;
#endif
  double w1_arr[] = {
    1, 1, 1, 1, 1, 1 };
  double w2_arr[] = {
    1, 1, 1, 1, 1 };

  zMat a;
  zVec b, w1, w2, ans;

  a = zMatCloneArray( a_arr, r, c );
  b = zVecCloneArray( b_arr, r );
  w1 = zVecCloneArray( w1_arr, c );
  w2 = zVecCloneArray( w2_arr, r );
  ans = zVecAlloc( c );

  printf( "MP solve (LU)\n" );
  zLESolveMP_LU( a, b, w1, w2, ans );
  printf( "A: " ); zMatWrite( a );
  printf( "b: " ); zVecWrite( b );
  printf( "x: " ); zVecWrite( ans );
  zMulMatVec( a, ans, b );
  printf( "A x: " ); zVecWrite( b );

  printf( "MP solve (QR)\n" );
  zLESolveMP( a, b, w1, w2, ans );
  printf( "A: " ); zMatWrite( a );
  printf( "b: " ); zVecWrite( b );
  printf( "x: " ); zVecWrite( ans );
  zMulMatVec( a, ans, b );
  printf( "A x: " ); zVecWrite( b );

  return 0;
}
