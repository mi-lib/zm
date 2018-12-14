#include <zm/zm_mat.h>

int main(void)
{
  zVec v1, v2;
  zMat m1, m2, m;
  double m1array[] = {
    2, 1,-1,
    1, 3, 2,
   -1,-1, 3 };
  double m2array[] = {
    1, 2, 1,
    2,-1,-2,
    3, 1,-1 };
  double v1array[] = { 1, 2, 3 };
  double v2array[] = { 2,-1,-2 };

  v1 = zVecCloneArray( v1array, 3 );
  v2 = zVecCloneArray( v2array, 3 );
  m1 = zMatCloneArray( m1array, 3, 3 );
  m2 = zMatCloneArray( m2array, 3, 3 );

  printf( "v1 = " ); zVecWrite( v1 );
  printf( "v2 = " ); zVecWrite( v2 );
  printf( "m1 = " ); zMatWrite( m1 );
  printf( "m2 = " ); zMatWrite( m2 );

  zMulMatVecDRC( m1, v1 );
  zMulVecMatDRC( v2, m1 );
  printf( "m1 * v1 = " ); zVecWrite( v1 );
  printf( "v2 * m1 = " ); zVecWrite( v2 );

  m = zMatAlloc( 3, 3 );
  zMulMatMat( m1, m2, m );
  printf( "m1 * m2 = " ); zMatWrite( m );
  zMatFree( m );

  m = zMatAlloc( 3, 3 );
  zMulMatTMat( m1, m2, m );
  printf( "m1^T * m2 = " ); zMatWrite( m );
  zMatFree( m );

  m = zMatAlloc( 3, 3 );
  zMulMatMatT( m1, m2, m );
  printf( "m1 * m2^T = " ); zMatWrite( m );
  zMatFree( m );

  zVecFree( v1 );
  zVecFree( v2 );
  zMatFree( m1 );
  zMatFree( m2 );
  return 0;
}
