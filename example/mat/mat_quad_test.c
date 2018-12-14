#include <zm/zm_mat.h>

#define N 200
#define M 300

int main(void)
{
  zMat a, wm, q1, q2, wat;
  zVec w;

  zRandInit();
  a = zMatAlloc( N, M );
  wm = zMatAllocSqr( M );
  q1 = zMatAllocSqr( N );
  q2 = zMatAllocSqr( N );
  wat = zMatAlloc( M, N );
  w = zVecAlloc( M );

  zVecRand( w, 0.1, 10 );
  zMatRand( a, -10, 10 );
  zMatDiag( wm, w );
#if 0
printf("A "); zMatWrite( a );
printf("w "); zVecWrite( w );
printf("W "); zMatWrite( wm );
#endif

  zMatQuad( a, w, q1 );
#if 0
printf("Q1 "); zMatWrite( q1 );
#endif
  zMulMatMatT( wm, a, wat );
#if 0
printf("WAT "); zMatWrite( wat );
#endif
  zMulMatMat( a, wat, q2 );
#if 0
printf("Q2 "); zMatWrite( q2 );
#endif

  printf( "%s\n", zBoolExpr( zMatIsEqual( q1, q2 ) ) );
  zMatSubDRC( q1, q2 );
  printf( "err = %g\n", zMatNorm(q1) );
  printf( "|max| = %g\n", zDataAbsMax(zMatArray(q1),zMatRowSize(q1)*zMatColSize(q1),NULL) );

  zMatFree( a );
  zMatFree( wm );
  zMatFree( q1 );
  zMatFree( q2 );
  zMatFree( wat );
  zVecFree( w );
  return 0;
}
