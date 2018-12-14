#include <zm/zm_le.h>

int main(int argc, char *argv[])
{
  int r = 30, c = 40;
  int rank;
  zMat a1, a2, l, q, tq;
  zIndex idx;

  zRandInit();
  a1 = zMatAlloc( r, c );
  zMatRandUniform( a1, -10, 10 );
  a2 = zMatClone( a1 );
  printf( ">>> LQ decomposition <<<\n" );
  rank = zLQDecompAlloc( a1, &l, &q, &idx );
  tq= zMatAllocSqr( _zMatRowSize(q) );

  zMatTouchup( l );
  zMatTouchup( q );
  printf( "rank=%d\n", rank );
  printf( "L=\n" );
  zMatImg( l );
  printf( "Q Q^T=\n" );
  zMulMatMatT( q, q, tq );
  zMatTouchup( tq );
  zMatWrite( tq );

  zMulMatMat( l, q, a2 );
  printf( "%s.\n", zMatIsEqual(a1,a2) ? "ok" : "error" );

  zMatFree( l );
  zMatFree( q );
  zMatFree( tq );
  zMatFree( a1 );
  zMatFree( a2 );
  zIndexFree( idx );
  return 0;
}
