#include <zm/zm_le.h>

int main(int argc, char *argv[])
{
  int rank;
  zMat a1, q, r, tq, a2;
  zIndex idx;

  zRandInit();
  a1 = zMatAlloc( 20, 30 );
  a2 = zMatClone( a1 );
  zMatRandUniform( a1, -10, 10 );

  q = zMatAlloc( _zMatRowSize(a1), zMin(_zMatRowSize(a1),_zMatColSize(a1)) );
  r = zMatAlloc( _zMatColSize(q), _zMatColSize(a1) );
  idx = zIndexCreate( _zMatRowSize(a1) );

  printf( ">>> QR decomposition <<<\n" );
  rank = zQRDecomp( a1, q, r, idx );
  zMatTouchup( q );
  zMatTouchup( r );
  tq= zMatAllocSqr( _zMatRowSize(q) );
  printf( "rank=%d\n", rank );
  printf( "R=\n" );
  zMatImg( r );
  printf( "Q^T Q=\n" );
  zMulMatTMat( q, q, tq );
  zMatTouchup( tq );
  zMatWrite( tq );

  zMulMatMat( q, r, a2 );
  printf( "%s.\n", zMatIsEqual(a1,a2) ? "ok" : "error" );

  zMatFree( q );
  zMatFree( r );
  zMatFree( tq );
  zMatFree( a1 );
  zMatFree( a2 );
  return 0;
}
