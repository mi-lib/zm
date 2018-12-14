#include <zm/zm_le.h>

int main(int argc, char *argv[])
{
  unsigned s;
  zMat m, ma, testm;

  zRandInit();

  s = argc > 1 ? atoi( argv[1] ) : 10;
  m = zMatAllocSqr( s );
  zMatRandUniform( m, -10, 10 );
  ma = zMatAllocSqr( s );
  testm = zMatAllocSqr( s );
  printf( "matrix: " ); zMatWrite( m );

  zMatAdj( m, ma );
  printf( "adjoint matrix: " ); zMatWrite( ma );

  zMulMatMat( m, ma, testm );
  zMatTouchup( testm );
  printf( "org x adj: " ); zMatWrite( testm );
  printf( "det = %.10g\n", zMatDet( m ) );

  zMatFree( m );
  zMatFree( ma );
  zMatFree( testm );
  return 0;
}
