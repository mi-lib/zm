#include <zm/zm_le.h>

#define N 10

int main(int argc, char *argv[])
{
  zMat m1, m2, m;
  int n;
  clock_t c1, c2;

  n = argc > 1 ? atoi( argv[1] ) : N;

  m1 = zMatAllocSqr( n );
  m2 = zMatAllocSqr( n );
  m  = zMatAllocSqr( n );
  zMatRand( m1, -10, 10 );

  printf( "m2 = m1^-1\n" );
  c1 = clock();
  zMatInv( m1, m2 );
  c2 = clock();
  zMulMatMat( m1, m2, m );
  zMatTouchup( m );
  zMatWrite( m );
  printf( "err=%g, clk = %ld\n", zMatSqrNorm(m)/n, c2-c1 );

  zMatFree( m1 );
  zMatFree( m2 );
  zMatFree( m );
  return 0;
}
