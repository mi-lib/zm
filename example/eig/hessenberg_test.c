#include <zm/zm_eig.h>

#define N 10

int main(void)
{
  zMat m, tm, p, tmp, mc;

  zRandInit();
  m = zMatAllocSqr( N );
  zMatRandUniform( m, -10, 10 );

  tm = zMatAllocSqr( N );
  p  = zMatAllocSqr( N );
  printf( "original matrix(M)\n" );
  zMatImg( m );

  zHess( m, tm, p );
  printf( "Hessenberg matrix(H)\n" );
  zMatImg( tm );

  tmp = zMatAllocSqr( N );
  mc = zMatAllocSqr( N );
  zMulMatMatT( tm, p, tmp );
  zMulMatMat( p, tmp, mc );
  zMatSubDRC( mc, m );
  printf( ">>ensurance<< || PHP^T - M || = %.16f\n", zMatNorm(mc) );

  zMatFree( m );
  zMatFree( tm );
  zMatFree( p );
  zMatFree( tmp );
  zMatFree( mc );
  return 0;
}
