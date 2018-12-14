#include <zm/zm_le.h>

#define N 100
#define MAX (1.0e12)

int main(int argc, char *argv[])
{
  zMat a;
  zVec b, x, ans, err;

  a = zMatAllocSqr( N );
  b = zVecAlloc( N );
  x = zVecAlloc( N );
  ans = zVecAlloc( N );
  err = zVecAlloc( N );

  zRandInit();
  zMatRand( a, -MAX, MAX );
  zVecRand( b, -MAX, MAX );

  printf( "LU decomposition method.\n" );
  zLESolveLU( a, b, x );
  zLEResidual( a, b, x, err );
  printf( "error = %.12g\n", zVecNorm( err ) );

  printf( "Gauss's elimination method.\n" );
  zLESolveGauss( a, b, x );
  zLEResidual( a, b, x, err );
  printf( "error = %.12g\n", zVecNorm( err ) );

  printf( "LU decomp. and residual iteration method.\n" );
  zLESolveRI( a, b, x );
  zLEResidual( a, b, x, err );
  printf( "error = %.12g\n", zVecNorm( err ) );

  zMatFree( a );
  zVecFree( b );
  zVecFree( x );
  zVecFree( ans );
  zVecFree( err );
  return 0;
}
