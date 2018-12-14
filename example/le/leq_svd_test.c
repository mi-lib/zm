#include <zm/zm_le.h>

#define M 300
#define N 400

#define D_MAX 10.0

void test(zMat a, zVec b, zVec x, zVec _b)
{
  register int i, j;

  for( i=0; i<M; i++ ){
    for( j=0; j<N; j++ )
      zMatSetElem( a, i, j, zRandF(-D_MAX,D_MAX) );
    zVecSetElem( b, i, zRandF(-D_MAX,D_MAX) );
  }
  zLESolveMP_SVD( a, b, x );
  zMulMatVec( a, x, _b );
  zVecSub( b, _b, _b );
  printf( "err = %g\n", zVecNorm( _b ) );
}

int main(int argc, char *argv[])
{
  zMat a;
  zVec b, x, _b;

  zRandInit();
  a = zMatAlloc( M, N );
  b = zVecAlloc( M );
  _b= zVecAlloc( M );
  x = zVecAlloc( N );

  test( a, b, x, _b );

  zMatFree( a );
  zVecFree( b );
  zVecFree( _b );
  zVecFree( x );
  return 0;
}
