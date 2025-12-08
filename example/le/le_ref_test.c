#include <zm/zm_le.h>

void test(zMat a, zVec b, zVec w, zVec d, zVec x, zVec y)
{
  zVec tmp;

  tmp = zVecAlloc( zVecSizeNC(d) );
  zLESolveRefNormMin( a, b, w, d, x );
  zMulMatVec( a, x, y );
  zVecSubDRC( y, b );
  zVecSub( d, x ,tmp );
  zVecAmp( tmp, w, x );
  printf( "|residual|=%.10g |error|=%.10g\n", zVecInnerProd(tmp,x), zVecNorm(y) );
  zVecFree( tmp );
}

#define N 100
#define M 300
int main(void)
{
  zMat a;
  zVec b, w, x, d, y;
  int i;

  zRandInit();
  a = zMatAlloc( N, M ); zMatRandUniform( a, -10, 10 );
  b = zVecAlloc( N ); zVecRandUniform( b, -10, 10 );
  w = zVecAlloc( M );
  d = zVecAlloc( M );
  x = zVecAlloc( M );
  y = zVecAlloc( N );

  zLESolveNormMin( a, b, NULL, d );
  for( i=0; i<zVecSizeNC(d); i++ )
    zVecElemNC(d,i) += zRandF(-10.0,10.0);

  zVecSetAll( w, 1.0 );
  test( a, b, w, d, x, y );

  zVecLinSpace( w, 1.0, zVecSize(w) );
  test( a, b, w, d, x, y );

  zMatFree( a );
  zVecFree( b );
  zVecFree( w );
  zVecFree( x );
  zVecFree( d );
  zVecFree( y );
  return 0;
}
