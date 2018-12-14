#include <zm/zm_eig.h>

#define N 200
#define M 100

int main(int argc, char *argv[])
{
  zMat m, u, v, l, tmp1, tmp2;
  zVec sv;
  int i, rank;

  zRandInit();
  m = zMatAlloc( N, M );
  u = zMatAllocSqr( N );
  v = zMatAlloc( N, M );
  sv = zVecAlloc( N );
  zMatRandUniform( m, -10, 10 );

  printf( "maximum singular value = %.10g\n", zSVMax( m ) );
  printf( "minimum singular value = %.10g\n", zSVMin( m ) );

  rank = zSVD( m, sv, u, v );
  zVecWrite( sv );

  printf( ">>ensurance\n" );
  l = zMatAlloc( N, rank );
  tmp1 = zMatAlloc( N, rank );
  tmp2 = zMatAlloc( N, M );
  for( i=0; i<rank; i++ )
    zMatSetElem( l, i, i, zVecElem(sv,i) );
  zMulMatMat( u, l, tmp1 );
  zMulMatMat( tmp1, v, tmp2 );
  zMatSubDRC( tmp2, m );
  printf( "|| ULV - A || = %g\n", zMatNorm(tmp2) );
  printf( " ... %s.\n", zMatIsTiny(tmp2) ? "OK" : "maybe a bug" );

  zMatFree( m );
  zMatFree( u );
  zMatFree( v );
  zMatFree( l );
  zMatFree( tmp1 );
  zMatFree( tmp2 );
  zVecFree( sv );
  return 0;
}
