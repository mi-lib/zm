#include <zm/zm_mat.h>

#define N 3
#define M 3

int main(void)
{
  zMat m;
  zVec v;
  int i, j;

  zRandInit();
  m = zMatAlloc( N, M );
  v = zVecAlloc( M );
  zMatRand( m, -10, 10 );
  zMatWrite( m );
  printf( "||m||_2 = %g\n", zMatNorm(m) );
  printf( "||m||_inf = %g\n", zMatInfNorm(m) );
  for( i=0; i<N; i++ ){
    for( j=0; j<M; j++ )
      zVecSetElem( v, j, zSgn(zMatElem(m,i,j)) );
    printf( "[%d] %g\n", i, zRawVecInnerProd(zMatRowArray(m,i),zVecArray(v),M) );
  }
  zMatFree( m );
  zVecFree( v );
  return 0;
}
