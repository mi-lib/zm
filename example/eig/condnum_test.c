#include <zm/zm_eig.h>

zMat hilbert_mat(zMat m)
{
  register int i, j;

  for( i=0; i<_zMatRowSize(m); i++ )
    for( j=0; j<_zMatColSize(m); j++ )
      zMatElem(m,i,j) = 1.0 / ( i + j + 1 );
  return m;
}

#define N 10

int main(int argc, char *argv[])
{
  zMat m;
  int i;

  for( i=3; i<=N; i++ ){
    m = zMatAllocSqr( i );
    hilbert_mat( m );
    printf( "dim=%2d, smax=%.10g, smin=%1.10g, kappa=%1.10g\n", i, zSVMax(m), zSVMin(m), zMatCondNum(m) );
    zMatFree( m );
  }
  return 0;
}
