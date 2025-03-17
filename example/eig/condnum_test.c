#include <zm/zm_mat_eig.h>

zMat hilbert_mat(zMat m)
{
  int i, j;

  for( i=0; i<zMatRowSizeNC(m); i++ )
    for( j=0; j<zMatColSizeNC(m); j++ )
      zMatElemNC(m,i,j) = 1.0 / ( i + j + 1 );
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
    printf( "dim=%2d, smax=%.10g, smin=%1.10g, kappa=%1.10g\n", i, zMatSingularValueMax(m), zMatSingularValueMin(m), zMatCondNum(m) );
    zMatFree( m );
  }
  return 0;
}
