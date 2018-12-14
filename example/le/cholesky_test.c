#include <zm/zm_le.h>

void test_mat(zMat m, zMat r, int n, double degratio)
{
  int i, j;

  for( i=0; i<n; i++ ){
    if( zRandF(0,1) > degratio )
      for( j=0; j<=i; j++ )
        zMatSetElem( r, i, j, zRandF(-10,10) );
    else
      for( j=0; j<=i; j++ )
        zMatSetElem( r, i, j, zMatElem(r,0,j)+zMatElem(r,i-1,j) );
  }
  zMulMatMatT( r, r, m );
}

#define N 100

int main(int argc, char *argv[])
{
  int i, j, k, n=10, rank;
  zMat m, l, r;
  zIndex index;

  m = zMatAllocSqr( n );
  r = zMatAllocSqr( n );
  zRandInit();
  for( i=j=k=0; i<N; i++ ){
    test_mat( m, r, n, 0.5 );
    if( ( rank = zCholeskyDecompAlloc( m, &l, &index ) ) < 0 ) continue;
    j++;
    zMulMatMatT( l, l, r );
    zMatSubDRC( r, m );
    printf( "rank=%d, error=%.16g ... %s.\n", rank, zMatNorm(r), zMatIsTiny(r) ? "ok" : "failure" );
    if( zMatIsTiny(r) ) k++;
  }
  printf( "success rate = %d/%d\n", k, j );

  zMatFreeAO( 3, m, r, l );
  zIndexFree( index );
  return 0;
}
