#include <zm/zm_le.h>

void test_lu(int r, int c, int n)
{
  zMat m, mc, l, u;
  zIndex index;
  int i, count_fail = 0;

  m = zMatAlloc( r, c );
  mc = zMatAlloc( r, c );
  l = zMatAlloc( r, r );
  u = zMatAlloc( r, c );
  index = zIndexAlloc( zMatColSizeNC(l) );
  for( i=0; i<n; i++ ){
    zMatRandUniform( m, -10, 10 );
    zMatDecompLU( m, l, u, index );
    zMulMatMat( l, u, mc );
    zMatSubDRC( m, mc );
    if( !zMatIsTiny( m ) ) count_fail++;
  }
  printf( "failure rate (%d x %d): LU decomposition %d/%d\n", r, c, count_fail, n );
  zMatFreeAtOnce( 4, m, mc, l, u );
  zIndexFree( index );
}

void test_cholesky(int r, int c, int n)
{
  int i;
  zMat m, mc, l, s;
  zIndex index;
  int count_fail = 0;

  m = zMatAllocSqr( r );
  mc = zMatAllocSqr( r );
  l = zMatAllocSqr( r );
  s = zMatAlloc( r, c );
  index = zIndexAlloc( r );
  for( i=0; i<n; i++ ){
    zMatRandUniform( s, -10, 10 );
    zMulMatMatT( s, s, m );
    zMatDecompCholesky( m, l, index );
    zMulMatMatT( l, l, mc );
    zMatSubDRC( m, mc );
    if( !zMatIsTiny( m ) ) count_fail++;
  }
  printf( "failure rate (%d x %d) (%d x %d): Cholesky decomposition %d/%d\n", r, c, c, r, count_fail, n );
  zMatFreeAtOnce( 4, m, mc, l, s );
  zIndexFree( index );
}

bool test_qqt(zMat m)
{
  register int i, j;

  if( !zMatIsSqr( m ) ){
    ZRUNERROR( ZM_ERR_MAT_NOTSQR );
    return false;
  }
  for( i=0; i<zMatRowSizeNC(m); i++ )
    for( j=0; j<zMatColSizeNC(m); j++ )
      if( !( i == j && zIsTiny( zMatElemNC(m,i,j) - 1 ) ) && !zIsTiny( zMatElemNC(m,i,j) ) )
        return false;
  return true;
}

void test_lq(int r, int c, int n)
{
  zMat m, mc, l, q, e;
  zIndex index;
  int i, rank, count_fail = 0;

  m = zMatAlloc( r, c );
  mc = zMatAlloc( r, c );
  l = zMatAlloc( r, r );
  q = zMatAlloc( r, c );
  e = zMatAllocSqr( r );
  index = zIndexAlloc( zMatColSizeNC(l) );
  for( i=0; i<n; i++ ){
    zMatRandUniform( m, -10, 10 );
    zMatSetColSize( l, r );
    zMatSetRowSize( q, r );
    rank = zMatDecompLQ( m, l, q, index );
    zMatColReg( l, rank );
    zMatRowReg( q, rank );
    zMulMatMat( l, q, mc );
    zMatSubDRC( m, mc );
    zMatSetRowSize( e, rank );
    zMatSetColSize( e, rank );
    zMulMatMatT( q, q, e );
    if( !zMatIsTiny( m ) || !test_qqt( e ) ) count_fail++;
  }
  printf( "failure rate (%d x %d): LQ decomposition %d/%d\n", r, c, count_fail, n );
  zMatFreeAtOnce( 5, m, mc, l, q, e );
  zIndexFree( index );
}

#define N 100

int main(void)
{
  zRandInit();
  test_cholesky( 8, 5, N );
  test_cholesky( 8, 8, N );
  test_lu( 5, 8, N );
  test_lu( 8, 5, N );
  test_lu( 8, 8, N );
  test_lq( 5, 8, N );
  test_lq( 8, 5, N );
  test_lq( 8, 8, N );
  return 0;
}
