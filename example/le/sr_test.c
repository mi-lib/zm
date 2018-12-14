#include <zm/zm_le.h>
#include <zm/zm_opt.h>
#include <liw/liw_debug.h>

zVec zLESolveSRCG(zMat a, zVec b, zVec wn, zVec we, zVec ans)
{
  register int i, j;
  zMat m1, m2;
  zVec v;

  m1 = zMatAlloc( _zMatColSize(a), _zMatRowSize(a) );
  m2 = zMatAllocSqr( zVecSizeNC(ans) );
  v = zVecAlloc( zVecSizeNC(ans) );
  for( i=0; i<_zMatRowSize(a); i++ )
    for( j=0; j<_zMatColSize(a); j++ )
      zMatSetElem( m1, j, i, zMatElem(a,i,j)*zVecElem(we,i) );
  zMulMatMatNC( m1, a, m2 );
  for( i=0; i<_zMatRowSize(m2); i++ )
    zMatElem(m2,i,i) += zVecElem( wn, i );
  zMulMatVecNC( m1, b, v );
  zVecRevDRC( v );
  zCGSolve( m2, v, ans, 0 );
  return ans;
}

void set_weight(zVec w, double v)
{
  register int i;

  for( i=0; i<zVecSize(w); i++ ) zVecSetElem( w, i, v );
}

#define ROW 20
#define COL 10

int main(void)
{
  double v;
  zMat a;
  zVec b, wn, we, x1, x2;
  clock_t clk;

  a = zMatAlloc( ROW, COL );
  zMatRand( a, -10, 10 );
  b = zVecAlloc( ROW );
  zVecRand( b, -10, 10 );
  wn = zVecAlloc( COL );
  we = zVecAlloc( ROW );
  zVecSetAll( we, 1.0 );
  x1 = zVecAlloc( COL );
  x2 = zVecAlloc( COL );

  zMatWrite( a );
  zVecWrite( b );
  zVecWrite( we );
  for( v=1; v>0.0000000001; v*=0.1 ){
    set_weight( wn, v );
    MEASURE_EXEC_CLOCK( &clk, zLESolveSR( a, b, wn, we, x1 ) );
    printf( "(direct solution) ... %ld\n", clk );
    zVecWrite( x1 );

    MEASURE_EXEC_CLOCK( &clk, zLESolveSRCG( a, b, wn, we, x2 ) );
    printf( "(indirect solution with CG method) ... %ld\n", clk );
    zVecWrite( x2 );

    if( zVecDist( x1, x2 ) > zTOL )
      ZRUNWARN( "possibly a bug inside of zLESolveSR" );
    getchar();
  }
  zLESolveNormMin( a, b, wn, x1 );
  zVecWrite( x1 );
  return 0;
}
