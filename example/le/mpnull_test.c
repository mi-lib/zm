#include <zm/zm_le.h>

#define N 30
#define M 50
#define R 20

void mat_deg(zMat m, int rank)
{
  int i, r1, r2;
  double s1, s2;

  for( i=rank; i<_zMatRowSize(m); i++ ){
    r1 = zRandI( 0, rank-1 );
    r2 = zRandI( 0, rank-1 );
    s1 = zRandF( -10, 10 );
    s2 = zRandF( -10, 10 );
    zRawVecLS( zMatRowBuf(m,i), _zMatColSize(m), 2, s1, zMatRowBuf(m,r1), s2, zMatRowBuf(m,r2) );
  }
}

int main(void)
{
  int rank;
  zMat a, mp, mn;
  zVec v, u, e;

  zRandInit();
  a = zMatAlloc( N, M );
  zMatRandUniform( a, -10, 10 );
  mat_deg( a, R );
  mp = zMatAlloc( M, N );
  mn = zMatAlloc( M, M );

  rank = zMPInvNull( a, mp, mn );
  zMatWrite( mp );
  zMatWrite( mn );

  /* assertion */
  v = zVecAlloc( M );
  u = zVecAlloc( M );
  e = zVecAlloc( N );
  zVecRandUniform( v, -10, 10 );
  zMulMatVec( mn, v, u );
  zMulMatVec( a, u, e );
  zVecWrite( u );
  zVecWrite( e );
  printf( "rank=%d, %s <%g,%g>\n", rank, zBoolExpr( zVecIsTiny(e) ), zVecMin(e,NULL), zVecMax(e,NULL) );

  zMatFreeAO( 3, a, mp, mn );
  zVecFreeAO( 3, v, u, e );
  return 0;
}
