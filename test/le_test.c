#include <zm/zm.h>

#define N 40
#define M 80
#define R 20

void mat_deg(zMat m, int rank)
{
  int i, r1, r2;
  double s1, s2;

  for( i=rank; i<zMatRowSizeNC(m); i++ ){
    r1 = zRandI( 0, rank-1 );
    r2 = zRandI( 0, rank-1 );
    s1 = zRandF( -10, 10 );
    s2 = zRandF( -10, 10 );
    zRawVecLinearSum( zMatRowBufNC(m,i), zMatColSizeNC(m), 2, s1, zMatRowBufNC(m,r1), s2, zMatRowBufNC(m,r2) );
  }
}

#define TOL (1.0e-10)

void assert_mpinv(void)
{
  zMat mp, m1, m2, ma, mb;
  zVec v1, v2, v3;
  int i;
  bool result1, result2, result3;

  mp = zMatAlloc( M, N );
  m1 = zMatAlloc( N, M );
  m2 = zMatAlloc( N, R );
  ma = zMatAlloc( M, R );
  mb = zMatAlloc( M, R );
  v1 = zVecAlloc( N );
  v2 = zVecAlloc( M );
  v3 = zVecAlloc( M );
  result1 = result2 = result3 = true;
  for( i=0; i<N; i++ ){
    zMatRandUniform( m1, -10, 10 );
    zMatRandUniform( m2, -10, 10 );
    zMatGetCol( m2, 0, v1 );

    zMPInv( m1, mp );
    zMulMatMat( mp, m2, ma );
    zMulMPInvMatMat( m1, m2, mb );
    if( !zMatEqual(ma,mb,TOL) ) result1 = false;

    zLESolveMP( m1, v1, NULL, NULL, v2 );
    zMulMatVec( mp, v1, v3 );
    if( !zVecEqual(v2,v3,TOL) ) result2 = false;
  }
  zAssert( zMPInv + zMulMPInvMatMat, result1 );
  zAssert( zLESolveMP, result2 );

  zMatFreeAtOnce( 5, m1, m2, ma, mb, mp );
  zVecFreeAtOnce( 3, v1, v2, v3 );
}

void assert_mpnull(void)
{
  zMat a, mp, mn;
  zVec v, u, e;

  a = zMatAlloc( N, M );
  zMatRandUniform( a, -10, 10 );
  mat_deg( a, R );
  mp = zMatAlloc( M, N );
  mn = zMatAlloc( M, M );
  zMPInvNull( a, mp, mn );

  v = zVecAlloc( M );
  u = zVecAlloc( M );
  e = zVecAlloc( N );
  zVecRandUniform( v, -10, 10 );
  zMulMatVec( mn, v, u );
  zMulMatVec( a, u, e );
  zAssert( zMPInvNull, zVecIsTol(e,TOL) );
  zMatFreeAtOnce( 3, a, mp, mn );
  zVecFreeAtOnce( 3, v, u, e );
}

int main(void)
{
  zRandInit();
  assert_mpinv();
  assert_mpnull();
  return EXIT_SUCCESS;
}
