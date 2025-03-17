#include <zm/zm.h>

#define N 40
#define M 80
#define R 20

void assert_mat_inv(void)
{
  zMat m1, m2, m;
  const int size = 10, testnum = 100;
  bool result = true;
  int k;

  m1 = zMatAllocSqr( size );
  m2 = zMatAllocSqr( size );
  m  = zMatAllocSqr( size );
  for( k=0; k<testnum; k++ ){
    zMatRandUniform( m1, -10, 10 );
    zMatInv( m1, m2 );
    zMulMatMat( m1, m2, m );
    if( !zMatIsIdent( m ) ){
      zMatPrint( m );
      result = false;
    }
  }
  zMatFreeAtOnce( 3, m1, m2, m );
  zAssert( zMatInv, result );
}

void assert_mat_mul_inv(void)
{
  zMat m1, m2, m3, m;
  const int rowsize = 10, colsize = 8, testnum = 100;
  bool result = true;
  int k;
  const double tol = 1.0e-10;

  m1 = zMatAlloc( rowsize, colsize );
  m2 = zMatAlloc( rowsize, colsize );
  m3 = zMatAlloc( rowsize, colsize );
  m = zMatAllocSqr( rowsize );
  for( k=0; k<testnum; k++ ){
    zMatRandUniform( m, -10, 10 );
    zMatRandUniform( m1,-10, 10 );
    zMulMatMat( m, m1, m2 );
    zMulInvMatMat( m, m2, m3 );
    zMatSub( m1, m3, m2 );
    if( !zMatIsTol( m2, tol ) ){
      eprintf( "case #%d: maximum error = %g\n", k, zDataAbsMax(zMatBuf(m2),zMatRowSizeNC(m2)*zMatColSizeNC(m2),NULL) );
      result = false;
    }
  }
  zMatFreeAtOnce( 4, m, m1, m2, m3 );
  zAssert( zMulInvMatMat, result );
}

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

void assert_mpinv(void)
{
  zMat mp, m1, m2, ma, mb;
  zVec v1, v2, v3;
  int i;
  bool result1, result2, result3;
  const double tol = 1.0e-10;

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
    if( !zMatEqual(ma,mb,tol) ) result1 = false;

    zLESolveMP( m1, v1, NULL, NULL, v2 );
    zMulMatVec( mp, v1, v3 );
    if( !zVecEqual(v2,v3,tol) ) result2 = false;
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
  const double tol = 1.0e-10;

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
  zAssert( zMPInvNull, zVecIsTol(e,tol) );
  zMatFreeAtOnce( 3, a, mp, mn );
  zVecFreeAtOnce( 3, v, u, e );
}

void assert_tridiagonal_equation(void)
{
  zVec a, b, c, d, ans;
  int i;
  double err;
  const int dim = 100;
  bool result = true;

  a = zVecAlloc( dim );
  b = zVecAlloc( dim );
  c = zVecAlloc( dim );
  d = zVecAlloc( dim );
  ans = zVecAlloc( dim );
  zVecRandUniform( a, -10, 10 );
  zVecRandUniform( b, -10, 10 );
  zVecRandUniform( c, -10, 10 );
  zVecRandUniform( d, -10, 10 );
  zLETridiagSolve( a, b, c, d, ans );

  for( i=0; i<zVecSizeNC(ans); i++ ){
    err = zVecElem(b,i)*zVecElem(ans,i);
    if( i > 0 )
      err += zVecElem(a,i)*zVecElem(ans,i-1);
    if( i < zVecSizeNC(a)-1 )
      err += zVecElem(c,i)*zVecElem(ans,i+1);
    err -= zVecElem(d,i);
    if( !zIsTiny( err ) ){
      eprintf( "residual #%d is not tiny: %g\n", i, err );
      result = false;
    }
  }
  zVecFreeAtOnce( 5, a, b, c, d, ans );
  zAssert( zLETridiagSolve, result );
}

void assert_lyapnov_equation(void)
{
  zMat a, b, x, tmp1, tmp2;
  const int size = 10, testnum = 100;
  const double tol = 1.0e-9;
  int k;
  bool result = true;

  a = zMatAllocSqr( size );
  b = zMatAllocSqr( size );
  x = zMatAllocSqr( size );
  tmp1 = zMatAllocSqr( size );
  tmp2 = zMatAllocSqr( size );
  for( k=0; k<testnum; k++ ){
    zMatRandUniform( a, -10, 10 );
    zMatRandUniform( b, -10, 10 );
    zLyapnovSolve( a, b, x );
    zMulMatMat( x, a, tmp1 );
    zMulMatTMat( a, x, tmp2 );
    zMatAddDRC( tmp1, tmp2 );
    zMatSubDRC( tmp1, b );
    if( !zMatIsTol( tmp1, tol ) ){
      eprintf( "case #%d: maximum error = %g\n", k, zDataAbsMax(zMatBuf(tmp1),zMatRowSizeNC(tmp1)*zMatColSizeNC(tmp1),NULL) );
      result = false;
    }
  }
  zMatFreeAtOnce( 5, a, b, x, tmp1, tmp2 );
  zAssert( zLyapnovSolve, result );
}

int main(void)
{
  zRandInit();
  assert_mat_inv();
  assert_mat_mul_inv();
  assert_mpinv();
  assert_mpnull();
  assert_tridiagonal_equation();
  assert_lyapnov_equation();
  return EXIT_SUCCESS;
}
