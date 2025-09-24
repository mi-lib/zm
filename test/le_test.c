#include <zm/zm.h>

#define N 40
#define M 80
#define R 20

void assert_mat_mat_sweep_out(void)
{
  zMat m1, m2, m3, m;
  zIndex index;
  bool result = true;
  const int size = 10, testnum = 1000;
  const double tol = 1.0e-6;
  int i, k;

  m  = zMatAllocSqr( size );
  m1 = zMatAllocSqr( size );
  m2 = zMatAllocSqr( zMatRowSizeNC(m1) );
  m3 = zMatAlloc( zMatRowSizeNC(m), zMatColSizeNC(m2) );
  index = zIndexCreate( zMatRowSizeNC(m1) );
  for( k=0; k<testnum; k++ ){
    zMatRandUniform( m1, -5, 5 );
    zMatCopy( m1, m );
    zMatIdent( m2 );
    for( i=0; i<zMatRowSizeNC(m1); i++ )
      zMatMatSweepOut( m1, m2, zIndexElemNC(index,i), i );
    zMulMatMat( m, m2, m3 );
    if( !zMatIsIdent( m3, tol ) ){
      eprintf( "maximum error = %g\n", zMatElemAbsMax( m3, NULL ) );
      result = false;
    }
  }
  zMatFreeAtOnce( 4, m, m1, m2, m3 );
  zIndexFree( index );
  zAssert( zMatMatSweepOut, result );
}

void assert_mat_vec_sweep_out(void)
{
  zMat m;
  zVec v, x;
  zIndex index;
  bool result = true;
  const double tol = 1.0e-6;
  const int size = 10, testnum = 1000;
  int i, k;

  m = zMatAllocSqr( size );
  v = zVecAlloc( size );
  x = zVecAlloc( size );
  index = zIndexCreate( size );
  for( k=0; k<testnum; k++ ){
    zMatRandUniform( m, -5, 5 );
    zVecRandUniform( v, -5, 5 );
    zLESolveGauss( m, v, x );
    for( i=0; i<zMatRowSize(m); i++ ){
      zMatPivoting( m, index, i, zIndexElem(index,i) );
      zMatVecSweepOut( m, v, zIndexElem(index,i), i );
    }
    zVecReorderDRC( v, index );
    zVecSubDRC( v, x );
    if( !zVecIsTol( v, tol ) ){
      eprintf( "maximum error = %g\n", zVecElemAbsMax( v, NULL ) );
      result = false;
    }
  }
  zMatFree( m );
  zVecFree( v );
  zIndexFree( index );
  zAssert( zMatVecSweepOut, result );
}

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
    if( !zMatIsIdent( m, zTOL ) ){
      zMatFPrint( stderr, m );
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
      eprintf( "case #%d: maximum error = %g\n", k, zMatElemAbsMax(m2,NULL) );
      result = false;
    }
  }
  zMatFreeAtOnce( 4, m, m1, m2, m3 );
  zAssert( zMulInvMatMat, result );
}

zMat make_mpinv_testcase1(void)
{
  /* regular case */
  double marray[] = {
    1, 2,
   -2,-4 };
  return zMatCloneArray( marray, 2, 2 );
}

zMat make_mpinv_testcase2(void)
{
  /* regular case (row excess) */
  double marray[] = {
    2, 1, 1,
    1, 1, 0 };
  return zMatCloneArray( marray, 2, 3 );
}

zMat make_mpinv_testcase3(void)
{
  /* regular case (row excess) */
  double marray[] = {
    4, 5,-2, 4,
    5, 3, 4, 3,
   -3, 1, 2,-4,
    8, 2, 1, 1,
   -1,-3, 2, 1,
    1, 3, 0,-2 };
  return zMatCloneArray( marray, 6, 4 );
}

zMat make_mpinv_testcase4(void)
{
  /* singular case (column excess) */
  double marray[] = {
    4, 5,-2, 4, 5, 3,
    4, 3,-3, 1, 2,-4,
    8, 2, 1, 1,-1,-3,
    4, 5,-2, 4, 5, 3,
    0, 2, 1, 3, 3, 7 };
  return zMatCloneArray( marray, 5, 6 );
}

zMat make_mpinv_testcase_regular_rand(void)
{
  /* singular case (column excess) */
  double marray[] = {
    4, 5,-2, 4, 5, 3,
    4, 3,-3, 1, 2,-4,
    8, 2, 1, 1,-1,-3,
    4, 5,-2, 4, 5, 3,
    0, 2, 1, 3, 3, 7 };
  return zMatCloneArray( marray, 5, 6 );
}

bool assert_mat_mpinv_one(zMat (* testcase)(void))
{
  zMat a, ai, b, d, e;
  bool result = true;

  a = testcase();
  ai = zMatAlloc( zMatColSizeNC(a), zMatRowSizeNC(a) );
  zMatMPInv( a, ai );

  b = zMatAllocSqr( zMatRowSizeNC(a) );
  d = zMatAlloc( zMatRowSizeNC(a), zMatColSizeNC(a) );
  e = zMatAlloc( zMatRowSizeNC(a), zMatColSizeNC(a) );
  zMulMatMat( a, ai, b );
  zMulMatMat( b, a, d );
  zMatSub( a, d, e );
  if( !zMatIsTiny( e ) ){
    printf( "|| A - AA^+A || = %.10f\n", zMatNorm(e) );
    result = false;
  }
  zMatFreeAtOnce( 3, b, d, e );

  b = zMatAllocSqr( zMatColSizeNC(a) );
  d = zMatAlloc( zMatColSizeNC(a), zMatRowSizeNC(a) );
  e = zMatAlloc( zMatColSizeNC(a), zMatRowSizeNC(a) );
  zMulMatMat( ai, a, b );
  zMulMatMat( b, ai, d );
  zMatSub( ai, d, e );
  if( !zMatIsTiny( e ) ){
    printf( "|| A^+ - A^+AA^+ || = %.10f\n", zMatNorm(e) );
    result = false;
  }
  zMatFreeAtOnce( 5, b, d, e, a, ai );
  return result;
}

void assert_mat_mpinv(void)
{
  zAssert( zMatMPInv,
    assert_mat_mpinv_one( make_mpinv_testcase1 ) &&
    assert_mat_mpinv_one( make_mpinv_testcase2 ) &&
    assert_mat_mpinv_one( make_mpinv_testcase3 ) &&
    assert_mat_mpinv_one( make_mpinv_testcase4 ) );
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

void assert_mat_mpinv_rand(void)
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
    mat_deg( m1, R );
    zMatRandUniform( m2, -10, 10 );
    zMatGetCol( m2, 0, v1 );

    zMatMPInv( m1, mp );
    zMulMatMat( mp, m2, ma );
    zMulMPInvMatMat( m1, m2, mb );
    if( !zMatEqual(ma,mb,tol) ) result1 = false;

    zLESolveMP( m1, v1, NULL, NULL, v2 );
    zMulMatVec( mp, v1, v3 );
    if( !zVecEqual(v2,v3,tol) ) result2 = false;
  }
  zAssert( zMatMPInv + zMulMPInvMatMat (random test), result1 );
  zAssert( zLESolveMP (random test), result2 );

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
  zMatMPInvNull( a, mp, mn );

  v = zVecAlloc( M );
  u = zVecAlloc( M );
  e = zVecAlloc( N );
  zVecRandUniform( v, -10, 10 );
  zMulMatVec( mn, v, u );
  zMulMatVec( a, u, e );
  zAssert( zMatMPInvNull, zVecIsTol(e,tol) );
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
  const double tol = 1.0e-10;

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
    if( !zIsTol( err, tol ) ){
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
    zLELyapnovSolve( a, b, x );
    zMulMatMat( x, a, tmp1 );
    zMulMatTMat( a, x, tmp2 );
    zMatAddDRC( tmp1, tmp2 );
    zMatSubDRC( tmp1, b );
    if( !zMatIsTol( tmp1, tol ) ){
      eprintf( "case #%d: maximum error = %g\n", k, zMatElemAbsMax(tmp1,NULL) );
      result = false;
    }
  }
  zMatFreeAtOnce( 5, a, b, x, tmp1, tmp2 );
  zAssert( zLELyapnovSolve, result );
}

void assert_mat_det_adj(void)
{
  zMat m, minv1, minv2;
  double det;
  const int size = 10, testnum = 100;
  int k;
  bool result = true;

  m = zMatAllocSqr( size );
  minv1 = zMatAllocSqr( size );
  minv2 = zMatAllocSqr( size );
  for( k=0; k<testnum; k++ ){
    zMatRandUniform( m, -10, 10 );
    zMatInv( m, minv1 );
    det = zMatDet( m );
    zMatAdj( m, minv2 );
    zMatDivDRC( minv2, det );
    if( !zMatEqual( minv1, minv2, zTOL ) ){
      zMatSubDRC( minv1, minv2 );
      eprintf( "case #%d: maximum error = %g\n", k, zMatElemAbsMax(minv1,NULL) );
      result = false;
    }
  }
  zMatFreeAtOnce( 3, m, minv1, minv2 );
  zAssert( zMatAdj + zMatDet, result );
}

int main(void)
{
  zRandInit();
  assert_mat_mat_sweep_out();
  assert_mat_vec_sweep_out();
  assert_mat_inv();
  assert_mat_mul_inv();
  assert_mat_mpinv();
  assert_mat_mpinv_rand();
  assert_mpnull();
  assert_tridiagonal_equation();
  assert_lyapnov_equation();
  assert_mat_det_adj();
  return EXIT_SUCCESS;
}
