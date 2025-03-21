#include <zm/zm.h>

#define N 100

bool check_mat_hessenberg(zMat m)
{
  int i, j, i1;

  for( i=2; i<zMatRowSizeNC(m); i++ ){
    i1 = i - 1;
    for( j=0; j<i1; j++ )
      if( !zIsTiny( zMatElemNC(m,i,j) ) ) return false;
  }
  return true;
}

void assert_mat_hessenberg(void)
{
  zMat m, tm, p, tmp, mc;
  const int size = 20;
  bool result = true;
  int i;

  m = zMatAllocSqr( size );
  tm = zMatAllocSqr( size );
  p  = zMatAllocSqr( size );
  tmp = zMatAllocSqr( size );
  mc = zMatAllocSqr( size );

  for( i=0; i<N; i++ ){
    zMatRandUniform( m, -10, 10 );
    zMatToHessenberg( m, tm, p );
    if( !check_mat_hessenberg( tm ) ){
      eprintf( "not a hessenberg matrix.\n" );
      zMatImg( tm );
      result = false;
    }
    zMulMatMatT( tm, p, tmp );
    zMulMatMat( p, tmp, mc );
    zMatSubDRC( mc, m );
    if( !zMatIsTiny( mc ) ){
      eprintf( "error matrix norm = %g\n", zMatNorm(mc) );
      result = false;
    }
  }
  zMatFreeAtOnce( 5, m, tm, p, tmp, mc );
  zAssert( zMatToHessenberg, result );
}

bool assert_mat_eig_one(double a[], int n)
{
  int i;
  zCVec eigval, eigvec, vt;
  zMat ma;
  zCMat cma, eigbase;
  zComplex z;
  bool result = true;
  const double tol = 1.0e-7;

  ma = zMatCloneArray( a, n, n );
  cma = zCMatAllocSqr( n );
  eigval = zCVecAlloc( n );
  eigbase = zCMatAllocSqr( n );
  eigvec = zCVecAlloc( n );
  vt = zCVecAlloc( n );
  zMatEig( ma, eigval, eigbase, 0 );
  /* confirmation */
  zMatToCMat( ma, cma );
  for( i=0; i<n; i++ ){
    zCMatGetCol( eigbase, i, eigvec );
    zCMulMatVec( cma, eigvec, vt );
    zComplexRev( zCVecElemNC(eigval,i), &z );
    zCVecCatDRC( vt, &z, eigvec );
    if( !zCVecIsTol( vt, tol ) ){
      eprintf( " error norm = %g\n", zCVecNorm(vt) );
      result = false;
    }
  }
  zMatFree( ma );
  zCMatFree( cma );
  zCMatFree( eigbase );
  zCVecFree( eigval );
  zCVecFree( eigvec );
  zCVecFree( vt );
  return result;
}

bool assert_mat_eig_testcase1(void)
{
  int n = 3;
  double a[] = {
    0.0,-6.0, -1.0,
    6.0, 2.0,-16.0,
   -5.0,20.0,-10.0,
  };
  /* eigv = -3.0710, -2.4645+17.6008i, -2.4645-17.6008i
     eigm =
       -0.8326  0.2003-0.1394i   0.2003+0.1394i
       -0.3553 -0.2110-0.6447i  -0.2110+0.6447i
       -0.4248 -0.6930          -0.6930
   */
  return assert_mat_eig_one( a, n );
}

bool assert_mat_eig_testcase2(void)
{
  int n = 3;
  double a[] = {
    6.0, 12.0, 19.0,
   -9.0,-20.0,-33.0,
    4.0,  9.0, 15.0,
  };
  /* eigv = -1, 1, 1
     eigm =
       -0.4741   -0.4082   -0.4082
        0.8127    0.8165    0.8165
       -0.3386   -0.4082   -0.4082
   */
  return assert_mat_eig_one( a, n );
}

bool assert_mat_eig_testcase3(void)
{
  int n = 3;
  double a[] = {
    3.0, 0.0, 0.0,
    0.0, 2.0,-5.0,
    0.0, 1.0,-2.0,
  };
  /* eigv = 3, i, -i */
  return assert_mat_eig_one( a, n );
}

bool assert_mat_eig_testcase4(void)
{
  int n = 3;
  double a[] = {
    0, 1, 1,
   -4, 4, 2,
    4,-3,-1,
  };
  /* eigv = 0, 1, 2 */
  return assert_mat_eig_one( a, n );
}

bool assert_mat_eig_testcase5(void)
{
  int n = 3;
  double a[] = {
    1, 0,-1,
   -2, 1, 3,
    2, 1, 2,
  };
  /* eigv = 3, w, w^2 */
  return assert_mat_eig_one( a, n );
}

bool assert_mat_eig_testcase6(void)
{
  int n = 3;
  double a[] = {
    1, 4,-4,
   -1,-3, 2,
    0, 2,-1,
  };
  /* eigv =-3,-1, 1
     eigm =
       -0.8165   0.8165   0.8944
        0.4082  -0.4082   0
       -0.4082  -0.4082   0.4472
   */
  return assert_mat_eig_one( a, n );
}

bool assert_mat_eig_testcase7(void)
{
  int n = 3;
  double a[] = {
    2, 1, 0,
    0, 1,-1,
    0, 2, 4,
  };
  /* eigv = 2, 2, 3,
     eigm =
        1  1 -0.4082
        0  0 -0.4082
        0  0  0.8165
   */
  return assert_mat_eig_one( a, n );
}

bool assert_mat_eig_testcase8(void)
{
  int n = 4;
  double a[] = {
    2, 1, 5, 1,
    1, 3, 7, 0,
    0, 0, 2, 1,
    2, 4, 1, 4,
  };
  /* eigv = 7.03608, 1.38087, 1.29152-2.62195i, 1.29152+2.62195i */
  return assert_mat_eig_one( a, n );
}

#define ASSERT_MAT_EIG_TESTCASE(n) do{\
  if( !assert_mat_eig_testcase##n() ){\
    eprintf( "failure in testcase "#n "\n" );\
    result = false;\
  }\
} while(0)

void assert_mat_eig(void)
{
  bool result = true;

  ASSERT_MAT_EIG_TESTCASE( 1 );
  ASSERT_MAT_EIG_TESTCASE( 2 );
  ASSERT_MAT_EIG_TESTCASE( 3 );
  ASSERT_MAT_EIG_TESTCASE( 4 );
  ASSERT_MAT_EIG_TESTCASE( 5 );
  ASSERT_MAT_EIG_TESTCASE( 6 );
  ASSERT_MAT_EIG_TESTCASE( 7 );
  ASSERT_MAT_EIG_TESTCASE( 8 );
  zAssert( zMatEig, result );
}

void assert_mat_sym_eig_power(void)
{
  zMat l, m, eigbase;
  zVec v, eigval, error;
  double s, s_ref;
  const int size = 10;
  const double tol = 1.0e-8;
  bool result_max = true, result_min = true;
  int i;

  l = zMatAllocSqr( size );
  m = zMatAllocSqr( size );
  eigbase = zMatAllocSqr( size );
  v = zVecAlloc( size );
  eigval = zVecAlloc( size );
  error = zVecAlloc( size );
  for( i=0; i<N; i++ ){
    zMatRandUniform( l, -10, 10 );
    zMulMatTMat( l, l, m );
    zMatSymEigJacobi( m, eigval, eigbase );
    s = zMatEigPower( m, v, 0 );
    zMulMatVec( m, v, error );
    zVecCatDRC( error, -s, v );
    if( !zEqual( s, ( s_ref = zVecElemMax(eigval,NULL) ), tol ) ){
      eprintf( "mistamch the greatest eigenvalue: %g / %g\n", s, s_ref );
      result_max = false;
    }
    s = zMatEigPowerInv( m, v, 0 );
    zMulMatVec( m, v, error );
    zVecCatDRC( error, -s, v );
    if( !zEqual( s, ( s_ref = zVecElemMin(eigval,NULL) ), tol ) ){
      eprintf( "mistamch the least eigenvalue: %g / %g\n", s, s_ref );
      result_min = false;
    }
  }
  zMatFreeAtOnce( 3, l, m, eigbase );
  zVecFreeAtOnce( 3, v, eigval, error );
  zAssert( zMatEigPower, result_max );
  zAssert( zMatEigPowerInv, result_min );
}

bool check_mat_sym_eig(zMat m, zVec eigval, zMat eigbase, int c)
{
  zVec v, e;
  int i;
  bool result = true;

  v = zVecAlloc( zVecSizeNC(eigval) );
  e = zVecAlloc( zVecSizeNC(eigval) );
  for( i=0; i<zVecSizeNC(eigval); i++ ){
    zMatGetCol( eigbase, i, v );
    zMulMatVec( m, v, e );
    zVecCatDRC( e, -zVecElem(eigval,i), v );
    if( !zVecIsTiny( e ) ){
      printf( "eig#%d: eig.val. = %g,  error = %g\n", i, zVecElem(eigval,i), zVecNorm(e) );
      result = false;
    }
  }
  if( c > 0 ) eprintf( "computation time(clk) = %d\n", c );
  zVecFree( v );
  zVecFree( e );
  return result;
}

void assert_mat_sym_eig(void)
{
#define R 1.41421356237
  double array[] = {
    1, R, R, 2,
    R,-R,-1, R,
    R,-1, R, R,
    2, R, R,-3,
  };
  zMat m, eigbase;
  zVec eigval;
  const int n = 4;
  bool result;

  m = zMatCloneArray( array, n, n );
  eigbase = zMatAllocSqr( n );
  eigval = zVecAlloc( n );
  zMatSymEigJacobi( m, eigval, eigbase );
  result = check_mat_sym_eig( m, eigval, eigbase, -1 );
  zMatFree( m );
  zMatFree( eigbase );
  zVecFree( eigval );
  zAssert( zMatSymEigJacobi, result );
}

void assert_mat_sym_eig_random(void)
{
  int n, i, j;
  zMat m, eigbase;
  zVec eigval;
  clock_t c1, c2;
  bool result_bisec, result_jacobi;

  n = N;
  m = zMatAllocSqr( n );
  eigbase = zMatAllocSqr( n );
  eigval = zVecAlloc( n );
  for( i=0; i<n; i++ )
    for( j=i; j<n; j++ ){
      zMatSetElem( m, i, j, zRandF(-10,10) );
      zMatSetElem( m, j, i, zMatElem(m,i,j) );
    }

  c1 = clock();
  zMatSymEigBisec( m, eigval, eigbase );
  c2 = clock();
  result_bisec = check_mat_sym_eig( m, eigval, eigbase, c2-c1 );

  c1 = clock();
  zMatSymEigJacobi( m, eigval, eigbase );
  c2 = clock();
  result_jacobi = check_mat_sym_eig( m, eigval, eigbase, c2-c1 );

  zMatFree( m );
  zMatFree( eigbase );
  zVecFree( eigval );

  zAssert( zMatSymEigBisec (random), result_bisec );
  zAssert( zMatSymEigJacobi (random), result_jacobi );
}

bool assert_mat_svd_one(double a[], int n, int m)
{
  zMat ma, u, v, s, tmp1, tmp2;
  zVec sv;
  int i, rank;
  bool result = true;

  ma = zMatCloneArray( a, n, m );
  u = zMatAllocSqr( n );
  v = zMatAlloc( n, m );
  sv = zVecAlloc( n );

  rank = zMatSVD( ma, u, sv, v );
/*
  zMatPrint( ma );
  zVecPrint( sv );
  zMatPrint( u );
  zMatPrint( v );
*/
  s = zMatAlloc( n, rank );
  tmp1 = zMatAlloc( n, rank );
  tmp2 = zMatAlloc( n, m );
  for( i=0; i<rank; i++ )
    zMatSetElem( s, i, i, zVecElem(sv,i) );
  zMulMatMat( u, s, tmp1 );
  zMulMatMat( tmp1, v, tmp2 );
  if( !zMatEqual( tmp2, ma, zTOL ) ){
    zMatSubDRC( tmp2, ma );
    eprintf( "|| USV - A || = %g\n", zMatNorm(tmp2) );
    result = false;
  }
  zMatFreeAtOnce( 6, ma, u, v, s, tmp1, tmp2 );
  zVecFree( sv );
  return result;
}

bool assert_mat_svd_testcase1(void)
{
  double a[] = {
    2, 3,-2,
   -3,14,-7,
   -5,19,-9,
  };
  int n = 3, m = 3;
  /* ans =
     | 0.11 -0.97  0.21 || 27.02 0    0    || -0.21  0.88 -0.43 |
     | 0.59  0.11 -0.80 ||  0    2.82 0    ||  0.96 -0.10  0.26 |
     | 0.80  0.22  0.56 ||  0    0    0.16 ||  0.19  0.46  0.87 |
   */
  return assert_mat_svd_one( a, n, m );
}

bool assert_mat_svd_testcase2(void)
{
  double a[] = {
    9, 4,
    6, 8,
    2, 7,
  };
  int n = 3, m = 2;
  /* ans =
     | -0.6105  0.7174  0.3355 || 14.9359 0      || -0.6925 -0.7214 |
     | -0.6646 -0.2336 -0.7098 ||  0      5.1883 ||  0.7214 -0.6925 |
     | -0.4308 -0.6563  0.6194 ||  0      0      |
   */
  return assert_mat_svd_one( a, n, m );
}

bool assert_mat_svd_testcase3(void)
{
  double a[] = {
    1, 2,
    3, 4,
    5, 6,
    7, 8,
  };
  int n = 4, m = 2;
  /* ans =
 | 0.1525  0.8226 -0.3945   -0.3800 ||
 | 0.3499  0.4214  0.2428    0.8007 || 14.2691 0      || 0.6414 -0.7672 |
 | 0.5474  0.0201  0.6979   -0.4614 ||  0      0.6268 || 0.7672  0.6414 |
 | 0.7448 -0.3812 -0.5462    0.0407 ||  0      0      |
   */
  return assert_mat_svd_one( a, n, m );
}

bool assert_mat_svd_testcase4(void)
{
  double a[] = {
    1, 2, 3, 4, 5,
    5, 6, 7, 8, 9,
    9,10,11,12,13,
  };
  int n = 3, m = 5;
  return assert_mat_svd_one( a, n, m );
}

void assert_mat_svd(void)
{
  bool result = true;

  if( !assert_mat_svd_testcase1() ) result = false;
  if( !assert_mat_svd_testcase2() ) result = false;
  if( !assert_mat_svd_testcase3() ) result = false;
  if( !assert_mat_svd_testcase4() ) result = false;
  zAssert( zMatSVD, result );
}

void assert_mat_svd_random(void)
{
  zMat a;
  bool result = true;
  int i;

  for( i=0; i<N; i++ ){
    a = zMatAlloc( zRandI(1,40), zRandI(1,40) );
    zMatRandUniform( a, -10, 10 );
    if( !assert_mat_svd_one( zMatBuf(a), zMatRowSize(a), zMatColSize(a) ) ) result = false;
    zMatFree( a );
  }
  zAssert( zMatSVD (random case), result );
}

void assert_mat_svd_minmax(void)
{
  zMat m, u, v;
  zVec sv;
  int i, rank;
  const int rowsize = 20, colsize = 10;
  double s_max, s_min;
  bool result = true;

  m = zMatAlloc( rowsize, colsize );
  u = zMatAllocSqr( rowsize );
  v = zMatAlloc( rowsize, colsize );
  sv = zVecAlloc( rowsize );
  for( i=0; i<N; i++ ){
    zMatRandUniform( m, -10, 10 );
    s_max = zMatSingularValueMax( m );
    s_min = zMatSingularValueMin( m );
    rank = zMatSVD( m, u, sv, v );
    if( !zEqual( s_max, zVecElemNC(sv,0), zTOL ) ){
      eprintf( "mismatch the greatest singular value: %g / %g\n", s_max, zVecElemNC(sv,0) );
      result = false;
    }
    if( !zEqual( s_min, zVecElemNC(sv,rank-1), zTOL ) ){
      eprintf( "mismatch the least singular value: %g / %g\n", s_min, zVecElemNC(sv,rank-1) );
      result = false;
    }
  }
  zMatFree( m );
  zMatFree( u );
  zMatFree( v );
  zVecFree( sv );
  zAssert( zMatSingularValueMax & zMatSingularValueMin, result );
}

int main(int argc, char *argv[])
{
  zRandInit();
  assert_mat_hessenberg();
  assert_mat_eig();
  assert_mat_sym_eig_power();
  assert_mat_sym_eig();
  assert_mat_sym_eig_random();
  assert_mat_svd();
  assert_mat_svd_random();
  assert_mat_svd_minmax();
  return 0;
}
