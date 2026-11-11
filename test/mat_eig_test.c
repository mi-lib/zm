#include <zm/zm.h>
#include <float.h>

#define N 100

/* tolerance of the normwise relative residual in the unit of n eps */
#define TOL_EIG 10.0

bool check_mat_hessenberg(zMat h)
{
  int i, j;

  for( i=2; i<zMatRowSizeNC(h); i++ )
    for( j=0; j<i-1; j++ )
      if( zMatElemNC(h,i,j) != 0 ) return false;
  return true;
}

/* check if p h p^T = m and p^T p = I with relative errors. */
bool check_mat_hessenberg_decomp(zMat m, zMat h, zMat p)
{
  zMat tmp, mc;
  int n;
  double norm;
  bool result = true;

  n = zMatRowSizeNC(m);
  tmp = zMatAllocSqr( n );
  mc = zMatAllocSqr( n );
  if( !check_mat_hessenberg( h ) ){
    eprintf( "not a Hessenberg matrix.\n" );
    result = false;
  }
  zMulMatMatT( h, p, tmp );
  zMulMatMat( p, tmp, mc );
  zMatSubDRC( mc, m );
  norm = zMatNorm( m );
  if( zMatNorm( mc ) > TOL_EIG * n * DBL_EPSILON * norm ){
    eprintf( "relative error = %g\n", zMatNorm(mc) / norm );
    result = false;
  }
  zMulMatTMat( p, p, tmp );
  zMatSubDRC( tmp, zMatIdentNC( mc ) );
  if( zMatNorm( tmp ) > TOL_EIG * n * DBL_EPSILON ){
    eprintf( "orthogonality error = %g\n", zMatNorm(tmp) );
    result = false;
  }
  zMatFreeAtOnce( 2, tmp, mc );
  return result;
}

void assert_mat_hessenberg(void)
{
  zMat m, h, p;
  const int size = 20;
  double scale[] = { 1.0, 1.0e-10, 1.0e-13, 1.0e-100, 1.0e10, 1.0e100 };
  bool result = true, result_scale = true;
  int i, k;

  m = zMatAllocSqr( size );
  h = zMatAllocSqr( size );
  p = zMatAllocSqr( size );
  for( i=0; i<N; i++ ){
    zMatRandUniform( m, -10, 10 );
    zMatToHessenberg( m, h, p );
    if( !check_mat_hessenberg_decomp( m, h, p ) ) result = false;
    for( k=0; k<sizeof(scale)/sizeof(double); k++ ){ /* scaled matrices */
      zMatMulDRC( m, scale[k] );
      if( !zMatToHessenberg( m, h, p ) || !check_mat_hessenberg_decomp( m, h, p ) ){
        eprintf( "failed for a matrix scaled by %g\n", scale[k] );
        result_scale = false;
      }
      zMatRandUniform( m, -10, 10 );
    }
  }
  zMatFreeAtOnce( 3, m, h, p );
  zAssert( zMatToHessenberg, result );
  zAssert( zMatToHessenberg (scaled matrix), result_scale );
}

/* check the normwise relative residuals of eigenpairs and the normalization of eigenvectors. */
bool check_mat_eig(zMat m, zCVec eigval, zCMat eigbase, double tol)
{
  int i, n;
  zCMat cm;
  zCVec eigvec, r;
  zComplex z;
  double norm;
  bool result = true;

  n = zMatRowSizeNC(m);
  cm = zCMatAllocSqr( n );
  eigvec = zCVecAlloc( n );
  r = zCVecAlloc( n );
  zMatToCMat( m, cm );
  norm = zMatNorm( m );
  for( i=0; i<n; i++ ){
    zCMatGetCol( eigbase, i, eigvec );
    if( !zIsTol( zCVecNorm(eigvec) - 1, TOL_EIG * n * DBL_EPSILON ) ){
      eprintf( "eigenvector #%d not normalized: norm = %g\n", i, zCVecNorm(eigvec) );
      result = false;
    }
    zCMulMatVec( cm, eigvec, r );
    zComplexRev( zCVecElemNC(eigval,i), &z );
    zCVecCatDRC( r, &z, eigvec );
    if( zCVecNorm( r ) > tol * n * DBL_EPSILON * norm ){
      eprintf( "eigenpair #%d: relative residual = %g\n", i, zCVecNorm(r) / norm );
      result = false;
    }
  }
  zCMatFree( cm );
  zCVecFree( eigvec );
  zCVecFree( r );
  return result;
}

bool assert_mat_eig_one(double a[], int n)
{
  zCVec eigval;
  zMat ma;
  zCMat eigbase;
  bool result;

  ma = zMatCloneArray( a, n, n );
  eigval = zCVecAlloc( n );
  eigbase = zCMatAllocSqr( n );
  result = zMatEig( ma, eigval, eigbase, 0 ) && check_mat_eig( ma, eigval, eigbase, TOL_EIG );
  zMatFree( ma );
  zCMatFree( eigbase );
  zCVecFree( eigval );
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

/* check if eigenvalues coincide with expected ones, regardless of the order. Each eigenvalue is
 * matched with the nearest one among the expected eigenvalues that are not matched yet. */
bool check_mat_eig_val(zCVec eigval, zCVec expected, double tol)
{
  int i, j, jmin, n;
  bool *used, result = true;
  double d, dmin;
  zComplex diff;

  n = zCVecSizeNC(eigval);
  used = zAlloc( bool, n );
  for( i=0; i<n; i++ ){
    for( dmin=HUGE_VAL, jmin=-1, j=0; j<n; j++ ){
      if( used[j] ) continue;
      zComplexSub( zCVecElemNC(eigval,i), zCVecElemNC(expected,j), &diff );
      if( ( d = zComplexAbs( &diff ) ) < dmin ){
        dmin = d;
        jmin = j;
      }
    }
    used[jmin] = true;
    if( dmin > tol ){
      eprintf( "eigenvalue %g%+gi mismatches %g%+gi\n",
        zCVecElemNC(eigval,i)->re, zCVecElemNC(eigval,i)->im, zCVecElemNC(expected,jmin)->re, zCVecElemNC(expected,jmin)->im );
      result = false;
    }
  }
  zFree( used );
  return result;
}

void assert_mat_eig_scale(void)
{
  double a[] = {
    1, 2, 3,
    4, 5, 6,
    7, 8,10,
  };
  const int n = 3;
  zMat m;
  zCVec eigval;
  zCMat eigbase;
  double s;
  int i;
  bool result = true;

  m = zMatCloneArray( a, n, n );
  eigval = zCVecAlloc( n );
  eigbase = zCMatAllocSqr( n );
  for( s=1.0e-300; s<1.0e300; s*=1.0e20 ){
    zMatCopyArray( a, n, n, m );
    zMatMulDRC( m, s );
    if( !zMatEig( m, eigval, eigbase, 0 ) ){
      eprintf( "failed for a matrix scaled by %g\n", s );
      result = false;
      continue;
    }
    /* eigenvectors are independent of the scale, and eigenvalues are proportional to it. */
    zMatCopyArray( a, n, n, m );
    for( i=0; i<n; i++ ){
      zCVecElemNC(eigval,i)->re /= s;
      zCVecElemNC(eigval,i)->im /= s;
    }
    if( !check_mat_eig( m, eigval, eigbase, TOL_EIG ) ){
      eprintf( "inaccurate for a matrix scaled by %g\n", s );
      result = false;
    }
  }
  zMatFree( m );
  zCVecFree( eigval );
  zCMatFree( eigbase );
  zAssert( zMatEig (scaled matrix), result );
}

/* matrices that are difficult for the double-shift QR method. */
bool assert_mat_eig_special_one(zMat m, zCVec expected, double tol)
{
  zCVec eigval;
  zCMat eigbase;
  bool result;

  eigval = zCVecAlloc( zMatRowSizeNC(m) );
  eigbase = zCMatAllocSqr( zMatRowSizeNC(m) );
  result = zMatEig( m, eigval, eigbase, 0 ) && check_mat_eig( m, eigval, eigbase, TOL_EIG );
  if( result && expected ) result = check_mat_eig_val( eigval, expected, tol );
  zCVecFree( eigval );
  zCMatFree( eigbase );
  zMatFree( m );
  zCVecFree( expected );
  return result;
}

void assert_mat_eig_special(void)
{
  zMat m;
  zCVec expected;
  const int n = 8;
  double coeff[] = { 1, -21, 175, -735, 1624, -1764, 720 }; /* (x-1)(x-2)...(x-6) */
  int i;
  bool result_cyclic, result_companion, result_jordan, result_trivial = true;

  /* cyclic permutation matrix, for which the shifts stagnate without exceptional shifts */
  m = zMatAllocSqr( n );
  expected = zCVecAlloc( n );
  for( i=0; i<n; i++ ){
    zMatSetElemNC( m, (i+1)%n, i, 1 );
    zComplexCreatePolar( zCVecElemNC(expected,i), 1, 2*zPI*i/n );
  }
  result_cyclic = assert_mat_eig_special_one( m, expected, TOL_EIG * n * DBL_EPSILON );
  /* companion matrix */
  m = zMatAllocSqr( 6 );
  expected = zCVecAlloc( 6 );
  for( i=0; i<6; i++ ){
    zMatSetElemNC( m, 0, i, -coeff[i+1] );
    if( i > 0 ) zMatSetElemNC( m, i, i-1, 1 );
    zComplexCreate( zCVecElemNC(expected,i), i+1, 0 );
  }
  result_companion = assert_mat_eig_special_one( m, expected, 1.0e-9 );
  /* Jordan block, a defective matrix */
  m = zMatAllocSqr( 4 );
  for( i=0; i<4; i++ ){
    zMatSetElemNC( m, i, i, 2 );
    if( i > 0 ) zMatSetElemNC( m, i-1, i, 1 );
  }
  result_jordan = assert_mat_eig_special_one( m, NULL, 0 );
  /* zero and identity matrices */
  if( !assert_mat_eig_special_one( zMatAllocSqr( 5 ), NULL, 0 ) ) result_trivial = false;
  if( !assert_mat_eig_special_one( zMatIdentNC( zMatAllocSqr( 5 ) ), NULL, 0 ) ) result_trivial = false;
  zAssert( zMatEig (cyclic permutation matrix), result_cyclic );
  zAssert( zMatEig (companion matrix), result_companion );
  zAssert( zMatEig (Jordan block), result_jordan );
  zAssert( zMatEig (zero and identity matrices), result_trivial );
}

void assert_mat_eig_random(void)
{
  zMat m, d, q, qinv;
  zCVec eigval, eigval_only;
  zCMat eigbase;
  int i, j, n;
  bool result = true, result_val = true, result_cluster = true;

  for( i=0; i<N; i++ ){
    n = zRandI( 1, 30 );
    m = zMatAllocSqr( n );
    eigval = zCVecAlloc( n );
    eigval_only = zCVecAlloc( n );
    eigbase = zCMatAllocSqr( n );
    zMatRandUniform( m, -10, 10 );
    if( !zMatEigDQR( m, eigval, eigbase, 0 ) || !check_mat_eig( m, eigval, eigbase, TOL_EIG ) )
      result = false;
    /* eigenvalues computed without eigenvectors */
    if( !zMatEigDQR( m, eigval_only, NULL, 0 ) ||
        !check_mat_eig_val( eigval_only, eigval, TOL_EIG * n * DBL_EPSILON * zMatNorm(m) ) )
      result_val = false;
    zMatFree( m );
    zCVecFree( eigval );
    zCVecFree( eigval_only );
    zCMatFree( eigbase );
  }
  for( i=0; i<N; i++ ){ /* clusters of multiple eigenvalues */
    n = zRandI( 2, 15 );
    m = zMatAllocSqr( n );
    d = zMatAllocSqr( n );
    q = zMatAllocSqr( n );
    qinv = zMatAllocSqr( n );
    eigval = zCVecAlloc( n );
    eigbase = zCMatAllocSqr( n );
    for( j=0; j<n; j++ ) zMatSetElemNC( d, j, j, zRandI(-2,2) );
    do{
      zMatRandUniform( q, -1, 1 );
    } while( !zMatInv( q, qinv ) );
    zMulMatMat( q, d, m );
    zMulMatMat( m, qinv, d );
    if( !zMatEigDQR( d, eigval, eigbase, 0 ) || !check_mat_eig( d, eigval, eigbase, 10*TOL_EIG ) )
      result_cluster = false;
    zMatFreeAtOnce( 4, m, d, q, qinv );
    zCVecFree( eigval );
    zCMatFree( eigbase );
  }
  zAssert( zMatEigDQR (random), result );
  zAssert( zMatEigDQR (eigenvalues only), result_val );
  zAssert( zMatEigDQR (multiple eigenvalues), result_cluster );
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
  assert_mat_eig_scale();
  assert_mat_eig_special();
  assert_mat_eig_random();
  assert_mat_sym_eig_power();
  assert_mat_sym_eig();
  assert_mat_sym_eig_random();
  assert_mat_svd();
  assert_mat_svd_random();
  assert_mat_svd_minmax();
  return 0;
}
