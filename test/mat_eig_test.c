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
    if( !( zCVecNorm( r ) <= tol * n * DBL_EPSILON * norm ) ){ /* negated to reject NaN */
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
    if( jmin < 0 ){ /* the eigenvalue is NaN */
      eprintf( "eigenvalue %g%+gi is invalid\n", zCVecElemNC(eigval,i)->re, zCVecElemNC(eigval,i)->im );
      result = false;
      continue;
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

/* eigenvectors are independent of the scale, and eigenvalues are proportional to it. */
bool assert_mat_eig_scale_one(double a[], int n, double s, zMat m, zCVec eigval, zCMat eigbase)
{
  int i;

  zMatCopyArray( a, n, n, m );
  zMatMulDRC( m, s );
  if( !zMatEig( m, eigval, eigbase, 0 ) ){
    eprintf( "failed for a matrix scaled by %g\n", s );
    return false;
  }
  zMatCopyArray( a, n, n, m );
  for( i=0; i<n; i++ ){
    zCVecElemNC(eigval,i)->re /= s;
    zCVecElemNC(eigval,i)->im /= s;
  }
  if( !check_mat_eig( m, eigval, eigbase, TOL_EIG ) ){
    eprintf( "inaccurate for a matrix scaled by %g\n", s );
    return false;
  }
  return true;
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
  bool result = true;

  m = zMatAllocSqr( n );
  eigval = zCVecAlloc( n );
  eigbase = zCMatAllocSqr( n );
  for( s=1.0e-300; s<1.0e300; s*=1.0e20 )
    if( !assert_mat_eig_scale_one( a, n, s, m, eigval, eigbase ) ) result = false;
  /* the maximum component 1e308 is greater than 2^1023, while the eigenvalues do not overflow. */
  if( !assert_mat_eig_scale_one( a, n, 1.0e307, m, eigval, eigbase ) ) result = false;
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

/* check if m eigbase = eigbase diag{eigval} and eigbase^T eigbase = I with relative errors. */
bool check_mat_sym_eig_rel(zMat m, zVec eigval, zMat eigbase, double tol)
{
  int i, n;
  zMat vtv;
  zVec v, e;
  double norm;
  bool result = true;

  n = zMatRowSizeNC(m);
  vtv = zMatAllocSqr( n );
  v = zVecAlloc( n );
  e = zVecAlloc( n );
  norm = zMatNorm( m );
  for( i=0; i<n; i++ ){
    zMatGetCol( eigbase, i, v );
    zMulMatVec( m, v, e );
    zVecCatDRC( e, -zVecElem(eigval,i), v );
    if( !( zVecNorm( e ) <= tol * n * DBL_EPSILON * norm ) ){ /* negated to reject NaN */
      eprintf( "eigenpair #%d: relative residual = %g\n", i, zVecNorm(e) / norm );
      result = false;
    }
  }
  zMulMatTMat( eigbase, eigbase, vtv );
  for( i=0; i<n; i++ ) zMatElemNC(vtv,i,i) -= 1;
  if( !( zMatNorm( vtv ) <= tol * n * DBL_EPSILON ) ){
    eprintf( "eigenvectors not orthonormal: |V^T V - I| = %g\n", zMatNorm(vtv) );
    result = false;
  }
  zMatFree( vtv );
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
  result = zMatSymEigJacobi( m, eigval, eigbase ) && check_mat_sym_eig_rel( m, eigval, eigbase, TOL_EIG );
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
  result_jacobi = zMatSymEigJacobi( m, eigval, eigbase );
  c2 = clock();
  eprintf( "computation time(clk) = %d\n", (int)( c2 - c1 ) );
  if( result_jacobi ) result_jacobi = check_mat_sym_eig_rel( m, eigval, eigbase, TOL_EIG );

  zMatFree( m );
  zMatFree( eigbase );
  zVecFree( eigval );

  zAssert( zMatSymEigBisec (random), result_bisec );
  zAssert( zMatSymEigJacobi (random), result_jacobi );
}

/* eigenvectors are independent of the scale, and eigenvalues are proportional to it. */
bool assert_mat_sym_eig_jacobi_scale_one(double a[], int n, double s, zMat m, zVec eigval, zMat eigbase)
{
  zMatCopyArray( a, n, n, m );
  zMatMulDRC( m, s );
  if( !zMatSymEigJacobi( m, eigval, eigbase ) ){
    eprintf( "failed for a matrix scaled by %g\n", s );
    return false;
  }
  zMatCopyArray( a, n, n, m );
  zVecDivDRC( eigval, s );
  if( !check_mat_sym_eig_rel( m, eigval, eigbase, TOL_EIG ) ){
    eprintf( "inaccurate for a matrix scaled by %g\n", s );
    return false;
  }
  return true;
}

void assert_mat_sym_eig_jacobi_scale(void)
{
  double a[] = {
    4, 1,-2, 2,
    1, 2, 0, 1,
   -2, 0, 3,-2,
    2, 1,-2,-1,
  };
  const int n = 4;
  zMat m, eigbase;
  zVec eigval;
  double s;
  bool result = true;

  m = zMatAllocSqr( n );
  eigbase = zMatAllocSqr( n );
  eigval = zVecAlloc( n );
  for( s=1.0e-300; s<1.0e300; s*=1.0e20 )
    if( !assert_mat_sym_eig_jacobi_scale_one( a, n, s, m, eigval, eigbase ) ) result = false;
  /* the maximum component 1e308 is greater than 2^1023, while the eigenvalues do not overflow. */
  if( !assert_mat_sym_eig_jacobi_scale_one( a, n, 2.5e307, m, eigval, eigbase ) ) result = false;
  zMatFree( m );
  zMatFree( eigbase );
  zVecFree( eigval );
  zAssert( zMatSymEigJacobi (scaled matrix), result );
}

/* expected is sorted and freed in this function. */
bool assert_mat_sym_eig_jacobi_special_one(zMat m, zVec expected)
{
  zMat eigbase;
  zVec eigval;
  int i;
  bool result;

  eigbase = zMatAllocSqr( zMatRowSizeNC(m) );
  eigval = zVecAlloc( zMatRowSizeNC(m) );
  result = zMatSymEigJacobi( m, eigval, eigbase ) && check_mat_sym_eig_rel( m, eigval, eigbase, TOL_EIG );
  if( result && expected ){
    zDataSort( zVecBuf(eigval), zVecSizeNC(eigval) );
    zDataSort( zVecBuf(expected), zVecSizeNC(expected) );
    for( i=0; i<zVecSizeNC(eigval); i++ )
      if( !( fabs( zVecElemNC(eigval,i) - zVecElemNC(expected,i) ) <=
             TOL_EIG * zVecSizeNC(eigval) * DBL_EPSILON * zMax( zMatNorm(m), 1 ) ) ){
        eprintf( "eigenvalue #%d: %.17g / %.17g\n", i, zVecElemNC(eigval,i), zVecElemNC(expected,i) );
        result = false;
      }
  }
  zMatFree( m );
  zMatFree( eigbase );
  zVecFree( eigval );
  zVecFree( expected );
  return result;
}

void assert_mat_sym_eig_jacobi_special(void)
{
  zMat m;
  zVec e;
  int n, i, j;
  bool result = true;

  for( n=1; n<=20; n++ ){
    /* zero matrix */
    if( !assert_mat_sym_eig_jacobi_special_one( zMatAllocSqr( n ), zVecAlloc( n ) ) ){
      eprintf( "failed for the zero matrix (n=%d)\n", n );
      result = false;
    }
    /* matrix of ones, the eigenvalues of which are n and zeros */
    m = zMatAllocSqr( n );
    e = zVecAlloc( n );
    for( i=0; i<n; i++ )
      for( j=0; j<n; j++ ) zMatSetElemNC( m, i, j, 1 );
    zVecSetElemNC( e, 0, n );
    if( !assert_mat_sym_eig_jacobi_special_one( m, e ) ){
      eprintf( "failed for the matrix of ones (n=%d)\n", n );
      result = false;
    }
    /* tridiagonal matrix of 2 and -1, the eigenvalues of which are 2 - 2 cos(k pi/(n+1)) */
    m = zMatAllocSqr( n );
    e = zVecAlloc( n );
    for( i=0; i<n; i++ ){
      zMatSetElemNC( m, i, i, 2 );
      if( i > 0 ) zMatSetElemNC( m, i, i-1, -1 );
      if( i < n-1 ) zMatSetElemNC( m, i, i+1, -1 );
      zVecSetElemNC( e, i, 2 - 2*cos( (i+1)*zPI/(n+1) ) );
    }
    if( !assert_mat_sym_eig_jacobi_special_one( m, e ) ){
      eprintf( "failed for the tridiagonal matrix (n=%d)\n", n );
      result = false;
    }
    /* identity matrix with off-diagonal components much smaller than zTOL */
    m = zMatAllocSqr( n );
    e = zVecAlloc( n );
    for( i=0; i<n; i++ )
      for( j=0; j<n; j++ ) zMatSetElemNC( m, i, j, i == j ? 1 : 1.0e-20 );
    zVecSetAll( e, 1 - 1.0e-20 );
    zVecSetElemNC( e, 0, 1 + ( n - 1 ) * 1.0e-20 );
    if( !assert_mat_sym_eig_jacobi_special_one( m, e ) ){
      eprintf( "failed for the identity matrix with tiny perturbation (n=%d)\n", n );
      result = false;
    }
    /* components of very different magnitudes */
    m = zMatAllocSqr( n );
    for( i=0; i<n; i++ )
      for( j=0; j<n; j++ ) zMatSetElemNC( m, i, j, i == j ? ( i % 2 ? 1.0e100 : -1.0e-100 ) : 1.0e-75 );
    if( !assert_mat_sym_eig_jacobi_special_one( m, NULL ) ){
      eprintf( "failed for the matrix of very different components (n=%d)\n", n );
      result = false;
    }
  }
  zAssert( zMatSymEigJacobi (special matrix), result );
}

/* check if m = u diag{sv} v, u^T u = I and v v^T = I for the first rank columns of u and rows of v,
 * and the singular values are in descending order, with relative errors. */
bool check_mat_svd(zMat m, zMat u, zVec sv, zMat v, int rank, double tol)
{
  int r, c, i, j, l;
  double norm, d, e, f;
  bool result = true;

  r = zMatRowSizeNC(m);
  c = zMatColSizeNC(m);
  if( rank < 0 || rank > zMin( r, c ) || zMatRowSizeNC(v) != rank ){
    eprintf( "invalid rank %d or row size of v %d\n", rank, zMatRowSizeNC(v) );
    return false;
  }
  for( i=0; i<zVecSizeNC(sv); i++ )
    if( i < rank ? !( zVecElemNC(sv,i) > 0 ) || ( i > 0 && zVecElemNC(sv,i) > zVecElemNC(sv,i-1) ) :
                   zVecElemNC(sv,i) != 0 ){
      eprintf( "invalid singular value #%d = %g\n", i, zVecElemNC(sv,i) );
      result = false;
    }
  /* the residual is computed with the matrix scaled by the greatest singular value to avoid overflow */
  f = rank > 0 ? zVecElemNC(sv,0) : 1;
  for( norm=0, e=0, i=0; i<r; i++ )
    for( j=0; j<c; j++ ){
      norm += _zSqr( zMatElemNC(m,i,j) / f );
      for( d=zMatElemNC(m,i,j)/f, l=0; l<rank; l++ )
        d -= zMatElemNC(u,i,l) * ( zVecElemNC(sv,l) / f ) * zMatElemNC(v,l,j);
      e += d * d;
    }
  norm = sqrt( norm );
  if( !( sqrt( e ) <= tol * zMax( r, c ) * DBL_EPSILON * norm ) ){
    eprintf( "relative residual |m - u s v| / |m| = %g\n", sqrt( e ) / norm );
    result = false;
  }
  for( i=0; i<rank; i++ )
    for( j=0; j<rank; j++ ){
      for( d=( i == j ? -1 : 0 ), l=0; l<r; l++ )
        d += zMatElemNC(u,l,i) * zMatElemNC(u,l,j);
      if( !( fabs( d ) <= tol * zMax( r, c ) * DBL_EPSILON ) ){
        eprintf( "left singular vectors #%d and #%d not orthonormal: %g\n", i, j, d );
        result = false;
      }
      for( d=( i == j ? -1 : 0 ), l=0; l<c; l++ )
        d += zMatElemNC(v,i,l) * zMatElemNC(v,j,l);
      if( !( fabs( d ) <= tol * zMax( r, c ) * DBL_EPSILON ) ){
        eprintf( "right singular vectors #%d and #%d not orthonormal: %g\n", i, j, d );
        result = false;
      }
    }
  return result;
}

bool assert_mat_svd_one(double a[], int n, int m)
{
  zMat ma, u, v;
  zVec sv;
  int rank;
  bool result;

  ma = zMatCloneArray( a, n, m );
  u = zMatAllocSqr( n );
  v = zMatAlloc( n, m );
  sv = zVecAlloc( n );
  rank = zMatSVD( ma, u, sv, v );
  result = check_mat_svd( ma, u, sv, v, rank, TOL_EIG );
  zMatFreeAtOnce( 3, ma, u, v );
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

void assert_mat_svd_scale(void)
{
  double a[] = {
    2, 3,-2,
   -3,14,-7,
   -5,19,-9,
  };
  const int n = 3;
  zMat m, u, v;
  zVec sv;
  double s;
  bool result = true;

  m = zMatAllocSqr( n );
  u = zMatAllocSqr( n );
  v = zMatAllocSqr( n );
  sv = zVecAlloc( n );
  /* the maximum component 1.9e308 at the last is greater than 2^1023, while the singular values do
   * not overflow. */
  for( s=1.0e-300; s<1.0e300; s = s < 1.0e280 ? s*1.0e20 : 1.0e307 ){
    zMatCopyArray( a, n, n, m );
    zMatMulDRC( m, s );
    zMatSetRowSize( v, n );
    if( zMatSVD( m, u, sv, v ) != n ||
        !check_mat_svd( m, u, sv, v, n, TOL_EIG ) ){
      eprintf( "failed for a matrix scaled by %g\n", s );
      result = false;
    }
    if( s == 1.0e307 ) break;
  }
  zMatFree( m );
  zMatFree( u );
  zMatFree( v );
  zVecFree( sv );
  zAssert( zMatSVD (scaled matrix), result );
}

/* a random orthogonal matrix as the eigenvectors of a random symmetric matrix */
void rand_orthogonal(zMat q)
{
  zMat s;
  zVec e;
  int i, j;

  s = zMatAllocSqr( zMatRowSizeNC(q) );
  e = zVecAlloc( zMatRowSizeNC(q) );
  for( i=0; i<zMatRowSizeNC(q); i++ )
    for( j=0; j<=i; j++ )
      zMatSetElemNC( s, i, j, ( zMatElemNC(s,j,i) = zRandF(-1,1) ) );
  zMatSymEigJacobi( s, e, q );
  zMatFree( s );
  zVecFree( e );
}

/* matrices of known ranks */
void assert_mat_svd_rank(void)
{
  zMat m, u, v, u0, v0;
  zVec sv;
  int k, r, c, q, i, j, l, rank;
  double sigma[40], d, scale;
  bool result = true;

  for( k=0; k<N; k++ ){
    r = zRandI( 1, 20 );
    c = zRandI( 1, 20 );
    q = zRandI( 0, zMin( r, c ) );
    scale = pow( 10, zRandF( -200, 200 ) );
    m = zMatAlloc( r, c );
    u = zMatAllocSqr( r );
    v = zMatAlloc( r, c );
    sv = zVecAlloc( r );
    u0 = zMatAllocSqr( r );
    v0 = zMatAllocSqr( c );
    rand_orthogonal( u0 );
    rand_orthogonal( v0 );
    /* singular values from scale to scale x 1e-8, which are distinguished from zero */
    for( i=0; i<q; i++ ) sigma[i] = scale * pow( 10, -8.0*i/q );
    for( i=0; i<r; i++ )
      for( j=0; j<c; j++ ){
        for( d=0, l=0; l<q; l++ ) d += zMatElemNC(u0,i,l) * sigma[l] * zMatElemNC(v0,j,l);
        zMatSetElemNC( m, i, j, d );
      }
    if( ( rank = zMatSVD( m, u, sv, v ) ) != q ){
      eprintf( "rank of a %dx%d matrix = %d (expected %d)\n", r, c, rank, q );
      result = false;
    } else
    if( !check_mat_svd( m, u, sv, v, rank, TOL_EIG ) ){
      result = false;
    } else
      for( i=0; i<q; i++ )
        if( !( fabs( zVecElemNC(sv,i) - sigma[i] ) <= TOL_EIG * zMax( r, c ) * DBL_EPSILON * sigma[0] ) ){
          eprintf( "singular value #%d = %.17g (expected %.17g)\n", i, zVecElemNC(sv,i), sigma[i] );
          result = false;
        }
    zMatFreeAtOnce( 5, m, u, v, u0, v0 );
    zVecFree( sv );
  }
  for( r=1; r<=10; r++ )
    for( c=1; c<=10; c++ ){
      m = zMatAlloc( r, c );
      u = zMatAllocSqr( r );
      v = zMatAlloc( r, c );
      sv = zVecAlloc( r );
      /* zero matrix */
      if( ( rank = zMatSVD( m, u, sv, v ) ) != 0 || !check_mat_svd( m, u, sv, v, rank, TOL_EIG ) ){
        eprintf( "failed for the %dx%d zero matrix\n", r, c );
        result = false;
      }
      /* matrix of ones */
      for( i=0; i<r; i++ )
        for( j=0; j<c; j++ ) zMatSetElemNC( m, i, j, 1 );
      zMatSetRowSize( v, r );
      if( ( rank = zMatSVD( m, u, sv, v ) ) != 1 || !check_mat_svd( m, u, sv, v, rank, TOL_EIG ) ||
          !( fabs( zVecElemNC(sv,0) - sqrt(r*c) ) <= TOL_EIG * zMax( r, c ) * DBL_EPSILON * sqrt(r*c) ) ){
        eprintf( "failed for the %dx%d matrix of ones\n", r, c );
        result = false;
      }
      /* identity matrix multiplied by a small value */
      zMatZero( m );
      for( i=0; i<zMin(r,c); i++ ) zMatSetElemNC( m, i, i, 1.0e-3 );
      zMatSetRowSize( v, r );
      if( ( rank = zMatSVD( m, u, sv, v ) ) != zMin(r,c) || !check_mat_svd( m, u, sv, v, rank, TOL_EIG ) ){
        eprintf( "failed for the %dx%d identity matrix multiplied by 1e-3\n", r, c );
        result = false;
      }
      zMatFreeAtOnce( 3, m, u, v );
      zVecFree( sv );
    }
  zAssert( zMatSVD (rank-deficient matrix), result );
}

void assert_mat_svd_invalid(void)
{
  zMat m, u, v;
  zVec sv;
  bool result = true;

  m = zMatAlloc( 3, 4 );
  u = zMatAllocSqr( 3 );
  v = zMatAlloc( 3, 4 );
  sv = zVecAlloc( 4 );
  if( zMatSVD( m, u, sv, v ) != -1 ) result = false; /* size mismatch of sv */
  zVecFree( sv );
  sv = zVecAlloc( 3 );
  zMatFree( u );
  u = zMatAllocSqr( 4 );
  if( zMatSVD( m, u, sv, v ) != -1 ) result = false; /* size mismatch of u */
  zMatFreeAtOnce( 3, m, u, v );
  zVecFree( sv );
  zAssert( zMatSVD (invalid size), result );
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
  assert_mat_sym_eig_jacobi_scale();
  assert_mat_sym_eig_jacobi_special();
  assert_mat_svd();
  assert_mat_svd_random();
  assert_mat_svd_scale();
  assert_mat_svd_rank();
  assert_mat_svd_invalid();
  assert_mat_svd_minmax();
  return 0;
}
