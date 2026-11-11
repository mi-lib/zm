#include <zm/zm.h>

/* normwise backward error of an approximate solution x of a x = b in the infinity norm, namely,
 * |a x - b| / ( |a||x| + |b| ). It does not depend on the condition number of a nor the scale of the
 * problem, and is in the order of the machine epsilon for a backward stable solver. */
double le_backward_error(const zMat a, const zVec x, const zVec b, zVec r)
{
  zMulMatVec( a, x, r );
  zVecSubDRC( r, b );
  return zVecInfNorm( r ) / ( zMatInfNorm( a ) * zVecInfNorm( x ) + zVecInfNorm( b ) );
}

/* normwise backward error of an approximate solution x of m x = b in the infinity norm, where r is
 * a workspace with the same size of b. See le_backward_error() for vectors. */
double mat_backward_error(const zMat m, const zMat x, const zMat b, zMat r)
{
  zMulMatMat( m, x, r );
  zMatSubDRC( r, b );
  return zMatInfNorm( r ) / ( zMatInfNorm( m ) * zMatInfNorm( x ) + zMatInfNorm( b ) );
}

/* linear equation solver */

void assert_mat_mat_sweep_out(void)
{
  zMat m1, m2, m3, m4, m;
  zIndex index;
  bool result = true;
  const int size = 10, testnum = 1000;
  const double tol = 1.0e-8;
  int i, k;

  m  = zMatAllocSqr( size );
  m1 = zMatAllocSqr( size );
  m2 = zMatAllocSqr( zMatRowSizeNC(m1) );
  m3 = zMatAlloc( zMatRowSizeNC(m), zMatColSizeNC(m2) );
  m4 = zMatAllocSqr( zMatRowSizeNC(m2) );
  index = zIndexCreate( zMatRowSizeNC(m1) );
  for( k=0; k<testnum; k++ ){
    zMatRandUniform( m1, -5, 5 );
    zMatCopy( m1, m );
    zMatIdent( m2 );
    for( i=0; i<zMatRowSizeNC(m1); i++ ){
      zMatPivoting( m1, index, i, i );
      zMatMatSweepOut( m1, m2, zIndexElemNC(index,i), i );
    }
    for( i=0; i<zMatRowSizeNC(m2); i++ ) /* reorder rows along with the pivoting index */
      zRawVecCopy( zMatRowBufNC(m2,zIndexElemNC(index,i)), zMatRowBufNC(m4,i), zMatColSizeNC(m2) );
    zMulMatMat( m, m4, m3 );
    if( !zMatIsIdent( m3, tol ) ){
      eprintf( "maximum error = %g\n", zMatElemAbsMax( m3, NULL, NULL ) );
      result = false;
    }
  }
  zMatFreeAtOnce( 5, m, m1, m2, m3, m4 );
  zIndexFree( index );
  zAssert( zMatMatSweepOut, result );
}

void assert_mat_vec_sweep_out(void)
{
  zMat m;
  zVec v, x;
  zIndex index;
  bool result = true;
  const double tol = 1.0e-8;
  const int size = 10, testnum = 1000;
  int i, k;
  double xmax;

  m = zMatAllocSqr( size );
  v = zVecAlloc( size );
  x = zVecAlloc( size );
  index = zIndexCreate( size );
  for( k=0; k<testnum; k++ ){
    zMatRandUniform( m, -5, 5 );
    zVecRandUniform( v, -5, 5 );
    zLESolveGauss( m, v, x );
    for( i=0; i<zMatRowSize(m); i++ ){
      zMatPivoting( m, index, i, i );
      zMatVecSweepOut( m, v, zIndexElem(index,i), i );
    }
    zVecReorderDRC( v, index );
    xmax = zVecElemAbsMax( x, NULL );
    zVecSubDRC( v, x );
    if( !zVecIsTol( v, tol * zMax( xmax, 1.0 ) ) ){ /* relative error to the solution */
      eprintf( "maximum error = %g (maximum absolute component of solution = %g)\n", zVecElemAbsMax( v, NULL ), xmax );
      result = false;
    }
  }
  zMatFree( m );
  zVecFreeAtOnce( 2, v, x );
  zIndexFree( index );
  zAssert( zMatVecSweepOut, result );
}

void assert_le(void)
{
  zMat a;
  zVec b, x, ans, r;
  zMat l, u; /* for LU decomposition */
  zMat adj, inv; /* Cramer's method */
  zIndex index;
  const double tol = 100; /* tolerance of the error of Cramer's method relative to kappa(a) eps |x|, since it is not backward stable */
  const double tol_be = 1.0e-14; /* tolerance of the normwise backward error */
  const int size = 10;
  const int n = 1000;
  int i;
  int count_cramer, count_gauss, count_lu, count_ri;
  double kappa;

  a = zMatAllocSqr( size );
  adj = zMatAllocSqr( size );
  inv = zMatAllocSqr( size );
  l = zMatAllocSqr( size );
  u = zMatAllocSqr( size );
  b = zVecAlloc( size );
  x = zVecAlloc( size );
  ans = zVecAlloc( size );
  r = zVecAlloc( size );
  index = zIndexAlloc( size );
  if( !a || !adj || !inv || !l || !u || !b || !x || !ans || !r || !index ) goto TERMINATE;
  count_cramer = count_gauss = count_lu = count_ri = 0;
  for( i=0; i<n; i++ ){
    /* generate a problem */
    zMatRandUniform( a, -10, 10 );
    zVecRandUniform( ans,-10, 10 );
    zMulMatVec( a, ans, b );
    /* Cramer's method (adjoint matrix / determinant) */
    zMatAdj( a, adj );
    zMatDivDRC( adj, zMatDet(a) );
    zMulMatVec( adj, b, x );
    /* the forward error is evaluated relative to the condition number of a,
       which is computed with zMatInv() tested in assert_mat_inv() beforehand. */
    kappa = zMatInv( a, inv ) ? zMatInfNorm( a ) * zMatInfNorm( inv ) : HUGE_VAL;
    if( zVecInfNorm( zVecSubDRC( x, ans ) ) <= tol * kappa * DBL_EPSILON * zVecInfNorm( ans ) ) count_cramer++;
    else{
      eprintf( "error (Cramer's method) " );
      zVecFPrint( stderr, x );
    }
    /* Gauss's elimination method */
    zLESolveGauss( a, b, x );
    if( le_backward_error( a, x, b, r ) <= tol_be ) count_gauss++;
    else{
      eprintf( "error (Gauss's elimination method) " );
      zVecFPrint( stderr, x );
    }
    /* LU decomposition method */
    zMatDecompLU( a, l, u, index );
    zLESolveLU( l, u, b, x, index );
    if( le_backward_error( a, x, b, r ) <= tol_be ) count_lu++;
    else{
      eprintf( "error (LU decomposition method) " );
      zVecFPrint( stderr, x );
    }
    /* LU decomposition + residual iteration method */
    zLESolveLURI( a, b, x );
    if( le_backward_error( a, x, b, r ) <= tol_be ) count_ri++;
    else{
      eprintf( "error (LU decomposition + residual iteration method) " );
      zVecFPrint( stderr, x );
    }
  }
 TERMINATE:
  zMatFreeAtOnce( 5, a, l, u, adj, inv );
  zVecFreeAtOnce( 4, b, x, ans, r );
  zIndexFree( index );

  eprintf( "success rate (Cramer) %d/%d ", count_cramer, n ); zAssert( zMatAdj + zMatDet, count_cramer == n );
  eprintf( "success rate (Gauss)  %d/%d ", count_gauss, n );  zAssert( zLESolveGauss,     count_gauss == n );
  eprintf( "success rate (LU)     %d/%d ", count_lu, n );     zAssert( zLESolveLU,        count_lu == n );
  eprintf( "success rate (LU+RI)  %d/%d ", count_ri, n );     zAssert( zLESolveLURI,      count_ri == n );
}

void assert_le_mat_capacity(void)
{
  zMat a;
  zVec b, x, ans, r;
  zIndex index;
  const double tol = 1.0e-14; /* tolerance of the normwise backward error */
  const int rowcapacity = 15, colcapacity = 13;
  const int size = 10;
  const int n = 1000;
  int i;
  int count_gauss = 0;

  a = zMatAlloc( rowcapacity, colcapacity );
  b = zVecAlloc( colcapacity );
  x = zVecAlloc( rowcapacity );
  ans = zVecAlloc( size );
  r = zVecAlloc( size );
  index = zIndexAlloc( size*2 );
  if( !a || !b || !x || !ans || !r || !index ) goto TERMINATE;
  zMatSetSize( a, size, size );
  zVecSetSize( b, size );
  zVecSetSize( x, size );
  zIndexSetSize( index, size );
  for( i=0; i<n; i++ ){
    /* generate a problem */
    zMatRandUniform( a, -10, 10 );
    zVecRandUniform( ans,-10, 10 );
    zMulMatVec( a, ans, b );
    /* Gauss's elimination method */
    zLESolveGauss( a, b, x );
    if( le_backward_error( a, x, b, r ) <= tol ) count_gauss++;
    else{
      eprintf( "error (Gauss's elimination method) " );
      zVecFPrint( stderr, x );
    }
  }
 TERMINATE:
  zMatFree( a );
  zVecFreeAtOnce( 4, b, x, ans, r );
  zIndexFree( index );
  zAssert( zLESolveGauss (more-capacity case), count_gauss == n );
}

void assert_le_gauss_seidel(void)
{
  zMat a, a_org;
  zVec b, x, ans;
  const double tol = ( 1.0e-9 );
  const int size = 10;
  const int n = 1000;
  int i, j;
  int count_gs = 0;

  a = zMatAllocSqr( size );
  a_org = zMatAllocSqr( size );
  b = zVecAlloc( size );
  x = zVecAlloc( size );
  ans = zVecAlloc( size );
  if( !a || !a_org || !b || !x || !ans ) goto TERMINATE;
  for( i=0; i<n; i++ ){
    /* generate a problem */
    zMatRandUniform( a_org, -10, 10 );
    zMulMatTMat( a_org, a_org, a );
    for( j=0; j<size; j++ ) zMatElemNC(a,j,j) += 200;
    zVecRandUniform( ans,-10, 10 );
    zMulMatVec( a, ans, b );
    /* Gauss-Seidel's method */
    zVecZero( x );
    zLESolveGaussSeidel( a, b, x );
    if( zVecIsTol( zVecSubDRC( x, ans ), tol ) ) count_gs++;
    else{
      eprintf( "error (Gauss-Seidel's method) " );
      zVecFPrint( stderr, x );
    }
  }
 TERMINATE:
  zMatFreeAtOnce( 2, a, a_org );
  zVecFreeAtOnce( 3, b, x, ans );
  eprintf( "success rate (Gauss-Seidel) %d/%d ", count_gs, n );
  zAssert( zLESolveGaussSeidel, count_gs == n );
}

void assert_le_gauss_seidel_scale(void)
{
  const double scale[] = { 1.0e-12, 1.0e-6, 1.0, 1.0e6, 1.0e12 };
  const double tol = 1.0e-10; /* tolerance of the relative error */
  const int size = 10, testnum = 100;
  zMat a, a_org;
  zVec b, x, ans;
  int i, j, k;
  double err;
  bool result = true;

  a = zMatAllocSqr( size );
  a_org = zMatAllocSqr( size );
  b = zVecAlloc( size );
  x = zVecAlloc( size );
  ans = zVecAlloc( size );
  /* the convergence has to be judged relatively to the scale of the solution */
  for( i=0; i<sizeof(scale)/sizeof(double); i++ )
    for( k=0; k<testnum; k++ ){
      zMatRandUniform( a_org, -10, 10 );
      zMulMatTMat( a_org, a_org, a );
      for( j=0; j<size; j++ ) zMatElemNC(a,j,j) += 200;
      zVecRandUniform( ans, -10, 10 );
      zVecMulDRC( ans, scale[i] );
      zMulMatVec( a, ans, b );
      zVecZero( x );
      zLESolveGaussSeidel( a, b, x );
      zVecSubDRC( x, ans );
      if( ( err = zVecInfNorm( x ) / zVecInfNorm( ans ) ) > tol ){
        eprintf( "scale %g, case #%d: relative error = %g\n", scale[i], k, err );
        result = false;
      }
    }
  zMatFreeAtOnce( 2, a, a_org );
  zVecFreeAtOnce( 3, b, x, ans );
  zAssert( zLESolveGaussSeidel (scaled solutions), result );
}

/* generalized linear equation solver (benchmarking) */

void generate_equation_general(zMat a, zVec b, zVec w, zVec w2, zVec x, zVec _b)
{
  zMatRandUniform( a, -10, 10 );
  zVecRandUniform( b, -10, 10 );
  zVecSetAll( w, 1.0 );
  zVecSetAll( w2, 1000.0 );
}

void generate_equation_illposed(zMat a, zVec b, zVec w, zVec w2, zVec x, zVec _b)
{
  int i, j, n;
  double e;

  for( j=0; j<zMatColSizeNC(a); j++ )
    zMatSetElem( a, 0, j, zRandF(-10,10) );
  zVecSetElem( b, 0, zRandF(-10,10) );
  zVecSetElem( w2, 0, 1.0e5 );
  for( i=1; i<zMatRowSizeNC(a); i++ ){
    zRawVecCopy( zMatRowBufNC(a,0), zMatRowBufNC(a,i), zMatColSizeNC(a) );
    zVecSetElem( b, i, zVecElem(b,0) );
    zVecSetElem( w2, i, zVecElem(w2,0) );
  }
  n = 0.8 * zMatMinSize( a );
  for( i=0; i<n; i++ ){
    e = zRandF(-0.01,0.01);
    zMatElemNC(a,i,i) += e;
    zVecElemNC(b,i  ) += e;
  }
  zVecSetAll( w, 1.0 );
}

#define NUM_METHOD 4

void try_le_gen_one(zMat a, zVec b, zVec w, zVec w2, zVec x, zVec _b, double error[NUM_METHOD], long delta_clock[NUM_METHOD])
{
  clock_t c1, c2;

  c1 = clock();
  zLESolveMP_LQ( a, b, w, w2, x );
  c2 = clock();
  delta_clock[0] = c2 - c1;
  error[0] = le_backward_error( a, x, b, _b );
  eprintf( "MP(LQ):%g %ld ", error[0], delta_clock[0] );

  c1 = clock();
  zLESolveMP_LU( a, b, w, w2, x );
  c2 = clock();
  delta_clock[1] = c2 - c1;
  error[1] = le_backward_error( a, x, b, _b );
  eprintf( "MP(LU):%g %ld ", error[1], delta_clock[1] );

  c1 = clock();
  zLESolveMP_SVD( a, b, x );
  c2 = clock();
  delta_clock[2] = c2 - c1;
  error[2] = le_backward_error( a, x, b, _b );
  eprintf( "MP(SVD):%g %ld ", error[2], delta_clock[2] );

  c1 = clock();
  zLESolveSR( a, b, w, w2, x );
  c2 = clock();
  delta_clock[3] = c2 - c1;
  error[3] = le_backward_error( a, x, b, _b );
  eprintf( "SR:%g %ld ", error[3], delta_clock[3] );

  eprintf( "\n" );
}

#define NUM_TRIAL 10

void eval_le_gen(const double error[NUM_TRIAL][NUM_METHOD], const double tol[NUM_METHOD], bool result[NUM_METHOD])
{
  int i, j;

  for( j=0; j<NUM_METHOD; j++ ) result[j] = true;
  for( i=0; i<NUM_TRIAL; i++ ){
    for( j=0; j<NUM_METHOD; j++ ){
      if( !zIsTol( error[i][j], tol[j] ) ) result[j] = false;
    }
  }
}

void assert_le_gen(void)
{
  /* tolerances of the normwise backward errors for MP(LQ), MP(LU), MP(SVD), and SR. Those for MP(SVD)
   * in ill-conditioned cases and SR are loose, since the former is based on the eigenvalues of a^T a,
   * and the latter does not intend to solve the equation exactly. In ill-conditioned cases, a row
   * which is close to linearly dependent on the others within ZM_LQ_DECOMP_TOL is discarded by LQ
   * decomposition, and MP(LU) loses accuracy in solving the normal equations of the factors. */
  const double tol_general[] = { 1.0e-14, 1.0e-12, 1.0e-13, 1.0e-6 };
  const double tol_illcond[] = { ZM_LQ_DECOMP_TOL, 1.0e-11, 1.0e-5, 1.0e-5 };
  const int rowsize = 150;
  const int colsize = 200;
  zMat a;
  zVec b, w, x, _b, w2;
  int trial;
  double error_general[NUM_TRIAL][NUM_METHOD];
  double error_illcond[NUM_TRIAL][NUM_METHOD];
  long delta_clock_general[NUM_TRIAL][NUM_METHOD];
  long delta_clock_illcond[NUM_TRIAL][NUM_METHOD];
  bool result_general[NUM_METHOD];
  bool result_illcond[NUM_METHOD];

  a = zMatAlloc( rowsize, colsize );
  b = zVecAlloc( rowsize );
  _b= zVecAlloc( rowsize );
  w = zVecAlloc( colsize );
  w2= zVecAlloc( rowsize );
  x = zVecAlloc( colsize );

  eprintf( "(general test cases)\n" );
  for( trial=0; trial<NUM_TRIAL; trial++ ){
    generate_equation_general( a, b, w, w2, x, _b );
    try_le_gen_one( a, b, w, w2, x, _b, error_general[trial], delta_clock_general[trial] );
  }
  eval_le_gen( error_general, tol_general, result_general );
  eprintf( "(ill-conditioned test cases)\n" );
  for( trial=0; trial<NUM_TRIAL; trial++ ){
    generate_equation_illposed( a, b, w, w2, x, _b );
    try_le_gen_one( a, b, w, w2, x, _b, error_illcond[trial], delta_clock_illcond[trial] );
  }
  eval_le_gen( error_illcond, tol_illcond, result_illcond );

  zMatFree( a );
  zVecFreeAtOnce( 5, b, _b, w, w2, x );
  zAssert( zLESolveMP_LQ, result_general[0] );
  zAssert( zLESolveMP_LU, result_general[1] );
  zAssert( zLESolveMP_SVD, result_general[2] );
  zAssert( zLESolveSR, result_general[3] );
  zAssert( zLESolveMP_LQ (ill-conditioned cases), result_illcond[0] );
  zAssert( zLESolveMP_LU (ill-conditioned cases), result_illcond[1] );
  zAssert( zLESolveMP_SVD (ill-conditioned cases), result_illcond[2] );
  zAssert( zLESolveSR (ill-conditioned cases), result_illcond[3] );
}

/* minimum-norm least-square solution of (f1 f2) x = b, where f1 and f2 are column-full-rank and row-full-rank. */
zVec le_gen_rowlong_ref(const zMat f1, const zMat f2, const zVec b, zVec x)
{
  zMat g;
  zVec y, z;

  g = zMatAllocSqr( zMatColSizeNC(f1) );
  y = zVecAlloc( zMatColSizeNC(f1) );
  z = zVecAlloc( zMatColSizeNC(f1) );
  zMulMatTMat( f1, f1, g );
  zMulMatTVec( f1, b, y );
  zLESolveGauss( g, y, z ); /* z = ( f1^T f1 )^-1 f1^T b */
  zMulMatMatT( f2, f2, g );
  zLESolveGauss( g, z, y ); /* y = ( f2 f2^T )^-1 z */
  zMulMatTVec( f2, y, x );
  zMatFree( g );
  zVecFreeAtOnce( 2, y, z );
  return x;
}

/* generalized linear equation solver for row-long (overdetermined) matrices. */
void assert_le_gen_rowlong(void)
{
  const int rowsize = 8, colsize = 5, testnum = 1000;
  const double tol = 1.0e-8;
  zMat a, f1, f2;
  zVec b, x_ref, x_lq;
  int i, rank;
  bool result = true;

  a = zMatAlloc( rowsize, colsize );
  b = zVecAlloc( rowsize );
  x_ref = zVecAlloc( colsize );
  x_lq = zVecAlloc( colsize );
  for( i=0; i<testnum; i++ ){
    rank = zRandI( 1, colsize );
    f1 = zMatAlloc( rowsize, rank );
    f2 = zMatAlloc( rank, colsize );
    zMatRandUniform( f1, -1, 1 );
    zMatRandUniform( f2, -1, 1 );
    zMulMatMat( f1, f2, a );
    zVecRandUniform( b, -1, 1 );
    le_gen_rowlong_ref( f1, f2, b, x_ref );
    zMatFreeAtOnce( 2, f1, f2 );
    zVecZero( x_lq );
    zLESolveMP_LQ( a, b, NULL, NULL, x_lq );
    if( !zVecEqual( x_lq, x_ref, tol * zMax( zVecElemAbsMax( x_ref, NULL ), 1.0 ) ) ) /* relative error */
      result = false;
  }
  zMatFree( a );
  zVecFreeAtOnce( 3, b, x_ref, x_lq );
  zAssert( zLESolveMP_LQ (row-long cases), result );
}

void assert_le_gen_aux(void)
{
  zMat a;
  zVec b, w, x1, x2, e1, e2, w2, aux;
  clock_t c1, c2, c3;
  long delta_clock_mp, delta_clock_mp_aux, delta_clock_sr, delta_clock_sr_aux;
  double e_max_mp, e_max_mp_aux, e_max_sr, e_max_sr_aux;
  const int rowsize = 150;
  const int colsize = 200;
  const double tol_sr = 1.0e-6;
  const double tol_mp = 1.0e-10;

  a = zMatAlloc( rowsize, colsize );
  b = zVecAlloc( rowsize );
  e1 = zVecAlloc( rowsize );
  e2 = zVecAlloc( rowsize );
  w = zVecAlloc( colsize );
  w2= zVecAlloc( rowsize );
  x1 = zVecAlloc( colsize );
  x2 = zVecAlloc( colsize );
  aux = zVecAlloc( colsize );

  zMatRandUniform( a, -10, 10 );
  zVecRandUniform( b, -10, 10 );
  zVecSetAll( w2, 1.0e6 );
  zVecSetAll( w, 1.0 );
  zVecRandUniform( aux, -10, 10 );

  c1 = clock();
  zLESolveMP( a, b, w, w2, x1 );
  c2 = clock();
  zLESolveMPAux( a, b, w, w2, x2, aux );
  c3 = clock();
  delta_clock_mp = c2 - c1;
  delta_clock_mp_aux = c3 - c2;
  zMulMatVec( a, x1, e1 ); zVecSubDRC( e1, b );
  zMulMatVec( a, x2, e2 ); zVecSubDRC( e2, b );
  e_max_mp = zVecElemAbsMax( e1, NULL );
  e_max_mp_aux = zVecElemAbsMax( e2, NULL );

  c1 = clock();
  zLESolveSR( a, b, w, w2, x1 );
  c2 = clock();
  zLESolveSRAux( a, b, w, w2, x2, aux );
  c3 = clock();
  delta_clock_sr = c2 - c1;
  delta_clock_sr_aux = c3 - c2;
  zMulMatVec( a, x1, e1 ); zVecSubDRC( e1, b );
  zMulMatVec( a, x2, e2 ); zVecSubDRC( e2, b );
  e_max_sr = zVecElemAbsMax( e1, NULL );
  e_max_sr_aux = zVecElemAbsMax( e2, NULL );

  zMatFree( a );
  zVecFreeAtOnce( 8, b, w, w2, x1, x2, e1, e2, aux );
  eprintf( "MP:%g %ld %g %ld ", e_max_mp, delta_clock_mp, e_max_mp_aux, delta_clock_mp_aux );
  zAssert( zLESolveMP / zLESolveMPAux, zIsTol( e_max_mp, tol_mp ) && zIsTol( e_max_mp_aux, tol_mp ) );

  eprintf( "SR:%g %ld %g %ld ", e_max_sr, delta_clock_sr, e_max_sr_aux, delta_clock_sr_aux );
  zAssert( zLESolveSR / zLESolveSRAux, zIsTol( e_max_sr, tol_sr ) && zIsTol( e_max_sr_aux, tol_sr ) );
}

void assert_le_mpnull(void)
{
  zMat a, mn;
  zVec b, w, w2, x1, x2, aux, dx, e1, e2;
  double e_max_mp, e_max_mp_aux;
  const int rowsize = 150;
  const int colsize = 200;
  const double tol_mp = 1.0e-10;

  a = zMatAlloc( rowsize, colsize );
  mn = zMatAllocSqr( colsize );
  b = zVecAlloc( rowsize );
  w = zVecAlloc( colsize );
  w2= zVecAlloc( rowsize );
  x1 = zVecAlloc( colsize );
  x2 = zVecAlloc( colsize );
  aux = zVecAlloc( colsize );
  dx = zVecAlloc( colsize );
  e1 = zVecAlloc( rowsize );
  e2 = zVecAlloc( rowsize );

  zMatRandUniform( a, -10, 10 );
  zVecRandUniform( b, -10, 10 );
  zVecSetAll( w2, 1.0e6 );
  zVecSetAll( w, 1.0 );
  zVecRandUniform( aux, -10, 10 );

  zLESolveMP( a, b, w, w2, x1 );
  zLESolveMPNull( a, b, w, w2, x2, mn );

  zMulMatVec( a, x1, e1 ); zVecSubDRC( e1, b );
  zMulMatVec( mn, aux, dx );
  zVecAddDRC( x2, dx );
  zMulMatVec( a, x2, e2 ); zVecSubDRC( e2, b );

  e_max_mp = zVecElemAbsMax( e1, NULL );
  e_max_mp_aux = zVecElemAbsMax( e2, NULL );

  zMatFreeAtOnce( 2, a, mn );
  zVecFreeAtOnce( 9, b, w, w2, x1, x2, e1, e2, aux, dx );
  zAssert( zLESolveMPNull, zIsTol( e_max_mp, tol_mp ) && zIsTol( e_max_mp_aux, tol_mp ) );
}

/* special linear equation solver */

/* check if the minimum-weighted-norm solution x of a x = b referring ref satisfies the optimality,
 * namely, if w^-1 (x - ref) is in the row space of a. q is a row-orthonormal matrix which spans the
 * row space of a, and d, p and v are workspaces. It returns the relative error of the projection. */
double le_ref_norm_min_optimality(const zMat q, const zVec x, const zVec ref, const zVec w, zVec d, zVec p, zVec v)
{
  zVecSub( x, ref, d );
  zVecDemDRC( d, w );
  zMulMatVec( q, d, v );
  zMulMatTVec( q, v, p );
  zVecSubDRC( p, d );
  return zVecInfNorm( p ) / zVecInfNorm( d );
}

void assert_le_ref_norm_min(void)
{
  const int rowsize = 30;
  const int colsize = 50;
  const int testnum = 1000;
  const double tol_err = 1.0e-14; /* tolerance of the normwise backward error */
  const double tol_opt = 1.0e-13; /* tolerance of the relative error of the optimality */
  zMat a, l, q;
  zVec b, ref, w, x, r, d, p, v;
  double err, opt, scale;
  int i, k;
  bool result = true;

  a = zMatAlloc( rowsize, colsize );
  l = zMatAllocSqr( rowsize );
  q = zMatAlloc( rowsize, colsize );
  b = zVecAlloc( rowsize );
  ref = zVecAlloc( colsize );
  w = zVecAlloc( colsize );
  x = zVecAlloc( colsize );
  r = zVecAlloc( rowsize );
  d = zVecAlloc( colsize );
  p = zVecAlloc( colsize );
  v = zVecAlloc( rowsize );
  for( k=0; k<testnum; k++ ){
    zMatRandUniform( a, -10, 10 );
    zVecRandUniform( b, -10, 10 );
    zVecRandUniform( ref, -10, 10 );
    zVecRandUniform( w, 0.1, 10 );
    if( k % 2 == 1 ) /* unbalanced scales of rows */
      for( i=0; i<rowsize; i++ ){
        scale = pow( 10, zRandF( -4, 4 ) );
        zRawVecMulDRC( zMatRowBufNC(a,i), scale, colsize );
        zVecElemNC(b,i) *= scale;
      }
    zLESolveRefNormMin( a, b, w, ref, x );
    zMatDecompLQ( a, l, q ); /* a is row-full-rank */
    err = le_backward_error( a, x, b, r );
    opt = le_ref_norm_min_optimality( q, x, ref, w, d, p, v );
    if( err > tol_err || opt > tol_opt ){
      eprintf( "case #%d: backward error = %g, optimality error = %g\n", k, err, opt );
      result = false;
    }
  }
  zMatFreeAtOnce( 3, a, l, q );
  zVecFreeAtOnce( 8, b, ref, w, x, r, d, p, v );
  zAssert( zLESolveRefNormMin, result );
}

/* normwise backward error of an approximate solution x of a tridiagonal equation */
double tridiag_backward_error(const zVec a, const zVec b, const zVec c, const zVec d, const zVec x)
{
  int i, n;
  double r, s, r_max = 0, m_norm = 0;

  n = zVecSizeNC(x);
  for( i=0; i<n; i++ ){
    r = zVecElemNC(b,i)*zVecElemNC(x,i) - zVecElemNC(d,i);
    s = fabs( zVecElemNC(b,i) );
    if( i > 0 ){
      r += zVecElemNC(a,i)*zVecElemNC(x,i-1);
      s += fabs( zVecElemNC(a,i) );
    }
    if( i < n-1 ){
      r += zVecElemNC(c,i)*zVecElemNC(x,i+1);
      s += fabs( zVecElemNC(c,i) );
    }
    r_max = zMax( r_max, fabs( r ) );
    m_norm = zMax( m_norm, s );
  }
  return r_max / ( m_norm * zVecInfNorm( x ) + zVecInfNorm( d ) );
}

void assert_tridiagonal_equation(void)
{
  zVec a, b, c, d, ans;
  int k;
  double err;
  const int dim = 100, testnum = 1000;
  bool result = true;
  const double tol = 1.0e-14; /* tolerance of the normwise backward error */

  a = zVecAlloc( dim );
  b = zVecAlloc( dim );
  c = zVecAlloc( dim );
  d = zVecAlloc( dim );
  ans = zVecAlloc( dim );
  /* random matrices, which are not necessarily diagonally dominant */
  for( k=0; k<testnum; k++ ){
    zVecRandUniform( a, -10, 10 );
    zVecRandUniform( b, -10, 10 );
    zVecRandUniform( c, -10, 10 );
    zVecRandUniform( d, -10, 10 );
    err = zLETridiagSolve( a, b, c, d, ans ) ? tridiag_backward_error( a, b, c, d, ans ) : HUGE_VAL;
    if( err > tol ){
      eprintf( "case #%d: backward error = %g\n", k, err );
      result = false;
    }
  }
  /* a matrix whose diagonal components are all zero, which requires pivoting */
  zVecSetAll( a, 1 );
  zVecZero( b );
  zVecSetAll( c, 1 );
  zVecRandUniform( d, -10, 10 );
  err = zLETridiagSolve( a, b, c, d, ans ) ? tridiag_backward_error( a, b, c, d, ans ) : HUGE_VAL;
  if( err > tol ){
    eprintf( "zero diagonal: backward error = %g\n", err );
    result = false;
  }
  zVecFreeAtOnce( 5, a, b, c, d, ans );
  zAssert( zLETridiagSolve, result );
}

void assert_lyapnov_equation(void)
{
  zMat a, b, x, tmp1, tmp2;
  const int size = 10, testnum = 100;
  const double tol = 1.0e-14; /* tolerance of the normwise backward error */
  int k;
  double err;
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
    /* normwise backward error of x a + a^T x = b */
    if( ( err = zMatInfNorm( tmp1 ) / ( 2 * zMatInfNorm( a ) * zMatInfNorm( x ) + zMatInfNorm( b ) ) ) > tol ){
      eprintf( "case #%d: backward error = %g\n", k, err );
      result = false;
    }
  }
  zMatFreeAtOnce( 5, a, b, x, tmp1, tmp2 );
  zAssert( zLELyapnovSolve, result );
}

void assert_lyapnov_equation_capacity(void)
{
  zMat a, b, x, xc;
  const int size = 10, testnum = 100;
  int i, k;
  bool result = true;

  a = zMatAllocSqr( size );
  b = zMatAllocSqr( size );
  x = zMatAllocSqr( size );
  xc = zMatAlloc( size, size*2 );
  zMatSetColSize( xc, size ); /* the column capacity is larger than the column size */
  for( k=0; k<testnum; k++ ){
    zMatRandUniform( a, -10, 10 );
    zMatRandUniform( b, -10, 10 );
    zLELyapnovSolve( a, b, x );
    zLELyapnovSolve( a, b, xc );
    for( i=0; i<size; i++ )
      if( memcmp( zMatRowBufNC(x,i), zMatRowBufNC(xc,i), sizeof(double)*size ) != 0 ){
        eprintf( "case #%d: the solution differs at the %dth row\n", k, i );
        result = false;
        break;
      }
  }
  zMatFreeAtOnce( 4, a, b, x, xc );
  zAssert( zLELyapnovSolve (more-capacity case), result );
}

/* matrix inversion */

void assert_mat_inv(void)
{
  zMat m1, m2, m;
  const int size = 10, testnum = 100;
  const double tol = 1.0e-14; /* tolerance of the normwise backward error */
  bool result = true;
  int i, k;
  double err;

  m1 = zMatAllocSqr( size );
  m2 = zMatAllocSqr( size );
  m  = zMatAllocSqr( size );
  for( k=0; k<testnum; k++ ){
    zMatRandUniform( m1, -10, 10 );
    zMatInv( m1, m2 );
    zMulMatMat( m1, m2, m );
    for( i=0; i<size; i++ ) zMatElemNC(m,i,i) -= 1.0;
    /* normwise backward error of m1 m2 = I */
    if( ( err = zMatInfNorm( m ) / ( zMatInfNorm( m1 ) * zMatInfNorm( m2 ) ) ) > tol ){
      eprintf( "case #%d: backward error = %g\n", k, err );
      result = false;
    }
  }
  zMatFreeAtOnce( 3, m1, m2, m );
  zAssert( zMatInv, result );
}

void assert_mat_mul_inv(void)
{
  zMat m1, m2, m3, m, r;
  const int rowsize = 10;
  const int colsize = 8;
  const int testnum = 100;
  bool result = true;
  int k;
  const double tol = 1.0e-14; /* tolerance of the normwise backward error */
  double err;

  m1 = zMatAlloc( rowsize, colsize );
  m2 = zMatAlloc( rowsize, colsize );
  m3 = zMatAlloc( rowsize, colsize );
  r = zMatAlloc( rowsize, colsize );
  m = zMatAllocSqr( rowsize );
  for( k=0; k<testnum; k++ ){
    zMatRandUniform( m, -10, 10 );
    zMatRandUniform( m1,-10, 10 );
    zMulMatMat( m, m1, m2 );
    zMulInvMatMat( m, m2, m3 );
    if( ( err = mat_backward_error( m, m3, m2, r ) ) > tol ){
      eprintf( "case #%d: backward error = %g\n", k, err );
      result = false;
    }
  }
  zMatFreeAtOnce( 5, m, m1, m2, m3, r );
  zAssert( zMulInvMatMat, result );
}

void assert_mat_mul_inv_capacity(void)
{
  zMat m1, m2, m3, m, r;
  const int colcapacity = 20;
  const int rowsize = 10;
  const int colsize = 8;
  const int testnum = 100;
  bool result = true;
  int k;
  const double tol = 1.0e-14; /* tolerance of the normwise backward error */
  double err;

  m1 = zMatAlloc( rowsize, colcapacity );
  m2 = zMatAlloc( rowsize, colsize );
  m3 = zMatAlloc( rowsize, colsize );
  r = zMatAlloc( rowsize, colsize );
  m = zMatAlloc( rowsize, colcapacity );
  zMatSetColSize( m1, colsize );
  zMatSetColSize( m, rowsize );
  for( k=0; k<testnum; k++ ){
    zMatRandUniform( m, -10, 10 );
    zMatRandUniform( m1,-10, 10 );
    zMulMatMat( m, m1, m2 );
    zMulInvMatMat( m, m2, m3 );
    if( ( err = mat_backward_error( m, m3, m2, r ) ) > tol ){
      eprintf( "case #%d: backward error = %g\n", k, err );
      result = false;
    }
  }
  zMatFreeAtOnce( 5, m, m1, m2, m3, r );
  zAssert( zMulInvMatMat (more-capacity case), result );
}

void assert_mat_mul_inv_mat_vec(void)
{
  zMat m1, m2, mat_ans, mat_org, mat_r;
  zVec v, vec_ans, vec_org, vec_r;
  const int size = 10;
  const int n = 100;
  int i;
  const double tol = 1.0e-14; /* tolerance of the normwise backward error */
  double err;
  bool result = true;

  m1 = zMatAllocSqr( size );
  m2 = zMatAllocSqr( size );
  mat_ans = zMatAllocSqr( size );
  mat_org = zMatAllocSqr( size );
  mat_r = zMatAllocSqr( size );
  v = zVecAlloc( size );
  vec_ans = zVecAlloc( size );
  vec_org = zVecAlloc( size );
  vec_r = zVecAlloc( size );
  if( !m1 || !m2 || !mat_ans || !mat_org || !mat_r || !v || !vec_ans || !vec_org || !vec_r ){
    ZALLOCERROR();
    goto TERMINATE;
  }
  for( i=0; i<n; i++ ){
    zMatRandUniform( mat_org, -10, 10 );
    zMatRandUniform( m1, -10, 10 );
    zMulMatMat( m1, mat_org, m2 );
    zVecRandUniform( vec_org, -10, 10 );
    zMulMatVec( m1, vec_org, v );
    zMulInvMatMatAndVec( m1, m2, v, mat_ans, vec_ans );
    if( ( err = le_backward_error( m1, vec_ans, v, vec_r ) ) > tol ){
      eprintf( "vector: backward error = %g\n", err );
      result = false;
    }
    if( ( err = mat_backward_error( m1, mat_ans, m2, mat_r ) ) > tol ){
      eprintf( "matrix: backward error = %g\n", err );
      result = false;
    }
  }

 TERMINATE:
  zMatFreeAtOnce( 5, m1, m2, mat_ans, mat_org, mat_r );
  zVecFreeAtOnce( 4, v, vec_ans, vec_org, vec_r );
  zAssert( zMulInvMatMatAndVec, result );
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
  zVec v;
  int i;
  bool result;
  const int rowsize = 40;
  const int colsize = 80;
  const int rank = 20;
  const int n = 10;
  const double tol = 1.0e-10;

  m1 = zMatAlloc( rowsize, colsize );
  m2 = zMatAlloc( rowsize, rank );
  ma = zMatAlloc( colsize, rank );
  mb = zMatAlloc( colsize, rank );
  mp = zMatAlloc( colsize, rowsize );
  v  = zVecAlloc( rowsize );
  result = true;
  for( i=0; i<n; i++ ){
    zMatRandUniform( m1, -10, 10 );
    mat_deg( m1, rank );
    zMatRandUniform( m2, -10, 10 );
    zMatGetCol( m2, 0, v );

    zMatMPInv( m1, mp );
    zMulMatMat( mp, m2, ma );
    zMulMPInvMatMat( m1, m2, mb );
    if( !zMatEqual(ma,mb,tol) ) result = false;
  }
  zAssert( zMatMPInv + zMulMPInvMatMat (random test), result );

  zMatFreeAtOnce( 5, m1, m2, ma, mb, mp );
  zVecFree( v );
}

void assert_mpnull(void)
{
  zMat a, mp, mn;
  zVec v, u, e;
  const int rowsize = 40;
  const int colsize = 80;
  const int rank = 20;
  const double tol = 1.0e-10;

  a = zMatAlloc( rowsize, colsize );
  zMatRandUniform( a, -10, 10 );
  mat_deg( a, rank );
  mp = zMatAlloc( colsize, rowsize );
  mn = zMatAlloc( colsize, colsize );
  zMatMPInvNull( a, mp, mn );

  v = zVecAlloc( colsize );
  u = zVecAlloc( colsize );
  e = zVecAlloc( rowsize );
  zVecRandUniform( v, -10, 10 );
  zMulMatVec( mn, v, u );
  zMulMatVec( a, u, e );
  zAssert( zMatMPInvNull, zVecIsTol(e,tol) );
  zMatFreeAtOnce( 3, a, mp, mn );
  zVecFreeAtOnce( 3, v, u, e );
}

void assert_mat_det_adj(void)
{
  zMat m, minv1, minv2;
  double det;
  const int size = 10, testnum = 100;
  int k;
  const double tol = 1.0e-10;
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
    if( !zMatEqual( minv1, minv2, tol ) ){
      zMatSubDRC( minv1, minv2 );
      eprintf( "case #%d: maximum error = %g\n", k, zMatElemAbsMax( minv1, NULL, NULL ) );
      result = false;
    }
  }
  zMatFreeAtOnce( 3, m, minv1, minv2 );
  zAssert( zMatAdj + zMatDet (matrix), result );
}

void assert_mat_det_scale(void)
{
  /* powers of two, by which matrices are scaled without round-off errors */
  const double scale[] = { 0x1p-10, 0x1p-7, 0x1p-3, 1.0, 0x1p10 };
  const int size = 10, testnum = 100;
  const double tol = 1.0e-12; /* tolerance of the relative error */
  zMat m, ms;
  double det, dets, s;
  int i, k;
  bool result = true;

  m = zMatAllocSqr( size );
  ms = zMatAllocSqr( size );
  for( i=0; i<sizeof(scale)/sizeof(double); i++ ){
    /* det( s I ) = s^n, which is not regarded as zero even if it is tiny */
    zMatIdent( m );
    zMatMulDRC( m, scale[i] );
    s = pow( scale[i], size );
    if( fabs( ( det = zMatDet( m ) ) - s ) > tol * s ){
      eprintf( "det(%g I) = %g (expected %g)\n", scale[i], det, s );
      result = false;
    }
    /* det( s a ) = s^n det( a ) */
    for( k=0; k<testnum; k++ ){
      zMatRandUniform( m, -10, 10 );
      zMatMul( m, scale[i], ms );
      det = zMatDet( m );
      dets = zMatDet( ms );
      if( fabs( dets - s * det ) > tol * fabs( s * det ) ){
        eprintf( "det(%g a) = %g, %g^n det(a) = %g\n", scale[i], dets, scale[i], s * det );
        result = false;
      }
    }
  }
  /* a singular matrix with two identical rows, whose determinant is negligible compared with the
   * Hadamard's bound, namely, the product of the norms of the row vectors */
  zMatRandUniform( m, -10, 10 );
  zRawVecCopy( zMatRowBufNC(m,0), zMatRowBufNC(m,size-1), size );
  for( s=1, i=0; i<size; i++ ) s *= zRawVecNorm( zMatRowBufNC(m,i), size );
  if( fabs( det = zMatDet( m ) ) > tol * s ){
    eprintf( "det of a singular matrix = %g (Hadamard's bound = %g)\n", det, s );
    result = false;
  }
  zMatFreeAtOnce( 2, m, ms );
  zAssert( zMatDet (scaled and singular matrices), result );
}

int main(void)
{
  zRandInit();
  assert_mat_mat_sweep_out();
  assert_mat_inv(); /* zMatInv() is used in assert_mat_det_adj() and assert_le() */
  assert_mat_det_adj(); /* zMatAdj() and zMatDet() are used in assert_le() */
  assert_mat_det_scale();
  assert_le();
  assert_mat_vec_sweep_out(); /* zLESolveGauss() is used to make the reference */
  assert_le_mat_capacity();
  assert_le_gauss_seidel();
  assert_le_gauss_seidel_scale();
  assert_le_gen();
  assert_le_gen_rowlong();
  assert_le_gen_aux();
  assert_le_mpnull();
  assert_le_ref_norm_min();
  assert_tridiagonal_equation();
  assert_lyapnov_equation();
  assert_lyapnov_equation_capacity();
  assert_mat_mul_inv();
  assert_mat_mul_inv_capacity();
  assert_mat_mul_inv_mat_vec();
  assert_mat_mpinv();
  assert_mat_mpinv_rand();
  assert_mpnull();
  return EXIT_SUCCESS;
}
