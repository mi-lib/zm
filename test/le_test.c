#include <zm/zm.h>

/* linear equation solver */

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
      eprintf( "maximum error = %g\n", zMatElemAbsMax( m3, NULL, NULL ) );
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

void assert_le(void)
{
  zMat a;
  zVec b, x, ans;
  zMat l, u; /* for LU decomposition */
  zMat adj;  /* Cramel's method */
  zIndex index;
  const double tol = ( 1.0e-9 );
  const int size = 10;
  const int n = 1000;
  int i;
  int count_cramel, count_gauss, count_lu, count_ri;

  a = zMatAllocSqr( size );
  adj = zMatAllocSqr( size );
  l = zMatAllocSqr( size );
  u = zMatAllocSqr( size );
  b = zVecAlloc( size );
  x = zVecAlloc( size );
  ans = zVecAlloc( size );
  index = zIndexAlloc( size );
  if( !a || !adj || !l || !u || !b || !x || !ans || !index ) goto TERMINATE;
  count_cramel = count_gauss = count_lu = count_ri = 0;
  for( i=0; i<n; i++ ){
    /* generate a problem */
    zMatRandUniform( a, -10, 10 );
    zVecRandUniform( ans,-10, 10 );
    zMulMatVec( a, ans, b );
    /* Cramel's method (adjoint matrix / determinant) */
    zMatAdj( a, adj );
    zMatDivDRC( adj, zMatDet(a) );
    zMulMatVec( adj, b, x );
    if( zVecIsTol( zVecSubDRC( x, ans ), tol*10 /* relaxed tolerance */ ) ) count_cramel++;
    else{
      eprintf( "error (Cramel's method) " );
      zVecFPrint( stderr, x );
    }
    /* Gauss's elimination method */
    zLESolveGauss( a, b, x );
    if( zVecIsTol( zVecSubDRC( x, ans ), tol ) ) count_gauss++;
    else{
      eprintf( "error (Gauss's elimination method) " );
      zVecFPrint( stderr, x );
    }
    /* LU decomposition method */
    zMatDecompLU( a, l, u, index );
    zLESolveLU( l, u, b, x, index );
    if( zVecIsTol( zVecSubDRC( x, ans ), tol ) ) count_lu++;
    else{
      eprintf( "error (LU decomposition method) " );
      zVecFPrint( stderr, x );
    }
    /* LU decomposition + residual iteration method */
    zLESolveLURI( a, b, x );
    if( zVecIsTol( zVecSubDRC( x, ans ), tol ) ) count_ri++;
    else{
      eprintf( "error (LU decomposition + residual iteration method) " );
      zVecFPrint( stderr, x );
    }
  }
 TERMINATE:
  zMatFreeAtOnce( 4, a, l, u, adj );
  zVecFreeAtOnce( 3, b, x, ans );
  zIndexFree( index );

  eprintf( "success rate (Cramel) %d/%d ", count_cramel, n ); zAssert( zMatAdj + zMatDet, count_cramel == n );
  eprintf( "success rate (Gauss)  %d/%d ", count_gauss, n );  zAssert( zLESolveGauss,     count_gauss == n );
  eprintf( "success rate (LU)     %d/%d ", count_lu, n );     zAssert( zLESolveLU,        count_lu == n );
  eprintf( "success rate (LU+RI)  %d/%d ", count_ri, n );     zAssert( zLESolveLURI,      count_ri == n );
}

void assert_le_gauss_seidel(void)
{
  zMat a, a_org;
  zVec b, x, ans;
  const double tol = ( 1.0e-10 );
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
  zLESolveMPLQ( a, b, w, w2, x );
  c2 = clock();
  zMulMatVec( a, x, _b );
  delta_clock[0] = c2 - c1;
  error[0] = zVecElemAbsMax( zVecSubDRC( _b, b ), NULL );
  eprintf( "MP(LQ):%g %ld ", error[0], delta_clock[0] );

  c1 = clock();
  zLESolveMPLU( a, b, w, w2, x );
  c2 = clock();
  zMulMatVec( a, x, _b );
  delta_clock[1] = c2 - c1;
  error[1] = zVecElemAbsMax( zVecSubDRC( _b, b ), NULL );
  eprintf( "MP(LU):%g %ld ", error[1], delta_clock[1] );

  c1 = clock();
  zLESolveMPSVD( a, b, x );
  c2 = clock();
  zMulMatVec( a, x, _b );
  delta_clock[2] = c2 - c1;
  error[2] = zVecElemAbsMax( zVecSubDRC( _b, b ), NULL );
  eprintf( "MP(SVD):%g %ld ", error[2], delta_clock[2] );

  c1 = clock();
  zLESolveSR( a, b, w, w2, x );
  c2 = clock();
  zMulMatVec( a, x, _b );
  delta_clock[3] = c2 - c1;
  error[3] = zVecElemAbsMax( zVecSubDRC( _b, b ), NULL );
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
  const double tol_general[] = { 1.0e-10, 1.0e-9, 1.0e-11, 1.0e-3 };
  const double tol_illcond[] = { 1.0e-10, 1.0e-9, 1.0e-2,  1.0e-2 };
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
  zAssert( zLESolveMPLQ, result_general[0] );
  zAssert( zLESolveMPLU, result_general[1] );
  zAssert( zLESolveMPSVD, result_general[2] );
  zAssert( zLESolveSR, result_general[3] );
  zAssert( zLESolveMPLQ (ill-conditioned cases), result_illcond[0] );
  zAssert( zLESolveMPLU (ill-conditioned cases), result_illcond[1] );
  zAssert( zLESolveMPSVD (ill-conditioned cases), result_illcond[2] );
  zAssert( zLESolveSR (ill-conditioned cases), result_illcond[3] );
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

/* special linear equation solver */

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
      eprintf( "case #%d: maximum error = %g\n", k, zMatElemAbsMax( tmp1, NULL, NULL ) );
      result = false;
    }
  }
  zMatFreeAtOnce( 5, a, b, x, tmp1, tmp2 );
  zAssert( zLELyapnovSolve, result );
}

/* matrix inversion */

void assert_mat_inv(void)
{
  zMat m1, m2, m;
  const int size = 10, testnum = 100;
  const double tol = 1.0e-10;
  bool result = true;
  int k;

  m1 = zMatAllocSqr( size );
  m2 = zMatAllocSqr( size );
  m  = zMatAllocSqr( size );
  for( k=0; k<testnum; k++ ){
    zMatRandUniform( m1, -10, 10 );
    zMatInv( m1, m2 );
    zMulMatMat( m1, m2, m );
    if( !zMatIsIdent( m, tol ) ){
      zMatTouchup( m, tol );
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
  const int rowsize = 10;
  const int colsize = 8;
  const int testnum = 100;
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
      eprintf( "case #%d: maximum error = %g\n", k, zMatElemAbsMax( m2, NULL, NULL ) );
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
      eprintf( "case #%d: maximum error = %g\n", k, zMatElemAbsMax( minv1, NULL, NULL ) );
      result = false;
    }
  }
  zMatFreeAtOnce( 3, m, minv1, minv2 );
  zAssert( zMatAdj + zMatDet (matrix), result );
}

int main(void)
{
  zRandInit();
  assert_mat_mat_sweep_out();
  assert_mat_vec_sweep_out();
  assert_le();
  assert_le_gauss_seidel();
  assert_le_gen();
  assert_le_gen_aux();
  assert_tridiagonal_equation();
  assert_lyapnov_equation();
  assert_mat_inv();
  assert_mat_mul_inv();
  assert_mat_mpinv();
  assert_mat_mpinv_rand();
  assert_mpnull();
  assert_mat_det_adj();
  return EXIT_SUCCESS;
}
