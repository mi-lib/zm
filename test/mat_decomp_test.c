#include <zm/zm_le.h>

zMat generate_matrix_composition(zMat mat, int rank)
{
  zMat row_facto, col_facto;

  row_facto = zMatAlloc( zMatRowSize(mat), rank );
  col_facto = zMatAlloc( rank, zMatColSize(mat) );
  if( !row_facto || !col_facto ){
    ZALLOCERROR();
    mat = NULL;
    goto TERMINATE;
  }
  zMatRandUniform( row_facto, -10, 10 );
  zMatRandUniform( col_facto, -10, 10 );
  zMulMatMat( row_facto, col_facto, mat );
 TERMINATE:
  zMatFreeAtOnce( 2, row_facto, col_facto );
  return mat;
}

bool check_matrix_composition(zMat mat, zMat left_mat, zMat right_mat, double tol)
{
  zMat mat_check;
  bool result;

  if( !( mat_check = zMatAlloc( zMatRowSize(mat), zMatColSize(mat) ) ) ){
    ZALLOCERROR();
    return false;
  }
  zMulMatMat( left_mat, right_mat, mat_check );
  if( !( result = zMatEqual( mat, mat_check, tol ) ) ){
    eprintf( "matrix unrecovered.\n" );
    zMatSubDRC( mat_check, mat );
    zMatTouchup( mat_check, tol );
    eprintf( "LQ - A =" ); zMatFPrint( stderr, mat_check );
  }
  zMatFree( mat_check );
  return result;
}

bool check_matrix_orthogonality(zMat q, int rank, double tol)
{
  zMat qqt;
  bool result;

  qqt = zMatAllocSqr( rank );
  zMulMatMatT( q, q, qqt );
  if( !( result = zMatIsIdent( qqt, tol ) ) ){
    eprintf( "orthogonality unsatisfied.\n" );
    zMatTouchup( qqt, tol );
    eprintf( "Q Q^T =" ); zMatFPrint( stderr, qqt );
  }
  zMatFree( qqt );
  return result;
}

bool assert_mat_decomp_cholesky(int rowsize, int colsize, int n)
{
  int i;
  zMat m, mc, l, s;
  zIndex index;
  int count_success = 0;

  m = zMatAllocSqr( rowsize );
  mc = zMatAllocSqr( rowsize );
  l = zMatAllocSqr( rowsize );
  s = zMatAlloc( rowsize, colsize );
  index = zIndexAlloc( rowsize );
  for( i=0; i<n; i++ ){
    zMatResetSize( l );
    zMatRandUniform( s, -10, 10 );
    zMulMatMatT( s, s, m );
    zMatDecompCholesky( m, l, index );
    zMulMatMatT( l, l, mc );
    zMatSubDRC( m, mc );
    if( zMatIsTiny( m ) ) count_success++;
  }
  eprintf( "success rate (%d x %d).(%d x %d) %d/%d ", rowsize, colsize, colsize, rowsize, count_success, n );
  zMatFreeAtOnce( 4, m, mc, l, s );
  zIndexFree( index );
  return count_success == n;
}

bool assert_mat_decomp_lu(int rowsize, int colsize, int rank, int n)
{
  zMat mat, l, u;
  zIndex index;
  int i, count_success = 0;
  const double tol = 1.0e-11;

  mat = zMatAlloc( rowsize, colsize );
  l = zMatAlloc( rowsize, rowsize );
  u = zMatAlloc( rowsize, colsize );
  index = zIndexAlloc( zMatColSizeNC(l) );
  for( i=0; i<n; i++ ){
    zMatResetSize( l );
    zMatResetSize( u );
    generate_matrix_composition( mat, rank );
    zMatDecompLU( mat, l, u, index );
    if( check_matrix_composition( mat, l, u, tol ) ) count_success++;
  }
  eprintf( "success rate (%d x %d) rank=%d, %d/%d ", rowsize, colsize, rank, count_success, n );
  zMatFreeAtOnce( 3, mat, l, u );
  zIndexFree( index );
  return count_success == n;
}

bool assert_mat_decomp_lq_one(int rowsize, int colsize, int rank_max, int n)
{
  zMat mat, l, q;
  int i, size, rank, rank_result, count_success = 0;
  const double tol = 1.0e-9;

  mat = zMatAlloc( rowsize, colsize );
  size = zMatMinSize( mat );
  l = zMatAlloc( rowsize, size );
  q = zMatAlloc( size, colsize );
  if( !mat || !l || !q ) goto TERMINATE;
  for( i=0; i<n; i++ ){
    zMatResetSize( l );
    zMatResetSize( q );
    rank = zRandI( 1, rank_max );
    generate_matrix_composition( mat, rank );
    rank_result = zMatDecompLQ( mat, l, q );
    if( rank_result != rank )
      eprintf( "assigned rank = %d / detected rank = %d, forcibly truncate.\n", rank, rank_result );
    if( check_matrix_composition( mat, l, q, tol ) && check_matrix_orthogonality( q, rank, tol ) ) count_success++;
  }
 TERMINATE:
  eprintf( "success rate (%d x %d), %d/%d ", zMatRowSize(mat), zMatColSize(mat), count_success, n );
  zMatFreeAtOnce( 3, mat, l, q );
  return count_success == n;
}

/* zero matrix case. */
bool assert_mat_decomp_lq_zero(int rowsize, int colsize)
{
  zMat m, l, q;
  int size;
  bool result = false;

  size = zMin( rowsize, colsize );
  m = zMatAlloc( rowsize, colsize );
  l = zMatAlloc( rowsize, size );
  q = zMatAlloc( size, colsize );
  if( m && l && q ){
    zMatZero( m );
    if( zMatDecompLQ( m, l, q ) == 0 &&
        zMatRowSizeNC(l) == rowsize && zMatColSizeNC(l) == 0 &&
        zMatRowSizeNC(q) == 0 && zMatColSizeNC(q) == colsize ) result = true;
  }
  zMatFreeAtOnce( 3, m, l, q );
  return result;
}

/* skipped-row compaction in the Q-workspace. */
bool assert_mat_decomp_lq_dependent_rows(void)
{
  zMat m[2], l, q;
  int i, rank;
  bool result = false;
  const double tol = 1.0e-10;

  m[0] = zMatCreateList( 3, 3, /* square case to skip the second row. */
    1.0, 0.0, 0.0,
    2.0, 0.0, 0.0,
    0.0, 1.0, 0.0 );
  m[1] = zMatCreateList( 3, 4, /* column-long case to skip the second row. */
    1.0, 0.0, 0.0, 0.0,
    2.0, 0.0, 0.0, 0.0,
    0.0, 1.0, 0.0, 0.0 );
  l = zMatAllocSqr( 3 );
  q = zMatAlloc( 3, 4 );
  if( !m[0] || !m[1] || !l || !q ) goto TERMINATE;
  result = true;
  for( i=0; i<2; i++ ){
    zMatSetSize( l, zMatRowSizeNC( m[i] ), zMatRowSizeNC( m[i] ) );
    zMatSetSize( q, zMatRowSizeNC( m[i] ), zMatColSizeNC( m[i] ) );
    rank = zMatDecompLQ( m[i], l, q );
    if( rank != 2 ||
        !check_matrix_composition( m[i], l, q, tol ) ||
        !check_matrix_orthogonality( q, rank, tol ) ) result = false;
  }
 TERMINATE:
  zMatFreeAtOnce( 4, m[0], m[1], l, q );
  return result;
}

/* check if tiny errors are rejected in row/column compactions. */
bool assert_mat_decomp_lq_rank(void)
{
  zMat m[3], l, q;
  int i, size, rank;
  bool result = false;
  const double residual = 1.0e-12;

  m[0] = zMatCreateList( 3, 3, /* square case */
    1.0, 0.0, 0.0,
    0.0, 1.0, 0.0,
    0.0, 1.0, residual );
  m[1] = zMatCreateList( 3, 4, /* column-long case */
    1.0, 0.0, 0.0, 0.0,
    0.0, 1.0, 0.0, 0.0,
    0.0, 1.0, residual, 0.0 );
  m[2] = zMatCreateList( 4, 3, /* row-long case */
    1.0, 0.0, 0.0,
    0.0, 1.0, 0.0,
    0.0, 1.0, residual,
    0.0, 1.0, 0.0 );
  l = zMatAlloc( 4, 3 );
  q = zMatAlloc( 3, 4 );
  if( !m[0] || !m[1] || !m[2] || !l || !q ) goto TERMINATE;
  result = true;
  for( i=0; i<3; i++ ){
    size = zMatMinSize( m[i] );
    zMatSetSize( l, zMatRowSizeNC( m[i] ), size );
    zMatSetSize( q, size, zMatColSizeNC( m[i] ) );
    rank = zMatDecompLQ( m[i], l, q );
    if( rank != 2 || zMatColSizeNC(l) != rank || zMatRowSizeNC(q) != rank ) result = false;
  }
 TERMINATE:
  zMatFreeAtOnce( 5, m[0], m[1], m[2], l, q );
  return result;
}

/* check LQ decomposition of a matrix with the expected rank (negative value to skip rank check). */
bool check_mat_decomp_lq(zMat m, int rank_expected, double tol_comp, double tol_orth)
{
  zMat l, q;
  int size, rank;
  bool result = false;

  size = zMatMinSize( m );
  l = zMatAlloc( zMatRowSizeNC(m), size );
  q = zMatAlloc( size, zMatColSizeNC(m) );
  if( !l || !q ) goto TERMINATE;
  if( ( rank = zMatDecompLQ( m, l, q ) ) <= 0 ) goto TERMINATE;
  if( rank_expected >= 0 && rank != rank_expected ){
    eprintf( "assigned rank = %d / detected rank = %d\n", rank_expected, rank );
    goto TERMINATE;
  }
  result = check_matrix_composition( m, l, q, tol_comp ) && check_matrix_orthogonality( q, rank, tol_orth );
 TERMINATE:
  zMatFreeAtOnce( 2, l, q );
  return result;
}

/* check if rank detection is independent of the scale of the matrix. */
bool assert_mat_decomp_lq_scale(void)
{
  zMat m[2];
  const double scale[] = { 1.0e-12, 1.0, 1.0e8 };
  const int nscale = sizeof(scale) / sizeof(double);
  int i, k;
  bool result = false;

  for( k=0; k<nscale; k++ ){
    m[0] = zMatCreateList( 3, 4, /* column-long case: rank 2 (the second row is slightly independent) */
      1.0, 2.0, 3.0, 4.0,
      2.0, 4.0, 6.0, 8.000001,
      3.0, 6.0, 9.0, 12.0 );
    m[1] = zMatCreateList( 4, 3, /* row-long case: rank 2 (the second row is slightly independent) */
      1.0, 2.0, 3.0,
      2.0, 4.0, 6.000001,
      3.0, 6.0, 9.0,
      4.0, 8.0, 12.0 );
    if( !m[0] || !m[1] ) goto TERMINATE;
    for( i=0; i<2; i++ ){
      zMatMulDRC( m[i], scale[k] );
      if( !check_mat_decomp_lq( m[i], 2, 1.0e-12 * scale[k], 1.0e-7 ) ){
        eprintf( "failed for scale = %g (%d x %d)\n", scale[k], zMatRowSizeNC(m[i]), zMatColSizeNC(m[i]) );
        goto TERMINATE;
      }
    }
    zMatFreeAtOnce( 2, m[0], m[1] );
  }
  return true;
 TERMINATE:
  zMatFreeAtOnce( 2, m[0], m[1] );
  return result;
}

/* check accuracy of LQ decomposition of ill-conditioned (Hilbert) matrices. */
bool assert_mat_decomp_lq_illcond(void)
{
  zMat m, mt;
  int n, i, j;
  bool result = true;
  const double tol_comp = 1.0e-12, tol_orth = 1.0e-6;

  for( n=6; n<=8; n+=2 ){
    m = zMatAlloc( n+2, n );
    mt = zMatAlloc( n, n+2 );
    if( !m || !mt ){
      zMatFreeAtOnce( 2, m, mt );
      return false;
    }
    for( i=0; i<n+2; i++ )
      for( j=0; j<n; j++ ) zMatSetElemNC( m, i, j, 1.0/(i+j+1) );
    zMatTNC( m, mt );
    if( !check_mat_decomp_lq( m, n, tol_comp, tol_orth ) ){ /* row-long case */
      eprintf( "failed for %d x %d Hilbert matrix\n", n+2, n );
      result = false;
    }
    if( !check_mat_decomp_lq( mt, n, tol_comp, tol_orth ) ){ /* column-long case */
      eprintf( "failed for %d x %d Hilbert matrix\n", n, n+2 );
      result = false;
    }
    zMatFreeAtOnce( 2, m, mt );
  }
  return result;
}

bool assert_mat_decomp_lq_full_one(int rowsize, int colsize, int rank_max, int n)
{
  zMat mat, l, q;
  int i, rank, rank_result, count_success = 0;
  const double tol = 1.0e-10;

  mat = zMatAlloc( rowsize, colsize );
  l = zMatAlloc( zMatRowSize(mat), zMatColSize(mat) );
  q = zMatAllocSqr( zMatColSize(mat) );
  if( !mat || !l || !q ) goto TERMINATE;
  for( i=0; i<n; i++ ){
    zMatResetSize( l );
    zMatResetSize( q );
    rank = zRandI( 1, rank_max );
    generate_matrix_composition( mat, rank );
    rank_result = zMatDecompLQFull( mat, l, q );
    if( rank_result != rank ){
      eprintf( "assigned rank = %d / detected rank = %d, forcibly truncate.\n", rank, rank_result );
    }
    zMatColResize( l, rank );
    zMatRowResize( q, rank );
    if( check_matrix_composition( mat, l, q, tol ) && check_matrix_orthogonality( q, rank, tol ) ) count_success++;
  }
 TERMINATE:
  eprintf( "number of success (%d x %d), %d/%d ", zMatRowSize(mat), zMatColSize(mat), count_success, n );
  zMatFreeAtOnce( 3, mat, l, q );
  return count_success == n;
}

bool assert_mat_decomp_lq_nullspace(int rowsize, int colsize, int rank, int n)
{
  zMat mat, l, q, qnull, zero_check;
  int i, rank_result, count_success = 0;
  const double tol = 1.0e-10;

  mat = zMatAlloc( rowsize, colsize );
  l = zMatAlloc( zMatRowSize(mat), zMatColSize(mat) );
  q = zMatAllocSqr( zMatColSize(mat) );
  qnull = zMatAllocSqr( zMatColSize(mat) );
  if( !mat || !l || !q || !qnull ) goto TERMINATE;
  for( i=0; i<n; i++ ){
    zMatResetSize( l );
    zMatResetSize( q );
    generate_matrix_composition( mat, rank );
    rank_result = zMatDecompLQNull( mat, l, q, qnull );
    if( zMatColSizeNC(mat) == rank_result ){
      count_success++;
    } else{
      zero_check = zMatAlloc( rowsize, zMatColSizeNC(mat) - rank_result );
      zMulMatMat( mat, qnull, zero_check );
      if( zMatIsTol( zero_check, tol ) ){
        count_success++;
      } else{
        eprintf( "%g ", zMatElemAbsMax( zero_check, NULL, NULL ) );
      }
      zMatFree( zero_check );
    }
  }
 TERMINATE:
  eprintf( "number of success (%d x %d) rank=%d) %d/%d ", zMatRowSize(mat), zMatColSize(mat), rank, count_success, n );
  zMatFreeAtOnce( 4, mat, l, q, qnull );
  return count_success == n;
}

int main(void)
{
  const int size_large = 8;
  const int size_small = 5;
  const int rank = 4;
  const int n = 100;

  zRandInit();
  zAssert( zMatDecompCholesky (8x5), assert_mat_decomp_cholesky( size_large, size_small, n ) );
  zAssert( zMatDecompCholesky (8x8), assert_mat_decomp_cholesky( size_large, size_large, n ) );
  zAssert( zMatDecompLU (5x8), assert_mat_decomp_lu( size_small, size_large, rank, n ) );
  zAssert( zMatDecompLU (8x5), assert_mat_decomp_lu( size_large, size_small, rank, n ) );
  zAssert( zMatDecompLU (8x8), assert_mat_decomp_lu( size_large, size_large, rank, n ) );
  zAssert( zMatDecompLQ (5x8), assert_mat_decomp_lq_one( size_small, size_large, size_small, n ) );
  zAssert( zMatDecompLQ (8x5), assert_mat_decomp_lq_one( size_large, size_small, size_small, n ) );
  zAssert( zMatDecompLQ (8x8), assert_mat_decomp_lq_one( size_large, size_large, size_small, n ) );
  zAssert( zMatDecompLQ (row-long zero-matrix), assert_mat_decomp_lq_zero( size_large, size_small ) );
  zAssert( zMatDecompLQ (column-long zero-matrix), assert_mat_decomp_lq_zero( size_small, size_large ) );
  zAssert( zMatDecompLQ (row compactions), assert_mat_decomp_lq_dependent_rows() );
  zAssert( zMatDecompLQ (error rejection in compactions), assert_mat_decomp_lq_rank() );
  zAssert( zMatDecompLQ (scale-independent rank detection), assert_mat_decomp_lq_scale() );
  zAssert( zMatDecompLQ (ill-conditioned matrices), assert_mat_decomp_lq_illcond() );
  zAssert( zMatDecompLQFull (5x8), assert_mat_decomp_lq_full_one( size_small, size_large, size_small, n ) );
  zAssert( zMatDecompLQFull (8x5), assert_mat_decomp_lq_full_one( size_large, size_small, size_small, n ) );
  zAssert( zMatDecompLQFull (8x8), assert_mat_decomp_lq_full_one( size_large, size_large, size_small, n ) );
  zAssert( zMatDecompLQNull (5x8), assert_mat_decomp_lq_nullspace( size_small, size_large, rank, n ) );
  zAssert( zMatDecompLQNull (8x5), assert_mat_decomp_lq_nullspace( size_large, size_small, rank, n ) );
  return 0;
}
