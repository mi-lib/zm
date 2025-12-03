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

bool check_matrix_composition(zMat mat, zMat left_mat, zMat right_mat)
{
  zMat mat_check;
  bool result;
  const double tol = 1.0e-8;

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

bool check_matrix_orthogonality(zMat q, int rank)
{
  zMat qqt;
  bool result;
  const double tol = 1.0e-9;

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

  mat = zMatAlloc( rowsize, colsize );
  l = zMatAlloc( rowsize, rowsize );
  u = zMatAlloc( rowsize, colsize );
  index = zIndexAlloc( zMatColSizeNC(l) );
  for( i=0; i<n; i++ ){
    generate_matrix_composition( mat, rank );
    zMatDecompLU( mat, l, u, index );
    if( check_matrix_composition( mat, l, u ) ) count_success++;
  }
  eprintf( "success rate (%d x %d) rank=%d, %d/%d ", rowsize, colsize, rank, count_success, n );
  zMatFreeAtOnce( 3, mat, l, u );
  zIndexFree( index );
  return count_success == n;
}

bool assert_mat_decomp_lq_one(int rowsize, int colsize, int rank, int n)
{
  zMat mat, l, q;
  int i, size, rank_result, count_success = 0;

  mat = zMatAlloc( rowsize, colsize );
  size = zMatMinSize( mat );
  l = zMatAlloc( rowsize, size );
  q = zMatAlloc( size, colsize );
  if( !mat || !l || !q ) goto TERMINATE;
  for( i=0; i<n; i++ ){
    zMatResetSize( l );
    zMatResetSize( q );
    generate_matrix_composition( mat, rank );
    rank_result = zMatDecompLQAndResize( mat, l, q );
    if( rank_result != rank ){
      eprintf( "assigned rank = %d / detected rank = %d, forcibly truncate.\n", rank, rank_result );
      zMatColResize( l, rank );
      zMatRowResize( q, rank );
    }
    if( check_matrix_composition( mat, l, q ) && check_matrix_orthogonality( q, rank ) ) count_success++;
  }
 TERMINATE:
  eprintf( "success rate (%d x %d) rank=%d, %d/%d ", zMatRowSize(mat), zMatColSize(mat), rank, count_success, n );
  zMatFreeAtOnce( 3, mat, l, q );
  return count_success == n;
}

bool assert_mat_decomp_lq_householder_one(int rowsize, int colsize, int rank, int n)
{
  zMat mat, l, q;
  int i, size, rank_result, count_success = 0;

  mat = zMatAlloc( rowsize, colsize );
  size = zMatMinSize( mat );
  l = zMatAlloc( zMatRowSize(mat), size );
  q = zMatAlloc( size, zMatColSize(mat) );
  if( !mat || !l || !q ) goto TERMINATE;
  for( i=0; i<n; i++ ){
    zMatResetSize( l );
    zMatResetSize( q );
    generate_matrix_composition( mat, rank );
    rank_result = zMatDecompLQ_Householder( mat, l, q );
    if( rank_result != rank ){
      eprintf( "assigned rank = %d / detected rank = %d, forcibly truncate.\n", rank, rank_result );
    }
    zMatColResize( l, rank );
    zMatRowResize( q, rank );
    if( check_matrix_composition( mat, l, q ) && check_matrix_orthogonality( q, rank ) ) count_success++;
  }
 TERMINATE:
  zMatFreeAtOnce( 3, mat, l, q );
  eprintf( "number of success (%d x %d) rank=%d) %d/%d ", zMatRowSize(mat), zMatColSize(mat), rank, count_success, n );
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
  zAssert( zMatDecompLU (5x8), assert_mat_decomp_lu( 5, 8, rank, n ) );
  zAssert( zMatDecompLU (8x5), assert_mat_decomp_lu( 8, 5, rank, n ) );
  zAssert( zMatDecompLU (8x8), assert_mat_decomp_lu( 8, 8, rank, n ) );
  zAssert( zMatDecompLQ (5x8), assert_mat_decomp_lq_one( 5, 8, rank, n ) );
  zAssert( zMatDecompLQ (8x5), assert_mat_decomp_lq_one( 8, 5, rank, n ) );
  zAssert( zMatDecompLQ (8x8), assert_mat_decomp_lq_one( 8, 8, rank, n ) );
  zAssert( zMatDecompLQ_Householder (5x8), assert_mat_decomp_lq_householder_one( 5, 8, 3, n ) );
  zAssert( zMatDecompLQ_Householder (8x5), assert_mat_decomp_lq_householder_one( 8, 5, 3, n ) );
  zAssert( zMatDecompLQ_Householder (8x8), assert_mat_decomp_lq_householder_one( 8, 8, 5, n ) );
  return 0;
}
