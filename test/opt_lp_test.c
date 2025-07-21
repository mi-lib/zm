#include <zm/zm_opt.h>

bool generate_lp(int n, int m, zVec *c, zMat *a, zVec *b, zVec *ans)
{
  zVec ans_dual, l;
  zIndex index;
  bool ret = false;
  int i;

  if( m >= n ){
    eprintf( "cannot generate a solvable problem with the number of equalities %d > the number of variables %d\n", m, n );
    return false;
  }
  *c = zVecAlloc( n );
  *a = zMatAlloc( m, n );
  *b = zVecAlloc( m );
  *ans = zVecAlloc( n );
  ans_dual = zVecAlloc( n );
  l = zVecAlloc( m );
  index = zIndexCreate( n );
  if( !*c || !*a || !*b || !*ans || !ans_dual || !l || !index )
    goto TERMINATE;

  zIndexShuffle( index, 0 );
  for( i=0; i<m; i++ ){
    zVecSetElemNC( *ans, zIndexElemNC(index,i), zRandF( 0, 10 ) );
    zVecSetElemNC( ans_dual, zIndexElemNC(index,i), 0 );
  }
  for( ; i<n; i++ ){
    zVecSetElemNC( *ans, zIndexElemNC(index,i), 0 );
    zVecSetElemNC( ans_dual, zIndexElemNC(index,i), zRandF( 0, 10 ) );
  }
  zMatRandUniform( *a, -10, 10 ); /* NOTE: this does not strictly guarantee that a is row-fullrank. */
  zMulMatVec( *a, *ans, *b );
  zVecRandUniform( l, -5, 5 ); /* lambda vector */
  zMulMatTVec( *a, l, *c );
  zVecAddDRC( *c, ans_dual );
  ret = true;

 TERMINATE:
  zVecFreeAtOnce( 2, ans_dual, l );
  zIndexFree( index );
  if( !ret ){
    zVecFree( *c );
    zMatFree( *a );
    zVecFree( *b );
    zVecFree( *ans );
  }
  return ret;
}

int check_answer(zVec x, zVec ans, double tol, bool *result)
{
  if( !zVecEqual( x, ans, tol ) ){
    *result = false;
    eprintf( ">> failure case detected <<\n" );
    zVecFPrint( stderr, ans );
    zVecFPrint( stderr, x );
    return 0;
  }
  return 1;
}

void assert_lp_simplex(void)
{
  zMat a;
  zVec c, b, ans, x;
  double cost;
  const int n = 10, m = 6;
  const double tol = 1.0e-8;
  int num_trial = 100;
  bool result = true;

  for( ; num_trial>0; num_trial-- ){
    if( !generate_lp( n, m, &c, &a, &b, &ans ) ) return;
    x = zVecAlloc( n );
    zLPSolveSimplex( a, b, c, x, &cost );
    check_answer( x, ans, tol, &result );
    zMatFree( a );
    zVecFreeAtOnce( 4, c, b, x, ans );
  }
  zAssert( zLPSolveSimplex (regular case), result );
}

int main(void)
{
  zRandInit();
  assert_lp_simplex();
  return 0;
}
