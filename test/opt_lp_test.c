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

bool test_lp_megiddo_dyer(int num_test, int num_constraint)
{
  zLPMegiddoDyer lp;
  double x1, x2, cost;
  double x1_ans, x2_ans, cost_ans;
  int i, j;
  double theta;
  bool result1, result2;

  for( i=0; i<num_test; i++ ){
    zLPMegiddoDyerInit( &lp );
    zLPMegiddoDyerSetCostCoefficient( &lp, zRandF(-5,5), zRandF(-5,5) );
    for( j=0; j<num_constraint; j++ ){
      theta = zRandF( -zPI, zPI );
      zLPMegiddoDyerAddConstraint( &lp, cos(theta), sin(theta), 5 + 5*cos(theta) );
    }
    result1 = zLPMegiddoDyerSolveDirect( &lp, &x1_ans, &x2_ans, &cost_ans );
    result2 = zLPMegiddoDyerSolve( &lp, &x1, &x2, &cost );
    zLPMegiddoDyerDestroy( &lp );
    if( !result1 && !result2 ) continue;
    if( !zEqual( x1, x1_ans, zTOL ) || !zEqual( x2, x2_ans, zTOL ) ){
      eprintf( "solution         : (x1*,x2*)=(%g,%g) (cost = %g)\n", x1, x2, cost );
      eprintf( "solution (direct): (x1*,x2*)=(%g,%g) (cost = %g)\n", x1_ans, x2_ans, cost_ans );
      return false;
    }
  }
  return true;
}

void assert_lp_megiddo_dyer(void)
{
  bool result;

  result = test_lp_megiddo_dyer( 100, 20 );
  zAssert( zLPMegiddoDyer (regular cases), result );
  result = test_lp_megiddo_dyer( 10, 5 );
  zAssert( zLPMegiddoDyer (likely unsolvable cases), result );
}

int main(void)
{
  zRandInit();
  assert_lp_simplex();
  assert_lp_megiddo_dyer();
  return 0;
}
