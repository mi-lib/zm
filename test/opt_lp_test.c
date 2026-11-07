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
  const double tol = 1.0e-10;

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
    if( !zEqual( x1, x1_ans, tol ) || !zEqual( x2, x2_ans, tol ) ){
      eprintf( "solution         : (x1*,x2*)=(%g,%g) (cost = %g)\n", x1, x2, cost );
      eprintf( "solution (direct): (x1*,x2*)=(%g,%g) (cost = %g)\n", x1_ans, x2_ans, cost_ans );
      return false;
    }
  }
  return true;
}

/* compare the solution of Megiddo-Dyer algorithm with that of simplex method, where the problem
 * with unbounded variables is converted to the standard form for the latter. */
bool test_lp_megiddo_dyer_simplex(int num_test, int num_constraint)
{
  zLPMegiddoDyer lp;
  zMat a, as;
  zVec b, c, x, cs, xs;
  double x1, x2, cost, x1_ans = 0, x2_ans = 0, cost_ans;
  int i, j;
  double theta;
  bool result1, result2, ret = false;
  const double tol = 1.0e-8; /* relative tolerance */

  a = zMatAlloc( num_constraint, 2 );
  b = zVecAlloc( num_constraint );
  c = zVecAlloc( 2 );
  x = zVecAlloc( 2 );
  if( !a || !b || !c || !x ) goto TERMINATE;
  for( i=0; i<num_test; i++ ){
    zLPMegiddoDyerInit( &lp );
    zVecSetElemNC( c, 0, zRandF(-5,5) );
    zVecSetElemNC( c, 1, zRandF(-5,5) );
    zLPMegiddoDyerSetCostCoefficient( &lp, zVecElemNC(c,0), zVecElemNC(c,1) );
    for( j=0; j<num_constraint; j++ ){
      theta = zRandF( -zPI, zPI );
      zMatSetElemNC( a, j, 0, cos(theta) );
      zMatSetElemNC( a, j, 1, sin(theta) );
      zVecSetElemNC( b, j, 5 + 5*cos(theta) );
      zLPMegiddoDyerAddConstraint( &lp, zMatElemNC(a,j,0), zMatElemNC(a,j,1), zVecElemNC(b,j) );
    }
    result1 = zLPMegiddoDyerSolve( &lp, &x1, &x2, &cost );
    zLPMegiddoDyerDestroy( &lp );
    if( !zLPUnb2Std( a, c, x, &as, &cs, &xs ) ) goto TERMINATE;
    if( ( result2 = zLPSolveSimplex( as, b, cs, xs, &cost_ans ) ) ){
      x1_ans = zVecElemNC(xs,0) - zVecElemNC(xs,2);
      x2_ans = zVecElemNC(xs,1) - zVecElemNC(xs,3);
    }
    zMatFree( as );
    zVecFreeAtOnce( 2, cs, xs );
    if( !result1 && !result2 ) continue;
    if( result1 != result2 ||
        sqrt( zSqr(x1-x1_ans) + zSqr(x2-x2_ans) ) > tol * ( 1 + sqrt( zSqr(x1) + zSqr(x2) ) ) ){
      eprintf( "solution          : (x1*,x2*)=(%g,%g) (cost = %g)\n", x1, x2, cost );
      eprintf( "solution (simplex): (x1*,x2*)=(%g,%g) (cost = %g)\n", x1_ans, x2_ans, cost_ans );
      goto TERMINATE;
    }
  }
  ret = true;
 TERMINATE:
  zMatFree( a );
  zVecFreeAtOnce( 3, b, c, x );
  return ret;
}

void assert_lp_megiddo_dyer(void)
{
  bool result;

  result = test_lp_megiddo_dyer( 100, 20 );
  zAssert( zLPMegiddoDyer (regular cases), result );
  result = test_lp_megiddo_dyer( 10, 5 );
  zAssert( zLPMegiddoDyer (likely unsolvable cases), result );
  result = test_lp_megiddo_dyer_simplex( 100, 20 );
  zAssert( zLPMegiddoDyer (comparison with simplex method), result );
  result = test_lp_megiddo_dyer_simplex( 10, 5 );
  zAssert( zLPMegiddoDyer (comparison with simplex method in likely unsolvable cases), result );
}

int main(void)
{
  zRandInit();
  assert_lp_simplex();
  assert_lp_megiddo_dyer();
  return 0;
}
