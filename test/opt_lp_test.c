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

/* generate a problem with redundant equality constraints from a regular problem with m independent
 * equality constraints, where nr redundant constraints, which are linear combinations of the independent
 * ones, are inserted at random positions, and the constraints are scaled by s. */
bool generate_lp_redundant(int n, int m, int nr, double s, zVec *c, zMat *a, zVec *b, zVec *ans)
{
  zMat a0;
  zVec b0;
  zIndex index;
  int i, j, k, p, q;
  double cp, cq;

  if( !generate_lp( n, m, c, &a0, &b0, ans ) ) return false;
  *a = zMatAlloc( m+nr, n );
  *b = zVecAlloc( m+nr );
  index = zIndexCreate( m+nr );
  if( !*a || !*b || !index ){
    zMatFree( *a );
    zVecFreeAtOnce( 3, *b, *c, *ans );
    zMatFree( a0 );
    zVecFree( b0 );
    zIndexFree( index );
    return false;
  }
  zIndexShuffle( index, 0 );
  for( i=0; i<m; i++ )
    for( j=0; j<n; j++ )
      zMatSetElemNC( *a, zIndexElemNC(index,i), j, zMatElemNC(a0,i,j) );
  for( ; i<m+nr; i++ ){
    p = zRandI( 0, m-1 );
    q = zRandI( 0, m-1 );
    cp = zRandF( -2, 2 );
    cq = zRandI( 0, 1 ) ? zRandF( -2, 2 ) : 0;
    k = zIndexElemNC(index,i);
    for( j=0; j<n; j++ )
      zMatSetElemNC( *a, k, j, cp*zMatElemNC(a0,p,j) + cq*zMatElemNC(a0,q,j) );
  }
  zMatMulDRC( *a, s );
  zMulMatVec( *a, *ans, *b );
  zMatFree( a0 );
  zVecFree( b0 );
  zIndexFree( index );
  return true;
}

/* check if x satisfies a x = b and x >= 0. */
bool check_lp_feasible(zMat a, zVec b, zVec x, double tol)
{
  zVec r;
  bool ret;

  if( !( r = zVecAlloc( zVecSizeNC(b) ) ) ) return false;
  zMulMatVec( a, x, r );
  zVecSubDRC( r, b );
  ret = zVecNorm( r ) <= tol * ( 1 + zVecNorm( b ) ) && zVecElemMin( x, NULL ) >= -tol * ( 1 + zVecElemAbsMax( x, NULL ) );
  zVecFree( r );
  return ret;
}

void assert_lp_simplex_redundant(void)
{
  zMat a;
  zVec c, b, ans, x;
  double cost, s;
  const int n = 10, m = 5;
  const double tol = 1.0e-8;
  int num_trial = 100;
  bool result_solve = true, result_feasible = true;

  for( ; num_trial>0; num_trial-- ){
    s = pow( 10, zRandF( -6, 6 ) );
    if( !generate_lp_redundant( n, m, zRandI( 1, 3 ), s, &c, &a, &b, &ans ) ) return;
    if( !( x = zVecAlloc( n ) ) ) break;
    if( !zLPSolveSimplex( a, b, c, x, &cost ) ){
      eprintf( ">> failure case detected (scale = %g) <<\n", s );
      result_solve = false;
    } else
      check_answer( x, ans, tol, &result_solve );
    if( !zLPFeasibleBase( a, b, x ) || !check_lp_feasible( a, b, x, tol ) ){
      eprintf( ">> infeasible base detected (scale = %g) <<\n", s );
      result_feasible = false;
    }
    zMatFree( a );
    zVecFreeAtOnce( 4, c, b, x, ans );
  }
  zAssert( zLPSolveSimplex (redundant equality constraints), result_solve );
  zAssert( zLPFeasibleBase (redundant equality constraints), result_feasible );
}

/* generate an infeasible problem from a regular problem with m equality constraints, where a copy of
 * the first constraint with a different constant is appended. */
bool generate_lp_inconsistent(int n, int m, zVec *c, zMat *a, zVec *b)
{
  zMat a0;
  zVec b0, ans;
  int i, j;

  if( !generate_lp( n, m, c, &a0, &b0, &ans ) ) return false;
  *a = zMatAlloc( m+1, n );
  *b = zVecAlloc( m+1 );
  if( *a && *b ){
    for( i=0; i<m; i++ ){
      for( j=0; j<n; j++ )
        zMatSetElemNC( *a, i, j, zMatElemNC(a0,i,j) );
      zVecSetElemNC( *b, i, zVecElemNC(b0,i) );
    }
    for( j=0; j<n; j++ )
      zMatSetElemNC( *a, m, j, zMatElemNC(a0,0,j) );
    zVecSetElemNC( *b, m, zVecElemNC(b0,0) + 1 + fabs( zVecElemNC(b0,0) ) );
  }
  zMatFree( a0 );
  zVecFreeAtOnce( 2, b0, ans );
  if( !*a || !*b ){
    zMatFree( *a );
    zVecFreeAtOnce( 2, *b, *c );
    return false;
  }
  return true;
}

/* an estimate of the error of the answer of zLPSolvePDIP_PC() for a problem with the answer ans, which
 * is based on the termination condition that the norm of the residual of the constraints and the duality
 * gap are not greater than zTOL times |b| and |c||x| + |b||y|, respectively. Let B and N be the columns
 * of a corresponding to the basic and the non-basic variables, respectively. The error of the non-basic
 * variables dxN is bounded by the duality gap divided by the dual slack variables z_N = c_N - N^T y, where
 * y satisfies B^T y = c_B, and the error of the basic variables is bounded by B^+ ( rp - N dxN ). */
double lp_pdip_err_bound(zMat a, zVec b, zVec c, zVec ans)
{
  zMat bm, nm, bp, bpn;
  zVec y, cb, zn;
  int i, j, kb, kn, nb;
  double gap, dxn = 0, bound = HUGE_VAL;

  for( nb=0, j=0; j<zVecSizeNC(ans); j++ )
    if( zVecElemNC(ans,j) > 0 ) nb++;
  bm = zMatAlloc( zMatRowSizeNC(a), nb );
  nm = zMatAlloc( zMatRowSizeNC(a), zVecSizeNC(ans)-nb );
  bp = zMatAlloc( nb, zMatRowSizeNC(a) );
  bpn = zMatAlloc( nb, zVecSizeNC(ans)-nb );
  y = zVecAlloc( zMatRowSizeNC(a) );
  cb = zVecAlloc( nb );
  zn = zVecAlloc( zVecSizeNC(ans)-nb );
  if( !bm || !nm || !bp || !bpn || !y || !cb || !zn ) goto TERMINATE;
  for( kb=kn=0, j=0; j<zVecSizeNC(ans); j++ )
    if( zVecElemNC(ans,j) > 0 ){
      for( i=0; i<zMatRowSizeNC(a); i++ ) zMatSetElemNC( bm, i, kb, zMatElemNC(a,i,j) );
      zVecSetElemNC( cb, kb++, zVecElemNC(c,j) );
    } else{
      for( i=0; i<zMatRowSizeNC(a); i++ ) zMatSetElemNC( nm, i, kn, zMatElemNC(a,i,j) );
      zVecSetElemNC( zn, kn++, zVecElemNC(c,j) );
    }
  if( zMatMPInv( bm, bp ) < 0 ) goto TERMINATE;
  zMulMatTVec( bp, cb, y ); /* B^T y = c_B */
  for( j=0; j<zVecSizeNC(zn); j++ )
    for( i=0; i<zMatRowSizeNC(nm); i++ )
      zVecElemNC(zn,j) -= zMatElemNC(nm,i,j) * zVecElemNC(y,i);
  zMulMatMat( bp, nm, bpn );
  gap = zTOL * ( zVecNorm(c) * zVecNorm(ans) + zVecNorm(b) * zVecNorm(y) );
  for( j=0; j<zVecSizeNC(zn); j++ ) dxn += 1 / zSqr( zVecElemNC(zn,j) );
  dxn = gap * sqrt( dxn );
  bound = dxn * ( 1 + zMatNorm(bpn) ) + zMatNorm(bp) * zTOL * zVecNorm(b);
 TERMINATE:
  zMatFreeAtOnce( 4, bm, nm, bp, bpn );
  zVecFreeAtOnce( 3, y, cb, zn );
  return bound;
}

/* solve a problem by zLPSolvePDIP_PC() and check if the error of the answer is within the estimate
 * by lp_pdip_err_bound() with a safety factor. */
bool check_lp_pdip(zMat a, zVec b, zVec c, zVec ans, zVec x)
{
  double cost, err, bound;
  const double safety = 10;

  if( !zLPSolvePDIP_PC( a, b, c, x, &cost ) ){
    eprintf( ">> failure case detected (unsolved) <<\n" );
    return false;
  }
  if( ( err = zVecDist( x, ans ) ) > safety * ( bound = lp_pdip_err_bound( a, b, c, ans ) ) ){
    eprintf( ">> failure case detected (error %g > estimate %g) <<\n", err, bound );
    return false;
  }
  return true;
}

void assert_lp_pdip(void)
{
  zMat a;
  zVec c, b, ans, x;
  double cost;
  const int n = 10, m = 5, num_trial = 100;
  int i;
  bool result_regular = true, result_scale = true, result_redundant = true, result_unsolvable = true;

  for( i=0; i<num_trial; i++ ){
    if( !generate_lp( n, m, &c, &a, &b, &ans ) ) return;
    if( !( x = zVecAlloc( n ) ) ) break;
    if( !check_lp_pdip( a, b, c, ans, x ) ) result_regular = false;
    /* scaled cost and constraints */
    zVecMulDRC( c, pow( 10, zRandF( -6, 6 ) ) );
    zMatMulDRC( a, pow( 10, zRandF( -6, 6 ) ) );
    zVecMulDRC( ans, pow( 10, zRandF( -6, 6 ) ) );
    zMulMatVec( a, ans, b );
    if( !check_lp_pdip( a, b, c, ans, x ) ) result_scale = false;
    zMatFree( a );
    zVecFreeAtOnce( 4, c, b, x, ans );
  }
  zAssert( zLPSolvePDIP_PC (regular case), result_regular );
  zAssert( zLPSolvePDIP_PC (scaled case), result_scale );

  for( i=0; i<num_trial; i++ ){
    if( !generate_lp_redundant( n, m, zRandI( 1, 3 ), pow( 10, zRandF( -6, 6 ) ), &c, &a, &b, &ans ) ) return;
    if( !( x = zVecAlloc( n ) ) ) break;
    if( !check_lp_pdip( a, b, c, ans, x ) ) result_redundant = false;
    zMatFree( a );
    zVecFreeAtOnce( 4, c, b, x, ans );
  }
  zAssert( zLPSolvePDIP_PC (redundant equality constraints), result_redundant );

  /* inconsistent redundant constraints, where the last row is a copy of the first one */
  for( i=0; i<num_trial; i++ ){
    if( !generate_lp_inconsistent( n, m, &c, &a, &b ) ) return;
    if( !( x = zVecAlloc( n ) ) ) break;
    zEchoOff();
    if( zLPSolvePDIP_PC( a, b, c, x, &cost ) ){
      eprintf( ">> infeasible problem solved <<\n" );
      result_unsolvable = false;
    }
    zEchoOn();
    zMatFree( a );
    zVecFreeAtOnce( 3, c, b, x );
  }

  /* unbounded problem: min. -x2 subject to x1 - x2 = 1, x >= 0 */
  a = zMatCreateList( 1, 2, 1.0, -1.0 );
  b = zVecCreateList( 1, 1.0 );
  c = zVecCreateList( 2, 0.0, -1.0 );
  x = zVecAlloc( 2 );
  zEchoOff();
  if( zLPSolvePDIP_PC( a, b, c, x, &cost ) ) result_unsolvable = false;
  /* infeasible problem: min. x1 + x2 subject to x1 + x2 = -1, x >= 0 */
  zMatSetElemNC( a, 0, 1, 1.0 );
  zVecSetElemNC( b, 0, -1.0 );
  zVecSetElemNC( c, 0, 1.0 );
  zVecSetElemNC( c, 1, 1.0 );
  if( zLPSolvePDIP_PC( a, b, c, x, &cost ) ) result_unsolvable = false;
  zEchoOn();
  zMatFree( a );
  zVecFreeAtOnce( 3, b, c, x );
  zAssert( zLPSolvePDIP_PC (infeasible and unbounded cases), result_unsolvable );
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
  assert_lp_simplex_redundant();
  assert_lp_pdip();
  assert_lp_megiddo_dyer();
  return 0;
}
