#include <zm/zm_opt.h>

bool generate_qp(int n, int m, double xmin, double xmax, zMat *q, zVec *c, zMat *a, zVec *b, zVec *ans)
{
  zMat r;
  zVec l, d;
  bool ret = false;
  int i;

  *q = zMatAllocSqr( n );
  *c = zVecAlloc( n );
  *a = zMatAlloc( m, n );
  *b = zVecAlloc( m );
  *ans = zVecAlloc( n );
  r = zMatAllocSqr( n );
  l = zVecAlloc( m );
  d = zVecAlloc( n );
  if( !*q || !*c || !*a || !*b || !*ans || !r || !l || !d ) goto TERMINATE;

  zVecRandUniform( *ans, xmin, xmax ); /* optimum answer */
  for( i=0; i<zVecSizeNC(l); i++ ) /* lambda vector */
    zVecSetElemNC( l, i, zRandI(0,1) == 0 ? zRandF( 1.0, 5.0 ) : 0 );

  zMatRandUniform( r, -10, 10 ); /* Q matrix */
  zMulMatTMat( r, r, *q );
  for( i=0; i<n; i++ )
    zMatElemNC(*q,i,i) += 1.0; /* to guarantee positive-definiteness of the matrix. */
  zMatRandUniform( *a, -10, 10 ); /* A matrix */
  zMulMatTVec( *a, l, *c ); /* c vector */
  zMulMatVec( *q, *ans, d );
  zVecSubDRC( *c, d );
  zMulMatVec( *a, *ans, *b ); /* b vector */
  for( i=0; i<zVecSizeNC(l); i++ )
    if( zVecElemNC(l,i) == 0 ) zVecElemNC(*b,i) -= zRandF( 1.0, 10.0 );
  ret = true;

 TERMINATE:
  zMatFree( r );
  zVecFree( l );
  zVecFree( d );
  if( !ret ){
    zMatFree( *q );
    zVecFree( *c );
    zMatFree( *a );
    zVecFree( *b );
    zVecFree( *ans );
  }
  return ret;
}

void check_answer(zVec x, zVec ans, double tol, bool *result)
{
  if( !zVecEqual( x, ans, tol ) ){
    *result = false;
    eprintf( ">> failure case detected <<\n" );
    zVecSubDRC( x, ans );
    zVecFPrint( stderr, x );
  }
}

void assert_qp_asm(void)
{
  zMat q, a;
  zVec c, b, ans, x;
  double cost;
  const int n = 10, m = 10;
  const double tol = 1.0e-8;
  int num_trial = 10;
  bool result = true;

  for( ; num_trial>0; num_trial-- ){
    if( !generate_qp( n, m, -10, 10, &q, &c, &a, &b, &ans ) ) return;
    x = zVecAlloc( n );
    zQPSolveASM( q, c, a, b, x, &cost );
    check_answer( x, ans, tol, &result );
    zMatFreeAtOnce( 2, q, a );
    zVecFreeAtOnce( 4, c, b, x, ans );
  }
  zAssert( zQPSolveASM (regular case), result );
}

/* a random matrix with orthonormal rows (the row size has to be not greater than the column size). */
bool mat_rand_orthonormal_row(zMat m)
{
  zMat r, l;
  zIndex index;
  bool ret;

  r = zMatAlloc( zMatRowSizeNC(m), zMatColSizeNC(m) );
  l = zMatAllocSqr( zMatRowSizeNC(m) );
  index = zIndexCreate( zMatRowSizeNC(m) );
  if( ( ret = r && l && index ) ){
    zMatRandUniform( r, -1, 1 );
    ret = zMatDecompLQ( r, l, m, index ) == zMatRowSizeNC(m);
  }
  zMatFreeAtOnce( 2, r, l );
  zIndexFree( index );
  return ret;
}

/* generate a problem with a known solution, where the condition numbers of q and of the active
 * constraints are kq and ka, respectively, nd constraints linearly dependent on the active ones
 * are also active, and the solution is scaled by s. */
bool generate_qp_illcond(int n, int m, int na, int nd, double kq, double ka, double s, zMat q, zVec c, zMat a, zVec b, zVec ans)
{
  zMat u, f, g;
  zVec l, d;
  int i, j;
  bool ret = false;

  u = zMatAllocSqr( n );
  f = zMatAllocSqr( na );
  g = zMatAlloc( na, n );
  l = zVecAlloc( m );
  d = zVecAlloc( n );
  if( !u || !f || !g || !l || !d ) goto TERMINATE;
  if( !mat_rand_orthonormal_row( u ) || !mat_rand_orthonormal_row( f ) || !mat_rand_orthonormal_row( g ) ) goto TERMINATE;
  /* q = u^T diag(1...1/kq) u */
  for( i=0; i<n; i++ )
    zVecSetElemNC( d, i, pow( kq, -(double)i/(n-1) ) );
  for( i=0; i<n; i++ )
    for( j=0; j<n; j++ )
      zMatSetElemNC( q, i, j, 0 );
  for( i=0; i<n; i++ )
    for( j=0; j<n; j++ )
      zRawVecCatDRC( zMatRowBufNC(q,j), zVecElemNC(d,i)*zMatElemNC(u,i,j), zMatRowBufNC(u,i), n );
  for( i=0; i<n; i++ ) /* to be exactly symmetric */
    for( j=0; j<i; j++ )
      zMatElemNC(q,i,j) = zMatElemNC(q,j,i) = 0.5 * ( zMatElemNC(q,i,j) + zMatElemNC(q,j,i) );
  /* active constraints f diag(1...1/ka) g, and constraints dependent on them */
  zMatRandUniform( a, -1, 1 );
  for( i=0; i<na; i++ )
    zRawVecMulDRC( zMatRowBufNC(g,i), pow( ka, -(double)i/(na-1) ), n );
  for( i=0; i<na; i++ ){
    zRawVecZero( zMatRowBufNC(a,i), n );
    for( j=0; j<na; j++ )
      zRawVecCatDRC( zMatRowBufNC(a,i), zMatElemNC(f,i,j), zMatRowBufNC(g,j), n );
  }
  for( i=na; i<na+nd; i++ ){
    zRawVecMul( zMatRowBufNC(a,zRandI(0,na-1)), zRandF(-1,1), zMatRowBufNC(a,i), n );
    zRawVecCatDRC( zMatRowBufNC(a,i), zRandF(-1,1), zMatRowBufNC(a,zRandI(0,na-1)), n );
  }
  /* the solution, Lagrange multipliers and the constant vectors */
  zVecRandUniform( ans, -10*s, 10*s );
  zMulMatVec( a, ans, b );
  for( i=0; i<m; i++ ){
    if( i < na + nd )
      zVecSetElemNC( l, i, zRandF( 1, 5 ) * s );
    else{
      zVecSetElemNC( l, i, 0 );
      zVecElemNC(b,i) -= zRandF( 1, 10 ) * s * zRawVecNorm( zMatRowBufNC(a,i), n );
    }
  }
  zMulMatTVec( a, l, c );
  zMulMatVec( q, ans, d );
  zVecSubDRC( c, d );
  for( i=m-1; i>0; i-- ){ /* shuffle the constraints */
    j = zRandI( 0, i );
    zMatSwapRow( a, i, j );
    zVecSwap( b, i, j );
  }
  ret = true;
 TERMINATE:
  zMatFreeAtOnce( 3, u, f, g );
  zVecFreeAtOnce( 2, l, d );
  return ret;
}

/* relative error of the solution. */
bool check_answer_rel(zVec x, zVec ans, double tol)
{
  zVecSubDRC( x, ans );
  if( zVecNorm( x ) <= tol * zVecNorm( ans ) ) return true;
  eprintf( ">> failure case detected (relative error = %g) <<\n", zVecNorm( x ) / zVecNorm( ans ) );
  return false;
}

void assert_qp_asm_illcond(void)
{
  zMat q, a;
  zVec c, b, ans, x;
  const int n = 30, m = 60, na = 15;
  const double tol_illcond = 1.0e-7, tol_degenerate = 1.0e-8;
  double s;
  int num_trial = 20;
  bool result_illcond = true, result_degenerate = true;

  q = zMatAllocSqr( n );
  a = zMatAlloc( m, n );
  c = zVecAlloc( n );
  b = zVecAlloc( m );
  ans = zVecAlloc( n );
  x = zVecAlloc( n );
  if( !q || !a || !c || !b || !ans || !x ) goto TERMINATE;
  for( ; num_trial>0; num_trial-- ){
    s = pow( 10, zRandF( -8, 8 ) );
    /* ill-conditioned quadratic term and active constraints */
    if( !generate_qp_illcond( n, m, na, 0, 1.0e6, 1.0e6, s, q, c, a, b, ans ) ) break;
    if( !zQPSolveASM( q, c, a, b, x, NULL ) || !check_answer_rel( x, ans, tol_illcond ) )
      result_illcond = false;
    /* active constraints linearly dependent on each other */
    if( !generate_qp_illcond( n, m, na, 5, 1.0e4, 1.0e4, s, q, c, a, b, ans ) ) break;
    if( !zQPSolveASM( q, c, a, b, x, NULL ) || !check_answer_rel( x, ans, tol_degenerate ) )
      result_degenerate = false;
  }
 TERMINATE:
  zMatFreeAtOnce( 2, q, a );
  zVecFreeAtOnce( 4, c, b, ans, x );
  zAssert( zQPSolveASM (ill-conditioned case), result_illcond );
  zAssert( zQPSolveASM (linearly dependent active constraints), result_degenerate );
}

void assert_qp_lemke(void)
{
  zMat q, a;
  zVec c, b, ans, x;
  double cost;
  const int n = 100, m = 50;
  const double tol = 1.0e-8;
  int num_trial = 10;
  bool result = true;

  for( ; num_trial>0; num_trial-- ){
    if( !generate_qp( n, m, 0, 10, &q, &c, &a, &b, &ans ) ) return;
    x = zVecAlloc( n );
    zQPSolveLemke( q, c, a, b, x, &cost );
    check_answer( x, ans, tol, &result );
    zMatFreeAtOnce( 2, q, a );
    zVecFreeAtOnce( 4, c, b, x, ans );
  }
  zAssert( zQPSolveLemke (regular case), result );
}

int main(void)
{
  zRandInit();
  assert_qp_asm();
  assert_qp_asm_illcond();
  assert_qp_lemke();
  return 0;
}
