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
  assert_qp_lemke();
  return 0;
}
