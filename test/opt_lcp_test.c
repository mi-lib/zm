#include <zm/zm_opt.h>

bool generate_lcp(int n, zMat *m, zVec *p, zVec *w_ans, zVec *z_ans)
{
  zMat r;
  zVec mz;
  int i;
  bool ret = false;

  *m = zMatAllocSqr( n );
  *p = zVecAlloc( n );
  *w_ans = zVecAlloc( n );
  *z_ans = zVecAlloc( n );
  mz = zVecAlloc( n );
  r = zMatAllocSqr( n );
  if( !*m || !*p || !*w_ans || !*z_ans || !mz || !r ) goto TERMINATE;
  for( i=0; i<n; i++ ){
    if( zRandI(0,1) == 0 ){
      zVecSetElemNC( *w_ans, i, 0 );
      zVecSetElemNC( *z_ans, i, zRandF( 0, 10 ) );
    } else{
      zVecSetElemNC( *w_ans, i, zRandF( 0, 10 ) );
      zVecSetElemNC( *z_ans, i, 0 );
    }
  }
  zMatRandUniform( r, -10, 10 );
  zMulMatTMat( r, r, *m );
  for( i=0; i<n; i++ )
    zMatElemNC(*m,i,i) += 1.0;
  zMatRandUniform( r, -1, 1 );
  zMatAddDRC( *m, r );
  zMatTDRC( r );
  zMatSubDRC( *m, r );

  zMulMatVec( *m, *z_ans, mz );
  zVecSub( *w_ans, mz, *p );

  ret = true;
 TERMINATE:
  zVecFree( mz );
  if( !ret ){
    zMatFree( *m );
    zVecFreeAtOnce( 3, *p, *w_ans, *z_ans );
  }
  return ret;
}

void check_answer(zVec w, zVec z, zVec w_ans, zVec z_ans, double tol, bool *result)
{
  if( !zVecEqual( w, w_ans, tol ) || !zVecEqual( z, z_ans, tol ) ){
    *result = false;
    eprintf( ">> failure case detected <<\n" );
    eprintf( "w:  " ); zVecFPrint( stderr, w );
    eprintf( "w*: " ); zVecFPrint( stderr, w_ans );
    eprintf( "z:  " ); zVecFPrint( stderr, z );
    eprintf( "z*: " ); zVecFPrint( stderr, z_ans );
  }
}

void assert_lcp_lemke(void)
{
  zMat m;
  zVec p, w, z, w_ans, z_ans;
  const int n = 100;
  const double tol = 1.0e-8;
  int num_trial = 10;
  bool result = true;

  for( ; num_trial>0; num_trial-- ){
    if( !generate_lcp( n, &m, &p, &w_ans, &z_ans ) ) return;
    w = zVecAlloc( n );
    z = zVecAlloc( n );
    zLCPSolveLemke( m, p, w, z );
    check_answer( w, z, w_ans, z_ans, tol, &result );
    zMatFree( m );
    zVecFreeAtOnce( 5, p, w, z, w_ans, z_ans );
  }
  zAssert( zLCPSolveLemke (regular case), result );
}

int main(void)
{
  zRandInit();
  assert_lcp_lemke();
  return 0;
}
