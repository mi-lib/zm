#include <zm/zm_eig.h>

void check(zMat m, zVec v, double s)
{
  zVec e;

  eprintf( "%.16g\n", s );
  e = zVecAlloc( zVecSizeNC(v) );
  zMulMatVec( m, v, e );
  zVecCatNCDRC( e, -s, v );
  eprintf( " %g ... %s.\n", zVecAbsMax(e,NULL), zVecIsTiny(e) ? "OK" : "maybe a bug" );
}

#define N 100

int main(int argc, char *argv[])
{
  zMat l, m, r;
  zVec v, e;
  double s;

  zRandInit();
  l = zMatAllocSqr( N );
  m = zMatAllocSqr( N );
  r = zMatAllocSqr( N );
  v = zVecAlloc( N );
  e = zVecAlloc( N );
  zMatRandUniform( l, -10, 10 );
  zMulMatTMat( l, l, m );

  zEigSymJacobi( m, e, r );
  eprintf( "maximal eigenvalue = %.16g\n", zVecMax(e,NULL) );
  eprintf( "minimal eigenvalue = %.16g\n", zVecMin(e,NULL) );

  eprintf( ">> maximal eigenvalue by Power method\n" );
  s = zEigPower( m, v, 0 );
  check( m, v, s );

  eprintf( ">> minimal eigenvalue by Power method\n" );
  s = zEigPowerInv( m, v, 0 );
  check( m, v, s );

  zMatFree( l );
  zMatFree( m );
  zMatFree( r );
  zVecFree( v );
  zVecFree( e );
  return 0;
}
