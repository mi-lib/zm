#include <zm/zm_eig.h>
#include <zm/zm_rand.h>

#define N 10
int main(int argc, char *argv[])
{
  int n;
  int i;
  zComplex *z;
  zMat ma;
  zCMat cma;
  zCVec *ve, vt;

  zRandInit();
  n = argc > 1 ? atoi( argv[1] ) : N;
  ma = zMatAllocSqr( n );
  zMatRandUniform( ma, -10, 10 );
  z = zAlloc( zComplex, n );
  ve = zAlloc( zCVec, n );
  for( i=0; i<n; i++ )
    ve[i] = zCVecAlloc( n );

  zEigSystem( ma, z, ve, 0 );

  /* ensurance */
  cma = zCMatAlloc( n, n );
  zMat2CMat( ma, cma );
  vt = zCVecAlloc( n );
  for( i=0; i<n; i++ ){
    printf( "eig.#%d val= ", i );
    zComplexWrite( &z[i] );
    printf( "\n" );
    printf( "eig.#%d vec.\n", i );
    zCVecWrite( ve[i] );
    zCMulMatVec( cma, ve[i], vt );
    zComplexRev( &z[i], &z[i] );
    zCVecCatDRC( vt, &z[i], ve[i] );
    printf( " (err) = %g\n", zCVecNorm(vt) );
  }

  zMatFree( ma );
  zCMatFree( cma );
  zCVecFree( vt );
  for( i=0; i<n; i++ )
    zCVecFree( ve[i] );
  zFree( ve );
  zFree( z );
  return 0;
}
