#include <zm/zm_eig.h>
#include <zm/zm_rand.h>

void check(zMat m, zVec s, zMat r, int c)
{
  zVec v, e;
  register int i;

  v = zVecAlloc( zVecSizeNC(s) );
  e = zVecAlloc( zVecSizeNC(s) );
  for( i=0; i<zVecSizeNC(s); i++ ){
    zMatGetCol( r, i, v );
    zMulMatVec( m, v, e );
    zVecCatDRC( e, -zVecElem(s,i), v );
    printf( "eig#%d: eig.val.=%g,  err = %g\n", i, zVecElem(s,i), zVecNorm(e) );
  }
  eprintf( "computation time(clk) = %d\n", c );
  zVecFree( v );
  zVecFree( e );
}

#define N 100

int main(int argc, char *argv[])
{
  int n, i, j;
  zMat m, r;
  zVec s;
  clock_t c1, c2;

  n = argc > 1 ? atoi( argv[1] ) : N;
  m = zMatAllocSqr( n );
  r = zMatAllocSqr( n );
  s = zVecAlloc( n );
  zRandInitMT( NULL );
  for( i=0; i<n; i++ )
    for( j=i; j<n; j++ ){
      zMatSetElem( m, i, j, zRandMTF(NULL,-10,10) );
      zMatSetElem( m, j, i, zMatElem(m,i,j) );
    }

  c1 = clock();
  zEigSymBisec( m, s, r );
  c2 = clock();
  check( m, s, r, c2-c1 );

  c1 = clock();
  zEigSymJacobi( m, s, r );
  c2 = clock();
  check( m, s, r, c2-c1 );

  zMatFree( m );
  zMatFree( r );
  zVecFree( s );
  return 0;
}
