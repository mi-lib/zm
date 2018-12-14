#include <zm/zm_eig.h>

void test(zMat m, zVec eig, zMat r, int n)
{
  register int i;
  zVec e, ev;

  e = zVecAlloc( _zMatRowSize(m) );
  ev= zVecAlloc( _zMatRowSize(r) );
  for( i=0; i<zVecSizeNC(eig); i++ ){
    zMatGetCol( r, i, ev );
    printf( "eig#%d:\n", i );
    printf( " eig-value = %g\n", zVecElem(eig,i) );
    printf( " eig-vec   = " ); zVecWrite(ev);
    zMulMatVec( m, ev, e );
    zVecCatDRC( e, -zVecElem(eig,i), ev );
    printf( " err = %g\n", zVecNorm(e) );
  }
  zVecFree( e );
  zVecFree( ev );
}

#define R 1.41421356237
#define N 4

int main(void)
{
  double array[] = {
    1, R, R, 2,
    R,-R,-1, R,
    R,-1, R, R,
    2, R, R,-3,
  };
  zMat m, r;
  zVec eig;

  m = zMatCloneArray( array, N, N );
  r = zMatAllocSqr( N );
  eig = zVecAlloc( N );

  printf( "original matrix\n" );
  zMatWrite( m );

  zEigSymJacobi( m, eig, r );
  test( m, eig, r, N );

  zMatFree( m );
  zMatFree( r );
  zVecFree( eig );
  return 0;
}
