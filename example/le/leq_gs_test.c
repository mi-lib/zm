#include <zm/zm_le.h>

#define N 10

int main(void)
{
  zMat a;
  zVec b, x;
  zIndex index;
  register int i;

  zRandInit();
  a = zMatAllocSqr( N );
  b = zVecAlloc( N );
  x = zVecAlloc( N );
  index = zIndexCreate( N );

  zMatRand( a, -10, 10 );
  zVecRand( b, -10, 10 );
  for( i=0; i<N; i++ ){
    zPivoting( a, index, i, i );
    zMatSetElem( a, zIndexElem(index,i), i, zRandF(2*N,3*N) );
  }

  zMatImg( a );

  printf( "Gauss's elimination method.\n" );
  zLESolveGauss( a, b, x );
  zVecWrite( x );

  printf( "Gauss-Seidel's method.\n" );
  zVecClear( x );
  zLESolveGS( a, b, x );
  zVecWrite( x );
  printf( "(confirmation) A x\n" );
  zVecWrite( b );
  zMulMatVec( a, x, b );
  zVecWrite( b );

  zMatFree( a );
  zVecFree( b );
  zVecFree( x );
  zIndexFree( index );
  return 0;
}
