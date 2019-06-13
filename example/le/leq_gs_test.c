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

  zMatRandUniform( a, -10, 10 );
  zVecRandUniform( b, -10, 10 );
  for( i=0; i<N; i++ ){
    zPivoting( a, index, i, i );
    zMatSetElem( a, zIndexElem(index,i), i, zRandF(2*N,3*N) );
  }

  zMatImg( a );

  printf( "Gauss's elimination method.\n" );
  zLESolveGauss( a, b, x );
  zVecPrint( x );

  printf( "Gauss-Seidel's method.\n" );
  zVecClear( x );
  zLESolveGS( a, b, x );
  zVecPrint( x );
  printf( "(confirmation) A x\n" );
  zVecPrint( b );
  zMulMatVec( a, x, b );
  zVecPrint( b );

  zMatFree( a );
  zVecFree( b );
  zVecFree( x );
  zIndexFree( index );
  return 0;
}
