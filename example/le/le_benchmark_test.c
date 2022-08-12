#include <zm/zm_le.h>

#define N  1000
#define SIZE 10
#define TOL (1.0e-10)

int main(void)
{
  zMat a, l, u, aa;
  zVec b, x, ans;
  zIndex index;
  int count_cramel, count_gauss, count_lu, count_ri, count_gs;
  register int i;

  zRandInit();
  a = zMatAllocSqr( SIZE );
  aa = zMatAllocSqr( SIZE );
  l = zMatAllocSqr( SIZE );
  u = zMatAllocSqr( SIZE );
  b = zVecAlloc( SIZE );
  x = zVecAlloc( SIZE );
  ans = zVecAlloc( SIZE );
  index = zIndexAlloc( SIZE );
  count_cramel = count_gauss = count_lu = count_ri = count_gs = 0;
  for( i=0; i<N; i++ ){
    /* generate a problem */
    zMatRandUniform( a,  -100, 100 );
    zVecRandUniform( ans,-10, 10 );
    zMulMatVec( a, ans, b );
    /* Cramel's method (adjoint matrix / determinant) */
    zMatAdj( a, aa );
    zMatDivDRC( aa, zMatDet(a) );
    zMulMatVec( aa, b, x );
    if( !zVecIsEqual( x, ans, TOL ) ){
      count_cramel++;
    }
    /* Gauss's elimination method */
    zLESolveGauss( a, b, x );
    if( !zVecIsEqual( x, ans, TOL ) ){
      count_gauss++;
    }
    /* LU decomposition method */
    zMatDecompLU( a, l, u, index );
    zLESolveLU( l, u, b, x, index );
    if( !zVecIsEqual( x, ans, TOL ) ){
      count_lu++;
    }
    /* Residual iteration method */
    zLESolveRI( a, b, x );
    if( !zVecIsEqual( x, ans, TOL ) ){
      count_ri++;
    }
    /* Gauss-Seidel's method */
    zLESolveGS( a, b, x );
    if( !zVecIsEqual( x, ans, TOL ) ){
      count_gs++;
    }
  }
  printf( "matrix size: %d x %d\n", SIZE, SIZE );
  printf( "failure rate: Cramel %d/%d, Gauss %d/%d, LU decomp %d/%d, Residual Iter. %d/%d, Gauss-Seidel %d/%d\n", count_cramel, N, count_gauss, N, count_lu, N, count_ri, N, count_gs, N );
  zMatFreeAO( 4, a, l, u, aa );
  zVecFreeAO( 3, b, x, ans );
  zIndexFree( index );
  return 0;
}
