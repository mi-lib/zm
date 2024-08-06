#include <zm/zm_le.h>

#define N 5

int main(void)
{
  zMat a, b, m;

  a = zMatAllocSqr( N );
  b = zMatAllocSqr( N );
  m = zMatAllocSqr( N );
  zMatRandUniform( a, -1e10, 1e10 );
  zMatInv( a, b );
  zMulMatMat( a, b, m );
  zMatPrint( m );

  zMatInvHotelling( a, b, zTOL, 0 );
  zMulMatMat( a, b, m );
  zMatPrint( m );

  zMatFreeAtOnce( 3, a, b, m );
  return 0;
}
