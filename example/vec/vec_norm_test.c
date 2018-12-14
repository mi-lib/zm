#include <zm/zm_vec.h>

#define N 5

int main(void)
{
  zVec v;

  zRandInit();
  v = zVecAlloc( N );
  zVecRandUniform( v, -10, 10 );
  zVecWrite( v );
  printf( "||v||_2 = %g\n", zVecNorm(v) );
  printf( "||v||_inf = %g\n", zVecInfNorm(v) );
  zVecFree( v );
  return 0;
}
