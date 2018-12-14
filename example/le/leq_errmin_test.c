#include <zm/zm_le.h>

void test(zMat a, zVec b, zVec w, zVec x, zVec _b)
{
  zMatWrite( a );
  zVecWrite( b );
  zLESolveErrorMin( a, b, w, x );
  printf( "ans: " );
  zVecWrite( x );
  zMulMatVec( a, x, _b );
  zVecSub( b, _b, _b );
  printf( "err=%g\n", zVecNorm( _b ) );
}

#define TEST 1

int main(int argc, char *argv[])
{
#if TEST == 1
  double aarr[] = {
    1, 1, 1,
  };
  double barr[] = {
    3, 3, 3,
  };
  int n=3, m=1;
  double we = 100.0;
#else
  double aarr[] = {
    1, 1, 1,
  };
  double barr[] = {
    2, 2, 2,
  };
  int n=3, m=1;
  double we = 1.0;
#endif
  zMat a;
  zVec b, w, x, _b;

  a = zMatCloneArray( aarr, n, m );
  b = zVecCloneArray( barr, n );
  _b= zVecAlloc( n );
  w = zVecAlloc( n );
  x = zVecAlloc( m );
  zVecSetAll( w, we );

  test( a, b, w, x, _b );

  zMatFree( a );
  zVecFree( b );
  zVecFree( _b );
  zVecFree( w );
  zVecFree( x );
  return 0;
}
