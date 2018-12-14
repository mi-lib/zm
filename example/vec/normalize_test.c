#include <zm/zm_vec.h>

int main(void)
{
  double norm;
  zVec v;

  v = zVecCreateList( 5, 1.0, 2.0, 3.0, 4.0, 5.0 );
  printf( "org: " ); zVecWrite( v );
  printf( "norm = %g\n", ( norm = zVecNorm( v ) ) );
  printf( "normalized: " ); zVecWrite( zVecNormalizeDRC(v) );
  printf( "norm = %g\n", zVecNorm( v ) );
  printf( "recover: " ); zVecWrite( zVecMulDRC(v,norm) );
  return 0;
}
