#include <zm/zm_vec.h>

int main(void)
{
  zVec data;

  data = zVecCreateList( 5, 1.0, 3.0, 4.0, 0.0, 2.0 );
  printf( "[data]\n" );
  zVecWrite( data );
  printf( "average  = %f\n", zVecAve( data ) );
  printf( "variance = %f\n", zVecVar( data ) );
  zVecFree( data );
  return 0;
}
