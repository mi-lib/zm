#include <zm/zm_vec.h>

int main(void)
{
  zVec v1, v2;

  v1 = zVecCreateList( 3, 0.25, 0.25*sqrt(3.0), 0.5 );
  v2 = zVecCreateList( 3,-0.25,-0.25*sqrt(3.0),-0.5 );

  zVecWrite( v1 );
  zVecWrite( v2 );
  printf( "dist = %f\n", zVecDist(v1,v2) );
  return 0;
}
