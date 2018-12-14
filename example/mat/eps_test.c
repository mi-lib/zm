#include <zm/zm_mat.h>

int main()
{
  zVec v;
  zMat m;

  v = zVecCreateList( 3, 1.0,-2.0, 3.0 );
  m = zMatCreateList( 2, 2, 2.0,-3.0, 4.0,-5.0 );
  zVecWrite( v );
  zMatWrite( m );
  while( 1 ){
    zVecWrite( v );
    if( zVecIsTiny( v ) ) break;
    zVecMulDRC( v, 0.1 );
  }
  while( 1 ){
    zMatWrite( m );
    if( zMatIsTiny( m ) ) break;
    zMatMulDRC( m, 0.1 );
  }
  return 0;
}
