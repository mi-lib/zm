#include <zm/zm_mat.h>

int main(void)
{
  zVec v;
  zMat m;

  v = zVecAlloc( 5 );
  zVecSetElemList( v, 1.0, 2.0, 3.0, 4.0, 5.0 );
  zVecWrite( v );
  zVecFree( v );

  v = zVecCreateList( 5,-1.0,-2.0,-3.0,-4.0,-5.0 );
  zVecWrite( v );
  zVecFree( v );

  m = zMatAlloc( 2, 3 );
  zMatSetElemList( m, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0 );
  zMatWrite( m );
  zMatFree( m );

  m = zMatCreateList( 2, 3, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0 );
  zMatWrite( m );
  zMatFree( m );
  return 0;
}
