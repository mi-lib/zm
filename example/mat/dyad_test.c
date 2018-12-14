#include <zm/zm_mat.h>

int main(void)
{
  zVec v1, v2;
  zMat m, dyad;

  v1 = zVecCreateList( 4, 1.0, 2.0, -2.0, -1.0 );
  v2 = zVecCreateList( 3, 1.0, 2.0, 3.0 );
  dyad = zMatAlloc( zVecSize(v1), zVecSize(v2) );
  m = zMatAlloc( zVecSize(v1), zVecSize(v2) );
  zVecWrite( v1 );
  zVecWrite( v2 );

  zVecDyad( v1, v2, dyad );
  zMatWrite( dyad );
  zMatClear( m );
  zMatAddDyad( m, v1, v2 );
  zMatWrite( m );
  zMatClear( m );
  zMatCatDyad( m, 2, v1, v2 );
  zMatWrite( m );

  zVecFree( v1 );
  zVecFree( v2 );
  zMatFree( m );
  zMatFree( dyad );
  return 0;
}
