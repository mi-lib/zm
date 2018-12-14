#include <zm/zm_pex.h>

int main(void)
{
  zPex p;
  zVec fact;

  fact = zVecCreateList( 3, 1.0, 2.0, 3.0 );
  p = zPexExp( fact );

  zVecWrite( fact );
  zPexExprX( p );

  zPexFree( p );
  zVecFree( fact );
  return 0;
}
