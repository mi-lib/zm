#include <zm/zm_ode.h>

/* sample: unit circle function */
zVec dp(double t, zVec p, void *dummy, zVec v)
{
  zVecSetElem( v, 0,-zVecElem(p,1) );
  zVecSetElem( v, 1, zVecElem(p,0) );
  return v;
}

#define DT 0.5
#define T  10.0

int main(void)
{
  zODE ode;
  zVec x;
  double t;

  zODEAssign( &ode, Adams, NULL, NULL );
  zODEInit( &ode, 2, 4, dp );
  x = zVecCreateList( 2, 1.0, 0.0 );
  zVecDataWrite( x );
  for( t=0; t<T; t+=DT ){
    zODEUpdate( &ode, t, x, DT, NULL );
    zVecDataWrite( x );
  }
  zVecFree( x );
  zODEDestroy( &ode );
  return 0;
}
