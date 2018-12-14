#include <zm/zm_ode.h>

/* sample: unit circle function */
zVec ddp(double t, zVec x, zVec dx, void *dummy, zVec ddx)
{
  zVecSetElem( ddx, 0,-zVecElem(x,0) );
  return ddx;
}

#define DT 0.1
#define T  1000.0

void output2(zVec x, zVec dx)
{
  printf( "%f %f\n", zVecElem(x,0), zVecElem(dx,0) );
}

int main(void)
{
  zODE2 ode;
  zVec x, dx;
  double t;

  zODE2Assign( &ode, Sympl, NULL, NULL, NULL, NULL );
  zODE2Init( &ode, 1, 0, ddp );
  x = zVecCreateList( 1, 1.0 );
  dx = zVecAlloc( 1 );
  output2( x, dx );
  for( t=0; t<T; t+=DT ){
    zODE2Update( &ode, t, x, dx, DT, NULL );
    output2( x, dx );
  }
  zVecFree( x );
  zVecFree( dx );
  zODE2Destroy( &ode );
  return 0;
}
