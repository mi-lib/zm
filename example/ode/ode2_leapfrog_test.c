#include <zm/zm_ode.h>

/* sample: unit circle function */
zVec ddp(double t, zVec x, zVec v, void *dummy, zVec a)
{
  zVecSetElem( a, 0,-zVecElem(x,0) );
  return a;
}

#define DT 0.1
#define T  1000.0

void output2(zVec x, zVec v)
{
  printf( "%f %f\n", zVecElem(x,0), zVecElem(v,0) );
}

int main(void)
{
  zODE2 ode;
  zVec x, v;
  double t;

  zODE2Assign( &ode, Leapfrog, NULL, NULL, NULL, NULL );
  zODE2Init( &ode, 1, 0, ddp );
  x = zVecCreateList( 1, 1.0 );
  v = zVecAlloc( 1 );
  zODE2InitHistLeapfrog( &ode, x, v, DT );
  output2( x, v );
  for( t=0; t<T; t+=DT ){
    zODE2Update( &ode, t, x, v, DT, NULL );
    output2( x, v );
  }
  zVecFree( x );
  zVecFree( v );
  zODE2Destroy( &ode );
  return 0;
}
