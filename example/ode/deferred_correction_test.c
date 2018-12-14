#include <zm/zm_ode.h>

/* sample: unit circle function */
zVec dp(double t, zVec p, void *dummy, zVec v)
{
  zVecSetElem( v, 0,-zVecElem(p,1) );
  zVecSetElem( v, 1, zVecElem(p,0) );
  return v;
}

void output(zVec x1, zVec x2)
{
  printf( "%.10g %.10g %.10g %.10g\n",
    zVecElem(x1,0), zVecElem(x1,1), zVecElem(x2,0), zVecElem(x2,1) );
}

#define DT 0.5
#define T  10.0

int main(void)
{
  zODE ode, ode_dc;
  zVec x, x_dc;
  double t;

  zODEAssign( &ode, Radau, NULL, NULL );
  zODEAssign( &ode_dc, Radau, NULL, NULL );
  zODEInit( &ode, 2, 0, dp );
  zODEInitDC( &ode_dc, 2, 0, dp );
  x = zVecCreateList( 2, 1.0, 0.0 );
  x_dc = zVecClone( x );
  output( x, x_dc );
  for( t=0; t<T; t+=DT ){
    zODEUpdate( &ode, t, x, DT, NULL );
    zODEUpdateDC( &ode_dc, t, x_dc, DT, NULL );
    output( x, x_dc );
  }
  zVecFree( x );
  zVecFree( x_dc );
  zODEDestroy( &ode );
  zODEDestroyDC( &ode_dc );
  return 0;
}
