#include <zm/zm_ode.h>

/* sample: unit circle function */
zVec dp(double t, zVec p, void *dummy, zVec v)
{
  zVecSetElem( v, 0,-zVecElem(p,1) );
  zVecSetElem( v, 1, zVecElem(p,0) );
  return v;
}

#define DT 0.1
#define T  10.0

int main(int argc, char *argv[])
{
  zODE ode;
  zVec x;
  double t;

  if( argc > 1 ){
    switch( atoi(argv[1]) ){
    case 1:  zODEAssign( &ode, CK45, NULL, NULL ); break;
    case 2:  zODEAssign( &ode, DP45, NULL, NULL ); break;
    default: zODEAssign( &ode, RKF45, NULL, NULL ); break;
    }
  } else
    zODEAssign( &ode, RKF45, NULL, NULL );
  zODEInit( &ode, 2, 0, dp );
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
