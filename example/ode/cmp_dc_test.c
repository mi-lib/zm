#include <zm/zm_ode.h>

/* sample: unit circle function */
zVec dp(double t, zVec p, void *dummy, zVec v)
{
  zVecSetElem( v, 0,-zVecElem(p,1) );
  zVecSetElem( v, 1, zVecElem(p,0) );
  return v;
}

#define DT 0.5
#define T  100.0

#define ODES 11

void output_all(zVec *x)
{
  int i;

  for( i=0; i<ODES; i++ )
    printf( "%f ", zVecNorm(x[i])-1 );
  zEndl();
}

int main(void)
{
  zODE ode[ODES];
  zVec x[ODES];
  double t;
  int i;

  zODEAssign( &ode[0], Euler, NULL, NULL );
  zODEAssign( &ode[1], Heun, NULL, NULL );
  zODEAssign( &ode[2], RK4, NULL, NULL );
  zODEAssign( &ode[3], RKG, NULL, NULL );
  zODEAssign( &ode[4], RKF45, NULL, NULL );
  zODEAssign( &ode[5], Adams, NULL, NULL );
  zODEAssign( &ode[6], BEuler, NULL, NULL );
  zODEAssign( &ode[7], Gauss, NULL, NULL );
  zODEAssign( &ode[8], Radau, NULL, NULL );
  zODEAssign( &ode[9], Gear, NULL, NULL );
  zODEAssign( &ode[10],TR, NULL, NULL );

  zODEInitDC( &ode[0], 2, 0, dp ); /* Euler method */
  zODEInitDC( &ode[1], 2, 0, dp ); /* Heun method */
  zODEInitDC( &ode[2], 2, 0, dp ); /* classical Runge-Kutta method */
  zODEInitDC( &ode[3], 2, 0, dp ); /* Runge-Kutta-Gill method */
  zODEInitDC( &ode[4], 2, 0, dp ); /* Runge-Kutta_Fehlberg method */
  zODEInitDC( &ode[5], 2, 5, dp ); /* PC on AB/AM formula */
  zODEInitDC( &ode[6], 2, 0, dp ); /* backward Euler method */
  zODEInitDC( &ode[7], 2, 0, dp ); /* Gauss method */
  zODEInitDC( &ode[8], 2, 0, dp ); /* Radau method */
  zODEInitDC( &ode[9], 2, 3, dp ); /* Gear method */
  zODEInitDC( &ode[10],2, 0, dp ); /* trapezoidal method */

  x[0] = zVecCreateList( 2, 1.0, 0.0 );
  for( i=1; i<ODES; i++ )
    x[i] = zVecClone( x[0] );
  zODEInitHist_Gear( &ode[9], x[9] );

  output_all( x );
  for( t=0; t<T; t+=DT ){
    eprintf( "%f/%f\n", t, T );
    for( i=0; i<ODES; i++ )
      zODEUpdateDC( &ode[i], t, x[i], DT, NULL );
    output_all( x );
  }

  eprintf( "Euler: err= %f\n", zVecNorm(x[0])-1 );
  eprintf( "Heun : err= %f\n", zVecNorm(x[1])-1 );
  eprintf( "RK4  : err= %f\n", zVecNorm(x[2])-1 );
  eprintf( "RKG  : err= %f\n", zVecNorm(x[3])-1 );
  eprintf( "RKF45: err= %f\n", zVecNorm(x[4])-1 );
  eprintf( "Adams: err= %f\n", zVecNorm(x[5])-1 );
  eprintf( "BEuler:err= %f\n", zVecNorm(x[6])-1 );
  eprintf( "Gauss: err= %f\n", zVecNorm(x[7])-1 );
  eprintf( "Radau: err= %f\n", zVecNorm(x[8])-1 );
  eprintf( "Gear : err= %f\n", zVecNorm(x[9])-1 );
  eprintf( "TR   : err= %f\n", zVecNorm(x[10])-1 );
  for( i=0; i<ODES; i++ ){
    zVecFree( x[i] );
    zODEDestroyDC( &ode[i] );
  }
  return 0;
}
