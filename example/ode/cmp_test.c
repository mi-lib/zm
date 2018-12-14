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

  zODEInit( &ode[0], 2, 0, dp ); /* Euler method */
  zODEInit( &ode[1], 2, 0, dp ); /* Heun method */
  zODEInit( &ode[2], 2, 0, dp ); /* classical Runge-Kutta method */
  zODEInit( &ode[3], 2, 0, dp ); /* Runge-Kutta-Gill method */
  zODEInit( &ode[4], 2, 0, dp ); /* Runge-Kutta_Fehlberg method */
  zODEInit( &ode[5], 2, 5, dp ); /* PC on AB/AM formula */
  zODEInit( &ode[6], 2, 0, dp ); /* backward Euler method */
  zODEInit( &ode[7], 2, 0, dp ); /* Gauss method */
  zODEInit( &ode[8], 2, 0, dp ); /* Radau method */
  zODEInit( &ode[9], 2, 3, dp ); /* Gear method */
  zODEInit( &ode[10],2, 0, dp ); /* trapezoidal method */

  x[0] = zVecCreateList( 2, 1.0, 0.0 );
  for( i=1; i<ODES; i++ )
    x[i] = zVecClone( x[0] );
  zODEInitHist_Gear( &ode[9], x[9] );

  output_all( x );
  for( t=0; t<T; t+=DT ){
    for( i=0; i<ODES; i++ )
      zODEUpdate( &ode[i], t, x[i], DT, NULL );
    output_all( x );
  }
  for( i=0; i<ODES; i++ ){
    zVecFree( x[i] );
    zODEDestroy( &ode[i] );
  }
  return 0;
}
