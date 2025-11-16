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

#define N_METHOD 13

void output_all(FILE *fp, zVec *x)
{
  int i;

  for( i=0; i<N_METHOD; i++ )
    fprintf( fp, "%g %g %g ", zVecElemNC(x[i],0), zVecElemNC(x[i],1), zVecNorm(x[i])-1 );
  fprintf( fp, "\n" );
}

int main(int argc, char *argv[])
{
  zODE ode[N_METHOD];
  zVec x[N_METHOD];
  double t;
  int i;
  FILE *fp;

  /* explicit Runge-Kutta methods */
  zODEAssign( &ode[0], Euler,  NULL, NULL ); zODEInit( &ode[0], 2, 0, dp ); /* Euler method */
  zODEAssign( &ode[1], Heun,   NULL, NULL ); zODEInit( &ode[1], 2, 0, dp ); /* Heun method */
  zODEAssign( &ode[2], RK4,    NULL, NULL ); zODEInit( &ode[2], 2, 0, dp ); /* classical Runge-Kutta method */
  zODEAssign( &ode[3], RKG,    NULL, NULL ); zODEInit( &ode[3], 2, 0, dp ); /* Runge-Kutta-Gill method */
  /* embedded Runge-Kutta methods */
  zODEAssign( &ode[4], RKF45,  NULL, NULL ); zODEInit( &ode[4], 2, 0, dp ); /* Runge-Kutta_Fehlberg method */
  zODEAssign( &ode[5], CK45,   NULL, NULL ); zODEInit( &ode[5], 2, 0, dp ); /* Cash-Karp method */
  zODEAssign( &ode[6], DP45,   NULL, NULL ); zODEInit( &ode[6], 2, 0, dp ); /* Dormand-Prince method */
  /* predicter-corrector methods */
  zODEAssign( &ode[7], Adams,  NULL, NULL ); zODEInit( &ode[7], 2, 5, dp ); /* PC on AB/AM formula */
  /* implicit Runge-Kutta methods */
  zODEAssign( &ode[8], BEuler, NULL, NULL ); zODEInit( &ode[8], 2, 0, dp ); /* backward Euler method */
  zODEAssign( &ode[9], Gauss,  NULL, NULL ); zODEInit( &ode[9], 2, 0, dp ); /* Gauss method */
  zODEAssign( &ode[10],Radau,  NULL, NULL ); zODEInit( &ode[10],2, 0, dp ); /* Radau method */
  /* other implicit methods */
  zODEAssign( &ode[11],Gear,   NULL, NULL ); zODEInit( &ode[11],2, 3, dp ); /* Gear method */
  zODEAssign( &ode[12],TR,     NULL, NULL ); zODEInit( &ode[12],2, 0, dp ); /* trapezoidal method */

  x[0] = zVecCreateList( 2, 1.0, 0.0 );
  for( i=1; i<N_METHOD; i++ )
    x[i] = zVecClone( x[0] );
  zODEInitHist_Gear( &ode[11], x[11] );

  fp = fopen( "result.dat", "w" );
  output_all( fp, x );
  for( t=0; t<T; t+=DT ){
    for( i=0; i<N_METHOD; i++ )
      zODEUpdate( &ode[i], t, x[i], DT, NULL );
    output_all( fp, x );
  }
  fclose( fp );
  for( i=0; i<N_METHOD; i++ ){
    zVecFree( x[i] );
    zODEDestroy( &ode[i] );
  }
  eprintf( "[methods] 0:Euler, 1:Heun, 2:RK4, 3:RKG, 4:RKF45, 5:CK45, 6:DP45, 7:Adams, 8:Backward-Euler, 9:Gauss(BK4), 10:Radau, 11:Gear, 12:Trapezoidal\n" );
  return 0;
}
