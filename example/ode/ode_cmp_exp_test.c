#include <zm/zm_ode.h>

#define A  100.0
#define A0 10.0

/* sample: exponential function */
zVec dp(double t, zVec p, void *dummy, zVec dp)
{
  zVecSetElem( dp, 0, -A*zVecElem(p,0) );
  return dp;
}

#define DT   0.01
#define STEP 20

#define N_METHOD 13

void output(FILE *fp, double t, zVec *x)
{
  int i;
  double a;

  a = A0*exp(-A*t);
  fprintf( fp, "%g ", a );
  for( i=0; i<N_METHOD; i++ )
    fprintf( fp, "%g %g ", zVecElem(x[i],0), zVecElem(x[i],0)-a );
  fprintf( fp, "\n" );
}

int main(void)
{
  zODE ode[N_METHOD];
  zVec x[N_METHOD];
  double t;
  int i, j;
  FILE *fp;

  /* explicit Runge-Kutta methods */
  zODEAssign( &ode[0], Euler,  NULL, NULL ); zODEInit( &ode[0], 1, 0, dp ); /* Euler method */
  zODEAssign( &ode[1], Heun,   NULL, NULL ); zODEInit( &ode[1], 1, 0, dp ); /* Heun method */
  zODEAssign( &ode[2], RK4,    NULL, NULL ); zODEInit( &ode[2], 1, 0, dp ); /* classical Runge-Kutta method */
  zODEAssign( &ode[3], RKG,    NULL, NULL ); zODEInit( &ode[3], 1, 0, dp ); /* Runge-Kutta-Gill method */
  /* embedded Runge-Kutta methods */
  zODEAssign( &ode[4], RKF45,  NULL, NULL ); zODEInit( &ode[4], 1, 0, dp ); /* Runge-Kutta_Fehlberg method */
  zODEAssign( &ode[5], CK45,   NULL, NULL ); zODEInit( &ode[5], 1, 0, dp ); /* Cash-Karp method */
  zODEAssign( &ode[6], DP45,   NULL, NULL ); zODEInit( &ode[6], 1, 0, dp ); /* Dormand-Prince method */
  /* predicter-corrector methods */
  zODEAssign( &ode[7], Adams,  NULL, NULL ); zODEInit( &ode[7], 1, 5, dp ); /* PC on AB/AM formula */
  /* implicit Runge-Kutta methods */
  zODEAssign( &ode[8], BEuler, NULL, NULL ); zODEInit( &ode[8], 1, 0, dp ); /* backward Euler method */
  zODEAssign( &ode[9], Gauss,  NULL, NULL ); zODEInit( &ode[9], 1, 0, dp ); /* Gauss method */
  zODEAssign( &ode[10],Radau,  NULL, NULL ); zODEInit( &ode[10],1, 0, dp ); /* Radau method */
  /* other implicit methods */
  zODEAssign( &ode[11],Gear,   NULL, NULL ); zODEInit( &ode[11],1, 3, dp ); /* Gear method */
  zODEAssign( &ode[12],TR,     NULL, NULL ); zODEInit( &ode[12],1, 0, dp ); /* trapezoidal method */

  x[0] = zVecCreateList( 1, A0 );
  for( i=1; i<N_METHOD; i++ )
    x[i] = zVecClone( x[0] );
  zODEInitHist_Gear( &ode[11], x[11] );

  fp = fopen( "result.dat", "w" );
  output( fp, 0, x );
  for( i=0; i<STEP; i++ ){
    t = i * DT;
    for( j=0; j<N_METHOD; j++ )
      zODEUpdate( &ode[j], t, x[j], DT, NULL );
    output( fp, t+DT, x );
  }
  fclose( fp );
  for( i=0; i<N_METHOD; i++ ){
    zVecFree( x[i] );
    zODEDestroy( &ode[i] );
  }
  eprintf( "[methods] 0:Euler, 1:Heun, 2:RK4, 3:RKG, 4:RKF45, 5:CK45, 6:DP45, 7:Adams, 8:Backward-Euler, 9:Gauss(BK4), 10:Radau, 11:Gear, 12:Trapezoidal\n" );
  return 0;
}
