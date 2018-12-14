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

#define ODES 10

void output(double t, zVec *x)
{
  register int i;
  double a;

  a = A0*exp(-A*t);
#if 0
  /* output error */
  for( i=0; i<ODES; i++ )
    printf( "%f ", zVecElem(x[i],0)-a );
#else
  /* output value */
  printf( "%f ", a );
  for( i=0; i<ODES; i++ )
    printf( "%f ", zVecElem(x[i],0) );
#endif
  zEndl();
}

int main(void)
{
  zODE ode[ODES];
  zVec x[ODES];
  double t;
  register int i, j;

  zODEAssign(&ode[0],Euler,NULL,NULL);  zODEInit(&ode[0],1,0,dp);
  zODEAssign(&ode[1],Heun,NULL,NULL);   zODEInit(&ode[1],1,0,dp);
  zODEAssign(&ode[2],RK4,NULL,NULL);    zODEInit(&ode[2],1,0,dp);
  zODEAssign(&ode[3],RKG,NULL,NULL);    zODEInit(&ode[3],1,0,dp);
  zODEAssign(&ode[4],RKF45,NULL,NULL);  zODEInit(&ode[4],1,0,dp);
  zODEAssign(&ode[5],Adams,NULL,NULL);  zODEInit(&ode[5],1,3,dp);
  zODEAssign(&ode[6],BEuler,NULL,NULL); zODEInit(&ode[6],1,0,dp);
  zODEAssign(&ode[7],Gauss,NULL,NULL);  zODEInit(&ode[7],1,0,dp);
  zODEAssign(&ode[8],Radau,NULL,NULL);  zODEInit(&ode[8],1,0,dp);
  zODEAssign(&ode[9],Gear,NULL,NULL);   zODEInit(&ode[9],1,3,dp);

  x[0] = zVecCreateList( 1, A0 );
  for( i=1; i<ODES; i++ )
    x[i] = zVecClone( x[0] );
  zODEInitHist_Gear( &ode[9], x[9] );

  output( 0, x );
  for( i=0; i<STEP; i++ ){
    t = i * DT;
    for( j=0; j<ODES; j++ )
      zODEUpdate( &ode[j], t, x[j], DT, NULL );
    output( t+DT, x );
  }
  for( i=0; i<ODES; i++ ){
    zVecFree( x[i] );
    zODEDestroy( &ode[i] );
  }
  return 0;
}
