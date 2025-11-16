#include <zm/zm_ode.h>

/* sample: unit circle function */
zVec ddp(double t, zVec x, zVec dx, void *dummy, zVec ddx)
{
  zVecSetElem( ddx, 0,-zVecElem(x,0) );
  return ddx;
}

#define DT 0.1
#define T  1000.0

#define N_METHOD 3

void output(FILE *fp, zVec x[], zVec dx[])
{
  int i;

  for( i=0; i<N_METHOD; i++ )
    fprintf( fp, "%g %g ", zVecElem(x[i],0), zVecElem(dx[i],0) );
  fprintf( fp, "\n" );
}

int main(void)
{
  zODE2 ode[N_METHOD];
  zVec x[N_METHOD], dx[N_METHOD];
  double t;
  int i;
  FILE *fp;

  zODE2Assign( &ode[0], Regular, NULL, NULL, NULL, NULL );  zODE2AssignRegular( &ode[0], RK4 );
  zODE2Assign( &ode[1], Leapfrog, NULL, NULL, NULL, NULL );
  zODE2Assign( &ode[2], Sympl, NULL, NULL, NULL, NULL );
  for( i=0; i<N_METHOD; i++ ){
    zODE2Init( &ode[i], 1, 0, ddp );
    x[i] = zVecCreateList( 1, 1.0 );
    dx[i] = zVecCreateList( 1, 0.0 );
  }
  zODE2InitHistLeapfrog( &ode[1], x[1], dx[1], DT );

  fp = fopen( "result.dat", "w" );
  output( fp, x, dx );
  for( t=0; t<T; t+=DT ){
    for( i=0; i<N_METHOD; i++ )
      zODE2Update( &ode[i], t, x[i], dx[i], DT, NULL );
    output( fp, x, dx );
  }
  for( i=0; i<N_METHOD; i++ ){
    zVecFree( x[i] );
    zVecFree( dx[i] );
    zODE2Destroy( &ode[i] );
  }
  eprintf( "[methods] 0:RK4, 1:Leapfrog, 2:Symplectic\n" );
  return 0;
}
