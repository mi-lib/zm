#include <zm/zm_mca.h>
#include <zm/zm_rand.h>

void gen_vec(zVecList *vl, int np, int nc, double xmin, double ymin, double xmax, double ymax)
{
  register int i, j;
  zVec vc;
  double xc, yc, rmax;
  double r, t;

  zListInit( vl );
  vc = zVecAlloc( 2 );
  rmax = sqrt( zSqr( xmax - xmin ) + zSqr( ymax - ymin ) ) / ( 2 * nc );
  for( i=0; i<nc; i++ ){
    xc = zRandF( xmin, xmax );
    yc = zRandF( ymin, ymax );
    for( j=0; j<np; j++ ){
      r = zRandFND(0,0.5*rmax);
      t = zRandF(0,zPIx2);
      zVecElem(vc,0) = xc + r * cos(t);
      zVecElem(vc,1) = yc + r * sin(t);
      zVecListInsertHead( vl, vc, true );
    }
  }
  zVecFree( vc );
}

#define NP 1000
#define NC 10

int main(int argc, char *argv[])
{
  zMCluster mc;
  zClusterListCell *vcc;
  zVecList points;
  FILE *fp;
  char filename[BUFSIZ];
  int np, nc, i = 0;

  zRandInit();
  np = argc > 1 ? atoi( argv[1] ) : NP;
  nc = argc > 2 ? atoi( argv[2] ) : NC;
  gen_vec( &points, np, nc, 0, 0, 10, 10 );
  zMClusterInit( &mc, 2, NULL, 2, NULL );
  printf( "X-means completed in %d times of iteration.\n",
    zMClusterXMeansBIC( &mc, &points, NULL, NULL ) );

  zListForEach( &mc.cl, vcc ){
    sprintf( filename, "%d", i++ );
    fp = fopen( filename, "w" );
    zClusterDataFWrite( fp, &vcc->data );
    fclose( fp );
  }
  zMClusterDestroy( &mc );
  zVecListDestroy( &points, true );
  return 0;
}
