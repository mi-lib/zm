#include <zm/zm_mca.h>

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
      r = zRandF(0,rmax);
      t = zRandF(0,zPIx2);
      zVecSetElem( vc, 0, xc + r * cos(t) );
      zVecSetElem( vc, 1, yc + r * sin(t) );
      zVecListInsertHead( vl, vc );
    }
  }
  zVecFree( vc );
}

#define NP 1000
#define NC 6

int main(int argc, char *argv[])
{
  zMCluster mc;
  zVecList points;
  int np, nc;
  double score;

  zRandInit();
  nc = argc > 1 ? atoi( argv[1] ) : NC;
  np = argc > 2 ? atoi( argv[2] ) : NP;
  gen_vec( &points, np, nc, 0, 0, 10, 10 );
  zMClusterInit( &mc, 2 );
  printf( "K-means completed in %d times of iteration.\n", zMClusterKMeans( &mc, &points, nc ) );
  zMClusterDataPrintFile( &mc, "" );
  score = zMClusterMeanSilhouette( &mc );
  zMClusterSilhouettePrintFile( &mc, "s" );
  printf( "mean silhouette = %.10g\n", score );

  zMClusterDestroy( &mc );
  zVecListDestroy( &points );
  return 0;
}
