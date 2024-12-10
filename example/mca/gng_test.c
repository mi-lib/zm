#include <zm/zm_gng.h>

void gen_vec(zVecList *vl, int np, int nc, double xmin, double ymin, double xmax, double ymax, double rmax)
{
  int i, j;
  zVec vc;
  double xc, yc, r, s, c;

  zListInit( vl );
  vc = zVecAlloc( 2 );
  for( i=0; i<nc; i++ ){
    xc = zRandF(xmin,xmax);
    yc = zRandF(ymin,ymax);
    for( j=0; j<np; j++ ){
      r = zRandF(0,rmax);
      zSinCos( zRandF(-zPI,zPI), &s, &c );
      zVecSetElemList( vc, xc + r * c, yc + r * s );
      zVecListInsertHead( vl, vc );
    }
  }
  zVecFree( vc );
}

void vec_output(zVecList *points)
{
  FILE *fp;
  zVecListCell *vp;

  fp = fopen( "src", "w" );
  zListForEach( points, vp )
    zVecValueFPrint( fp, vp->data );
  fclose( fp );
}

#define NP 300
#define NC 6
#define R  1.5

int main(int argc, char *argv[])
{
  zGNG gng;
  zVecList points;
  int np, nc, i;
  FILE *fp1, *fp2;
  char filename[BUFSIZ];

  zRandInit();
  np = argc > 1 ? atoi( argv[1] ) : NP;
  nc = argc > 2 ? atoi( argv[2] ) : NC;
  gen_vec( &points, np, nc, R, R, 10-R, 10-R, R );
  vec_output( &points );

  zGNGInit( &gng, 2, NULL, &points );
  zGNGSetBatchTrialSize( &gng, 10 );

  fp2 = fopen( "n", "w" );
  for( i=0; i<zListSize(&points); i++ ){
    zGNGUpdate( &gng, &points );
    if( i % 10 == 0 ){
      sprintf( filename, "%04d", i );
      fp1 = fopen( filename, "w" );
      zGNGFWrite( fp1, &gng );
      fclose( fp1 );
    }
    fprintf( fp2, "%d %d\n", zListSize(&gng.unitlist), zListSize(&gng.edgelist) );
  }
  fclose( fp2 );

  zGNGDestroy( &gng );
  zVecListDestroy( &points );
  return 0;
}
