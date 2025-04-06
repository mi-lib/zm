#include <zm/zm_mva.h>
#include <zm/zm_rand.h>

void gen_vec(zVecList *vl, int np, int nc, double xmin, double ymin, double zmin, double xmax, double ymax, double zmax)
{
  int i, j;
  zVec vc;
  double xc, yc, zc, a, b, c, e;
  double x, y, z;

  zListInit( vl );
  vc = zVecAlloc( 3 );
  for( i=0; i<nc; i++ ){
    xc = zRandF(xmin,xmax);
    yc = zRandF(ymin,ymax);
    zc = zRandF(zmin,zmax);
    a = zRandF(-1,1);
    b = zRandF(-1,1);
    c = zc - a*xc - b*yc;
    for( j=0; j<np; j++ ){
      x = zRandF( xmin, xmax );
      y = zRandF( ymin, ymax );
      z = a*x + b*y + c;
      e = zRandND0( NULL );
      zVecSetElem( vc, 0, x + e * a );
      zVecSetElem( vc, 1, y + e * b );
      zVecSetElem( vc, 2, z + e );
      zVecListInsertHead( vl, vc );
    }
  }
  zVecFree( vc );
}

void vec_output(zVecList *points)
{
  FILE *fp;

  fp = fopen( "src", "w" );
  zVecListFPrint( fp, points );
  fclose( fp );
}

void gmm_output(zGMM *gmm)
{
  FILE *fp;
  zGMMListCell *gc;

  fp = fopen( "gmm", "w" );
  fprintf( fp, "set isosamples 20\n" );
  fprintf( fp, "unset key\n" );
  fprintf( fp, "splot 'src'" );
  zListForEach( &gmm->glist, gc ){
    fprintf( fp, ", %g*x+%g*y+%g",
      zVecElem(gc->data.core,0),
      zVecElem(gc->data.core,1),
      zVecElem(gc->data.core,2) );
  }
  fprintf( fp, "\n" );
  fclose( fp );
}

#define NP 1000
#define NC 3

int main(int argc, char *argv[])
{
  zGMM gmm;
  zVecList points;
  int np, nc;

  zRandInit();
  np = argc > 1 ? atoi( argv[1] ) : NP;
  nc = argc > 2 ? atoi( argv[2] ) : NC;
  gen_vec( &points, np, nc, 0, 0, 0, 10, 10, 10 );
  vec_output( &points );

  zGMMInit( &gmm, nc, 3 );
  zGMMSetMethodLS( &gmm, NULL );
  zGMMCreateEM( &gmm, &points );
  gmm_output( &gmm );
  zGMMDestroy( &gmm );
  zVecListDestroy( &points );
  return 0;
}
