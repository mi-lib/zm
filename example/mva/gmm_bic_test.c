#include <zm/zm_mva.h>
#include <zm/zm_rand.h>

void gen_vec(zVecList *vl, int np, int nc, double xmin, double ymin, double xmax, double ymax)
{
  int i, j;
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
      r = zRandND( NULL, 0, 0.5*rmax );
      t = zRandF( 0, zPIx2 );
      zVecSetElemNC( vc, 0, xc + r * cos(t) );
      zVecSetElemNC( vc, 1, yc + r * sin(t) );
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
  zListForEach( points, vp ){
    fprintf( fp, "%g %g 0\n", zVecElem(vp->data,0), zVecElem(vp->data,1) );
  }
  fclose( fp );
}

void gmm_output(zGMM *gmm)
{
  FILE *fp;
  char filename[BUFSIZ];
  zGMMListCell *gc;
  double x0, y0;

  sprintf( filename, "gmm%d", zListSize(&gmm->glist) );
  fp = fopen( filename, "w" );
  fprintf( fp, "set isosamples 30\n" );
  fprintf( fp, "set contour base\n" );
  fprintf( fp, "unset key\n" );
  fprintf( fp, "splot 'src', " );
  zListForEach( &gmm->glist, gc ){
    x0 = zVecElem(gc->data.core,0);
    y0 = zVecElem(gc->data.core,1);
    fprintf( fp, " +%g/sqrt(2*pi*%g)*exp(-0.5*(%g*(x-%g)**2+2*%g*(x-%g)*(y-%g)+%g*(y-%g)**2))",
      gc->data.weight,
      gc->data._cov_det,
      zMatElem(gc->data._cov_inv,0,0),x0,
      zMatElem(gc->data._cov_inv,0,1),x0,y0,
      zMatElem(gc->data._cov_inv,1,1),y0 );
  }
  fprintf( fp, "\n" );
  fclose( fp );
}

#define NP 1000
#define NC 4

int main(int argc, char *argv[])
{
  zGMM gmm;
  zVecList points;
  int np, nc, k, i;

  zRandInit();
  np = argc > 1 ? atoi( argv[1] ) : NP;
  nc = argc > 2 ? atoi( argv[2] ) : NC;
  gen_vec( &points, np, nc, 0, 0, 10, 10 );
  vec_output( &points );

  for( i=0; i<6; i++ ){
    k = nc + i - 2;
    zGMMInit( &gmm, k, 2 );
    zGMMCreateEM( &gmm, &points );
    gmm_output( &gmm );
    printf( "%d: BIC = %g, AIC = %g\n", k, zGMMBIC( &gmm, &points ), zGMMAIC( &gmm ) );
    zGMMDestroy( &gmm );
  }
  zVecListDestroy( &points );
  return 0;
}
