#include <zm/zm_mca.h>

zVec errorLSM(zClusterMethod *cm, zVec p, zVec core, void *dummy, zVec err)
{
  double xm, e;

  xm = zVecElem(p,zVecSizeNC(p)-1);
  zVecSetElem( p, zVecSizeNC(p)-1, 1 );
  e = zVecInnerProd(core,p) - xm;
  zVecSetElem( p, zVecSizeNC(p)-1, xm );
  zVecSetElem( err, 0, e );
  return err;
}

zVec coreLSM(zClusterMethod *cm, zVecList *pl, void *dummy, zVec core)
{
  zVecListCell *vc;
  zMat c;
  zVec b;
  double xm;

  c = zMatAllocSqr( zVecSizeNC(zListTail(pl)->data) );
  b = zVecAlloc( zVecSizeNC(zListTail(pl)->data) );
  if( c == NULL || b == NULL ){
    ZALLOCERROR();
    core = NULL;
    goto TERMINATE;
  }
  zListForEach( pl, vc ){
    xm = zVecElem(vc->data,zVecSizeNC(vc->data)-1);
    zVecSetElem( vc->data, zVecSizeNC(vc->data)-1, 1 );
    zMatAddDyadNC( c, vc->data, vc->data );
    zVecCatNCDRC( b, xm, vc->data );
    zVecSetElem( vc->data, zVecSizeNC(vc->data)-1, xm );
  }
  zLESolveGauss( c, b, core );
 TERMINATE:
  zMatFree( c );
  zVecFree( b );
  return core;
}




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
      e = zRandF( -1, 1 );
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
  zVecListCell *vp;

  fp = fopen( "src", "w" );
  zListForEach( points, vp )
    zVecValueFPrint( fp, vp->data );
  fclose( fp );
}

void plane_output(zMCluster *mc)
{
  FILE *fp;
  zClusterListCell *c;

  fp = fopen( "mc", "w" );
  fprintf( fp, "set isosamples 20\n" );
  fprintf( fp, "unset key\n" );
  fprintf( fp, "splot 'src'" );
  zListForEach( zMClusterClusterList(mc), c ){
    fprintf( fp, ", %g*x+%g*y+%g",
      zVecElem(c->data.core,0),
      zVecElem(c->data.core,1),
      zVecElem(c->data.core,2) );
  }
  fprintf( fp, "\n" );
  fclose( fp );
}

#define NP 1000
#define NC 3

int main(int argc, char *argv[])
{
  zMCluster mc;
  zVecList points;
  int np, nc;

  zRandInit();
  nc = argc > 1 ? atoi( argv[1] ) : NC;
  np = argc > 2 ? atoi( argv[2] ) : NP;
  gen_vec( &points, np, nc, 0, 0, 0, 10, 10, 10 );
  vec_output( &points );

  zMClusterInit( &mc, 3 );
  zMClusterSetErrorFunc( &mc, 1, errorLSM, NULL );
  zMClusterSetCoreFunc( &mc, 3, coreLSM, NULL );
  printf( "K-means completed in %d times of iteration.\n", zMClusterKMeans( &mc, &points, nc ) );
  plane_output( &mc );

  zMClusterDestroy( &mc );
  zVecListDestroy( &points );
  return 0;
}
