#include <zm/zm_mca.h>

void gen_vec(zVecList *vl, int n, double xc, double yc, double r1, double r2, double theta)
{
  int i;
  zVec vc;
  double s, c, st, ct, xt, yt, d;

  zListInit( vl );
  zSinCos( theta, &s, &c );
  for( i=0; i<n; i++ ){
    d = zRandF(0.0,1.0);
    zSinCos( zRandF(-zPI,zPI), &st, &ct );
    xt = r1 * d * ct;
    yt = r2 * d * st;
    vc = zVecAlloc( 2 );
    zVecSetElem( vc, 0, xc + xt * c - yt * s );
    zVecSetElem( vc, 1, yc + xt * s + yt * c );
    zVecAddrListInsertHead( vl, vc );
  }
}

#define N   1000
#define DIM    2

int main(int argc, char *argv[])
{
  zVecList src, dest;
  zVec p, mean;
  zMat cov;
  int n;
  FILE *fp;

  zRandInit();
  n = argc > 1 ? atoi( argv[1] ) : N;
  gen_vec( &src, n, 10, 10, 10, 5, zDeg2Rad(30) );
  fp = fopen( "s", "w" );
  zVecListFPrint( fp, &src );
  fclose( fp );

  p = zVecAlloc( DIM );
  mean = zVecAlloc( DIM );
  cov = zMatAllocSqr( DIM );
  zVecListMeanCov( &src, mean, cov );

  zVecListGenRandND( &dest, N, mean, cov );
  fp = fopen( "d", "w" );
  zVecListFPrint( fp, &dest );
  fclose( fp );

  zVecListDestroy( &src );
  zVecListDestroy( &dest );
  zVecFree( p );
  zVecFree( mean );
  zMatFree( cov );
  return 0;
}
