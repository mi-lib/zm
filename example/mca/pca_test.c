#include <zm/zm_mca.h>

void gen_vec(zVecList *vl, int n, double xc, double yc, double r1, double r2, double theta)
{
  register int i;
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
    zVecElem(vc,0) = xc + xt * c - yt * s;
    zVecElem(vc,1) = yc + xt * s + yt * c;
    zVecListInsertHead( vl, vc, false );
  }
}

#define N   1000
#define DIM    2
#define CR     0.9

int main(int argc, char *argv[])
{
  zVecList points;
  zVecListCell *pc;
  zVec p, ave, score;
  zMat loading;
  double cr;
  int n;
  FILE *fp;

  zRandInit();
  n = argc > 1 ? atoi( argv[1] ) : N;
  cr = argc > 2 ? atof( argv[2] ) : CR;
  gen_vec( &points, n, 10, 10, 10, 5, zDeg2Rad(30) );
  fp = fopen( "s", "w" );
  zVecListFWrite( fp, &points );
  fclose( fp );

  p = zVecAlloc( DIM );
  ave = zVecAlloc( DIM );
  score = zVecAlloc( DIM );
  loading = zMatAllocSqr( DIM );
  printf( "executing PCA. the number of PC = %d\n",
    zPCA( &points, cr, ave, score, loading ) );

  fp = fopen( "a", "w" );
  zListForEach( &points, pc ){
    zVecSub( pc->data, ave, p );
    zMulMatTVec( loading, p, score );
    zVecDataFWrite( fp, score );
  }
  fclose( fp );

  zVecListDestroy( &points, true );
  zVecFreeAO( 3, p, ave, score );
  zMatFree( loading );
  return 0;
}
