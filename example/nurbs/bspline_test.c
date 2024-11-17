#include <zm/zm_nurbs.h>

#define STEP 100
#define DIM  3

void output_src(zSeq *seq)
{
  FILE *fp;
  zSeqCell *cp;

  fp = fopen( "src", "w" );
  zListForEach( seq, cp )
    fprintf( fp, "%f %f\n", zVecElemNC(cp->data.v,0), zVecElemNC(cp->data.v,1) );
  fclose( fp );
}

int main(int argc, char *argv[])
{
  zBSpline bspline;
  zSeq seq;
  double xp[] = { 0.0, 2.0,-2.0, 0.0, 2.0,-2.0, 0.0 };
  double yp[] = { 0.0, 2.0, 3.0, 6.0, 7.5, 9.0,10.0 };
  double knot[] = { 0, 0, 0, 0, 1, 2, 3, 4, 4, 4, 4 };
  zVec v;
  double t;
  int num, i;

  num = sizeof(xp) / sizeof(double);
  zListInit( &seq );
  for( i=num-1; i>=0; i-- ){
    v = zVecCreateList( 2, xp[i], yp[i] );
    zSeqEnqueue( &seq, v, 1.0 );
  }
  output_src( &seq );
  zBSplineCreate( &bspline, &seq, DIM );
  zSeqFree( &seq );
  for( i=0; i<zBSplineKnotNum(&bspline); i++ )
    zBSplineSetKnot( &bspline, i, knot[i] );

  v = zVecAlloc( 2 );
  zBSplineSetSlice( &bspline, STEP );
  for( i=0; i<=zBSplineSlice(&bspline); i++ ){
    t = zBSplineKnotSlice( &bspline, i );
    printf( "%1.2f ", t );
    zBSplineVec( &bspline, t, v );
    printf( "%g %g ", zVecElemNC(v,0), zVecElemNC(v,1) );
    zBSplineVecDiff( &bspline, t, 0, v );
    printf( "%g %g ", zVecElemNC(v,0), zVecElemNC(v,1) );
    zBSplineVecDiff( &bspline, t, 1, v );
    printf( "%g %g ", zVecElemNC(v,0), zVecElemNC(v,1) );
    zBSplineVecDiff( &bspline, t, 2, v );
    printf( "%g %g\n", zVecElemNC(v,0), zVecElemNC(v,1) );
  }
  zVecFree( v );
  zBSplineDestroy( &bspline );
  return 0;
}
