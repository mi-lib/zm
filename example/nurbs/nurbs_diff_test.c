/* 2014. 8.25 Thanks to Mr. Ken'ya Tanaka */
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
  zNURBS nurbs;
  zSeq seq;
  double xp[] = { 0.0, 2.0,-2.0, 0.0, 2.0,-2.0, 0.0 };
  double yp[] = { 0.0, 2.0, 3.0, 6.0, 7.5, 9.0,10.0 };
  double weight[] = { 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0 };
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
  zNURBSCreate( &nurbs, &seq, DIM, STEP );
  zSeqFree( &seq );
  for( i=0; i<num ; i++ )
    zNURBSSetWeight( &nurbs, i, weight[i] );
  for( i=0; i<zNURBSKnotNum(&nurbs); i++ )
    zNURBSSetKnot( &nurbs, i, knot[i] );

  v = zVecAlloc( 2 );
  for( i=0; i<=STEP; i++ ){
    t = zNURBSKnotSlice( &nurbs, i );
    printf( "%1.2f ", t );
    zNURBSVec( &nurbs, t, v );
    printf( "%g %g ", zVecElemNC(v,0), zVecElemNC(v,1) );
    zNURBSVecDiff( &nurbs, t, 0, v );
    printf( "%g %g ", zVecElemNC(v,0), zVecElemNC(v,1) );
    zNURBSVecDiff( &nurbs, t, 1, v );
    printf( "%g %g ", zVecElemNC(v,0), zVecElemNC(v,1) );
    zNURBSVecDiff( &nurbs, t, 2, v );
    printf( "%g %g\n", zVecElemNC(v,0), zVecElemNC(v,1) );
  }
  zVecFree( v );
  zNURBSDestroy( &nurbs );
  return 0;
}
