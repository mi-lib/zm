#include <zm/zm_nurbs.h>

#define STEP 100

void test_weight(zNURBS *nurbs, int w)
{
  double t;
  int i;
  FILE *fp;
  char filename[BUFSIZ];
  zVec v;

  sprintf( filename, "w%d", w );
  fp = fopen( filename, "w" );
  zNURBSSetWeight( nurbs, 4, (double)w );
  v = zVecAlloc( zVecSizeNC( zNURBSCP(nurbs,0) ) );
  zNURBSSetSlice( nurbs, STEP );
  for( i=0; i<=zNURBSSlice(nurbs); i++ ){
    t = zNURBSKnotSlice( nurbs, i );
    if( zNURBSVec( nurbs, t, v ) )
      zVecValueFPrint( fp, v );
  }
  zVecFree( v );
  fclose( fp );
}

void output_src(zSeq *seq)
{
  FILE *fp;
  zSeqCell *cp;

  fp = fopen( "src", "w" );
  zListForEach( seq, cp )
    fprintf( fp, "%f %f\n", zVecElemNC(cp->data.v,0), zVecElemNC(cp->data.v,1) );
  fclose( fp );
}

#define ORDER 3

int main(int argc, char *argv[])
{
  zNURBS nurbs;
  zSeq seq;
  zVec v;
  int order, num, i;
  /* example data array */
  double xp[] = { 2.0, 3.0, 5.0, 4.0, 5.0, 7.0 };
  double yp[] = { 3.0,-1.0,-2.0, 0.0, 4.0, 1.5 };

  order = argc > 1 ? atoi(argv[1]) : ORDER;
  /* creation of x-values and y-values vector */
  num = sizeof(xp) / sizeof(double);
  zListInit( &seq );
  for( i=0; i<num; i++ ){
    v = zVecCreateList( 2, xp[i], yp[i] );
    zSeqEnqueue( &seq, v, 1.0 );
  }
  output_src( &seq );

  /* creation of spline interpolator */
  if( zNURBSCreate( &nurbs, &seq, order ) ){
    for( i=0; i<5; i++ )
      test_weight( &nurbs, i );
    zNURBSDestroy( &nurbs );
  }
  zSeqFree( &seq );
  return 0;  
}
