#include <zm/zm_nurbs.h>

#define STEP 100

void test_weight(zNURBS1 *nurbs, int w)
{
  double t;
  int i;
  FILE *fp;
  char filename[BUFSIZ];
  zVec v;

  sprintf( filename, "w%d", w );
  fp = fopen( filename, "w" );
  zNURBS1Weight(nurbs,4) = (double)w;
  v = zVecAlloc( zVecSizeNC( zNURBS1CP(nurbs,0) ) );
  for( i=0; i<=STEP; i++ ){
    t = ( zNURBS1KnotE(nurbs)-zNURBS1Knot0(nurbs) ) * i / STEP + zNURBS1Knot0(nurbs);
    if( zNURBS1Vec( nurbs, t, v ) )
      zVecDataFWrite( fp, v );
  }
  zVecFree( v );
  fclose( fp );
}

void output_src(zSeq *seq)
{
  FILE *fp;
  zSeqListCell *cp;

  fp = fopen( "src", "w" );
  zListForEach( seq, cp )
    fprintf( fp, "%f %f\n", zVecElem(cp->data.v,0), zVecElem(cp->data.v,1) );
  fclose( fp );
}

#define DIM 3

int main(int argc, char *argv[])
{
  zNURBS1 nurbs;
  zSeq seq;
  zVec v;
  int dim, num, i;
  /* example data array */
  double xp[] = { 2.0, 3.0, 5.0, 4.0, 5.0, 7.0 };
  double yp[] = { 3.0,-1.0,-2.0, 0.0, 4.0, 1.5 };

  dim = argc > 1 ? atoi(argv[1]) : DIM;
  /* creation of x-values and y-values vector */
  num = sizeof(xp) / sizeof(double);
  zListInit( &seq );
  for( i=0; i<num; i++ ){
    v = zVecCreateList( 2, xp[i], yp[i] );
    zSeqEnqueue( &seq, v, 1.0 );
  }
  output_src( &seq );

  /* creation of spline interpolator */
  if( zNURBS1Create( &nurbs, &seq, dim ) ){
    zNURBS1KnotFWrite( stdout, &nurbs );
    for( i=0; i<5; i++ )
      test_weight( &nurbs, i );
    zNURBS1Destroy( &nurbs );
  }
  zSeqFree( &seq );
  return 0;  
}
