#include <zm/zm_nurbs.h>

#define STEP 100
#define N     50

void nn_test(zNURBS *nurbs)
{
  register int i;
  FILE *fp;
  zVec v, n, nn;

  v = zVecAlloc( zVecSizeNC( zNURBSCP(nurbs,0) ) );
  n = zVecAlloc( zVecSizeNC( zNURBSCP(nurbs,0) ) );
  nn = zVecAlloc( zVecSizeNC( zNURBSCP(nurbs,0) ) );

  fp = fopen( "p", "w" );
  for( i=0; i<=STEP; i++ ){
    if( zNURBSVec( nurbs, (double)i/STEP, v ) )
      zVecDataFPrint( fp, v );
  }
  fclose( fp );
  fp = fopen( "nn", "w" );
  for( i=0; i<N; i++ ){
    zVecSetElemList( v, zRandF(-1,7), zRandF(-3,5) );
    zNURBSVecNN( nurbs, v, nn );
    zVecDataFPrint( fp, v );
    zVecDataFPrint( fp, nn );
    fprintf( fp, "\n" );
  }
  fclose( fp );
  zVecFreeAtOnce( 3, v, n, nn );
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

#define DIM 3

int main(int argc, char *argv[])
{
  zNURBS nurbs;
  zSeq seq;
  zVec v;
  int num, i;
  /* example data array */
  double xp[] = { 2.0, 3.0, 5.0, 4.0, 5.0, 7.0 };
  double yp[] = { 3.0,-1.0,-2.0, 0.0, 4.0, 1.5 };

  zRandInit();
  /* creation of x-values and y-values vector */
  num = sizeof(xp) / sizeof(double);
  zListInit( &seq );
  for( i=0; i<num; i++ ){
    v = zVecCreateList( 2, xp[i], yp[i] );
    zSeqEnqueue( &seq, v, 1.0 );
  }
  output_src( &seq );

  /* creation of spline interpolator */
  if( zNURBSCreate( &nurbs, &seq, DIM, 0 ) ){
    zNURBSKnotNormalize( &nurbs );
    nn_test( &nurbs );
    zNURBSDestroy( &nurbs );
  }
  zSeqFree( &seq );
  return 0;  
}
