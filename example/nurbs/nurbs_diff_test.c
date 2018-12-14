/* 2014. 8.25 Thanks to Mr. Ken'ya Tanaka */
#include <zm/zm_nurbs.h>

#define STEP 100
#define DIM  3

int main(int argc, char *argv[])
{
  zNURBS nurbs;
  zSeq seq;
  zVec v;
  double t;
  int num, i;
  /* example data array */
  double xp[] = { 0.0, 2.0,-2.0, 0.0, 2.0,-2.0, 0.0 };
  double yp[] = { 0.0, 2.0, 3.0, 6.0, 7.5, 9.0,10.0 };
  double weight[] = { 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0 };
  double knot[] = { 0, 0, 0, 0, 1, 1, 1, 2, 2, 2, 2 };

  /* creation of x-values and y-values vector */
  num = sizeof(xp) / sizeof(double);
  zListInit( &seq );
  for( i=num-1; i>=0; i-- ){
    v = zVecCreateList( 2, xp[i], yp[i] );
    zSeqEnqueue( &seq, v, 1.0 );
  }
  zNURBSCreate( &nurbs, &seq, DIM );
  zSeqFree( &seq );

  for( i=0; i<num ; i++ )
    zNURBSWeight(&nurbs, i) = weight[i];

  zVecClear( nurbs.knot );
  for( i=0; i<zVecSize(nurbs.knot); i++ )
    zVecSetElem( nurbs.knot, i, knot[i] );

  /* creation of spline interpolator */
  v = zVecAlloc( 2 );
  for( i=0; i<=STEP; i++ ){
    t = ( zNURBSKnotE(&nurbs)-zNURBSKnot0(&nurbs) ) * i / STEP + zNURBSKnot0(&nurbs);
    printf( "%1.2f ", t );
    zNURBSVec( &nurbs, t, v );
    printf( "%g %g ", zVecElem(v,0), zVecElem(v,1) );
    zNURBSVecDiff( &nurbs, t, v, 0 );
    printf( "%g %g ", zVecElem(v,0), zVecElem(v,1) );
    zNURBSVecDiff( &nurbs, t, v, 1);
    printf( "%g %g ", zVecElem(v,0), zVecElem(v,1) );
    zNURBSVecDiff( &nurbs, t, v, 2);
    printf( "%g %g\n", zVecElem(v,0), zVecElem(v,1) );
  }
  zVecFree( v );
  zNURBSDestroy( &nurbs );
  return 0;
}
