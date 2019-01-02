/* 2014. 8.25 Thanks to Mr. Ken'ya Tanaka */
#include <zm/zm_nurbs.h>

#define STEP 100
#define DIM  3

int main(int argc, char *argv[])
{
  zNURBS1 nurbs;
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
  zNURBS1Create( &nurbs, &seq, DIM );
  zSeqFree( &seq );

  for( i=0; i<num ; i++ )
    zNURBS1Weight(&nurbs, i) = weight[i];

  zVecClear( nurbs.knot );
  for( i=0; i<zVecSize(nurbs.knot); i++ )
    zVecSetElem( nurbs.knot, i, knot[i] );

  /* creation of spline interpolator */
  v = zVecAlloc( 2 );
  for( i=0; i<=STEP; i++ ){
    t = ( zNURBS1KnotE(&nurbs)-zNURBS1Knot0(&nurbs) ) * i / STEP + zNURBS1Knot0(&nurbs);
    printf( "%1.2f ", t );
    zNURBS1Vec( &nurbs, t, v );
    printf( "%g %g ", zVecElem(v,0), zVecElem(v,1) );
    zNURBS1VecDiff( &nurbs, t, v, 0 );
    printf( "%g %g ", zVecElem(v,0), zVecElem(v,1) );
    zNURBS1VecDiff( &nurbs, t, v, 1);
    printf( "%g %g ", zVecElem(v,0), zVecElem(v,1) );
    zNURBS1VecDiff( &nurbs, t, v, 2);
    printf( "%g %g\n", zVecElem(v,0), zVecElem(v,1) );
  }
  zVecFree( v );
  zNURBS1Destroy( &nurbs );
  return 0;
}
