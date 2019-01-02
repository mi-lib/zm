#include <zm/zm_nurbs.h>

#define STEP 100

int main(int argc, char *argv[])
{
  zNURBS1 nurbs;
  zSeq seq;
  zVec v;
  int num, i;
  /* example data array */
  double xp[] = { 1.0, 1.0, 0.0, -1.0, -1.0, -1.0,  0.0,  1.0, 1.0 };
  double yp[] = { 0.0, 1.0, 1.0,  1.0,  0.0, -1.0, -1.0, -1.0, 0.0 };

  /* creation of x-values and y-values vector */
  num = sizeof(xp) / sizeof(double);
  zListInit( &seq );
  for( i=0; i<num; i++ ){
    v = zVecCreateList( 2, xp[i], yp[i] );
    zSeqEnqueue( &seq, v, 1.0 /* dummy */ );
  }
  zNURBS1Create( &nurbs, &seq, 2 );
  zSeqFree( &seq );

  zNURBS1Knot(&nurbs,0) = 0.0;
  zNURBS1Knot(&nurbs,1) = 0.0;
  zNURBS1Knot(&nurbs,2) = 0.0;
  zNURBS1Knot(&nurbs,3) = zPI_2;
  zNURBS1Knot(&nurbs,4) = zPI_2;
  zNURBS1Knot(&nurbs,5) = zPI;
  zNURBS1Knot(&nurbs,6) = zPI;
  zNURBS1Knot(&nurbs,7) = 1.5 * zPI;
  zNURBS1Knot(&nurbs,8) = 1.5 * zPI;
  zNURBS1Knot(&nurbs,9) = zPIx2;
  zNURBS1Knot(&nurbs,10) = zPIx2;
  zNURBS1Knot(&nurbs,11) = zPIx2;

  zNURBS1Weight(&nurbs,0) = 1.0;
  zNURBS1Weight(&nurbs,1) = sqrt(2) * 0.5;
  zNURBS1Weight(&nurbs,2) = 1.0;
  zNURBS1Weight(&nurbs,3) = sqrt(2) * 0.5;
  zNURBS1Weight(&nurbs,4) = 1.0;
  zNURBS1Weight(&nurbs,5) = sqrt(2) * 0.5;
  zNURBS1Weight(&nurbs,6) = 1.0;
  zNURBS1Weight(&nurbs,7) = sqrt(2) * 0.5;
  zNURBS1Weight(&nurbs,8) = 1.0;

  /* creation of spline interpolator */
  v = zVecAlloc( 2 );
  for( i=0; i<STEP; i++ ){
    zNURBS1Vec( &nurbs, zPIx2 * i / STEP, v );
    printf( "%g %g ", zVecElem(v,0), zVecElem(v,1) );
    zVecNormalizeDRC( v );
    printf( "%g %g\n", zVecElem(v,0), zVecElem(v,1) );
  }
  zVecFree( v );
  zNURBS1Destroy( &nurbs );
  return 0;  
}
