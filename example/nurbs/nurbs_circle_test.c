#include <zm/zm_nurbs.h>

#define STEP 100

int main(int argc, char *argv[])
{
  zNURBS nurbs;
  zSeq seq;
  zVec v;
  int num, i;
  double xp[] = { 1.0, 1.0, 0.0, -1.0, -1.0, -1.0,  0.0,  1.0, 1.0 };
  double yp[] = { 0.0, 1.0, 1.0,  1.0,  0.0, -1.0, -1.0, -1.0, 0.0 };

  /* create xy-values */
  num = sizeof(xp) / sizeof(double);
  zListInit( &seq );
  for( i=0; i<num; i++ ){
    v = zVecCreateList( 2, xp[i], yp[i] );
    zSeqEnqueue( &seq, v, 1.0 /* dummy */ );
  }
  zNURBSCreate( &nurbs, &seq, 2 );
  zSeqFree( &seq );

  zNURBSSetKnot( &nurbs, 0, 0.0 );
  zNURBSSetKnot( &nurbs, 1, 0.0 );
  zNURBSSetKnot( &nurbs, 2, 0.0 );
  zNURBSSetKnot( &nurbs, 3, zPI_2 );
  zNURBSSetKnot( &nurbs, 4, zPI_2 );
  zNURBSSetKnot( &nurbs, 5, zPI );
  zNURBSSetKnot( &nurbs, 6, zPI );
  zNURBSSetKnot( &nurbs, 7, 1.5 * zPI );
  zNURBSSetKnot( &nurbs, 8, 1.5 * zPI );
  zNURBSSetKnot( &nurbs, 9, zPIx2 );
  zNURBSSetKnot( &nurbs,10, zPIx2 );
  zNURBSSetKnot( &nurbs,11, zPIx2 );

  zNURBSSetWeight( &nurbs, 0, 1.0 );
  zNURBSSetWeight( &nurbs, 1, sqrt(2) * 0.5 );
  zNURBSSetWeight( &nurbs, 2, 1.0 );
  zNURBSSetWeight( &nurbs, 3, sqrt(2) * 0.5 );
  zNURBSSetWeight( &nurbs, 4, 1.0 );
  zNURBSSetWeight( &nurbs, 5, sqrt(2) * 0.5 );
  zNURBSSetWeight( &nurbs, 6, 1.0 );
  zNURBSSetWeight( &nurbs, 7, sqrt(2) * 0.5 );
  zNURBSSetWeight( &nurbs, 8, 1.0 );

  /* create a spline interpolator */
  v = zVecAlloc( 2 );
  for( i=0; i<STEP; i++ ){
    zNURBSVec( &nurbs, zPIx2 * i / STEP, v );
    printf( "%g %g ", zVecElem(v,0), zVecElem(v,1) );
    zVecNormalizeDRC( v );
    printf( "%g %g\n", zVecElem(v,0), zVecElem(v,1) );
  }
  zVecFree( v );
  zNURBSDestroy( &nurbs );
  return 0;  
}
