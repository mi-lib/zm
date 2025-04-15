#include <zm/zm_nurbs.h>

#define TEST_NURBS_FILENAME "test_nurbs.ztk"

void create_test_nurbs(zNURBS *nurbs)
{
  zSeq seq;
  zVec v;
  int n, i;
  double xp[] = { 2.0, 3.0, 5.0, 4.0, 5.0, 7.0 };
  double yp[] = { 3.0,-1.0,-2.0, 0.0, 4.0, 1.5 };
  FILE *fp;

  n = sizeof(xp) / sizeof(double);
  zListInit( &seq );
  for( i=0; i<n; i++ ){
    v = zVecCreateList( 2, xp[i], yp[i] );
    zSeqEnqueue( &seq, v, 1.0 );
  }
  zNURBSCreate( nurbs, &seq, 3 );
  zNURBSSetWeight( nurbs, 4, 2.0 );
  zSeqFree( &seq );
  fp = fopen( TEST_NURBS_FILENAME, "w" );
  zNURBSFPrintZTK( fp, nurbs );
  fclose( fp );
}

void read_test_nurbs(zNURBS *nurbs)
{
  ZTK ztk;

  ZTKInit( &ztk );
  ZTKParse( &ztk, TEST_NURBS_FILENAME );
  zNURBSFromZTK( nurbs, &ztk );
  ZTKDestroy( &ztk );
}

void output_nurbs(zNURBS *nurbs)
{
  double t;
  zVec v;
  const int slice = 100;
  int i;
  FILE *fp;

  fp = fopen( "d", "w" );
  v = zVecAlloc( zVecSizeNC( zNURBSCP(nurbs,0) ) );
  zNURBSSetSlice( nurbs, slice );
  for( i=0; i<=zNURBSSlice(nurbs); i++ ){
    t = zNURBSKnotSlice( nurbs, i );
    if( zNURBSVec( nurbs, t, v ) )
      zVecValueFPrint( fp, v );
  }
  zVecFree( v );
  fclose( fp );
}

int main(int argc, char *argv[])
{
  zNURBS nurbs;

  if( argc > 1 ){
    read_test_nurbs( &nurbs );
  } else{
    create_test_nurbs( &nurbs );
  }
  output_nurbs( &nurbs );
  zNURBSDestroy( &nurbs );
  return 0;  
}
