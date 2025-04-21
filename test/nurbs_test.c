#include <zm/zm_nurbs.h>

void assert_nurbs_clone(void)
{
  zNURBS src, dest;
  zSeq seq;
  zVec vs, vd;
  int num, i, k;
  const int testnum = 100;
  double xp[] = { 2.0, 3.0, 5.0, 4.0, 5.0, 7.0 };
  double yp[] = { 3.0,-1.0,-2.0, 0.0, 4.0, 1.5 };
  bool result = true;

  num = sizeof(xp) / sizeof(double);
  zListInit( &seq );
  for( i=0; i<num; i++ ){
    vs = zVecCreateList( 2, xp[i], yp[i] );
    zSeqEnqueue( &seq, vs, 1.0 );
  }
  zNURBSCreate( &src, &seq, 3 );
  zSeqFree( &seq );
  zNURBSClone( &src, &dest );

  vs = zVecAlloc( 2 );
  vd = zVecAlloc( 2 );
  for( k=0; k<testnum; k++ ){
    i = zRandI( 0, zNURBSSlice(&src) );
    zNURBSVec( &src,  zNURBSKnotSlice(&src ,i), vs );
    zNURBSVec( &dest, zNURBSKnotSlice(&dest,i), vd );
    if( !zVecMatch( vs, vd ) ) result = false;
  }
  zVecFree( vs );
  zVecFree( vd );
  zNURBSDestroy( &dest );
  zNURBSDestroy( &src );
  zAssert( zNURBSClone, result );
}

int main(int argc, char *argv[])
{
  zRandInit();
  assert_nurbs_clone();
  return 0;  
}
