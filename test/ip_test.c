#include <zm/zm_ip.h>

void assert_spline_cubicequation(void)
{
  zSeq seq;
  zIP ip;
  zVec v, a, b, c, d, v_cubic;
  const int test_num = 100;
  int point_num, i, k;
  double tp[] = { 0, 2, 3, 5, 6 };
  double vp[] = { 0, 3, 1,-2, 4 };
  double t;
  bool result = true;

  point_num = sizeof(tp) / sizeof(double);
  zSeqInit( &seq );
  for( t=0, i=0; i<point_num; i++ ){
    v = zVecCreateList( 1, vp[i] );
    zSeqEnqueue( &seq, v, tp[i]-t );
    t = tp[i];
  }
  /* cubic spline interpolator */
  v = zVecAlloc( 1 );
  zVecZero( v );
  zIPCreateSpline( &ip, &seq, ZSPLINE_FIX_EDGE, v, ZSPLINE_FIX_EDGE, v );
  a = zVecAlloc( 1 );
  b = zVecAlloc( 1 );
  c = zVecAlloc( 1 );
  d = zVecAlloc( 1 );
  v_cubic = zVecAlloc( 1 );
  for( k=0; k<test_num; k++ ){
    i = zRandI( 0, point_num-2 );
    zIPSplineCoeff( &ip, i, a, b, c, d );
    t = zRandF( tp[i], tp[i+1] );
    zIPVec( &ip, t, v );
    t -= tp[i];
    zVecZero( v_cubic );
    zVecCatNCDRC( v_cubic, t*t*t, a );
    zVecCatNCDRC( v_cubic, t*t,   b );
    zVecCatNCDRC( v_cubic, t,     c );
    zVecAddNCDRC( v_cubic,        d );
    if( !zVecEqual( v, v_cubic, zTOL ) ){
      eprintf( "k=%d, i=%d\n", k, i );
      zVecFPrint( stderr, v );
      zVecFPrint( stderr, v_cubic );
      result = false;
    }
  }
  zVecFreeAtOnce( 6, v, a, b, c, d, v_cubic );
  zIPDestroy( &ip );
  zSeqFree( &seq );
  zAssert( zIPSplineCoeff, result );
}

int main(void)
{
  zRandInit();
  assert_spline_cubicequation();
  return 0;
}
