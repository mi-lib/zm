#include <zm/zm_ip.h>

void usage(void)
{
  fprintf( stderr, "usage: spline_test <etype1> <etype2>\n" );
  fprintf( stderr, " etype1/2 ... 'fix' or 'free'\n" );
  exit( 1 );
}

/* one can choose each edge type by command line argument.
 */
int etype(char *tp)
{
  if( !strcmp( tp, "fix" ) )
    return ZSPLINE_FIX_EDGE;
  if( !strcmp( tp, "free" ) )
    return ZSPLINE_FREE_EDGE;
  return ZSPLINE_INVALID;
}

int main(int argc, char *argv[])
{
  zSeq seq;
  zIP ip;
  zVec v;
  int point_num;
  double t, tmax;
  register int i;
  /* example data array */
  double tp[] = { 0, 2, 3, 5, 6 };
  double vp1[] = { 0, 3, 1,-2, 4 };
  double vp2[] = { 1, 2,-3, 4, 2 };
  int e1, e2;

  if( argc < 3 ) usage();
  e1 = etype( argv[1] );
  e2 = etype( argv[2] );

  /* creation of x-values and y-values vector */
  point_num = sizeof(tp) / sizeof(double);
  zSeqInit( &seq );
  for( t=0, i=0; i<point_num; i++ ){
    v = zVecCreateList( 2, vp1[i], vp2[i] );
    zSeqEnqueue( &seq, v, tp[i]-t );
    t = tp[i];
  }
  /* creation of spline interpolator */
  v = zVecAlloc( 2 );
  zIPCreateSpline( &ip, &seq, e1, v, e2, v );

  /* value, velocity, acceleration */
  tmax = tp[point_num-1];
  v = zVecAlloc( 2 );
  for( i=0; ; i++ ){
    t = tp[0] + 0.1 * i;
    zIPVec( &ip, t, v );
    printf( "%f %f %f ", t, zVecElem(v,0), zVecElem(v,1) );
    zIPVel( &ip, t, v );
    printf( "%f %f ", zVecElem(v,0), zVecElem(v,1) );
    zIPAcc( &ip, t, v );
    printf( "%f %f\n", zVecElem(v,0), zVecElem(v,1) );
    if( t > tmax ) break;
  }

  /* destruction of instances */
  zIPDestroy( &ip );
  zVecFree( v );
  return 0;  
}
