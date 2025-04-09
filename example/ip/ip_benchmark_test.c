#include <zm/zm_ip.h>

#define TEST 0
#define DT   0.01

int enqueue(zSeq *seq, int point_num, double tp[], double vp[])
{
  zVec v;
  double t;
  int i;
  FILE *fp;

  fp = fopen( "ip.dat", "w" );
  zSeqInit( seq );
  for( t=0, i=0; i<point_num; i++ ){
    fprintf( fp, "%.10g %.10g\n", tp[i], vp[i] );
    v = zVecCreateList( 1, vp[i] );
    zSeqEnqueue( seq, v, tp[i]-t );
    t = tp[i];
  }
  fclose( fp );
  return point_num;
}

void output(zIP *ip, int point_num, double tp[], double dt, const char *filename)
{
  double t, tmax;
  zVec v;
  int i;
  FILE *fp;

  fp = fopen( filename, "w" );
  tmax = tp[point_num-1];
  v = zVecAlloc( 1 );
  for( i=0; ; i++ ){
    t = tp[0] + dt * i;
    zIPVec( ip, t, v ); fprintf( fp, "%g %g ", t, zVecElem(v,0) ); /* value */
    zIPVel( ip, t, v ); fprintf( fp, "%g ",  zVecElem(v,0) ); /* velocity */
    zIPAcc( ip, t, v ); fprintf( fp, "%g\n", zVecElem(v,0) ); /* acceleration */
    if( t >= tmax ) break;
  }
  zVecFree( v );
  fclose( fp );
}

void output_plotscript_one(char *type, double tmin, double tmax)
{
  printf( "set output 'ip_%s_p.png'\n", type );
  printf( "plot [%g:%g] 0, '%s' u 1:2 w l, 'ip.dat' w p pt 7 lt 8\n", tmin, tmax, type );
  printf( "set output 'ip_%s_v.png'\n", type );
  printf( "plot [%g:%g] 0, '%s' u 1:3 w l\n", tmin, tmax, type );
  printf( "set output 'ip_%s_a.png'\n", type );
  printf( "plot [%g:%g] 0, '%s' u 1:4 w l\n", tmin, tmax, type );
}

void output_plotscript(double tmin, double tmax)
{
  printf( "unset key\n" );
  printf( "set terminal png\n" );
  output_plotscript_one( "lag", tmin, tmax );
  output_plotscript_one( "sp2", tmin, tmax );
  output_plotscript_one( "aki", tmin, tmax );
  output_plotscript_one( "mak", tmin, tmax );
  output_plotscript_one( "pch", tmin, tmax );
  printf( "!montage -tile 5x3 -geometry 60%% ip_lag_p.png ip_sp2_p.png ip_aki_p.png ip_mak_p.png ip_pch_p.png ip_lag_v.png ip_sp2_v.png ip_aki_v.png ip_mak_v.png ip_pch_v.png ip_lag_a.png ip_sp2_a.png ip_aki_a.png ip_mak_a.png ip_pch_a.png ip.png\n" );
}

int main(int argc, char *argv[])
{
  zSeq seq;
  zIP ip;
  zVec v;
  int point_num;
  /* example data array */
#if TEST == 1
  double tp[] = { 0, 2, 3, 5, 6 };
  double vp[] = { 0, 3, 1,-2, 4 };
#elif TEST == 2
  double tp[] = { 0, 2, 3, 5, 6 };
  double vp[] = { 1, 2,-3, 4, 2 };
#elif TEST == 3
  /* Cauchy-Lorentz function */
  double tp[] = { -1.0, -0.75, -0.5, -0.25, 0, 0.25, 0.5, 0.75, 1.0 };
  double vp[] = { 1.0/26.0, 16.0/241.0, 4.0/29.0, 16.0/41.0, 1.0, 16.0/41.0, 4.0/29.0, 16.0/241.0, 1.0/26.0 };
#elif TEST == 4
  /* modified Akima testcase */
  double tp[] = { 1, 2, 3, 4, 5, 6, 7, 8 };
  double vp[] = {-1,-1,-1, 0, 1, 1, 1, 1 };
#else
  /* PCHIP testcase */
  double tp[]  = { 0, 2, 3, 5, 6, 8, 9, 11, 12, 14, 15 };
  double vp[] = { 10, 10, 10, 10, 10, 10, 10.5, 15, 50, 60, 85 };
#endif

  point_num = enqueue( &seq, sizeof(tp) / sizeof(double), tp, vp );

  /* linear interpolator */
  zIPCreateLinear( &ip, &seq );
  output( &ip, point_num, tp, DT, "lin" );
  zIPDestroy( &ip );

  /* Lagrange interpolator */
  zIPCreateLagrange( &ip, &seq );
  output( &ip, point_num, tp, DT, "lag" );
  zIPDestroy( &ip );

  /* Chebyshev interpolator */
  zIPCreateChebyshev( &ip, &seq );
  output( &ip, point_num, tp, DT, "che" );
  zIPDestroy( &ip );

  /* cubic spline interpolator */
  v = zVecAlloc( 1 );
  zIPCreateSpline( &ip, &seq, ZSPLINE_FIX_EDGE, v, ZSPLINE_FIX_EDGE, v );
  zVecFree( v );
  output( &ip, point_num, tp, DT, "sp1" );
  zIPDestroy( &ip );

  zIPCreateSpline( &ip, &seq, ZSPLINE_FREE_EDGE, NULL, ZSPLINE_FREE_EDGE, NULL );
  output( &ip, point_num, tp, DT, "sp2" );
  zIPDestroy( &ip );

  /* Akima spline interpolator */
  zIPCreateAkima( &ip, &seq );
  output( &ip, point_num, tp, DT, "aki" );
  zIPDestroy( &ip );

  /* Akima spline interpolator */
  zIPCreateModifiedAkima( &ip, &seq );
  output( &ip, point_num, tp, DT, "mak" );
  zIPDestroy( &ip );

  /* PCHIP */
  zIPCreatePCHIP( &ip, &seq );
  output( &ip, point_num, tp, DT, "pch" );
  zIPDestroy( &ip );

  zSeqFree( &seq );

  output_plotscript( tp[0], tp[point_num-1] );
  return 0;  
}
