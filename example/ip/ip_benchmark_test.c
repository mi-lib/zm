#include <zm/zm_ip.h>

#define TEST 1
#define DT   0.1

int enqueue(zSeq *seq, int point_num, double tp[], double vp[])
{
  zVec v;
  double t;
  register int i;
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
  register int i;
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

void output_plotscript(void)
{
  printf( "set terminal png\n" );
  printf( "set output 'ip.png'\n" );
  printf( "unset key\n" );
  printf( "set multiplot\n" );
  printf( "set size 0.5, 0.5\n" );
  printf( "set origin 0, 0\n" );
  printf( "plot [0:6][-4:5] 0, 'lag' u 1:2 w l, 'lin' u 1:2 w l lt 4, 'ip.dat' w p pt 7 lt 8\n" );
  printf( "plot [0:6][-4:5] 0, 'che' u 1:2 w l, 'ip.dat' w p pt 7 lt 8\n" );
  printf( "set origin 0, 0.5\n" );
  printf( "plot [0:6][-4:5] 0, 'sp1' u 1:2 w l, 'ip.dat' w p pt 7 lt 8\n" );
  printf( "set origin 0.5, 0.5\n" );
  printf( "plot [0:6][-4:5] 0, 'sp2' u 1:2 w l, 'ip.dat' w p pt 7 lt 8\n" );
  printf( "set origin 0.5, 0\n" );
  printf( "plot [0:6][-4:5] 0, 'aki' u 1:2 w l, 'ip.dat' w p pt 7 lt 8\n" );
  printf( "unset multiplot\n" );
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
#else
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

  /* PCHIP */
  zIPCreatePCHIP( &ip, &seq );
  output( &ip, point_num, tp, DT, "pch" );
  zIPDestroy( &ip );

  zSeqFree( &seq );

  output_plotscript();
  return 0;  
}
