#include <zm/zm_ip.h>

int main(int argc, char *argv[])
{
  zSeq seq;
  zIPData dat;
  int point_num;
  double t, tmax;
  register int i;
  /* example data array */
  double tp[] = { 0, 1, 3, 6 };

  /* creation of x-values and y-values vector */
  point_num = sizeof(tp) / sizeof(double);
  zSeqInit( &seq );
  for( t=0, i=0; i<point_num; i++ ){
    zSeqEnqueue( &seq, zVecAlloc(1), tp[i]-t );
    t = tp[i];
  }
  /* creation of spline interpolator */
  zIPDataAlloc( &dat, &seq );
  tmax = tp[point_num-1];

  for( i=0; ; i++ ){
    t = tp[0] + 0.1 * i;
    printf( "%g %d\n", t, zIPSeg(&dat,t) );
    if( t > tmax ) break;
  }

  /* destruction of instances */
  zIPDataFree( &dat );
  return 0;  
}
