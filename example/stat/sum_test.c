#include <zm/zm_stat.h>

double simple_sum(double *data, int num)
{
  register int i;
  double result;

  for( result=0, i=0; i<num; i++ )
    result += data[i];
  return result;
}

int cmp1(const void *v1, const void *v2)
{
  double d1, d2;

  d1 = *(double *)v1;
  d2 = *(double *)v2;
  return d1 >= d2 ? 1 : -1;
}

int cmp2(const void *v1, const void *v2)
{
  double d1, d2;

  d1 = *(double *)v1;
  d2 = *(double *)v2;
  return d1 <= d2 ? 1 : -1;
}

#define N 100000
#define N2     1

int main(void)
{
  double dat[N], s1, s2, s3, s4;
  int i;

  zRandInit();
  for( i=0; i<N2; i++ )
    dat[i] = zRandF(0,1e16);
  for( i=N2; i<N; i++ )
    dat[i] = zRandF(0,1e-5);

  s2 = zDataSum( dat, N );
  s3 = simple_sum( dat, N ); /* stupid sum */

  qsort( dat, N, sizeof(double), cmp1 );
  s1 = simple_sum( dat, N ); /* probably the most accurate sum */
  qsort( dat, N, sizeof(double), cmp2 );
  s4 = simple_sum( dat, N ); /* probably the worst sum */

  printf( "%.16f\n",       s1 );
  printf( "%.16f %.16f\n", s2, s1-s2 );
  printf( "%.16f %.16f\n", s3, s1-s3 );
  printf( "%.16f %.16f\n", s4, s1-s4 );

  return 0;
}
