#include <zm/zm_stat.h>
#include <zm/zm_rand.h>

#define N 1000

int main(int argc, char *argv[])
{
  zHistogram histogram;
  double data[N];
  int i;

  zRandInit();
  for( i=0; i<N; i++ )
    data[i] = zRandND( NULL, 10.0, 5.0 );
  zHistogramCreateAuto( &histogram, data, N, argc > 1 ? atoi(argv[1]) : 20 );
  zHistogramFPrint( stdout, &histogram );
  zHistogramDestroy( &histogram );
  return 0;
}
