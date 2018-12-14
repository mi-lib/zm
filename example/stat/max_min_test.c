#include <zm/zm_stat.h>

#define N 20
int main(void)
{
  double data[N];
  register int i;

  zRandInit();
  for( i=0; i<N; i++ ){
    data[i] = zRandI( -100, 100 );
    printf( "%d\t", (int)data[i] );
    if( i % 5 == 4 ) zEndl();
  }
  printf( "max = %f\n", zDataMax( data, N, NULL ) );
  printf( "min = %f\n", zDataMin( data, N, NULL ) );
  printf( "|max| = %f\n", zDataAbsMax( data, N, NULL ) );
  printf( "|min| = %f\n", zDataAbsMin( data, N, NULL ) );
  return 0;
}
