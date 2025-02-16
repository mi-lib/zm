#include <zm/zm.h>

#define NUM 100

void assert_histogram(void)
{
  double data[] = { 0, 1, 1, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 5, 5, 6 };
  zHistogram histogram;
  bool result = true;

  zHistogramCreateAuto( &histogram, data, sizeof(data)/sizeof(double), 7 );
  if( zHistogramBinSize( &histogram ) != 7 ||
      zHistogramBinCount( &histogram, 0 ) != 1 ||
      zHistogramBinCount( &histogram, 1 ) != 2 ||
      zHistogramBinCount( &histogram, 2 ) != 3 ||
      zHistogramBinCount( &histogram, 3 ) != 4 ||
      zHistogramBinCount( &histogram, 4 ) != 3 ||
      zHistogramBinCount( &histogram, 5 ) != 2 ||
      zHistogramBinCount( &histogram, 6 ) != 1 ) result = false;
  zHistogramDestroy( &histogram );
  zAssert( zHistogramCreateAuto, result );
}

void assert_is_included(void)
{
  double data[NUM];
  int i;
  bool result = true;

  for( i=0; i<NUM; i++ )
    data[i] = zRandF(-10,10);
  for( i=0; i<NUM; i++ )
    if( !zDataIsIncluded( data, NUM, data[i], zTOL ) ) result = false;
  zAssert( zDataIsIncluded, result );
}

int main(int argc, char *argv[])
{
  zRandInit();
  assert_is_included();
  assert_histogram();
  return 0;
}
