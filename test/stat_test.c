#include <zm/zm.h>

#define NUM 100

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
  return 0;
}
