#include <zm/zm.h>

#define NUM 100

bool assert_is_included(void)
{
  double data[NUM];
  register int i;

  for( i=0; i<NUM; i++ )
    data[i] = zRandF(-10,10);
  for( i=0; i<NUM; i++ )
    if( !zDataIsIncluded( data, NUM, data[i] ) ) return false;
  return true;
}

int main(int argc, char *argv[])
{
  zRandInit();
  zAssert( zDataIsIncluded, assert_is_included() );
  return 0;
}
