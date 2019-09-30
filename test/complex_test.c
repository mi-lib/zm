#include <zm/zm.h>

#define NUM 100

bool assert_isincluded(void)
{
  register int i;
  zComplex c[NUM];

  for( i=0; i<NUM; i++ )
    zComplexCreate( &c[i], zRandF(-10,10), zRandF(-10,10) );
  for( i=0; i<NUM; i++ )
    if( !zComplexValIsIncluded( c, NUM, &c[i], zTOL ) ) return false;
  return true;
}

int main(int argc, char *argv[])
{
  zRandInit();
  zAssert( zComplexValIsIncluded, assert_isincluded() );
  return 0;
}
