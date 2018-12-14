#include <zm/zm_vec.h>

int main(void)
{
  zVec v;

  zRandInit();
  v = zVecAlloc( 10 );
  zVecRandUniform( v, -10, 10 );
  printf( "(before shifted)\n" ); zVecWrite( v );
  zVecShift( v, -zVecMin(v,NULL) );
  printf( "(after shifted)\n" );  zVecWrite( v );
  zVecFree( v );
  return 0;
}
