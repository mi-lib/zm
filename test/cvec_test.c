#include <zm/zm.h>

#define TEST_VEC_SIZE 1000

void assert_misc(void)
{
  zCVec v1, v2, v3;

  v1 = zCVecAlloc( TEST_VEC_SIZE );
  zCVecRandUniform( v1, -10, -10, 10, 10 );
  v2 = zCVecClone( v1 );
  v3 = zCVecAlloc( TEST_VEC_SIZE );
  zCVecRandUniform( v3, -10, -10, 10, 10 );
  zAssert( zCVecIsEqual, zCVecIsEqual(v1,v2) && !zCVecIsEqual(v1,v3) );

  zCVecFree( v1 );
  zCVecFree( v2 );
  zCVecFree( v3 );
}

int main(void)
{
  zRandInit();
  assert_misc();

  return EXIT_SUCCESS;
}
