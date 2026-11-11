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
  zAssert( zCVecEqual, zCVecEqual(v1,v2,zTOL) && !zCVecEqual(v1,v3,zTOL) );

  zCVecFree( v1 );
  zCVecFree( v2 );
  zCVecFree( v3 );
}

void assert_cmul_cdiv_drc(void)
{
  zCVec v1, v2, v3;
  zComplex z;

  v1 = zCVecAlloc( TEST_VEC_SIZE );
  v2 = zCVecAlloc( TEST_VEC_SIZE );
  zCVecRandUniform( v1, -10, -10, 10, 10 );
  zComplexCreate( &z, zRandF(1,10), zRandF(-10,10) );
  v3 = zCVecClone( v1 );
  zCVecCMul( v1, &z, v2 );
  zCVecCMulDRC( v3, &z );
  zAssert( zCVecCMulDRC, zCVecEqual( v2, v3, zTOL ) );
  zCVecCopy( v1, v3 );
  zCVecCDiv( v1, &z, v2 );
  zCVecCDivDRC( v3, &z );
  zAssert( zCVecCDivDRC, zCVecEqual( v2, v3, zTOL ) );
  zCVecFree( v1 );
  zCVecFree( v2 );
  zCVecFree( v3 );
}

int main(void)
{
  zRandInit();
  assert_misc();
  assert_cmul_cdiv_drc();
  return EXIT_SUCCESS;
}
