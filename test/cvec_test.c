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

/* division by a complex number with a very small or large absolute value */
void assert_cdiv_scale(void)
{
  zCVec v1, v2;
  zComplex z, r;
  double scale[] = { 1.0e-7, 1.0e-200, 1.0e200 };
  int i, j;
  bool result = true, result_zero;

  v1 = zCVecAlloc( TEST_VEC_SIZE );
  v2 = zCVecAlloc( TEST_VEC_SIZE );
  zCVecRandUniform( v1, -10, -10, 10, 10 );
  for( j=0; j<sizeof(scale)/sizeof(double); j++ ){
    zComplexCreate( &z, zRandF(-10,10)*scale[j], zRandF(-10,10)*scale[j] );
    if( !zCVecCDiv( v1, &z, v2 ) ){
      result = false;
      continue;
    }
    for( i=0; i<TEST_VEC_SIZE; i++ ){ /* v2 z = v1 */
      zComplexCMul( zCVecElemNC(v2,i), &z, &r );
      zComplexSubDRC( &r, zCVecElemNC(v1,i) );
      if( zComplexAbs( &r ) > 10 * DBL_EPSILON * zComplexAbs( zCVecElemNC(v1,i) ) ) result = false;
    }
  }
  zComplexCreate( &z, 0, 0 );
  result_zero = zCVecCDiv( v1, &z, v2 ) == NULL;
  zCVecFree( v1 );
  zCVecFree( v2 );
  zAssert( zCVecCDiv (very small or large divisor), result );
  zAssert( zCVecCDiv (division by zero), result_zero );
}

int main(void)
{
  zRandInit();
  assert_misc();
  assert_cmul_cdiv_drc();
  assert_cdiv_scale();
  return EXIT_SUCCESS;
}
