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

/* division by a tiny real number, and normalization with tolerances */
void assert_div_normalize_tol(void)
{
  zCVec v1, v2;
  double k[] = { 1.0e-13, 4.0e-309 };
  int i, j;
  bool result_div = true, result_zero, result_normalize;

  v1 = zCVecAlloc( TEST_VEC_SIZE );
  v2 = zCVecAlloc( TEST_VEC_SIZE );
  for( j=0; j<sizeof(k)/sizeof(double); j++ ){
    zCVecRandUniform( v1, -1, -1, 1, 1 );
    if( k[j] < 1.0e-300 ) zCVecMulDRC( v1, 1.0e-300 ); /* avoid overflow of quotients */
    if( !zCVecDiv( v1, k[j], v2 ) ){
      result_div = false;
      continue;
    }
    for( i=0; i<TEST_VEC_SIZE; i++ )
      if( fabs( zCVecElemNC(v2,i)->re - zCVecElemNC(v1,i)->re/k[j] ) > DBL_EPSILON * fabs( zCVecElemNC(v1,i)->re/k[j] ) ||
          fabs( zCVecElemNC(v2,i)->im - zCVecElemNC(v1,i)->im/k[j] ) > DBL_EPSILON * fabs( zCVecElemNC(v1,i)->im/k[j] ) )
        result_div = false;
  }
  result_zero = zCVecDiv( v1, 0, v2 ) == NULL;
  zCVecRandUniform( v1, -1, -1, 1, 1 );
  zCVecMulDRC( v1, 1.0e-13 / zCVecNorm( v1 ) );
  result_normalize = !zCVecNormalize( v1, 0, v2 ) && zCVecIsTol( v2, 0 ) &&
    zCVecNormalize( v1, 1.0e-20, v2 ) && zIsTol( zCVecNorm( v2 ) - 1, 10*DBL_EPSILON ) &&
    zCVecNormalizeDRC( v1, 1.0e-20 ) && zIsTol( zCVecNorm( v1 ) - 1, 10*DBL_EPSILON );
  zCVecFree( v1 );
  zCVecFree( v2 );
  zAssert( zCVecDiv (very small divisor), result_div );
  zAssert( zCVecDiv (division by zero), result_zero );
  zAssert( zCVecNormalize + zCVecNormalizeDRC (tolerance), result_normalize );
}

int main(void)
{
  zRandInit();
  assert_misc();
  assert_cmul_cdiv_drc();
  assert_cdiv_scale();
  assert_div_normalize_tol();
  return EXIT_SUCCESS;
}
