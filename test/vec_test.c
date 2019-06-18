#include <zm/zm.h>

#define TEST_VEC_SIZE 10

void assert_get_put(void)
{
  const int size = TEST_VEC_SIZE;
  zVec test_vec1, test_vec2, test_vec3;

  test_vec1 = zVecAlloc( size );
  test_vec2 = zVecAlloc( size );
  test_vec3 = zVecAlloc( 3 );
  zVecRandUniform( test_vec1, -10, 10 );
  zVecRandUniform( test_vec2, -10, 10 );

  zVecGet( test_vec1, 2, test_vec3 );
  zAssert( zVecGet, memcmp( zVecBufNC(test_vec1)+2, zVecBufNC(test_vec3), sizeof(double)*3 ) == 0 );
  zVecPut( test_vec2, 2, test_vec3 );
  zAssert( zVecPut, memcmp( zVecBufNC(test_vec2)+2, zVecBufNC(test_vec3), sizeof(double)*3 ) == 0 );
  zVecCopy( test_vec1, test_vec2 );
  zVecSwap( test_vec2, 3, 6 );
  zAssert( zVecSwap, zVecElemNC(test_vec1,3) == zVecElemNC(test_vec2,6) && zVecElemNC(test_vec1,6) == zVecElemNC(test_vec2,3) );
  zVecFreeAO( 3, test_vec1, test_vec2, test_vec3 );
}

void assert_misc(void)
{
  const int size = TEST_VEC_SIZE;
  zVec test_vec1, test_vec2;
  double val1, val2, dval;
  register int i;
  bool result;

  test_vec1 = zVecAlloc( size );
  test_vec2 = zVecAlloc( size );
  val1 = zRandF( -10, 10 );
  val2 = zRandF( -10, 10 );
  zVecSetAll( test_vec1, val1 );
  for( result=true, i=0; i<size; i++ )
    if( zVecElemNC(test_vec1,i) != val1 ) result = false;
  zAssert( zVecSetAll, result );
  zVecLinSpace( test_vec1, val1, val2 );
  dval = ( val2 - val1 ) / ( size - 1 );
  for( result=true, i=1; i<size; i++ )
    if( !zIsTiny( zVecElemNC(test_vec1,i) - zVecElemNC(test_vec1,i-1) - dval ) ) result = false;
  zAssert( zVecLinSpace, result );
  zVecCopy( test_vec1, test_vec2 );
  zVecShift( test_vec2, dval );
  for( result=true, i=0; i<size; i++ )
    if( !zIsTiny( zVecElemNC(test_vec2,i) - zVecElemNC(test_vec1,i) - dval ) ) result = false;
  zAssert( zVecShift, result );
  zVecFreeAO( 2, test_vec1, test_vec2 );
}

void assert_arith(void)
{
  const int size = TEST_VEC_SIZE;
  zVec test_vec1, test_vec2, test_vec3;
  double k;
  register int i;
  bool result;

  test_vec1 = zVecAlloc( size );
  test_vec2 = zVecAlloc( size );
  test_vec3 = zVecAlloc( size );
  zVecRandUniform( test_vec1, -10, 10 );
  zVecRandUniform( test_vec2, -10, 10 );
  k = zRandF( -10, 10 );

  zVecAdd( test_vec1, test_vec2, test_vec3 );
  for( result=true, i=0; i<size; i++ )
    if( !zIsTiny( zVecElemNC(test_vec1,i)+zVecElemNC(test_vec2,i)-zVecElemNC(test_vec3,i) ) ) result = false;
  zAssert( zVecAdd, result );
  zVecSub( test_vec1, test_vec2, test_vec3 );
  for( result=true, i=0; i<size; i++ )
    if( !zIsTiny( zVecElemNC(test_vec1,i)-zVecElemNC(test_vec2,i)-zVecElemNC(test_vec3,i) ) ) result = false;
  zAssert( zVecSub, result );
  zVecRev( test_vec1, test_vec3 );
  for( result=true, i=0; i<size; i++ )
    if( !zIsTiny( zVecElemNC(test_vec1,i)+zVecElemNC(test_vec3,i) ) ) result = false;
  zAssert( zVecRev, result );
  zVecMul( test_vec1, k, test_vec3 );
  for( result=true, i=0; i<size; i++ )
    if( !zIsTiny( zVecElemNC(test_vec1,i)*k-zVecElemNC(test_vec3,i) ) ) result = false;
  zAssert( zVecMul, result );
  zVecDiv( test_vec1, k, test_vec3 );
  for( result=true, i=0; i<size; i++ )
    if( !zIsTiny( zVecElemNC(test_vec1,i)/k-zVecElemNC(test_vec3,i) ) ) result = false;
  zAssert( zVecDiv, result );
  zVecAmp( test_vec1, test_vec2, test_vec3 );
  for( result=true, i=0; i<size; i++ )
    if( !zIsTiny( zVecElemNC(test_vec1,i)*zVecElemNC(test_vec2,i)-zVecElemNC(test_vec3,i) ) ) result = false;
  zAssert( zVecAmp, result );
  zVecDem( test_vec1, test_vec2, test_vec3 );
  for( result=true, i=0; i<size; i++ )
    if( !zIsTiny( zVecElemNC(test_vec1,i)/zVecElemNC(test_vec2,i)-zVecElemNC(test_vec3,i) ) ) result = false;
  zAssert( zVecDem, result );
  zVecCat( test_vec1, k, test_vec2, test_vec3 );
  for( result=true, i=0; i<size; i++ )
    if( !zIsTiny( zVecElemNC(test_vec1,i)+k*zVecElemNC(test_vec2,i)-zVecElemNC(test_vec3,i) ) ) result = false;
  zAssert( zVecCat, result );
  zVecFreeAO( 3, test_vec1, test_vec2, test_vec3 );
}

void assert_normalize(void)
{
  const int size = TEST_VEC_SIZE;
  zVec test_vec1, test_vec2;
  double norm;
  register int i;
  bool result;

  test_vec1 = zVecAlloc( size );
  test_vec2 = zVecAlloc( size );
  zVecRandUniform( test_vec1, -10, 10 );
  norm = zVecNorm( test_vec1 );
  zVecNormalize( test_vec1, test_vec2 );
  zAssert( zVecNormalize, zIsTiny( zVecNorm(test_vec2) - 1 ) );
  for( result=true, i=0; i<size; i++ )
    if( !zIsTiny( zVecElemNC(test_vec1,i)/zVecElemNC(test_vec2,i) - norm ) ) result = false;
  zAssert( zVecNorm, result );
  zVecFreeAO( 2, test_vec1, test_vec2 );
}

int main(void)
{
  zRandInit();
  assert_get_put();
  assert_misc();
  assert_arith();
  assert_normalize();

  return EXIT_SUCCESS;
}
