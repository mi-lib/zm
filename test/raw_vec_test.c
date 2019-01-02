#include <zm/zm.h>

#define TEST_VEC_SIZE    10
#define TEST_VEC_BUFSIZ 100

void assert_get_put(void)
{
  double test_vec1[TEST_VEC_BUFSIZ];
  double test_vec2[TEST_VEC_BUFSIZ];
  double test_vec3[TEST_VEC_BUFSIZ];

  zRawVecRandUniform( test_vec1, TEST_VEC_SIZE, -10, 10 );
  zRawVecRandUniform( test_vec2, TEST_VEC_SIZE, -10, 10 );
  zRawVecGet( test_vec1, 2, test_vec3, 3 );
  zAssert( zRawVecGet, memcmp( test_vec1+2, test_vec3, sizeof(double)*3 ) == 0 );
  zRawVecPut( test_vec2, 2, test_vec3, 3 );
  zAssert( zRawVecPut, memcmp( test_vec2+2, test_vec3, sizeof(double)*3 ) == 0 );
  zRawVecCopy( test_vec1, test_vec2, TEST_VEC_SIZE );
  zRawVecSwap( test_vec2, 3, 6 );
  zAssert( zRawVecSwap, test_vec1[3] == test_vec2[6] && test_vec1[6] == test_vec2[3] );
}

void assert_misc(void)
{
  double test_vec1[TEST_VEC_BUFSIZ];
  double test_vec2[TEST_VEC_BUFSIZ];
  double val1, val2, dval;
  register int i;
  bool result;

  val1 = zRandF( -10, 10 );
  val2 = zRandF( -10, 10 );
  zRawVecSetAll( test_vec1, TEST_VEC_SIZE, val1 );
  for( result=true, i=0; i<TEST_VEC_SIZE; i++ )
    if( test_vec1[i] != val1 ) result = false;
  zAssert( zRawVecSetAll, result );
  zRawVecLinSpace( test_vec1, TEST_VEC_SIZE, val1, val2 );
  dval = ( val2 - val1 ) / ( TEST_VEC_SIZE - 1 );
  for( result=true, i=1; i<TEST_VEC_SIZE; i++ )
    if( !zIsTiny( test_vec1[i] - test_vec1[i-1] - dval ) ) result = false;
  zAssert( zRawVecLinSpace, result );
  zRawVecCopy( test_vec1, test_vec2, TEST_VEC_SIZE );
  zRawVecShift( test_vec2, TEST_VEC_SIZE, dval );
  for( result=true, i=0; i<TEST_VEC_SIZE; i++ )
    if( !zIsTiny( test_vec2[i] - test_vec1[i] - dval ) ) result = false;
  zAssert( zRawVecShift, result );
}

void assert_arith(void)
{
  double test_vec1[TEST_VEC_BUFSIZ];
  double test_vec2[TEST_VEC_BUFSIZ];
  double test_vec3[TEST_VEC_BUFSIZ];
  double k;
  register int i;
  bool result;

  zRawVecRandUniform( test_vec1, TEST_VEC_SIZE, -10, 10 );
  zRawVecRandUniform( test_vec2, TEST_VEC_SIZE, -10, 10 );
  k = zRandF( -10, 10 );

  zRawVecAdd( test_vec1, test_vec2, test_vec3, TEST_VEC_SIZE );
  for( result=true, i=0; i<TEST_VEC_SIZE; i++ )
    if( !zIsTiny( test_vec1[i]+test_vec2[i]-test_vec3[i] ) ) result = false;
  zAssert( zRawVecAdd, result );
  zRawVecSub( test_vec1, test_vec2, test_vec3, TEST_VEC_SIZE );
  for( result=true, i=0; i<TEST_VEC_SIZE; i++ )
    if( !zIsTiny( test_vec1[i]-test_vec2[i]-test_vec3[i] ) ) result = false;
  zAssert( zRawVecSub, result );
  zRawVecRev( test_vec1, test_vec3, TEST_VEC_SIZE );
  for( result=true, i=0; i<TEST_VEC_SIZE; i++ )
    if( !zIsTiny( test_vec1[i]+test_vec3[i] ) ) result = false;
  zAssert( zRawVecRev, result );
  zRawVecMul( test_vec1, k, test_vec3, TEST_VEC_SIZE );
  for( result=true, i=0; i<TEST_VEC_SIZE; i++ )
    if( !zIsTiny( test_vec1[i]*k-test_vec3[i] ) ) result = false;
  zAssert( zRawVecMul, result );
  zRawVecDiv( test_vec1, k, test_vec3, TEST_VEC_SIZE );
  for( result=true, i=0; i<TEST_VEC_SIZE; i++ )
    if( !zIsTiny( test_vec1[i]/k-test_vec3[i] ) ) result = false;
  zAssert( zRawVecDiv, result );
  zRawVecAmp( test_vec1, test_vec2, test_vec3, TEST_VEC_SIZE );
  for( result=true, i=0; i<TEST_VEC_SIZE; i++ )
    if( !zIsTiny( test_vec1[i]*test_vec2[i]-test_vec3[i] ) ) result = false;
  zAssert( zRawVecAmp, result );
  zRawVecDem( test_vec1, test_vec2, test_vec3, TEST_VEC_SIZE );
  for( result=true, i=0; i<TEST_VEC_SIZE; i++ )
    if( !zIsTiny( test_vec1[i]/test_vec2[i]-test_vec3[i] ) ) result = false;
  zAssert( zRawVecDem, result );
  zRawVecCat( test_vec1, k, test_vec2, test_vec3, TEST_VEC_SIZE );
  for( result=true, i=0; i<TEST_VEC_SIZE; i++ )
    if( !zIsTiny( test_vec1[i]+k*test_vec2[i]-test_vec3[i] ) ) result = false;
  zAssert( zRawVecCat, result );
}

void assert_normalize(void)
{
  double test_vec1[TEST_VEC_BUFSIZ];
  double test_vec2[TEST_VEC_BUFSIZ];
  double norm;
  register int i;
  bool result;

  zRawVecRandUniform( test_vec1, TEST_VEC_SIZE, -10, 10 );
  norm = zRawVecNorm( test_vec1, TEST_VEC_SIZE );
  zRawVecNormalize( test_vec1, TEST_VEC_SIZE, test_vec2 );
  zAssert( zRawVecNormalize, zIsTiny( zRawVecNorm(test_vec2,TEST_VEC_SIZE) - 1 ) );
  for( result=true, i=0; i<TEST_VEC_SIZE; i++ )
    if( !zIsTiny( test_vec1[i]/test_vec2[i] - norm ) ) result = false;
  zAssert( zRawVecNorm, result );
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
