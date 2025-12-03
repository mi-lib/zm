#include <zm/zm.h>

#define TEST_VEC_SIZE    10
#define TEST_VEC_BUFSIZ 100

void assert_raw_vec_get_put(void)
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

void assert_raw_vec_misc(void)
{
  double test_vec1[TEST_VEC_BUFSIZ];
  double test_vec2[TEST_VEC_BUFSIZ];
  double val1, val2, dval;
  int i;
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
  zRawVecShift( test_vec1, TEST_VEC_SIZE, dval, test_vec2 );
  for( result=true, i=0; i<TEST_VEC_SIZE; i++ )
    if( !zIsTiny( test_vec2[i] - test_vec1[i] - dval ) ) result = false;
  zAssert( zRawVecShift, result );
  zRawVecCopy( test_vec1, test_vec2, TEST_VEC_SIZE );
  zRawVecShiftDRC( test_vec2, TEST_VEC_SIZE, dval );
  for( result=true, i=0; i<TEST_VEC_SIZE; i++ )
    if( !zIsTiny( test_vec2[i] - test_vec1[i] - dval ) ) result = false;
  zAssert( zRawVecShiftDRC, result );
}

void assert_raw_vec_equal(void)
{
  double *v1, *v2;
  const int size = 10;
  int i;
  bool result1, result2;

  v1 = zAlloc( double, size );
  v2 = zAlloc( double, size );
  if( !v1 || !v2 ){
    ZALLOCERROR();
    return;
  }
  zRawVecRandUniform( v1, size, -10, 10 );
  zRawVecCopy( v1, v2, size );
  for( i=0; i<size; i++ )
    *( v2 + i ) += zTOL * 0.5;
  result1 = zRawVecEqual( v1, v2, size, zTOL );
  for( i=0; i<size; i++ )
    *( v2 + i ) += zTOL * 10;
  result2 = zRawVecEqual( v1, v2, size, zTOL );
  free( v1 );
  free( v2 );
  zAssert( zRawVecEqual (positive case), result1 );
  zAssert( zRawVecEqual (negative case), !result2 );
}

void assert_raw_vec_match(void)
{
  double *v1, *v2;
  const int size = 10;
  bool result;

  v1 = zAlloc( double, size );
  v2 = zAlloc( double, size );
  if( !v1 || !v2 ){
    ZALLOCERROR();
    return;
  }
  zRawVecRandUniform( v1, size, -10, 10 );
  zRawVecCopy( v1, v2, size );
  result = zRawVecMatch( v1, v2, size );
  free( v1 );
  free( v2 );
  zAssert( zRawVecCopy + zRawVecMatch, result );
}

void assert_raw_vec_arith(void)
{
  double test_vec1[TEST_VEC_BUFSIZ];
  double test_vec2[TEST_VEC_BUFSIZ];
  double test_vec3[TEST_VEC_BUFSIZ];
  double k;
  int i;
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

void assert_raw_vec_innerprod(void)
{
  const double v1[]  = { 3,-2, 5 };
  const double v2[]  = {-6, 1,-4 };
  const double v3[]  = { 5, 3,-6 };
  const double vec[] = { 2,-1, 3 };
  const int size = 3;

  zAssert( zRawVecInnerProd,
    zEqual( zRawVecInnerProd( v1, vec, size ), 23, 0 ) &&
    zEqual( zRawVecInnerProd( v2, vec, size ),-25, 0 ) &&
    zEqual( zRawVecInnerProd( v3, vec, size ),-11, 0 ) );
}

void assert_raw_vec_normalize(void)
{
  double test_vec1[TEST_VEC_BUFSIZ];
  double test_vec2[TEST_VEC_BUFSIZ];
  double norm;
  int i;
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
  assert_raw_vec_get_put();
  assert_raw_vec_equal();
  assert_raw_vec_match();
  assert_raw_vec_misc();
  assert_raw_vec_arith();
  assert_raw_vec_innerprod();
  assert_raw_vec_normalize();
  return EXIT_SUCCESS;
}
