#include <zm/zm.h>

#define TEST_VEC_SIZE 10
#define N 1000
#define TOL  (1.0e-10)

void assert_clone(void)
{
  const int size = TEST_VEC_SIZE;
  zVec src, dest;

  src = zVecAlloc( size );
  zVecRandUniform( src, -10, 10 );
  dest = zVecClone( src );
  zAssert( zVecClone, zVecMatch( src, dest ) );
  zVecFree( src );
  zVecFree( dest );
}

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
  zVecFreeAtOnce( 3, test_vec1, test_vec2, test_vec3 );
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
  zVecFreeAtOnce( 2, test_vec1, test_vec2 );
}

void assert_shift(void)
{
  const int size = TEST_VEC_SIZE;
  zVec src, dest, error;
  double shift;
  int i, j;
  bool result1, result2;

  src = zVecAlloc( size );
  dest = zVecAlloc( size );
  error = zVecAlloc( size );
  for( result1=result2=true, i=0; i<N; i++ ){
    shift = zRandF( -10, 10 );
    zVecRandUniform( src, -10, 10 );
    zVecShift( src, shift, dest );
    zVecSub( dest, src, error );
    for( j=0; j<zVecSizeNC(error); j++ )
      if( !zEqual( zVecElemNC(error,j), shift, zTOL ) ) result1 = false;
    zVecShiftDRC( dest, -shift );
    if( !zVecEqual( src, dest, zTOL ) ) result2 = false;
  }
  zVecFreeAtOnce( 3, src, dest, error );
  zAssert( zVecShift, result1 );
  zAssert( zVecShiftDRC, result2 );
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
  zVecFreeAtOnce( 3, test_vec1, test_vec2, test_vec3 );
}

void assert_scale(void)
{
  const int size = TEST_VEC_SIZE;
  zVec min, max, src, dest;
  double elem_min, elem_max;
  int i, j;
  bool result1, result2;

  min = zVecAlloc( size );
  max = zVecAlloc( size );
  src = zVecAlloc( size );
  dest = zVecAlloc( size );
  for( result1=true, i=0; i<N; i++ ){
    zVecRandUniform( min, -10, 10 );
    zVecRandUniform( max, -10, 10 );
    for( j=0; j<zVecSizeNC(min); j++ )
      if( zVecElemNC(min,j) > zVecElemNC(max,j) )
        zSwap( double, zVecElemNC(min,j), zVecElemNC(max,j) );
    zVecRandUniform( src, 0, 1 );
    zVecScale( src, min, max, dest );
    zVecSubDRC( dest, min );
    zVecSubDRC( max, min );
    zVecDemDRC( dest, max );
    if( !zVecEqual( src, dest, TOL ) ){
      zVecSubDRC( src, dest );
      zVecPrint( src );
      result1 = false;
    }
  }
  for( result2=true, i=0; i<N; i++ ){
    elem_min = zRandF( -10,  0 );
    elem_max = zRandF(   0, 10 );
    zVecRandUniform( src, 0, 1 );
    zVecScaleUniform( src, elem_min, elem_max, dest );
    zVecShiftDRC( dest, -elem_min );
    zVecDivDRC( dest, ( elem_max - elem_min ) );
    if( !zVecEqual( src, dest, TOL ) ){
      zVecSubDRC( src, dest );
      zVecPrint( src );
      result2 = false;
    }
  }
  zVecFreeAtOnce( 4, min, max, src, dest );
  zAssert( zVecScale, result1 );
  zAssert( zVecScaleUniform, result2 );
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
  zVecFreeAtOnce( 2, test_vec1, test_vec2 );
}

#define N_NODE 1000
#define DIM 3

void assert_nearest_neighbor(void)
{
  zVecList list;
  zVecTree tree, *node;
  zVec v, nn;
  int i;
  double dmin1, dmin2;
  bool result;

  zListInit( &list );
  zVecTreeInit( &tree, DIM );
  v = zVecAlloc( DIM );
  for( i=0; i<N_NODE; i++ ){
    zVecRandUniform( v, -10, 10 );
    zVecTreeAdd( &tree, v );
    zVecListInsertHead( &list, v );
  }
  zVecRandUniform( v, -10, 10 );

  dmin1 = zVecTreeNN( &tree, v, &node );
  dmin2 = zVecListNN( &list, v, &nn );
  result = ( dmin1 == dmin2 ) && zVecEqual( node->v, nn, 0 ) ? true : false;

  zVecFree( v );
  zVecListDestroy( &list );
  zVecTreeDestroy( &tree );
  zAssert( zVecListNN + zVecTreeNN, result );
}

int main(void)
{
  zRandInit();
  assert_get_put();
  assert_misc();
  assert_shift();
  assert_arith();
  assert_scale();
  assert_normalize();

  assert_nearest_neighbor();

  return EXIT_SUCCESS;
}
