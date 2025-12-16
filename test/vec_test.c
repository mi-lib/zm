#include <zm/zm.h>

#define TEST_VEC_SIZE 10
#define N 1000
#define TOL  (1.0e-10)

void assert_vec_clone(void)
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

void assert_vec_get_put(void)
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

void assert_vec_resize(void)
{
  zVec src, dest, part;
  const int size = 5;
  const int regsize = 3;
  const int n = 100;
  int i;
  bool result = true;

  src  = zVecAlloc( size );
  dest = zVecAlloc( size );
  part = zVecAlloc( regsize );
  for( i=0; i<n; i++ ){
    zVecRandUniform( src, -10, 10 );
    zVecResetSize( dest );
    zVecCopy( src, dest );
    zVecGet( src, 0, part );
    zVecResize( dest, regsize );
    if( !zVecEqual( dest, part, zTOL ) ) result = false;
  }
  zVecFreeAtOnce( 3, src, dest, part );
  zAssert( zVecResetSize + zVecResize, result );
}

void assert_vec_misc(void)
{
  const int size = TEST_VEC_SIZE;
  zVec test_vec1, test_vec2;
  double val1, val2, dval;
  int i;
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

void assert_vec_shift(void)
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

void assert_vec_arith(void)
{
  const int size = TEST_VEC_SIZE;
  zVec test_vec1, test_vec2, test_vec3;
  double k;
  int i;
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

void assert_vec_scale(void)
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

void assert_vec_normalize(void)
{
  const int size = TEST_VEC_SIZE;
  zVec test_vec1, test_vec2;
  double norm;
  int i;
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

void assert_vec_orthogonalize(void)
{
  zVec v, n, o;
  const int size = 10, testnum = 100;
  int i;
  bool result1 = true, result2 = true;

  v = zVecAlloc( size );
  n = zVecAlloc( size );
  o = zVecAlloc( size );
  for( i=0; i<testnum; i++ ){
    zVecRandUniform( v, -10, 10 );
    zVecRandUniform( n, -10, 10 );
    zVecOrthogonalize( v, n, o );
    if( !zIsTiny( zVecInnerProd( n, o ) ) ){
      eprintf( "[%d] error = %g\n", i, zVecInnerProd( n, o ) );
      result1 = false;
    }
    zVecOrthogonalize( v, n, n );
    if( !zVecEqual( o, n, zTOL ) ) result2 = false;
  }
  zVecFreeAtOnce( 3, v, n, o );
  zAssert( zVecProj + zVecOrthogonalize, result1 );
  zAssert( zVecOrthogonalize (override), result2 );
}

double vec_edge_dist2(const zVec v, const zVec v1, const zVec v2)
{
  zVec e, n, o;
  double d = -1;

  e = zVecAlloc( zVecSizeNC(v) );
  n = zVecAlloc( zVecSizeNC(v) );
  o = zVecAlloc( zVecSizeNC(v) );
  if( e && n && o ){
    zVecSubNC( v2, v1, n );
    zVecSubNC( v,  v1, e );
    d = zVecNorm( zVecOrthogonalize( e, n, o ) );
  }
  zVecFree( e );
  zVecFree( n );
  zVecFree( o );
  return d;
}

void assert_vec_edge_dist(void)
{
  zVec v, v1, v2;
  double d1, d2;
  const int size = 10, testnum = 10;
  int i;
  bool result = true;

  v  = zVecAlloc( size );
  v1 = zVecAlloc( size );
  v2 = zVecAlloc( size );
  for( i=0; i<testnum; i++ ){
    zVecRandUniform( v,  -10, 10 );
    zVecRandUniform( v1, -10, 10 );
    zVecRandUniform( v2, -10, 10 );
    if( !zEqual( ( d1 = zVecEdgeDist( v, v1, v2 ) ), ( d2 = vec_edge_dist2( v, v1, v2 ) ), zTOL ) ){
      eprintf( "%g  / %g\n", d1, d2 );
      result = false;
    }
  }
  zVecFreeAtOnce( 3, v, v1, v2 );
  zAssert( zVecEdgeDist, result );
}

#define N_NODE 1000
#define DIM 3

void assert_vec_nearest_neighbor(void)
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
  assert_vec_get_put();
  assert_vec_resize();
  assert_vec_misc();
  assert_vec_shift();
  assert_vec_arith();
  assert_vec_scale();
  assert_vec_normalize();
  assert_vec_orthogonalize();
  assert_vec_edge_dist();

  assert_vec_nearest_neighbor();

  return EXIT_SUCCESS;
}
