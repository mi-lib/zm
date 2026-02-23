#include <zm/zm_vec.h>

void assert_vec_setlist()
{
  zVec v = new zVecStruct {3};
  v->setList( 1.0, 2.0, 3.0 );
  zAssert( C++::zVec (setList), (*v)[0] == 1.0 && (*v)[1] == 2.0 && (*v)[2] == 3.0 );
  delete v;
}

void assert_vec_clone()
{
  const int size = 10;
  zVec src  = new zVecStruct{size};
  zVec dest = new zVecStruct{size};
  src->randUniform( -10, 10 );
  dest = src->clone();
  zAssert( C++::zVec (clone), zVecMatch( src, dest ) );
  delete src;
  delete dest;
}

void assert_vec_get_put()
{
  const int size = 10;
  zVec test_vec1 = new zVecStruct{ size };
  zVec test_vec2 = new zVecStruct{ size };
  zVec test_vec3 = new zVecStruct{ 3 };
  test_vec1->randUniform( -10, 10 );
  test_vec2->randUniform( -10, 10 );
  test_vec1->get( 2, test_vec3 );
  zAssert( C++::zVec (get), memcmp( test_vec1->buf+2, test_vec3->buf, sizeof(double)*3 ) == 0 );
  test_vec2->put( 2, test_vec3 );
  zAssert( C++::zVec (put), memcmp( test_vec2->buf+2, test_vec3->buf, sizeof(double)*3 ) == 0 );
  test_vec2->copy( test_vec1 );
  test_vec2->swap( 3, 6 );
  zAssert( C++::zVec (swap), test_vec1->buf[3] == test_vec2->buf[6] && test_vec1->buf[6] == test_vec2->buf[3] );
  delete test_vec1;
  delete test_vec2;
  delete test_vec3;
}

void assert_vec_resize()
{
  const int size = 5;
  const int regsize = 3;
  const int n = 100;
  zVec src  = new zVecStruct{ size };
  zVec dest = new zVecStruct{ size };
  zVec part = new zVecStruct{ regsize };
  bool result = true;
  for(int i=0; i<n; i++ ){
    src->randUniform( -10, 10 );
    dest->resetSize();
    dest->copy( src );
    src->get( 0, part );
    dest->resize( regsize );
    if( !zVecEqual( dest, part, zTOL ) ) result = false;
  }
  delete src;
  delete dest;
  delete part;
  zAssert( C++::zVec (resetSize + resize), result );
}

void assert_vec_misc()
{
  const int size = 10;
  zVec test_vec1 = new zVecStruct{ size };
  zVec test_vec2 = new zVecStruct{ size };
  double val1 = zRandF( -10, 10 );
  double val2 = zRandF( -10, 10 );
  test_vec1->setAll( val1 );
  bool result = true;
  for(int i=0; i<size; i++ )
    if( test_vec1->buf[i] != val1 ) result = false;
  zAssert( C++::zVec (setAll), result );
  test_vec1->linspace( val1, val2 );
  double dval = ( val2 - val1 ) / ( size - 1 );
  result = true;
  for(int i=1; i<size; i++ )
    if( !zIsTiny( test_vec1->buf[i] - test_vec1->buf[i-1] - dval ) ) result = false;
  zAssert( C++::zVec (linspace), result );
  delete test_vec1;
  delete test_vec2;
}

void assert_vec_shift()
{
  const int size = 10;
  const int n = 100;
  zVec src   = new zVecStruct{ size };
  zVec dest  = new zVecStruct{ size };
  zVec error = new zVecStruct{ size };
  double shift;
  bool result = true;
  for(int i=0; i<n; i++ ){
    src->randUniform( -10, 10 );
    dest->copy( src );
    dest->shift( ( shift = zRandF( -10, 10 ) ) );
    zVecSub( dest, src, error );
    for(int j=0; j<error->size; j++ )
      if( !zEqual( error->buf[j], shift, zTOL ) ) result = false;
  }
  zVecFreeAtOnce( 3, src, dest, error );
  zAssert( C++::zVec (shift), result );
}

void assert_vec_scale()
{
  const int size = 10;
  const int n = 100;
  zVec min  = new zVecStruct{ size };
  zVec max  = new zVecStruct{ size };
  zVec src  = new zVecStruct{ size };
  zVec dest = new zVecStruct{ size };
  double elem_min, elem_max;
  bool result = true;
  for(int i=0; i<n; i++ ){
    min->randUniform( -10, 10 );
    max->randUniform( -10, 10 );
    for(int j=0; j<min->size; j++ )
      if( min->buf[j] > max->buf[j] )
        zSwap( double, min->buf[j], max->buf[j] );
    src->randUniform( 0, 1 );
    dest->copy( src );
    dest->scale( min, max );
    zVecSubDRC( dest, min );
    zVecSubDRC( max, min );
    zVecDemDRC( dest, max );
    if( !zVecEqual( src, dest, zTOL ) ){
      zVecSubDRC( src, dest );
      zVecPrint( src );
      result = false;
    }
  }
  zAssert( C++::zVec (scale), result );
  result = true;
  for(int i=0; i<n; i++ ){
    elem_min = zRandF( -10,  0 );
    elem_max = zRandF(   0, 10 );
    src->randUniform( 0, 1 );
    dest->copy( src );
    dest->scaleUniform( elem_min, elem_max );
    dest->shift( -elem_min );
    zVecDivDRC( dest, ( elem_max - elem_min ) );
    if( !zVecEqual( src, dest, zTOL ) ){
      zVecSubDRC( src, dest );
      zVecPrint( src );
      result = false;
    }
  }
  zAssert( C++::zVec (scaleUniform), result );
  zVecFreeAtOnce( 4, min, max, src, dest );
}

void assert_vec_normalize()
{
  const int size = 10;
  zVec test_vec1 = new zVecStruct{ size };
  zVec test_vec2 = new zVecStruct{ size };
  test_vec1->randUniform( -10, 10 );
  double norm = test_vec1->norm();
  test_vec2->copy( test_vec1 );
  test_vec2->normalize();
  zAssert( C++::zVec (normalize), zIsTiny( test_vec2->norm() - 1 ) );
  bool result = true;
  for(int i=0; i<size; i++ )
    if( !zIsTiny( test_vec1->buf[i]/test_vec2->buf[i] - norm ) ) result = false;
  zAssert( C++::zVec (norm), result );
  zVecFreeAtOnce( 2, test_vec1, test_vec2 );
}

int main()
{
  assert_vec_setlist();
  assert_vec_clone();
  assert_vec_get_put();
  assert_vec_resize();
  assert_vec_misc();
  assert_vec_shift();
  assert_vec_scale();
  assert_vec_normalize();
  return 0;
}
