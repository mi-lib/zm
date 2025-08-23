/* ZM - Z's Mathematics Toolbox
 * Copyright (C) 1998 Tomomichi Sugihara (Zhidao)
 *
 * zm_vec_array - vector array class.
 */

#include <zm/zm_vec.h>

/* allocate a vector array. */
bool zVecArrayAlloc(zVecArray *array, int dim, int num)
{
  int i;

  zArrayAlloc( array, zVec, num );
  for( i=0; i<num; i++ )
    if( !( *zArrayElem(array,i) = zVecAlloc( dim ) ) ){
      ZALLOCERROR();
      zVecArrayFree( array );
      return false;
    }
  return true;
}

/* free vector array. */
void zVecArrayFree(zVecArray *array)
{
  int i;

  if( zArrayBuf(array) )
    for( i=0; i<zIndexSizeNC(array); i++ )
      zVecFree( *zArrayElem(array,i) );
  zArrayFree( array );
}

/* fill the array with the same vector. */
zVec zVecArrayFill(zVecArray *array, zVec v)
{
  int i;

  if( !zVecSizeEqual(*zArrayHead(array),v) ){
    ZRUNERROR( ZM_ERR_VEC_SIZEMISMATCH );
    return NULL;
  }
  for( i=0; i<zIndexSizeNC(array); i++ )
    zVecCopyNC( v, *zArrayElem(array,i) );
  return v;
}

/* create zVec from a zVecArray */
zVec *zVecCreateFromzVecArray(const zVecArray *src_array, zVec *dest_vec)
{
  int i, vec_size, array_size;

  vec_size = zVecSize( *zArrayElemNC(src_array, 0) );
  array_size = zArraySize(src_array);
  *dest_vec = zVecAlloc( vec_size * array_size );
  for( i=0; i<array_size; i++ )
    zVecPutNC( *dest_vec, i * vec_size, *zArrayElemNC( src_array, i ) );
  return dest_vec;
}

/* create zVecArray from a zVec */
zVecArray *zVecArrayCreateFromzVec(const zVec src_vec, const int vec_size, const int array_size, zVecArray *dest_array)
{
  int i;

  if( vec_size * array_size != zVecSize(src_vec) ){
    ZRUNERROR( ZM_ERR_VEC_SIZEMISMATCH );
    return NULL;
  }
  zVecArrayAlloc( dest_array, vec_size, array_size );
  for( i=0; i<array_size; i++ )
    zVecGetNC( src_vec, vec_size * i, *zArrayElemNC(dest_array, i) );
  return dest_array;
}

