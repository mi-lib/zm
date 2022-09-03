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
    for( i=0; i<zArraySize(array); i++ )
      zVecFree( *zArrayElem(array,i) );
  zArrayFree( array );
}

/* fill the array with the same vector. */
zVec zVecArrayFill(zVecArray *array, zVec v)
{
  int i;

  if( !zVecSizeIsEqual(*zArrayHead(array),v) ){
    ZRUNERROR( ZM_ERR_SIZMIS_VEC );
    return NULL;
  }
  for( i=0; i<zArraySize(array); i++ )
    zVecCopyNC( v, *zArrayElem(array,i) );
  return v;
}
