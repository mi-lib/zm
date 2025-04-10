/* ZM - Z's Mathematics Toolbox
 * Copyright (C) 1998 Tomomichi Sugihara (Zhidao)
 *
 * zm_vec_ring - vector ring class.
 */

#include <zm/zm_vec.h>

/* allocate a vector ring. */
bool zVecRingAlloc(zVecRing *ring, int dim, int num)
{
  int i;

  zRingAlloc( ring, zVec, num );
  for( i=0; i<num; i++ )
    if( !( *zRingElem(ring,i) = zVecAlloc( dim ) ) ){
      ZALLOCERROR();
      zVecRingFree( ring );
      return false;
    }
  return true;
}

/* free a vector ring. */
void zVecRingFree(zVecRing *ring)
{
  int i;

  if( zRingBuf(ring) )
    for( i=0; i<zRingSize(ring); i++ )
      zVecFree( *zRingElem(ring,i) );
  zRingFree( ring );
}

/* fill a ring with the same vector. */
zVec zVecRingFill(zVecRing *ring, zVec v)
{
  int i;

  if( !zVecSizeEqual(*zRingHead(ring),v) ){
    ZRUNERROR( ZM_ERR_VEC_SIZEMISMATCH );
    return NULL;
  }
  for( i=0; i<zRingSize(ring); i++ )
    zVecCopyNC( v, *zRingElem(ring,i) );
  return v;
}

/* concatenate vectors of a ring. */
zVec zVecRingCat(zVec v, const zVec c, const zVecRing *ring)
{
  int i;

  if( !zVecSizeEqual(*zRingHead(ring),v) ){
    ZRUNERROR( ZM_ERR_VEC_SIZEMISMATCH );
    return NULL;
  }
  for( i=0; i<zVecSizeNC(c); i++ )
    zVecCatNCDRC( v, zVecElem(c,i), *zRingElem(ring,i) );
  return v;
}

/* linear sum of vectors of a ring. */
zVec zVecRingLinearSum(zVec v, const zVec c, const zVecRing *ring)
{
  zVecZero( v );
  return zVecRingCat( v, c, ring );
}
