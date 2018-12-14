/* ZM - Z's Mathematics Toolbox
 * Copyright (C) 1998 Tomomichi Sugihara (Zhidao)
 *
 * zm_vec_ring - vector ring class.
 */

#include <zm/zm_vec.h>

/* zVecRingAlloc
 * - allocate vector ring.
 */
bool zVecRingAlloc(zVecRing *ring, int dim, int num)
{
  register int i;

  zRingAlloc( ring, zVec, num );
  for( i=0; i<num; i++ )
    if( !( *zRingElem(ring,i) = zVecAlloc( dim ) ) ){
      ZALLOCERROR();
      zVecRingFree( ring );
      return false;
    }
  return true;
}

/* zVecRingFree
 * - free vector ring.
 */
void zVecRingFree(zVecRing *ring)
{
  register int i;

  if( zRingBuf(ring) )
    for( i=0; i<zRingNum(ring); i++ )
      zVecFree( *zRingElem(ring,i) );
  zRingFree( ring );
}

/* zVecRingFill
 * - fill the ring with the same vector.
 */
zVec zVecRingFill(zVecRing *ring, zVec v)
{
  register int i;

  if( !zVecSizeIsEqual(*zRingHead(ring),v) ){
    ZRUNERROR( ZM_ERR_SIZMIS_VEC );
    return NULL;
  }
  for( i=0; i<zRingNum(ring); i++ )
    zVecCopyNC( v, *zRingElem(ring,i) );
  return v;
}

/* zVecRingCat
 * - concatenate vector ring.
 */
zVec zVecRingCat(zVec v, zVec c, zVecRing *ring)
{
  register int i;

  if( !zVecSizeIsEqual(*zRingHead(ring),v) ){
    ZRUNERROR( ZM_ERR_SIZMIS_VEC );
    return NULL;
  }
  for( i=0; i<zVecSizeNC(c); i++ )
    zVecCatNCDRC( v, zVecElem(c,i), *zRingElem(ring,i) );
  return v;
}

/* zVecRingLS
 * - linear sum of vector ring.
 */
zVec zVecRingLS(zVec v, zVec c, zVecRing *ring)
{
  zVecClear( v );
  return zVecRingCat( v, c, ring );
}
