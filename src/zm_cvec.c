/* ZM - Z's Mathematics Toolbox
 * Copyright (C) 1998 Tomomichi Sugihara (Zhidao)
 *
 * zm_cvec - complex vector class.
 */

#include <zm/zm_cvec.h>

/* ********************************************************** */
/* CLASS: zCVec
 * double precision floating point value vector class
 * ********************************************************** */

/* zCVecAlloc
 * - allocate vector.
 */
zCVec zCVecAlloc(int size)
{
  zCVec v;

  if( !( v = zAlloc( zCVecStruct, 1 ) ) ){
    ZALLOCERROR();
    return NULL;
  }
  if( !( v->elem = zAlloc( zComplex, size ) ) ){
    ZALLOCERROR();
    free( v );
    return NULL;
  }
  zCVecSetSize( v, size );
  return v;
}

/* zCVecFree
 * - free vector.
 */
void zCVecFree(zCVec v)
{
  if( !v ) return;
  zFree( v->elem );
  free( v );
}

/* zCVecClear
 * - cleanup vector.
 */
zCVec zCVecClear(zCVec v)
{
  register int i;

  for( i=0; i<_zCVecSize(v); i++ )
    zCVecSetElem( v, i, ZCOMPLEXZERO );
  return v;
}

/* zCVecCopyNC
 * - copy vector without checking size consistency.
 */
zCVec zCVecCopyNC(zCVec src, zCVec dest)
{
  register int i;

  for( i=0; i<_zCVecSize(dest); i++ )
    zComplexCopy( zCVecElem(src,i), zCVecElem(dest,i) );
  return dest;
}

/* zCVecCopy
 * - copy vector.
 */
zCVec zCVecCopy(zCVec src, zCVec dest)
{
  return zCVecSizeIsEqual( src, dest ) ?
    zCVecCopyNC( src, dest ) : NULL;
}

/* zCVecClone
 * - clone vector.
 */
zCVec zCVecClone(zCVec src)
{
  zCVec dest;

  if( ( dest = zCVecAlloc( zCVecSize(src) ) ) )
    zCVecCopyNC( src, dest );
  return dest;
}

/* zVec2CVec
 * - convert real vector to complex vector.
 */
zCVec zVec2CVec(zVec v, zCVec cv)
{
  register int i;

  for( i=0; i<_zCVecSize(cv); i++ )
    zComplexCreate( zCVecElem(cv,i), zVecElem(v,i), 0 );
  return cv;
}

/* zCVecIsEqual
 * - see if two vectors are equal.
 */
bool zCVecIsEqual(zCVec v1, zCVec v2)
{
  register int i;
  zComplex e;

  if( !zCVecSizeIsEqual( v1, v2 ) ) return false;
  for( i=0; i<_zCVecSize(v1); i++ ){
    zComplexSub( zCVecElem(v1,i), zCVecElem(v2,i), &e );
    if( !zComplexIsTiny( &e ) ) return false;
  }
  return true;
}

/* zCVecIsTol
 * - test for tiny vector.
 */
bool zCVecIsTol(zCVec v, double tol)
{
  register int i;

  for( i=0; i<_zCVecSize(v); i++ )
    if( !zComplexIsTol( zCVecElem(v,i), tol ) ) return false;
  return true;
}

/* zCVecAddNC
 * - add two vectors without checking size consistency.
 */
zCVec zCVecAddNC(zCVec v1, zCVec v2, zCVec v)
{
  register int i;

  for( i=0; i<_zCVecSize(v); i++ )
    zComplexAdd( zCVecElem(v1,i), zCVecElem(v2,i), zCVecElem(v,i) );
  return v;
}

/* zCVecSubNC
 * - subtract vectors without checking size consistency.
 */
zCVec zCVecSubNC(zCVec v1, zCVec v2, zCVec v)
{
  register int i;

  for( i=0; i<_zCVecSize(v); i++ )
    zComplexSub( zCVecElem(v1,i), zCVecElem(v2,i), zCVecElem(v,i) );
  return v;
}

/* zCVecRevNC
 * - reverse vector without checking size consistency.
 */
zCVec zCVecRevNC(zCVec v1, zCVec v)
{
  register int i;

  for( i=0; i<_zCVecSize(v); i++ )
    zComplexRev( zCVecElem(v1,i), zCVecElem(v,i) );
  return v;
}

/* zCVecMulNC
 * - multiply vector by double value without checking
 *   size consistency.
 */
zCVec zCVecMulNC(zCVec v1, zComplex *z, zCVec v)
{
  register int i;

  for( i=0; i<_zCVecSize(v); i++ )
    zComplexCMul( zCVecElem(v1,i), z, zCVecElem(v,i) );
  return v;
}

/* zCVecDivNC
 * - divide vector by double value without checking size
 *   consistency.
 */
zCVec zCVecDivNC(zCVec v1, zComplex *z, zCVec v)
{
  register int i;
  double r;
  zComplex dz;

  r = zComplexSqrAbs( z );
  zComplexConj( z, &dz );
  zComplexDiv( &dz, r, &dz );
  for( i=0; i<_zCVecSize(v); i++ )
    zComplexCMul( zCVecElem(v1,i), &dz, zCVecElem(v,i) );
  return v;
}

/* zCVecCatNC
 * - concatenate vector by a complex value without checking size
 *   consistency.
 */
zCVec zCVecCatNC(zCVec v1, zComplex *z, zCVec v2, zCVec v)
{
  register int i;
  zComplex dz;

  for( i=0; i<_zCVecSize(v); i++ ){
    zComplexCMul( zCVecElem(v2,i), z, &dz );
    zComplexAdd( zCVecElem(v1,i), &dz, zCVecElem(v,i) );
  }
  return v;
}

#define __z_cvec_size_check_2(v1,v2) \
  if( !zCVecSizeIsEqual(v1,v2) ){\
    ZRUNERROR( ZM_ERR_SIZMIS_VEC );\
    return NULL;\
  }
#define __z_cvec_size_check_3(v1,v2,v) \
  if( !zCVecSizeIsEqual(v1,v2) || !zCVecSizeIsEqual(v1,v) ){\
    ZRUNERROR( ZM_ERR_SIZMIS_VEC );\
    return NULL;\
  }

/* zCVecAdd
 * - add two vectors.
 */
zCVec zCVecAdd(zCVec v1, zCVec v2, zCVec v)
{
  __z_cvec_size_check_3( v1, v2, v );
  return zCVecAddNC( v1, v2, v );
}

/* zCVecSub
 * - subtract vectors.
 */
zCVec zCVecSub(zCVec v1, zCVec v2, zCVec v)
{
  __z_cvec_size_check_3( v1, v2, v );
  return zCVecSubNC( v1, v2, v );
}

/* zCVecRev
 * - reverse vector.
 */
zCVec zCVecRev(zCVec v1, zCVec v)
{
  __z_cvec_size_check_2( v1, v );
  return zCVecRevNC( v1, v );
}

/* zCVecMul
 * - multiply vector by value.
 */
zCVec zCVecMul(zCVec v1, zComplex *z, zCVec v)
{
  __z_cvec_size_check_2( v1, v );
  return zCVecMulNC( v1, z, v );
}

/* zCVecDiv
 * - divide vector by value.
 */
zCVec zCVecDiv(zCVec v1, zComplex *z, zCVec v)
{
  __z_cvec_size_check_2( v1, v );
  if( zComplexIsTiny( z ) ){
    ZRUNERROR( ZM_ERR_ZERODIV );
    return NULL;
  }
  return zCVecDivNC( v1, z, v );
}

/* zCVecCat
 * - concatenate vector by a complex value.
 */
zCVec zCVecCat(zCVec v1, zComplex *z, zCVec v2, zCVec v)
{
  __z_cvec_size_check_3( v1, v2, v );
  return zCVecCatNC( v1, z, v2, v );
}

/* zCVecInnerProdNC
 * - inner products of vector without checking size consistency.
 */
zComplex *zCVecInnerProdNC(zCVec v1, zCVec v2, zComplex *z)
{
  register int i;
  zComplex c;

  zComplexClear( z );
  for( i=0; i<_zCVecSize(v1); i++ ){
    zComplexCMul( zCVecElem(v1,i), zCVecElem(v2,i), &c );
    zComplexAdd( z, &c, z );
  }
  return z;
}

/* zCVecInnerProd
 * - inner products of vector.
 */
zComplex *zCVecInnerProd(zCVec v1, zCVec v2, zComplex *z)
{
  if( !zCVecSizeIsEqual(v1,v2) ){
    ZRUNERROR( ZM_ERR_SIZMIS_VEC );
    return 0;
  }
  return zCVecInnerProdNC( v1, v2, z );
}

/* zCVecSqrNorm
 * - squared norm of vector.
 */
double zCVecSqrNorm(zCVec v)
{
  zComplex z;

  return zComplexAbs( zCVecInnerProdNC( v, v, &z ) );
}

/* zCVecNormalize
 * - normalization of vector.
 */
zCVec zCVecNormalize(zCVec src, zCVec dest)
{
  register int i;
  double r;

  if( zIsTiny( ( r = zCVecNorm( src ) ) ) ){
    ZRUNWARN( ZM_ERR_ZERONORM );
    return NULL;
  }
  for( i=0; i<_zCVecSize(dest); i++ )
    zComplexDiv( zCVecElem(src,i), r, zCVecElem(dest,i) );
  return dest;
}

/* zCVecFWrite
 * - output of vector to file.
 */
void zCVecFWrite(FILE *fp, zCVec v)
{
  register int i;

  if( !v )
    fprintf( fp, "(null vector)\n" );
  else{
    fprintf( fp, "%d (\n", _zCVecSize(v) );
    for( i=0; i<_zCVecSize(v); i++ ){
      fprintf( fp, "  " );
      zComplexFWrite( fp, zCVecElem(v,i) );
      fprintf( fp, "\n" );
    }
    fprintf( fp, ")\n" );
  }
}
