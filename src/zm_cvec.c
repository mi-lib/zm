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

/* allocate a complex vector. */
zCVec zCVecAlloc(int size)
{
  zCVec v;

  if( !( v = zAlloc( zCVecStruct, 1 ) ) ){
    ZALLOCERROR();
    return NULL;
  }
  if( !( zArrayBuf(v) = zAlloc( zComplex, size ) ) ){
    ZALLOCERROR();
    free( v );
    return NULL;
  }
  zCVecSetSizeNC( v, size );
  return v;
}

/* free a complex vector. */
void zCVecFree(zCVec v)
{
  if( !v ) return;
  zFree( zArrayBuf(v) );
  free( v );
}

/* zero a complex vector. */
zCVec zCVecZero(zCVec v)
{
  register int i;

  for( i=0; i<zCVecSizeNC(v); i++ )
    zComplexZero( zCVecElemNC( v, i ) );
  return v;
}

/* copy a complex vector without checking size consistency. */
zCVec zCVecCopyNC(zCVec src, zCVec dest)
{
  register int i;

  for( i=0; i<zCVecSizeNC(dest); i++ )
    zComplexCopy( zCVecElemNC(src,i), zCVecElemNC(dest,i) );
  return dest;
}

/* copy a complex vector. */
zCVec zCVecCopy(zCVec src, zCVec dest)
{
  return zCVecSizeIsEqual( src, dest ) ?
    zCVecCopyNC( src, dest ) : NULL;
}

/* clone a complex vector. */
zCVec zCVecClone(zCVec src)
{
  zCVec dest;

  if( ( dest = zCVecAlloc( zCVecSizeNC(src) ) ) )
    zCVecCopyNC( src, dest );
  return dest;
}

/* convert a real vector to a complex vector. */
zCVec zVec2CVec(zVec v, zCVec cv)
{
  register int i;

  for( i=0; i<zCVecSizeNC(cv); i++ )
    zComplexCreate( zCVecElemNC(cv,i), zVecElemNC(v,i), 0 );
  return cv;
}

/* check if two complex vectors are equal. */
bool zCVecIsEqual(zCVec v1, zCVec v2)
{
  register int i;
  zComplex e;

  if( !zCVecSizeIsEqual( v1, v2 ) ) return false;
  for( i=0; i<zCVecSizeNC(v1); i++ ){
    zComplexSub( zCVecElemNC(v1,i), zCVecElemNC(v2,i), &e );
    if( !zComplexIsTiny( &e ) ) return false;
  }
  return true;
}

/* check if a complex vector is tiny. */
bool zCVecIsTol(zCVec v, double tol)
{
  register int i;

  for( i=0; i<zCVecSizeNC(v); i++ )
    if( !zComplexIsTol( zCVecElemNC(v,i), tol ) ) return false;
  return true;
}

/* add two complex vectors without checking size consistency. */
zCVec zCVecAddNC(zCVec v1, zCVec v2, zCVec v)
{
  register int i;

  for( i=0; i<zCVecSizeNC(v); i++ )
    zComplexAdd( zCVecElemNC(v1,i), zCVecElemNC(v2,i), zCVecElemNC(v,i) );
  return v;
}

/* subtract a complex vector from another without checking size consistency. */
zCVec zCVecSubNC(zCVec v1, zCVec v2, zCVec v)
{
  register int i;

  for( i=0; i<zCVecSizeNC(v); i++ )
    zComplexSub( zCVecElemNC(v1,i), zCVecElemNC(v2,i), zCVecElemNC(v,i) );
  return v;
}

/* reverse a complex vector without checking size consistency. */
zCVec zCVecRevNC(zCVec v1, zCVec v)
{
  register int i;

  for( i=0; i<zCVecSizeNC(v); i++ )
    zComplexRev( zCVecElemNC(v1,i), zCVecElemNC(v,i) );
  return v;
}

/* multiply a complex vector by a complex scalar value
 * without checking size consistency. */
zCVec zCVecMulNC(zCVec v1, zComplex *z, zCVec v)
{
  register int i;

  for( i=0; i<zCVecSizeNC(v); i++ )
    zComplexCMul( zCVecElemNC(v1,i), z, zCVecElemNC(v,i) );
  return v;
}

/* divide a complex vector by a complex scalar value
 * without checking size consistency. */
zCVec zCVecDivNC(zCVec v1, zComplex *z, zCVec v)
{
  register int i;
  double r;
  zComplex dz;

  r = zComplexSqrAbs( z );
  zComplexConj( z, &dz );
  zComplexDiv( &dz, r, &dz );
  for( i=0; i<zCVecSizeNC(v); i++ )
    zComplexCMul( zCVecElemNC(v1,i), &dz, zCVecElemNC(v,i) );
  return v;
}

/* concatenate a complex vector with another multiplied by a
 * complex scalar value without checking size consistency. */
zCVec zCVecCatNC(zCVec v1, zComplex *z, zCVec v2, zCVec v)
{
  register int i;
  zComplex dz;

  for( i=0; i<zCVecSizeNC(v); i++ ){
    zComplexCMul( zCVecElemNC(v2,i), z, &dz );
    zComplexAdd( zCVecElemNC(v1,i), &dz, zCVecElemNC(v,i) );
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

/* add two complex vectors. */
zCVec zCVecAdd(zCVec v1, zCVec v2, zCVec v)
{
  __z_cvec_size_check_3( v1, v2, v );
  return zCVecAddNC( v1, v2, v );
}

/* subtract a complex vector from another. */
zCVec zCVecSub(zCVec v1, zCVec v2, zCVec v)
{
  __z_cvec_size_check_3( v1, v2, v );
  return zCVecSubNC( v1, v2, v );
}

/* reverse a complex vector. */
zCVec zCVecRev(zCVec v1, zCVec v)
{
  __z_cvec_size_check_2( v1, v );
  return zCVecRevNC( v1, v );
}

/* multiply a complex vector by a complex scalar value. */
zCVec zCVecMul(zCVec v1, zComplex *z, zCVec v)
{
  __z_cvec_size_check_2( v1, v );
  return zCVecMulNC( v1, z, v );
}

/* divide a complex vector by a complex scalar value. */
zCVec zCVecDiv(zCVec v1, zComplex *z, zCVec v)
{
  __z_cvec_size_check_2( v1, v );
  if( zComplexIsTiny( z ) ){
    ZRUNERROR( ZM_ERR_ZERODIV );
    return NULL;
  }
  return zCVecDivNC( v1, z, v );
}

/* concatenate a complex vector with another multiplied by a complex scalar value. */
zCVec zCVecCat(zCVec v1, zComplex *z, zCVec v2, zCVec v)
{
  __z_cvec_size_check_3( v1, v2, v );
  return zCVecCatNC( v1, z, v2, v );
}

/* inner product of two complex vectors without checking size consistency. */
zComplex *zCVecInnerProdNC(zCVec v1, zCVec v2, zComplex *z)
{
  register int i;
  zComplex c;

  zComplexZero( z );
  for( i=0; i<zCVecSizeNC(v1); i++ ){
    zComplexCMul( zCVecElemNC(v1,i), zCVecElemNC(v2,i), &c );
    zComplexAdd( z, &c, z );
  }
  return z;
}

/* inner product of two complex vector. */
zComplex *zCVecInnerProd(zCVec v1, zCVec v2, zComplex *z)
{
  if( !zCVecSizeIsEqual(v1,v2) ){
    ZRUNERROR( ZM_ERR_SIZMIS_VEC );
    return 0;
  }
  return zCVecInnerProdNC( v1, v2, z );
}

/* squared norm of a complex vector. */
double zCVecSqrNorm(zCVec v)
{
  zComplex z;

  return zComplexAbs( zCVecInnerProdNC( v, v, &z ) );
}

/* normalize a complex vector. */
zCVec zCVecNormalize(zCVec src, zCVec dest)
{
  register int i;
  double r;

  if( zIsTiny( ( r = zCVecNorm( src ) ) ) ){
    ZRUNWARN( ZM_ERR_ZERONORM );
    return NULL;
  }
  for( i=0; i<zCVecSizeNC(dest); i++ )
    zComplexDiv( zCVecElemNC(src,i), r, zCVecElemNC(dest,i) );
  return dest;
}

/* print a complex vector to a file. */
void zCVecFPrint(FILE *fp, zCVec v)
{
  register int i;

  if( !v )
    fprintf( fp, "(null vector)\n" );
  else{
    fprintf( fp, "%d (\n", zCVecSizeNC(v) );
    for( i=0; i<zCVecSizeNC(v); i++ ){
      fprintf( fp, "  " );
      zComplexFPrint( fp, zCVecElemNC(v,i) );
      fprintf( fp, "\n" );
    }
    fprintf( fp, ")\n" );
  }
}
